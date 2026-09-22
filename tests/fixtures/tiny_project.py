#!/usr/bin/env python3
"""A tiny synthetic project, and the expected answers for Stages 2 to 5.

    tiny_project.py build   <dest> <repo> [seq_type]
    tiny_project.py verify  <dest>
    tiny_project.py qualify <dest>

`build` writes a reference and Stage 1 outputs for 8 samples. Stages 0 and 1
need salmon and fastp, so their products are written here directly; everything
downstream then runs for real against them. `seq_type` defaults to tagseq, which
is the branch that takes countsFromAbundance="no"; pass `rnaseq` for the bulk
branch, which takes "lengthScaledTPM".

`verify` recomputes the expected gene counts from the fixture itself, by the
route tximport is not taking (parse the GTF, sum NumReads per gene), and checks
what the stages produced against it. An assertion that just re-ran the pipeline's
own arithmetic would confirm nothing.

`qualify` is V04, "numerical agreement with an independent trusted route within
declared tolerances", for the bulk branch. It reimplements lengthScaledTPM here,
in Python, from the definition, and compares cell by cell.

Deterministic: one seed, fixed below, so a regenerated fixture is the same
fixture. Nothing here is random at check time.
"""
import csv
import json
import math
import random
import re
import shutil
import sys
from pathlib import Path

SEED = 20260922
N_FLAT = 180             # genes with no treatment effect
N_DE = 10                # genes up, and the same number down
BASE = 500               # reads per gene per sample before effect and noise
EFFECT = 8.0             # fold change carried by the DE genes
NOISE = 0.08             # within-group spread

CONTROL = ["S01", "S02", "S03", "S04"]
TREATED = ["S05", "S06", "S07", "S08"]
SAMPLES = CONTROL + TREATED

# Two samples are unusable, one by each route the QC policy has. Both must end
# up excluded from the count matrix, with a reason, and still appear in the
# disposition table.
LOW_MAPPING = "S03"      # mapping rate under sample_qc.mapping_rate_min
LOW_DEPTH = "S07"        # mapping fine, too few reads land on genes

CONTRAST = "treatment__treated_vs_control"
BARRED_TX = "T002|ENSEMBL|extra"    # a quant.sf name tximport must cut at "|"

# G001's two isoforms are deliberately far apart in length, and which one is
# used switches between the groups at constant read depth. That is the only
# configuration in which lengthScaledTPM differs from raw counts at all: for a
# single-transcript gene whose effective length is the same in every sample, the
# correction cancels exactly. Without this the bulk branch and the tag branch
# would agree everywhere and neither would be tested.
SHORT_ISOFORM, LONG_ISOFORM = "T001a", "T001b"
EFFLEN = {SHORT_ISOFORM: 500.0, LONG_ISOFORM: 5000.0, "T002": 1200.0}
SWITCH = {"control": (0.80, 0.20), "treated": (0.20, 0.80)}


def genes():
    """(gene_id, [transcript ids], treated-group multiplier)."""
    out = [
        # Two transcripts, so gene counts must be their sum. This is the one
        # thing tx2gene aggregation can get wrong without looking wrong.
        ("G001", [SHORT_ISOFORM, LONG_ISOFORM], 1.0),
        # Its quant.sf name carries a "|" suffix. If ignoreAfterBar does not
        # apply, this gene goes missing rather than erroring.
        ("G002", ["T002"], 1.0),
    ]
    out += [(f"UP{i:03d}", [f"UP{i:03d}_t1"], EFFECT) for i in range(N_DE)]
    out += [(f"DN{i:03d}", [f"DN{i:03d}_t1"], 1.0 / EFFECT) for i in range(N_DE)]
    out += [(f"G{i + 100:03d}", [f"G{i + 100:03d}_t1"], 1.0) for i in range(N_FLAT)]
    return out


def efflen(tx, i):
    """Effective length. Spread across the single-transcript genes, so nothing
    downstream can be right only because every length was the same number."""
    return EFFLEN.get(tx, 400.0 + 25.0 * (i % 37))


def quant_name(tx):
    return BARRED_TX if tx == "T002" else tx


def read_quant(path):
    """quant.sf -> {tx_id: (num_reads, tpm, effective_length)}, bar stripped."""
    out = {}
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            out[row["Name"].split("|")[0]] = (float(row["NumReads"]),
                                              float(row["TPM"]),
                                              float(row["EffectiveLength"]))
    return out


def read_tx2gene(dest):
    """GTF -> {tx: gene}, by plain text, never through txdbmaker."""
    out = {}
    for line in (dest / "reference" / "annotation.gtf").read_text().splitlines():
        attrs = line.split("\t")[8]
        out[re.search(r'transcript_id "([^"]+)"', attrs).group(1)] = \
            re.search(r'gene_id "([^"]+)"', attrs).group(1)
    return out


# --------------------------------------------------------------------------
def build(dest: Path, repo: Path, seq_type: str = "tagseq"):
    rng = random.Random(SEED)
    gs = genes()
    (dest / "config").mkdir(parents=True, exist_ok=True)
    ref = dest / "reference"
    ref.mkdir(exist_ok=True)

    lengths = {}
    gtf, fa = [], []
    for i, (gid, txs, _) in enumerate(gs):
        for j, tx in enumerate(txs):
            lengths[tx] = efflen(tx, i)
            span = int(lengths[tx]) + 150            # genomic span > effective
            start = 1 + j * 20000
            gtf.append(f"chr1\tfixture\texon\t{start}\t{start + span - 1}\t.\t+\t."
                       f'\tgene_id "{gid}"; transcript_id "{tx}";')
            fa.append(f">{quant_name(tx)}\n" + "ACGT" * (span // 4))
    (ref / "annotation.gtf").write_text("\n".join(gtf) + "\n")
    (ref / "transcriptome.fa").write_text("\n".join(fa) + "\n")

    # Stand-ins for what Stage 0 publishes, so Snakemake does not schedule a
    # reference build to satisfy inputs of a Stage 1 whose outputs exist.
    (ref / "salmon_idx").mkdir(exist_ok=True)
    (ref / "salmon_idx" / "info.json").write_text('{"fixture": true}\n')
    (ref / "reference.lock.json").write_text('{"fixture": true}\n')

    for s in SAMPLES:
        d = dest / "results" / "quant" / s
        d.mkdir(parents=True, exist_ok=True)
        group = "treated" if s in TREATED else "control"

        reads = {}
        for gid, txs, mult in gs:
            eff = mult if s in TREATED else 1.0
            for k, tx in enumerate(txs):
                # G001 keeps its total depth and moves it between a short and a
                # long isoform. Raw counts stay flat; the molecule count does
                # not, which is the bias lengthScaledTPM exists to remove.
                share = SWITCH[group][k] if gid == "G001" else 1.0
                reads[tx] = max(1, int(round(BASE * eff * share
                                             * (1.0 + rng.uniform(-NOISE, NOISE)))))

        # TPM as Salmon defines it, so the abundance column is consistent with
        # the counts and lengths rather than a placeholder. tximport reads this
        # column as abundance, and lengthScaledTPM is computed from it.
        rate = {tx: n / lengths[tx] for tx, n in reads.items()}
        denom = sum(rate.values())
        rows = ["Name\tLength\tEffectiveLength\tTPM\tNumReads"]
        for tx in reads:
            rows.append(f"{quant_name(tx)}\t{int(lengths[tx]) + 150}"
                        f"\t{lengths[tx]:.1f}\t{1e6 * rate[tx] / denom:.6f}"
                        f"\t{reads[tx]}")
        (d / "quant.sf").write_text("\n".join(rows) + "\n")

        # Depth is quoted against the floor for this assay: sample_qc asks for
        # 3M reads on genes for tagseq and 20M for bulk. A single set of numbers
        # would pass one branch and exclude the whole cohort on the other.
        deep, shallow = (5_000_000, 2_000_000) if seq_type == "tagseq" \
                        else (30_000_000, 12_000_000)
        map_rate, n_reads = {LOW_MAPPING: (0.45, deep),
                             LOW_DEPTH: (0.82, shallow)}.get(s, (0.88, deep))
        (d / "metrics.json").write_text(json.dumps({
            "sample": s, "seq_type": seq_type,
            "num_reads_processed": n_reads, "mapping_rate": map_rate,
            "min_map_rate_threshold": 0.60, "flagged_low_mapping": map_rate < 0.60,
            "detected_libtype": "SF", "expected_libtype": "",
            "libtype_mismatch": False, "salmon_bias_flags": "",
            "salmon_version": "fixture", "fastp_version": "fixture",
            "fastp_report": "",
        }, indent=2) + "\n")

    sheet = ["sample_id\tfastq_url\ttreatment"] + [
        f"{s}\t/dev/null/{s}.fq.gz\t" + ("treated" if s in TREATED else "control")
        for s in SAMPLES]
    (dest / "config" / "samples.tsv").write_text("\n".join(sheet) + "\n")

    # seq_type tagseq makes countsFromAbundance "no", the branch whose answer is
    # exactly predictable: gene counts are summed NumReads. rnaseq makes it
    # "lengthScaledTPM", which section 5.7 qualifies numerically.
    # A placeholder rather than an f-string: the YAML flow mappings below are
    # full of braces, and every one would have to be doubled.
    (dest / "config" / "config.yaml").write_text("""\
project: {name: fixture, description: tiny synthetic cohort, output_dir: "results/"}
organism:
  common_name: fixture
  scientific_name: Fixtura synthetica
  tax_id: 1
  kegg_code: null
  orgdb_package: org.Fsynthetica.eg.db
reference: {accession: FIXTURE_1.0, assembly_name: FIX1, cache_dir: null}
samples: {sheet: config/samples.tsv, seq_type: SEQ_TYPE}
model: {fixed_effects: "~ treatment", random_effects: null, primary_factor: treatment}
contrasts:
  - id: "treatment"
    type: "pairwise"
    factor: "treatment"
    reverse: true
sample_policy: {exclude_failing_qc: true}
downstream: {run_go: false, run_kegg: false, run_wgcna: true}
orgdb: {strategy: skip}
hpc: {delete_fastq_after_quant: false, samples_in_flight: null}
""".replace("SEQ_TYPE", seq_type))
    shutil.copy(repo / "config" / "thresholds.yaml",
                dest / "config" / "thresholds.yaml")


# --------------------------------------------------------------------------
def verify(dest: Path):
    fails = []

    def ck(cond, msg):
        print(("  ok   " if cond else "  FAIL ") + msg)
        if not cond:
            fails.append(msg)

    def tsv(p):
        with open(dest / p) as f:
            return list(csv.DictReader(f, delimiter="\t"))

    def js(p):
        return json.loads((dest / p).read_text())

    # Expectation rebuilt from the fixture, independently of tximport.
    tx2gene = {}
    for line in (dest / "reference" / "annotation.gtf").read_text().splitlines():
        attrs = line.split("\t")[8]
        tx2gene[re.search(r'transcript_id "([^"]+)"', attrs).group(1)] = \
            re.search(r'gene_id "([^"]+)"', attrs).group(1)

    exp = {}
    for s in SAMPLES:
        for row in tsv(f"results/quant/{s}/quant.sf"):
            g = tx2gene[row["Name"].split("|")[0]]
            exp.setdefault(g, {}).setdefault(s, 0)
            exp[g][s] += int(row["NumReads"])

    with open(dest / "results" / "counts.tsv") as f:
        r = csv.DictReader(f, delimiter="\t")
        cols = [c for c in r.fieldnames if c != "gene_id"]
        got = {row["gene_id"]: {c: int(row[c]) for c in cols} for row in r}

    print("== Stage 2: gene-level aggregation ==")
    kept = [s for s in SAMPLES if s not in (LOW_MAPPING, LOW_DEPTH)]
    ck(cols == kept, f"cohort after QC exclusion is {cols}")
    ck(set(got) == set(exp), f"all {len(exp)} genes present, none invented")
    bad = [(g, c) for g in got for c in cols if got[g][c] != exp[g][c]]
    ck(not bad, f"every gene x sample count equals the summed NumReads "
                f"({len(got) * len(cols)} cells)"
                + ("" if not bad else f"; first mismatch {bad[0]}"))

    q = {r["Name"]: int(r["NumReads"]) for r in tsv("results/quant/S01/quant.sf")}
    ck(got["G001"]["S01"] == q["T001a"] + q["T001b"] and q["T001a"] != q["T001b"],
       f"G001 is the SUM of its transcripts ({q['T001a']}+{q['T001b']}"
       f"={got['G001']['S01']}), not either one")
    ck("G002" in got and got["G002"]["S01"] == q[BARRED_TX],
       "G002 survived the '|' in its quant.sf name (ignoreAfterBar)")

    print("== Stage 2: QC disposition ==")
    disp = {r["sample_id"]: r for r in tsv("results/sample_disposition.tsv")}
    ck(len(disp) == len(SAMPLES),
       f"all {len(SAMPLES)} samples have a disposition row, excluded ones included")
    ck(disp[LOW_MAPPING]["status"] == "flagged"
       and "mapping rate" in disp[LOW_MAPPING]["flag"]
       and disp[LOW_MAPPING]["include"] == "FALSE",
       f"{LOW_MAPPING} excluded on mapping rate: {disp[LOW_MAPPING]['flag']}")
    ck(disp[LOW_DEPTH]["status"] == "flagged"
       and "mapped reads" in disp[LOW_DEPTH]["flag"]
       and disp[LOW_DEPTH]["include"] == "FALSE",
       f"{LOW_DEPTH} excluded on depth: {disp[LOW_DEPTH]['flag']}")
    ck(all(disp[s]["status"] == "ok" for s in cols),
       f"the {len(cols)} retained samples are all status=ok")

    print("== Stage 2: count provenance (CT01-CT03) ==")
    prov = js("results/counts_provenance.json")
    ck(prov["counts_from_abundance"] == "no", "tagseq took countsFromAbundance='no'")
    ck(prov["length_correction_applied"] is False, "no length correction was applied")
    ck(prov["further_length_correction_permitted"] is False,
       "a further correction is refused (CT03)")

    print("== Stage 3: differential expression ==")
    de = {r["gene_id"]: r for r in
          tsv(f"results/DE_Results/{CONTRAST}_DE_analysis.tsv")}
    ck(len(de) == len(exp), f"{len(de)} genes tested")
    up = [float(de[f"UP{i:03d}"]["logFC"]) for i in range(N_DE)]
    dn = [float(de[f"DN{i:03d}"]["logFC"]) for i in range(N_DE)]
    flat = [float(de[g]["logFC"]) for g in de if g.startswith("G")]
    # The sign is the point: a flipped contrast would put these the other way
    # round and every biological conclusion drawn from it would invert.
    ck(all(x > 2.5 for x in up),
       f"all {N_DE} UP genes have logFC > 2.5 (min {min(up):.2f}, "
       f"expected ~{math.log2(EFFECT):.2f})")
    ck(all(x < -2.5 for x in dn),
       f"all {N_DE} DN genes have logFC < -2.5 (max {max(dn):.2f})")
    ck(max(abs(x) for x in flat) < 1.0,
       f"the {len(flat)} flat genes stay under |logFC| 1.0 "
       f"(max {max(abs(x) for x in flat):.2f})")
    sig = [g for g in de if g.startswith("UP") and float(de[g]["adj.P.Val"]) < 0.05]
    ck(len(sig) == N_DE, f"{len(sig)}/{N_DE} UP genes are significant at FDR 5%")
    d = {r["contrast_name"]: r for r in
         tsv("results/DE_Results/contrast_directions.tsv")}
    ck(d[CONTRAST]["numerator_level"] == "treated"
       and d[CONTRAST]["denominator_level"] == "control",
       "the contrast name matches its recorded direction (treated - control)")

    print("== Stage 4: enrichment status taxonomy ==")
    st = tsv("results/enrichment_status.tsv")
    ck(len(st) == 2 and {r["method"] for r in st} == {"GO", "KEGG"},
       "one row per method")
    # Both are switched off in the config. Master plan 11 keeps a module the
    # user turned off apart from one that errored; reporting these as failed
    # made a policy skip end the stage.
    ck(all(r["status"] == "skipped_by_policy" for r in st),
       "both are skipped_by_policy, not failed (master plan 11)")
    em = js("metrics/enrichment.json")
    ck(em["n_failed"] == 0 and em["n_skipped"] == 2,
       f"metrics agree: {em['n_failed']} failed, {em['n_skipped']} skipped")

    print("== Stage 5: WGCNA gate ==")
    wm = js("metrics/wgcna.json")
    ck(wm.get("skipped") is True
       and f"n_samples={len(cols)}" in wm["reason"] and "min_samples=15" in wm["reason"],
       f"skipped with a stated reason: {wm['reason']}")

    print(f"\n{len(fails)} failed")
    return 1 if fails else 0


# --------------------------------------------------------------------------
# V04: "numerical agreement with an independent trusted route within declared
# tolerances" (master plan 12). The route under test is Snakemake -> 02_aggregate.R
# -> tx2gene built by txdbmaker -> tximport. The route here is this file: the
# GTF parsed as text, and lengthScaledTPM evaluated from its definition, in a
# different language. Agreement means the pipeline wired tximport up correctly
# (right mapping, right countsFromAbundance, right columns, right rounding);
# disagreement is a wiring defect, because both sides read the same quant.sf.
#
# TOLERANCE, declared: every cell must be within 0.5 of the independent value.
# 02_aggregate.R writes round(), so half a count is the entire difference
# rounding is allowed to make, and nothing else may differ at all.
TOLERANCE = 0.5 + 1e-6


def length_scaled_tpm(dest, samples):
    """counts[gene][sample], by the definition rather than by calling tximport.

        abundance[g,s] = sum of transcript TPM
        counts[g,s]    = sum of transcript NumReads
        length[g,s]    = sum(TPM * efflen) / abundance          (TPM-weighted)
        new[g,s]       = abundance[g,s] * mean over samples of length[g,s]
        out[g,s]       = new[g,s] * colsum(counts)[s] / colsum(new)[s]

    The final rescale is what makes each column sum back to the raw read total,
    which is the invariant asserted separately below.
    """
    tx2gene = read_tx2gene(dest)
    abundance, counts, wlength = {}, {}, {}
    for s in samples:
        for tx, (n, tpm, ln) in read_quant(
                dest / "results" / "quant" / s / "quant.sf").items():
            g = tx2gene[tx]
            abundance.setdefault(g, {}).setdefault(s, 0.0)
            counts.setdefault(g, {}).setdefault(s, 0.0)
            wlength.setdefault(g, {}).setdefault(s, 0.0)
            abundance[g][s] += tpm
            counts[g][s] += n
            wlength[g][s] += tpm * ln

    # Every gene in this fixture has non-zero abundance in every sample, so
    # tximport's replaceMissingLength() fallback is never reached. It is
    # therefore not qualified here.
    for g in abundance:
        for s in samples:
            if abundance[g][s] <= 0:
                raise AssertionError(
                    f"{g}/{s} has zero abundance; this fixture does not cover "
                    "tximport's missing-length fallback")

    mean_len = {g: sum(wlength[g][s] / abundance[g][s] for s in samples) / len(samples)
                for g in abundance}
    new = {g: {s: abundance[g][s] * mean_len[g] for s in samples} for g in abundance}
    scale = {s: sum(counts[g][s] for g in counts) / sum(new[g][s] for g in new)
             for s in samples}
    return ({g: {s: new[g][s] * scale[s] for s in samples} for g in new},
            counts)


def qualify(dest: Path):
    fails = []

    def ck(cond, msg):
        print(("  ok   " if cond else "  FAIL ") + msg)
        if not cond:
            fails.append(msg)

    prov = json.loads((dest / "results" / "counts_provenance.json").read_text())
    print("== the bulk branch was taken (CT01, CT02) ==")
    ck(prov["counts_from_abundance"] == "lengthScaledTPM",
       f"countsFromAbundance is {prov['counts_from_abundance']!r}")
    ck(prov["length_correction_applied"] is True, "a length correction was applied")
    ck(prov["further_length_correction_permitted"] is False,
       "a second correction is still refused (CT02)")

    with open(dest / "results" / "counts.tsv") as f:
        r = csv.DictReader(f, delimiter="\t")
        cols = [c for c in r.fieldnames if c != "gene_id"]
        got = {row["gene_id"]: {c: float(row[c]) for c in cols} for row in r}

    # The independent route sees every sample; the pipeline's matrix has had the
    # QC-excluded ones removed, so compare on the columns that survived.
    exp, raw = length_scaled_tpm(dest, SAMPLES)

    print(f"== V04: agreement with an independent route (tolerance {TOLERANCE}) ==")
    ck(set(got) == set(exp), f"both routes produced the same {len(exp)} genes")
    diffs = [(g, c, abs(got[g][c] - exp[g][c])) for g in got for c in cols]
    worst = max(diffs, key=lambda d: d[2])
    ck(worst[2] <= TOLERANCE,
       f"every one of {len(diffs)} cells agrees within {TOLERANCE}; "
       f"worst is {worst[0]}/{worst[1]} at {worst[2]:.4f}")
    # The declared tolerance is what rounding may cost. In practice the two
    # routes agree to floating point and the whole of the residual is round(),
    # so assert that too: it is the stronger claim, and if it ever weakens to
    # merely "within 0.5" something has changed that is worth looking at. R and
    # Python both round half to even, so this is a fair comparison.
    exact = sum(got[g][c] == round(exp[g][c]) for g in got for c in cols)
    ck(exact == len(diffs),
       f"{exact}/{len(diffs)} cells equal round(independent) exactly, so the "
       "residual above is rounding and nothing else")

    print("== the correction is real, and is applied where it should be ==")
    # Column sums are preserved by the final rescale. This is a property of the
    # definition, not something read off tximport's output.
    for c in cols[:1] + cols[-1:]:
        pipeline_sum = sum(got[g][c] for g in got)
        raw_sum = sum(raw[g][c] for g in raw)
        ck(abs(pipeline_sum - raw_sum) <= len(got) * 0.5,
           f"{c}: library total is preserved ({pipeline_sum:.0f} vs raw "
           f"{raw_sum:.0f}, within rounding)")

    # G001 holds its read depth and moves it from a short isoform to a long one.
    # Raw counts cannot see that; lengthScaledTPM must.
    ctl, trt = [s for s in cols if s in CONTROL], [s for s in cols if s in TREATED]
    raw_ratio = (sum(raw["G001"][s] for s in trt) / len(trt)) / \
                (sum(raw["G001"][s] for s in ctl) / len(ctl))
    cor_ratio = (sum(got["G001"][s] for s in trt) / len(trt)) / \
                (sum(got["G001"][s] for s in ctl) / len(ctl))
    ck(0.9 < raw_ratio < 1.1,
       f"G001 raw reads are flat across groups (ratio {raw_ratio:.2f})")
    ck(cor_ratio < 0.5,
       f"after correction G001 is much lower in treated (ratio {cor_ratio:.2f}): "
       "the same reads from a longer isoform are fewer molecules")

    # A single-transcript gene whose effective length is the same in every
    # sample is only touched by the per-column rescale, so all of them move by
    # one common factor. A per-gene length error would break this.
    single = [g for g in got if g not in ("G001",)]
    ratios = [got[g][cols[0]] / raw[g][cols[0]] for g in single]
    ck(max(ratios) - min(ratios) < 0.01,
       f"the {len(single)} single-transcript genes all move by one common "
       f"factor ({min(ratios):.4f} to {max(ratios):.4f})")

    print(f"\n{len(fails)} failed")
    return 1 if fails else 0


if __name__ == "__main__":
    if sys.argv[1] == "build":
        build(Path(sys.argv[2]), Path(sys.argv[3]),
              sys.argv[4] if len(sys.argv) > 4 else "tagseq")
    elif sys.argv[1] == "verify":
        sys.exit(verify(Path(sys.argv[2])))
    elif sys.argv[1] == "qualify":
        sys.exit(qualify(Path(sys.argv[2])))
    else:
        sys.exit(f"usage: {sys.argv[0]} build <dest> <repo> [seq_type] "
                 f"| verify <dest> | qualify <dest>")

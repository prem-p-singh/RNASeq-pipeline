#!/usr/bin/env python3
"""A tiny synthetic project, and the expected answers for Stages 2 to 5.

    tiny_project.py build  <dest> <repo>
    tiny_project.py verify <dest>

`build` writes a reference and Stage 1 outputs for 8 samples. Stages 0 and 1
need salmon and fastp, so their products are written here directly; everything
downstream then runs for real against them.

`verify` recomputes the expected gene counts from the fixture itself, by the
route tximport is not taking (parse the GTF, sum NumReads per gene), and checks
what the stages produced against it. An assertion that just re-ran the pipeline's
own arithmetic would confirm nothing.

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


def genes():
    """(gene_id, [transcript ids], treated-group multiplier)."""
    out = [
        # Two transcripts, so gene counts must be their sum. This is the one
        # thing tx2gene aggregation can get wrong without looking wrong.
        ("G001", ["T001a", "T001b"], 1.0),
        # Its quant.sf name carries a "|" suffix. If ignoreAfterBar does not
        # apply, this gene goes missing rather than erroring.
        ("G002", ["T002"], 1.0),
    ]
    out += [(f"UP{i:03d}", [f"UP{i:03d}_t1"], EFFECT) for i in range(N_DE)]
    out += [(f"DN{i:03d}", [f"DN{i:03d}_t1"], 1.0 / EFFECT) for i in range(N_DE)]
    out += [(f"G{i + 100:03d}", [f"G{i + 100:03d}_t1"], 1.0) for i in range(N_FLAT)]
    return out


def quant_name(tx):
    return BARRED_TX if tx == "T002" else tx


# --------------------------------------------------------------------------
def build(dest: Path, repo: Path):
    rng = random.Random(SEED)
    gs = genes()
    (dest / "config").mkdir(parents=True, exist_ok=True)
    ref = dest / "reference"
    ref.mkdir(exist_ok=True)

    gtf, fa = [], []
    for gid, txs, _ in gs:
        for j, tx in enumerate(txs):
            start = 1 + j * 2000
            gtf.append(f"chr1\tfixture\texon\t{start}\t{start + 999}\t.\t+\t."
                       f'\tgene_id "{gid}"; transcript_id "{tx}";')
            fa.append(f">{quant_name(tx)}\n" + "ACGT" * 250)
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
        rows = ["Name\tLength\tEffectiveLength\tTPM\tNumReads"]
        for gid, txs, mult in gs:
            eff = mult if s in TREATED else 1.0
            for tx in txs:
                # Split the two-transcript gene unevenly, so its sum is not the
                # same number as either part.
                share = 0.75 if tx.endswith("a") else 0.25 if tx.endswith("b") else 1.0
                n = max(1, int(round(BASE * eff * share
                                     * (1.0 + rng.uniform(-NOISE, NOISE)))))
                rows.append(f"{quant_name(tx)}\t1000\t850.0\t0.0\t{n}")
        (d / "quant.sf").write_text("\n".join(rows) + "\n")

        rate, reads = {LOW_MAPPING: (0.45, 5_000_000),
                       LOW_DEPTH: (0.82, 2_000_000)}.get(s, (0.88, 5_000_000))
        (d / "metrics.json").write_text(json.dumps({
            "sample": s, "seq_type": "tagseq",
            "num_reads_processed": reads, "mapping_rate": rate,
            "min_map_rate_threshold": 0.60, "flagged_low_mapping": rate < 0.60,
            "detected_libtype": "SF", "expected_libtype": "",
            "libtype_mismatch": False, "salmon_bias_flags": "",
            "salmon_version": "fixture", "fastp_version": "fixture",
            "fastp_report": "",
        }, indent=2) + "\n")

    sheet = ["sample_id\tfastq_url\ttreatment"] + [
        f"{s}\t/dev/null/{s}.fq.gz\t" + ("treated" if s in TREATED else "control")
        for s in SAMPLES]
    (dest / "config" / "samples.tsv").write_text("\n".join(sheet) + "\n")

    # seq_type tagseq makes countsFromAbundance "no", which is the branch whose
    # answer is exactly predictable: gene counts are summed NumReads.
    (dest / "config" / "config.yaml").write_text("""\
project: {name: fixture, description: tiny synthetic cohort, output_dir: "results/"}
organism:
  common_name: fixture
  scientific_name: Fixtura synthetica
  tax_id: 1
  kegg_code: null
  orgdb_package: org.Fsynthetica.eg.db
reference: {accession: FIXTURE_1.0, assembly_name: FIX1, cache_dir: null}
samples: {sheet: config/samples.tsv, seq_type: tagseq}
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
""")
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


if __name__ == "__main__":
    if sys.argv[1] == "build":
        build(Path(sys.argv[2]), Path(sys.argv[3]))
    elif sys.argv[1] == "verify":
        sys.exit(verify(Path(sys.argv[2])))
    else:
        sys.exit(f"usage: {sys.argv[0]} build <dest> <repo> | verify <dest>")

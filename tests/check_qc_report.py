#!/usr/bin/env python3
"""Self-check for workflow/scripts/qc_report.py.

The things that must not break:
  - the sample list comes from the sheet, so a stale directory from an earlier
    run on a different cohort never reaches the report
  - a sample listed in the sheet but missing from disk appears, as NA
  - an absent metric is NA in the table and draws no bar, and is never 0
  - a metric that genuinely IS 0 is still reported and drawn as 0
  - build_gene_matrix.clean_desc strips the species it is given, not a
    hardcoded one (that script shares read_sample_ids with qc_report)

Imports the module directly, which works because matplotlib is imported inside
make_charts. Run:  python3 tests/check_qc_report.py
"""
import importlib.util, json, os, tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
spec = importlib.util.spec_from_file_location(
    "qc_report", ROOT / "workflow" / "scripts" / "qc_report.py")
qc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(qc)

TMP = Path(tempfile.mkdtemp())


def write_sample(quant, sid, *, fastp=True, salmon=True, reads=1_000_000,
                 q30=0.95, mapped=900_000, decoy=50_000, processed=1_000_000):
    d = quant / sid
    (d / "aux_info").mkdir(parents=True, exist_ok=True)
    (d / "quant.sf").write_text("Name\tLength\tEffectiveLength\tTPM\tNumReads\n")
    if fastp:
        (d / "fastp.json").write_text(json.dumps({"summary": {
            "before_filtering": {"total_reads": reads, "q30_rate": q30},
            "after_filtering": {"total_reads": reads, "q30_rate": q30}}}))
    if salmon:
        (d / "aux_info" / "meta_info.json").write_text(json.dumps({
            "num_processed": processed, "num_mapped": mapped,
            "num_decoy_fragments": decoy, "percent_mapped": 90.0,
            "frag_length_mean": 210.0}))


def sheet(path, ids):
    path.write_text("sample_id\tgroup\n" + "".join(f"{i}\ta\n" for i in ids))
    return path


# --- read_sample_ids --------------------------------------------------
sh = sheet(TMP / "samples.tsv", ["S1", "S2", "S3"])
assert qc.read_sample_ids(sh) == ["S1", "S2", "S3"]

# comments and blank lines must not become sample IDs
sh2 = TMP / "commented.tsv"
sh2.write_text("# a note\n\nsample_id\tgroup\nS1\ta\n\n# trailing\nS2\tb\n")
assert qc.read_sample_ids(sh2) == ["S1", "S2"], qc.read_sample_ids(sh2)

bad = TMP / "nohdr.tsv"
bad.write_text("name\tgroup\nS1\ta\n")
try:
    qc.read_sample_ids(bad)
    raise AssertionError("a sheet without sample_id must be rejected")
except SystemExit as e:
    assert "sample_id" in str(e)

# --- stale directories are ignored ------------------------------------
quant = TMP / "quant"
for sid in ["S1", "S2"]:
    write_sample(quant, sid)
write_sample(quant, "OLD_COHORT_SAMPLE")        # left over from a previous run

recs = qc.collect(quant, ["S1", "S2"])
assert [r["sample"] for r in recs] == ["S1", "S2"], [r["sample"] for r in recs]
assert all(r["has_quant"] for r in recs)

# --- a sheet sample with no directory is present, as NA ---------------
recs = qc.collect(quant, ["S1", "S2", "NEVER_RAN"])
missing = [r for r in recs if r["sample"] == "NEVER_RAN"][0]
assert missing["has_quant"] is False
assert missing["reads_before"] is None and missing["num_mapped"] is None
assert qc.cell(missing["reads_before"]) == "NA"

# --- absent vs genuinely-zero metrics ---------------------------------
write_sample(quant, "NO_FASTP", fastp=False)          # salmon only
write_sample(quant, "ZERO_MAPPED", mapped=0, decoy=0) # real zeros
recs = qc.collect(quant, ["S1", "NO_FASTP", "ZERO_MAPPED"])
by = {r["sample"]: r for r in recs}

assert by["NO_FASTP"]["reads_before"] is None, "absent fastp must not become 0"
assert by["NO_FASTP"]["num_mapped"] == 900_000, "salmon metrics must still load"
assert by["ZERO_MAPPED"]["num_mapped"] == 0, "a real zero must survive as 0"

# present() skips the absent one and keeps the real zero
xs, ys = qc.present(recs, "reads_before", 1e6)
assert 1 not in xs, "a bar was placed for a sample with no fastp metrics"
assert xs == [0, 2], xs
xs, ys = qc.present(recs, "num_mapped", 1e6)
assert xs == [0, 1, 2] and ys[2] == 0.0, (xs, ys)

# --- TSV writes NA, not 0 ---------------------------------------------
out = TMP / "out"
out.mkdir()
qc.write_tsv(recs, out)
lines = (out / "alignment_summary.tsv").read_text().splitlines()
hdr = lines[0].split("\t")
row = dict(zip(hdr, lines[2].split("\t")))
assert row["sample"] == "NO_FASTP"
assert row["reads_before"] == "NA", row
assert row["q30_before"] == "NA", row
row0 = dict(zip(hdr, lines[3].split("\t")))
assert row0["sample"] == "ZERO_MAPPED" and row0["num_mapped"] == "0", row0

# --- tick labels flag the gaps ----------------------------------------
labels = qc.tick_labels(recs)
assert labels[0] == "S1", labels
assert "no fastp" in labels[1], labels
assert "no salmon" not in labels[1], labels

# --- identity comes from the run, with honest fallbacks ---------------
name, species, ref, title = qc.identity("MyRun", "Vitis vinifera", "GCF_1.2", "ASM9")
assert (name, species) == ("MyRun", "Vitis vinifera")
assert ref == "GCF_1.2 (ASM9)", ref
assert title == "MyRun — Vitis vinifera: cleaning & alignment QC", title

name, species, ref, title = qc.identity("", "", "", "")
assert name == "RNA-Seq run" and species == "" and ref == "unspecified"
assert title == "RNA-Seq run: cleaning & alignment QC", title

# the old hardcoded identity must not be reachable from any code path
src = (ROOT / "workflow" / "scripts" / "qc_report.py").read_text()
for baked in ("NCGM-4408", "Ascochyta", "GCF_004011695"):
    assert baked not in src, f"{baked} is still hardcoded in qc_report.py"

# --- build_gene_matrix.py, which shares read_sample_ids with this module ---
spec2 = importlib.util.spec_from_file_location(
    "build_gene_matrix", ROOT / "workflow" / "scripts" / "build_gene_matrix.py")
bgm = importlib.util.module_from_spec(spec2)
spec2.loader.exec_module(bgm)

# clean_desc used to strip one hardcoded binomial, so every other genome kept
# its species name in the description column.
assert bgm.clean_desc("Vitis vinifera chalcone synthase (CHS1), mRNA",
                      "Vitis vinifera") == "chalcone synthase"
assert bgm.clean_desc("PREDICTED: Vitis vinifera chalcone synthase, "
                      "transcript variant X2, mRNA",
                      "Vitis vinifera") == "chalcone synthase"
assert bgm.clean_desc("Ascochyta rabiei hypothetical protein, partial mRNA",
                      "Ascochyta rabiei") == "hypothetical protein"
# no species given: nothing is guessed away
assert bgm.clean_desc("Vitis vinifera chalcone synthase, mRNA",
                      "") == "Vitis vinifera chalcone synthase"
assert "Ascochyta" not in (ROOT / "workflow" / "scripts" / "build_gene_matrix.py").read_text()

# the cross-script import build_gene_matrix does at runtime must resolve
import sys
sys.path.insert(0, str(ROOT / "workflow" / "scripts"))
from qc_report import read_sample_ids as _rsi
assert _rsi(sh) == ["S1", "S2", "S3"]

print("check_qc_report.py: all assertions passed")

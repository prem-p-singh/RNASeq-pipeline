#!/usr/bin/env python3
"""Self-check for reference provenance (P0 / R17c, rules RF01, RF02, D11).

The defect this guards against: the QC report stated "decoy-aware salmon"
unconditionally while retrieve.smk built a transcriptome-only index, so every
report made a reference-construction claim the run did not support.

What must not break:
  - a transcriptome-only index is recorded as such, never as decoy-aware
  - a decoy-aware index is recorded as decoy-aware only when decoys were given
  - file identity (size and sha256) is recorded for the bundle files
  - the report's reference line is derived from the lock and from nothing else
  - with no lock, the report says construction is unrecorded and makes no
    decoy claim
  - the literal string "decoy-aware" cannot be produced for a decoy-free bundle
  - retrieve.smk and qc_report.py contain no unconditional decoy assertion

Run:  python3 tests/check_reference_lock.py
"""
import hashlib
import importlib.util
import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
LOCK_PY = ROOT / "workflow" / "scripts" / "reference_lock.py"
TMP = Path(tempfile.mkdtemp())

spec = importlib.util.spec_from_file_location("reference_lock", LOCK_PY)
rl = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rl)

spec2 = importlib.util.spec_from_file_location(
    "qc_report", ROOT / "workflow" / "scripts" / "qc_report.py")
qc = importlib.util.module_from_spec(spec2)
spec2.loader.exec_module(qc)


def make_bundle(name, *, decoys=False):
    d = TMP / name
    (d / "salmon_idx").mkdir(parents=True, exist_ok=True)
    fa = d / "transcriptome.fa"
    fa.write_text(">t1\nACGTACGTACGT\n>t2\nTTTTGGGGCCCC\n")
    gtf = d / "annotation.gtf"
    gtf.write_text('chr1\ttest\texon\t1\t12\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n')
    (d / "salmon_idx" / "info.json").write_text('{"index_version": 5, "k": 31}')
    decoy_file = None
    if decoys:
        decoy_file = d / "decoys.txt"
        decoy_file.write_text("chr1\n")
    return d, fa, gtf, decoy_file


def run_lock(d, fa, gtf, decoy_file):
    out = d / "reference.lock.json"
    cmd = [sys.executable, str(LOCK_PY), "--out", str(out),
           "--transcriptome", str(fa), "--gtf", str(gtf),
           "--index", str(d / "salmon_idx"),
           "--organism", "Vitis vinifera", "--tax-id", "29760",
           "--accession", "GCF_1.1", "--assembly-name", "ASM1",
           "--transcriptome-url", "https://example/tx.fa.gz",
           "--gtf-url", "https://example/ann.gtf.gz",
           "--kmer", "31"]
    if decoy_file:
        cmd += ["--decoys", str(decoy_file)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return json.loads(out.read_text()), out


# --- 1. transcriptome-only bundle is never called decoy-aware ----------
d, fa, gtf, _ = make_bundle("plain")
lock, lock_path = run_lock(d, fa, gtf, None)

decoys = lock["index"]["decoys"]
assert decoys["used"] is False, decoys
assert decoys["source"] is None, decoys
assert decoys["label"] == "transcriptome-only (no decoys)", decoys
assert "decoy-aware" not in json.dumps(lock).replace('"decoy-aware"', ""), (
    "a decoy-free bundle must not contain the phrase decoy-aware")
assert lock["index"]["construction"] == "salmon index -t <transcriptome>", lock["index"]

# file identity is real, not asserted
tx = lock["files"]["transcriptome"]
assert tx["present"] and tx["bytes"] == fa.stat().st_size, tx
assert tx["sha256"] == hashlib.sha256(fa.read_bytes()).hexdigest(), tx
assert lock["files"]["annotation"]["present"], lock["files"]

# identifiers carried through
assert lock["assembly"]["accession"] == "GCF_1.1"
assert lock["organism"]["tax_id"] == 29760
assert lock["sources"]["gtf_url"] == "https://example/ann.gtf.gz"

# --- 2. the report's reference line comes from the lock ---------------
_, desc = qc.reference_description(str(lock_path))
assert "GCF_1.1 (ASM1)" in desc, desc
assert "transcriptome-only (no decoys)" in desc, desc
assert "decoy-aware" not in desc, f"report would claim decoy-aware: {desc}"
assert "mapping-based" in desc, desc

# --- 3. no lock means unrecorded, and no decoy claim ------------------
_, desc_missing = qc.reference_description(str(TMP / "does_not_exist.json"))
assert "unrecorded" in desc_missing, desc_missing
assert "decoy" not in desc_missing.replace("decoy status unrecorded", ""), desc_missing

bad = TMP / "corrupt.json"
bad.write_text("{not json")
_, desc_bad = qc.reference_description(str(bad))
assert "unrecorded" in desc_bad, desc_bad

# --- 4. a genuinely decoy-aware bundle IS labelled so -----------------
d2, fa2, gtf2, decoy_file = make_bundle("decoyed", decoys=True)
lock2, lock2_path = run_lock(d2, fa2, gtf2, decoy_file)
assert lock2["index"]["decoys"]["used"] is True, lock2["index"]["decoys"]
assert lock2["index"]["decoys"]["label"] == "decoy-aware", lock2["index"]["decoys"]
assert lock2["index"]["decoys"]["source"] == str(decoy_file)
_, desc2 = qc.reference_description(str(lock2_path))
assert "decoy-aware" in desc2, desc2

# --- 5. the sources of the original defect stay clean -----------------
qc_src = (ROOT / "workflow" / "scripts" / "qc_report.py").read_text()
# strip comments and docstrings so the explanation of the fix is not a hit
code_only = "\n".join(
    l for l in qc_src.splitlines() if not l.strip().startswith("#"))
code_only = re.sub(r'""".*?"""', "", code_only, flags=re.S)
assert "decoy-aware" not in code_only, (
    "qc_report.py asserts decoy-aware outside a comment; it must come from the lock")

smk = (ROOT / "workflow" / "rules" / "retrieve.smk").read_text()
index_cmd = [l for l in smk.splitlines() if "salmon index" in l and "-t" in l]
assert index_cmd, "salmon index command not found in retrieve.smk"
for line in index_cmd:
    if " -d " in line or "--decoys" in line:
        continue
    # Transcriptome-only build: nothing may claim otherwise.
    assert "decoy" not in line.lower(), line

print("check_reference_lock.py: all assertions passed")

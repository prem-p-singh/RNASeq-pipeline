#!/usr/bin/env python3
"""Canonical lanes reach preparation intact; invalid mates cannot publish success."""
from copy import deepcopy
import gzip
import json
from pathlib import Path
import sys
import tempfile

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
import metadata as md
from prepare_reads import prepare

with tempfile.TemporaryDirectory() as directory:
    root = Path(directory)
    tables = {"samples": [{"sample_id": "001", "condition": "control"}],
              "libraries": [{"sample_id": "001", "library_id": "lib1", "layout": "paired",
                             "assay_family": "bulk", "umi": "no", "strandedness": "reverse"}],
              "reads": []}
    for run in ("run_b", "run_a"):
        for mate in (1, 2):
            path = root / f"{run} read {mate}.fq.gz"
            with gzip.open(path, "wb") as handle:
                for i in range(2):
                    handle.write(f"@{run}_{i}/{mate}\nACGT\n+\nIIII\n".encode())
            tables["reads"].append({"read_unit_id": f"{run}_{mate}", "library_id": "lib1",
                                    "role": f"R{mate}", "run": run, "lane": "1", "uri": str(path)})
    md.write_tables(tables, root / "metadata")
    cfg = {"samples": {"metadata_dir": "metadata", "sheet": "metadata/samples.tsv", "seq_type": "rnaseq_paired"}}
    manifest = md.execution_inputs(root, cfg)["001"]
    assert [u["run"] for u in manifest["units"]] == ["run_a", "run_b"]
    ledger = prepare(manifest, root / "output")
    assert ledger["fragments"] == 4 and len(ledger["units"]) == 2
    with gzip.open(root / "output/prepared_R1.fastq.gz", "rt") as handle:
        merged = handle.read()
    assert merged.index("run_a") < merged.index("run_b")
    assert all(Path(r["uri"]).exists() for r in tables["reads"])
    assert len(ledger["units"][0]["inputs"][0]["sha256"]) == 64

    def rejected(candidate, message):
        try:
            prepare(candidate, root / "bad")
        except ValueError as exc:
            assert message in str(exc), str(exc)
        else:
            raise AssertionError("Invalid reads accepted")
        assert not (root / "bad/read_preparation.json").exists()
        assert not (root / "bad/prepared_R1.fastq.gz").exists()

    bad = deepcopy(manifest)
    bad["units"][0]["r1"]["checksum"] = "0" * 64
    rejected(bad, "Checksum mismatch")
    collision = root / "collision"
    collision.mkdir()
    source_ledger = collision / "read_preparation.json"
    source_ledger.write_bytes(b"user-owned data")
    bad = deepcopy(manifest)
    bad["units"][0]["r1"]["uri"] = str(source_ledger)
    try:
        prepare(bad, collision)
    except ValueError as exc:
        assert "collision" in str(exc)
    else:
        raise AssertionError("Source/output collision accepted")
    assert source_ledger.read_bytes() == b"user-owned data"
    bad = deepcopy(manifest)
    bad["units"][0]["r2"] = None
    rejected(bad, "mate structure")
    bad = deepcopy(manifest)
    bad["units"][1]["r1"]["uri"] = bad["units"][0]["r1"]["uri"]
    rejected(bad, "Repeated input")
    with gzip.open(manifest["units"][1]["r2"]["uri"], "wb") as handle:
        handle.write(b"@wrong/2\nACGT\n+\nIIII\n")
    rejected(manifest, "Mismatched mate")
    tables["libraries"][0]["umi"] = "yes"
    md.write_tables(tables, root / "metadata")
    try:
        md.execution_inputs(root, cfg)
    except ValueError as exc:
        assert "umi=no" in str(exc)
    else:
        raise AssertionError("UMI library accepted")

print("check_read_units.py: canonical run/lane preparation and failure checks passed")

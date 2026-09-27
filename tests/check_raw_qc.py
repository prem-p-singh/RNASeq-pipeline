#!/usr/bin/env python3
"""Raw QC needs no reference/model; failed libraries cannot erase earlier reports."""
import gzip
import json
import shutil
from pathlib import Path
import subprocess
import sys
import tempfile
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from metadata import write_tables
from raw_qc import run

if "--real" in sys.argv and not shutil.which("fastp"):
    print("needs locked fastp 0.23.4")
    sys.exit(77)

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / "config").mkdir()
    # Deliberately omit references and the DE model: QC must not need them.
    (root / "config/config.yaml").write_text(
        'samples:\n  metadata_dir: metadata\n  sheet: metadata/samples.tsv\n  seq_type: rnaseq_paired\n')
    tables = {"samples": [], "libraries": [], "reads": []}
    for sid in ("001", "NA"):
        tables["samples"].append({"sample_id": sid})
        tables["libraries"].append({"sample_id": sid, "library_id": sid, "layout": "paired",
            "assay_family": "bulk", "umi": "no", "strandedness": "unknown"})
        for mate in (1, 2):
            src = root / f"{sid}_{mate}.fq.gz"
            with gzip.open(src, "wt") as f:
                f.write(f'@read/{mate}\nACGT\n+\nIIII\n')
            tables["reads"].append({"read_unit_id": f"{sid}_{mate}", "library_id": sid,
                "role": f"R{mate}", "uri": str(src)})
    write_tables(tables, root / "metadata")
    originals = {r["uri"]: Path(r["uri"]).read_bytes() for r in tables["reads"]}

    def fastp(cmd, **kwargs):
        if cmd == ["fastp", "--version"]:
            return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="fastp 0.23.4")
        assert "--in2" in cmd and "--disable_trim_poly_g" in cmd
        assert "--disable_adapter_trimming" in cmd and "--disable_quality_filtering" in cmd
        assert "--disable_length_filtering" in cmd and "--overrepresentation_analysis" in cmd
        out = Path(cmd[cmd.index("--json") + 1])
        if out.parent.name == "NA":
            raise subprocess.CalledProcessError(1, cmd)
        out.write_text(json.dumps({"summary": {phase: {"total_reads": 2}
            for phase in ("before_filtering", "after_filtering")}}))
        Path(cmd[cmd.index("--html") + 1]).write_text("<html>QC</html>")
        return subprocess.CompletedProcess(cmd, 0)

    with patch("raw_qc.subprocess.run", side_effect=fastp):
        try:
            run(root, root / "qc")
        except subprocess.CalledProcessError:
            pass
        else:
            raise AssertionError("Tool failure was swallowed")
        data = json.loads((root / "qc/manifest.json").read_text())
        assert data["samples"]["001"]["status"] == "succeeded"
        assert data["samples"]["NA"]["status"] == "failed"
        assert (root / "qc/001/fastp.html").exists()
        assert (root / "qc/index.html").exists()
        assert not list((root / "qc").glob("*/.reads-*"))
        try:
            run(root, root / "qc")
        except FileExistsError:
            pass
        else:
            raise AssertionError("Existing reports overwritten")
    assert all(Path(p).read_bytes() == data for p, data in originals.items())
    if "--real" in sys.argv:
        for layout in ("paired", "single"):
            if layout == "single":
                for lib in tables["libraries"]:
                    lib["layout"] = "single"
                tables["reads"] = [r for r in tables["reads"] if r["role"] == "R1"]
                write_tables(tables, root / "metadata")
                cfg = root / "config/config.yaml"
                cfg.write_text(cfg.read_text().replace("rnaseq_paired", "rnaseq_single"))
            result = run(root, root / f"real-{layout}")
            assert all(s["status"] == "succeeded" for s in result["samples"].values())
        assert all(Path(p).read_bytes() == data for p, data in originals.items())
        print("Real fastp: SE/PE count conservation and reports passed without references or model")

print("check_raw_qc.py: independent QC failure/provenance/ownership contract passed (mocked fastp)")

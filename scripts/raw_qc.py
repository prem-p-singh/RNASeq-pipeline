#!/usr/bin/env python3
"""Independent, observation-only fastp QC for canonical non-UMI bulk projects."""
import argparse
import html
import json
from pathlib import Path
import subprocess
import tempfile

import yaml

from metadata import execution_inputs
from prepare_reads import prepare


def run(project, output, threads=2):
    project, output = Path(project).resolve(), Path(output).resolve()
    if not 1 <= threads <= 16:
        raise ValueError("threads must be between 1 and 16")
    cfg = yaml.safe_load((project / "config/config.yaml").read_text())
    inputs = execution_inputs(project, cfg)
    if inputs is None:
        raise ValueError("Independent raw QC requires canonical Samples/Libraries/Reads tables")
    version = subprocess.run(["fastp", "--version"], check=True, capture_output=True, text=True)
    version = (version.stdout + version.stderr).strip()
    if version.split()[-1:] != ["0.23.4"]:
        raise ValueError(f"Use the locked fastp 0.23.4 runtime; found {version!r}")
    # A new output directory prevents stale success and protects existing reports.
    output.mkdir(parents=True, exist_ok=False)
    report = {"schema_version": 1, "mode": "raw_qc_only", "fastp_version": version,
              "processing": "observation only; no trimming, filtering or deduplication",
              "unavailable": ["contamination screening", "adapter detection", "alignment metrics"],
              "samples": {sid: {"status": "pending"} for sid in inputs}}

    def publish():
        (output / "manifest.json").write_text(json.dumps(report, indent=2) + "\n")
        rows = []
        for sid, item in report["samples"].items():
            link = f'<a href="{sid}/fastp.html">report</a>' if item["status"] == "succeeded" else ""
            rows.append(f'<tr><td>{html.escape(sid)}</td><td>{item["status"]}</td><td>{link}</td></tr>')
        (output / "index.html").write_text('<!doctype html><meta charset="utf-8"><title>Raw-read QC</title>'
            '<h1>Raw-read QC</h1><p>Observation only. No reference, quantification or inference required. '
            'Contamination, adapter detection and alignment metrics are unavailable.</p>'
            '<table><tr><th>Sample</th><th>Status</th><th>QC</th></tr>' + ''.join(rows) + '</table>\n')

    publish()
    for sid, manifest in inputs.items():
        item = report["samples"][sid]
        dest = output / sid
        dest.mkdir()
        try:
            with tempfile.TemporaryDirectory(prefix=".reads-", dir=dest) as tmp:
                ledger = prepare(manifest, tmp)
                (dest / "read_preparation.json").write_text(json.dumps(ledger, indent=2) + "\n")
                command = ["fastp", "--in1", str(Path(tmp) / "prepared_R1.fastq.gz"),
                           "--json", str(dest / "fastp.json"), "--html", str(dest / "fastp.html"),
                           "--thread", str(threads), "--disable_adapter_trimming",
                           "--disable_quality_filtering", "--disable_length_filtering",
                           "--disable_trim_poly_g", "--overrepresentation_analysis"]
                if manifest["layout"] == "paired":
                    command.extend(["--in2", str(Path(tmp) / "prepared_R2.fastq.gz")])
                item["command"] = command
                with (dest / "fastp.log").open("w") as log:
                    subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
                data = json.loads((dest / "fastp.json").read_text())
                expected = ledger["fragments"] * (2 if manifest["layout"] == "paired" else 1)
                for phase in ("before_filtering", "after_filtering"):
                    if data["summary"][phase]["total_reads"] != expected:
                        raise ValueError("fastp read counts disagree with validated input fragments")
                if not (dest / "fastp.html").stat().st_size:
                    raise ValueError("fastp HTML report is empty")
                item.update(status="succeeded", fragments=ledger["fragments"], reads=expected)
        except Exception as exc:
            item.update(status="failed", error=str(exc))
            raise
        finally:
            publish()
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-d", "--project-dir", required=True)
    parser.add_argument("--out", required=True, help="New output directory; existing directories are refused")
    parser.add_argument("--threads", type=int, default=2)
    args = parser.parse_args()
    run(args.project_dir, args.out, args.threads)

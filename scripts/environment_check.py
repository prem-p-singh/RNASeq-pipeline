#!/usr/bin/env python3
"""Verify an activated Linux environment against the release's explicit lock.

Uses installed Conda metadata, executable versions and actual R namespace loads.
This verifies package identity, not a byte-for-byte audit of installed files.
"""
import argparse
import hashlib
import importlib
import json
import os
from pathlib import Path
import platform
import shutil
import socket
import subprocess
import sys
from urllib.parse import urldefrag

ROOT = Path(__file__).resolve().parent.parent


def verify(prefix, lock):
    prefix, lock = Path(prefix).resolve(), Path(lock)
    issues, versions = [], {}
    if platform.system() != "Linux" or platform.machine() != "x86_64":
        issues.append("The release runtime supports linux-64; use the Linux deployment.")
    if not lock.is_file():
        return {"status": "failed", "issues": [f"Missing explicit lock: {lock}"]}
    expected = {}
    for line in lock.read_text().splitlines():
        if line.startswith("https://"):
            url, digest = urldefrag(line)
            if len(digest) != 64:
                issues.append(f"Lock entry has no SHA256: {url}")
            expected[Path(url).name] = digest
    installed = {}
    for record in (prefix / "conda-meta").glob("*.json"):
        data = json.loads(record.read_text())
        installed[data.get("fn", "")] = data
    if not expected:
        issues.append("Lock has no package records")
    for filename, digest in expected.items():
        rec = installed.get(filename)
        if rec is None or rec.get("sha256") != digest:
            issues.append(f"Missing or mismatched package: {filename}")
    for filename in installed.keys() - expected.keys():
        issues.append(f"Unrecorded Conda package: {filename}")
    for name in ("yaml", "pandas", "openpyxl", "matplotlib"):
        try:
            mod = importlib.import_module(name)
            if not Path(mod.__file__).resolve().is_relative_to(prefix):
                issues.append(f"Python module outside environment: {name}")
        except Exception as exc:
            issues.append(f"Cannot import {name}: {exc}")
    for tool in ("python3", "Rscript", "snakemake", "salmon", "fastp", "multiqc"):
        path = shutil.which(tool)
        if not path or not Path(path).resolve().is_relative_to(prefix):
            issues.append(f"Executable absent or outside environment: {tool}")
            continue
        run = subprocess.run([path, "--version"], capture_output=True, text=True)
        versions[tool] = (run.stdout + run.stderr).strip()
        if run.returncode:
            issues.append(f"Version probe failed: {tool}")
    if shutil.which("Rscript"):
        probe = '''p <- c("tximport","txdbmaker","GenomicFeatures","AnnotationDbi",
          "AnnotationForge","edgeR","limma","variancePartition","clusterProfiler",
          "GOSemSim","WGCNA","emmeans","jsonlite","dplyr","readr","tibble")
          root <- normalizePath(commandArgs(TRUE)[1])
          for (n in p) {
            if (!requireNamespace(n, quietly=TRUE)) stop("Missing R package: ", n)
            if (!startsWith(normalizePath(find.package(n)), paste0(root, "/")))
              stop("R package outside environment: ", n)
          }
          cat(paste(p, vapply(p, function(n) as.character(packageVersion(n)), ""), collapse="\\n"))'''
        run = subprocess.run(["Rscript", "--vanilla", "-e", probe, str(prefix)],
                             capture_output=True, text=True)
        versions["R_packages"] = run.stdout.strip()
        if run.returncode:
            issues.append(f"R namespace verification failed: {run.stderr.strip()}")
    return {"schema_version": 1, "status": "failed" if issues else "verified",
            "lock_sha256": hashlib.sha256(lock.read_bytes()).hexdigest(),
            "prefix": str(prefix), "platform": platform.platform(),
            "hostname": socket.gethostname(), "slurm_job_id": os.getenv("SLURM_JOB_ID"),
            "package_count": len(installed), "versions": versions, "issues": issues}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", default=sys.prefix)
    parser.add_argument("--lock", type=Path, default=ROOT / "environments/linux-64.explicit.txt")
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    report = verify(args.prefix, args.lock)
    payload = json.dumps(report, indent=2) + "\n"
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(payload)
    print(payload)
    sys.exit(report["status"] != "verified")

#!/usr/bin/env python3
"""Explicit non-UMI bulk preprocessing and read-loss provenance."""
import argparse
import hashlib
import itertools
import json
import re
from pathlib import Path
import subprocess

from prepare_reads import records

DEFAULTS = dict(poly_g="auto", adapter_r1="", adapter_r2="", minimum_length=15,
                qualified_quality_phred=15, unqualified_percent_limit=40)


def policy_args(policy, paired):
    unknown = set(policy) - set(DEFAULTS)
    if unknown:
        raise ValueError(f"Unknown preprocessing settings: {sorted(unknown)}")
    p = DEFAULTS | policy
    if p["poly_g"] not in ("auto", "on", "off"):
        raise ValueError("preprocessing.poly_g must be auto, on or off")
    for name, low, high in (("minimum_length", 1, 10000),
                            ("qualified_quality_phred", 0, 93),
                            ("unqualified_percent_limit", 0, 100)):
        if type(p[name]) is not int or not low <= p[name] <= high:
            raise ValueError(f"preprocessing.{name} must be an integer in {low}..{high}")
    for name in ("adapter_r1", "adapter_r2"):
        if not isinstance(p[name], str) or set(p[name].upper()) - set("ACGT"):
            raise ValueError(f"preprocessing.{name} must be an A/C/G/T sequence or empty")
    if p["adapter_r2"] and (not paired or not p["adapter_r1"]):
        raise ValueError("adapter_r2 requires paired input and an explicit adapter_r1")
    args = ["--length_required", str(p["minimum_length"]),
            "--qualified_quality_phred", str(p["qualified_quality_phred"]),
            "--unqualified_percent_limit", str(p["unqualified_percent_limit"]),
            "--overrepresentation_analysis"]
    if p["poly_g"] != "auto":
        args.append("--trim_poly_g" if p["poly_g"] == "on" else "--disable_trim_poly_g")
    if p["adapter_r1"]:
        args.extend(["--adapter_sequence", p["adapter_r1"]])
        if paired:
            args.extend(["--adapter_sequence_r2", p["adapter_r2"] or p["adapter_r1"]])
    elif paired:
        args.append("--detect_adapter_for_pe")
    return p, args


def process(r1, r2, outdir, sample, threads, policy):
    if not re.fullmatch(r"[A-Za-z0-9_][A-Za-z0-9_.-]*", sample):
        raise ValueError("Invalid sample identifier")
    if type(threads) is not int or threads < 1:
        raise ValueError("threads must be a positive integer")
    paired = bool(r2)
    p, flags = policy_args(policy, paired)
    out = Path(outdir)
    sources = [Path(r1)] + ([Path(r2)] if paired else [])
    trims = [out / f"{sample}_R1.trim.fastq.gz", out / f"{sample}_R2.trim.fastq.gz"] if paired else [out / f"{sample}.trim.fastq.gz"]
    reserved = trims + [out / f"fastp.{ext}" for ext in ("json", "html", "log")] + [out / "preprocessing.json"]
    if {x.resolve() for x in sources} & {x.resolve() for x in reserved}:
        raise ValueError("Preprocessing output/source collision")
    out.mkdir(parents=True, exist_ok=True)
    (out / "preprocessing.json").unlink(missing_ok=True)
    inputs = []
    for path in sources:
        with path.open("rb") as f:
            digest = hashlib.file_digest(f, "sha256").hexdigest()
        inputs.append(dict(path=str(path), bytes=path.stat().st_size, sha256=digest))
    readers = [records(path) for path in sources]
    fragments = 0
    try:
        for unit in itertools.zip_longest(*readers):
            if any(x is None for x in unit) or (paired and unit[0][0] != unit[1][0]):
                raise ValueError("Mismatched input mate IDs/counts")
            fragments += 1
    finally:
        for reader in readers:
            reader.close()
    if not fragments:
        raise ValueError("Empty FASTQ input")
    version = subprocess.run(["fastp", "--version"], check=True, capture_output=True, text=True)
    version = (version.stdout + version.stderr).strip()
    cmd = ["fastp", "-i", str(sources[0]), "-o", str(trims[0]),
           "--json", str(out / "fastp.json"), "--html", str(out / "fastp.html"),
           "--thread", str(threads)] + flags
    if paired:
        cmd.extend(["-I", str(sources[1]), "-O", str(trims[1])])
    with (out / "fastp.log").open("w") as log:
        subprocess.run(cmd, check=True, stdout=log, stderr=subprocess.STDOUT)
    data = json.loads((out / "fastp.json").read_text())
    before = data["summary"]["before_filtering"]["total_reads"]
    after = data["summary"]["after_filtering"]["total_reads"]
    mates = 2 if paired else 1
    if before != fragments * mates or not 0 <= after <= before or after % mates:
        raise ValueError("fastp read counts disagree with validated fragment units")
    if not all(path.is_file() for path in trims) or not (out / "fastp.html").stat().st_size:
        raise ValueError("Missing preprocessing output")
    loss = (before - after) / before
    ledger = dict(schema_version=1, sample=sample, layout="paired" if paired else "single",
                  fastp_version=version, policy=p, command=cmd, input_files=inputs,
                  fragments_before=fragments, fragments_after=after // mates,
                  reads_before=before, reads_after=after, fraction_lost=loss,
                  filtering_result=data.get("filtering_result", {}),
                  warnings=["More than 10% of reads lost: review preprocessing and library quality"] if loss > 0.10 else [],
                  contamination_screening="not_performed")
    (out / "preprocessing.json").write_text(json.dumps(ledger, indent=2) + "\n")
    return ledger


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--r1", required=True)
    ap.add_argument("--r2", default="")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--policy-json", default="{}")
    a = ap.parse_args()
    process(a.r1, a.r2, a.outdir, a.sample, a.threads, json.loads(a.policy_json))

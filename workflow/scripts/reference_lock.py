#!/usr/bin/env python3
"""reference_lock.py — record what the reference bundle actually is.

Implements RF01 ("lock organism, assembly, annotation release, file hashes,
identifiers, tool/index parameters and construction steps") and the
reference.lock.json contract in master plan section 12.1.

The immediate reason this exists: the QC report used to state "decoy-aware
salmon" unconditionally while retrieve.smk built a transcriptome-only index.
RF02 requires that a transcriptome-only reference is recorded as such and never
labelled decoy-aware, and D11 requires recording the actual mapping mode,
decoys and construction rather than a prose description. Nothing may describe
the reference from anything other than this record.

Usage:
    reference_lock.py --out reference.lock.json \
        --transcriptome reference/transcriptome.fa \
        --gtf reference/annotation.gtf \
        --index reference/salmon_idx \
        --organism "Vitis vinifera" --tax-id 29760 \
        --accession GCF_030704535.1 --assembly-name ASM3070453v1 \
        --transcriptome-url URL --gtf-url URL \
        --kmer 31
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
from datetime import datetime, timezone


def sha256(path: str, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def file_record(path: str) -> dict:
    """Identity of one bundle file. Absent is recorded, never guessed."""
    if not path or not os.path.isfile(path):
        return {"path": path or None, "present": False,
                "bytes": None, "sha256": None}
    return {"path": path, "present": True,
            "bytes": os.path.getsize(path), "sha256": sha256(path)}


def salmon_version() -> str | None:
    """The version on PATH now. NOT necessarily the one that built the index.

    WORKING_PLAN 3.1: querying the available Salmon "need not identify the
    version that built a reused index". Callers pass --salmon-version, captured
    in the shell that ran `salmon index`; this is only the fallback, and the
    record marks which of the two it holds.
    """
    try:
        out = subprocess.run(["salmon", "--version"], capture_output=True,
                             text=True, timeout=60)
    except (OSError, subprocess.SubprocessError):
        return None
    text = (out.stdout or out.stderr or "").strip().splitlines()
    return text[0].replace('"', "").strip() if text else None


def index_record(index_dir: str, kmer: int, decoy_file: str | None,
                 built_version: str | None = None) -> dict:
    """What was actually built.

    `decoys` is derived from whether a decoy file was supplied to
    `salmon index -d`, not from an assumption about the route. A
    transcriptome-only index reports decoys.used == False and a
    `construction` string that says so, which is what any report must print.
    """
    used = bool(decoy_file)
    info_json = os.path.join(index_dir, "info.json") if index_dir else ""
    info = None
    if info_json and os.path.isfile(info_json):
        try:
            with open(info_json) as fh:
                info = json.load(fh)
        except (OSError, ValueError):
            info = None

    # Build-time capture is authoritative; the ambient query is a labelled
    # fallback so a reader can tell the difference.
    version = built_version or salmon_version()
    provenance = "build_time" if built_version else "observed_at_record_time"

    return {
        "tool": "salmon",
        "tool_version": version,
        "tool_version_provenance": provenance,
        "path": index_dir or None,
        "present": bool(index_dir and os.path.isdir(index_dir)),
        "kmer": kmer,
        # Salmon's own index metadata, when it wrote any. Not synthesised.
        "salmon_info": info,
        "decoys": {
            "used": used,
            "source": decoy_file if used else None,
            # RF02: the exact wording a report is allowed to use.
            "label": "decoy-aware" if used else "transcriptome-only (no decoys)",
        },
        # D11: salmon distinguishes mapping-based from alignment-based use.
        # This bundle is built for mapping-based quantification against the
        # transcriptome; it is not an alignment-based route.
        "quantification_mode": "mapping-based",
        "construction": (
            "salmon index -t <transcriptome> -d <decoys>"
            if used else "salmon index -t <transcriptome>"
        ),
    }


def build_lock(args) -> dict:
    return {
        "schema_version": 1,
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "organism": {
            "scientific_name": args.organism or None,
            "tax_id": args.tax_id,
        },
        "assembly": {
            "accession": args.accession or None,
            "assembly_name": args.assembly_name or None,
        },
        "sources": {
            "transcriptome_fasta_url": args.transcriptome_url or None,
            "gtf_url": args.gtf_url or None,
            "genome_fasta_url": args.genome_url or None,
        },
        "files": {
            "transcriptome": file_record(args.transcriptome),
            "annotation": file_record(args.gtf),
        },
        "index": index_record(args.index, args.kmer, args.decoys,
                              args.salmon_version or None),
        # RF09: what a reader must match before reusing this entry.
        "cache": {"key_inputs": json.loads(args.cache_key_inputs)
                  if args.cache_key_inputs else None},
    }


def describe(lock: dict) -> str:
    """One line a report may print. Derived only from the lock."""
    asm = lock.get("assembly", {}) or {}
    idx = lock.get("index", {}) or {}
    acc = asm.get("accession") or "unspecified accession"
    name = asm.get("assembly_name")
    ref = f"{acc} ({name})" if name else acc
    decoys = (idx.get("decoys") or {}).get("label") or "construction unrecorded"
    tool = idx.get("tool") or "unknown quantifier"
    ver = (idx.get("tool_version") or "").strip()
    # salmon --version prints "salmon 1.10.3", so joining tool and version
    # naively gave "salmon salmon 1.10.3".
    if ver.lower().startswith(tool.lower()):
        tool_s = ver
    else:
        tool_s = f"{tool} {ver}" if ver else tool
    mode = idx.get("quantification_mode") or "unrecorded mode"
    return f"{ref}, {tool_s}, {mode}, {decoys}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--transcriptome", default="")
    ap.add_argument("--gtf", default="")
    ap.add_argument("--index", default="")
    ap.add_argument("--organism", default="")
    ap.add_argument("--tax-id", type=int, default=None)
    ap.add_argument("--accession", default="")
    ap.add_argument("--assembly-name", default="")
    ap.add_argument("--transcriptome-url", default="")
    ap.add_argument("--gtf-url", default="")
    ap.add_argument("--genome-url", default="")
    ap.add_argument("--kmer", type=int, default=31)
    ap.add_argument("--salmon-version", default="",
                    help="version captured in the shell that built the index; "
                         "authoritative, unlike the version on PATH later")
    ap.add_argument("--cache-key-inputs", default="",
                    help="JSON of the construction parameters this entry was "
                         "built from, for RF09 reuse checks")
    ap.add_argument("--decoys", default="",
                    help="decoy file passed to `salmon index -d`; empty means "
                         "no decoys were used, which is recorded as such")
    args = ap.parse_args()

    lock = build_lock(args)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(lock, fh, indent=2, sort_keys=False)
        fh.write("\n")
    print(f"reference_lock: {args.out}")
    print(f"  {describe(lock)}")


if __name__ == "__main__":
    main()

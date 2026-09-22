#!/usr/bin/env python3
"""build_salmon_index.py — build a Salmon index safely into a shared cache.

Implements the cache half of R03 and R17d. Master plan 5.3 requires a shared
cache to use "content/version identity, an ordinary filesystem build lock, a
temporary build directory, validation, and atomic publication. Readers never
consume an incomplete cache entry."

Previously `salmon index` wrote straight into the shared path, so an
interrupted build left a partial index at the name every project reads, and two
runs on the same accession raced each other.

Order here: take the build lock; if another run published a valid entry while
we waited, use it; otherwise build into a staging directory, capture the Salmon
version from the process that actually built it, write reference.lock.json
inside the staging directory, and only then rename it into place.

The Salmon version matters. WORKING_PLAN 3.1 notes that querying the available
Salmon "need not identify the version that built a reused index", so it is read
here, in the same invocation as the build, and recorded as build_time.
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "scripts"))
import reference_cache  # noqa: E402


def salmon_version() -> str:
    out = subprocess.run(["salmon", "--version"], capture_output=True, text=True)
    text = (out.stdout or out.stderr or "").strip().splitlines()
    return text[0].replace('"', "").strip() if text else "unknown"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--transcriptome", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--final-index", required=True)
    ap.add_argument("--staging", required=True)
    ap.add_argument("--lock-out", required=True)
    ap.add_argument("--helper", required=True, help="reference_lock.py")
    ap.add_argument("--kmer", type=int, default=31)
    ap.add_argument("--decoys", default="")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--key-inputs", default="")
    for opt in ("organism", "accession", "assembly-name", "transcriptome-url",
                "gtf-url", "genome-url"):
        ap.add_argument(f"--{opt}", default="")
    ap.add_argument("--tax-id", type=int, default=0)
    a = ap.parse_args()

    final_index = Path(a.final_index)
    staging = Path(a.staging)
    entry = final_index.parent               # the cache entry directory
    lock_out = Path(a.lock_out)

    with reference_cache.BuildLock(final_index):
        # Another run may have finished while we waited for the lock. Its
        # output is as valid as ours would be, so do not rebuild.
        if final_index.is_dir() and (final_index / "info.json").is_file() \
                and lock_out.is_file():
            print(f"build_salmon_index: {final_index} was published by another "
                  f"run while waiting for the lock; reusing it")
            return

        if staging.exists():
            shutil.rmtree(staging)
        staging.parent.mkdir(parents=True, exist_ok=True)

        cmd = ["salmon", "index", "-t", a.transcriptome, "-i", str(staging),
               "-k", str(a.kmer), "--threads", str(a.threads)]
        if a.decoys:
            cmd += ["-d", a.decoys]
        print("build_salmon_index: " + " ".join(cmd), flush=True)
        subprocess.run(cmd, check=True)

        if not (staging / "info.json").is_file():
            raise SystemExit(
                f"salmon index produced no info.json in {staging}; refusing to "
                f"publish an incomplete entry")

        # Written inside the staging directory, so it is published with the
        # index in the same rename and can never describe a different build.
        staged_lock = staging / "reference.lock.json"
        lock_cmd = [
            sys.executable, a.helper,
            "--out", str(staged_lock),
            "--transcriptome", a.transcriptome,
            "--gtf", a.gtf,
            "--index", str(staging),
            "--kmer", str(a.kmer),
            "--decoys", a.decoys,
            "--salmon-version", salmon_version(),
            "--organism", a.organism,
            "--tax-id", str(a.tax_id),
            "--accession", a.accession,
            "--assembly-name", a.assembly_name,
            "--transcriptome-url", a.transcriptome_url,
            "--gtf-url", a.gtf_url,
            "--genome-url", a.genome_url,
        ]
        if a.key_inputs:
            lock_cmd += ["--cache-key-inputs", a.key_inputs]
        subprocess.run(lock_cmd, check=True)

        # The lock records the staging path; rewrite it to the published one so
        # the entry does not describe a directory that no longer exists.
        rec = json.loads(staged_lock.read_text())
        rec["index"]["path"] = str(final_index)
        staged_lock.write_text(json.dumps(rec, indent=2) + "\n")

        entry.mkdir(parents=True, exist_ok=True)
        shutil.move(str(staged_lock), str(lock_out))
        result = reference_cache.publish(staging, final_index)
        print(f"build_salmon_index: {result} -> {final_index}")


if __name__ == "__main__":
    main()

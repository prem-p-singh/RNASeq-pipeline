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
import re
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "scripts"))
import reference_cache  # noqa: E402
from reference_lock import file_record


def salmon_version() -> str:
    out = subprocess.run(["salmon", "--version"], capture_output=True, text=True)
    text = (out.stdout or out.stderr or "").strip().splitlines()
    return text[0].replace('"', "").strip() if text else "unknown"


def check_identifiers(fasta, gtf, decoys=""):
    """Every quantifiable FASTA target must map unambiguously to one GTF gene."""
    names = []
    with open(fasta) as f:
        for line in f:
            if line.startswith(">"):
                names.append(line[1:].split()[0].split("|")[0])
    if len(names) != len(set(names)):
        raise ValueError("Duplicate transcript identifiers in reference FASTA")
    mapping = {}
    with open(gtf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9:
                raise ValueError("Malformed GTF row")
            tx = re.search(r'transcript_id "([^"]+)"', fields[8])
            gene = re.search(r'gene_id "([^"]+)"', fields[8])
            if tx and gene:
                previous = mapping.setdefault(tx[1], gene[1])
                if previous != gene[1]:
                    raise ValueError(f"Transcript {tx[1]} maps to multiple genes")
    excluded = set(Path(decoys).read_text().splitlines()) if decoys else set()
    targets = set(names) - excluded
    absent = targets - mapping.keys()
    if not targets or absent:
        raise ValueError(f"FASTA/GTF mismatch: {len(absent)} unmapped targets; examples: {sorted(absent)[:10]}")
    return {"quantifiable_targets": len(targets), "mapped_targets": len(targets),
            "annotation_transcripts": len(mapping), "decoys": len(excluded)}


def prepare_genome_decoys(transcriptome, genome, combined, names):
    """Append genomic sequences after transcripts, as required by Salmon."""
    seen = set()
    decoys = []
    with open(combined, "w") as out:
        for source, is_genome in ((transcriptome, False), (genome, True)):
            count = 0
            with open(source) as inp:
                for line in inp:
                    if line.startswith(">"):
                        identifier = line[1:].split()[0]
                        if "|" in identifier or identifier in seen:
                            raise ValueError(f"Ambiguous/duplicate FASTA identifier: {identifier}")
                        seen.add(identifier)
                        count += 1
                        if is_genome:
                            decoys.append(identifier)
                    out.write(line.rstrip("\n") + "\n")
            if not count:
                raise ValueError(f"Empty FASTA reference: {source}")
    Path(names).write_text("\n".join(decoys) + "\n")


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
    ap.add_argument("--genome", default="")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--key-inputs", default="")
    for opt in ("organism", "accession", "assembly-name", "transcriptome-url",
                "gtf-url", "genome-url"):
        ap.add_argument(f"--{opt}", default="")
    ap.add_argument("--tax-id", type=int, default=0)
    a = ap.parse_args()
    compatibility = check_identifiers(a.transcriptome, a.gtf, "" if a.decoys == "genome" else a.decoys)
    if a.decoys == "genome" and not a.genome:
        raise ValueError("Genome decoys require --genome")

    final_index = Path(a.final_index)
    staging = Path(a.staging)
    entry = final_index.parent               # the cache entry directory
    lock_out = Path(a.lock_out)

    with reference_cache.BuildLock(final_index):
        # Another run may have finished while we waited for the lock. Its
        # output is as valid as ours would be, so do not rebuild.
        if final_index.is_dir() and (final_index / "info.json").is_file() \
                and lock_out.is_file():
            ref = {"accession": a.accession, "transcriptome_fasta_url": a.transcriptome_url,
                   "gtf_url": a.gtf_url, "genome_fasta_url": a.genome_url}
            problems = reference_cache.verify_entry(entry, ref, a.kmer, a.decoys or None)
            if problems:
                raise RuntimeError("Existing reference is incompatible: " + "; ".join(problems))
            print(f"build_salmon_index: {final_index} was published by another "
                  f"run while waiting for the lock; reusing it")
            return

        # A directory without the completion lock is an interrupted build.
        # We hold its build lock and only remove this owned index, never inputs.
        if final_index.exists():
            shutil.rmtree(final_index)
        if staging.exists():
            shutil.rmtree(staging)
        staging.parent.mkdir(parents=True, exist_ok=True)

        fasta, decoy_file = a.transcriptome, a.decoys
        if a.decoys == "genome":
            fasta, decoy_file = str(entry / "gentrome.fa"), str(entry / "decoys.txt")
            prepare_genome_decoys(a.transcriptome, a.genome, fasta, decoy_file)
            compatibility = check_identifiers(fasta, a.gtf, decoy_file)
        cmd = ["salmon", "index", "-t", fasta, "-i", str(staging),
               "-k", str(a.kmer), "--threads", str(a.threads)]
        if a.decoys:
            cmd += ["-d", decoy_file]
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
            "--decoys", decoy_file,
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
        rec["identifier_compatibility"] = compatibility
        if a.decoys == "genome":
            rec["files"]["genome"] = file_record(a.genome)
            rec["files"]["decoy_names"] = file_record(decoy_file)
            rec["files"]["index_fasta"] = file_record(fasta)
        staged_lock.write_text(json.dumps(rec, indent=2) + "\n")

        entry.mkdir(parents=True, exist_ok=True)
        result = reference_cache.publish(staging, final_index)
        # The lock is the completion record. Publish it only after the index,
        # and retain an identical copy inside the index for recovery.
        temp_lock = lock_out.with_suffix(".json.tmp")
        temp_lock.write_text(json.dumps(rec, indent=2) + "\n")
        temp_lock.replace(lock_out)
        if a.decoys == "genome":
            Path(fasta).unlink()  # reproducible concatenation; retain source references
        print(f"build_salmon_index: {result} -> {final_index}")


if __name__ == "__main__":
    main()

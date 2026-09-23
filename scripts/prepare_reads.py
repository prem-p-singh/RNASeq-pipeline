#!/usr/bin/env python3
"""Validate and merge read units from one compatible library, in mate order."""
import argparse
from contextlib import ExitStack
import gzip
import hashlib
import itertools
import json
import os
from pathlib import Path
import subprocess
import tempfile

from metadata import local_path


def records(path):
    with gzip.open(path, "rb") as handle:
        while True:
            header = handle.readline()
            if not header:
                return
            seq, plus, qual = (handle.readline() for _ in range(3))
            if (not header.startswith(b"@") or not plus.startswith(b"+")
                    or not seq.strip() or len(seq.rstrip(b"\r\n")) != len(qual.rstrip(b"\r\n"))):
                raise ValueError(f"Malformed four-line FASTQ record in {path}")
            tokens = header[1:].split()
            if not tokens:
                raise ValueError(f"Missing FASTQ read ID in {path}")
            ident = tokens[0]
            ident = ident[:-2] if ident.endswith((b"/1", b"/2")) else ident
            yield ident, header + seq + plus + qual


def prepare(manifest, destination):
    paired = manifest["layout"] == "paired"
    if manifest["layout"] not in ("single", "paired") or not manifest.get("units"):
        raise ValueError("Expected a single/paired library with at least one read unit")
    dest = Path(destination).resolve()
    dest.mkdir(parents=True, exist_ok=True)
    outputs = [dest / "prepared_R1.fastq.gz"]
    if paired:
        outputs.append(dest / "prepared_R2.fastq.gz")
    reserved = {p.resolve() for p in outputs + [dest / "read_preparation.json"]}
    for unit in manifest["units"]:
        for role in ("r1", "r2"):
            record = unit.get(role)
            source = local_path(record["uri"]) if record else None
            if source is not None and source.resolve() in reserved:
                raise ValueError(f"Output/source collision: {source}")
    ledger = {"schema_version": 1, "library_id": manifest["library_id"],
              "layout": manifest["layout"], "units": [], "fragments": 0}
    # Keep partial results and downloads private; failure never publishes a ledger.
    (dest / "read_preparation.json").unlink(missing_ok=True)
    with tempfile.TemporaryDirectory(prefix=".prepare-", dir=dest) as tmp, ExitStack() as stack:
        tmp = Path(tmp)
        writers = [stack.enter_context(gzip.open(tmp / p.name, "wb")) for p in outputs]
        seen = set()
        for index, unit in enumerate(manifest["units"]):
            if not unit.get("r1") or bool(unit.get("r2")) != paired:
                raise ValueError(f"Read unit {index}: mate structure disagrees with library layout")
            paths, provenance = [], []
            for role in ("r1", "r2") if paired else ("r1",):
                record = unit[role]
                uri = record["uri"]
                path = local_path(uri)
                if path is None:
                    path = tmp / f"download-{index}-{role}.fq.gz"
                    if uri.startswith("s3://"):
                        subprocess.run(["aws", "s3", "cp", uri, str(path)], check=True)
                    elif uri.startswith(("https://", "http://", "ftp://")):
                        subprocess.run(["curl", "--fail", "--retry", "3", "-sSL", "-o", str(path), uri], check=True)
                    else:
                        raise ValueError(f"Unsupported read URI: {uri}")
                    identity = uri
                else:
                    path = path.resolve()
                    identity = str(path)
                if identity in seen or path in outputs or path == dest / "read_preparation.json":
                    raise ValueError(f"Repeated input or output/source collision: {uri}")
                seen.add(identity)
                with path.open("rb") as handle:
                    digest = hashlib.file_digest(handle, "sha256").hexdigest()
                expected = record.get("checksum", "").removeprefix("sha256:")
                if expected and expected.lower() != digest:
                    raise ValueError(f"Checksum mismatch for {uri}")
                if record.get("bytes") not in (None, "") and int(record["bytes"]) != path.stat().st_size:
                    raise ValueError(f"Byte-size mismatch for {uri}")
                paths.append(path)
                provenance.append({"read_unit_id": record["read_unit_id"], "role": role,
                                   "uri": uri, "sha256": digest, "bytes": path.stat().st_size})
            count = 0
            readers = [records(path) for path in paths]
            try:
                for reads in itertools.zip_longest(*readers):
                    if any(r is None for r in reads) or (paired and reads[0][0] != reads[1][0]):
                        raise ValueError(f"Mismatched mate IDs/counts in run={unit.get('run')} lane={unit.get('lane')}")
                    for writer, read in zip(writers, reads):
                        writer.write(read[1])
                    count += 1
            finally:
                for reader in readers:
                    reader.close()
            if not count:
                raise ValueError(f"Empty input in run={unit.get('run')} lane={unit.get('lane')}")
            ledger["units"].append({"run": unit.get("run"), "lane": unit.get("lane"),
                                     "fragments": count, "inputs": provenance})
            ledger["fragments"] += count
            for path in paths:
                if path.parent == tmp:
                    path.unlink()
        stack.close()
        for output in outputs:
            os.replace(tmp / output.name, output)
        (tmp / "ledger.json").write_text(json.dumps(ledger, indent=2) + "\n")
        os.replace(tmp / "ledger.json", dest / "read_preparation.json")
    return ledger


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest-json", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()
    prepare(json.loads(args.manifest_json), args.outdir)

#!/usr/bin/env python3
"""reference_cache.py — content identity, build locking, atomic publication.

Master plan 5.3: "Shared caches use content/version identity, an ordinary
filesystem build lock, a temporary build directory, validation, and atomic
publication. Assembly accession alone is insufficient: annotation, decoys, tool
version, and index parameters can differ. Readers never consume an incomplete
cache entry."

RF09: "Cache entry found | Reuse only if content, construction parameters and
tool compatibility match; accession name alone is insufficient."

The cache was keyed on accession alone, so two projects that agreed on the
genome but differed in annotation release, k-mer or decoy status would share
one directory and silently consume each other's index.

Split of responsibility, because not everything is knowable at the same time:

  cache_key()      construction parameters known when the DAG is built, so the
                   path itself distinguishes incompatible builds
  reference.lock   what actually built the entry, written during the build
  verify_entry()   read-time check that a found entry is complete and was
                   built from the parameters now being asked for

`flock` is absent on macOS, so the lock is a mkdir, which is atomic on every
POSIX filesystem. A stale lock is reported with its recorded owner rather than
broken automatically: master plan 9.3 says stale-lock recovery "requires
recorded ownership checks".
"""
from __future__ import annotations

import hashlib
import json
import os
import shutil
import socket
import time
from pathlib import Path

# Bumped when the meaning of a key input changes, so old entries do not appear
# compatible with a new construction scheme.
KEY_SCHEMA = 1


def key_inputs(reference: dict, kmer: int = 31, decoys: str | None = None) -> dict:
    """The construction parameters that make two indexes interchangeable.

    Accession is included for readability of the resulting path, but it is the
    URLs, k-mer and decoy status that actually determine what was built.
    """
    reference = reference or {}
    return {
        "key_schema": KEY_SCHEMA,
        "accession": reference.get("accession") or "unspecified_assembly",
        "transcriptome_fasta_url": reference.get("transcriptome_fasta_url") or "",
        "gtf_url": reference.get("gtf_url") or "",
        "kmer": int(kmer),
        "decoys": decoys or "none",
    }


def cache_key(reference: dict, kmer: int = 31, decoys: str | None = None) -> str:
    """Directory name for this exact construction: <accession>-<digest>.

    The accession stays in the name so the cache is readable by a human, and
    the digest is what actually separates incompatible builds.
    """
    inputs = key_inputs(reference, kmer, decoys)
    payload = json.dumps(inputs, sort_keys=True).encode()
    digest = hashlib.sha256(payload).hexdigest()[:12]
    acc = str(inputs["accession"]).replace("/", "_")
    return f"{acc}-{digest}"


# --------------------------------------------------------------- completeness
REQUIRED_ENTRY_PARTS = ("reference.lock.json", "transcriptome.fa",
                        "annotation.gtf", "salmon_idx/info.json")


def entry_problems(entry_dir) -> list[str]:
    """What is missing or empty in a cache entry. Empty list means usable."""
    entry = Path(entry_dir)
    out = []
    for rel in REQUIRED_ENTRY_PARTS:
        p = entry / rel
        if not p.exists():
            out.append(f"missing {rel}")
        elif p.is_file() and p.stat().st_size == 0:
            out.append(f"empty {rel}")
    return out


def verify_entry(entry_dir, reference: dict, kmer: int = 31,
                 decoys: str | None = None) -> list[str]:
    """RF09: reuse only if the entry is complete and matches the request.

    Returns reasons not to reuse. An entry whose lock records different
    construction parameters is refused even when the path matches, which is the
    case a hand-edited or partially-restored cache produces.
    """
    problems = entry_problems(entry_dir)
    if problems:
        return problems

    lock_path = Path(entry_dir) / "reference.lock.json"
    try:
        lock = json.loads(lock_path.read_text())
    except (OSError, ValueError) as e:
        return [f"reference.lock.json unreadable: {e}"]

    recorded = (lock.get("cache") or {}).get("key_inputs")
    if recorded is None:
        return ["reference.lock.json records no cache.key_inputs"]

    wanted = key_inputs(reference, kmer, decoys)
    diffs = [f"{k}: entry has {recorded.get(k)!r}, request wants {v!r}"
             for k, v in wanted.items() if recorded.get(k) != v]
    return diffs


# ---------------------------------------------------------------------- lock
class BuildLock:
    """Atomic mkdir lock with a recorded owner.

    Used as a context manager around a cache build. Waits for a concurrent
    build rather than duplicating it, and never removes another holder's lock.
    """

    def __init__(self, target_dir, timeout=3600, poll=2.0):
        self.path = Path(str(target_dir) + ".lock")
        self.timeout = timeout
        self.poll = poll
        self.acquired = False

    def owner(self) -> dict:
        try:
            return json.loads((self.path / "owner.json").read_text())
        except (OSError, ValueError):
            return {}

    def __enter__(self):
        deadline = time.time() + self.timeout
        while True:
            try:
                self.path.mkdir(parents=True)          # atomic
                (self.path / "owner.json").write_text(json.dumps({
                    "pid": os.getpid(),
                    "host": socket.gethostname(),
                    "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ",
                                                 time.gmtime()),
                }, indent=2))
                self.acquired = True
                return self
            except FileExistsError:
                if time.time() >= deadline:
                    o = self.owner()
                    raise TimeoutError(
                        f"reference cache build lock held for longer than "
                        f"{self.timeout}s: {self.path}\n"
                        f"  holder: pid={o.get('pid')} host={o.get('host')} "
                        f"since={o.get('started_utc')}\n"
                        f"  Not removed automatically. Confirm that build is "
                        f"dead, then delete the directory to recover.")
                time.sleep(self.poll)

    def __exit__(self, *exc):
        if self.acquired:
            shutil.rmtree(self.path, ignore_errors=True)
        return False


def publish(staging_dir, final_dir):
    """Move a completed build into place, or discard it if we lost the race.

    os.rename is atomic within a filesystem, so a reader sees either no entry
    or a complete one, never a half-written directory.
    """
    staging, final = Path(staging_dir), Path(final_dir)
    final.parent.mkdir(parents=True, exist_ok=True)
    try:
        os.rename(staging, final)
        return "published"
    except OSError:
        # Another build published first. Theirs is as valid as ours.
        if final.exists():
            shutil.rmtree(staging, ignore_errors=True)
            return "already_present"
        raise


def demo():
    import tempfile
    ref = {"accession": "GCF_1.1",
           "transcriptome_fasta_url": "https://x/rna.fa.gz",
           "gtf_url": "https://x/ann.gtf.gz"}

    # accession alone must not decide the key
    k1 = cache_key(ref, kmer=31)
    assert cache_key(ref, kmer=25) != k1, "k-mer must change the key"
    assert cache_key(dict(ref, gtf_url="https://x/other.gtf.gz")) != k1, \
        "annotation release must change the key"
    assert cache_key(ref, decoys="genome.fa") != k1, "decoys must change the key"
    assert cache_key(ref, kmer=31) == k1, "same inputs must give the same key"
    assert k1.startswith("GCF_1.1-"), k1

    d = Path(tempfile.mkdtemp())
    entry = d / k1
    (entry / "salmon_idx").mkdir(parents=True)
    assert entry_problems(entry), "incomplete entry must be refused"

    for rel, text in (("transcriptome.fa", ">t\nACGT\n"),
                      ("annotation.gtf", "c\tt\texon\t1\t4\t.\t+\t.\tgene_id \"g\";\n"),
                      ("salmon_idx/info.json", "{}")):
        (entry / rel).write_text(text)
    (entry / "reference.lock.json").write_text(json.dumps(
        {"cache": {"key_inputs": key_inputs(ref, 31)}}))
    assert entry_problems(entry) == [], entry_problems(entry)
    assert verify_entry(entry, ref, 31) == []
    assert verify_entry(entry, ref, 25), "a k-mer mismatch must refuse reuse"

    with BuildLock(entry, timeout=5) as lk:
        assert lk.path.is_dir()
        inner = BuildLock(entry, timeout=1, poll=0.1)
        try:
            inner.__enter__()
            raise AssertionError("a second holder acquired the same lock")
        except TimeoutError:
            pass
    assert not (Path(str(entry) + ".lock")).exists(), "lock not released"

    print("reference_cache.demo: ok")


if __name__ == "__main__":
    demo()

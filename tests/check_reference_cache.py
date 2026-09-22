#!/usr/bin/env python3
"""Self-check for the shared reference cache (P0 / R03, R17d; V19).

Master plan 5.3: "Shared caches use content/version identity, an ordinary
filesystem build lock, a temporary build directory, validation, and atomic
publication. Assembly accession alone is insufficient: annotation, decoys, tool
version, and index parameters can differ. Readers never consume an incomplete
cache entry." RF09 says the same for reuse.

The defect: the cache was keyed on accession alone, so two projects agreeing on
the genome but differing in annotation release, k-mer or decoy status shared one
directory and consumed each other's index. `salmon index` also wrote straight
into that shared path, so an interrupted build left a partial index where every
project reads, and two runs on the same accession raced.

What must not break:
  - the key separates builds that differ in any construction parameter
  - an incomplete entry is refused, not consumed
  - an entry built from other parameters is refused even at a matching path
  - two builders cannot hold the lock at once, and the loser waits
  - a stale lock is reported with its owner, never broken silently
  - publication is atomic, and losing the race discards the loser's work
  - the Snakefile keys the cache by content, not by accession

Run:  python3 tests/check_reference_cache.py
"""
import importlib.util
import json
import multiprocessing as mp
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
spec = importlib.util.spec_from_file_location(
    "reference_cache", ROOT / "scripts" / "reference_cache.py")
rc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rc)

REF = {"accession": "GCF_1.1",
       "transcriptome_fasta_url": "https://x/rna.fa.gz",
       "gtf_url": "https://x/ann.gtf.gz"}

# --- 1. content identity, not accession -------------------------------
base = rc.cache_key(REF, 31)
assert base.startswith("GCF_1.1-"), base
variants = {
    "different k-mer":     rc.cache_key(REF, 25),
    "different annotation": rc.cache_key(dict(REF, gtf_url="https://x/v2.gtf.gz"), 31),
    "different transcriptome": rc.cache_key(dict(REF, transcriptome_fasta_url="https://x/v2.fa.gz"), 31),
    "decoys added":        rc.cache_key(REF, 31, "genome.fa"),
}
for label, k in variants.items():
    assert k != base, f"{label} produced the same cache key"
assert len(set(variants.values())) == len(variants), "variants collided"
assert rc.cache_key(REF, 31) == base, "same inputs must be stable"

# two accessions that differ only by accession still differ
assert rc.cache_key(dict(REF, accession="GCF_2.2"), 31) != base


def make_entry(root, ref=REF, kmer=31, complete=True, key_inputs=None):
    e = Path(root) / rc.cache_key(ref, kmer)
    (e / "salmon_idx").mkdir(parents=True, exist_ok=True)
    (e / "transcriptome.fa").write_text(">t\nACGT\n")
    (e / "annotation.gtf").write_text('c\tt\texon\t1\t4\t.\t+\t.\tgene_id "g";\n')
    if complete:
        (e / "salmon_idx" / "info.json").write_text('{"index_version":5}')
        (e / "reference.lock.json").write_text(json.dumps(
            {"cache": {"key_inputs": key_inputs or rc.key_inputs(ref, kmer)}}))
    return e


# --- 2. an incomplete entry is refused --------------------------------
d = Path(tempfile.mkdtemp())
partial = make_entry(d, complete=False)
probs = rc.entry_problems(partial)
assert probs, "incomplete entry accepted"
assert any("reference.lock.json" in p for p in probs), probs
assert any("info.json" in p for p in probs), probs

# a zero-byte index sentinel is incomplete too
full = make_entry(Path(tempfile.mkdtemp()))
assert rc.entry_problems(full) == []
(full / "salmon_idx" / "info.json").write_text("")
assert any("empty" in p for p in rc.entry_problems(full)), rc.entry_problems(full)

# --- 3. reuse is refused when parameters disagree ---------------------
d2 = Path(tempfile.mkdtemp())
e = make_entry(d2)
assert rc.verify_entry(e, REF, 31) == []
diffs = rc.verify_entry(e, REF, 25)
assert diffs and any("kmer" in x for x in diffs), diffs

# an entry whose lock was written for other parameters, at a matching path
e2 = make_entry(Path(tempfile.mkdtemp()),
                key_inputs=rc.key_inputs(dict(REF, gtf_url="https://x/other.gtf.gz"), 31))
diffs = rc.verify_entry(e2, REF, 31)
assert diffs and any("gtf_url" in x for x in diffs), diffs

# a lock with no cache record cannot be verified, so it is not reused
e3 = make_entry(Path(tempfile.mkdtemp()))
(e3 / "reference.lock.json").write_text(json.dumps({"schema_version": 1}))
assert rc.verify_entry(e3, REF, 31), "lock without key_inputs was accepted"


# --- 4. the build lock actually excludes (R17d, V19) ------------------
def hold(target, seconds, started, done):
    with rc.BuildLock(target, timeout=30):
        started.set()
        time.sleep(seconds)
    done.set()


if __name__ == "__main__":
    d3 = Path(tempfile.mkdtemp())
    target = d3 / "entry"
    started, done = mp.Event(), mp.Event()
    p = mp.Process(target=hold, args=(target, 1.5, started, done))
    p.start()
    assert started.wait(10), "first holder never acquired the lock"

    # while held, a second builder must not get in
    t0 = time.time()
    try:
        with rc.BuildLock(target, timeout=0.4, poll=0.05):
            raise AssertionError("two builders held the same lock")
    except TimeoutError as e:
        msg = str(e)
    assert "held for longer than" in msg, msg
    assert "pid=" in msg and "host=" in msg, f"owner not reported: {msg}"
    assert "Not removed automatically" in msg, "a stale lock must not be broken silently"

    # and it becomes available once released, rather than staying stuck
    p.join(20)
    assert done.is_set(), "holder did not finish"
    with rc.BuildLock(target, timeout=5):
        pass
    assert not Path(str(target) + ".lock").exists(), "lock left behind"

    # --- 5. publication is atomic; the loser discards its work --------
    d4 = Path(tempfile.mkdtemp())
    staging, final = d4 / "staging", d4 / "final"
    staging.mkdir()
    (staging / "info.json").write_text("{}")
    assert rc.publish(staging, final) == "published"
    assert (final / "info.json").is_file() and not staging.exists()

    # a second builder finds the entry already there and drops its staging
    staging2 = d4 / "staging2"
    staging2.mkdir()
    (staging2 / "info.json").write_text("{}")
    assert rc.publish(staging2, final) == "already_present"
    assert not staging2.exists(), "loser left its staging directory behind"
    assert (final / "info.json").is_file(), "winner's entry was disturbed"

    # --- 6. the Snakefile keys by content -----------------------------
    snake = (ROOT / "Snakefile").read_text()
    assert "reference_cache" in snake, "Snakefile does not use the cache module"
    assert "CACHE_KEY" in snake and "cache_key(" in snake, snake[:0]
    assert "/ str(_accession)" not in snake, (
        "Snakefile still keys the cache on the bare accession")

    print("check_reference_cache.py: all assertions passed")

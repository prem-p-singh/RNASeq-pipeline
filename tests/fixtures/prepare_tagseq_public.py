#!/usr/bin/env python3
"""Fetch pinned public runs for one QuantSeq kit; does not run analysis.

Usage: prepare_tagseq_public.py KIT NEW_DESTINATION
  KIT is one of the keys in tagseq_public.json's "kits" object
  ("QuantSeq FWD", "QuantSeq REV", "QuantSeq FWD-UMI").
"""
import gzip
import hashlib
import json
import random
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts'))
from tagseq import check_spacer  # noqa: E402

MAX_READS = 2_000_000
SEED = 1


def digest(path, algorithm='sha256'):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, algorithm).hexdigest()


def fetch(url, target, expected_md5, expected_bytes):
    subprocess.run(['curl', '--fail', '--retry', '3', '-sSL', url, '-o', str(target)], check=True)
    if digest(target, 'md5') != expected_md5 or target.stat().st_size != expected_bytes:
        raise ValueError(f'Public source identity mismatch: {target.name}')
    return dict(url=url, bytes=target.stat().st_size, sha256=digest(target))


def downsample(src, dst, max_reads=MAX_READS, seed=SEED):
    """Seeded reservoir sampling (Algorithm R) over whole 4-line FASTQ records.
    Deterministic for a given seed; keeps <= max_reads records; input order is
    preserved when the source has <= max_reads records to begin with.
    """
    rng = random.Random(seed)
    reservoir = []
    seen = 0
    record = []
    with gzip.open(src, 'rt') as fin:
        for line in fin:
            record.append(line)
            if len(record) == 4:
                if seen < max_reads:
                    reservoir.append(record)
                else:
                    j = rng.randint(0, seen)
                    if j < max_reads:
                        reservoir[j] = record
                record = []
                seen += 1
    if record:
        raise ValueError(f'{src.name}: FASTQ record count not a multiple of 4')
    with gzip.open(dst, 'wt') as fout:
        for rec in reservoir:
            fout.writelines(rec)
    return dict(reads_in=seen, reads_kept=len(reservoir))


def prepare(kit, destination):
    manifest = json.loads(Path(__file__).with_name('tagseq_public.json').read_text())
    if kit not in manifest['kits']:
        raise ValueError(f'Unknown kit {kit!r}; choose one of {sorted(manifest["kits"])}')
    spec = manifest['kits'][kit]
    out = Path(destination).resolve()
    out.mkdir(parents=True, exist_ok=False)
    raw_dir = out / 'raw'
    raw_dir.mkdir()
    ledger = out / 'download_manifest.json'
    evidence = dict(status='preparing', kit=kit, study=spec['study'], downsample=manifest['downsample'], runs={})
    try:
        for run in spec['runs']:
            acc = run['accession']
            raw = raw_dir / f'{acc}.fastq.gz'
            source = fetch(run['fastq_url'], raw, run['fastq_md5'], run['fastq_bytes'])
            derived_path = out / f'{acc}.fastq.gz'
            counts = downsample(raw, derived_path)
            entry = dict(group=run['group'], sample_title=run['sample_title'],
                         source=source, derived=dict(sha256=digest(derived_path), **counts))
            if kit == 'QuantSeq FWD-UMI':
                entry['umi_spacer_fraction'] = check_spacer(derived_path)
            evidence['runs'][acc] = entry
        samples_tsv = out / 'samples.tsv'
        with samples_tsv.open('w') as f:
            f.write('sample_id\tbiological_unit_id\tcondition\tfastq_path\n')
            for run in spec['runs']:
                acc = run['accession']
                f.write(f"{acc}\t{acc}\t{run['group']}\t{out / (acc + '.fastq.gz')}\n")
        (out / 'profile.json').write_text(json.dumps(spec, indent=2) + '\n')
        evidence['status'] = 'downloaded_not_qualified'
    except Exception as exc:
        evidence.update(status='failed', error=str(exc))
        raise
    finally:
        ledger.write_text(json.dumps(evidence, indent=2) + '\n')
    return evidence


if __name__ == '__main__':
    if len(sys.argv) != 3:
        raise SystemExit('Usage: prepare_tagseq_public.py KIT NEW_DESTINATION')
    result = prepare(sys.argv[1], sys.argv[2])
    print(result['status'])

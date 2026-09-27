#!/usr/bin/env python3
"""Fetch pinned technical controls for the plant small-RNA route; does not run analysis."""
import csv
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys
import urllib.request

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from small_rna import mature_reference

ENA_TEMPLATE = 'https://www.ebi.ac.uk/ena/portal/api/filereport?accession={0}&result=read_run&fields=fastq_ftp,fastq_md5,fastq_bytes'


def digest(path, algorithm='sha256'):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, algorithm).hexdigest()


def fetch(url, target, expected, algorithm='sha256', size=None):
    subprocess.run(['curl', '--fail', '--retry', '3', '-sSL', url, '-o', str(target)], check=True)
    if digest(target, algorithm) != expected or (size is not None and target.stat().st_size != size):
        raise ValueError(f'Public source identity mismatch: {target.name}')
    return dict(url=url, bytes=target.stat().st_size, sha256=digest(target))


def fetch_record(url, target):
    """Download a source whose checksum is not pinned yet; record it into the manifest instead."""
    subprocess.run(['curl', '--fail', '--retry', '3', '-sSL', url, '-o', str(target)], check=True)
    return dict(url=url, bytes=target.stat().st_size, sha256=digest(target))


def ena_fastq_files(accession):
    """Resolve this run's FASTQ URL(s), md5(s) and byte size(s) from the ENA filereport API."""
    with urllib.request.urlopen(ENA_TEMPLATE.format(accession)) as resp:
        text = resp.read().decode()
    rows = list(csv.DictReader(io.StringIO(text), delimiter='\t'))
    if len(rows) != 1 or not rows[0].get('fastq_ftp'):
        raise ValueError(f'{accession}: ENA filereport returned no single-run fastq_ftp')
    ftps, md5s = rows[0]['fastq_ftp'].split(';'), rows[0]['fastq_md5'].split(';')
    sizes = [int(x) for x in rows[0]['fastq_bytes'].split(';')]
    if not (len(ftps) == len(md5s) == len(sizes)):
        raise ValueError(f'{accession}: ENA filereport fastq_ftp/md5/bytes counts disagree')
    return list(zip(ftps, md5s, sizes))


def prepare(destination):
    manifest = json.loads(Path(__file__).with_name('plant_small_rna_public.json').read_text())
    out = Path(destination).resolve()
    out.mkdir(parents=True, exist_ok=False)
    evidence = dict(status='preparing', biological_inference_eligible=False, sources={}, runs={})
    ledger = out / 'download_manifest.json'
    try:
        evidence['sources']['genome.fa'] = fetch_record(manifest['genome_source']['url'], out / 'genome.fa')
        prefix = manifest['mature_mirna_source']['prefix']
        raw = out / 'mature.fa.all'
        evidence['sources']['mature.fa.all'] = fetch_record(manifest['mature_mirna_source']['url'], raw)
        n = mature_reference(raw, prefix, out / 'ath_mature.fa')
        evidence['reference_derivation'] = dict(source='miRBase mature.fa', prefix=prefix, sequences=n,
                                                 fasta_sha256=digest(out / 'ath_mature.fa'))
        for run in manifest['runs']:
            acc = run['accession']
            files = ena_fastq_files(acc)
            names = [f'{acc}.fastq.gz'] if len(files) == 1 else [f'{acc}_{i + 1}.fastq.gz' for i in range(len(files))]
            for name, (ftp, md5, size) in zip(names, files):
                url = ftp if ftp.startswith('http') else 'https://' + ftp
                evidence['sources'][name] = fetch(url, out / name, md5, 'md5', size)
            evidence['runs'][acc] = names
        (out / 'profile.json').write_text(json.dumps(manifest, indent=2) + '\n')
        evidence['status'] = 'downloaded_not_qualified'
    except Exception as exc:
        evidence.update(status='failed', error=str(exc))
        raise
    finally:
        ledger.write_text(json.dumps(evidence, indent=2) + '\n')
    return evidence


if __name__ == '__main__':
    if sys.argv[1:] == ['--help']:
        print('Usage: prepare_plant_small_rna_public.py NEW_DESTINATION')
        sys.exit(0)
    if len(sys.argv) != 2:
        raise SystemExit('Usage: prepare_plant_small_rna_public.py NEW_DESTINATION')
    result = prepare(sys.argv[1])
    print(result['status'])

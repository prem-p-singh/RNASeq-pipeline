#!/usr/bin/env python3
"""Fetch pinned technical controls for small-RNA development; does not run analysis."""
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import openpyxl


def digest(path, algorithm='sha256'):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, algorithm).hexdigest()


def fetch(url, target, expected, algorithm='sha256', size=None):
    subprocess.run(['curl', '--fail', '--retry', '3', '-sSL', url, '-o', str(target)], check=True)
    if digest(target, algorithm) != expected or (size is not None and target.stat().st_size != size):
        raise ValueError(f'Public source identity mismatch: {target.name}')
    return dict(url=url, bytes=target.stat().st_size, sha256=digest(target))


def prepare(destination):
    manifest = json.loads(Path(__file__).with_name('small_rna_public.json').read_text())
    out = Path(destination).resolve()
    out.mkdir(parents=True, exist_ok=False)
    evidence = dict(status='preparing', biological_inference_eligible=False, sources={})
    ledger = out/'download_manifest.json'
    try:
        for source, filename in zip(manifest['sources'], ['study.soft.gz', 'supplement.xlsx']):
            evidence['sources'][filename] = fetch(source['url'], out/filename, source['sha256'])
        wb = openpyxl.load_workbook(out/'supplement.xlsx', read_only=True, data_only=True)
        try:
            rows = list(wb['Table S2'].values)[3:]
        finally:
            wb.close()
        if len(rows) != 334:
            raise ValueError('Unexpected published reference-table size')
        references, aliases = {}, {}
        for ident, alias, sequence, length, kind, modification, a, b in rows:
            if not ident or any(c.isspace() for c in ident) or set(sequence)-set('ACGT') or len(sequence) != length:
                raise ValueError('Invalid published sequence record')
            record = (sequence, float(a.removesuffix('x')), float(b.removesuffix('x')))
            if ident in references and references[ident] != record:
                raise ValueError('Conflicting sequence or concentration under one published ID')
            references[ident] = record
            aliases.setdefault(ident, []).append(dict(alias=alias, kind=kind, modification=modification))
        if len(references) != 331 or len({v[0] for v in references.values()}) != 331:
            raise ValueError('Unexpected published duplicate structure')
        (out/'reference.fa').write_text(''.join(f'>{key}\n{value[0]}\n' for key,value in references.items()))
        with (out/'expected_ratios.tsv').open('w') as handle:
            writer=csv.writer(handle,delimiter='\t');writer.writerow(['feature_id','length','pool_A','pool_B','expected_A_over_B'])
            for key,(sequence,a,b) in references.items():
                writer.writerow([key,len(sequence),a,b,a/b])
        evidence['reference_derivation'] = dict(source_table='Table S2', input_rows=334, unique_ids_and_sequences=331,
            duplicate_records={k:v for k,v in aliases.items() if len(v)>1},
            explanation='Three repeated IDs have identical sequence and ratio but differing alias/biotype labels. '
                        'Each sequence is represented once; all source labels are preserved here. '
                        'No biological identity is inferred from these conflicting labels.',
            fasta_sha256=digest(out/'reference.fa'), ratio_sha256=digest(out/'expected_ratios.tsv'))
        for run in manifest['runs']:
            key=run['accession']
            evidence['sources'][key+'.fastq.gz']=fetch(run['fastq_url'],out/(key+'.fastq.gz'),
                run['fastq_md5'],'md5',run['fastq_bytes'])
            evidence['sources'][key+'.published_counts.txt.gz']=fetch(run['native_counts_url'],
                out/(key+'.published_counts.txt.gz'),run['native_counts_sha256'])
        (out/'profile.json').write_text(json.dumps(manifest,indent=2)+'\n')
        evidence['status']='downloaded_not_qualified'
    except Exception as exc:
        evidence.update(status='failed',error=str(exc));raise
    finally:
        ledger.write_text(json.dumps(evidence,indent=2)+'\n')
    return evidence


if __name__ == '__main__':
    if len(sys.argv)!=2:
        raise SystemExit('Usage: prepare_small_rna_public.py NEW_DESTINATION')
    result=prepare(sys.argv[1]);print(result['status'])

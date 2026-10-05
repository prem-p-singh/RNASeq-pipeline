#!/usr/bin/env python3
"""Fetch the pinned public runs for one dual-organism study; does not run analysis.

Usage: prepare_dual_public.py STUDY NEW_DESTINATION
  STUDY is one of the keys in dual_public.json's "studies" object ("virus", "cellular").
"""
import csv
import gzip
import io
import json
import sys
import urllib.request
from pathlib import Path

from prepare_tagseq_public import digest, downsample, fetch

ENA_REPORT = 'https://www.ebi.ac.uk/ena/portal/api/filereport?accession={0}&result=read_run&fields=fastq_ftp,fastq_md5,fastq_bytes'


def ena_md5(run):
    """MD5 values from the ENA file report, for a run whose manifest entry has none.

    The file report must name the pinned URLs and byte sizes, or the run is not the one that was selected.
    """
    with urllib.request.urlopen(ENA_REPORT.format(run['accession'])) as response:
        rows = list(csv.DictReader(io.StringIO(response.read().decode()), delimiter='\t'))
    if len(rows) != 1:
        raise ValueError(f"{run['accession']}: the ENA file report did not return one run")
    urls = ['https://' + u for u in rows[0]['fastq_ftp'].split(';')]
    sizes = [int(x) for x in rows[0]['fastq_bytes'].split(';')]
    if urls != run['fastq_url'] or sizes != run['fastq_bytes']:
        raise ValueError(f"{run['accession']}: the ENA file report differs from the pinned URLs or sizes: {urls} {sizes}")
    return rows[0]['fastq_md5'].split(';')


def as_is(src, dst):
    """No read cap: move the verified download into place unchanged and count its records."""
    with gzip.open(src, 'rb') as f:
        lines = sum(1 for _ in f)
    if lines % 4:
        raise ValueError(f'{src.name}: FASTQ line count not a multiple of 4')
    src.replace(dst)
    return dict(reads_in=lines // 4, reads_kept=lines // 4)


def exon_annotation(source_gtf, genome, out_gtf, out_fasta):
    """One exon row and one transcript sequence per gene row of a GTF that has no exon rows.

    For a RefSeq virus annotation (gene, CDS and codon rows only). Returns the number of genes.
    """
    contigs, name = {}, None
    with gzip.open(genome, 'rt') as f:
        for line in f:
            if line.startswith('>'):
                name = line[1:].split()[0]
                contigs[name] = []
            else:
                contigs[name].append(line.strip())
    contigs = {k: ''.join(v) for k, v in contigs.items()}
    complement = str.maketrans('ACGTacgt', 'TGCAtgca')
    rows, sequences = [], []
    with gzip.open(source_gtf, 'rt') as f:
        for line in f:
            fields = line.rstrip('\n').split('\t')
            if line.startswith('#') or len(fields) != 9:
                continue
            if fields[2] == 'exon':
                raise ValueError(f'{source_gtf} already has exon rows; use it as it is')
            if fields[2] != 'gene':
                continue
            gene = fields[8].split('gene_id "')[1].split('"')[0]
            start, end = int(fields[3]), int(fields[4])
            sequence = contigs[fields[0]][start - 1:end]
            if len(sequence) != end - start + 1:
                raise ValueError(f'{gene}: {fields[0]}:{start}-{end} is outside the genome sequence')
            if fields[6] == '-':
                sequence = sequence.translate(complement)[::-1]
            rows.append('\t'.join(fields[:2] + ['exon'] + fields[3:8] + [f'gene_id "{gene}"; transcript_id "{gene}_tx";']) + '\n')
            sequences.append(f'>{gene}_tx\n{sequence}\n')
    if not rows:
        raise ValueError(f'No gene rows in {source_gtf}')
    Path(out_gtf).write_text(''.join(rows))
    Path(out_fasta).write_text(''.join(sequences))
    return len(rows)


def prepare(study, destination):
    manifest = json.loads(Path(__file__).with_name('dual_public.json').read_text())
    if study not in manifest['studies']:
        raise ValueError(f'Unknown study {study!r}; choose one of {sorted(manifest["studies"])}')
    spec = manifest['studies'][study]
    out = Path(destination).resolve()
    out.mkdir(parents=True, exist_ok=False)
    raw_dir = out / 'raw'
    raw_dir.mkdir()
    ledger = out / 'download_manifest.json'
    evidence = dict(status='preparing', key=study, study=spec['study'], downsample=spec['downsample'], references={}, runs={})
    try:
        source_gtf = spec['pathogen'].get('source_gtf')
        if source_gtf:
            files = {}
            for name in ('genome_fasta', 'source_gtf'):
                pinned = spec['pathogen'][name]
                files[name] = raw_dir / pinned['url'].rsplit('/', 1)[1]
                evidence['references'][name] = fetch(pinned['url'], files[name], pinned['md5'], pinned['bytes'])
            genes = exon_annotation(files['source_gtf'], files['genome_fasta'], out / 'pathogen.exons.gtf', out / 'pathogen.transcripts.fa')
            evidence['references']['derived'] = dict(genes=genes, gtf_sha256=digest(out / 'pathogen.exons.gtf'),
                                                     transcripts_sha256=digest(out / 'pathogen.transcripts.fa'))
        rows = []
        for run in spec['runs']:
            acc = run['accession']
            md5s = run['fastq_md5'] or ena_md5(run)
            kept, entry = [], dict(group=run['group'], title=run['title'], files=[])
            if not len(run['fastq_url']) == len(md5s) == len(run['fastq_bytes']):
                raise ValueError(f"{acc}: {len(run['fastq_url'])} URL(s), {len(md5s)} MD5 value(s) and {len(run['fastq_bytes'])} size(s)")
            cap = spec['downsample']['max_reads_per_run']
            for url, md5, size in zip(run['fastq_url'], md5s, run['fastq_bytes']):
                name = url.rsplit('/', 1)[1]
                source = fetch(url, raw_dir / name, md5, size)
                # The same seed keeps the same record numbers in both mates of a pair.
                counts = as_is(raw_dir / name, out / name) if cap is None else downsample(raw_dir / name, out / name, cap, spec['downsample']['seed'])
                entry['files'].append(dict(source=dict(source, md5=md5), derived=dict(sha256=digest(out / name), **counts)))
                kept.append(str(out / name))
            if len({f['derived']['reads_in'] for f in entry['files']}) != 1:
                raise ValueError(f'{acc}: the mate files hold different numbers of records')
            evidence['runs'][acc] = entry
            rows.append([acc, run['group'], *kept, *([''] * (2 - len(kept)))])
        with (out / 'samples.tsv').open('w') as f:
            f.write('sample_id\tcondition\tfastq_path\tfastq_path_r2\n')
            f.writelines('\t'.join(row) + '\n' for row in rows)
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
        raise SystemExit('Usage: prepare_dual_public.py STUDY NEW_DESTINATION')
    result = prepare(sys.argv[1], sys.argv[2])
    print(result['status'])

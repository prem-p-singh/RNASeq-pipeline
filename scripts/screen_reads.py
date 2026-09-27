#!/usr/bin/env python3
"""Diagnostic screening of non-UMI bulk reads against an explicitly declared panel."""
import argparse
import csv
import hashlib
import html
import itertools
import json
from pathlib import Path
import random
import re
import shutil
import subprocess
import tempfile

from prepare_reads import records

SUFFIXES = ('1', '2', '3', '4', 'rev.1', 'rev.2')
DEFAULTS = dict(enabled=False, fragments=100000, seed=1, references=[])


def validate(policy):
    if not isinstance(policy, dict) or set(policy) - set(DEFAULTS):
        raise ValueError('Unknown screening settings')
    p = DEFAULTS | policy
    if type(p['enabled']) is not bool:
        raise ValueError('screening.enabled must be boolean')
    for name in ('fragments', 'seed'):
        if type(p[name]) is not int or p[name] < (1 if name == 'fragments' else 0):
            raise ValueError(f'screening.{name} must be a positive count or nonnegative seed')
    if not isinstance(p['references'], list):
        raise ValueError('screening.references must be a list')
    names = set()
    for ref in p['references']:
        if not isinstance(ref, dict) or set(ref) != {'name', 'role', 'index'}:
            raise ValueError('Each screen reference needs exactly name, role and index')
        if not isinstance(ref['name'], str) or not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*', ref['name']) or ref['name'] in names:
            raise ValueError('Screen reference names must be unique alphanumeric identifiers')
        names.add(ref['name'])
        if ref['role'] not in ('expected', 'possible_contaminant'):
            raise ValueError('Screen reference role must be expected or possible_contaminant')
        if not isinstance(ref['index'], str) or not ref['index'].strip():
            raise ValueError('Screen reference index must be a Bowtie2 basename')
    if p['enabled'] and not any(r['role'] == 'expected' for r in p['references']):
        raise ValueError('Screening requires an explicitly declared expected reference')
    return p


def index_files(basename):
    for ext in ('bt2', 'bt2l'):
        files = [Path(f'{basename}.{suffix}.{ext}').expanduser().resolve() for suffix in SUFFIXES]
        if all(p.is_file() and p.stat().st_size for p in files):
            return files
    raise ValueError(f'Missing or incomplete Bowtie2 index: {basename}')


def identity(path):
    with path.open('rb') as handle:
        digest = hashlib.file_digest(handle, 'sha256').hexdigest()
    return dict(path=str(path), bytes=path.stat().st_size, sha256=digest)


def sample_reads(sources, limit, seed, destination):
    """Reservoir sampling keeps mate pairs together; tool classifies each mate separately."""
    rng = random.Random(seed)
    reservoir, total = [], 0
    readers = [records(path) for path in sources]
    try:
        for unit in itertools.zip_longest(*readers):
            if any(x is None for x in unit) or len({x[0] for x in unit}) != 1:
                raise ValueError('Mismatched screen input mate IDs/counts')
            total += 1
            if total <= limit:
                reservoir.append((total, unit))
            else:
                slot = rng.randrange(total)
                if slot < limit:
                    reservoir[slot] = (total, unit)
    finally:
        for reader in readers:
            reader.close()
    if not total:
        raise ValueError('Empty screening input')
    reservoir.sort(key=lambda x: x[0])
    for mate in range(len(sources)):
        with (destination / f'R{mate+1}.fastq').open('wb') as handle:
            for _, unit in reservoir:
                handle.write(unit[mate][1])
    return total, len(reservoir)


def parse_report(path, panel, count):
    lines = path.read_text().splitlines()
    table = [line for line in lines if line and not line.startswith(('#Fastq', '%'))]
    rows = list(csv.DictReader(table, delimiter='\t'))
    if len(rows) != len(panel) or {r['Genome'] for r in rows} != {r['name'] for r in panel}:
        raise ValueError('Screen output does not match declared reference panel')
    result = []
    roles = {r['name']: r['role'] for r in panel}
    for row in rows:
        total = int(row['#Reads_processed'])
        unmapped = int(row['#Unmapped'])
        unique_panel = int(row['#One_hit_one_genome']) + int(row['#Multiple_hits_one_genome'])
        shared_panel = int(row['#One_hit_multiple_genomes']) + int(row['Multiple_hits_multiple_genomes'])
        if total != count or min(unmapped, unique_panel, shared_panel) < 0 or unmapped + unique_panel + shared_panel != total:
            raise ValueError('Screen read-count conservation failed')
        result.append(dict(reference=row['Genome'], role=roles[row['Genome']], reads=total,
                           unmapped=unmapped, exclusive_to_reference=unique_panel,
                           shared_with_other_references=shared_panel))
    return result


def screen(sources, outdir, sample, policy, threads=2):
    p = validate(policy)
    if type(threads) is not int or threads < 1:
        raise ValueError('Screen threads must be a positive integer')
    out = Path(outdir).resolve()
    sources = [Path(x).resolve() for x in sources]
    if len(sources) not in (1, 2):
        raise ValueError('Screening supports one or two bulk read files')
    panels = [(ref, index_files(ref['index'])) for ref in p['references']] if p['enabled'] else []
    if any(path.is_relative_to(out) for path in sources + [f for _, files in panels for f in files]):
        raise ValueError('Screen output must not contain source reads or reference indexes')
    out.mkdir(parents=True, exist_ok=True)
    for name in ('screening.json', 'screening.tsv', 'screening.html', 'screen.log', 'R1_screen.txt', 'R2_screen.txt'):
        (out / name).unlink(missing_ok=True)
    ledger = dict(schema_version=1, sample=sample, status='not_performed', policy=p,
                  input_stage='after_preprocessing', classification_unit='read', reads_removed=0,
                  interpretation='Matches are diagnostic, not organism abundance or proof of contamination. '
                  'Shared sequences remain ambiguous. Missing panel organisms cannot be excluded. '
                  'Bowtie2 is not splice-aware; prefer expected transcriptome indexes for RNA libraries.',
                  results=[])
    if p['enabled']:
        versions = {}
        for tool, required in [('fastq_screen', '0.16.0'), ('bowtie2', '2.5.4')]:
            exe = shutil.which(tool)
            if not exe or not re.fullmatch(r'[A-Za-z0-9_./+-]+', exe):
                raise ValueError(f'{tool} must be installed at a shell-safe path')
            run = subprocess.run([exe, '--version'], check=True, capture_output=True, text=True)
            versions[tool] = (run.stdout + run.stderr).strip()
            if required not in versions[tool].splitlines()[0]:
                raise ValueError(f'Expected {tool} {required}')
        ledger.update(input_files=[identity(path) for path in sources], versions=versions,
                      references=[dict(**ref, files=[identity(f) for f in files]) for ref, files in panels])
        # FastQ Screen constructs shell commands internally. Only generated safe
        # names reach it; original paths (including spaces) are never interpolated.
        with tempfile.TemporaryDirectory(prefix='rnaseq-screen-', dir='/tmp') as temp:
            work = Path(temp)
            config = []
            for i, (ref, files) in enumerate(panels):
                for suffix, path in zip(SUFFIXES, files):
                    (work / f'panel{i}.{suffix}.{path.suffix[1:]}').symlink_to(path)
                config.append(f'DATABASE {ref["name"]} {work}/panel{i}\n')
            (work / 'screen.conf').write_text(''.join(config))
            total, selected = sample_reads(sources, p['fragments'], p['seed'], work)
            ledger.update(fragments_available=total, fragments_sampled=selected,
                          fraction_sampled=selected/total, sampling='uniform reservoir without replacement')
            cmd = ['fastq_screen', '--conf', str(work / 'screen.conf'), '--aligner', 'bowtie2',
                   '--threads', str(threads), '--subset', '0', '--outdir', str(work)]
            cmd += [str(work / f'R{i+1}.fastq') for i in range(len(sources))]
            ledger['command'] = cmd
            with (out / 'screen.log').open('w') as log:
                subprocess.run(cmd, check=True, cwd=work, stdout=log, stderr=subprocess.STDOUT)
            for i in range(len(sources)):
                report = work / f'R{i+1}_screen.txt'
                rows = parse_report(report, p['references'], selected)
                ledger['results'].extend(dict(mate=f'R{i+1}', **row) for row in rows)
                shutil.copyfile(report, out / report.name)
            if ledger['input_files'] != [identity(path) for path in sources] or any(
                ref['files'] != [identity(f) for f in files] for ref, (_, files) in zip(ledger['references'], panels)):
                raise ValueError('Screen input or reference changed during analysis')
        ledger['status'] = 'completed'
    rows = ledger['results']
    columns = ['mate', 'reference', 'role', 'reads', 'unmapped', 'exclusive_to_reference', 'shared_with_other_references']
    with (out / 'screening.tsv').open('w') as handle:
        writer = csv.DictWriter(handle, columns, delimiter='\t'); writer.writeheader(); writer.writerows(rows)
    rendered = '<tr>' + ''.join(f'<th>{html.escape(c)}</th>' for c in columns) + '</tr>'
    for row in rows:
        rendered += '<tr>' + ''.join(f'<td>{html.escape(str(row[c]))}</td>' for c in columns) + '</tr>'
    (out / 'screening.html').write_text('<!doctype html><meta charset="utf-8"><title>Diagnostic screening</title>'
        f'<h1>{html.escape(str(sample))}: {ledger["status"]}</h1><p>{ledger["interpretation"]}</p>'
        f'<p>Reads removed: 0. Sampling: {ledger.get("fragments_sampled", 0)} fragments.</p><table border="1">{rendered}</table>')
    (out / 'screening.json').write_text(json.dumps(ledger, indent=2) + '\n')
    return ledger


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reads', nargs='+', required=True)
    parser.add_argument('--outdir', required=True)
    parser.add_argument('--sample', required=True)
    parser.add_argument('--policy-json', required=True)
    parser.add_argument('--threads', type=int, default=2)
    args = parser.parse_args()
    screen(args.reads, args.outdir, args.sample, json.loads(args.policy_json), args.threads)

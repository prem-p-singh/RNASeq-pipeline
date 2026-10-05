#!/usr/bin/env python3
"""Host and pathogen in one library: combined reference, read labels and splitting."""
import argparse
import copy
import csv
import gzip
import io
import itertools
import json
import os
import re
import subprocess
import tempfile
from collections import Counter
from pathlib import Path

import yaml

import config_resolve
from metadata import read_tsv
from prepare_reads import records

REPO = Path(__file__).resolve().parent.parent
ROLES = ('host', 'pathogen')
LABELS = ('host', 'pathogen', 'ambiguous', 'unassigned')
STAR_EXIT_GRACE = 5  # seconds for a STAR whose output ended to exit
DEFAULTS = dict(multimap_limit=20, ambiguous_warn_fraction=0.05)


def validate(cfg):
    d = dict(DEFAULTS, **(cfg.get('dual_organism') or {}))
    if cfg.get('samples', {}).get('seq_type') in ('tagseq', 'srnaseq'):
        raise ValueError(f"dual_organism does not support seq_type {cfg['samples']['seq_type']}")
    for role in ROLES:
        o = dict(d.get(role) or {})
        if not isinstance(o.get('label'), str) or not re.fullmatch(r'[A-Za-z0-9_.-]+', o['label']):
            raise ValueError(f'dual_organism.{role}.label must be letters, digits, _ . - (no | or spaces)')
        if not o.get('genome_fasta_url'):
            raise ValueError(f'dual_organism.{role}.genome_fasta_url is required')
        if o.get('objectives') is None:
            o['objectives'] = ['gene_expression', 'differential_expression'] if role == 'host' else ['gene_expression']
        d[role] = o
    if d['host']['label'] == d['pathogen']['label']:
        raise ValueError('dual_organism host and pathogen labels must differ')
    for role in ROLES:
        if not d[role].get('gtf_url'):
            raise ValueError(f'dual_organism.{role}.gtf_url is required')
    # Each organism is quantified on its own afterwards; Salmon needs that organism's transcripts.
    quantifier = (cfg.get('analysis') or {}).get('quantifier', 'auto')
    for role in ROLES:
        salmon = quantifier == 'salmon' or (quantifier == 'auto' and 'alignment' not in d[role]['objectives'])
        if salmon and set(d[role]['objectives']) != {'qc'} and not d[role].get('transcriptome_fasta_url'):
            raise ValueError(f'dual_organism.{role}.transcriptome_fasta_url is required on the Salmon route')
    if d['pathogen'].get('min_assigned_reads') is None:
        d['pathogen']['min_assigned_reads'] = 10000
    if type(d['multimap_limit']) is not int or d['multimap_limit'] < 1:
        raise ValueError('dual_organism.multimap_limit must be a positive integer')
    return d


def combine(host_fa, pathogen_fa, host_gtf, pathogen_gtf, out_fa, out_gtf, host_label, pathogen_label):
    counts = {}
    with open(out_fa, 'w') as out:
        for role, fa, lab in (('host', host_fa, host_label), ('pathogen', pathogen_fa, pathogen_label)):
            seen = set()
            with open(fa) as f:
                for line in f:
                    if line.startswith('>'):
                        tokens = line[1:].split()
                        if not tokens:
                            raise ValueError(f'Empty FASTA header in {fa}')
                        if tokens[0] in seen:
                            raise ValueError(f'Duplicate contig {tokens[0]} in the {role} genome')
                        seen.add(tokens[0])
                        out.write(f'>{lab}|{tokens[0]}\n')
                    else:
                        out.write(line if line.endswith('\n') else line + '\n')
            counts[f'{role}_contigs'] = len(seen)
    with open(out_gtf, 'w') as out:
        for gtf, lab in ((host_gtf, host_label), (pathogen_gtf, pathogen_label)):
            if gtf is None:
                continue
            with open(gtf) as f:
                for line in f:
                    line = line if line.endswith('\n') else line + '\n'
                    out.write(line if line.startswith('#') or not line.strip() else f'{lab}|{line}')
    return counts


def label(organisms):
    if organisms == {'host'}:
        return 'host'
    if organisms == {'pathogen'}:
        return 'pathogen'
    return 'ambiguous' if organisms else 'unassigned'


def assign(sam_lines, fastqs, outs, host_label, pathogen_label):
    """One pass over STAR SAM (alignments of a read adjacent, input order) and the FASTQ(s).

    Outputs are written to <name>.tmp and renamed into place only after every check passes.
    """
    roles = {host_label: 'host', pathogen_label: 'pathogen'}
    counts = dict(reads=0, host=0, pathogen=0, ambiguous=0, unassigned=0, too_many_loci=0)
    groups = itertools.groupby((l for l in sam_lines if not l.startswith('@')), key=lambda l: l.split('\t', 1)[0])
    readers = [records(p) for p in fastqs]
    temps = {role: [Path(f'{p}.tmp') for p in paths] for role, paths in outs.items()}
    handles = {role: [] for role in outs}
    try:
        for role, paths in temps.items():
            handles[role] = []
            for p in paths:
                handles[role].append(gzip.open(p, 'wb'))
        for unit in zip(*readers):
            counts['reads'] += 1
            try:
                name, lines = next(groups)
            except StopIteration:
                raise ValueError(f'conservation failure: no alignment record for read {counts["reads"]}') from None
            if any(name != u[0].decode() for u in unit):
                raise ValueError(f'STAR output order mismatch at read {counts["reads"]}: {name}')
            organisms, too_many = set(), False
            for line in lines:
                f = line.rstrip('\n').split('\t')
                if int(f[1]) & 4:
                    too_many = too_many or 'uT:A:3' in f[11:]
                else:
                    prefix = f[2].split('|', 1)[0]
                    if prefix not in roles:
                        raise ValueError(f'Unexpected contig {f[2]} for read {name}: prefix is neither organism label')
                    organisms.add(roles[prefix])
            lab = label(organisms)
            counts[lab] += 1
            counts['too_many_loci'] += lab == 'unassigned' and too_many
            if lab in handles:
                for h, rec in zip(handles[lab], unit):
                    h.write(rec[1])
        if next(groups, None) is not None:
            raise ValueError('conservation failure: alignment records beyond the last read')
        if counts['host'] + counts['pathogen'] + counts['ambiguous'] + counts['unassigned'] != counts['reads']:
            raise ValueError('conservation failure in dual-organism assignment')
        for hs in handles.values():
            for h in hs:
                h.close()
        for role, paths in outs.items():
            for p, t in zip(paths, temps[role]):
                t.replace(p)
    except BaseException:
        for hs in handles.values():
            for h in hs:
                h.close()
        for ts in temps.values():
            for t in ts:
                t.unlink(missing_ok=True)
        raise
    finally:
        for r in readers:
            r.close()
    return counts


def align(index, reads, outdir, sample, host_label, pathogen_label, limit=20, threads=4):
    """Align once to the combined genome, split the reads and write assignment.json last."""
    import star_counts
    index, outdir = Path(index).resolve(), Path(outdir).resolve()
    lock = json.loads((index / 'reference.lock.json').read_text())
    star = star_counts.version('STAR')
    if star != lock['request']['STAR']:
        raise ValueError('STAR version differs from index builder')
    reads = [str(Path(r).resolve()) for r in reads]
    inputs = [star_counts.identity(r) for r in reads]
    outdir.mkdir(parents=True, exist_ok=True)
    outs = {role: [outdir / f'{role}{mate}.fastq.gz' for mate in (('_R1', '_R2') if len(reads) == 2 else ('',))]
            for role in ROLES}
    record = outdir / 'assignment.json'
    for path in [record] + [p for paths in outs.values() for p in paths]:
        path.unlink(missing_ok=True)
    with tempfile.TemporaryDirectory(prefix='.star-', dir=outdir) as tmp:
        cmd = ['STAR', '--runThreadN', str(threads), '--genomeDir', str(index),
               '--readFilesIn', *reads, '--readFilesCommand', 'zcat',
               '--outFilterMultimapNmax', str(limit), '--outSAMmultNmax', '-1',
               '--outSAMtype', 'SAM', '--outStd', 'SAM', '--outSAMorder', 'PairedKeepInputOrder',
               '--outSAMunmapped', 'Within', '--outSAMattributes', 'NH', 'HI', 'AS',
               '--outFileNamePrefix', tmp + '/']
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
        try:
            counts = assign(proc.stdout, reads, outs, host_label, pathogen_label)
        except BaseException as exc:
            # A STAR that failed has ended the stream and exits at once: report it, not the symptom.
            # One still running after the grace period did not fail; stop it and keep the stream error.
            try:
                failed = proc.wait(timeout=STAR_EXIT_GRACE) != 0
            except subprocess.TimeoutExpired:
                proc.terminate()
                proc.wait()
                failed = False
            proc.stdout.close()
            if failed:
                raise RuntimeError(f'STAR exited with status {proc.returncode}') from exc
            raise
        proc.stdout.close()
        if proc.wait() != 0:
            for path in (p for paths in outs.values() for p in paths):
                path.unlink(missing_ok=True)
            raise RuntimeError(f'STAR exited with status {proc.returncode}')
    result = dict(schema_version=1, sample=sample, counts=counts, labels=dict(host=host_label, pathogen=pathogen_label),
                  multimap_limit=limit, versions={'STAR': star}, command=cmd, input_files=inputs,
                  reference_lock_sha256=star_counts.identity(index / 'reference.lock.json')['sha256'])
    staged = Path(f'{record}.tmp')
    staged.write_text(json.dumps(result, indent=2) + '\n')
    staged.replace(record)
    return result


def table(records, warn):
    lines = ['\t'.join(['sample', 'reads', *LABELS, 'too_many_loci', *(f'{k}_fraction' for k in LABELS), 'warnings'])]
    for r in records:
        c = r['counts']
        fractions = [c[k] / c['reads'] if c['reads'] else 0.0 for k in LABELS]
        warnings = f'ambiguous fraction {fractions[2]:.4f} above {warn}' if fractions[2] > warn else ''
        lines.append('\t'.join([r['sample'], *(str(c[k]) for k in ('reads', *LABELS, 'too_many_loci')),
                                *(f'{x:.6f}' for x in fractions), warnings]))
    return '\n'.join(lines) + '\n'


def _from_child(path):
    """A parent-project path as seen from <project>/dual/<role>."""
    path = Path(path).expanduser()
    return str(path if path.is_absolute() else Path('../..') / path)


def child_config(cfg, role, objectives=None):
    """The single-organism project config for one organism of a dual project.

    Reads nothing from the project, so the parent can check before any read is assigned exactly
    what children() writes afterwards.
    """
    full = config_resolve.resolve(REPO, [('project', cfg)])[0]
    organism = validate(full)[role]
    c = {k: copy.deepcopy(v) for k, v in cfg.items() if k != 'dual_organism'}
    for block in ('project', 'analysis', 'reference', 'samples', 'thresholds'):
        c[block] = dict(c.get(block) or {})
    c['project']['output_dir'] = _from_child(Path(full['project']['output_dir']) / role) + '/'
    c['analysis']['objectives'] = objectives or organism['objectives']
    c['preprocessing'] = {'passthrough': True}
    c['screening'] = {'enabled': False}
    c['samples'].update(sheet='config/samples.tsv', metadata_dir=None)
    # A decoy-name file and the gene-information table describe the host and are named relative to the parent.
    host = role == 'host'
    decoys, table_path = full['reference']['decoys'], full['reference']['annotation_tsv']['path']
    if decoys not in (None, '', 'genome'):
        decoys = _from_child(decoys) if host else None
    c['reference'].update({k: organism.get(k) for k in ('genome_fasta_url', 'gtf_url', 'transcriptome_fasta_url')},
                          accession=organism['label'], decoys=decoys,
                          annotation_tsv=dict(c['reference'].get('annotation_tsv') or {},
                                              path=_from_child(table_path) if host and table_path else None))
    if not host:
        # Salmon 1.10.3 can crash in position-bias correction on a very small transcript set.
        c['analysis']['salmon_position_bias'] = False
        # Reads are pre-assigned: inclusion is min_assigned_reads, not the mapping-rate or depth floor.
        c['downstream'] = dict(run_go=False, run_kegg=False, run_wgcna=False)
        c['thresholds']['sample_qc'] = dict(c['thresholds'].get('sample_qc') or {},
                                            mapping_rate_min=0, min_reads_on_genes_rnaseq=0)
    if not c['thresholds']:
        del c['thresholds']
    return c


def _write(path, text):
    """Leave an unchanged file alone, so a relaunch does not look like a new configuration."""
    if not path.is_file() or path.read_text() != text:
        path.write_text(text)


def children(project, cfg, assignment_tsv=None):
    """Write one single-organism project per organism under <project>/dual/ from a finished assignment.

    cfg is the parent project config as authored. The pathogen child takes the samples with at least
    min_assigned_reads pathogen reads; with none, no pathogen project is written.
    """
    project = Path(project).resolve()
    full = config_resolve.resolve(REPO, [('project', cfg)])[0]
    dual = validate(full)
    out = project / Path(full['project']['output_dir']).expanduser()
    sheet = read_tsv(project / full['samples']['sheet'])
    if not sheet:
        raise ValueError(f"No samples in {project / full['samples']['sheet']}")
    with open(assignment_tsv or out / 'assignment.tsv', newline='') as f:
        assigned = {r['sample']: int(r['pathogen']) for r in csv.DictReader(f, delimiter='\t')}
    missing = [r['sample_id'] for r in sheet if r['sample_id'] not in assigned]
    if missing:
        raise ValueError(f"The assignment table has no row for sample(s): {', '.join(missing)}")
    floor = dual['pathogen']['min_assigned_reads']
    low = [r['sample_id'] for r in sheet if assigned[r['sample_id']] < floor]
    objectives = {role: dual[role]['objectives'] for role in ROLES}
    de, reason = 'not_requested', ''
    if set(objectives['pathogen']) & {'differential_expression', 'enrichment'}:
        factor = full['model']['primary_factor']
        kept = Counter(r.get(factor) for r in sheet if r['sample_id'] not in low)
        short = sorted(level for level in {r.get(factor) for r in sheet} - {None} if kept[level] < 2)
        if factor not in sheet[0]:
            reason = f'primary factor {factor!r} is not a column of the sample sheet'
        elif short:
            reason = f"fewer than 2 quantified samples in {factor} level(s): {', '.join(short)}"
        de = 'not_estimable' if reason else 'requested'
        if reason:
            objectives['pathogen'] = ['gene_expression']
    status = {'host': {'samples': [r['sample_id'] for r in sheet]},
              'pathogen': {'status': 'quantified', 'samples': [r['sample_id'] for r in sheet if r['sample_id'] not in low],
                           'not_quantified': low, 'min_assigned_reads': floor, 'de': de, 'reason': reason}}
    if not status['pathogen']['samples']:
        status['pathogen'].update(status='skipped', de='skipped', reason=f'no sample has at least {floor} pathogen reads')
    mates = ('_R1', '_R2') if full['samples']['seq_type'] == 'rnaseq_paired' else ('',)
    for role in ROLES:
        child = project / 'dual' / role
        if not status[role]['samples']:
            (child / 'config/config.yaml').unlink(missing_ok=True)
            continue
        rows = []
        for r in sheet:
            if r['sample_id'] not in status[role]['samples']:
                continue
            reads = [out / 'quant' / r['sample_id'] / 'dual' / f'{role}{mate}.fastq.gz' for mate in mates]
            for path in reads:
                if not path.is_file():
                    raise ValueError(f'Split reads not found: {path}')
            rows.append({'sample_id': r['sample_id'], **dict(zip(('fastq_url', 'fastq_url_r2'), map(str, reads))),
                         **{k: v for k, v in r.items() if k not in ('sample_id', 'fastq_url', 'fastq_url_r2')}})
        table_text = io.StringIO()
        writer = csv.DictWriter(table_text, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
        (child / 'config').mkdir(parents=True, exist_ok=True)
        _write(child / 'config/config.yaml', yaml.safe_dump(child_config(cfg, role, objectives[role]), sort_keys=False))
        _write(child / 'config/samples.tsv', table_text.getvalue())
        # The launcher reads tier sizes from the project's own copy.
        _write(child / 'config/thresholds.yaml', (project / 'config/thresholds.yaml').read_text())
    (out / 'pathogen').mkdir(parents=True, exist_ok=True)
    with open(out / 'pathogen/samples_not_quantified.tsv', 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t', lineterminator='\n')
        writer.writerow(['sample', 'pathogen_reads', 'reason'])
        writer.writerows([s, assigned[s], f'below min_assigned_reads {floor}'] for s in low)
    (out / 'pathogen/status.json').write_text(json.dumps(status['pathogen'], indent=2) + '\n')
    return status


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest='command', required=True)
    c = sub.add_parser('combine')
    for name in ('host', 'pathogen', 'host-gtf', 'pathogen-gtf', 'out-fa', 'out-gtf'):
        c.add_argument(f'--{name}', required=True)
    a = sub.add_parser('assign')
    a.add_argument('--index', required=True); a.add_argument('--reads', nargs='+', required=True)
    a.add_argument('--outdir', required=True); a.add_argument('--sample', required=True)
    a.add_argument('--multimap-limit', type=int, default=20); a.add_argument('--threads', type=int, default=4)
    for p in (c, a):
        p.add_argument('--host-label', required=True); p.add_argument('--pathogen-label', required=True)
    t = sub.add_parser('table')
    t.add_argument('--warn', type=float, required=True); t.add_argument('--out', required=True)
    t.add_argument('records', nargs='+')
    k = sub.add_parser('children')
    k.add_argument('--project', required=True); k.add_argument('--configfile', default='config/config.yaml')
    args = ap.parse_args()
    if args.command == 'combine':
        with tempfile.TemporaryDirectory(prefix='.combine-', dir=Path(args.out_fa).parent) as tmp:
            fa, gtf = Path(tmp) / 'combined.fa', Path(tmp) / 'combined.gtf'
            combine(args.host, args.pathogen, args.host_gtf, args.pathogen_gtf, fa, gtf,
                    args.host_label, args.pathogen_label)
            os.replace(fa, args.out_fa)
            os.replace(gtf, args.out_gtf)
    elif args.command == 'assign':
        align(args.index, args.reads, args.outdir, args.sample, args.host_label, args.pathogen_label,
              args.multimap_limit, args.threads)
    elif args.command == 'children':
        made = children(args.project, yaml.safe_load((Path(args.project) / args.configfile).read_text()))
        p = made['pathogen']
        print(f"dual-organism children: host {len(made['host']['samples'])} sample(s); pathogen {p['status']}, "
              f"{len(p['samples'])} sample(s), {len(p['not_quantified'])} not quantified, "
              f"differential expression {p['de']}" + (f" ({p['reason']})" if p['reason'] else ''))
    else:
        Path(args.out).write_text(table([json.loads(Path(r).read_text()) for r in args.records], args.warn))

#!/usr/bin/env python3
"""Host and pathogen in one library: combined reference, read labels and splitting."""
import gzip
import itertools
import re
from pathlib import Path

from prepare_reads import records

ROLES = ('host', 'pathogen')
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
        o.setdefault('objectives', ['gene_expression', 'differential_expression'] if role == 'host' else ['gene_expression'])
        d[role] = o
    if d['host']['label'] == d['pathogen']['label']:
        raise ValueError('dual_organism host and pathogen labels must differ')
    if not d['host'].get('gtf_url'):
        raise ValueError('dual_organism.host.gtf_url is required')
    d['pathogen'].setdefault('min_assigned_reads', 10000)
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

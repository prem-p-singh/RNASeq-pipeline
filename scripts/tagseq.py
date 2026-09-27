#!/usr/bin/env python3
"""Named QuantSeq kit profiles; unknown kits block the TAG-seq route."""
import csv
import gzip
from pathlib import Path
import yaml

KITS = yaml.safe_load((Path(__file__).resolve().parent.parent / 'config/tagseq_kits.yaml').read_text())['kits']


def load_kit(name):
    if name not in KITS:
        raise ValueError(f'Unknown TAG-seq kit {name!r}; choose one of {sorted(KITS)}')
    return dict(KITS[name])


def kit_policy(kit):
    return dict(adapter_r1=kit['adapter_r1'], trim_front=kit['trim_front'], umi_len=kit['umi_len'],
                umi_skip=kit['umi_skip'], poly_x=kit['poly_x'], minimum_length=20)


def merge_policy(preprocessing, kit):
    """The kit wins over defaults; a config that sets a kit key to another value is refused."""
    from preprocessing import DEFAULTS
    policy = kit_policy(kit)
    clashes = sorted(k for k, v in policy.items()
                     if k in preprocessing and preprocessing[k] not in (v, DEFAULTS[k]))
    if clashes:
        raise ValueError('The TAG-seq kit sets ' + ', '.join(f'preprocessing.{k}' for k in clashes)
                         + '; remove them from the config or drop analysis.kit')
    return preprocessing | policy


def check_spacer(fastq_gz, umi_len=6, spacer='TATA', reads=10000, minimum=0.8):
    """Declared FWD-UMI geometry must be visible in the reads before any trimming."""
    seen = hits = 0
    with gzip.open(fastq_gz, 'rt') as handle:
        for i, line in enumerate(handle):
            if i % 4 == 1:
                seen += 1
                hits += line[umi_len:umi_len + len(spacer)] == spacer
                if seen == reads:
                    break
    fraction = hits / seen if seen else 0.0
    if fraction < minimum:
        raise ValueError(f'{fraction:.1%} of reads carry {spacer} after a {umi_len} nt UMI; '
                         f'expected at least {minimum:.0%}. Check the declared kit.')
    return fraction


def check_strand(kit, libtypes):
    """Declared strands ({sample: libtype}) must agree with the kit; blank/unknown defer to it."""
    expected = {'forward': 'SF', 'reverse': 'SR'}[kit['strandedness']]
    bad = sorted(s for s, t in libtypes.items() if t in ('U', 'SF', 'SR', 'IU', 'ISF', 'ISR') and t != expected)
    if bad:
        raise ValueError(f'expected_libtype conflicts with the TAG-seq kit strand {expected} for: {", ".join(bad)}')


def umi_table(genes_tsv, umi_tsv, out):
    """umi_tools per-gene molecules in STAR gene order; unobserved annotated genes are 0."""
    with open(genes_tsv) as f:
        genes = [r[0] for r in list(csv.reader(f, delimiter='\t'))[1:]]
    with open(umi_tsv) as f:
        rows = list(csv.reader(f, delimiter='\t'))
    if not rows or rows[0] != ['gene', 'count']:
        raise ValueError(f'Unexpected umi_tools header {rows[:1]}; expected gene, count')
    counts = {g: int(n) for g, n in rows[1:]}
    outside = sorted(set(counts) - set(genes))
    if outside:
        raise ValueError(f'{len(outside)} UMI genes outside the annotation, e.g. {outside[:3]}')
    with open(out, 'w') as f:
        f.write('gene_id\tcount\n' + ''.join(f'{g}\t{counts.get(g, 0)}\n' for g in genes))

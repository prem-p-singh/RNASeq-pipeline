#!/usr/bin/env python3
"""Named QuantSeq kit profiles; unknown kits block the TAG-seq route."""
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

#!/usr/bin/env python3
"""Named QuantSeq kit profiles; unknown kits block the TAG-seq route."""
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

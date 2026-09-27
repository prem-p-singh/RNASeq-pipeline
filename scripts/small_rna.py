#!/usr/bin/env python3
"""Small-RNA kit profiles, mature-miRNA references and read counting."""
from pathlib import Path
import yaml

TABLE = yaml.safe_load((Path(__file__).resolve().parent.parent / 'config/small_rna_kits.yaml').read_text())


def load_kit(name):
    if name in TABLE.get('blocked', {}):
        raise ValueError(f'{name}: {TABLE["blocked"][name]}')
    if name not in TABLE['kits']:
        raise ValueError(f'Unknown small-RNA kit {name!r}; choose one of {sorted(TABLE["kits"])}')
    return dict(TABLE['kits'][name])


def kit_policy(kit):
    return dict(adapter_r1=kit['adapter_r1'], minimum_length=kit['minimum_length'],
                length_limit=kit['length_limit'], poly_g='on')

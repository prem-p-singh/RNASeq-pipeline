#!/usr/bin/env python3
"""Small-RNA kit adapters and size selection come from one cited table."""
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from small_rna import load_kit, kit_policy
from preprocessing import policy_args

truseq = load_kit('TruSeq Small RNA')
assert truseq['adapter_r1'] == 'TGGAATTCTCGGGTGCCAAGG'
p, args = policy_args(kit_policy(truseq), paired=False)
assert args[args.index('--length_required') + 1] == '18' and args[args.index('--length_limit') + 1] == '30'
assert '--trim_poly_g' in args
try:
    load_kit('QIAseq miRNA')
except ValueError as exc:
    assert 'UMI' in str(exc)
else:
    raise AssertionError('UMI small-RNA kit accepted')
print('small RNA kit table passed')

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

from preprocessing import loss_summary
fake = {'summary': {'before_filtering': {'total_reads': 100}, 'after_filtering': {'total_reads': 70}},
        'filtering_result': {'passed_filter_reads': 70, 'too_short_reads': 10, 'too_long_reads': 20,
                             'low_quality_reads': 0, 'too_many_N_reads': 0}}
s = loss_summary(fake)
assert s['filtering_result']['too_long_reads'] == 20
assert sum(v for k, v in s['filtering_result'].items() if k != 'passed_filter_reads') == 30
print('small RNA loss summary passed')

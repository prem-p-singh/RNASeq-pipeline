#!/usr/bin/env python3
"""Kit facts come from one cited table; unknown kits block."""
import gzip
import sys
import tempfile
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from tagseq import load_kit, KITS, kit_policy, check_spacer
from preprocessing import policy_args

fwd, rev, umi = (load_kit(k) for k in ('QuantSeq FWD', 'QuantSeq REV', 'QuantSeq FWD-UMI'))
assert fwd['strandedness'] == 'forward' and fwd['umi_len'] == 0
assert rev['strandedness'] == 'reverse' and rev['umi_len'] == 0
assert umi['strandedness'] == 'forward' and umi['umi_len'] == 6 and umi['umi_skip'] == 4
assert all(k['status'] == 'experimental' and k['sources'] for k in KITS.values())
try:
    load_kit('QuantSeq generic')
except ValueError as exc:
    assert 'Unknown TAG-seq kit' in str(exc)
else:
    raise AssertionError('unknown kit accepted')
print('TAG-seq kit table passed')

p, args = policy_args(kit_policy(load_kit('QuantSeq FWD-UMI')), paired=False)
assert args[args.index('--umi_len') + 1] == '6' and args[args.index('--umi_skip') + 1] == '4'
assert '--umi' in args and args[args.index('--umi_loc') + 1] == 'read1'
assert args[args.index('--trim_front1') + 1] == '12' and '--trim_poly_x' in args
p, args = policy_args(kit_policy(load_kit('QuantSeq REV')), paired=False)
assert '--umi' not in args and '--trim_front1' not in args
try:
    policy_args(kit_policy(load_kit('QuantSeq FWD')), paired=True)
except ValueError as exc:
    assert 'single-end' in str(exc)
else:
    raise AssertionError('paired TAG-seq accepted')
print('TAG-seq preprocessing arguments passed')

with tempfile.TemporaryDirectory() as t:
    good, bad = Path(t, 'good.fq.gz'), Path(t, 'bad.fq.gz')
    with gzip.open(good, 'wt') as g, gzip.open(bad, 'wt') as b:
        for i in range(200):
            g.write(f'@r{i}\nACGTAC' + 'TATA' + 'G' * 40 + '\n+\n' + 'I' * 50 + '\n')
            b.write(f'@r{i}\nACGTAC' + 'CCCC' + 'G' * 40 + '\n+\n' + 'I' * 50 + '\n')
    assert check_spacer(good) == 1.0
    try:
        check_spacer(bad)
    except ValueError as exc:
        assert '0.0%' in str(exc) and 'TATA' in str(exc)
    else:
        raise AssertionError('missing spacer accepted')
print('UMI spacer check passed')

from tagseq import merge_policy
from preprocessing import DEFAULTS
merged = merge_policy(dict(DEFAULTS, qualified_quality_phred=20), umi)
assert merged['umi_len'] == 6 and merged['trim_front'] == 12 and merged['qualified_quality_phred'] == 20
assert merge_policy(dict(DEFAULTS, trim_front=12), fwd)['trim_front'] == 12  # same value is not a conflict
try:
    merge_policy(dict(DEFAULTS, trim_front=5), fwd)
except ValueError as exc:
    assert 'preprocessing.trim_front' in str(exc) and 'kit' in str(exc), exc
else:
    raise AssertionError('conflicting kit key accepted')
print('TAG-seq kit policy merge passed')

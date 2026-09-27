#!/usr/bin/env python3
"""Kit facts come from one cited table; unknown kits block."""
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from tagseq import load_kit, KITS

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

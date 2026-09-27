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
fake = {'summary': {'before_filtering': {'total_reads': 100}, 'after_filtering': {'total_reads': 70, 'read1_mean_length': 22}},
        'filtering_result': {'passed_filter_reads': 70, 'too_short_reads': 10, 'too_long_reads': 20,
                             'low_quality_reads': 0, 'too_many_N_reads': 0}}
s = loss_summary(fake)
assert s['filtering_result']['too_long_reads'] == 20
assert sum(v for k, v in s['filtering_result'].items() if k != 'passed_filter_reads') == 30
assert s['size_selection'] == {'too_short_reads': 10, 'too_long_reads': 20, 'mean_length_after': 22}
print('small RNA loss summary passed')

import json
import tempfile
from small_rna import mature_reference, count_alignments, count_metrics
with tempfile.TemporaryDirectory() as t:
    t = Path(t)
    (t / 'mature.fa').write_text('>hsa-miR-1 MIMAT1\nUGGAAUGUAAAGAAGUAUGUAU\n>mmu-miR-1 MIMAT2\nUGGAAUGUAAAGAAGUAUGUAU\n')
    assert mature_reference(t / 'mature.fa', 'hsa', t / 'ref.fa') == 1
    assert (t / 'ref.fa').read_text() == '>hsa-miR-1\nTGGAATGTAAAGAAGTATGTAT\n'
    try:
        mature_reference(t / 'mature.fa', 'vvx', t / 'none.fa')
    except ValueError as exc:
        assert 'vvx' in str(exc)
    else:
        raise AssertionError('empty species reference accepted')
    # SAM from bowtie -a --best --strata: r1 unique, r2 two equal hits, r3 unaligned.
    # Bowtie 1 SAM has no NH tag: multi-mapping is detected from repeated read names.
    (t / 'a.sam').write_text('@HD\tVN:1.0\n'
        'r1\t0\thsa-miR-1\t1\t255\t22M\t*\t0\t0\tA\tI\n'
        'r2\t0\thsa-miR-2a\t1\t255\t22M\t*\t0\t0\tA\tI\n'
        'r2\t0\thsa-miR-2b\t1\t255\t22M\t*\t0\t0\tA\tI\n'
        'r3\t4\t*\t0\t0\t*\t*\t0\t0\tA\tI\n')
    totals = count_alignments(t / 'a.sam', t / 'counts.tsv')
    assert totals == {'unique': 1, 'multi_mapped': 1, 'unaligned': 1}, totals
    assert (t / 'counts.tsv').read_text() == 'gene_id\tcount\nhsa-miR-1\t1\n'
    # Aggregation reads STAR-shaped metrics: mapping_rate = (unique + multi) / reads.
    m = count_metrics(t / 'a.sam', t / 'c2.tsv', t / 'metrics.json', 'S1')
    assert json.loads((t / 'metrics.json').read_text()) == m
    assert (m['sample'], m['num_reads_processed'], m['num_assigned']) == ('S1', 3, 1)
    assert abs(m['mapping_rate'] - 2 / 3) < 1e-12
print('animal miRNA reference and counting passed')

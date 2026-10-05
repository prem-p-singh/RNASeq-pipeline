#!/usr/bin/env python3
"""Dual-organism labels, combined reference and stream splitting."""
import gzip, sys, tempfile
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from dual_organism import validate, combine, label, assign

ok = {'dual_organism': {'enabled': True,
      'host': {'label': 'grape', 'genome_fasta_url': 'file:///h.fa', 'gtf_url': 'file:///h.gtf'},
      'pathogen': {'label': 'GRBV', 'genome_fasta_url': 'file:///p.fa'}},
      'samples': {'seq_type': 'rnaseq_single'}}
v = validate(ok)
assert v['multimap_limit'] == 20 and v['ambiguous_warn_fraction'] == 0.05
assert v['pathogen']['min_assigned_reads'] == 10000 and v['pathogen']['objectives'] == ['gene_expression']
for path, value, word in ((('host', 'label'), 'a|b', 'label'), (('pathogen', 'label'), 'grape', 'differ'),
                          (('pathogen', 'genome_fasta_url'), None, 'genome_fasta_url')):
    bad = {**ok, 'dual_organism': {**ok['dual_organism'], path[0]: {**ok['dual_organism'][path[0]], path[1]: value}}}
    try: validate(bad)
    except ValueError as exc: assert word in str(exc), exc
    else: raise AssertionError(path)
try: validate({**ok, 'samples': {'seq_type': 'tagseq'}})
except ValueError as exc: assert 'tagseq' in str(exc)
else: raise AssertionError('tagseq accepted')

assert label({'host'}, False) == 'host' and label({'pathogen'}, False) == 'pathogen'
assert label({'host', 'pathogen'}, False) == 'ambiguous' and label(set(), False) == 'unassigned'
assert label(set(), True) == 'unassigned'

with tempfile.TemporaryDirectory() as t:
    t = Path(t)
    (t / 'h.fa').write_text('>chr1 host\nACGT\n>chr2\nAAAA\n'); (t / 'p.fa').write_text('>chr1\nGGGG\n')
    (t / 'h.gtf').write_text('chr1\tx\texon\t1\t4\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n')
    n = combine(t / 'h.fa', t / 'p.fa', t / 'h.gtf', None, t / 'c.fa', t / 'c.gtf', 'grape', 'GRBV')
    assert n == {'host_contigs': 2, 'pathogen_contigs': 1}
    assert (t / 'c.fa').read_text() == '>grape|chr1\nACGT\n>grape|chr2\nAAAA\n>GRBV|chr1\nGGGG\n'
    assert (t / 'c.gtf').read_text().startswith('grape|chr1\tx\texon')
    (t / 'dup.fa').write_text('>chr1\nA\n>chr1\nC\n')
    try: combine(t / 'dup.fa', t / 'p.fa', t / 'h.gtf', None, t / 'd.fa', t / 'd.gtf', 'grape', 'GRBV')
    except ValueError as exc: assert 'Duplicate' in str(exc)
    else: raise AssertionError('duplicate contig accepted')

    def fq(path, names):
        with gzip.open(path, 'wt') as f:
            for n in names: f.write(f'@{n} comment\nACGT\n+\nIIII\n')
    names = ['r1/1', 'r2/1', 'r3/1', 'r4/1', 'r5/1']
    fq(t / 'in.fq.gz', names)
    def sam(read, flag, ref, tag=''):
        return f'{read}\t{flag}\t{ref}\t1\t255\t4M\t*\t0\t0\tACGT\tIIII\tNH:i:1{tag}\n'
    lines = ['@HD\tVN:1.4\n', sam('r1', 0, 'grape|chr1'), sam('r2', 0, 'GRBV|chr1'),
             sam('r3', 0, 'grape|chr1'), sam('r3', 256, 'GRBV|chr1'),
             sam('r4', 4, '*', '\tuT:A:1'), sam('r5', 4, '*', '\tuT:A:3')]
    outs = {'host': [t / 'host.fq.gz'], 'pathogen': [t / 'path.fq.gz']}
    c = assign(lines, [t / 'in.fq.gz'], outs, 'grape', 'GRBV')
    assert c == {'reads': 5, 'host': 1, 'pathogen': 1, 'ambiguous': 1, 'unassigned': 2, 'too_many_loci': 1}, c
    assert gzip.open(t / 'host.fq.gz', 'rt').read().startswith('@r1/1 comment\n')
    assert gzip.open(t / 'path.fq.gz', 'rt').read().count('@') == 1
    # Out-of-order alignments must fail loudly, never mislabel.
    try: assign(['@HD\n', sam('r2', 0, 'GRBV|chr1'), sam('r1', 0, 'grape|chr1')], [t / 'in.fq.gz'], outs, 'grape', 'GRBV')
    except ValueError as exc: assert 'order' in str(exc)
    else: raise AssertionError('out-of-order SAM accepted')
    # A missing read is a conservation failure.
    try: assign(['@HD\n', sam('r1', 0, 'grape|chr1')], [t / 'in.fq.gz'], outs, 'grape', 'GRBV')
    except ValueError as exc: assert 'conservation' in str(exc) or 'order' in str(exc)
    else: raise AssertionError('missing reads accepted')
print('dual-organism labels, reference and splitting passed')

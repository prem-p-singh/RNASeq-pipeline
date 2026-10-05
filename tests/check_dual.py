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
# Null fields from the config template mean unset.
nulls = {**ok, 'dual_organism': {**ok['dual_organism'], 'pathogen': {**ok['dual_organism']['pathogen'], 'objectives': None, 'min_assigned_reads': None, 'gtf_url': None}}}
assert validate(nulls)['pathogen'] == v['pathogen'] | {'gtf_url': None}
for path, value, word in ((('host', 'label'), 'a|b', 'label'), (('pathogen', 'label'), 'grape', 'differ'),
                          (('pathogen', 'genome_fasta_url'), None, 'genome_fasta_url')):
    bad = {**ok, 'dual_organism': {**ok['dual_organism'], path[0]: {**ok['dual_organism'][path[0]], path[1]: value}}}
    try: validate(bad)
    except ValueError as exc: assert word in str(exc), exc
    else: raise AssertionError(path)
try: validate({**ok, 'samples': {'seq_type': 'tagseq'}})
except ValueError as exc: assert 'tagseq' in str(exc)
else: raise AssertionError('tagseq accepted')

assert label({'host'}) == 'host' and label({'pathogen'}) == 'pathogen'
assert label({'host', 'pathogen'}) == 'ambiguous' and label(set()) == 'unassigned'

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
    # Hardening: paired input, mate mismatch, malformed references, no partial output.
    fq(t / 'a1.fq.gz', ['r1/1', 'r2/1']); fq(t / 'a2.fq.gz', ['r1/2', 'r2/2']); fq(t / 'z2.fq.gz', ['zz/2', 'r2/2'])
    po = {'host': [t / 'h1.fq.gz', t / 'h2.fq.gz'], 'pathogen': [t / 'p1.fq.gz', t / 'p2.fq.gz']}
    pair = [t / 'a1.fq.gz', t / 'a2.fq.gz']
    c = assign(['@HD\n', sam('r1', 0, 'grape|chr1'), sam('r1', 0, 'grape|chr2'), sam('r2', 0, 'grape|chr1'), sam('r2', 4, '*')],
               pair, po, 'grape', 'GRBV')
    assert c['host'] == 2, c
    assert all(gzip.open(p, 'rt').read().count('@') == 2 for p in po['host'])
    for p in po['host'] + po['pathogen']: p.unlink(missing_ok=True)
    c = assign(['@HD\n', sam('r1', 0, 'grape|chr1'), sam('r1', 0, 'GRBV|chr1'), sam('r2', 0, 'GRBV|chr1'), sam('r2', 4, '*')],
               pair, po, 'grape', 'GRBV')
    assert c['ambiguous'] == 1 and c['pathogen'] == 1 and c['host'] == 0, c
    assert gzip.open(po['pathogen'][0], 'rt').read().count('@') == 1
    assert gzip.open(po['host'][0], 'rt').read() == ''
    def none_left(ps):
        return not any(p.exists() or Path(f'{p}.tmp').exists() for ps_ in ps.values() for p in ps_)
    for p in po['host'] + po['pathogen']: p.unlink()
    for bad, word in (([sam('r1', 0, 'grape|chr1'), sam('r2', 0, 'GRBV|chr1')], 'order'),
                      ([sam('r1', 0, 'grape|chr1'), sam('r2', 0, 'other|chr1')], 'neither')):
        fqs = [t / 'a1.fq.gz', t / 'z2.fq.gz'] if word == 'order' else pair
        try: assign(['@HD\n'] + bad, fqs, po, 'grape', 'GRBV')
        except ValueError as exc: assert word in str(exc), exc
        else: raise AssertionError(word)
        assert none_left(po)
    (t / 'pd.fa').write_text('>a\nA\n>a\nC\n')
    try: combine(t / 'h.fa', t / 'pd.fa', t / 'h.gtf', None, t / 'e.fa', t / 'e.gtf', 'grape', 'GRBV')
    except ValueError as exc: assert 'Duplicate' in str(exc) and 'pathogen' in str(exc)
    else: raise AssertionError('pathogen duplicate accepted')
    (t / 'nn.fa').write_text('>a\nACGT'); (t / 'nn.gtf').write_text('a\tx\texon\t1\t4\t.\t+\t.\tgene_id "g";')
    (t / 'pn.fa').write_text('>b\nGG\n'); (t / 'pn.gtf').write_text('b\tx\texon\t1\t2\t.\t+\t.\tgene_id "p";\n')
    combine(t / 'nn.fa', t / 'pn.fa', t / 'nn.gtf', t / 'pn.gtf', t / 'n.fa', t / 'n.gtf', 'grape', 'GRBV')
    assert (t / 'n.fa').read_text() == '>grape|a\nACGT\n>GRBV|b\nGG\n'
    assert [l[:7] for l in (t / 'n.gtf').read_text().splitlines()] == ['grape|a', 'GRBV|b\t']
    (t / 'eh.fa').write_text('>\nA\n')
    try: combine(t / 'eh.fa', t / 'p.fa', t / 'h.gtf', None, t / 'f.fa', t / 'f.gtf', 'grape', 'GRBV')
    except ValueError as exc: assert 'eh.fa' in str(exc)
    else: raise AssertionError('empty header accepted')
print('dual-organism labels, reference and splitting passed')

# Pass-through preprocessing for already-processed reads (Task 2)
sys.path.insert(0, str(ROOT / 'scripts'))
from preprocessing import policy_args
p, args = policy_args({'passthrough': True}, paired=False)
assert args == ['-A', '-Q', '-L', '-G'], args
p, args = policy_args({'passthrough': True}, paired=True)
assert args == ['-A', '-Q', '-L', '-G'], args
assert '-A' not in policy_args({}, paired=False)[1]
try: policy_args({'passthrough': True, 'adapter_r1': 'ACGT'}, paired=False)
except ValueError as exc: assert 'passthrough' in str(exc)
else: raise AssertionError('passthrough with trimming settings accepted')
print('pass-through preprocessing passed')

# Cohort assignment table (Task 3)
from dual_organism import table
HEADER = ('sample\treads\thost\tpathogen\tambiguous\tunassigned\ttoo_many_loci\thost_fraction\tpathogen_fraction\t'
          'ambiguous_fraction\tunassigned_fraction\twarnings')
rows = table([{'sample': 's1', 'counts': dict(reads=100, host=80, pathogen=10, ambiguous=6, unassigned=4, too_many_loci=1)},
              {'sample': 's2', 'counts': dict(reads=0, host=0, pathogen=0, ambiguous=0, unassigned=0, too_many_loci=0)}], 0.05).splitlines()
assert rows[0] == HEADER and len(rows) == 3, rows
s1, s2 = (dict(zip(HEADER.split('\t'), r.split('\t'))) for r in rows[1:])
assert s1['host_fraction'] == '0.800000' and 'ambiguous fraction' in s1['warnings'], s1
assert s2['warnings'] == '' and s2['host_fraction'] == '0.000000', s2
print('dual-organism assignment table passed')

# STAR wiring, with a stand-in STAR: record written last, STAR failure reported, nothing partial left.
import json, os, stat
import dual_organism
from dual_organism import align
dual_organism.STAR_EXIT_GRACE = 1
with tempfile.TemporaryDirectory() as t:
    t = Path(t)
    (t / 'bin').mkdir(); (t / 'idx').mkdir()
    (t / 'bin/STAR').write_text('#!/bin/sh\n[ "$1" = --version ] && { echo 2.7.11b; exit 0; }\ncat "$FAKE_SAM"\n'
                                '[ "$FAKE_RC" = hang ] && exec sleep 60\n[ "$FAKE_RC" = kill ] && kill -9 $$\nexit "$FAKE_RC"\n')
    (t / 'bin/STAR').chmod(stat.S_IRWXU)
    os.environ['PATH'] = f"{t / 'bin'}{os.pathsep}{os.environ['PATH']}"
    (t / 'idx/reference.lock.json').write_text(json.dumps({'request': {'STAR': '2.7.11b'}}))
    with gzip.open(t / 'in.fq.gz', 'wt') as f:
        for n in ('r1', 'r2'): f.write(f'@{n}\nACGT\n+\nIIII\n')
    (t / 'ok.sam').write_text('@HD\tVN:1.4\n' + sam('r1', 0, 'grape|chr1') + sam('r2', 0, 'GRBV|chr1'))
    (t / 'short.sam').write_text(sam('r1', 0, 'grape|chr1'))
    os.environ.update(FAKE_SAM=str(t / 'ok.sam'), FAKE_RC='0')
    r = align(t / 'idx', [t / 'in.fq.gz'], t / 'out', 's1', 'grape', 'GRBV', limit=7, threads=2)
    assert r['counts']['host'] == r['counts']['pathogen'] == 1 and r['sample'] == 's1', r
    assert r['command'][r['command'].index('--outFilterMultimapNmax') + 1] == '7'
    assert json.loads((t / 'out/assignment.json').read_text()) == r
    assert sorted(p.name for p in (t / 'out').iterdir()) == ['assignment.json', 'host.fastq.gz', 'pathogen.fastq.gz']
    for sam_file, rc, status in (('ok.sam', '3', 3), ('short.sam', '3', 3), ('ok.sam', 'kill', -9), ('short.sam', 'kill', -9)):
        os.environ.update(FAKE_SAM=str(t / sam_file), FAKE_RC=rc)
        try: align(t / 'idx', [t / 'in.fq.gz'], t / 'out', 's1', 'grape', 'GRBV')
        except RuntimeError as exc: assert f'STAR exited with status {status}' in str(exc), exc
        else: raise AssertionError('STAR failure accepted')
        assert not list((t / 'out').iterdir()), list((t / 'out').iterdir())
    # A STAR that is still running when the stream is wrong is stopped, and the stream error is reported.
    (t / 'swap.sam').write_text(sam('r2', 0, 'GRBV|chr1') + sam('r1', 0, 'grape|chr1'))
    os.environ.update(FAKE_SAM=str(t / 'swap.sam'), FAKE_RC='hang')
    try: align(t / 'idx', [t / 'in.fq.gz'], t / 'out', 's1', 'grape', 'GRBV')
    except ValueError as exc: assert 'order' in str(exc), exc
    else: raise AssertionError('out-of-order STAR output accepted')
    assert not list((t / 'out').iterdir()), list((t / 'out').iterdir())
    (t / 'idx/reference.lock.json').write_text(json.dumps({'request': {'STAR': '2.7.10a'}}))
    try: align(t / 'idx', [t / 'in.fq.gz'], t / 'out', 's1', 'grape', 'GRBV')
    except ValueError as exc: assert 'version' in str(exc)
    else: raise AssertionError('STAR version mismatch accepted')
print('dual-organism STAR wiring passed')

#!/usr/bin/env python3
"""Dual-organism labels, combined reference and stream splitting."""
import gzip, sys, tempfile
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from dual_organism import validate, combine, label, assign

ok = {'dual_organism': {'enabled': True,
      'host': {'label': 'grape', 'genome_fasta_url': 'file:///h.fa', 'gtf_url': 'file:///h.gtf'},
      'pathogen': {'label': 'GRBV', 'genome_fasta_url': 'file:///p.fa', 'gtf_url': 'file:///p.gtf'}},
      'samples': {'seq_type': 'rnaseq_single'}, 'analysis': {'quantifier': 'star'}}
v = validate(ok)
# Salmon quantifies against transcripts: each organism needs its own.
for quantifier in ('salmon', 'auto'):
    try: validate({**ok, 'analysis': {'quantifier': quantifier}})
    except ValueError as exc: assert 'dual_organism.host.transcriptome_fasta_url' in str(exc), exc
    else: raise AssertionError('Salmon route without transcripts accepted')
assert v['multimap_limit'] == 20 and v['ambiguous_warn_fraction'] == 0.05
assert v['pathogen']['min_assigned_reads'] == 10000 and v['pathogen']['objectives'] == ['gene_expression']
# Null fields from the config template mean unset.
nulls = {**ok, 'dual_organism': {**ok['dual_organism'], 'pathogen': {**ok['dual_organism']['pathogen'], 'objectives': None, 'min_assigned_reads': None}}}
assert validate(nulls)['pathogen'] == v['pathogen']
for path, value, word in ((('host', 'label'), 'a|b', 'label'), (('pathogen', 'label'), 'grape', 'differ'),
                          (('pathogen', 'genome_fasta_url'), None, 'genome_fasta_url'),
                          # Each organism is counted against its own annotation afterwards.
                          (('host', 'gtf_url'), None, 'dual_organism.host.gtf_url'),
                          (('pathogen', 'gtf_url'), None, 'dual_organism.pathogen.gtf_url')):
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

# Child projects (Task 4): one ordinary single-organism project per organism, built from the assignment table.
import copy, csv, subprocess, yaml
import config_resolve
from recommend import recommend
from dual_organism import children, child_config
PARENT = {'project': {'name': 'dual', 'output_dir': 'results/'},
          'analysis': {'quantifier': 'star'},
          'screening': {'enabled': True, 'references': [{'name': 'Host', 'role': 'expected', 'fasta': 'genome'}]},
          'reference': {'star_overhang': 49, 'decoys': 'config/decoys.txt', 'annotation_tsv': {'path': 'reference/annotation_info.tsv'}},
          'samples': {'sheet': 'config/samples.tsv', 'seq_type': 'rnaseq_paired', 'expected_libtype': 'ISR'},
          'model': {'fixed_effects': '~ treatment', 'random_effects': None, 'primary_factor': 'treatment'},
          'contrasts': [{'id': 'treatment', 'type': 'pairwise', 'factor': 'treatment'}],
          'thresholds': {'sample_qc': {'duplication_rate_warn': 0.5}},
          'dual_organism': {'enabled': True,
                            'host': {'label': 'grape', 'genome_fasta_url': 'file:///h.fa', 'gtf_url': 'file:///h.gtf',
                                     'transcriptome_fasta_url': 'file:///h.rna.fa'},
                            'pathogen': {'label': 'GRBV', 'genome_fasta_url': 'file:///p.fa', 'gtf_url': 'file:///p.gtf',
                                         'objectives': ['gene_expression', 'differential_expression']}}}
MOCK, INFECTED = ['m1', 'm2', 'm3'], ['i1', 'i2', 'i3']
SNAPSHOT = copy.deepcopy(PARENT)


def parent(t, pathogen_reads):
    """A finished parent run: sheet, assignment table and split files."""
    for d in ('config', 'results', 'reference'):
        (t / d).mkdir()
    (t / 'reference/annotation_info.tsv').write_text('Gene.stable.ID\n')
    (t / 'config/decoys.txt').write_text('chr1\n')
    (t / 'config/thresholds.yaml').write_text((ROOT / 'config/thresholds.yaml').read_text())
    (t / 'config/samples.tsv').write_text('# cohort\nsample_id\tfastq_url\tfastq_url_r2\ttreatment\n' + ''.join(
        f'{s}\thttps://x/{s}_1.fq.gz\thttps://x/{s}_2.fq.gz\t{"mock" if s in MOCK else "infected"}\n' for s in MOCK + INFECTED))
    (t / 'results/assignment.tsv').write_text(table(
        [{'sample': s, 'counts': dict(reads=10 ** 6, host=10 ** 6 - 100 - pathogen_reads[s], pathogen=pathogen_reads[s],
                                      ambiguous=60, unassigned=40, too_many_loci=0)} for s in MOCK + INFECTED], 0.05))
    for s in MOCK + INFECTED:
        (t / 'results/quant' / s / 'dual').mkdir(parents=True)
        for name in ('host_R1', 'host_R2', 'pathogen_R1', 'pathogen_R2'):
            (t / 'results/quant' / s / 'dual' / f'{name}.fastq.gz').touch()
    return t / 'results/assignment.tsv'


def child(t, role):
    cfg = yaml.safe_load((t / 'dual' / role / 'config/config.yaml').read_text())
    resolved, _, issues = config_resolve.resolve(ROOT, [('project', cfg)])
    assert not config_resolve.blocking(issues, ROOT), issues
    assert recommend(resolved)['status'] == 'selected', recommend(resolved)['issues']
    assert recommend(resolved)['route'] == 'star_counts' and not resolved['dual_organism']['enabled']
    with open(t / 'dual' / role / 'config/samples.tsv', newline='') as f:
        return cfg, resolved, list(csv.DictReader(f, delimiter='\t'))


with tempfile.TemporaryDirectory() as t:
    t = Path(t).resolve()
    tsv = parent(t, {s: 12 if s in MOCK else 50000 for s in MOCK + INFECTED})
    status = children(t, PARENT, tsv)
    assert status['host'] == {'samples': MOCK + INFECTED}, status
    assert status['pathogen']['samples'] == INFECTED and status['pathogen']['not_quantified'] == MOCK, status
    assert status['pathogen']['status'] == 'quantified' and status['pathogen']['de'] == 'not_estimable', status
    assert 'mock' in status['pathogen']['reason'] and 'treatment' in status['pathogen']['reason'], status
    assert json.loads((t / 'results/pathogen/status.json').read_text()) == status['pathogen']
    with open(t / 'results/pathogen/samples_not_quantified.tsv', newline='') as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    assert [r['sample'] for r in rows] == MOCK and all(r['pathogen_reads'] == '12' for r in rows), rows
    assert all('below min_assigned_reads' in r['reason'] for r in rows), rows

    cfg, resolved, sheet = child(t, 'host')
    assert [r['sample_id'] for r in sheet] == MOCK + INFECTED
    assert sheet[0]['fastq_url'] == str(t / 'results/quant/m1/dual/host_R1.fastq.gz'), sheet[0]
    assert sheet[0]['fastq_url_r2'] == str(t / 'results/quant/m1/dual/host_R2.fastq.gz') and sheet[0]['treatment'] == 'mock'
    assert cfg['preprocessing'] == {'passthrough': True} and cfg['project']['output_dir'] == '../../results/host/'
    assert cfg['screening'] == {'enabled': False} and cfg['samples']['metadata_dir'] is None
    assert cfg['reference'] | {'annotation_tsv': 0} == {
        'star_overhang': 49, 'annotation_tsv': 0, 'genome_fasta_url': 'file:///h.fa', 'gtf_url': 'file:///h.gtf',
        'transcriptome_fasta_url': 'file:///h.rna.fa', 'accession': 'grape', 'decoys': '../../config/decoys.txt'}, cfg['reference']
    # Host files named relative to the parent project are reached from dual/host.
    assert cfg['reference']['annotation_tsv'] == {'path': '../../reference/annotation_info.tsv'}
    assert (t / 'dual/host' / cfg['reference']['decoys']).is_file() and (t / 'dual/host' / cfg['reference']['annotation_tsv']['path']).is_file()
    # The parent check and the written project come from one function.
    assert cfg == child_config(PARENT, 'host')
    absolute = {**PARENT, 'project': {'output_dir': '/data/out'}, 'reference': {'decoys': 'genome', 'annotation_tsv': {'path': '/data/info.tsv'}}}
    made = child_config(absolute, 'host')
    assert (made['project']['output_dir'], made['reference']['decoys'], made['reference']['annotation_tsv']['path']) == ('/data/out/host/', 'genome', '/data/info.tsv')
    assert child_config(absolute, 'pathogen')['reference']['decoys'] == 'genome'
    unchanged = [(p, p.stat().st_mtime_ns) for p in (t / 'dual/host/config').iterdir()]
    assert cfg['analysis'] == {'quantifier': 'star', 'objectives': ['gene_expression', 'differential_expression']}
    assert resolved['thresholds']['sample_qc']['mapping_rate_min'] == 0.60 and resolved['downstream']['run_go'] is True
    assert (t / 'dual/host/config/thresholds.yaml').read_text() == (t / 'config/thresholds.yaml').read_text()

    cfg, resolved, sheet = child(t, 'pathogen')
    assert [r['sample_id'] for r in sheet] == INFECTED
    assert sheet[0]['fastq_url'] == str(t / 'results/quant/i1/dual/pathogen_R1.fastq.gz'), sheet[0]
    assert cfg['analysis']['objectives'] == ['gene_expression'] and cfg['project']['output_dir'] == '../../results/pathogen/'
    assert cfg['downstream'] == {'run_go': False, 'run_kegg': False, 'run_wgcna': False}
    assert cfg['reference']['accession'] == 'GRBV' and cfg['reference']['transcriptome_fasta_url'] is None
    # The decoy list and the gene-information table describe the host.
    assert cfg['reference']['decoys'] is None and cfg['reference']['annotation_tsv'] == {'path': None}
    # Pathogen reads are pre-assigned: inclusion is min_assigned_reads, not the mapping-rate or depth floor.
    assert cfg['thresholds']['sample_qc'] == {'duplication_rate_warn': 0.5, 'mapping_rate_min': 0, 'min_reads_on_genes_rnaseq': 0}
    assert PARENT == SNAPSHOT, 'children changed the parent config'

    # Enough quantified samples in every level: differential expression stays requested.
    (t / 'results/assignment.tsv').write_text((t / 'results/assignment.tsv').read_text().replace('\t12\t', '\t20000\t'))
    status = children(t, PARENT, tsv)
    assert status['pathogen']['de'] == 'requested' and status['pathogen']['not_quantified'] == [], status
    # A relaunch leaves an unchanged child project untouched.
    assert unchanged == [(p, p.stat().st_mtime_ns) for p, _ in unchanged]
    assert child(t, 'pathogen')[0]['analysis']['objectives'] == ['gene_expression', 'differential_expression']
    nofactor = children(t, {**PARENT, 'model': {**PARENT['model'], 'primary_factor': 'dose'}}, tsv)['pathogen']
    assert nofactor['de'] == 'not_estimable' and 'dose' in nofactor['reason'], nofactor
    counts_only = {**PARENT, 'dual_organism': {**PARENT['dual_organism'], 'pathogen': {**PARENT['dual_organism']['pathogen'], 'objectives': None}}}
    assert children(t, counts_only, tsv)['pathogen']['de'] == 'not_requested'

    # A sample the assignment table does not cover, or a missing split file, stops the run.
    (t / 'results/quant/i3/dual/host_R2.fastq.gz').unlink()
    try: children(t, PARENT, tsv)
    except ValueError as exc: assert 'host_R2.fastq.gz' in str(exc), exc
    else: raise AssertionError('missing split file accepted')
    (t / 'results/assignment.tsv').write_text(''.join((t / 'results/assignment.tsv').read_text().splitlines(True)[:-1]))
    try: children(t, PARENT, tsv)
    except ValueError as exc: assert 'i3' in str(exc), exc
    else: raise AssertionError('sample missing from the assignment table accepted')

with tempfile.TemporaryDirectory() as t:
    t = Path(t).resolve()
    tsv = parent(t, {s: 12 for s in MOCK + INFECTED})
    status = children(t, PARENT, tsv)
    assert status['pathogen']['status'] == 'skipped' and status['pathogen']['samples'] == [], status
    assert status['pathogen']['de'] == 'skipped' and '10000' in status['pathogen']['reason'], status
    assert status['pathogen']['not_quantified'] == MOCK + INFECTED and not (t / 'dual/pathogen').exists()
    assert json.loads((t / 'results/pathogen/status.json').read_text())['status'] == 'skipped'
    assert len(child(t, 'host')[2]) == 6
    # A pathogen project left by an earlier run is not launched again once no sample qualifies.
    low = (t / 'results/assignment.tsv').read_text()
    (t / 'results/assignment.tsv').write_text(low.replace('\t12\t', '\t20000\t'))
    assert children(t, PARENT)['pathogen']['status'] == 'quantified' and (t / 'dual/pathogen/config/config.yaml').is_file()
    (t / 'results/assignment.tsv').write_text(low)
    assert children(t, PARENT)['pathogen']['status'] == 'skipped' and not (t / 'dual/pathogen/config/config.yaml').exists()
print('dual-organism child projects passed')

# Launcher: the released single launch is unchanged; a dual project runs the parent, then each child.
with tempfile.TemporaryDirectory() as t:
    t = Path(t).resolve()
    (t / 'bin').mkdir(); (t / 'p').mkdir()
    (t / 'bin/snakemake').write_text('#!/bin/sh\necho "$PWD|$*" >> "$LOG"\n[ "$(basename "$PWD")" != "$FAIL" ] || exit 7\n')
    (t / 'bin/snakemake').chmod(stat.S_IRWXU)
    parent(t / 'p', {s: 12 if s in MOCK else 50000 for s in MOCK + INFECTED})
    counts = {**PARENT, 'dual_organism': {**PARENT['dual_organism'], 'host': {**PARENT['dual_organism']['host'], 'objectives': ['gene_expression']}}}
    (t / 'p/config/config.yaml').write_text(yaml.safe_dump(counts))
    (t / 'p/gates').mkdir()
    tail = (ROOT / 'submit.sh').read_text().split('# --- Launch ---', 1)[1]
    head = (f'set -euo pipefail\nREPO={ROOT}\nPROJDIR={t / "p"}\nCONFIG=config/config.yaml\nEXECUTOR=local\nCORES=2\n'
            'MAX_CONC=3\nDRY_RUN=$DRY\nPASSTHRU=(--keep-going)\ncd "$PROJDIR"\n#')
    args = f'--snakefile {ROOT}/Snakefile --executor local --cores 2 --scheduler greedy --configfile config/config.yaml --config repo_dir={ROOT} --jobs 3'

    def submit(route, fail='', dry=''):
        (t / 'p/gates/recommendation.json').write_text(json.dumps({'route': route}))
        (t / 'log').write_text('')
        env = dict(os.environ, PATH=f"{t / 'bin'}{os.pathsep}{os.environ['PATH']}", LOG=str(t / 'log'), FAIL=fail, DRY=dry)
        rc = subprocess.run(['bash', '-c', head + tail], env=env, capture_output=True, text=True).returncode
        return rc, (t / 'log').read_text().splitlines()

    assert submit('star_counts') == (0, [f'{t}/p|{args} --keep-going']), submit('star_counts')
    assert submit('star_counts', fail='p') == (7, [f'{t}/p|{args} --keep-going'])
    assert not (t / 'p/dual').exists()
    assert submit('dual_assignment', dry='--dry-run') == (0, [f'{t}/p|{args} --dry-run --keep-going'])
    assert submit('dual_assignment', fail='p') == (7, [f'{t}/p|{args} --keep-going']) and not (t / 'p/dual').exists()
    assert submit('dual_assignment') == (0, [f'{t}/p{d}|{args} --keep-going' for d in ('', '/dual/host', '/dual/pathogen')])
    assert submit('dual_assignment', fail='host') == (7, [f'{t}/p{d}|{args} --keep-going' for d in ('', '/dual/host')])
    assert submit('dual_assignment', fail='pathogen')[0] == 7
print('dual-organism launcher passed')

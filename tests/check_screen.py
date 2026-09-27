#!/usr/bin/env python3
"""Known exclusive/shared/unmapped reads, mate sampling, and no read deletion."""
import gzip
import json
from pathlib import Path
import random
import shutil
import subprocess
import sys
import tempfile
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from screen_reads import screen, sample_reads, validate, index_files
for bad in ({'enabled': 'true'}, {'enabled': True}, {'fragments': 0}, {'seed': -1},
            {'references': [{'name':'X', 'index':'x', 'role':'unexpected'}]}):
    try: validate(bad)
    except ValueError: pass
    else: raise AssertionError(bad)
with tempfile.TemporaryDirectory() as temp:
    p = Path(temp); source = p / 'source with spaces.fq.gz'
    with gzip.open(source, 'wt') as f:
        for i in range(50): f.write(f'@{i}\nACGTACGT\n+\nIIIIIIII\n')
    a = p/'a'; b=p/'b'; a.mkdir(); b.mkdir()
    assert sample_reads([source,source], 7, 9, a) == (50,7)
    sample_reads([source,source], 7, 9, b)
    assert (a/'R1.fastq').read_bytes() == (b/'R1.fastq').read_bytes() == (a/'R2.fastq').read_bytes()
    mismatched = p/'mismatched.fq.gz'
    with gzip.open(mismatched, 'wt') as f: f.write('@wrong\nACGTACGT\n+\nIIIIIIII\n')
    try: sample_reads([source,mismatched], 7, 9, a)
    except ValueError: pass
    else: raise AssertionError('mismatched mates accepted')
    before = source.read_bytes()
    result = screen([source], p/'disabled', 'sample', {})
    assert result['status']=='not_performed' and source.read_bytes()==before
    try: screen([source],p,'sample',{})
    except ValueError: pass
    else: raise AssertionError('source/output collision allowed')
if '--real' not in sys.argv:
    print('Screen policy, sampling, disabled status and source protection passed'); sys.exit(0)
if any(not shutil.which(t) for t in ('fastq_screen','bowtie2','bowtie2-build')):
    print('Needs pinned screening module');sys.exit(77)
with tempfile.TemporaryDirectory(prefix='screen with spaces ') as temp:
    p=Path(temp); rng=random.Random(345)
    seq=lambda: ''.join(rng.choices('ACGT',k=180))
    host,other,shared,none = [seq() for _ in range(4)]
    refs=[]
    for name,unique,role in [('Host',host,'expected'),('Other',other,'possible_contaminant')]:
        fasta=p/f'{name}.fa'; fasta.write_text(f'>unique\n{unique}\n>shared\n{shared}\n')
        idx=p/name
        subprocess.run(['bowtie2-build',str(fasta),str(idx)],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
        assert len(index_files(str(idx)))==6
        refs.append(dict(name=name,role=role,index=str(idx)))
    source=p/'reads.fq.gz'
    with gzip.open(source,'wt') as f:
        for i,s in enumerate([host]*3+[other]*2+[shared]*4+[none]):
            f.write(f'@{i}\n{s[:100]}\n+\n'+('I'*100)+'\n')
    before=source.read_bytes()
    policy=dict(enabled=True,fragments=100,seed=7,references=refs)
    for reads,label in [([source],'single'),([source,source],'paired')]:
        result=screen(reads,p/label,label,policy)
        assert result['fragments_sampled']==result['fragments_available']==10
        assert result['reads_removed']==0 and source.read_bytes()==before
        for row in result['results']:
            assert row['reads']==10 and row['shared_with_other_references']==4,row
            assert row['exclusive_to_reference']==(3 if row['reference']=='Host' else 2),row
        assert len(result['results'])==2*len(reads)
        assert json.loads((p/label/'screening.json').read_text())['status']=='completed'
print('Real screen SE/PE exclusive/shared controls and original read preservation passed')

#!/usr/bin/env python3
"""Known gene/junction/intron fragments through real STAR and featureCounts."""
import csv
import gzip
from pathlib import Path
import random
import shutil
import sys
import tempfile
ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT/'scripts'))
from star_counts import build_index,quantify
if any(not shutil.which(tool) for tool in ('STAR','featureCounts','samtools')):
    print('Needs STAR, featureCounts and samtools');sys.exit(77)
rng=random.Random(42)
def sequence(n):return ''.join(rng.choices('ACGT',k=n))
def reverse(s):return s.translate(str.maketrans('ACGT','TGCA'))[::-1]
with tempfile.TemporaryDirectory() as tmp:
    p=Path(tmp);genome=p/'genome.fa';gtf=p/'genes.gtf'
    fasta=[];annotation=[];fragments=[]
    for gene in ('A','B'):
        first,last,intron=sequence(300),sequence(300),'GT'+sequence(196)+'AG'
        fasta.append(f'>chr{gene}\n{first+intron+last}\n')
        for start,end in [(1,300),(501,800)]:
            annotation.append(f'chr{gene}\tfixture\texon\t{start}\t{end}\t.\t+\t.\tgene_id "{gene}"; transcript_id "{gene}_tx";\n')
        transcript=first+last
        fragments += [(f'{gene}_{i}',transcript[215+i:415+i]) for i in range(8)]
        if gene=='A':fragments.append(('intronic',intron[15:175]))
    fragments.append(('unmapped','A'*200))
    genome.write_text(''.join(fasta));gtf.write_text(''.join(annotation))
    r1,r2=p/'R1.fq.gz',p/'R2.fq.gz'
    with gzip.open(r1,'wt') as a,gzip.open(r2,'wt') as b:
        for name,frag in fragments:
            a.write(f'@{name}/1\n{frag[:100]}\n+\n'+('I'*100)+'\n')
            b.write(f'@{name}/2\n{reverse(frag[-100:])}\n+\n'+('I'*100)+'\n')
    idx=p/'index'
    idx.mkdir()  # Snakemake pre-creates parents of nested declared outputs.
    record=build_index(genome,gtf,idx,2,99)
    assert build_index(genome,gtf,idx,2,99)==record
    for layout,mate,strand in [('paired',r2,1),('single','',1),('wrong-strand',r2,2)]:
        (p/layout).mkdir()
        result=quantify(r1,mate,idx,gtf,p/layout,strand,2)
        with (p/layout/'gene_counts.tsv').open() as f:
            counts={row['gene_id']:int(row['count']) for row in csv.DictReader(f,delimiter='\t')}
        expected={'A':8,'B':8} if strand==1 else {'A':0,'B':0}
        assert counts==expected,(layout,counts)
        assert result['num_reads_processed']==18
        assert result['num_assigned']==sum(expected.values())
        assert result['assignment_summary']['Unassigned_NoFeatures']>=1
        assert (p/layout/'Aligned.sortedByCoord.out.bam.bai').is_file()
    with (idx/'SAindex').open('ab') as f:f.write(b'tamper')
    try:build_index(genome,gtf,idx,2,99)
    except ValueError as exc:assert 'integrity' in str(exc)
    else:raise AssertionError('tampered index reused')
print('STAR junction mapping, SE/PE fragment counts, strand and intron controls passed')

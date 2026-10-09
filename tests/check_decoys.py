#!/usr/bin/env python3
"""Build a real decoy index and verify reference identity and cache reuse."""
import json
from pathlib import Path
import random
import shutil
import subprocess
import sys
import tempfile
ROOT = Path(__file__).resolve().parent.parent
sys.path[:0] = [str(ROOT/'workflow/scripts'), str(ROOT/'scripts')]
from build_salmon_index import prepare_genome_decoys
from reference_cache import key_inputs, verify_entry
with tempfile.TemporaryDirectory() as tmp:
    p = Path(tmp); rng = random.Random(21)
    tx, genome, gtf = p/'transcriptome.fa', p/'genome.fa', p/'annotation.gtf'
    tx.write_text('>tx1\n'+''.join(rng.choices('ACGT',k=2000))+'\n')
    genome.write_text('>chr1\n'+''.join(rng.choices('ACGT',k=3000))+'\n')
    gtf.write_text('chr1\ttest\texon\t1\t2000\t.\t+\t.\tgene_id "gene1"; transcript_id "tx1";\n')
    combined, names = p/'gentrome.fa', p/'decoys.txt'
    prepare_genome_decoys(tx,genome,combined,names)
    assert combined.read_text() == tx.read_text()+genome.read_text()
    assert names.read_text() == 'chr1\n'
    genome.write_text(tx.read_text())
    try: prepare_genome_decoys(tx,genome,combined,names)
    except ValueError: pass
    else: raise AssertionError('colliding genomic names accepted')
    genome.write_text('>chr1\n'+''.join(rng.choices('ACGT',k=3000))+'\n')
    if not shutil.which('salmon'):
        print('Salmon required for real index qualification'); sys.exit(77)
    ref={'accession':'fixture', 'transcriptome_fasta_url':tx.as_uri(), 'gtf_url':gtf.as_uri(), 'genome_fasta_url':genome.as_uri()}
    cmd=[sys.executable,str(ROOT/'workflow/scripts/build_salmon_index.py'),
         '--transcriptome',str(tx),'--gtf',str(gtf),'--genome',str(genome),
         '--final-index',str(p/'salmon_idx'),'--staging',str(p/'.building'),
         '--lock-out',str(p/'reference.lock.json'),'--helper',str(ROOT/'workflow/scripts/reference_lock.py'),
         '--decoys','genome','--threads','2','--accession','fixture',
         '--transcriptome-url',tx.as_uri(),'--gtf-url',gtf.as_uri(),'--genome-url',genome.as_uri(),
         '--key-inputs',json.dumps(key_inputs(ref,31,'genome'))]
    subprocess.run(cmd,check=True,capture_output=True)
    lock=json.loads((p/'reference.lock.json').read_text())
    assert lock['index']['decoys']['used']
    assert lock['identifier_compatibility']['decoys']==1
    assert not verify_entry(p,ref,31,'genome')
    subprocess.run(cmd,check=True,capture_output=True)
    genome.write_text(genome.read_text()+'A\n')
    assert verify_entry(p,ref,31,'genome')
print('Real genome-decoy index, reuse and tamper rejection passed')

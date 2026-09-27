#!/usr/bin/env python3
"""STAR genomic alignment and featureCounts gene counting for non-UMI bulk reads.

Consumes preprocessed gzip FASTQ. Strandedness is explicit: 0 unstranded,
1 forward, 2 reverse. Paired libraries count fragments, not individual mates.
"""
import argparse
import csv
import hashlib
import json
import math
import shutil
from pathlib import Path
import subprocess
import tempfile

from reference_cache import BuildLock, publish


def cache_key(reference):
    inputs = {k: reference.get(k) for k in ("genome_fasta_url", "gtf_url", "star_overhang")}
    inputs["STAR"] = "2.7.11b"
    return "star-" + hashlib.sha256(json.dumps(inputs,sort_keys=True).encode()).hexdigest()[:16]


def verify_index(destination, genome=None, annotation=None):
    dest = Path(destination)
    saved = json.loads((dest/'reference.lock.json').read_text())
    for key, path in (("genome", genome), ("annotation", annotation)):
        if path is not None and identity(path)['sha256'] != saved['request'][key]['sha256']:
            raise ValueError(f'STAR indexed {key} differs from current reference')
    for name,digest in saved['index_sha256'].items():
        if not (dest/name).is_file() or identity(dest/name)['sha256'] != digest:
            raise ValueError(f'STAR index integrity failure: {name}')
    return saved


def identity(path):
    path = Path(path).resolve()
    with path.open('rb') as f:
        digest = hashlib.file_digest(f, 'sha256').hexdigest()
    return dict(path=str(path), bytes=path.stat().st_size, sha256=digest)


def version(tool, flag='--version'):
    p = subprocess.run([tool, flag], check=True, capture_output=True, text=True)
    return (p.stdout + p.stderr).strip().splitlines()[0]


def run(command, log):
    with Path(log).open('w') as f:
        subprocess.run(command, check=True, stdout=f, stderr=subprocess.STDOUT)


def build_index(genome, gtf, destination, threads=4, overhang=99):
    if threads < 1 or overhang < 1:
        raise ValueError('threads and splice-junction overhang must be positive')
    genome, gtf, dest = Path(genome).resolve(), Path(gtf).resolve(), Path(destination).resolve()
    names, bases = set(), 0
    with genome.open() as f:
        for line in f:
            if line.startswith('>'):
                name = line[1:].split()[0]
                if name in names:
                    raise ValueError(f'Duplicate genome contig: {name}')
                names.add(name)
            else:
                bases += len(line.strip())
    if not names or not bases:
        raise ValueError('Empty genome FASTA')
    exons = 0
    with gtf.open() as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            row = line.rstrip().split('\t')
            if len(row) != 9 or row[0] not in names:
                raise ValueError('GTF row/contig does not match genome FASTA')
            exons += row[2] == 'exon'
    if not exons:
        raise ValueError('Gene counting requires exon annotations')
    star_version = version('STAR')
    request = dict(genome=identity(genome), annotation=identity(gtf),
                   STAR=star_version, overhang=overhang)
    dest.parent.mkdir(parents=True, exist_ok=True)
    with BuildLock(dest):
        # Snakemake creates parent directories for the declared child outputs.
        if dest.is_dir() and not any(dest.iterdir()):
            dest.rmdir()
        if dest.exists():
            manifest = dest/'reference.lock.json'
            if not manifest.is_file():
                raise ValueError('Incomplete existing STAR index; use a new destination')
            saved = verify_index(dest)
            if saved['request'] != request:
                raise ValueError('Existing STAR index has different reference or build settings')
            return saved
        with tempfile.TemporaryDirectory(prefix='.star-',dir=dest.parent) as tmp:
            stage = Path(tmp)/'index'; stage.mkdir()
            cmd = ['STAR','--runMode','genomeGenerate','--runThreadN',str(threads),
                   '--genomeDir',str(stage),'--genomeFastaFiles',str(genome),
                   '--sjdbGTFfile',str(gtf),'--sjdbOverhang',str(overhang),
                   '--genomeSAindexNbases',str(min(14,max(1,int(math.log2(bases)/2-1))))]
            run(cmd,stage/'build.log')
            required = ['Genome','SA','SAindex','genomeParameters.txt']
            if not all((stage/name).is_file() and (stage/name).stat().st_size for name in required):
                raise ValueError('STAR did not produce a complete index')
            saved = dict(schema_version=1,request=request,command=cmd,
                         index_sha256={p.name:identity(p)['sha256'] for p in stage.iterdir()
                                       if p.is_file() and not p.name.endswith(('.log', '.out'))})
            (stage/'reference.lock.json').write_text(json.dumps(saved,indent=2)+'\n')
            publish(stage,dest)
    return saved


def quantify(r1, r2, index, gtf, destination, strand, threads=4, sample=None, replace_owned=False):
    if strand not in (0,1,2) or threads < 1:
        raise ValueError('Declare strand 0, 1 or 2 and a positive thread count')
    index, gtf, dest = Path(index).resolve(), Path(gtf).resolve(), Path(destination).resolve()
    if any(Path(x).resolve().is_relative_to(dest) for x in (r1,r2,index,gtf) if x):
        raise ValueError('Read/reference source cannot be inside the STAR output directory')
    if dest.exists() and (not dest.is_dir() or any(dest.iterdir())):
        marker = dest/'metrics.json'
        if not replace_owned or not marker.is_file() or json.loads(marker.read_text()).get('producer') != 'STAR_featureCounts':
            raise ValueError('Use a new STAR output directory; only completed owned outputs can be replaced')
    record = json.loads((index/'reference.lock.json').read_text())
    if record['request']['annotation']['sha256'] != identity(gtf)['sha256']:
        raise ValueError('Counting annotation differs from the indexed annotation')
    tool_versions = {'STAR':version('STAR'),'featureCounts':version('featureCounts','-v'),
                     'samtools':version('samtools')}
    if tool_versions['STAR'] != record['request']['STAR']:
        raise ValueError('STAR version differs from index builder')
    reads = [str(Path(r1).resolve())] + ([str(Path(r2).resolve())] if r2 else [])
    inputs = [identity(path) for path in reads]
    if dest.exists():
        shutil.rmtree(dest)
    dest.mkdir(parents=True)
    cmd = ['STAR','--runThreadN',str(threads),'--genomeDir',str(index),
           '--readFilesIn',*reads,'--readFilesCommand','zcat',
           '--twopassMode','Basic','--outFilterMultimapNmax','1','--outSAMtype','BAM','SortedByCoordinate',
           '--outSAMattributes','NH','HI','AS','nM','--outSAMmapqUnique','60',
           '--outSAMunmapped','Within','--limitBAMsortRAM','4000000000',
           '--outFileNamePrefix',str(dest)+'/']
    run(cmd,dest/'alignment.log')
    bam = dest/'Aligned.sortedByCoord.out.bam'
    subprocess.run(['samtools','quickcheck',str(bam)],check=True)
    subprocess.run(['samtools','index','-@',str(threads),str(bam)],check=True)
    counts = dest/'featureCounts.tsv'
    count_cmd = ['featureCounts','-T',str(threads),'-a',str(gtf),'-t','exon','-g','gene_id',
                 '-s',str(strand),'-o',str(counts)]
    if r2:
        count_cmd += ['-p','--countReadPairs','-B','-C']
    count_cmd.append(str(bam))
    run(count_cmd,dest/'counting.log')
    with counts.open() as f:
        rows = list(csv.DictReader((line for line in f if not line.startswith('#')),delimiter='\t'))
    if not rows or len(rows[0]) != 7 or len({r['Geneid'] for r in rows}) != len(rows):
        raise ValueError('Unexpected featureCounts gene table')
    count_column = list(rows[0])[-1]
    measured = [(r['Geneid'],int(r[count_column]),int(r['Length'])) for r in rows]
    if any(n < 0 or length < 1 for _,n,length in measured):
        raise ValueError('Invalid gene counts/lengths')
    with Path(str(counts)+'.summary').open() as f:
        summary = {r[0]:int(r[1]) for r in list(csv.reader(f,delimiter='\t'))[1:]}
    if sum(n for _,n,_ in measured) != summary['Assigned']:
        raise ValueError('Gene sum disagrees with assigned fragments/reads')
    stats = {}
    for line in (dest/'Log.final.out').read_text().splitlines():
        if '|' in line:
            key,value = line.split('|',1); stats[key.strip()] = value.strip()
    processed = int(stats['Number of input reads'])
    mapped = int(stats['Uniquely mapped reads number']) + int(stats['Number of reads mapped to multiple loci'])
    if processed < 1 or not 0 <= summary['Assigned'] <= processed or sum(summary.values()) != processed:
        raise ValueError('Alignment/counting unit conservation failed')
    for filename,column,label in [('gene_counts.tsv',1,'count'),('gene_lengths.tsv',2,'exonic_bases')]:
        with (dest/filename).open('w') as f:
            writer=csv.writer(f,delimiter='\t');writer.writerow(['gene_id',label])
            writer.writerows((row[0],row[column]) for row in measured)
    result=dict(schema_version=1,sample=sample or dest.name,num_mapped=mapped,producer='STAR_featureCounts',versions=tool_versions,
                count_unit='fragments' if r2 else 'reads',strandedness=strand,
                num_reads_processed=processed,num_assigned=summary['Assigned'],
                mapping_rate=mapped/processed,assignment_rate=summary['Assigned']/processed,
                assignment_summary=summary,input_files=inputs,
                reference_lock_sha256=identity(index/'reference.lock.json')['sha256'],
                commands=[cmd,count_cmd],multimapping='excluded',overlapping_genes='excluded',
                paired_policy='both mates mapped; chimeric pairs excluded' if r2 else None)
    (dest/'metrics.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


if __name__ == '__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    sub=ap.add_subparsers(dest='command',required=True)
    idx=sub.add_parser('index');idx.add_argument('--genome',required=True)
    idx.add_argument('--gtf',required=True);idx.add_argument('--out',required=True)
    idx.add_argument('--threads',type=int,default=4);idx.add_argument('--overhang',type=int,default=99)
    q=sub.add_parser('quant');q.add_argument('--r1',required=True);q.add_argument('--r2',default='')
    q.add_argument('--index',required=True);q.add_argument('--gtf',required=True)
    q.add_argument('--out',required=True);q.add_argument('--strand',type=int,choices=[0,1,2],required=True)
    q.add_argument('--threads',type=int,default=4)
    q.add_argument('--sample')
    q.add_argument('--replace-owned',action='store_true')
    a=ap.parse_args()
    if a.command=='index':build_index(a.genome,a.gtf,a.out,a.threads,a.overhang)
    else:quantify(a.r1,a.r2,a.index,a.gtf,a.out,a.strand,a.threads,a.sample,a.replace_owned)

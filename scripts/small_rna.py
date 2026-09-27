#!/usr/bin/env python3
"""Small-RNA kit profiles, mature-miRNA references and read counting."""
import csv
import json
import sys
from pathlib import Path
import yaml

TABLE = yaml.safe_load((Path(__file__).resolve().parent.parent / 'config/small_rna_kits.yaml').read_text())


def load_kit(name):
    if name in TABLE.get('blocked', {}):
        raise ValueError(f'{name}: {TABLE["blocked"][name]}')
    if name not in TABLE['kits']:
        raise ValueError(f'Unknown small-RNA kit {name!r}; choose one of {sorted(TABLE["kits"])}')
    return dict(TABLE['kits'][name])


def kit_policy(kit):
    return dict(adapter_r1=kit['adapter_r1'], minimum_length=kit['minimum_length'],
                length_limit=kit['length_limit'], poly_g='on')


def mature_reference(mature_fa, prefix, out):
    """Species mature miRNAs from miRBase mature.fa, RNA to DNA alphabet."""
    records, keep, name = [], False, None
    for line in Path(mature_fa).read_text().splitlines():
        if line.startswith('>'):
            name = line[1:].split()[0]
            keep = name.startswith(prefix + '-')
        elif keep:
            records.append((name, line.strip().upper().replace('U', 'T')))
    if not records:
        raise ValueError(f'No mature miRNAs with prefix {prefix!r}; check the miRBase species code')
    Path(out).write_text(''.join(f'>{n}\n{s}\n' for n, s in records))
    return len(records)


def count_alignments(sam, out_tsv):
    """Unique hits count per miRNA; a read aligned more than once is reported, never split."""
    hits, unaligned = {}, set()
    with open(sam) as fh:  # streamed: -a output for 20-30M reads does not fit in memory as text
        for line in fh:
            if line.startswith('@'):
                continue
            read, flag, ref = line.split('\t', 3)[:3]
            if int(flag) & 4:
                unaligned.add(read)
            else:
                hits.setdefault(read, []).append(ref)
    counts = {}
    for refs in hits.values():
        if len(refs) == 1:
            counts[refs[0]] = counts.get(refs[0], 0) + 1
    Path(out_tsv).write_text('gene_id\tcount\n' + ''.join(f'{g}\t{c}\n' for g, c in sorted(counts.items())))
    return dict(unique=sum(counts.values()), multi_mapped=sum(len(r) > 1 for r in hits.values()),
                unaligned=len(unaligned))


def count_metrics(sam, out_tsv, metrics_json, sample):
    """Counts plus a metrics.json in the shape 02_aggregate.R reads for STAR."""
    t = count_alignments(sam, out_tsv)
    reads = sum(t.values())
    if reads < 1:
        raise ValueError(f'{sam}: no reads in the alignment')
    mapped = t['unique'] + t['multi_mapped']
    m = dict(schema_version=1, sample=sample, producer='bowtie_mature_mirna', count_unit='reads',
             num_reads_processed=reads, num_mapped=mapped, num_assigned=t['unique'],
             mapping_rate=mapped / reads, assignment_rate=t['unique'] / reads,
             alignment_summary=t, multimapping='excluded')
    Path(metrics_json).write_text(json.dumps(m, indent=2) + '\n')
    return m


def shortstack_command(genome, reads, known, outdir, threads):
    """ShortStack 4 on pre-trimmed reads (fastp already trimmed, so never --autotrim)."""
    cmd = ['ShortStack', '--genomefile', str(genome), '--readfile', *map(str, reads),
           '--outdir', str(outdir), '--threads', str(threads)]
    return cmd + (['--known_miRNAs', str(known)] if known else [])


def shortstack_matrix(counts_txt, results_txt, out_counts, out_loci):
    """Counts.txt to gene_id + one column per sample (names from its header); loci table from Results.txt."""
    with open(counts_txt) as fh:
        rows = list(csv.reader(fh, delimiter='\t'))
    head = rows[0]
    if head[:3] != ['Coords', 'Name', 'MIRNA'] or len(head) < 4:
        raise ValueError(f'{counts_txt}: expected Coords, Name, MIRNA then one column per sample; got {head}')
    Path(out_counts).write_text('\t'.join(['gene_id'] + head[3:]) + '\n' +
                                ''.join('\t'.join([r[1]] + r[3:]) + '\n' for r in rows[1:]))
    keep = ['Name', 'Chrom', 'Start', 'End', 'Strand', 'DicerCall', 'MIRNA', 'Known_miRNAs']
    with open(results_txt) as fh:
        loci = list(csv.DictReader(fh, delimiter='\t'))
    Path(out_loci).write_text('\t'.join(keep) + '\n' + ''.join('\t'.join(r[k] for k in keep) + '\n' for r in loci))
    return len(rows) - 1


def matrix_metrics(matrix_tsv, preprocessing_json, sample, metrics_json):
    """Plant per-sample metrics in count_metrics' shape: reads in loci over reads after preprocessing."""
    with open(matrix_tsv) as fh:
        table = list(csv.DictReader(fh, delimiter='\t'))
    if not table or sample not in table[0]:
        raise ValueError(f'{matrix_tsv}: no column for sample {sample!r}; check ShortStack sample naming')
    in_loci = sum(int(r[sample]) for r in table)
    reads = json.loads(Path(preprocessing_json).read_text())['reads_after']
    if reads < 1:
        raise ValueError(f'{preprocessing_json}: no reads after preprocessing')
    m = dict(schema_version=1, sample=sample, producer='shortstack_loci', count_unit='reads',
             num_reads_processed=reads, num_mapped=in_loci, num_assigned=in_loci,
             mapping_rate=in_loci / reads, assignment_rate=in_loci / reads,
             multimapping='placed by ShortStack')
    Path(metrics_json).write_text(json.dumps(m, indent=2) + '\n')
    return m


if __name__ == '__main__':
    {'reference': mature_reference, 'count': count_metrics, 'shortstack_matrix': shortstack_matrix,
     'plant_metrics': matrix_metrics}[sys.argv[1]](*sys.argv[2:])

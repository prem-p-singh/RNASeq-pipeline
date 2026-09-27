#!/usr/bin/env python3
"""Small-RNA kit profiles, mature-miRNA references and read counting."""
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
    for line in Path(sam).read_text().splitlines():
        if line.startswith('@'):
            continue
        read, flag, ref = line.split('\t')[:3]
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


if __name__ == '__main__':
    {'reference': mature_reference, 'count': count_metrics}[sys.argv[1]](*sys.argv[2:])

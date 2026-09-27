#!/usr/bin/env python3
"""Synthetic small-RNA animal project for tests/check_small_rna_workflow.sh.

Usage: small_rna_workflow_project.py PROJECT_DIR REPO_DIR
Writes the project config, reads, a local mature-miRNA FASTA and expected.json.
"""
import csv
import gzip
import json
import sys
from pathlib import Path

import yaml

ADAPTER = 'TGGAATTCTCGGGTGCCAAGG'  # TruSeq Small RNA, config/small_rna_kits.yaml
FILLER = 'ACTCAAT'                 # 3' of the adapter; no G run for poly-G trimming
INSERTS = {
    'miR-1': 'TGAGGTAGTAGGTTGTATAGTT',
    'miR-3': 'TAGCTTATCAGACTGATGTTGA',
    'miR-2': 'TAGCAGCACGTAAATATTGGCG',   # fix-miR-2a and fix-miR-2b
    'oth': 'TAAAGTGCTTATAGTGCAGGTAG',    # oth-miR-1, filtered out by the prefix
    'random': 'GATCCATGACTTGCAGTCAAGC',  # in no record
    'short': 'TCAGTCAGTC',               # 10 nt, below minimum_length 18
}
DESIGN = {
    'S1': {'miR-1': 5, 'miR-3': 2, 'miR-2': 3, 'oth': 1, 'random': 1, 'short': 1},
    'S2': {'miR-1': 3, 'miR-3': 4, 'miR-2': 1, 'oth': 2, 'short': 1},
    'S3': {'miR-1': 6, 'miR-2': 2, 'random': 2, 'short': 1},
}


def rna(seq):
    return seq.replace('T', 'U')


def build(proj, repo):
    (proj / 'inputs').mkdir(parents=True)
    (proj / 'config').mkdir()
    m3 = rna(INSERTS['miR-3'])
    (proj / 'inputs/mature.fa').write_text(
        f">fix-miR-1 MIMAT9000001\n{rna(INSERTS['miR-1'])}\n"
        f">fix-miR-3 MIMAT9000003\n{m3[:12]}\n{m3[12:]}\n"
        f">fix-miR-2a MIMAT9000002\n{rna(INSERTS['miR-2'])}\n"
        f">fix-miR-2b MIMAT9000004\n{rna(INSERTS['miR-2'])}\n"
        f">oth-miR-1 MIMAT9000005\n{rna(INSERTS['oth'])}\n")
    rows, expected = [], {}
    for sid, design in DESIGN.items():
        reads = [(f'{sid}_{kind}_{i}', INSERTS[kind] + ADAPTER + FILLER)
                 for kind, n in design.items() for i in range(n)]
        path = proj / f'inputs/{sid}.fq.gz'
        with gzip.open(path, 'wt') as f:
            f.writelines(f'@{name}\n{seq}\n+\n{"I" * len(seq)}\n' for name, seq in reads)
        rows.append([sid, str(path), ''])
        counts = {f'fix-{k}': design[k] for k in ('miR-1', 'miR-3') if design.get(k)}
        expected[sid] = dict(
            counts=counts, reads_before=len(reads),
            alignment=dict(unique=sum(counts.values()), multi_mapped=design.get('miR-2', 0),
                           unaligned=design.get('oth', 0) + design.get('random', 0)))
    with open(proj / 'config/samples.tsv', 'w') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['sample_id', 'fastq_url', 'fastq_url_r2'])
        w.writerows(rows)
    cfg = {
        'project': {'name': 'small_rna_workflow', 'output_dir': 'results/'},
        'organism': {'scientific_name': 'Fixtura synthetica', 'tax_id': 1, 'kegg_code': None,
                     'orgdb_package': 'org.Fixture.eg.db', 'mirbase_prefix': 'fix'},
        'reference': {'accession': 'FIXTURE_SRNA_1', 'cache_dir': None,
                      'mature_mirna_url': (proj / 'inputs/mature.fa').as_uri()},
        'samples': {'seq_type': 'srnaseq', 'sheet': 'config/samples.tsv'},
        'analysis': {'objectives': ['gene_expression'], 'small_rna_profile': 'animal',
                     'small_rna_kit': 'TruSeq Small RNA'},
        'model': {'fixed_effects': '~ 1', 'random_effects': None, 'primary_factor': ''},
        'contrasts': [],
        'downstream': {'run_go': False, 'run_kegg': False, 'run_wgcna': False},
        'orgdb': {'strategy': 'skip'},
        'hpc': {'delete_fastq_after_quant': False},
        'thresholds': {'sample_qc': {'min_reads_on_genes_srnaseq': 1}},
    }
    (proj / 'config/config.yaml').write_text(yaml.safe_dump(cfg))
    (proj / 'config/thresholds.yaml').write_text((repo / 'config/thresholds.yaml').read_text())
    (proj / 'expected.json').write_text(json.dumps(expected, indent=2) + '\n')
    return expected


if __name__ == '__main__':
    print(json.dumps(build(Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve()), indent=2))

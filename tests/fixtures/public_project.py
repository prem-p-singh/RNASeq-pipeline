#!/usr/bin/env python3
"""Pinned nf-core/GSE110004 yeast smoke dataset, three biological units per group.

Uses accession-level biological metadata from the pinned dataset README, not the
artificial lane/group labels in nf-core's workflow-test samplesheet. Chromosome-I
filtered/downsampled reads test operation, not genome-wide biological conclusions.
"""
import csv
import hashlib
import json
import subprocess
from pathlib import Path
import sys
import yaml

REVISION = "626c8fab639062eade4b10747e919341cbf9b41a"
BASE = f"https://raw.githubusercontent.com/nf-core/test-datasets/{REVISION}"

dest, repo = Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve()
flags = sys.argv[3:]
if set(flags) - {'--star', '--screen'} or len(flags) != len(set(flags)):
    raise SystemExit('Usage: public_project.py PROJECT REPO [--star] [--screen]')
star_route = '--star' in flags
(dest / "config").mkdir(parents=True)
(dest / "inputs").mkdir()
# The upstream FASTA adds a GFP transgene absent from genes.gtf.gz. Build an
# explicitly documented yeast-only fixture; do not weaken production ID checks.
original = dest / "inputs/upstream_transcriptome.fasta"
source_url = BASE + "/reference/transcriptome.fasta"
subprocess.run(["curl", "--fail", "--retry", "3", "-sSL", source_url, "-o", str(original)], check=True)
records = original.read_text().split(">")[1:]
removed = [r.split()[0] for r in records if r.split()[0] == "Gfp_transgene_gene"]
assert removed == ["Gfp_transgene_gene"], removed
transcripts = dest / "inputs/yeast_transcriptome.fasta"
transcripts.write_text("".join(">" + r for r in records if r.split()[0] != "Gfp_transgene_gene"))
(dest / "inputs/reference_derivation.json").write_text(json.dumps({
    "source_url": source_url, "source_sha256": hashlib.sha256(original.read_bytes()).hexdigest(),
    "removed_unannotated_test_transgenes": removed,
    "derived_sha256": hashlib.sha256(transcripts.read_bytes()).hexdigest(),
}, indent=2) + "\n")
cfg = {
    "project": {"name": "GSE110004_chrI_smoke", "output_dir": "results/"},
    "organism": {"scientific_name": "Saccharomyces cerevisiae", "tax_id": 4932,
                 "kegg_code": "sce", "orgdb_package": "org.Sc.sgd.db"},
    "reference": {"accession": "R64-1-1_chrI", "cache_dir": None,
                  "transcriptome_fasta_url": transcripts.as_uri(),
                  "genome_fasta_url": None,
                  "gtf_url": BASE + "/reference/genes.gtf.gz",
                  "annotation_tsv": {"path": None}},
    "samples": {"seq_type": "rnaseq_paired", "sheet": "config/samples.tsv"},
    "model": {"fixed_effects": "~ treatment", "random_effects": None, "primary_factor": "treatment"},
    "contrasts": [{"id": "genotype", "type": "pairwise", "factor": "treatment", "reverse": True}],
    "downstream": {"run_go": False, "run_kegg": False, "run_wgcna": False},
    "orgdb": {"strategy": "skip"},
    "hpc": {"delete_fastq_after_quant": False},
    # Explicit smoke-test thresholds for subsampled chromosome-I data only.
    "thresholds": {"sample_qc": {"min_reads_on_genes_rnaseq": 1, "mapping_rate_min": 0.20}},
}
if star_route:
    # Upstream deliberately includes a not_in_genome contig in the GTF.
    # Derive a compatible fixture; keep production reference checks strict.
    genome = dest/'inputs/genome.fa'
    upstream_gtf = dest/'inputs/upstream_genes.gtf'
    for name, target in [('genome.fasta', genome), ('genes.gtf', upstream_gtf)]:
        subprocess.run(['curl', '--fail', '--retry', '3', '-sSL', BASE+'/reference/'+name,
                        '-o', str(target)], check=True)
    contigs = {line[1:].split()[0] for line in genome.read_text().splitlines() if line.startswith('>')}
    lines = upstream_gtf.read_text().splitlines(keepends=True)
    removed_lines = [line for line in lines if line.strip() and not line.startswith('#') and line.split('\t')[0] not in contigs]
    assert {line.split('\t')[0] for line in removed_lines} == {'not_in_genome'}
    annotation = dest/'inputs/yeast_genes.gtf'
    annotation.write_text(''.join(line for line in lines if line not in removed_lines))
    derivation = dest/'inputs/reference_derivation.json'
    record = json.loads(derivation.read_text())
    record['genomic_annotation'] = {
        'source_url': BASE+'/reference/genes.gtf',
        'source_sha256': hashlib.sha256(upstream_gtf.read_bytes()).hexdigest(),
        'removed_contigs': ['not_in_genome'], 'removed_rows': len(removed_lines),
        'derived_sha256': hashlib.sha256(annotation.read_bytes()).hexdigest(),
        'genome_source_url': BASE+'/reference/genome.fasta',
        'genome_sha256': hashlib.sha256(genome.read_bytes()).hexdigest(),
    }
    derivation.write_text(json.dumps(record, indent=2)+'\n')
    cfg['reference']['genome_fasta_url'] = genome.as_uri()
    cfg['reference']['gtf_url'] = annotation.as_uri()
    cfg['samples']['expected_libtype'] = 'ISR'
    cfg['analysis'] = {'quantifier': 'star', 'backend': 'edger_ql'}
if '--screen' in flags:
    index = dest/'inputs/yeast_screen'
    with (dest/'inputs/screen_index.log').open('w') as log:
        subprocess.run(['bowtie2-build', str(transcripts), str(index)], check=True,
                       stdout=log, stderr=subprocess.STDOUT)
    cfg['screening'] = {'enabled': True, 'fragments': 10000, 'seed': 1,
        'references': [{'name': 'Yeast', 'role': 'expected', 'index': str(index)}]}
(dest / "config/config.yaml").write_text(yaml.safe_dump(cfg))
(dest / "config/thresholds.yaml").write_text((repo / "config/thresholds.yaml").read_text())
with open(dest / "config/samples.tsv", "w") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["sample_id", "fastq_url", "fastq_url_r2", "treatment"])
    for i in range(6):
        acc = f"SRR{6357070+i}"
        w.writerow([acc, f"{BASE}/testdata/GSE110004/{acc}_1.fastq.gz",
                    f"{BASE}/testdata/GSE110004/{acc}_2.fastq.gz", "WT" if i < 3 else "Rap1_AID_uninduced"])
(dest / "DATA_SOURCE.md").write_text(
    f"Source: {BASE}/README.md\n\nGSE110004; SRR6357070-SRR6357075, six distinct biological replicates.\n"
    "Filtered to yeast chromosome I and downsampled by nf-core. Technical smoke test only; "
    "no claim of reproducing the paper's genome-wide differential expression.\n\n"
    "Reference derivation: upstream transcriptome includes Gfp_transgene_gene without a matching "
    "GTF record. This fixture removes that one artificial transgene, preserving all yeast records. "
    "Original and derived checksums are recorded in inputs/reference_derivation.json.\n")
if star_route:
    with (dest/'DATA_SOURCE.md').open('a') as f:
        f.write(f'\nSTAR uses the pinned genomic reference and reverse strandedness documented in '
                f'{BASE}/samplesheet/samplesheet.csv. Only strandedness is taken from that '
                'workflow-test sheet; biological identities remain accession-level README metadata.\n')
        f.write('The fixture GTF excludes upstream not_in_genome rows; the genome lacks that '
                'test contig. Original/derived hashes and removed-row count are recorded.\n')

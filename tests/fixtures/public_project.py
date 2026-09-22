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

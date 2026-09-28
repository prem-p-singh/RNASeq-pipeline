#!/usr/bin/env python3
"""Build a project from prepared small-RNA public data (see prepare_small_rna_public.py
and prepare_plant_small_rna_public.py). gene_expression only: these are technical
controls, not a biological design, so no DE model is fitted.
"""
import csv
import json
from pathlib import Path
import sys

import yaml

BASE_CFG = dict(
    project={"output_dir": "results/"},
    analysis={"objectives": ["gene_expression"], "small_rna_kit": "TruSeq Small RNA"},
    samples={"seq_type": "srnaseq", "sheet": "config/samples.tsv"},
    model={"fixed_effects": "~ 1", "random_effects": None, "primary_factor": ""},
    contrasts=[],
    downstream={"run_go": False, "run_kegg": False, "run_wgcna": False},
    orgdb={"strategy": "skip"},
    hpc={"delete_fastq_after_quant": False},
)


def _write(dest, cfg, rows, repo):
    (dest / "config").mkdir(parents=True)
    (dest / "config/config.yaml").write_text(yaml.safe_dump(cfg))
    (dest / "config/thresholds.yaml").write_text((repo / "config/thresholds.yaml").read_text())
    with open(dest / "config/samples.tsv", "w") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["sample_id", "fastq_url", "fastq_url_r2"])
        w.writerows(rows)


def build_animal(dest, data, repo):
    """Publication's 331 synthetic sequences (not miRBase); relabelled so the standard
    mature_mirna_reference rule's prefix filter accepts them unchanged."""
    manifest = json.loads((repo / "tests/fixtures/small_rna_public.json").read_text())
    reference = (data / "reference.fa").read_text()
    records = reference.split(">")[1:]
    relabelled = "".join(f">pub-{r}" for r in records)
    (dest / "inputs").mkdir(parents=True)
    (dest / "inputs/mature_source.fa").write_text(relabelled)
    cfg = dict(BASE_CFG)
    cfg["project"] = dict(cfg["project"], name="small_rna_animal_public")
    cfg["organism"] = {"mirbase_prefix": "pub"}
    cfg["reference"] = {"cache_dir": None, "mature_mirna_url": (dest / "inputs/mature_source.fa").as_uri()}
    cfg["analysis"] = dict(cfg["analysis"], small_rna_profile="animal")
    rows = [[run["pool"], (data / f"{run['accession']}.fastq.gz").as_uri(), ""] for run in manifest["runs"]]
    _write(dest, cfg, rows, repo)


def build_plant(dest, data, repo):
    """Genome + optional known-miRNA reference from the real TAIR10/miRBase sources;
    the standard mature_mirna_reference rule filters mature.fa.all to the ath prefix."""
    manifest = json.loads((repo / "tests/fixtures/plant_small_rna_public.json").read_text())
    (dest / "inputs").mkdir(parents=True)
    cfg = dict(BASE_CFG)
    cfg["project"] = dict(cfg["project"], name="small_rna_plant_public")
    cfg["organism"] = {"mirbase_prefix": manifest["mirbase_prefix"]}
    cfg["reference"] = {"cache_dir": None,
                         "genome_fasta_url": (data / "genome.fa.gz").as_uri(),
                         "mature_mirna_url": (data / "mature.fa.all").as_uri()}
    cfg["analysis"] = dict(cfg["analysis"], small_rna_profile="plant")
    rows = [[run["accession"], (data / f"{run['accession']}.fastq.gz").as_uri(), ""] for run in manifest["runs"]]
    _write(dest, cfg, rows, repo)


if __name__ == "__main__":
    mode, data_dir, dest, repo = sys.argv[1], Path(sys.argv[2]).resolve(), Path(sys.argv[3]).resolve(), Path(sys.argv[4]).resolve()
    {"--animal": build_animal, "--plant": build_plant}[mode](dest, data_dir, repo)

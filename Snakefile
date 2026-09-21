# =============================================================================
# Generic RNA-Seq Pipeline — Snakemake orchestrator
#
# Entry point for all stages. Reads config.yaml + samples.tsv and dispatches
# to per-sample (Stage 1) or aggregate (Stages 2-5) rules.
#
# Usage:
#   ./submit.sh                     # auto-selects tier, launches via SLURM
#   snakemake -n                    # dry-run, show DAG
#   snakemake --cores 4             # local run (small datasets only)
# =============================================================================

import os
import pandas as pd
from pathlib import Path

# --- Config --------------------------------------------------------------
configfile: "config/config.yaml"

# Merge threshold defaults into config
import yaml
with open("config/thresholds.yaml") as f:
    _thresh_defaults = yaml.safe_load(f)
config.setdefault("thresholds", {})
for k, v in _thresh_defaults.items():
    config["thresholds"].setdefault(k, v)

# --- Load sample sheet ---------------------------------------------------
samples = pd.read_csv(
    config["samples"]["sheet"],
    sep="\t",
    comment="#",
).set_index("sample_id", drop=False)

SAMPLES = samples.index.tolist()
N_SAMPLES = len(SAMPLES)

# --- Output roots --------------------------------------------------------
OUT = Path(config["project"]["output_dir"])
REF = Path("reference")
QUANT = OUT / "quant"
METRICS = Path("metrics")
GATES = Path("gates")

for d in (OUT, REF, QUANT, METRICS, GATES):
    d.mkdir(parents=True, exist_ok=True)

# --- Rule modules --------------------------------------------------------
include: "workflow/rules/common.smk"
include: "workflow/rules/retrieve.smk"
include: "workflow/rules/build_orgdb.smk"
include: "workflow/rules/per_sample.smk"
include: "workflow/rules/qc_report.smk"
include: "workflow/rules/aggregate.smk"

# --- Conditional target helpers -----------------------------------------
def all_targets():
    """Build the final target list. OrgDb is only required when enrichment runs."""
    targets = [
        # Stage 1 — one quant.sf per sample
        *expand(str(QUANT / "{sample}/quant.sf"), sample=SAMPLES),
        # Stage 1b — comparative QC report (before/after-clean charts + MultiQC)
        OUT / "qc_report/qc_report.done",
        # Stage 2 — aggregated counts
        OUT / "counts.tsv",
        # Stage 3 — DE results
        OUT / "de_done.flag",
        # Stage 4 — enrichment
        OUT / "enrichment_done.flag",
        # Stage 5 — WGCNA (may be skipped; sentinel always emitted)
        OUT / "wgcna_done.flag",
    ]
    # Only ask for the OrgDb when GO enrichment will actually consume it.
    if orgdb_required():
        targets.append(orgdb_sentinel())
    return targets


# --- Final targets -------------------------------------------------------
# default_target makes this THE target for a bare `snakemake` call. Without
# it, Snakemake uses whichever rule it saw first, and because the rule modules
# are included above, that is fetch_transcriptome -- a bare launch would stop
# after downloading the reference.
rule all:
    default_target: True
    input:
        all_targets(),

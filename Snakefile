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

# Resolve defaults through the one shared path (R17e, R05; master plan 6.3).
#
# This used to be a one-level `setdefault` over config/thresholds.yaml. Because
# config.template.yaml carried its own partial `thresholds.wgcna`, the richer
# definition in thresholds.yaml was never merged and 05_wgcna.R read NULL for
# both merge heights. WGCNA swallows that error and returns unmerged modules,
# so the run looked successful while the requested merging never happened.
#
# scripts/config_resolve.py deep-merges, so overriding one key keeps its
# siblings, and reports unknown keys and type conflicts rather than absorbing
# them. preflight surfaces those same findings as CFG001/CFG003 before launch;
# here they are a hard stop, because by this point a run is starting.
import sys as _sys
_repo = Path(config.get("repo_dir", ".")).resolve()
_sys.path.insert(0, str(_repo / "scripts"))
import config_resolve as _cr

_project_layer = {k: v for k, v in config.items() if k != "repo_dir"}
config, CONFIG_ORIGINS, _cfg_issues = _cr.resolve(_repo, [("project", _project_layer)])
config["repo_dir"] = str(_repo)          # injected by submit.sh, not user-authored

if _cfg_issues:
    raise RuntimeError(
        "Configuration problems (preflight reports the same findings):\n  "
        + "\n  ".join(f"{code}: {detail}" for code, detail in _cfg_issues))

# --- Load sample sheet ---------------------------------------------------
samples = pd.read_csv(
    config["samples"]["sheet"],
    sep="\t",
    comment="#",
).set_index("sample_id", drop=False)

SAMPLES = samples.index.tolist()
N_SAMPLES = len(SAMPLES)

# --- Where the workflow itself lives -------------------------------------
# Shell rules invoke helper scripts by path. Those live with the workflow, not
# in the project, and the working directory is the project, so a repo-relative
# path like "workflow/scripts/01_qc_quant.sh" would not resolve.
# submit.sh passes --config repo_dir=<checkout>; the "." default keeps a bare
# `snakemake` run from inside the checkout working as before.
# Passed as config rather than read from a Snakemake attribute so the value is
# explicit and survives into SLURM job invocations.
REPO_DIR = Path(config.get("repo_dir", ".")).resolve()

# --- Output roots --------------------------------------------------------
# All relative, so they land in whatever directory the run was launched from.
# submit.sh chdirs into the project first, which is what keeps two projects
# using this one checkout out of each other's state.
OUT = Path(config["project"]["output_dir"])
QUANT = OUT / "quant"
METRICS = Path("metrics")
GATES = Path("gates")

# Reference artifacts (transcriptome, GTF, salmon index) are derived purely from
# the assembly, so they can be shared between projects. Keying the cache by
# accession means two projects on the same genome reuse one index, while two on
# different genomes cannot overwrite each other. Set reference.cache_dir to null
# to keep them inside the project instead.
# Note: reference.annotation_tsv.path stays project-relative and is unaffected.
#
# ponytail: no lock on the shared cache. Two runs starting a build of the SAME
# accession at the same moment can race, because Snakemake's lock covers a
# working directory, not this cache. Add an flock around the index rules if that
# becomes real; sequential runs and different accessions are already safe.
_ref_cache = (config.get("reference") or {}).get("cache_dir")
_accession = (config.get("reference") or {}).get("accession") or "unspecified_assembly"
if _ref_cache:
    REF = Path(os.path.expanduser(str(_ref_cache))) / str(_accession)
else:
    REF = Path("reference")

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

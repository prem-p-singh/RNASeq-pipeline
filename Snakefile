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
import hashlib
import json
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

_cfg_blocking = _cr.blocking(_cfg_issues, _repo)
for _code, _detail in _cfg_issues:
    if (_code, _detail) not in _cfg_blocking:
        print(f"config warning {_code}: {_detail}")
if _cfg_blocking:
    raise RuntimeError(
        "Configuration problems (preflight reports the same findings):\n  "
        + "\n  ".join(f"{code}: {detail}" for code, detail in _cfg_blocking))

from recommend import recommend as _recommend, requires_inference
from preprocessing import policy_args
from screen_reads import validate as validate_screen, index_files
from small_rna import load_kit, kit_policy
import tagseq as _tagseq
RECOMMENDATION = _recommend(config)
if RECOMMENDATION["issues"]:
    raise RuntimeError("Unsupported analysis: " + "; ".join(RECOMMENDATION["issues"]))
# Small-RNA and named TAG-seq kits fix adapter, trimming and UMI settings; the kit table wins
# over the generic policy (a TAG-seq config that sets a kit key to another value is refused).
TAGSEQ_KIT = _tagseq.load_kit(config["analysis"]["kit"]) if RECOMMENDATION["route"] == "star_tagseq" else None
UMI_ROUTE = RECOMMENDATION["count_treatment"] == "umi_gene_counts"
PREPROCESSING = config["preprocessing"]
if config["samples"]["seq_type"] == "srnaseq":
    PREPROCESSING = PREPROCESSING | kit_policy(load_kit(config["analysis"]["small_rna_kit"]))
if TAGSEQ_KIT:
    PREPROCESSING = _tagseq.merge_policy(PREPROCESSING, TAGSEQ_KIT)
policy_args(PREPROCESSING, config["samples"]["seq_type"] == "rnaseq_paired")

SCREEN_POLICY = validate_screen(config["screening"])
# FASTA panel members (transcriptome, genome, path or URL) are indexed in the
# project by build_screen_index; the run policy then names only Bowtie2 indexes.
SCREEN_BUILT = {r["name"]: r["fasta"] for r in SCREEN_POLICY["references"] if "fasta" in r} if SCREEN_POLICY["enabled"] else {}
SCREEN_INDEXES = [str(f) for ref in SCREEN_POLICY["references"] if "index" in ref for f in index_files(ref["index"])] if SCREEN_POLICY["enabled"] else []
SCREEN_INDEXES += [f"screening_panel/{name}/built.txt" for name in SCREEN_BUILT]
SCREEN_RUN_POLICY = dict(SCREEN_POLICY, references=[
    dict(name=r["name"], role=r["role"], index=str(Path("screening_panel", r["name"], r["name"]).resolve()))
    if "fasta" in r else r for r in SCREEN_POLICY["references"]])

QC_ONLY = set(config["analysis"]["objectives"]) == {"qc"}
STAR_ROUTE = RECOMMENDATION["route"] in ("star_counts", "star_tagseq")
SMALL_RNA_ROUTE = RECOMMENDATION["route"] in ("small_rna_animal", "small_rna_plant")
PLANT_SRNA_ROUTE = RECOMMENDATION["route"] == "small_rna_plant"
# Plant counts are one cohort matrix (ShortStack clusters loci jointly); the per-sample product is its QC record.
QUANT_PRODUCT = ("star/umi_gene_counts.tsv" if UMI_ROUTE else "star/gene_counts.tsv" if STAR_ROUTE else "small_rna/metrics.json" if PLANT_SRNA_ROUTE
                 else "small_rna/counts.tsv" if SMALL_RNA_ROUTE else "quant.sf")

# --- Load sample sheet ---------------------------------------------------
import metadata as _metadata
CANONICAL_INPUTS = _metadata.execution_inputs(Path.cwd(), config)
if STAR_ROUTE and not TAGSEQ_KIT and CANONICAL_INPUTS is not None:
    for sample_id, library in CANONICAL_INPUTS.items():
        if library['expected_libtype'] not in ('U', 'SF', 'SR', 'IU', 'ISF', 'ISR'):
            raise ValueError(f'STAR requires known strandedness for {sample_id}')
samples = pd.read_csv(
    config["samples"]["sheet"],
    sep="\t",
    comment="#",
    converters={"sample_id": str},
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

# Shared references are keyed by URLs, index parameters and the pinned Salmon
# version. The completion lock records source checksums and is verified on reuse.
import reference_cache as _refcache

_ref_cfg = config.get("reference") or {}
_ref_cache = _ref_cfg.get("cache_dir")
# Index construction parameters. Kept here because they must be knowable when
# the DAG is built: they decide the cache path.
SALMON_KMER = int(_ref_cfg.get("kmer") or 31)
SALMON_DECOYS = _ref_cfg.get("decoys") or None

# RF09 and master plan 5.3: "Assembly accession alone is insufficient:
# annotation, decoys, tool version, and index parameters can differ." The key
# was the bare accession, so two projects agreeing on the genome but differing
# in annotation release, k-mer or decoy status shared one directory and
# consumed each other's index. The digest below covers all of them.
import star_counts as _star
if QC_ONLY:
    CACHE_KEY = "unused-qc"
elif STAR_ROUTE:
    CACHE_KEY = _star.cache_key(_ref_cfg)
elif SMALL_RNA_ROUTE:
    CACHE_KEY = "srna-" + hashlib.sha256(json.dumps({
        "mature_mirna_url": _ref_cfg.get("mature_mirna_url"),
        "genome_fasta_url": _ref_cfg.get("genome_fasta_url"),
        "mirbase_prefix": config["organism"].get("mirbase_prefix")}, sort_keys=True).encode()).hexdigest()[:16]
else:
    CACHE_KEY = _refcache.cache_key(_ref_cfg, SALMON_KMER, SALMON_DECOYS)
if _ref_cache:
    REF = Path(os.path.expanduser(str(_ref_cache))) / CACHE_KEY
else:
    REF = Path("reference")

# Refuse an entry that is incomplete or was built from other parameters,
# rather than reading whatever happens to sit at that path.
if not QC_ONLY and not STAR_ROUTE and not SMALL_RNA_ROUTE and _ref_cache and REF.exists():
    _entry_issues = _refcache.verify_entry(REF, _ref_cfg, SALMON_KMER, SALMON_DECOYS)
    _fatal = [i for i in _entry_issues if not i.startswith(("missing ", "empty "))]
    if _fatal:
        raise RuntimeError(
            "Reference cache entry does not match this run (RF09):\n  "
            + str(REF) + "\n  " + "\n  ".join(_fatal))
    if _entry_issues:
        print("Reference cache entry incomplete; it will be rebuilt:")
        for i in _entry_issues:
            print(f"  {i}")

if STAR_ROUTE and (REF / "star_idx/reference.lock.json").is_file():
    _star.verify_index(REF / "star_idx", REF / "genome.fa", REF / "annotation.gtf")

for d in (OUT, QUANT, METRICS, GATES) + (() if QC_ONLY else (REF,)):
    d.mkdir(parents=True, exist_ok=True)

_snapshot = json.dumps(config, sort_keys=True, indent=2) + "\n"
_snapshot_path = GATES / "resolved_config.json"
if not _snapshot_path.exists() or _snapshot_path.read_text() != _snapshot:
    _snapshot_path.write_text(_snapshot)
_lock_path = REPO_DIR / "environments/linux-64.explicit.txt"
RUNTIME_ID = hashlib.sha256(_lock_path.read_bytes()).hexdigest() if _lock_path.exists() else "unlocked"
if STAR_ROUTE:
    RUNTIME_ID += hashlib.sha256((REPO_DIR / "environments/star-linux-64.explicit.txt").read_bytes()).hexdigest()
if SMALL_RNA_ROUTE:
    RUNTIME_ID += hashlib.sha256((REPO_DIR / "environments/srna-linux-64.explicit.txt").read_bytes()).hexdigest()
if UMI_ROUTE:
    RUNTIME_ID += hashlib.sha256((REPO_DIR / "environments/umi-linux-64.explicit.txt").read_bytes()).hexdigest()


if SCREEN_POLICY["enabled"]:
    RUNTIME_ID += hashlib.sha256((REPO_DIR / "environments/screen-linux-64.explicit.txt").read_bytes()).hexdigest()

# --- Invalidate results a flag claims but the filesystem lacks ----------
# R09 / master plan 11: "A done flag alone never proves a result exists."
# Runs here, at parse time, so it happens BEFORE Snakemake resolves the DAG:
# clearing a stale flag and manifest makes the producing stage look incomplete,
# and it is rescheduled in this same invocation rather than the next one.
import artifacts as _artifacts

_stale = _artifacts.invalidate_stale(OUT)
if _stale:
    print("Stale results detected; affected stages will rerun:")
    print(_artifacts.format_report(_stale))

# --- Rule modules --------------------------------------------------------
include: "workflow/rules/common.smk"
include: "workflow/rules/retrieve.smk"
include: "workflow/rules/build_orgdb.smk"
include: "workflow/rules/per_sample.smk"
include: "workflow/rules/star.smk"
include: "workflow/rules/small_rna.smk"
include: "workflow/rules/qc_report.smk"
include: "workflow/rules/aggregate.smk"

# --- Conditional target helpers -----------------------------------------
def all_targets():
    """Build the final target list. OrgDb is only required when enrichment runs."""
    early = [OUT / "preprocessing_report/multiqc_report.html",
             *expand(str(QUANT / "{sample}/fastp.html"), sample=SAMPLES)]
    early += expand(str(QUANT / "{sample}/screening/screening.{ext}"), sample=SAMPLES, ext=["json", "tsv", "html"])
    if QC_ONLY:
        return early
    targets = early + [
        # Stage 1 — one quant.sf per sample
        *expand(str(QUANT / "{sample}" / QUANT_PRODUCT), sample=SAMPLES),
        *(expand(str(QUANT / "{sample}/read_preparation.json"), sample=SAMPLES) if CANONICAL_INPUTS is not None else []),
        # Stage 2 — aggregated counts
        OUT / "counts.tsv",
    ]
    # Stage 1b — comparative QC report (before/after-clean charts + MultiQC).
    # ponytail: the report reads Salmon/STAR logs only; small-RNA runs rely on
    # the preprocessing report and small_rna/metrics.json until it learns Bowtie.
    if not SMALL_RNA_ROUTE:
        targets += [OUT / "qc_report/qc_report.done",
                    OUT / "qc_report/qc_charts.html",
                    OUT / "qc_report/comparative_charts.png",
                    OUT / "qc_report/alignment_summary.tsv",
                    OUT / "qc_report/multiqc/multiqc_report.html"]
    if STAR_ROUTE:
        targets.extend(expand(str(QUANT / "{sample}/star/Aligned.sortedByCoord.out.bam{suffix}"),
                              sample=SAMPLES, suffix=["", ".bai"]))
    if requires_inference(config):
        targets.extend([OUT / "de_done.flag", OUT / "enrichment_done.flag", OUT / "wgcna_done.flag"])
    elif "coexpression" in config["analysis"]["objectives"]:
        targets.append(OUT / "wgcna_done.flag")
    # Only ask for the OrgDb when GO enrichment will actually consume it.
    if requires_inference(config) and orgdb_required():
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

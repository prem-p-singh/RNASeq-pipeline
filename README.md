<div align="center">
<img src="./assets/readme-banner.svg" alt="RNA-Seq pipeline: preflight, reference, quantify, differential expression, enrichment and WGCNA" width="100%" />
</div>

# RNASeq pipeline

A configuration-driven Snakemake workflow for annotated bulk RNA-seq gene expression on Linux x86_64. It checks study design before quantification, records sample exclusions and reference provenance, and produces differential-expression tables and QC reports. GO/KEGG enrichment and WGCNA are optional.

[![Linux release validation](https://github.com/prem-p-singh/RNASeq-pipeline/actions/workflows/release-validation.yml/badge.svg)](https://github.com/prem-p-singh/RNASeq-pipeline/actions/workflows/release-validation.yml)
[![Lightweight checks](https://github.com/prem-p-singh/RNASeq-pipeline/actions/workflows/checks.yml/badge.svg)](https://github.com/prem-p-singh/RNASeq-pipeline/actions/workflows/checks.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-285F47)](LICENSE)

## Release scope

This is a release candidate for the implemented bulk workflow. The broader multi-assay RNA-seq platform remains a development roadmap. See [release qualification and limitations](docs/RELEASE.md) before using results for a publication.

| Input/design | Implemented route | Qualification |
|---|---|---|
| Non-UMI bulk, paired-end or single-end | fastp → Salmon → tximport `lengthScaledTPM` → limma-voom | Linux end-to-end release tests |
| Bulk with a declared random effect | Same measurement route → dream | Provisional; design checks exist, but no end-to-end mixed-model qualification |
| 3′ TAG-seq, non-UMI | fastp → Salmon → counts without length correction | Provisional; no named kit qualified |
| GO enrichment | Signed-statistic GSEA with a compatible OrgDb | Local annotation fixture; verify mapping for your organism |
| KEGG enrichment | Signed-statistic GSEA | Requires current network data; live result qualification remains separate |
| WGCNA | Filtered log-CPM, declared parameters, memory-aware blocks | Synthetic module-recovery tests; no large-cohort memory guarantee |
| UMI, single-cell, spatial, small RNA, long reads, splicing, fusions, variants, external count import | Not implemented | Rejected rather than substituted with bulk analysis |

One sample-sheet row represents one biological measurement and one FASTQ or read pair. Multiple lanes must currently be merged externally with a recorded method. The separate metadata-table validators do not yet make multi-library execution available.

Salmon is the implemented quantifier, not a universal recommendation. [Executable recommendation rules](config/recommendation_rules.yaml) explain the supported selection and record rule IDs in `gates/recommendation.json`. The declared dependence structure selects limma-voom or dream. Increasing sample count does not silently change the scientific model.

## Environment

Use **Linux x86_64**, Bash, Conda, `curl`, `gzip` and `flock`. FARM supplies these prerequisites; macOS/ARM is not a qualified runtime. Analysis runs on allocated compute nodes, not the login node.

The [explicit lock](environments/linux-64.explicit.txt) records 579 exact package builds and SHA256 hashes. It includes Python 3.12.14, Snakemake 9.19.0, R 4.5.2, Salmon 1.10.3, fastp 0.23.4 and the analysis libraries. `environment.yml` is a maintainer solve specification, not the installation lock.

```bash
# FARM only; elsewhere make your existing Conda installation available.
module load conda/latest

# Choose a location with enough space, visible on every worker.
export RNASEQ_ENV_PREFIX=/absolute/shared/path/rnaseq-runtime
# Optional: place downloaded Conda packages on an appropriate filesystem.
export CONDA_PKGS_DIRS=/absolute/path/conda-packages

prefix=$(bash scripts/bootstrap.sh)
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$prefix"
```

Bootstrap downloads missing packages into a new environment, checks exact installed package metadata, imports the required Python/R libraries and checks executable versions and locations. A mismatched existing environment is refused; create a separate prefix to rebuild it. It does not install Conda or alter the operating system. Normal launches repeat verification. Annotation packages built for a project are kept in the configured OrgDb cache, outside the release runtime.

The environment, Conda download cache, references and dataset all need space. Choose paths according to actual filesystem capacity and quota. There is **no default 20 GB limit**.

## Create a project

Keep projects outside the source checkout. For reproducible reference selection, edit an explicit project configuration:

```bash
export REPO=/absolute/path/RNASeq_pipeline
export PROJECT=/absolute/shared/path/my_study
mkdir -p "$PROJECT/config"
cp "$REPO/config/config.template.yaml" "$PROJECT/config/config.yaml"
cp "$REPO/config/thresholds.yaml" "$PROJECT/config/thresholds.yaml"
cp "$REPO/config/samples.tsv.template" "$PROJECT/config/samples.tsv"
```

Edit all three files. Supply versioned transcriptome FASTA and matching GTF URLs, sample FASTQ paths/URLs, organism identity, assay/layout, model and contrasts. Use `rnaseq_single` or `rnaseq_paired` for ordinary bulk data. Review the default downstream settings; disable analyses you do not want. Review QC thresholds for your experiment rather than copying the deliberately permissive smoke-test thresholds.

Alternatively, copy [setup_inputs.template.yaml](scripts/setup_inputs.template.yaml) into your project, edit it, and run:

```bash
python3 "$REPO/scripts/setup.py" --project-dir "$PROJECT" "$PROJECT/setup_inputs.yaml"
```

Setup can resolve references from external services. Review the generated URLs and metadata before launch; a service's “latest” assembly is not a permanent reference identity. The interactive FARM helper is `bash scripts/new_project.sh PROJECT_NAME` and expects files in `~/new_project_inbox/PROJECT_NAME`.

## Plan, launch and resume

```bash
# Preflight and storage plan only, using an already available Python/R runtime:
bash "$REPO/submit.sh" -d "$PROJECT" --plan-only

# Verify/bootstrap the runtime and inspect the complete DAG:
bash "$REPO/submit.sh" -d "$PROJECT" --dry-run

# FARM/SLURM execution; use the same command to resume:
bash "$REPO/submit.sh" -d "$PROJECT"
```

The shipped profiles contain **UC Davis FARM** account, partition and QoS settings. Adapt those values for another cluster. `-p small|medium|large` selects a profile; per-rule resources and storage planning constrain its concurrency. Worker nodes must see the project, source checkout, environment, input and reference paths.

For a local Linux worker allocation, after bootstrapping and activating the runtime:

```bash
cd "$PROJECT"
python3 "$REPO/scripts/preflight.py" -d "$PROJECT"
python3 "$REPO/scripts/plan_resources.py" -d "$PROJECT" --out gates/resource_plan.json
snakemake --snakefile "$REPO/Snakefile" --configfile config/config.yaml \
  --config repo_dir="$REPO" --cores 4 --scheduler greedy
```

`--plan-only` is the read-mostly planning entry point. Snakemake DAG construction also writes resolved configuration and invalidates stale completion markers, so its dry run can repair bookkeeping even though it does not execute analysis rules.

The storage plan measures available input sizes, estimates unknown remote objects conservatively, accounts for retained downloads, references, environment and concurrent temporary files, then reduces concurrency if needed. More samples increase retained output; at fixed concurrency scratch peak need not grow with cohort size. Estimates are assumptions, not benchmark guarantees. Filesystem free space does not always reveal an account quota: provide a real remaining quota through the documented legacy `hpc.storage_budget_gb` setting when required. That optional constraint has no default. `hpc.samples_in_flight` can lower concurrency further.

## Outputs and provenance

```text
project/
  config/                         inputs, model and thresholds
  gates/                          preflight, recommendation, environment,
                                  resolved config, resource plan and decisions
  reference/                      local reference, index and reference.lock.json
  results/
    quant/<sample>/                quant.sf, fastp reports, metrics and input hashes
    qc_report/                     comparative HTML/PNG/TSV and MultiQC
    counts.tsv                     retained sample gene counts
    gene_lengths.tsv               effective lengths for retained samples
    gene_import.rds                 retained tximport data
    counts_provenance.json          count origin and length policy
    sample_disposition.tsv          inclusion/exclusion and reasons
    DE_Results/                    contrast-specific results
    de_manifest.tsv                authoritative contrast/output list
    GO_results/, KEGG_results/     optional enrichment tables
    enrichment_status.tsv          success, empty, skipped, unavailable or failed
    WGCNA/, wgcna_manifest.tsv      optional network outputs and status
  metrics/                        machine-readable stage summaries
  logs/                           per-rule logs
```

Shared reference caches use a digest of source URLs and index construction parameters, including the pinned Salmon version. The reference lock records actual source checksums and identifier compatibility. Use immutable/versioned URLs: a changed remote object at the same URL is not automatically discovered. Local input changes, analysis configuration and runtime-lock changes participate in rerun decisions. Missing manifest-listed results are rebuilt. Preserve the whole project and the exact source revision with your results.

User-supplied local FASTQs are never deleted. `hpc.delete_fastq_after_quant` controls cleanup of downloads and trimmed reads created by the workflow after validated quantification. Failed work retains its intermediates for diagnosis.

## Validation

```bash
# Full runtime required; any skip is a release failure.
bash tests/run_all.sh --strict
# Networked public-data smoke test, separate from the offline fixtures.
bash tests/check_public.sh
# FARM: clean environment installation and scientific tests on a worker.
sbatch scripts/validate_farm.sbatch
```

The strict suite covers configuration/design gates, count semantics, WGCNA behavior, reference/cache safety, output recovery, isolated projects, raw single/paired FASTQ processing and GO enrichment. The public fixture uses six biological samples from a pinned, downsampled GSE110004 yeast dataset. It tests workflow operation, not reproduction of the paper's genome-wide conclusions. See [release evidence](docs/RELEASE.md).

Developer mode (`bash tests/run_all.sh`) permits clearly reported dependency skips. The lightweight CI job is not the scientific release gate. Require the separate `release-validation` job before tagging a release.

## License and citation

Pipeline source is [MIT licensed](LICENSE). External tools and reference/data sources retain their own licenses. The sequencing handbook PDF and private planning documents are not redistributed. Citation metadata is in [CITATION.cff](CITATION.cff).

Maintained by Prem Pratap Singh, Department of Viticulture and Enology, UC Davis.

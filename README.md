<div align="center">

<img src="./assets/readme-banner.svg" alt="RNA-Seq pipeline: preflight, reference, quantification, differential expression, enrichment and co-expression" width="100%" />

<br />

[![self-checks](https://github.com/prem-p-singh/RNASeq-pipeline/actions/workflows/checks.yml/badge.svg)](https://github.com/prem-p-singh/RNASeq-pipeline/actions/workflows/checks.yml)
[![License](https://img.shields.io/badge/license-MIT-285F47.svg)](LICENSE)
[![Snakemake](https://img.shields.io/badge/snakemake-9.19-397C68.svg)](https://snakemake.readthedocs.io)
[![Python](https://img.shields.io/badge/python-3.12-46567D.svg)](https://www.python.org)
[![R](https://img.shields.io/badge/R-Bioconductor-9C482D.svg)](https://bioconductor.org)
[![Status](https://img.shields.io/badge/status-pre--release-7A5B16.svg)](#development-status)

</div>

## Overview

A configuration-driven Snakemake workflow for bulk RNA-seq and 3′ TAGseq differential
expression on SLURM clusters. Species, experimental design and downstream analyses are
set in configuration; no code is modified between projects.

The workflow validates the study design before committing compute. Model estimability,
metadata integrity and the requested assay/design/backend combination are checked at
preflight, so a design that cannot support the requested comparisons fails in seconds
rather than after the cohort has been quantified.

**Pipeline stages**

```
preflight → reference preparation → quantification → sample disposition
          → differential expression → enrichment / co-expression → report
```

## Development status

This is pre-release software. Stages 0 and 1 have been executed on real data; stages 2
through 5 have not been run end to end.

| Stage | Implementation | Executed on real data |
|---|---|---|
| 0. Reference retrieval and indexing | Complete | Yes |
| 1. Read QC and quantification | Complete | Yes |
| 2. Aggregation to gene counts | Complete | Not yet |
| 3. Differential expression | Complete | Not yet |
| 4. Functional enrichment | Complete | Not yet |
| 5. Co-expression networks | Complete | Not yet |

All five stages are implemented as parameterised refactors of the analyses in
[`examples/grape/`](examples/grape/). The blocking dependency is the runtime: the target
cluster's R module does not provide tximport, variancePartition, clusterProfiler or
WGCNA, so stages 2 through 5 require the conda environment defined in
[`environment.yml`](environment.yml). That environment has not yet been built, and
`environment.lock.yml` therefore does not exist. Downstream results should not be
treated as validated until both are done.

Quantitative agreement with an independently scripted reference analysis has not been
established. See [Validation](#validation).

## Features

**Design validation before execution.** Preflight checks configuration keys, sample
sheet integrity, scope, and model estimability (rank deficiency, residual degrees of
freedom, independent biological replicates), emitting 23 stable issue codes with
severity and remedy. It calls the same `workflow/scripts/_design.R` used by the
differential expression stage, so the two cannot disagree.

**Declarative contrasts.** Comparisons are specified as structured records and built
through `emmeans`, not by evaluating configuration strings as code. Contrast directions
are recorded alongside results.

**Assay-aware count handling.** Transcript abundance is length-corrected for bulk
RNA-seq and left uncorrected for 3′ tag counts, which do not scale with transcript
length.

**Explicit sample disposition.** Inclusion and exclusion decisions are recorded per
sample with the reason and the governing threshold, and excluded samples are dropped
before the count matrix is written.

**Truthful stage status.** Each analysis reports one of succeeded, empty, failed,
skipped by configuration, or unavailable prerequisites. A zero-discovery result is a
successful analysis; a missing annotation makes enrichment unavailable rather than
silently empty.

**Bounded resource use.** The launcher derives the peak per-sample working set from the
assay, caps concurrency to remain within a declared storage budget, and refuses runs
that cannot fit. Intermediates the pipeline created are released after quantification is
verified; source FASTQ files supplied by the user are never modified or deleted.

**Project isolation.** Each dataset runs in its own directory against a read-only
checkout. Reference bundles are cached by assembly accession and shared across projects
on the same genome.

## Scope

Supported combinations are declared in [`config/spec.yaml`](config/spec.yaml) and
enforced at preflight.

| Capability | Status | Conditions |
|---|---|---|
| Paired-end bulk RNA-seq | Supported | Both mates present for every library |
| Single-end bulk RNA-seq | Supported | Length-corrected abundance |
| 3′ TAGseq | Provisional | Counts not length-corrected; no kit validated to date |
| Independent groups, factorial designs | Supported | limma-voom; contrasts must be estimable |
| Repeated measures | Supported | dream; replication counted in subjects, not rows |
| Non-model organisms | Conditional | Requires compatible annotation and identifier coverage |
| Functional enrichment | Optional | Requires an OrgDb or a versioned gene-set table |
| Co-expression networks | Optional | Requires sufficient samples after filtering |
| Small RNA-seq, single-cell | Not supported | Distinct preprocessing and data models |
| UMI protocols | Not supported | No validated extraction or deduplication step |
| Multi-lane libraries | Not supported | Lanes must be merged before the workflow runs |
| External quantification import | Not supported | Reference and provenance validation not implemented |

## Requirements

- SLURM cluster, or a local machine for small validation runs
- conda or mamba
- Python 3.12, Snakemake 9.19, Salmon 1.10.3, fastp 0.23.4
- R with Bioconductor: tximport, edgeR, limma, variancePartition, clusterProfiler, WGCNA

All runtime dependencies are pinned in [`environment.yml`](environment.yml).

## Installation

```bash
git clone https://github.com/prem-p-singh/RNASeq-pipeline.git
cd RNASeq-pipeline
conda env create -f environment.yml
conda activate rnaseq-pipeline
```

Record the resolved environment before running an analysis intended for publication:

```bash
conda env export --no-builds > environment.lock.yml
```

## Usage

Each project occupies its own directory. The checkout is read-only at run time, so
multiple projects share a single installation without interfering.

```bash
PROJ=~/rnaseq_projects/my_study
mkdir -p "$PROJ"

# Configuration from seven inputs
cp scripts/setup_inputs.template.yaml "$PROJ/setup_inputs.yaml"
$EDITOR "$PROJ/setup_inputs.yaml"
python3 scripts/setup.py --project-dir "$PROJ" "$PROJ/setup_inputs.yaml"

# Validate before committing compute
python3 scripts/preflight.py -d "$PROJ"

# Execute; resumes from the last completed step if interrupted
./submit.sh -d "$PROJ"
```

`submit.sh` runs preflight itself and will not submit if it reports an error.

<details>
<summary><strong>Guided setup</strong></summary>

For interactive configuration on the cluster, upload a metadata spreadsheet and run:

```bash
scripts/start_new.sh <project_name> <metadata_file> [fastq_dir]
```

Sample-ID and factor columns are detected from the spreadsheet; remaining fields are
prompted. The project is created under `~/rnaseq_projects/`, overridable with
`RNASEQ_PROJECTS_ROOT`.

</details>

<details>
<summary><strong>Preflight output</strong></summary>

```text
assay=tagseq  design=independent_groups  backend=limma_voom  count_treatment=no
samples=6  model_vars=treatment  unit=sample
design: 2 coefficients, 4 residual df, 3 biological replicates per sample
stages: reference, quantification, aggregation, differential_expression, enrichment, report
  not planned: wgcna (downstream.run_wgcna is false)
```

Written to `gates/preflight_issues.tsv` (code, severity, scope, message, detail,
remedy) and `gates/preflight_plan.json` (sample and design counts, resolved settings,
and the planned status of every stage with a reason for each omission).

</details>

<details>
<summary><strong>Project layout</strong></summary>

```text
~/rnaseq_projects/<name>/
  config/      config.yaml, samples.tsv, thresholds.yaml
  results/     counts, DE_Results, Enrichment, WGCNA
  gates/       preflight_issues.tsv, preflight_plan.json, decisions.log
  metrics/     per-stage JSON consumed by subsequent stages
  logs/        per-rule logs

~/rnaseq_reference_cache/<accession>/
               transcriptome, GTF and Salmon index, shared across projects
               on the same assembly
```

</details>

## Repository layout

```text
Snakefile                 Workflow entry point
submit.sh                 Launcher: preflight, tier selection, storage cap
scripts/
  preflight.py            Pre-execution validation
  setup.py                Project configuration generator
  new_project.sh          Interactive setup
  presets/                Curated organism metadata by NCBI taxID
config/
  spec.yaml               Supported assays, designs, backends; issue codes
  config.template.yaml    Configuration reference and validation schema
  thresholds.yaml         Decision-gate cutoffs
workflow/
  rules/                  Snakemake modules, one per stage
  scripts/                Stage implementations and shared design helpers
tests/                    Self-checks
examples/grape/           Source analyses this workflow generalises
```

## Documentation

| Document | Contents |
|---|---|
| [DESIGN.md](DESIGN.md) | Architecture, stage boundaries, decision gates, storage model |
| [INPUTS.md](INPUTS.md) | The seven setup inputs and every file the wizard generates |
| [STRATEGY.md](STRATEGY.md) | Tier selection and concurrency behaviour |
| [HANDBOOK_INSIGHTS.md](HANDBOOK_INSIGHTS.md) | Methodological choices, with the reasoning and sources behind them |
| [config/config.template.yaml](config/config.template.yaml) | Every configuration key, annotated; also the validation schema |
| [config/thresholds.yaml](config/thresholds.yaml) | Default decision-gate cutoffs |
| [config/spec.yaml](config/spec.yaml) | Supported scope and the issue-code catalogue |

## Validation

The repository contains seven self-checks. There is no test framework; each file is a
script that exits non-zero on failure.

```bash
bash tests/run_all.sh
```

| Check | Covers |
|---|---|
| `check_de.R` | Contrast construction and direction, design estimability, sample disposition |
| `check_wgcna.R` | Block sizing against the memory allocation, recorded parameters |
| `check_setup_helpers.py` | Reference URL resolution, mate detection, sample-to-file matching |
| `check_qc_report.py` | Missing versus zero metrics, exclusion of stale sample directories |
| `check_preflight.py` | Issue codes, scope enforcement, estimability gate |
| `check_storage_policy.sh` | Intermediate cleanup, preservation of source reads |
| `check_project_isolation.sh` | Independence of concurrent project state |

Checks whose dependencies are absent report `SKIP` with the reason and are counted
separately from passes. Continuous integration installs base R but not the Bioconductor
stack, so `check_de.R` and `check_wgcna.R` are skipped there and the skip count is
reported.

These checks verify implementation behaviour. They do not establish numerical agreement
with an independent reference analysis, which remains outstanding.

## Citation

If this workflow contributes to published work, please cite the repository. Citation
metadata is provided in [`CITATION.cff`](CITATION.cff); GitHub renders it under
**Cite this repository**.

## License

Released under the MIT License. See [LICENSE](LICENSE).

## Author

**Prem Pratap Singh**, Postdoctoral Scholar, Department of Viticulture and Enology,
University of California, Davis.

[Website](https://www.prempsingh.com) ·
[Google Scholar](https://scholar.google.com/citations?user=UGFMZEYAAAAJ&hl=en) ·
[ORCID](https://orcid.org/0000-0001-7921-9379) ·
[LinkedIn](https://www.linkedin.com/in/prem-p-singh)

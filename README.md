# RNA-Seq Pipeline

A generic, HPC-aware, adaptive RNA-Seq analysis pipeline. Runs on any dataset with a sample sheet and a config file — no code changes needed to swap species, experimental design, or project.

Built for a **20 GB storage budget** on a SLURM cluster, with **storage-aware concurrency** (the launcher caps how many samples are in flight so the working set stays inside the budget) and **adaptive downstream decisions** (DE backend, batch correction, enrichment filtering, WGCNA skip) driven by metrics from prior stages.

---

## Quick start

Each project gets its own directory. The pipeline checkout is read-only at run
time, so several projects can share one copy without overwriting each other.

**Guided (recommended).** From your Mac, upload a metadata spreadsheet and start
the wizard on FARM:

```bash
scripts/start_new.sh <project_name> <metadata_file> [fastq_dir]
```

That uploads into `~/new_project_inbox/<project_name>/` and drops you into
`scripts/new_project.sh`, which detects the sample-ID and factor columns, asks a
few questions, then runs setup and launches. The project lands in
`~/rnaseq_projects/<project_name>/` (override with `RNASEQ_PROJECTS_ROOT`).

**By hand**, if you would rather drive it yourself:

```bash
PROJ=~/rnaseq_projects/my_study
mkdir -p "$PROJ"

# 1. Fill in 7 minimal fields
cp scripts/setup_inputs.template.yaml "$PROJ/setup_inputs.yaml"
$EDITOR "$PROJ/setup_inputs.yaml"     # project name, taxID, fastq source, metadata file, ...

# 2. Setup wizard: writes config.yaml, samples.tsv, thresholds.yaml and the
#    NCBI annotation table into the project, and resolves the OrgDb
python3 scripts/setup.py --project-dir "$PROJ" "$PROJ/setup_inputs.yaml"

# 3. Check the project before spending any compute on it (optional: submit.sh
#    runs this itself and refuses to submit if it finds errors)
python3 scripts/preflight.py -d "$PROJ"

# 4. Launch. The submitter counts samples, picks the SLURM tier, and caps
#    concurrency to stay inside the storage budget
./submit.sh -d "$PROJ"

# 5. Resume after any interruption: same command
./submit.sh -d "$PROJ"
```

The first run for a given genome fetches the reference and builds a Salmon index.
Those are cached in `~/rnaseq_reference_cache/<accession>/` and shared by every
project on that assembly, so it happens once, not once per project.
See [**INPUTS.md**](INPUTS.md) for the full list of what the wizard generates.

`intake_template.xlsx` is a **planning aid**: it lists every decision you will be
asked about so you can work them out beforehand. Nothing parses it; you type the
answers into the prompts. `rnaseq_pipeline.n8n.json` is likewise a **diagram** of
the stages for reading in n8n, not an executable workflow.

---

## Where to read next

| Doc | Read when you want to… |
|---|---|
| [**INPUTS.md**](INPUTS.md) | Use the setup wizard — 7 inputs, auto-generates everything external |
| [**DESIGN.md**](DESIGN.md) | Understand the architecture, stages, decision gates, and storage budget |
| [**STRATEGY.md**](STRATEGY.md) | Understand how `./submit.sh -d <project>` auto-picks small / medium / large SLURM tier based on sample count |
| [`config/config.template.yaml`](config/config.template.yaml) | See every knob you can set per project |
| [`config/thresholds.yaml`](config/thresholds.yaml) | See default decision-gate cutoffs (mapping rate, WGCNA power, etc.) |
| [`examples/grape/`](examples/grape/) | The worked example — original grape/GRBV project that this pipeline was refactored from |

---

## Directory map

```
RNASeq_pipeline/
├── README.md                ← you are here
├── DESIGN.md                ← architecture + decision gates
├── STRATEGY.md              ← three-tier auto-selection
│
├── Snakefile                ← Snakemake entry point
├── submit.sh                ← launcher (counts samples, picks tier)
│
├── scripts/
│   ├── preflight.py         ← validates a project before any compute
│   ├── setup.py             ← generates a project's config from 7 inputs
│   └── new_project.sh       ← guided setup on FARM
│
├── config/
│   ├── spec.yaml            ← supported assays/designs/backends + issue codes
│   ├── config.template.yaml ← also the schema: preflight rejects unknown keys
│   ├── thresholds.yaml
│   └── samples.tsv.template
│
├── profiles/
│   ├── small/               ← ≤ 20 samples, serial
│   ├── medium/              ← 21–200 samples, 20 concurrent
│   └── large/               ← 201+ samples, throttled
│
├── workflow/
│   ├── rules/
│   │   ├── common.smk       ← helpers (logging, metrics)
│   │   ├── retrieve.smk     ← Stage 0: fetch ref + build Salmon index
│   │   ├── per_sample.smk   ← Stage 1: one SLURM task per sample
│   │   └── aggregate.smk    ← Stages 2–5: DE, enrichment, WGCNA
│   └── scripts/
│       ├── 01_qc_quant.sh           ← fastp + salmon, releases intermediates
│       ├── 02_aggregate.R           ← tximport + sample disposition
│       ├── 03_de.R                  ← limma-voom or dream, declarative contrasts
│       ├── 04_enrichment.R          ← GO + KEGG, per-contrast status
│       └── 05_wgcna.R               ← WGCNA, blocks sized to the allocation
│
├── tests/                   ← runnable self-checks (see below)
└── examples/
    └── grape/               ← worked example (original .Rmd / .R sources)
```

Each run writes into its own project directory, not here:

```
~/rnaseq_projects/<name>/
├── config/                  ← config.yaml, samples.tsv, thresholds.yaml
├── results/                 ← counts, DE_Results, Enrichment, WGCNA
├── gates/
│   ├── preflight_issues.tsv ← codes, severity, remedy
│   ├── preflight_plan.json  ← counts, design facts, planned stages
│   └── decisions.log        ← every automatic decision this run made
├── metrics/                 ← per-stage *.json (drives the next stage)
└── logs/                    ← per-rule logs

~/rnaseq_reference_cache/<accession>/   ← transcriptome, GTF, salmon index
                                          (shared by projects on that assembly)
```

---

## Pipeline status

| Stage | Implemented | Executed on real data |
|---|---|---|
| 0 — fetch reference + index | yes (`workflow/rules/retrieve.smk`) | yes |
| 1 — per-sample QC + quant | yes (`workflow/scripts/01_qc_quant.sh`) | yes |
| 2 — aggregate counts | yes (`workflow/scripts/02_aggregate.R`) | **not yet** |
| 3 — differential expression | yes (`workflow/scripts/03_de.R`) | **not yet** |
| 4 — enrichment (GO + KEGG) | yes (`workflow/scripts/04_enrichment.R`) | **not yet** |
| 5 — WGCNA | yes (`workflow/scripts/05_wgcna.R`) | **not yet** |

All five stages are written, as parameterized refactors of the `examples/grape/`
files. Stages 2 to 5 have **not** been run end to end yet: FARM's `R/4.4.2`
module does not provide tximport, variancePartition, clusterProfiler or WGCNA, so
they need the conda environment in `environment.yml`. Creating that environment is
the next step before trusting any downstream output.

## Preflight

`submit.sh` runs `scripts/preflight.py` first and refuses to submit if it reports
an error. It checks, in seconds and before any data is fetched:

- configuration keys against `config/config.template.yaml`, so a typo is rejected
  rather than silently ignored
- the assay, design family and DE backend against `config/spec.yaml`, which is
  the machine-readable record of what this pipeline claims to support
- the sample sheet: duplicate ids, blank ids, model columns that do not exist,
  model variables with only one level
- **model estimability**: rank deficiency, residual degrees of freedom, and
  independent biological replicates. This previously only surfaced in Stage 3,
  after the whole cohort had been quantified.

It writes `gates/preflight_issues.tsv` (stable issue codes, severity, remedy) and
`gates/preflight_plan.json` (counts, design facts, and which stages will run and
why not). The estimability check calls the same `workflow/scripts/_design.R` that
Stage 3 uses, so preflight and the DE stage cannot reach different verdicts.

## Self-checks

No test framework; each file is a script that exits non-zero on failure.

```bash
Rscript tests/check_de.R                # contrasts, design estimability, sample disposition
Rscript tests/check_wgcna.R             # block sizing, recorded parameters
python3 tests/check_setup_helpers.py    # URL resolution, mate detection, sample matching
python3 tests/check_qc_report.py        # missing-vs-zero metrics, stale sample dirs
python3 tests/check_preflight.py        # issue codes, scope matrix, estimability gate
bash    tests/check_storage_policy.sh   # intermediate cleanup, source FASTQs preserved
bash    tests/check_project_isolation.sh  # two projects cannot touch each other's state
bash    tests/run_all.sh                # all of the above, one line each
```

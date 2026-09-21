# RNA-Seq Pipeline

A generic, HPC-aware, adaptive RNA-Seq analysis pipeline. Runs on any dataset with a sample sheet and a config file — no code changes needed to swap species, experimental design, or project.

Built for a **20 GB storage budget** on a SLURM cluster, with **per-sample fragmentation** (one FASTQ in flight at a time) and **adaptive downstream decisions** (DE backend, batch correction, enrichment filtering, WGCNA skip) driven by metrics from prior stages.

---

## Quick start

```bash
# 1. Fill in 7 minimal fields — everything else is auto-generated
cp scripts/setup_inputs.template.yaml setup_inputs.yaml
$EDITOR setup_inputs.yaml            # project name, taxID, fastq source, metadata file, ...

# 2. Run the setup wizard — produces config.yaml, samples.tsv, annotation,
#    Ensembl map, and OrgDb (either installs from Bioconductor or builds locally)
python scripts/setup.py setup_inputs.yaml

# 3. Launch — the submitter counts samples and picks the right SLURM profile
./submit.sh

# 4. Resume after any interruption (same command; Snakemake picks up where it stopped)
./submit.sh
```

The first run also fetches the reference genome and builds a Salmon index (one-time, ~1 GB).
See [**INPUTS.md**](INPUTS.md) for the full list of what the wizard generates.

---

## Where to read next

| Doc | Read when you want to… |
|---|---|
| [**INPUTS.md**](INPUTS.md) | Use the setup wizard — 7 inputs, auto-generates everything external |
| [**DESIGN.md**](DESIGN.md) | Understand the architecture, stages, decision gates, and storage budget |
| [**STRATEGY.md**](STRATEGY.md) | Understand how `./submit.sh` auto-picks small / medium / large SLURM tier based on sample count |
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
├── config/
│   ├── config.template.yaml
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
│       ├── 01_qc_quant.sh           ← fastp | salmon | delete FASTQ
│       ├── 02_aggregate.R           ← (pending) tximport
│       ├── 03_de.R                  ← (pending) adaptive DE
│       ├── 04_enrichment.R          ← (pending) GO + KEGG
│       └── 05_wgcna.R               ← (pending) adaptive WGCNA
│
├── examples/
│   └── grape/               ← worked example (original .Rmd / .R sources)
│
├── gates/                   ← runtime: decisions.log (every auto-decision)
└── metrics/                 ← runtime: per-stage *.json (drives next stage)
```

---

## Pipeline status

| Stage | Status |
|---|---|
| 0 — fetch reference + index | ✅ Implemented (`workflow/rules/retrieve.smk`) |
| 1 — per-sample QC + quant | ✅ Implemented (`workflow/scripts/01_qc_quant.sh`) |
| 2 — aggregate counts | ✅ Implemented (`workflow/scripts/02_aggregate.R`) |
| 3 — differential expression | ✅ Implemented (`workflow/scripts/03_de.R`) |
| 4 — enrichment (GO + KEGG) | ✅ Implemented (`workflow/scripts/04_enrichment.R`) |
| 5 — WGCNA | ✅ Implemented (`workflow/scripts/05_wgcna.R`) |

Pipeline is end-to-end implemented. Parameterized refactors of the `examples/grape/` files.

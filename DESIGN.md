# Generic RNA-Seq Pipeline — HPC-Aware, Adaptive, Config-Driven

A reusable RNA-Seq analysis pipeline that runs on any dataset, fits in a **20 GB HPC home quota**, processes samples in **per-sample fragments**, and **adapts downstream parameters** based on metrics from earlier stages.

The wine grape / GRBV study is the worked example; nothing in the pipeline is hardcoded to *Vitis vinifera*.

---

## 1. Design principles

1. **One `config.yaml` per project.** Species, paths, model formula, contrasts, thresholds — all in config. Scripts never hardcode biology.
2. **Fragmented execution.** Snakemake DAG with per-sample job arrays. Any interruption is resumable.
3. **Aggressive cleanup.** Raw FASTQ is deleted immediately after quantification; BAMs never exist (Salmon is pseudoalignment).
4. **Metrics → decisions.** Each stage emits a small `*_metrics.json`. The next stage reads it and picks parameters (or tools) accordingly. Every auto-decision is logged to `gates/decisions.log` for audit.
5. **Fail loud, fail early.** Decision gates warn when a dataset exits the sanity range (low mapping, confounded batch, etc.) rather than silently proceeding.

---

## 2. Storage budget (20 GB hard cap)

| Artifact | Size | Persistent? |
|---|---|---|
| Salmon index (plant/mammal) | ~1 GB | yes |
| GTF + annotation tables | ~300 MB | yes |
| One raw FASTQ in flight | 0.2–3 GB | transient |
| `quant.sf` per sample | ~10 MB | yes |
| DE / GO / KEGG / WGCNA outputs | < 500 MB total | yes |
| **Steady-state persistent** | **~3–4 GB** | |
| **Peak (processing 1 sample)** | **~5–7 GB** | |

Leaves ~12 GB headroom. **STAR is rejected** — its index alone is 15–25 GB. HISAT2 (~4 GB index) is a fallback option, BAMs still must not persist.

---

## 3. Architecture

```
┌─────────────────────────────────────────────────────────────┐
│  config/config.yaml          (species, paths, model, etc.)  │
│  config/thresholds.yaml      (all decision-gate cutoffs)    │
│  config/samples.tsv          (sample sheet)                 │
└─────────────────────────────────────────────────────────────┘
                            │
                            ▼
                     ┌──────────────┐
                     │  Snakefile   │   orchestrates DAG, SLURM submission
                     └──────────────┘
                            │
    ┌───────────────────────┼────────────────────────┐
    ▼                       ▼                        ▼
Stage 0              Stage 1–2                Stage 3–5
retrieve           per-sample jobs           aggregate R
(ref, anno,        (job array on             (runs once at end)
 samplesheet)       $SLURM_ARRAY_TASK_ID)

    │                       │                        │
    ▼                       ▼                        ▼
 metrics/             quant/{sample}/           DE_Results/
 retrieve.json        per-sample-metrics.json   GO_results/
                                                KEGG_results/
                                                WGCNA/
                                                gates/decisions.log
```

---

## 4. Stages (generic form)

### Stage 0 — Retrieve
- **Inputs:** URLs in `config.yaml` (genome, GTF, sample FASTQ manifest)
- **Process:** Build Salmon index once; pin versions; hash-verify downloads
- **Outputs:** `reference/salmon_idx/`, `reference/annotation.tsv`
- **Decision gate:** If genome > 3 GB or transcriptome > 200k seqs → warn, check `k` parameter
- **Cleanup:** Keep compressed FASTA; drop uncompressed

### Stage 1 — Per-sample QC + trim + quant  (**fragmented; one SLURM task per sample**)
- **Inputs:** One raw FASTQ URL + Salmon index
- **Process:** `wget` or `aws s3 cp` → `fastp` → stdin → `salmon quant` (no intermediate FASTQ written)
- **Outputs:** `quant/{sample}/quant.sf`, `quant/{sample}/metrics.json`
- **Decision gate:**
  - Mapping rate < 60% → **flag sample**, recommend alternate reference
  - Reads assigned to genes < project threshold → **flag sample**, downstream excludes it
- **Cleanup:** `rm` the FASTQ as soon as `quant.sf` validates

### Stage 2 — Aggregate counts
- **Inputs:** all `quant/*/quant.sf` + sample sheet
- **Process:** `tximport` to gene-level counts, apply annotation, export matrix
- **Outputs:** `counts.tsv`, `counts_metrics.json` (sample drop-outs, depth distribution)
- **Decision gate:**
  - Library size CV > 3× or any sample < median/10 → warn
  - # samples < 3/group → downstream picks **edgeR exactTest** instead of voom/dream

### Stage 3 — Differential expression  (*adapts to dataset shape*)
- **Inputs:** `counts.tsv`, `config.yaml`, `counts_metrics.json`
- **Process:**
  1. Normalize (TMM); `filterByExpr`
  2. MDS/PCA → compute `var_batch` vs `var_condition`; if batch explains > 30% of PC1–2 and isn't confounded with condition → **add batch covariate** automatically (log decision)
  3. Pick DE backend:
     - n/group ≥ 3 **and** has random effect in formula → `dream`
     - n/group ≥ 3 **no** random effect → `limma-voom`
     - n/group < 3 → `edgeR exactTest`
  4. Fit, extract contrasts via `emmeans` from `config.contrasts`
- **Outputs:** `DE_Results/{contrast}.tsv`, `de_metrics.json` (# sig genes per contrast)
- **Decision gate:**
  - Any contrast with 0 sig genes → downgrade to non-directional GSEA (`|t|`) for that contrast
  - Any contrast with > 5000 sig genes → flag likely confounder, recommend re-check batch

### Stage 4 — Functional enrichment  (*adapts to DE result shape*)
- **Inputs:** `DE_Results/*.tsv`, `config.organism` (OrgDb name + KEGG code)
- **Process:**
  - **GO GSEA:** rank by `|t|` (non-directional) — robust default from clusterProfiler best-practice
  - **KEGG GSEA:** rank by signed `t` (directional)
  - If `# significant terms > 20` → auto-apply `clusterProfiler::simplify(cutoff=0.7, by="p.adjust", select_fun=min)`
- **Outputs:** `GO_results/{contrast}.tsv`, `KEGG_results/{contrast}.tsv`, dotplots
- **Decision gate:**
  - Custom OrgDb not installed → emit helper: *"run `make_orgdb.R` with tax_id=X"*
  - Organism not in KEGG list → skip KEGG step, log skip

### Stage 5 — WGCNA  (*adaptive power selection*)
- **Inputs:** log2-CPM from Stage 3 (filtered matrix only)
- **Process:**
  - `pickSoftThreshold` over 1:60; pick lowest power with scale-free R² ≥ 0.85
  - Cap at 30 (function limit); if R² at 30 < 0.8 → **warn in report**
  - `blockwiseModules`; if > 50 modules → raise `mergeCutHeight` and rerun once
- **Outputs:** `WGCNA/MEs.csv`, `Module_genes.xlsx`, dendrogram, hub-gene heatmaps, per-module GO
- **Decision gate:** Skip WGCNA entirely if n_samples < 15 (unreliable below this threshold)

### Stage 6 — Report
- Knit each stage `.Rmd` to HTML
- Aggregate `gates/decisions.log` into an "Auto-decisions" section of the final report

---

## 5. Decision-gate matrix (summary)

| Gate | Metric | Threshold | Auto-action |
|---|---|---|---|
| Reference choice | Mapping rate | < 60% | Flag, suggest alt reference |
| Sample QC | Reads on genes | < 1M (TAGseq) / < 10M (RNA-Seq) | Exclude sample |
| Batch correction | PC1–2 var by batch | > 30% and not confounded | Add batch covariate |
| DE backend | n/group, random effect | see Stage 3 logic | Pick edgeR / voom / dream |
| Enrichment cleanup | # sig terms | > 20 | Apply `simplify()` |
| KEGG skip | Organism code | not in KEGG | Skip step |
| WGCNA skip | n samples | < 15 | Skip step |
| WGCNA power | R² at power 30 | < 0.8 | Warn, use power=30 |
| Module merge | # modules | > 50 | Re-run with higher `mergeCutHeight` |

All thresholds live in `config/thresholds.yaml` — editable per project.

---

## 6. Fragmentation — how it fits in 20 GB

**Snakemake job array, one task per sample:**

```bash
# submit.sh — launch as a SLURM job array sized to N samples
sbatch --array=1-$(wc -l < config/samples.tsv) run_sample.slurm
```

`run_sample.slurm` runs one sample end-to-end (download → trim → quant → delete FASTQ → emit metrics), then exits. Each task's peak disk is one FASTQ + a Salmon output dir; Snakemake only launches the next sample-task after the prior FASTQ is cleaned up.

Once every `quant.sf` exists, Snakemake fires the **single aggregate job** (Stages 2–6) on the tiny count matrix — that job is memory-bound, not storage-bound.

**Resume behavior:** Snakemake checkpoints on artifact presence. Killing the job mid-array and restarting re-queues only the un-finished samples.

---

## 7. Swapping the worked example

To run this on *any* new dataset, only `config/config.yaml`, `config/samples.tsv`, and — if the species is new to Bioconductor — a one-time `make_orgdb.R` invocation need to change. **No `.Rmd` or `.R` script needs editing.**

See [`config/config.template.yaml`](config/config.template.yaml) for the full schema.

---

## 8. Files this design adds/changes

```
config/
  config.template.yaml       ← new, copy per project
  thresholds.yaml            ← new, default decision-gate cutoffs
  samples.tsv                ← project-specific sample sheet
Snakefile                    ← new, orchestrator
workflow/
  rules/
    retrieve.smk             ← Stage 0
    per_sample.smk           ← Stage 1 (SLURM job array)
    aggregate.smk            ← Stages 2–6
  scripts/
    01_qc_quant.sh           ← wraps fastp | salmon
    02_aggregate.R           ← tximport + counts_metrics.json
    03_de.R                  ← parameterized DE (was analysis_20*.Rmd)
    04_enrichment.R          ← parameterized GO/KEGG (was GO_20*.Rmd, KEGG_20*.Rmd)
    05_wgcna.R               ← parameterized WGCNA (was WGCNA.Rmd)
gates/
  decisions.log              ← every auto-decision, with threshold + value
metrics/
  *.json                     ← per-stage metric dumps
```

Existing `analysis_2021.Rmd` etc. stay as the **project-specific** example — the generic `03_de.R` is a refactor of them with all grape-specific code moved to config.

---

## 9. Quick run (on HPC)

```bash
# one-time per project
cp config/config.template.yaml config/config.yaml
$EDITOR config/config.yaml                      # fill in species, paths, model
$EDITOR config/samples.tsv                       # paste sample manifest

# if non-Bioconductor organism
Rscript workflow/scripts/make_orgdb.R           # reads config.yaml

# dry-run the DAG
snakemake -n --configfile config/config.yaml

# submit to SLURM (fragmented)
snakemake --profile slurm --configfile config/config.yaml

# resume after interruption
snakemake --profile slurm --configfile config/config.yaml   # picks up where it left off
```

# Handbook insights — pipeline assessment

> Historical interpretation, retained for context. Current implementation and qualification are documented in [README.md](README.md) and [docs/RELEASE.md](docs/RELEASE.md). The earlier exclusions of RNA assay families and fixed-storage rationale below are superseded; they are not release guarantees.

Source: `~/Desktop/sequencing handbook.pdf` (Sequencing Techniques — A Working Reference, April 2026).

This doc captures what we took from the handbook and made part of the pipeline, what we examined and rejected, and what's still on the table for later. It's the "why" behind the choices that aren't obvious from the code alone.

---

## What the handbook validated about our current design

Reading the handbook end-to-end was reassuring — the pipeline is already aligned with most of its recommendations:

- **Salmon pseudoalignment** when DE is the primary goal (handbook p47) — our choice
- **`tximport` → DESeq2 / limma-voom / DREAM** chain (p47) — our flow
- **DREAM for repeated-measures / longitudinal / paired** (p47) — exactly what our auto-backend rule picks when a random effect is present in the model
- **Minimum 3 biological reps / group; 5+ much better; 12+ for subtle effects** (p60) — our `de_backend_rules` thresholds match this
- **Benjamini-Hochberg FDR at q < 0.05** (p60) — what `limma`/`dream` give us
- **Snakemake + conda + Singularity / Apptainer** as the reproducibility stack (p53) — our orchestration choice
- **Raw FASTQ on cold storage; BAMs deleted; counts/DMRs on fast storage + backed up** (p53) — exactly our 20 GB strategy

Reassurance, not work.

---

## What we implemented (the "Tier 1" bundle)

Five concrete changes shipped on 2026-05-29:

### 1. Salmon technical-correction flags (`--seqBias`, `--posBias`)

**Source:** handbook section 7.2.2 ("Salmon with `--gcBias`, `--seqBias`, `--posBias` for technical correction").

**Where:** `workflow/scripts/01_qc_quant.sh` — the Salmon invocation.

**What changed:** previously we ran only `--gcBias`. Now:
- Standard RNA-seq (`rnaseq_single`, `rnaseq_paired`) → `--gcBias --seqBias --posBias`
- TAGseq → `--gcBias` only

**Why the TAGseq exception:** 3′-tag libraries are intentionally 3′-end-biased — applying positional-bias correction would fight the assay's biology and corrupt counts. This is implied by the handbook's separate TAGseq treatment (section 3.2) and confirmed by Salmon documentation.

**Recorded in metrics:** `metrics.json` now includes `salmon_bias_flags` so the chosen flags are visible per sample.

---

### 2. Strandedness verification gate

**Source:** handbook section 4.1.3 and the failure-modes table (p61): *"RNA-seq counts ~50% expected → Strandedness misconfigured → Verify with `infer_experiment.py` or Salmon library type A."*

**Where:** `workflow/scripts/01_qc_quant.sh` + new config field `samples.expected_libtype`.

**What changed:**
- Salmon's auto-detected `library_types[0]` is now extracted from `aux_info/meta_info.json` and recorded as `detected_libtype` in `metrics.json`.
- If `config.samples.expected_libtype` is set (e.g. `"ISR"` for paired dUTP-stranded), the pipeline case-insensitively compares it to what Salmon detected. Mismatch → `WARN` line in `gates/decisions.log` and `libtype_mismatch: true` in `metrics.json`.
- Default behavior is unchanged when `expected_libtype` is `null` — Salmon auto-detect runs as before, just with telemetry.

**Why it matters:** the handbook flags this as one of the most common silent failure modes — a wrong stranded setting halves counts and quietly invalidates DE analysis. We now catch it.

**Common library-type codes** (paste into `expected_libtype`):
- `ISR` — paired-end, dUTP-stranded (NEBNext Ultra II Directional, TruSeq Stranded mRNA)
- `ISF` — paired-end, forward-stranded
- `IU` — paired-end, unstranded
- `SR` — single-end, stranded-reverse (Lexogen QuantSeq, most TAGseq)
- `SF` — single-end, forward-stranded
- `U` — single-end, unstranded

---

### 3. Effect-size filter alongside FDR

**Source:** handbook section 10.2.3 (p60): *"A statistically significant fold change of 1.05 with q = 0.001 across 30,000 cells is real but biologically uninteresting. Always filter on both."*

**Where:** `workflow/scripts/03_de.R` — both the voom/dream branch and the edgeR exactTest branch.

**What changed:** for each contrast, the pipeline now writes **two** files:
- `<contrast>_DE_analysis.tsv` — full results (unchanged)
- `<contrast>_DE_meaningful.tsv` — filtered on **both** `adj.P.Val < thresholds.de_meaningful.adjp_max` (default 0.05) **and** `|logFC| >= thresholds.de_meaningful.log2fc_min` (default 1.0 = 2-fold change)

**Why two files:** downstream users (collaborators, enrichment, publication tables) almost always want the meaningful subset. Keeping the full table available means you don't lose flexibility for re-thresholding later.

---

### 5. MultiQC report aggregation  →  superseded by `qc_report` rule

**Source:** handbook section 7.1 ("Universal early steps") and the curated tools index (p62).

**Originally planned:** new rule `workflow/rules/multiqc.smk` running MultiQC standalone.

**Replaced by:** Ritu contributed a richer implementation in `workflow/rules/qc_report.smk` + `workflow/scripts/qc_report.py`. Their rule does MultiQC **plus** generates:
- `comparative_charts.png` — before/after-clean read counts + alignment-fate stacked bars per sample
- `qc_charts.html` — self-contained HTML embedding the PNG + per-sample alignment table
- `alignment_summary.tsv` — flat table for ad-hoc analysis
- `multiqc/multiqc_report.html` — the full MultiQC report

Strictly more useful than MultiQC alone — kept their version, dropped the planned `multiqc.smk`.

**Dependency:** `multiqc` is still in `environment.yml` (required by `qc_report.smk`). Plus `matplotlib` for the PNG generator (add to `environment.yml` if missing).

---

### 7. Version-pinned conda environment

**Source:** handbook section 7.15 (p53): *"One environment per pipeline; version-pin every tool (`fastp=0.23.4`, not `fastp`)."*

**Where:** new `environment.yml` at the pipeline root.

**What changed:**
- Pinned exact versions for `python`, `snakemake`, `pandas`, `pyyaml`, `openpyxl`.
- `snakemake-executor-plugin-slurm` and `multiqc` left as floating minor versions (compatibility is robust; pinning these too aggressively breaks rebuilds).
- File includes the `conda env create -f environment.yml` recipe in the header.
- Header also documents the `conda env export --no-builds > environment.lock.yml` workflow for fully locked re-creation after first build.

**Why it matters:** an unpinned `conda create ... snakemake pandas` resolves to whatever's latest on the day you run it — which can break months later. The pinned file is the difference between "I ran this analysis in May 2026" and "I can re-run this analysis identically in May 2027."

---

## Real-world hardening pulled from Ritu (parallel work, merged in)

While I was implementing the handbook bundle locally, Ritu was actively running a 12-sample paired-end fungal RNA-seq dataset on the validation cluster. Her version of the pipeline made several improvements I hadn't anticipated. Merged into local on 2026-05-30:

- **Paired-end support** — `samples.tsv` gained a `fastq_url_r2` column; `common.smk` gained `get_fastq_url_r2()`; `per_sample.smk` passes `--url2`; `01_qc_quant.sh` runs fastp with `--detect_adapter_for_pe` and feeds R1+R2 to salmon. The pipeline now handles SE and PE in one rule.
- **No-streaming for `--gcBias`** — discovered the hard way that salmon's `--gcBias` re-reads the input, which deadlocks if it comes from a pipe. PE branch writes trimmed reads to disk before salmon reads them. SE branch still streams (single-pass `--gcBias` is OK there).
- **Realistic resources** — `mem_mb=16000`, `runtime=90`. Old defaults (6000/45) OOMed on real data with a decoy-aware Salmon index.
- **`--trim_poly_g`** — fastp removes the poly-G tails NextSeq/NovaSeq two-color chemistry introduces. Standard for any modern Illumina dataset.
- **Strict-mode-safe module load** — wraps `source /etc/profile.d/modules.sh` in `set +euo pipefail` to avoid the "conda deactivate: command not found" → exit 127 trap that kills the whole script under `set -e`.
- **Custom QC report system** — `qc_report.smk` + `qc_report.py` (item #5 above).
- **Gene matrix builder** — `workflow/scripts/build_gene_matrix.py`. Collapses transcript-level `quant.sf` files to gene-level counts + TPM matrices, attaches RefSeq descriptions from the transcriptome FASTA. Independent of `02_aggregate.R` (uses GTF directly instead of tximport). Useful as a quick downstream-ready output for collaborators who want a spreadsheet.

These earned changes are now part of the local source-of-truth. They didn't come from the handbook — they came from the discipline of actually running real data.

## Tier 1 changes also updated

Two changes adjacent to the main bundle, also from the handbook:

### Threshold defaults bumped to handbook ranges

In `config/thresholds.yaml`:
- `min_reads_on_genes_tagseq`: 1,000,000 → **3,000,000** (handbook p7: 3–10M typical for TAGseq)
- `min_reads_on_genes_rnaseq`: 10,000,000 → **20,000,000** (handbook p7: 20–40M for bulk RNA-seq DE)
- `duplication_rate_warn`: 0.80 → **0.30** (handbook Appendix A: > 30% = over-amplified)

The old defaults were "won't break things" floors; the new defaults are "what a healthy library looks like."

### New `de_meaningful` thresholds section
- `adjp_max: 0.05`
- `log2fc_min: 1.0`

Drives item #3.

---

## What we examined but deferred

These showed up in the handbook and are worth doing eventually — but not in this bundle.

### Dual host-pathogen quantification (handbook 8.2)

**The opportunity:** concatenate the host transcriptome with a pathogen transcriptome (e.g. GRBV) into one Salmon index, quantify both organisms in a single pass. Gives viral transcript abundance alongside host gene counts from the same library.

**Why deferred:** half-day of work and only relevant for plant pathology projects with a known viral target. Not every project needs it. Easy to add when a real use-case shows up — would extend `retrieve.smk` to fetch + concatenate FASTAs before indexing, and `config.yaml` would gain a `reference.pathogen_transcriptome_url` field.

### Power-analysis pre-flight (handbook 10.2.5)

**The opportunity:** before launching SLURM jobs, `submit.sh` could compute the smallest detectable fold change for the design (3v3 → 2-fold, 5v5 → 1.5-fold, 12v12 → 1.2-fold per handbook p60) and warn if under-powered for the expected biology.

**Why deferred:** the math is simple but the UX needs care — false alarms on small pilot studies would be annoying. Better to add once we know how often it'd actually fire.

### Adapter-dimer / over-amplification flag (handbook 10.3)

**The opportunity:** fastp already detects adapter content; we could add a gate that flags samples with `> 20%` adapter content or `> 30%` PCR duplicates as "library prep issue, re-check."

**Why deferred:** MultiQC (item #5) already surfaces this visually in the HTML report. Adding an automatic gate on top is value, but lower priority. Easy if/when needed.

### Plant functional-annotation databases (handbook 7.2.5)

**The opportunity:** add support for **PlantCyc**, **Mapman**, **PLAZA** alongside GO/KEGG. For grape specifically: **Grape Genomic Database (URGI)**, **VitisVea**, **PlantCyc Vinifera**.

**Why deferred:** clusterProfiler handles GO and KEGG natively; PlantCyc / Mapman require custom term-to-gene mappings that aren't packaged uniformly. Worth doing if/when we start interpreting plant-specific pathway results seriously.

---

## What we examined and rejected

Some handbook content doesn't apply to this pipeline and would just add complexity:

| Handbook section | Why we skipped it |
|---|---|
| ChIP-seq / ATAC-seq / CUT&RUN | Different assays — different alignment (bowtie2/bwa), different downstream (MACS2). Belongs in a sister pipeline if needed. |
| Hi-C / 3C / Micro-C | Same — entirely separate pipeline (HiC-Pro, Juicer). |
| Ribo-seq, CLIP-seq | Specialized assays with their own tool chains. |
| Single-cell, spatial transcriptomics | Different data structure (cell × gene matrix vs. sample × gene). Cell Ranger / Space Ranger / Seurat universe. |
| Long-read RNA-seq (Iso-Seq, ONT) | Different aligners (minimap2), different isoform tools (SQANTI3, IsoQuant). Add as a sister pipeline if needed. |
| Variant calling, methylation, metagenomics | Out of scope for an RNA-seq pipeline. |

We can absolutely build separate pipelines for any of these in the future. They just don't belong in this one.

---

## Reading the handbook before adding new techniques

The end-notes (p71) give the right next move when a new technique is needed:

> The right next step when a new technique is needed: find a recent application paper using that technique on tissue similar to yours, replicate their pipeline first, then optimize.

Worth keeping in mind whenever we extend the pipeline beyond what's here.

---

*Last updated: 2026-05-29.*

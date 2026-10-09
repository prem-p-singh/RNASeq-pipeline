# Changelog

## Unreleased

- Removed the n8n orchestration integration (`integrations/n8n/`) and its three checks. It was not part of the planned 2.1.0 scope. The release gate is 35 checks. Qualification of this source: pending.

## 2.1.0 — 2026-10-08

Scope: bulk completion (milestones M1–M3). Later assay families are planned for 3.0.0 and beyond.

- Workbook schema 5: analysis goal (differential expression, expression only, QC only), repeated measures with a subject column (random subject with dream, or fixed subject block for pairing), optional fixed-effects formula, optional Contrasts tab (pairwise with levels or linear weights), named library protocol, and contamination screening questions.
- A changed workbook reports the changed settings and affected stages; `submit.sh --accept-revision` applies it and keeps the previous configuration.
- Named non-UMI bulk kit profiles (TruSeq, NEBNext) set adapters and strand; conflicting answers are rejected. The sequencing platform sets poly-G trimming unless chosen explicitly.
- Screening panels accept FASTA sources (path, URL, project transcriptome or genome) that the workflow indexes automatically, with the build version and FASTA checksum recorded.
- New projects record the handbook DE default (DESeq2 for 3–12 biological units per group, limma-voom above 12); `auto` in existing configurations is unchanged. Two units per group needs an explicit method.
- Design checks compare an interaction (DESeq2), a donor-blocked pairing (edgeR QL) and a continuous covariate (limma-voom) against direct fits.

- Optional FastQ Screen/Bowtie2 diagnostic panels: seeded fragment sampling, separate mate/read counts, reference checksums and explicit disabled status; never removes reads or excludes samples. Isolated locked deployment preserves the core interpreter.

- Workbook schema 4 connects bulk quantifier/DE selection, complete custom reference sets, decoy profile and basic preprocessing to generated configuration; older schemas remain readable.

- Independent, observation-only raw QC for canonical non-UMI bulk projects using locked fastp 0.23.4; no reference or statistical design required.
- Durable per-sample reports and run status, input checksum/mate validation, exact command provenance and protected source reads. Existing output directories are refused to prevent stale results.
- Separate resumable preprocessing and quantification producers, explicit adapter/quality/poly-G settings, read-loss provenance and an early MultiQC report.
- QC-only, expression-only and coexpression objectives; local/SLURM launcher selection and direct Excel launch/resume with snapshot checks.
- Fixed-effect edgeR QL and DESeq2 Wald adapters, shared preflight contrast validation, explicit pair direction and named coefficient contrasts.
- Full-genome Salmon decoy construction, identifier checks and reference/cache integrity verification.
- STAR/featureCounts bulk route with declared strand, uniquely aligned fragment counts, indexed BAMs and a separately pinned tool module; public (38992390) and recovery (38992412) checks passed.
- Candidate Linux runtime adds DESeq2 1.50.2 while preserving the previous 579 package records.
- This branch combines the deployment-qualified source `82e70c5` with the 2.1.0 work. The combined source `0a857c6` passed the 38-check gate with zero skips and both public routes (job 39620298), SLURM-worker execution, resume and the shared reference cache (39620305) and the stock container check (39618728). See docs/RELEASE.md. The hosted release-validation workflow has not run on this source.
- Deployment qualification (M8) passed for this branch's bulk scope on source `82e70c5`: fresh installation, the 37-check gate with zero skips, both public routes, launcher-to-SLURM-worker execution, resume after interruption, two projects on one reference cache entry and a stock Ubuntu container. See docs/RELEASE.md for jobs, limits and the release matrix. The hosted release-validation workflow has not run on this source.
- SLURM profiles wait 120 seconds for worker outputs (30 or 60 before). On NFS mounts that cache directories for 60 seconds, a resumed or retried job was reported as missing its outputs.
- `scripts/qualify_deployment.sbatch` runs the distributed and container checks; `tests/check_public.sh --verify` checks a project that another launcher ran; `RNASEQ_VALIDATE_ROOT` relocates the validation runtime.
- M1–M7 remain incomplete; new backend qualification is in progress. M2 remains partial: named chemistry profiles and new assay support remain open; declared-panel screening does not establish panel sensitivity.
- Full regression: 35 checks, 0 failed, 0 skipped, plus Salmon (with screening) and STAR public routes (39010068). This run predates the n8n checks and the deployment qualification below. M4–M7 remain incomplete. M2 remains partial: named chemistry profiles and new assay support remain open; declared-panel screening does not establish panel sensitivity.

## 2.0.0 — 2026-09-23

V2 releases the currently implemented **annotated, non-UMI bulk RNA-seq workflow**. The broader multi-assay roadmap is not part of this release's supported scope. Workbook schema version 3 is independent of software version 2.0.0.

### Added

- Direct Excel setup, stable field IDs, source workbook/checksum and cell-level provenance.
- Samples, Libraries and Reads tabs with explicit ownership of inline versus external metadata.
- Executable canonical records: multiple sequencing runs/lanes within one library per sample, mate-ID/count and FASTQ validation, optional SHA256/size verification and retained read provenance.
- Recovery when canonical read provenance is missing; text IDs including leading zeros and literal `NA` survive setup and analysis.
- Generic SLURM validation and explicit SSH-host configuration for the optional uploader.
- README workflow diagram, complete installation/intake/run/resume guidance and capability boundaries.

### Changed and repaired

- Storage estimates warn without blocking launch or reducing concurrency. No default 20 GB limit; estimates include complete libraries and retained merged/trimmed reads.
- PCA batch diagnostics preserve the declared scientific model and record review suggestions.
- Complete canonical-table relationships are validated; duplicate roles are retained for diagnostics instead of overwritten.
- Empty shell arguments are quoted correctly in both legacy and canonical execution.
- Environment verification uses the existing 579-package Linux lock. No new scientific runtime dependencies were added in this upgrade.

### Qualification and limits

The strict suite has 22 checks, including real synthetic FASTQ-to-report processing, planted DE effects, quantitative count agreement, WGCNA, enrichment, recovery and project isolation. Release publication requires passing scientific GitHub CI, including the separate public-data smoke test. See [release evidence](docs/RELEASE.md).

TAG-seq has no qualified kit; dream lacks complete end-to-end mixed-model qualification; live KEGG/annotation services have separate limitations. UMI, small-RNA, cell/nucleus, spatial, long-read, dual-organism and specialized routes remain unimplemented. See [the capability table](README.md#what-you-can-run).

# Changelog

## Unreleased

- Add the n8n orchestration handoff and an explicit three-test integration gate. It targets the separately reviewed schema-7 development source; this schema-4 parent gains no workbook revision, TAG-seq/QuantSeq, or small-RNA routes. Missing release-smoke prerequisites are visible skips and block strict qualification. See integrations/n8n/README.md and VALIDATION.md.

- Optional FastQ Screen/Bowtie2 diagnostic panels: seeded fragment sampling, separate mate/read counts, reference checksums and explicit disabled status; never removes reads or excludes samples. Isolated locked deployment preserves the core interpreter.

- Workbook schema 4 connects bulk quantifier/DE selection, complete custom reference sets, decoy profile and basic preprocessing to generated configuration; older schemas remain readable.

- Independent, observation-only raw QC for canonical non-UMI bulk projects using locked fastp 0.23.4; no reference or statistical design required.
- Durable per-sample reports and run status, input checksum/mate validation, exact command provenance and protected source reads. Existing output directories are refused to prevent stale results.
- Separate resumable preprocessing and quantification producers, explicit adapter/quality/poly-G settings, read-loss provenance and an early MultiQC report.
- QC-only, expression-only and coexpression objectives; local/SLURM launcher selection and direct Excel launch/resume with snapshot checks.
- Fixed-effect edgeR QL and DESeq2 Wald adapters, shared preflight contrast validation, explicit pair direction and named coefficient contrasts.
- Full-genome Salmon decoy construction, identifier checks and reference/cache integrity verification.
- STAR/featureCounts bulk route with declared strand, uniquely aligned fragment counts, indexed BAMs and a separately pinned tool module; public and recovery qualification in progress.
- Candidate Linux runtime adds DESeq2 1.50.2 while preserving the previous 579 package records.
- Deployment qualification (M8) passed for this branch's bulk scope on source `82e70c5`: fresh installation, the 37-check gate with zero skips, both public routes, launcher-to-SLURM-worker execution, resume after interruption, two projects on one reference cache entry and a stock Ubuntu container. See docs/RELEASE.md for jobs, limits and the release matrix. The hosted release-validation workflow has not run on this source.
- SLURM profiles wait 120 seconds for worker outputs (30 or 60 before). On NFS mounts that cache directories for 60 seconds, a resumed or retried job was reported as missing its outputs.
- `scripts/qualify_deployment.sbatch` runs the distributed and container checks; `tests/check_public.sh --verify` checks a project that another launcher ran; `RNASEQ_VALIDATE_ROOT` relocates the validation runtime.
- M1–M7 remain incomplete; new backend qualification is in progress. M2 remains partial: named chemistry profiles and new assay support remain open; declared-panel screening does not establish panel sensitivity.

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

# Changelog

## Unreleased

- Independent, observation-only raw QC for canonical non-UMI bulk projects using locked fastp 0.23.4; no reference or statistical design required.
- Durable per-sample reports and run status, input checksum/mate validation, exact command provenance and protected source reads. Existing output directories are refused to prevent stale results.
- M2 remains partial: this does not implement contamination screening, chemistry profiles or new assay support.

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

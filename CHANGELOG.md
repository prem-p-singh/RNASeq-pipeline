# Changelog

## Unreleased

Nothing yet.

## 3.0.0 (2026-11-07)

Adds three library types beside bulk RNA-seq.

### Added

- Small RNA for animals (Bowtie 1 to the species' mature miRNAs, unique hits counted) and plants (ShortStack 4.1.0 loci on the genome, with DicerCall size classes). Start it from the workbook's sRNA-seq tab or the project configuration. Reads are counted, not molecules: no UMI deduplication, isomiR resolution or phased-siRNA calling. Named kits are in `config/small_rna_kits.yaml`; TruSeq Small RNA is validated and UMI small-RNA kits are refused.
- 3′ end-tag libraries for Lexogen QuantSeq FWD, REV and FWD-UMI (`config/tagseq_kits.yaml`), from the workbook's TAGseq tab or `analysis.kit`. The kit sets trimming and the STAR/featureCounts strand. FWD-UMI checks the TATA spacer, moves the UMI into the read name and reports umi_tools molecules per gene beside read counts. All three kits are validated.
- Host and pathogen in the same bulk libraries (dual organism), from `dual_organism` in the project configuration or the Bulk RNA-seq tab (workbook schema 8). STAR aligns each read to a combined genome and labels it host, pathogen, ambiguous or unassigned; `results/assignment.tsv` reports all four per sample. Salmon or STAR then runs once per organism into `results/host/` and `results/pathogen/`. Samples with fewer than `min_assigned_reads` pathogen reads (default 500) are listed in `results/pathogen/samples_not_quantified.tsv`. Both routes are validated.
- `analysis.salmon_position_bias` (default `true`) turns Salmon's `--posBias` off when `false`. Salmon 1.10.3 crashed intermittently in this correction when a transcript length class had no reads. The pathogen project of a dual-organism run always has it off.
- Two more locked runtimes: small RNA (`environments/srna-linux-64.explicit.txt`) and UMI (`environments/umi-linux-64.explicit.txt`). They are installed only when a project needs them.

### Checks

- The release gate has 41 checks: the 35 from 2.1 plus unit and workflow checks for small RNA, QuantSeq and dual organism.
- Public-data checks: `tests/check_small_rna_public.sh --animal|--plant`, `tests/check_tagseq_public.sh --fwd|--rev|--umi` and `tests/check_dual_public.sh --virus|--cellular`, each with acceptance criteria written into its fixture before any run.
- The hosted release-validation workflow builds all five runtimes.

### Changed

- Projects started before this version recompute from preprocessing on relaunch, because the preprocessing record and the quantification script changed. Results are unchanged.

## 2.1.2 (2026-10-17)

- Added a Debian package (`packaging/`, `.github/workflows/package.yml`) providing `rnaseq-pipeline` and `rnaseq-new-project` commands. No pipeline code changed.

## 2.1.1 (2026-10-08)

- Removed the n8n automation integration and its three checks. It was not part of the planned release. No pipeline code changed. All 35 checks passed with none skipped, with the public study run through both Salmon and STAR.

## 2.1.0 (2026-10-08)

Superseded by 2.1.1.

### Added

- Intake workbook: analysis goal (differential expression, expression only, QC only); paired and repeated measures with a subject column; optional fixed-effects formula; optional Contrasts tab; named library protocol; contamination screening questions; method and reference choices.
- A changed workbook reports the changed settings and affected stages; `submit.sh --accept-revision` applies it and keeps the previous configuration.
- `submit.sh` can launch and resume directly from the workbook, locally or through SLURM.
- STAR and featureCounts route with a declared strand, uniquely aligned fragment counts and indexed BAMs, in a separately locked tool module.
- edgeR QL and DESeq2 Wald beside limma-voom and dream, with shared contrast checks, explicit pair direction and named coefficient contrasts.
- New projects record a default DE method: DESeq2 for 3 to 12 biological units per group, limma-voom above 12. `auto` in existing configurations is unchanged.
- QC-only, expression-only and coexpression goals.
- Salmon index with full-genome decoys, identifier checks and cache integrity verification.
- Named bulk library protocols (TruSeq, NEBNext) set adapters and strand; the sequencing platform sets poly-G trimming unless chosen explicitly.
- Optional FastQ Screen and Bowtie2 contamination panels, built from FASTA sources and indexed by the workflow. Screening never removes reads or excludes samples.
- Stand-alone raw read QC command that needs no reference and no statistical design.
- Preprocessing and quantification run as separate, resumable steps, with read-loss records and an early MultiQC report.

### Changed

- The core runtime adds DESeq2 1.50.2 (580 packages).
- SLURM profiles wait 120 seconds for worker outputs. On shared filesystems that cache directory listings, a resumed job could be reported as missing its outputs.

### Checks

- `scripts/qualify_deployment.sbatch` checks SLURM-worker execution, resume, a shared reference cache and a stock container.
- Read screening reports matches; it does not establish how sensitive a panel is.

## 2.0.0 (2026-09-23)

V2 releases the currently implemented **annotated, non-UMI bulk RNA-seq workflow**. The broader multi-assay roadmap is not part of this release's supported scope.

### Added

- Direct Excel setup, stable field IDs, source workbook/checksum and cell-level provenance.
- Samples, Libraries and Reads tabs with explicit ownership of inline versus external metadata.
- Executable canonical records: multiple sequencing runs/lanes within one library per sample, mate-ID/count and FASTQ validation, optional SHA256/size verification and retained read provenance.
- Recovery when canonical read provenance is missing; text IDs including leading zeros and literal `NA` survive setup and analysis.
- Generic SLURM validation and explicit SSH-host configuration for the optional uploader.
- README workflow diagram, complete installation/intake/run/resume guidance and capability boundaries.

### Changed and repaired

- Storage estimates warn without blocking launch or reducing concurrency. Estimates include complete libraries and retained merged and trimmed reads.
- PCA batch diagnostics preserve the declared scientific model and record review suggestions.
- Complete canonical-table relationships are validated; duplicate roles are retained for diagnostics instead of overwritten.
- Empty shell arguments are quoted correctly in both legacy and canonical execution.
- Environment verification uses the existing 579-package Linux lock. No new scientific runtime dependencies were added in this upgrade.

### Qualification and limits

The strict suite has 22 checks, including real synthetic FASTQ-to-report processing, planted DE effects, quantitative count agreement, WGCNA, enrichment, recovery and project isolation. Release publication requires passing scientific GitHub CI, including the separate public-data smoke test. See [release evidence](docs/RELEASE.md).

TAG-seq has no qualified kit; dream lacks complete end-to-end mixed-model qualification; live KEGG/annotation services have separate limitations. UMI, small-RNA, cell/nucleus, spatial, long-read, dual-organism and specialized routes remain unimplemented. See [the capability table](README.md#what-you-can-run).

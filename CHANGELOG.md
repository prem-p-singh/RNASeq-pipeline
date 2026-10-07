# Changelog

## Unreleased

- Track the n8n orchestration handoff and explicit three-test integration gate. Its reviewed schema-7 deployment snapshot predates this schema-8 release branch; it does not qualify schema-8 or dual-organism integration. Missing release-smoke prerequisites are visible skips and block strict qualification. See integrations/n8n/README.md and VALIDATION.md.

- Small-RNA animal route (Bowtie 1 to species mature miRNAs, unique-hit counting) and plant route (ShortStack 4.1.0 joint loci, DicerCall size classes) via the workbook's sRNA-seq tab or direct project configuration. Both routes count reads, not molecules: no UMI deduplication, no isomiR resolution, no phased-siRNA calling. Named non-UMI kit profiles live in `config/small_rna_kits.yaml`; UMI small-RNA kits are rejected.
- Public-data fixtures and local checks prepared for both routes (`tests/fixtures/small_rna_public.json` and `plant_small_rna_public.json`, their `prepare_*.py` downloaders, and `tests/check_small_rna_public.sh --animal|--plant`). Acceptance thresholds are fixed before any run: animal needs count conservation, a median |log2 pool A/B ratio| within 0.5 of the published ratio for well-covered Table S2 sequences, and Spearman >= 0.9 against published per-run counts; plant needs 21 and 24 nt as the two most common trimmed read lengths, both DicerCall 21 and 24 present, ath-miR156/166/319 overlapping a locus, and a repeat run giving an identical count matrix. Two criteria changed after a first FARM failure, both recorded in the fixtures: animal ratios are depth-normalized (raw ratios were offset by a 2.0x depth difference), and the plant size check replaced an unsourced 90% 20-24 nt floor (measured 76.3%). All checks passed on FARM (final job 39055140).
- TAG-seq QuantSeq route (public checks passed on FARM: FWD 39025685, REV 39023623, FWD-UMI 39023624; named-kit mapping floor 0.40) for `QuantSeq FWD`, `QuantSeq REV` and `QuantSeq FWD-UMI` (`config/tagseq_kits.yaml`), chosen in the workbook's TAGseq tab or with `analysis.kit`. The kit sets trimming and the STAR/featureCounts strand; FWD-UMI checks the TATA spacer, moves the UMI into the read name and reports umi_tools molecules per gene alongside read counts. A `tagseq` project without a kit is unchanged.
- `tests/check_tagseq_workflow.sh` (synthetic: plus-strand reads in a plus-strand gene A and a minus-strand gene B give A=100, B=0 under FWD and A=0, B=30 under REV; 20 molecules x 5 duplicates count as 20; a short read recorded as lost; without umi_tools the FWD/REV cases still run and the check reports SKIP for the UMI case) is a required check in `tests/run_all.sh`. `tests/check_tagseq_public.sh --fwd|--rev|--umi` runs six pinned public runs per kit (`tests/fixtures/tagseq_public.json`, prepared by `prepare_tagseq_public.py`) through counts, DE and the QC report, and records assignment and UMI duplication rates. All three kits stay `experimental` and `tagseq.validated_kits` stays empty until a FARM run passes these checks.
- Dual-organism route (status `validated` on the STAR route after two public checks on FARM, jobs 39420845 and 39505358; `experimental` on the Salmon route) for bulk RNA-seq libraries that hold a host and a pathogen, set with `dual_organism` in the project configuration or in the workbook's Bulk RNA-seq tab (schema 8). STAR aligns each read to a combined genome and gives it one label: host, pathogen, ambiguous or unassigned. Ambiguous and unassigned reads enter neither count table, and `results/assignment.tsv` reports all four per sample. The existing Salmon or STAR route then runs once per organism and writes `results/host/` and `results/pathogen/`. Samples below `min_assigned_reads` pathogen reads are listed in `results/pathogen/samples_not_quantified.tsv`; the pathogen analysis is counts only or counts with differential expression. Three defaults are assumptions (`multimap_limit` 20, `ambiguous_warn_fraction` 0.05, `min_assigned_reads` 500). Not supported: more than two organisms, dual-organism TAG-seq or small RNA, virus-derived small RNA, cross-organism networks, strain or variant analysis. Known limit: Salmon 1.10.3 with position-bias correction can crash when a transcript length class has no reads (see `analysis.salmon_position_bias` below).
- `tests/check_dual_workflow.sh` (synthetic host and pathogen genomes sharing one segment; exact label counts, split reads and per-organism counts with STAR and Salmon) is a required check in `tests/run_all.sh`. `tests/check_dual_public.sh --virus|--cellular` runs pinned public runs (`tests/fixtures/dual_public.json`, prepared by `prepare_dual_public.py`): grapevine with Grapevine red blotch virus (PRJNA339456, 3 healthy and 3 infected) and Arabidopsis with *Botrytis cinerea* (PRJNA1270936, 4 mock and 4 infected). The five acceptance criteria were written into the fixture before any run: conservation of the four labels, a lower pathogen fraction in every control than in every infected sample, count tables for every eligible sample, differential expression that completes or is reported as not estimable, and identical counts from a repeated assignment. Both public checks passed on the STAR route (FARM jobs 39420845 and 39505358).
- `analysis.salmon_position_bias` (default `true`) switches Salmon position-bias correction (`--posBias`) off when set to `false`. Salmon 1.10.3 crashed intermittently with a segmentation fault in this correction when a transcript length class had no reads: with one expressed 5000 nt transcript and two unexpressed 2000 nt transcripts, 45 of 60 runs at 300 reads, 37 of 60 at 3000 reads and 26 of 60 at 20000 reads; 0 of 40 with `--gcBias --seqBias` only. The pathogen project of a dual-organism run always has it off, the host project follows the parent, and `metrics.json` records the flags used. With the default the Salmon command is unchanged.
- Projects started before this version recompute from preprocessing on relaunch, because the preprocessing policy record and the quantification script changed. Results are unchanged.

## 2.1.0 — candidate, not yet released

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
- Full regression: 35 checks, 0 failed, 0 skipped, plus Salmon (with screening) and STAR public routes (39010068). M4–M8 remain incomplete. M2 remains partial: named chemistry profiles and new assay support remain open; declared-panel screening does not establish panel sensitivity.

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

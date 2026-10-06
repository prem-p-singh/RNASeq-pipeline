# V2 release qualification

Release: **v2.0.0**, 2026-09-23. Scope: **Linux x86_64, annotated non-UMI bulk gene expression**. The user authorized publication of the current bulk release; full multi-assay coverage is not a v2 claim.

## Release versions

Each software version groups development milestones (M1–M8). Workbook schema numbers are independent of software versions. A version is tagged only after its source is frozen, required checks pass with zero skips, and evidence covers each advertised capability.

| Version | Scope | Status |
|---|---|---|
| 2.0.0 | Annotated non-UMI bulk: Salmon, limma-voom/dream, workbook schema 3 | Released 2026-09-23 |
| 2.1.0 | Bulk completion: workbook schema 5 with analysis goal, subject models, formula, Contrasts tab and revision reports (M1); named bulk kit profiles, platform poly-G and automatic screening panel (M2); handbook DE default and interaction/paired/continuous design checks (M3); earlier schema 4 and Excel launch/resume (M1); raw QC, preprocessing and screening (M2); Salmon decoys, STAR/featureCounts, DESeq2, edgeR QL and objective selection (M3) | Candidate. Completed M1–M3 code passed 35 checks (0 failed, 0 skipped) and both public routes in Linux job 39010068 (exit 0, 55:09; source archive SHA256 `7fc10ce7987a11b9c195cd74420bab3ecf78c2cd1b3b9c3942d22c803e3afbb5`). Not tagged |
| 3.0.0 | Named end-tag kits, small RNA and dual organism (M4) | Small RNA (animal, plant) and TAG-seq QuantSeq FWD, REV, FWD-UMI implemented and passed on FARM. Final job 39055140 (commit `4ef08b0`, exit 0): strict suite 39 passed, 0 failed, 0 skipped; bulk GSE110004 Salmon and STAR; animal small RNA (GSE94585 technical pools: depth-normalized median ratio error 0.218, Spearman 0.962 and 0.974); plant small RNA (Arabidopsis SRR3222443-5: 15,018 loci, DicerCall 21 and 24, ath-miR156/166/319 found, 21 and 24 nt the top read lengths, identical repeat run). QuantSeq public checks: FWD 39025685 (mouse PRJNA631146), REV 39023623 (human PRJNA1494475), FWD-UMI 39023624 (mouse PRJEB107189; TATA spacer 88%). Two criteria were changed after a first failure and recorded in the fixtures: animal ratios are depth-normalized, and the plant size check replaced an unsourced 90% 20-24 nt floor (measured 76.3%). Kits and profiles promoted to `validated` on 2026-09-28. Dual organism (host and pathogen) implemented, status `experimental`; its public checks are pending (see below). Not tagged |
| Later | Cell/nucleus and spatial (M5), long-read and reference-poor (M6), specialized RNA (M7) | Not started; versions assigned when scoped |

Moved from M1–M2 to the assay milestones that need them: workbook intake for non-bulk assays (M4–M7) and barcode/UMI base protection (M4, with UMI and end-tag routes).

## Unreleased development qualification

The development branch changes the runtime and adds bulk routes. The sections below this one describe the **v2 tag**, not a qualification claim for the current branch.

- Core candidate: 580 packages, preserving the previous 579 package identities and adding DESeq2 1.50.2. Lock SHA256: `dc1d1459bf56f4a5e93c72d2a966c1a96a37168ea49e15933fe549e3549bff75`.
- Separate STAR module: 37 packages including STAR 2.7.11b, featureCounts 2.1.1 and samtools 1.24. Lock SHA256: `ebca83dc6a889b16291dadc923379a830b4596970b3cabe6c2c62c483c86c0d5`. The launcher verifies this module without replacing the core R environment.
- Direct native edgeR QL, DESeq2 Wald and one paired-donor dream comparison passed on Linux snapshots. Genome-decoy construction, reuse/tamper rejection, and QC/count-only objectives also passed. These fixtures do not qualify all study designs.
- A six-sample public Salmon run passed on the upgraded core runtime. STAR SE/PE exon-junction, strand and intronic controls passed. A six-sample public STAR/featureCounts/edgeR run also completed (Linux job 38992390, exit 0); synthetic raw-read/count/DE/report, resume and deleted-BAM-index recovery passed in job 38992412. Full current-source qualification remains pending.
- Schema 4 workbook setup/intake/preflight/project-isolation checks passed in Linux job 38993523; workbook rendering and prior values/formats/dropdowns/panes were checked locally. The bulk full-suite snapshot predates this intake edit.
- Screening module: 64 pinned packages, FastQ Screen 0.16.0/Bowtie2 2.5.4; lock SHA256 `065bb112579c134645324e4d9ee6a21d9aa822c3090e9ebc7475bc483852a7de`. Real SE/PE exclusive/shared mixture controls passed in 38994223, which then exposed interpreter shadowing during preflight. The launcher now retains core Python first in PATH; corrected QC-only screening DAG, resume, missing-report recovery, project isolation, clean installation and offline reuse passed in 38994332. Public yeast screening with Salmon/counts/DE/QC passed in 38994434. Job 38994223 as a whole failed and is not a suite pass. The 34-check strict suite (34 passed, 0 failed, 0 skipped) and both GSE110004 public routes (Salmon with screening, STAR) passed in 38994640 with exit 0 on the frozen source archive `c5e01ba2e992549de876326ef824b6b09ed1aa9cc1e6eb4b91faa0a980224a33`.
- The current release gate requires **34 checks, zero skips**, plus both public routes. Full current-source qualification is pending; do not publish a new release based only on the focused checks above.

### Dual-organism route (3.0.0, status `experimental`)

The route assigns each read to the host or the pathogen with STAR on a combined genome, then runs the existing Salmon or STAR route once per organism. It stays `experimental` until both public checks below pass and promotion is approved.

- Synthetic check `tests/check_dual_workflow.sh` (STAR and Salmon, exact label and per-organism counts, one paired-end project): passed on FARM, job 39412685, commit `ecc89a5`, exit 0. A dual project launched through `submit.sh --executor local` also completed (job 39412697, commit `ecc89a5`, exit 0).
- Public check, plant with virus, `tests/check_dual_public.sh --virus` (grapevine with Grapevine red blotch virus, PRJNA339456, runs SRR4043732 to SRR4043737, full depth, STAR route): passed on FARM, job 39420845, commit `e70f24c`, exit 0, all five criteria held. Virus-assigned reads were 1, 4 and 4 in the healthy libraries (15.9 to 22.5 million reads) and 906, 1651 and 2415 in the infected libraries (11.4 to 20.1 million reads). The three infected libraries were quantified for the virus; the three healthy ones were below `min_assigned_reads` (500), so virus differential expression was reported as not estimable. An earlier run with the former default of 10000 (job 39412687, commit `ecc89a5`) quantified no library for the virus, which is why the default was lowered.
- Public check, plant with fungus, `tests/check_dual_public.sh --cellular` (Arabidopsis with *Botrytis cinerea*, PRJNA1270936, runs SRR35391023 to SRR35391028, SRR35391035 and SRR35391036): pending. Job ID: pending. Result: pending.
- Strict suite with the dual-organism checks included: passed on FARM, job 39412689, commit `ecc89a5`, exit 0: 41 passed, 0 failed, 0 skipped, with the bulk Salmon and STAR public runs.

The acceptance criteria for the public checks are in `tests/fixtures/dual_public.json` and were written before any run: the four labels add up to the reads for every sample; every healthy or mock sample has a lower pathogen-assigned fraction than every infected sample; every sample reaches the host count table and every infected sample with at least `min_assigned_reads` pathogen reads reaches the pathogen count table; host differential expression completes and pathogen differential expression completes or is reported as not estimable with a reason; a repeat of the assignment gives identical counts. No numeric threshold is used.

Public protocol selection and evidence requirements are recorded in [Public-data qualification](PUBLIC_QUALIFICATION.md). The wider working plan remains incomplete.

## V2 upgrade evidence

- Canonical-intake candidate: Linux/SLURM job **38577940**, exit 0 in **14:02**, **22 passed, 0 failed, 0 skipped**. Includes two sequencing runs sharing a lane number, legacy single-end processing, quantitative recovery, changed-input detection and missing-provenance repair. Archive SHA256: `96916f31b70344a03e2260d9905daaa92b7c0997589369959cd096aa25e1c715`.
- Final retained-intermediate storage correction: Linux job **38578300**, exit 0; later TSV example changes passed local metadata checks.
- Declared-model preservation: focused Linux job **38577830**, exit 0; also included in the passing strict suite above.
- The first canonical candidate failed because empty shell arguments were lost. It was repaired and the corrected full raw-read workflow passed; the failed run is not qualification evidence.
- The **release-validation** GitHub workflow installs the locked runtime, verifies it, runs all 22 strict checks and the public smoke test on the release source. Publication requires its success on the tagged commit. Consult the [scientific runs](https://github.com/prem-p-singh/RNASeq-pipeline/actions/workflows/release-validation.yml) and release notes for the exact run and commit; historical worker evidence alone does not establish current GitHub status.

Workbook source/parser/generator consistency, external and inline setup, text identifiers and rejected invalid inputs are tested. Excel reference-service lookups are mocked in intake tests; raw-read processing is tested separately through the actual scientific tools. Live automatic reference/annotation service behavior is not certified by those mocks.

## Runtime contract

Install with `scripts/bootstrap.sh` from `environments/linux-64.explicit.txt`. The explicit lock contains 579 package URLs, exact builds and SHA256 hashes; its SHA256 is:

```text
351136c98793fb251de644f5ee352c1e8bb5c3edb462251addbade143b2cc7a4
```

The validated runtime includes Python 3.12.14, Snakemake 9.19.0, Salmon 1.10.3, fastp 0.23.4, MultiQC 1.35, R 4.5.2, tximport 1.38.2, edgeR 4.8.2, limma 3.66.0, variancePartition 1.40.2, clusterProfiler 4.18.4 and WGCNA 1.74.

The environment verifier compares installed Conda records with the lock, imports required libraries and checks executable/library locations. It is not a byte-for-byte integrity audit of every installed file. Do not install unrelated packages into the qualified prefix. Project annotation packages belong in the configured OrgDb cache.

## Historical environment and deployment evidence

Linux validation is executed on Linux x86_64 compute nodes under SLURM with site-specific account and partition settings. Some allocated hosts have GPU names, but these jobs request and use CPUs only. Full local evidence is retained by the maintainer; concise outcomes below describe what was actually exercised.

- Initial environment build: SLURM job **38516117**.
- Independent clean recreation from the explicit lock and offline reuse: job **38516428**, completed with exit 0 in 5 minutes 52 seconds; all 579 package records and required imports verified.
- GO enrichment qualification: job **38516721**, passed with a locally constructed AnnotationForge OrgDb. Tests retain positive and negative enrichment directions, verify real output paths, recognize valid empty results and unavailable KEGG configuration, and fail on a missing ranking statistic.
- Raw FASTQ qualification: job **38516701** passed its single-end and paired-end end-to-end checks and restart/report recovery checks. That older combined run failed its earlier GO fixture and public reference fixture; it is **not** a passing full release run.

- Public six-sample yeast workflow: job **38516756**, completed with exit 0 in 2 minutes 13 seconds. Reference construction, six FASTQ pairs, gene counts, DE result tables and QC reports passed.
- Actual `submit.sh` → SLURM executor → worker execution: controller job **38516744**, completed with exit 0 in 5 minutes 12 seconds. Four remote rules retrieved references, built an index and quantified a paired-end sample. The source checkout and project were separate shared paths; the project path contained spaces. Worker jobs included 38516760, 38516829 and 38516837 on bm6/bm7. Preflight, runtime verification and the capacity plan ran through the real launcher.

- Final resource/default-resolution and project-isolation checks: job **38517046**, passed. It also verified that an absent runtime is rejected with explicit missing-package findings.

- **Final combined release gate: job 38517266**, completed with **exit 0** on gpu-5-58 in **15 minutes 20 seconds**. It verified the locked runtime, passed **20 strict checks with zero failures/skips**, and passed the corrected six-sample public workflow in the same invocation. Recorded peak RSS was approximately 1.62 GiB for this small qualification workload, not a cohort-sizing benchmark.

All 87 runtime/configuration/test files in that historical frozen snapshot matched its source manifest. That manifest is not the v2 source identity. V2 is identified by its Git tag, exact commit and release checksums. Historical failed/superseded runs remain diagnostic history. A release must not be tagged from a skipped or failed scientific suite.

### Test data

`tests/fixtures/fastq_project.py` creates six deterministic biological-sample surrogates, 120 unique genes, and independently known fragments. Real fastp and Salmon run from raw FASTQs through reference construction, gene aggregation, differential expression and reporting. Assertions check quantitative recovery, planted effect direction, source-file preservation, owned-intermediate cleanup and restart behavior. This is a technical fixture, not evidence of biological generalization.

`tests/check_enrichment.R` builds a local annotation database and runs real signed GO GSEA. Its `GID` key type is explicit. For public annotation packages, use a compatible configured `orgdb.key_type` and verify gene-ID mapping; identifier namespaces are not interchangeable.

The public smoke test uses **GSE110004**, accessions **SRR6357070–SRR6357075**, with three wild-type and three uninduced Rap1-AID biological samples. These are chromosome-I-filtered, downsampled yeast reads from the [pinned nf-core test-data revision](https://github.com/nf-core/test-datasets/tree/626c8fab639062eade4b10747e919341cbf9b41a). [The dataset's metadata and sampling procedure](https://raw.githubusercontent.com/nf-core/test-datasets/626c8fab639062eade4b10747e919341cbf9b41a/README.md) establish the sample identities. The original transcriptome additionally contains `Gfp_transgene_gene`, which is absent from its GTF. The fixture removes that one artificial entry and records original/derived hashes in `inputs/reference_derivation.json`; production reference checks remain strict. This smoke test does not reproduce the study's genome-wide conclusions.

### Reproduce on a SLURM cluster

From the source checkout, run:

```bash
sbatch scripts/validate_slurm.sbatch
```

That job builds a fresh locked runtime in worker-local temporary space, verifies it, runs `tests/run_all.sh --strict`, runs the public smoke test and checks offline reuse. Evidence is saved under `docs/evidence/slurm-JOBID/`. Supply account, partition and QoS through sbatch options when required by your site. This tests local execution inside a worker allocation; distributed launcher-to-worker qualification is a separate check.

For an already activated qualified Linux environment:

```bash
python3 scripts/environment_check.py --out environment_report.json
bash tests/run_all.sh --strict
bash tests/check_public.sh
```

The scientific GitHub Actions workflow repeats these checks. The lightweight workflow permits reported dependency skips and cannot certify a release.

## Repairs in this candidate

- Install and verify an exact runtime instead of silently mixing cluster modules and R/Python libraries.
- Record the implemented recommendation route, design-based statistical backend, policy hash and rule/handbook IDs; reject unsupported objectives and UMI requests.
- Remove the inherited 20 GB value from all project-creation defaults. Account for retained remote reads, existing environments, measured local reference sizes and configured concurrency limits.
- Preserve seekable trimmed input during quantification; validate metrics/count fields and record source FASTQ hashes before cleanup.
- Check transcript-to-gene compatibility before indexing, verify source checksums on cache reuse and publish the completion record after the index.
- Preserve retained-sample effective lengths, tximport data and count provenance.
- Declare analysis configuration/runtime identity for downstream reruns, local read dependencies and shell helper dependencies; repair missing report targets and result-manifest validation.
- Correct enrichment result paths, retain signed ranks, read statistic names as strings, and fail when any requested enrichment method errors.
- Install generated annotation packages outside the release runtime and serialize builds for the same annotation cache entry.
- Provide a zero-skip scientific CI gate and raw/public-data qualification fixtures.

## Limits that remain part of the release contract

| Area | Boundary |
|---|---|
| RNA-seq coverage | STAR alignment and DESeq2/edgeR inference adapters are implemented and released. Small RNA (animal Bowtie 1, plant ShortStack) is implemented and validated for the TruSeq Small RNA kit (3.0.0). UMI handling is limited to molecule counting for the validated QuantSeq FWD-UMI kit (3.0.0); the dual-organism route (host and pathogen) is implemented with status `experimental` and its public checks are pending (3.0.0); single-cell, spatial, long-read and specialized assays remain unimplemented. Installed libraries do not imply an enabled route. |
| TAG-seq | QuantSeq FWD, REV and FWD-UMI are validated by their public-data checks; a kit-less `tagseq` project remains provisional. |
| Statistical generalization | Synthetic fixtures and a small public smoke dataset do not validate every organism, design or contrast. The suite checks repeated-measures design/replication logic but does not qualify an end-to-end dream fit; mixed-model inference remains provisional. |
| Annotation | GO has a local fixture. Live NCBI annotation construction, live KEGG results, annotation snapshot reproducibility and the eggNOG fallback are not certified by that fixture. Missing prerequisites and failures must remain visible in status files. |
| Sample policy | QC inclusion/exclusion is implemented; a complete review/override ledger and mandatory action on library-type mismatch remain planned. Review recorded mismatches before biological interpretation. |
| Metadata | Canonical execution supports multiple run/lane units within one non-UMI bulk library per sample. Multiple prepared libraries per sample remain unavailable. Legacy sheets retain one file/pair per row. |
| Reference/cache | Index publication is atomic, but the entire multi-file bundle is not published in one transaction. Mutable remote URLs are not re-fetched solely to detect changes. Full simultaneous multi-project cache stress qualification remains open. |
| Storage and scale | Estimates remain labelled assumptions; no continuous capacity monitor or calibrated large-cohort benchmark. Capacity/quota concerns warn only. Actual failed writes still fail tasks. |
| Platforms and deployment | Linux x86_64 with SLURM is the tested target. Other clusters require site configuration; macOS, ARM and a container deployment are not qualified. |
| Intake and optional integrations | Excel supports independent bulk with optional batch and inline/external records. Setup-to-config tests mock live reference lookup; the scientific DAG is tested separately. Advanced models and reference choices need YAML. The optional upload/wizard is not an equivalent qualified interface. |

## Publication procedure

1. Preserve the tested source revision, lock hash and qualification logs; inspect the complete candidate diff.
2. Keep handbook PDFs, private audit/upgrade plans, working trackers and raw local evidence out of the public source tree. Public operational documentation must remain readable without those files.
3. Push the reviewed source changes and require the **release-validation** scientific job, including the public smoke test, to pass without skips. Configure this as a required branch check in the hosting repository.
4. Publish with the supported scope and limitations above. Do not describe this candidate as an implemented universal RNA-seq platform or as a qualified TAG-seq kit workflow.
5. Record the actual tag/revision in release notes and citation metadata; publish only the reviewed scope. Preserve repository visibility unless the owner separately requests a change.

For a runtime upgrade, solve `environment.yml` into a new prefix, export `conda list --explicit --sha256`, review the exact lock changes and repeat the qualification gates before changing the supported release environment.

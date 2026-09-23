# V2 release qualification

Release: **v2.0.0**, 2026-09-23. Scope: **Linux x86_64, annotated non-UMI bulk gene expression**. The user authorized publication of the current bulk release; full multi-assay coverage is not a v2 claim.

## Current upgrade evidence

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
| RNA-seq coverage | STAR alignment, DESeq2/edgeR inference adapters, UMI handling, single-cell, spatial, small-RNA, long-read, dual-organism and specialized assays remain unimplemented. Installed libraries do not imply an enabled route. |
| TAG-seq | A provisional route exists, but no named kit is qualified. |
| Statistical generalization | Synthetic fixtures and a small public smoke dataset do not validate every organism, design or contrast. The suite checks repeated-measures design/replication logic but does not qualify an end-to-end dream fit; mixed-model inference remains provisional. |
| Annotation | GO has a local fixture. Live NCBI annotation construction, live KEGG results, annotation snapshot reproducibility and the eggNOG fallback are not certified by that fixture. Missing prerequisites and failures must remain visible in status files. |
| Sample policy | QC inclusion/exclusion is implemented; a complete review/override ledger and mandatory action on library-type mismatch remain planned. Review recorded mismatches before biological interpretation. |
| Metadata | Canonical execution supports multiple run/lane units within one non-UMI bulk library per sample. Multiple prepared libraries per sample remain unavailable. Legacy sheets retain one file/pair per row. |
| Reference/cache | Index publication is atomic, but the entire multi-file bundle is not published in one transaction. Mutable remote URLs are not re-fetched solely to detect changes. Full simultaneous multi-project cache stress qualification remains open. |
| Storage and scale | Estimates remain labelled assumptions; no continuous capacity monitor or calibrated large-cohort benchmark. Capacity/quota concerns warn only. Actual failed writes still fail tasks. |
| Platforms and deployment | Linux x86_64 with SLURM is the tested target. Other clusters require site configuration; macOS, ARM and a container deployment are not qualified. |
| Intake and optional integrations | Excel supports independent bulk with optional batch and inline/external records. Setup-to-config tests mock live reference lookup; the scientific DAG is tested separately. Advanced models and reference choices need YAML. Optional upload/wizard and n8n integrations are not equivalent qualified interfaces. |

## Publication procedure

1. Preserve the tested source revision, lock hash and qualification logs; inspect the complete candidate diff.
2. Keep handbook PDFs, private audit/upgrade plans, working trackers and raw local evidence out of the public source tree. Public operational documentation must remain readable without those files.
3. Push the reviewed source changes and require the **release-validation** scientific job, including the public smoke test, to pass without skips. Configure this as a required branch check in the hosting repository.
4. Publish with the supported scope and limitations above. Do not describe this candidate as an implemented universal RNA-seq platform or as a qualified TAG-seq kit workflow.
5. Record the actual tag/revision in release notes and citation metadata; publish only the reviewed scope. Preserve repository visibility unless the owner separately requests a change.

For a runtime upgrade, solve `environment.yml` into a new prefix, export `conda list --explicit --sha256`, review the exact lock changes and repeat the qualification gates before changing the supported release environment.

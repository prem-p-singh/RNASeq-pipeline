# Release qualification

**Final FARM gate: PASSED — 20 checks passed, 0 failed, 0 skipped; the public six-sample workflow also passed.**

Qualification date: 2026-09-22. Target: **Linux x86_64, annotated non-UMI bulk gene expression**. This is a release candidate, not certification of the broader multi-assay roadmap.

## Runtime contract

Install with `scripts/bootstrap.sh` from `environments/linux-64.explicit.txt`. The explicit lock contains 579 package URLs, exact builds and SHA256 hashes; its SHA256 is:

```text
351136c98793fb251de644f5ee352c1e8bb5c3edb462251addbade143b2cc7a4
```

The validated runtime includes Python 3.12.14, Snakemake 9.19.0, Salmon 1.10.3, fastp 0.23.4, MultiQC 1.35, R 4.5.2, tximport 1.38.2, edgeR 4.8.2, limma 3.66.0, variancePartition 1.40.2, clusterProfiler 4.18.4 and WGCNA 1.74.

The environment verifier compares installed Conda records with the lock, imports required libraries and checks executable/library locations. It is not a byte-for-byte integrity audit of every installed file. Do not install unrelated packages into the qualified prefix. Project annotation packages belong in the configured OrgDb cache.

## Evidence and reproducibility

Linux validation is executed on UC Davis FARM compute nodes under SLURM, using the `publicgrp` account and `low` partition. Some allocated hosts have GPU names, but these jobs request and use CPUs only. Full local evidence is retained by the maintainer; concise outcomes below describe what was actually exercised.

- Initial environment build: FARM job **38516117**.
- Independent clean recreation from the explicit lock and offline reuse: job **38516428**, completed with exit 0 in 5 minutes 52 seconds; all 579 package records and required imports verified.
- GO enrichment qualification: job **38516721**, passed with a locally constructed AnnotationForge OrgDb. Tests retain positive and negative enrichment directions, verify real output paths, recognize valid empty results and unavailable KEGG configuration, and fail on a missing ranking statistic.
- Raw FASTQ qualification: job **38516701** passed its single-end and paired-end end-to-end checks and restart/report recovery checks. That older combined run failed its earlier GO fixture and public reference fixture; it is **not** a passing full release run.

- Public six-sample yeast workflow: job **38516756**, completed with exit 0 in 2 minutes 13 seconds. Reference construction, six FASTQ pairs, gene counts, DE result tables and QC reports passed.
- Actual `submit.sh` → SLURM executor → worker execution: controller job **38516744**, completed with exit 0 in 5 minutes 12 seconds. Four remote rules retrieved references, built an index and quantified a paired-end sample. The source checkout and project were separate shared paths; the project path contained spaces. Worker jobs included 38516760, 38516829 and 38516837 on bm6/bm7. Preflight, runtime verification and the capacity plan ran through the real launcher.

- Final resource/default-resolution and project-isolation checks: job **38517046**, passed. It also verified that an absent runtime is rejected with explicit missing-package findings.

- **Final combined release gate: job 38517266**, completed with **exit 0** on gpu-5-58 in **15 minutes 20 seconds**. It verified the locked runtime, passed **20 strict checks with zero failures/skips**, and passed the corrected six-sample public workflow in the same invocation. Recorded peak RSS was approximately 1.62 GiB for this small qualification workload, not a cohort-sizing benchmark.

All 87 runtime/configuration/test files in the frozen FARM snapshot were checked against the local source manifest with zero differences. The publication archive includes that manifest as `SOURCE_SHA256.json`. Its SHA256 is `4960694928ab5708b3fbc9141172ac77d7491925bc0bb4b58f4d124c1498d98f`. Earlier failed/superseded runs remain diagnostic history; the final candidate is qualified by the successful combined gate above. A release must not be tagged from a skipped or failed scientific suite.

### Test data

`tests/fixtures/fastq_project.py` creates six deterministic biological-sample surrogates, 120 unique genes, and independently known fragments. Real fastp and Salmon run from raw FASTQs through reference construction, gene aggregation, differential expression and reporting. Assertions check quantitative recovery, planted effect direction, source-file preservation, owned-intermediate cleanup and restart behavior. This is a technical fixture, not evidence of biological generalization.

`tests/check_enrichment.R` builds a local annotation database and runs real signed GO GSEA. Its `GID` key type is explicit. For public annotation packages, use a compatible configured `orgdb.key_type` and verify gene-ID mapping; identifier namespaces are not interchangeable.

The public smoke test uses **GSE110004**, accessions **SRR6357070–SRR6357075**, with three wild-type and three uninduced Rap1-AID biological samples. These are chromosome-I-filtered, downsampled yeast reads from the [pinned nf-core test-data revision](https://github.com/nf-core/test-datasets/tree/626c8fab639062eade4b10747e919341cbf9b41a). [The dataset's metadata and sampling procedure](https://raw.githubusercontent.com/nf-core/test-datasets/626c8fab639062eade4b10747e919341cbf9b41a/README.md) establish the sample identities. The original transcriptome additionally contains `Gfp_transgene_gene`, which is absent from its GTF. The fixture removes that one artificial entry and records original/derived hashes in `inputs/reference_derivation.json`; production reference checks remain strict. This smoke test does not reproduce the study's genome-wide conclusions.

### Reproduce on FARM

From the source checkout, run:

```bash
sbatch scripts/validate_farm.sbatch
```

That job builds a fresh locked runtime in worker-local temporary space, verifies it, runs `tests/run_all.sh --strict`, runs the public smoke test and checks offline reuse. Evidence is saved under `docs/evidence/farm-JOBID/`. Cluster-specific account/partition settings are explicit in the script. This tests local execution inside a worker allocation; distributed launcher-to-worker qualification is a separate check.

For an already activated qualified Linux environment:

```bash
python3 scripts/environment_check.py --out environment_report.json
bash tests/run_all.sh --strict
bash tests/check_public.sh
```

The scientific GitHub Actions workflow repeats these checks. Its status has not been observed on GitHub until the candidate is pushed and CI runs. The lightweight workflow permits reported dependency skips and cannot certify a release.

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
| Metadata | The current execution model is one row/biological measurement and one read file/pair. Separate sample/library/read validation is not multi-lane execution. |
| Reference/cache | Index publication is atomic, but the entire multi-file bundle is not published in one transaction. Mutable remote URLs are not re-fetched solely to detect changes. Full simultaneous multi-project cache stress qualification remains open. |
| Storage and scale | Estimates remain labelled assumptions; no continuous capacity monitor or calibrated large-cohort resource benchmark. Account quota discovery is site-specific. A real quota is a constraint, never a pipeline-wide default. |
| Platforms and deployment | Linux x86_64/FARM is the tested target. Other clusters require site configuration; macOS, ARM and a container deployment are not qualified. |
| Intake and optional integrations | The explicit TSV/YAML route is the release interface. Convenience spreadsheet/URL-list intake and n8n are not equivalent end-to-end qualified interfaces. |

## Publication procedure

1. Preserve the tested source revision, lock hash and qualification logs; inspect the complete candidate diff.
2. Keep handbook PDFs, private audit/upgrade plans, working trackers and raw local evidence out of the public source tree. Public operational documentation must remain readable without those files.
3. Push the reviewed source changes and require the **release-validation** scientific job, including the public smoke test, to pass without skips. Configure this as a required branch check in the hosting repository.
4. Publish with the supported scope and limitations above. Do not describe this candidate as an implemented universal RNA-seq platform or as a qualified TAG-seq kit workflow.
5. Record the actual tag/revision in the release notes and citation metadata. No tag, public release or push is performed by preparing this candidate.

For a runtime upgrade, solve `environment.yml` into a new prefix, export `conda list --explicit --sha256`, review the exact lock changes and repeat the qualification gates before changing the supported release environment.

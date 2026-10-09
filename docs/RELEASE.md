# V2 release qualification

Release: **v2.1.0**, qualified 2026-10-08 on source `0a857c6`; publication still requires the hosted workflow to pass on the tagged commit. Previous release: v2.0.0, 2026-09-23. Scope: **Linux x86_64, annotated non-UMI bulk gene expression**. The user authorized publication of the current bulk release; full multi-assay coverage is not a v2 claim.

## Release versions

Each software version groups development milestones (M1–M8). Workbook schema numbers are independent of software versions. A version is tagged only after its source is frozen, required checks pass with zero skips, and evidence covers each advertised capability.

| Version | Scope | Status |
|---|---|---|
| 2.0.0 | Annotated non-UMI bulk: Salmon, limma-voom/dream, workbook schema 3 | Released 2026-09-23 |
| 2.1.0 | Bulk completion: workbook schema 5 with analysis goal, subject models, formula, Contrasts tab and revision reports (M1); named bulk kit profiles, platform poly-G and automatic screening panel (M2); handbook DE default and interaction/paired/continuous design checks (M3); earlier schema 4 and Excel launch/resume (M1); raw QC, preprocessing and screening (M2); Salmon decoys, STAR/featureCounts, DESeq2, edgeR QL and objective selection (M3) | Qualified 2026-10-08 on source `0a857c6`: 38 checks, 0 failed, 0 skipped, both public routes, SLURM workers, resume, shared cache and a stock container (jobs 39620298, 39620305, 39618728). The 2.1.0 code alone had passed 35 checks in job 39010068. Not published: the hosted workflow has not run on this source |
| 3.0.0 | Named end-tag kits, small RNA and dual organism (M4) | Test data prepared only |
| Later | Cell/nucleus and spatial (M5), long-read and reference-poor (M6), specialized RNA (M7) | Not started; versions assigned when scoped |

Moved from M1–M2 to the assay milestones that need them: workbook intake for non-bulk assays (M4–M7) and barcode/UMI base protection (M4, with UMI and end-tag routes).

## 2.1.0 qualification

Version 2.1.0 changes the runtime and adds bulk routes. This section and [Deployment qualification](#deployment-qualification) record its evidence; the gate that qualifies the released source is under [Combined 2.1.0 source](#combined-210-source). The sections from [V2 upgrade evidence](#v2-upgrade-evidence) onward describe the **v2.0.0 tag**.

- Core candidate: 580 packages, preserving the previous 579 package identities and adding DESeq2 1.50.2. Lock SHA256: `dc1d1459bf56f4a5e93c72d2a966c1a96a37168ea49e15933fe549e3549bff75`.
- Separate STAR module: 37 packages including STAR 2.7.11b, featureCounts 2.1.1 and samtools 1.24. Lock SHA256: `ebca83dc6a889b16291dadc923379a830b4596970b3cabe6c2c62c483c86c0d5`. The launcher verifies this module without replacing the core R environment.
- Direct native edgeR QL, DESeq2 Wald and one paired-donor dream comparison passed on Linux snapshots. Genome-decoy construction, reuse/tamper rejection, and QC/count-only objectives also passed. These fixtures do not qualify all study designs.
- A six-sample public Salmon run passed on the upgraded core runtime. STAR SE/PE exon-junction, strand and intronic controls passed. A six-sample public STAR/featureCounts/edgeR run also completed (Linux job 38992390, exit 0); synthetic raw-read/count/DE/report, resume and deleted-BAM-index recovery passed in job 38992412. Full current-source qualification remains pending.
- Schema 4 workbook setup/intake/preflight/project-isolation checks passed in Linux job 38993523; workbook rendering and prior values/formats/dropdowns/panes were checked locally. The bulk full-suite snapshot predates this intake edit.
- Screening module: 64 pinned packages, FastQ Screen 0.16.0/Bowtie2 2.5.4; lock SHA256 `065bb112579c134645324e4d9ee6a21d9aa822c3090e9ebc7475bc483852a7de`. Real SE/PE exclusive/shared mixture controls passed in 38994223, which then exposed interpreter shadowing during preflight. The launcher now retains core Python first in PATH; corrected QC-only screening DAG, resume, missing-report recovery, project isolation, clean installation and offline reuse passed in 38994332. Public yeast screening with Salmon/counts/DE/QC passed in 38994434. Job 38994223 as a whole failed and is not a suite pass. The frozen-source 34-check strict suite and both public routes passed in Linux job 38994640: 34 passed, 0 failed, 0 skipped, exit 0, 1:22:45. Retained log: [screen-final-38994640.log](evidence/development-20260926/screen-final-38994640.log), SHA256 `ff1ccaccba982d272272c46899a8422180025ba2be7c9e169aab5931a9a3cc03`. The qualified commit is `7714b84`; subsequent commit `01190c1` adds documentation and an unexercised small-RNA fixture, with no production-code change. This qualification covers the tested fixtures only and does not close M1–M8; it does not qualify the current release gate.
- The current release gate requires **37 checks, zero skips**, plus both public routes. Linux job **39597462** passed it on source `82e70c5`: 37 passed, 0 failed, 0 skipped, then both public routes and offline reuse. See [Deployment qualification](#deployment-qualification) for the scope of that pass and for what still blocks publication.

## Deployment qualification

Source: commit `82e70c56572b7bcc64847eabd7380dd07ff44982` on `claude/m8-qualification`, three commits after `731f664`. `git archive` SHA256 `454609dd18fd3f649571f738f6fd35e0f146d109a61c36ef2ac666a0cf4981af`, verified on the cluster before extraction. Later commits on this branch change documentation only. All jobs ran on 2026-10-07 (Pacific) on Linux x86_64 under SLURM, with source, runtimes, projects and temporary files on one NFS filesystem.

| Check | Job | Result |
|---|---|---|
| Fresh installation and release gate | 39597462, exit 0, 1:30:03 | Core (580 packages, lock `dc1d1459...`), STAR (37, `ebca83dc...`) and screening (64, `065bb112...`) runtimes built from the locks in a new directory and verified. `tests/run_all.sh --strict`: **37 passed, 0 failed, 0 skipped**. GSE110004 Salmon with screening and GSE110004 STAR each passed six samples, counts, DE and QC report. All three runtimes were then reused with `CONDA_OFFLINE=true`. |
| Launcher to SLURM workers | 39597460, exit 0, 1:05:21 | `submit.sh --executor slurm` ran the six-sample public Salmon project from a controller job. 27 rule jobs completed as separate SLURM jobs; none ran in the controller allocation. |
| Interruption and resume | same job | The launcher received SIGINT at 16 of 28 steps with one sample quantified. A second `submit.sh` finished the project with 11 worker jobs. The finished sample kept its modification time, so it was not recomputed. |
| Two projects, one reference cache entry | same job | Two projects naming the same reference and `reference.cache_dir` were launched together, 27 worker jobs each. The cache held one entry and no build lock afterwards. Per-sample count totals differed by at most 0.023%; the 1% tolerance was written into the script before the run. |
| Stock Linux container | 39597461, exit 0, 0:26:18 | `condaforge/miniforge3:26.7.2-0` (Ubuntu 24.04.5, conda 26.7.2) under Apptainer 1.5.4 with `--containall --no-home`. Fresh 580-package core install verified, `tests/check_end_to_end.sh` and the public Salmon route passed, offline reuse passed. |

Log SHA256: `tests.log` `b0421898829713fa682551cfee315d9ef120899652478e8ebd40b7e887f218fc`, `public.log` `67aeeb28d5c384b5ffac9018f318c5a736c3f08199e0f468acee4a30aecc0fef`, `public-star.log` `31763084ca6b30a605bdd00e3a2e115219b3f5a1f9890acd00ef09c73eb9c1ca`, `distributed.log` `b0a185a4460709e353cb50b6e09242fc7a8230f28a7e3bb3e96bcde6584feb3c`, `restart-2.log` `ae82d1257c0fc89ee966a108c101cea534e3a76c84aca0af1398be3e279a560e`, `container.log` `0f89ede7d9827c7e7c66c0d260ebaa29703a93c24463ac10ff207863215f8e8b`.

The resume check found a defect. With `latency-wait: 30`, job 39597092 failed five quantification attempts on resume with `MissingOutputException` although the workers had written the files. The filesystem is mounted with `acdirmin=60`: after the controller removes the outputs of an interrupted or failed job, its NFS client can report the rewritten files as absent for 60 seconds or more. All three profiles now wait 120 seconds. Runs that never removed outputs did not show the failure.

Runs that are not evidence: 39595347 and 39597092 failed in the interruption step (the first on a signal sent to the wrong process by the qualification script, the second on the defect above). Jobs 39595476 and 39595346 passed the gate and the container check on `4f1f09d`, before the profile change, and were repeated on `82e70c5`.

Prerequisites of the three n8n checks were supplied from outside the checkout: Node.js 22.14.0, the reviewed schema-7 bundle (source SHA256 `c0927cf7...`, as recorded in `integrations/n8n/source.review.json`) and its verified small-RNA (89 packages) and UMI (119 packages) runtimes. Worker jobs took the SLURM account and partition from `SBATCH_ACCOUNT` and `SBATCH_PARTITION`; the profiles were not edited.

### Release matrix for this source

| Capability | Status | Evidence |
|---|---|---|
| Fresh locked installation and offline reuse, Linux x86_64 | Passed | 39597462; container 39597461 |
| Salmon bulk route, raw reads to report | Passed on synthetic fixtures and the downsampled public study | 39597462 |
| STAR/featureCounts bulk route | Passed on synthetic controls and the downsampled public study | 39597462 |
| limma-voom, edgeR QL, DESeq2 Wald, one paired-donor dream comparison | Passed on fixtures | 39597462 |
| Declared-panel read screening | Passed on known mixtures and the public study; panel sensitivity not established | 39597462 |
| `submit.sh` with the SLURM executor | Passed for the Salmon route, six samples | 39597460 |
| Resume after an interrupted launcher | Passed for SIGINT on the Salmon route. A killed controller with orphaned workers was not tested | 39597460 |
| Simultaneous projects on one reference cache | Passed for two projects. Larger numbers and the STAR index were not tested | 39597460 |
| Stock Ubuntu container userland | Passed for the core runtime and Salmon route. No container image is shipped; STAR, screening and SLURM submission from inside a container were not tested | 39597461 |
| n8n orchestration gate | Passed with external prerequisites | 39597462 |
| Hosted `release-validation` workflow | **Not run on this source.** It installs neither Node.js nor the reviewed bundle, so its strict step would stop on two skipped checks until the workflow supplies them | none |
| Second SLURM site, macOS, ARM | Not tested | none |
| Large cohorts, full-size genomes | Not tested | none |
| TAG-seq kits, small RNA, dual organism, single cell, spatial, long read, specialized assays | Not on this source | none |

This closes the deployment milestone for the bulk scope of this source. It does not complete the wider assay roadmap, and it is not a published release: publication still requires the hosted workflow to pass on the tagged commit.

### Combined 2.1.0 source

Source: commit `0a857c6` on `claude/release-2.1`: the deployment-qualified source above plus the 2.1.0 work (workbook schema 5, named bulk kit profiles, automatic screening panel, design checks; previously qualified alone in job 39010068 with 35 passed, 0 failed, 0 skipped). `git archive` SHA256 `4f1acc8c551bec5599aab235946178e3e7b5d0829ac0b0f2f20b6cca6b7ff92f`, verified on the cluster before extraction. The jobs in the table above do not cover this source; the jobs below do. All ran on 2026-10-08 (Pacific) on Linux x86_64 under SLURM with source, runtimes and projects on one NFS filesystem.

| Check | Job | Result |
|---|---|---|
| Fresh installation and release gate | 39620298, exit 0, 1:41:45 | Core, STAR and screening runtimes built from the locks in a new directory and verified. `tests/run_all.sh --strict`: **38 passed, 0 failed, 0 skipped**. GSE110004 Salmon with screening and GSE110004 STAR each passed six samples, counts, DE and QC report. The runtimes were then reused with `CONDA_OFFLINE=true`. |
| Launcher to SLURM workers | 39620305, exit 0, 1:24:08 | `submit.sh --executor slurm` ran the six-sample public Salmon project from a controller job; 27 rule jobs completed as separate SLURM jobs. The verified core runtime of job 39620298 was reused. |
| Interruption and resume | same job | The launcher was interrupted and a second `submit.sh` finished the project with 10 worker jobs. |
| Two projects, one reference cache entry | same job | Two projects naming the same reference and cache directory were launched together. Largest per-sample count-total difference 0.018%, inside the 1% tolerance written into the script before the run. |
| Stock Linux container | 39618728, exit 0, 0:25:46 | `condaforge/miniforge3:26.7.2-0` (Ubuntu 24.04.5, conda 26.7.2) under Apptainer: fresh core install, `tests/check_end_to_end.sh` and the public Salmon route passed. |

Log SHA256: `tests.log` `041e7220ecd597c082b79d499d6589b576828b88d648b8755248c51ec550e029`, `public.log` `24b2e3c1504d70242af769f70190359844b7c5f7270eac9a41a27c095df0e086`, `public-star.log` `579bf25615b591920a92e377fbd2de650bc4dd13e89ce5b37efd43ac82a0492a`, `distributed.log` `f646475f00743eb95be36494c0be3a770cd2e827fdb16899db99bc99e777e1e4`, `restart-2.log` `9cccebc9e95af33df1f201d899e9ba99d6b07f01fe587916c11a2fcbe32fd5d1`, `container.log` `ce941a75cd02888dcc02680be5a841ec31296ffc6c2e1dd7f6e14fec81576c9e`.

How these runs differ from the table above. Jobs 39620298 and 39620305 were submitted with `HOME` set to a directory on the shared filesystem, so that Conda's package cache could not fall back to the submitting account's home directory. The n8n prerequisites were again supplied from outside the checkout (Node.js 22.14.0, the reviewed schema-7 bundle, and its small-RNA and UMI runtimes). The release matrix above applies to this source with the gate count changed to 38; its limits are unchanged.

Runs that are not evidence: 39618727 started a fresh installation whose package downloads filled the submitting account's home quota; 39618726 was cancelled after it hit the same quota errors; 39620291 failed during installation with a Conda package-cache error, before any test.

The hosted `release-validation` workflow has not run on this source. This is not a published release.

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

That job builds a fresh locked runtime in worker-local temporary space, verifies it, runs `tests/run_all.sh --strict`, runs the public smoke test and checks offline reuse. Evidence is saved under `docs/evidence/slurm-JOBID/`. Supply account, partition and QoS through sbatch options when required by your site. Set `RNASEQ_VALIDATE_ROOT` to build the runtime on a shared filesystem instead of worker-local `/tmp`. This tests local execution inside a worker allocation.

For launcher-to-worker execution, resume and a shared reference cache, and separately for a stock container:

```bash
RNASEQ_QUALIFY_ROOT=/shared/dir sbatch scripts/qualify_deployment.sbatch
RNASEQ_QUALIFY_ROOT=/shared/dir sbatch scripts/qualify_deployment.sbatch container
```

`RNASEQ_QUALIFY_ROOT` must be visible to every worker. Evidence is saved under `docs/evidence/deploy-MODE-JOBID/`.

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
| Reference/cache | Index publication is atomic, but the entire multi-file bundle is not published in one transaction. Mutable remote URLs are not re-fetched solely to detect changes. Two simultaneous projects building one Salmon cache entry passed (job 39597460); larger numbers of projects and the STAR index remain untested. |
| Storage and scale | Estimates remain labelled assumptions; no continuous capacity monitor or calibrated large-cohort benchmark. Capacity/quota concerns warn only. Actual failed writes still fail tasks. |
| Platforms and deployment | Linux x86_64 with SLURM is the tested target, on one site. Other clusters require site configuration, and shared filesystems that cache directories for longer than 120 seconds need a larger `latency-wait`. A stock Ubuntu container userland passed installation and the Salmon route; no container image is shipped. macOS and ARM are not qualified. |
| Intake and optional integrations | Excel supports independent bulk with optional batch and inline/external records. Setup-to-config tests mock live reference lookup; the scientific DAG is tested separately. Advanced models and reference choices need YAML. The optional upload/wizard is not an equivalent qualified interface. |

## Publication procedure

1. Preserve the tested source revision, lock hash and qualification logs; inspect the complete candidate diff.
2. Keep handbook PDFs, private audit/upgrade plans, working trackers and raw local evidence out of the public source tree. Public operational documentation must remain readable without those files.
3. Push the reviewed source changes and require the **release-validation** scientific job, including the public smoke test, to pass without skips. Configure this as a required branch check in the hosting repository.
4. Publish with the supported scope and limitations above. Do not describe this candidate as an implemented universal RNA-seq platform or as a qualified TAG-seq kit workflow.
5. Record the actual tag/revision in release notes and citation metadata; publish only the reviewed scope. Preserve repository visibility unless the owner separately requests a change.

For a runtime upgrade, solve `environment.yml` into a new prefix, export `conda list --explicit --sha256`, review the exact lock changes and repeat the qualification gates before changing the supported release environment.

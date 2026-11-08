# Release qualification

Current release: **v3.0.0** (2026-11-07). It runs on **Linux x86_64** and covers annotated organisms in four library types: **bulk RNA-seq** without UMIs, **small RNA** (animal and plant), **Lexogen QuantSeq 3′ end-tag** (FWD, REV and FWD-UMI) and **host plus pathogen** in the same bulk libraries.

| Version | Date | What it is |
|---|---|---|
| 3.0.0 | 2026-11-07 | Current release. Adds small RNA, QuantSeq 3′ end-tag and host-plus-pathogen routes, and two locked runtimes |
| 2.1.2 | 2026-10-17 | Same pipeline code as 2.1.1, with a Debian package and apt repository. Superseded by 3.0.0 |
| 2.1.1 | 2026-10-08 | Same pipeline code as 2.1.0, without the n8n integration. Superseded by 2.1.2 |
| 2.1.0 | 2026-10-08 | Added STAR, edgeR and DESeq2, read screening and the larger intake workbook. Superseded by 2.1.1 |
| 2.0.0 | 2026-09-23 | First release: Salmon with limma-voom |

A version is released only after every required check passes with none skipped.

## What was tested for 3.0.0

The SLURM runs used source `12fec99` (`git archive` SHA256 `16ec5aa910424c9b67376150cf8e7b2dc71832687292d208ebe67c6176f86535`, checked on the cluster). The hosted run used `4a56555`, which changes only the two GitHub workflow files. Later commits up to the tag change only documentation and the README animation.

| Check | Where | Result |
|---|---|---|
| Install all five locked environments from scratch and verify them | SLURM cluster (job 40546128); [hosted workflow](https://github.com/prem-p-singh/RNASeq-pipeline/actions/runs/42451991881) | Passed: core 580, STAR 37, screening 64, small RNA 89 and UMI 119 packages verified |
| 41 checks, none skipped | Same two runs | 41 passed, 0 failed, 0 skipped |
| Six-sample public study (GSE110004) through Salmon with read screening, and through STAR | Same two runs | Counts, DE tables and QC report produced on both routes |
| Reuse of the installed environments without network access | SLURM cluster (job 40546128) | Passed for all five environments |
| Launcher sending each step to SLURM as its own job | SLURM cluster (job 40546129, 1 h 36 min) | 27 rule jobs, each a separate SLURM job, all completed |
| Resume after an interrupted launch | Same job | Interrupted at 16 of 28 steps; the second launch finished with 10 worker jobs, and the finished sample was not recomputed |
| Two projects sharing one reference cache entry | Same job | One cache entry; 27 worker jobs per project; per-sample count totals differed by at most 1.28e-04 |
| Stock Ubuntu container (`condaforge/miniforge3:26.7.2-0`, Ubuntu 24.04.5) | SLURM cluster (job 40546089) | Fresh core install, `check_end_to_end.sh` and the Salmon public study passed |

The 41 checks are the 35 from 2.1 plus a unit check and a synthetic workflow check for each new route: small RNA, QuantSeq and host plus pathogen.

Cluster logs, SHA256: `tests.log` `ae170f8c0000e8bf7d5dbe1a346f9c66f4e9834e5fde3857628810b9f8c4c637`, `public.log` `e550c617e636b3b69e9eb9156e36bc5fccaffc8e0318a55b43e57589bd4be942`, `public-star.log` `9792f80a49a4287e62b1f6b3654a3f62a74aabec907dbb14460ec3e757040d3e`, `distributed.log` `b8018c16e5ad01d1ce005e8778d7047f0bdb790919c8389e63aa73023e727ab9`, `container.log` `eefc6c5cbbe159036f59b48f39dca955076e9a4fc20f41b2c7ff03cf015cc04c`. The cluster jobs used Miniforge 26.7.2 and a separate home directory, so a personal `.condarc` could not add package caches or default packages to the clean installs.

### Public-data checks for the new routes

These ran on the cluster on earlier sources of this line and were not repeated on the 3.0.0 source. Acceptance criteria were written into each fixture before the first run.

| Route | Data | Job | Result |
|---|---|---|---|
| Small RNA, animal | GSE94585 technical pools | 39055140 (commit `4ef08b0`) | 238 sequences at ≥ 100 reads per pool; median log2 ratio difference 0.218; Spearman 0.962 and 0.974 against published counts |
| Small RNA, plant | SRR3222443 to SRR3222445 | 39055140 | 15018 loci; DicerCall 21 to 24 present; 21 and 24 nt the most common lengths; a repeat run gave an identical count matrix |
| QuantSeq FWD | six public runs | 39025685 | Counts, DE tables and QC report for all six; mapping rate 0.485 (floor 0.40) |
| QuantSeq REV | six public runs | 39023623 | Counts, DE tables and QC report for all six |
| QuantSeq FWD-UMI | six public runs | 39023624 | Counts, DE tables and QC report for all six, with molecule counts and duplication rates |
| Host plus pathogen, STAR | Grapevine with Grapevine red blotch virus, PRJNA339456 | 39420845 (commit `e70f24c`) | All five acceptance criteria held |
| Host plus pathogen, STAR | Arabidopsis with *Botrytis cinerea*, PRJNA1270936 | 39505358 (commit `7d53d81`) | All five acceptance criteria held |
| Host plus pathogen, Salmon | Both studies above | 39586469, 39586470 (commit `2c25c17`) | All five acceptance criteria held on each |

For the virus, whose genes overlap and whose libraries are unstranded, per-gene counts differed between the STAR and Salmon routes. The checks do not compare the routes. See [Options](OPTIONS.md#host-and-pathogen-in-the-same-libraries) for the known limits of the host-plus-pathogen route.

## What was tested for 2.1

| Check | Where | Result |
|---|---|---|
| Install the locked environments from scratch and verify them | SLURM cluster; hosted workflow | Passed |
| 35 checks, none skipped | SLURM cluster (job 39621454); [hosted workflow on the tagged commit](https://github.com/prem-p-singh/RNASeq-pipeline/actions/runs/37890271662) | 35 passed, 0 failed, 0 skipped |
| Six-sample public study through Salmon, with read screening | Same two runs | Counts, DE tables and QC report produced |
| Same public study through STAR and featureCounts | Same two runs | Counts, DE tables and QC report produced |
| Reuse of the installed environments without network access | SLURM cluster | Passed |
| Launcher sending each step to SLURM as its own job | SLURM cluster (job 39620305), on the 2.1.0 source | Passed for the Salmon route, six samples |
| Resume after an interrupted launch | Same job | Passed; finished work was not recomputed |
| Two projects sharing one reference cache entry | Same job | Passed |
| Stock Ubuntu container | SLURM cluster (job 39618728), on the 2.1.0 source | Core install and Salmon route passed |

The rows above record v2.1.1 unless a row names another version. The last four rows ran on the 2.1.0 source and were not repeated for 2.1.1, which differs only by the removed n8n files, the test runner and documentation.

v2.1.2 changes no pipeline code; it adds the Debian package. The hosted workflow passed on its [tagged commit](https://github.com/prem-p-singh/RNASeq-pipeline/actions/runs/39074182797). The SLURM runs were not repeated for 2.1.2.

The 35 checks cover input validation, read preparation, reference and cache safety, count agreement between routes, the direction of DE results, the three DE methods, paired designs, enrichment, WGCNA, storage planning and recovery after failures.

### Test data

- **Synthetic reads:** `tests/fixtures/fastq_project.py` builds six samples and 120 genes with known fragments and planted expression changes. Real fastp, Salmon and STAR run on them from raw FASTQ to report, and the checks compare the results with the known truth.
- **Public study:** GSE110004, runs SRR6357070 to SRR6357075, three wild-type and three Rap1-AID yeast samples, taken as chromosome-I, downsampled reads from a [pinned nf-core test-data revision](https://github.com/nf-core/test-datasets/tree/626c8fab639062eade4b10747e919341cbf9b41a). The transcriptome there contains one artificial entry, `Gfp_transgene_gene`, that is absent from its GTF; the test removes it and records the change. This run shows that the routes work on real reads. It does not reproduce the study's conclusions.

## Runtime

Install with `scripts/bootstrap.sh`. Each environment is an explicit lock file with exact builds and checksums.

| Environment | Lock file | Packages | SHA256 of the lock file |
|---|---|---|---|
| Core | `environments/linux-64.explicit.txt` | 580 | `dc1d1459bf56f4a5e93c72d2a966c1a96a37168ea49e15933fe549e3549bff75` |
| STAR module | `environments/star-linux-64.explicit.txt` | 37 | `ebca83dc6a889b16291dadc923379a830b4596970b3cabe6c2c62c483c86c0d5` |
| Screening module | `environments/screen-linux-64.explicit.txt` | 64 | `065bb112579c134645324e4d9ee6a21d9aa822c3090e9ebc7475bc483852a7de` |
| Small RNA module | `environments/srna-linux-64.explicit.txt` | 89 | `bb8f4a47cb20f77d632a9a3dbd0aef77e6a1dfd6df43ee3dbb47e70e975514c3` |
| UMI module | `environments/umi-linux-64.explicit.txt` | 119 | `70ef9ba037dbfac0713cb81f345a5df7489eef5b4301ddcc4e8b420210075873` |

Main tool versions: Python 3.12.14, Snakemake 9.19.0, R 4.5.2, Salmon 1.10.3, fastp 0.23.4, MultiQC 1.35, tximport 1.38.2, limma 3.66.0, edgeR 4.8.2, DESeq2 1.50.2, variancePartition 1.40.2, clusterProfiler 4.18.4, WGCNA 1.74, STAR 2.7.11b, featureCounts 2.1.1, samtools 1.24, FastQ Screen 0.16.0, Bowtie2 2.5.4, Bowtie 1.3.1, ShortStack 4.1.0, cutadapt 5.2 and umi_tools 1.1.6.

The verifier compares the installed packages with the lock, imports the required libraries and checks where each tool is found. It does not audit every installed file. Do not install other packages into these environments.

## Repeat the tests

On a SLURM cluster, from the source checkout:

```bash
sbatch scripts/validate_slurm.sbatch
```

This installs fresh environments, runs `tests/run_all.sh --strict`, runs the public study through both routes and checks reuse without network access. Add your site's account, partition and QoS options. Set `RNASEQ_VALIDATE_ROOT` to a shared directory if worker-local `/tmp` is too small.

For the SLURM-worker, resume and shared-cache checks, and separately for a container:

```bash
RNASEQ_QUALIFY_ROOT=/shared/dir sbatch scripts/qualify_deployment.sbatch
RNASEQ_QUALIFY_ROOT=/shared/dir sbatch scripts/qualify_deployment.sbatch container
```

In an already installed environment, follow the commands under [Validation](../README.md#validation) in the README.

## Known limits

| Area | Limit |
|---|---|
| Assays | Bulk RNA-seq without UMIs, small RNA with the TruSeq Small RNA kit, the three QuantSeq kits, and host plus pathogen in bulk libraries; one library per sample. Other assay types are listed under [Future scope](../README.md#future-scope) |
| Small RNA | Reads are counted, not molecules: no UMI deduplication, isomiR resolution or phased-siRNA calling. Public checks used technical controls |
| Host plus pathogen | Two organisms only; bulk libraries only; samples with no host reads are not supported. Three defaults are assumptions, not calibrated values. See [Options](OPTIONS.md#host-and-pathogen-in-the-same-libraries) |
| Evidence | The checks use synthetic reads, one small bulk public study and the public studies named above for the new routes. They do not prove correctness for every organism, design or contrast |
| Mixed models | A random subject effect with `dream` was checked on one paired-donor fixture. Treat mixed-model results as provisional |
| Read screening | Diagnostic only. How sensitive a given panel is has not been measured |
| Enrichment | GO is tested with a local annotation database. KEGG needs live network access and organism mapping, and its live results are not covered by the tests |
| Sample exclusion | Samples failing QC can be excluded by policy. A library-type mismatch is recorded but does not stop the run; review it before interpreting results |
| References | A remote file that changes at the same URL is not detected. Use versioned URLs |
| Storage | Storage estimates are approximate and only warn. A real failed write still stops the affected step |
| Scale | Large cohorts and full-size genomes have not been benchmarked |
| Platforms | Tested on Linux x86_64 with SLURM at one site. Other clusters need their own profile settings; filesystems that cache directory listings for more than 120 seconds need a larger `latency-wait`. No container image is shipped. macOS and ARM are not supported |

## Publishing a release

1. Freeze the source and keep the lock files and test logs for it.
2. Push the source and require the `release-validation` workflow to pass with no skipped check.
3. Tag the commit that passed, and record the version in `CITATION.cff` and the changelog.
4. Describe only what was tested, with the limits above.

To change the runtime, solve `environment.yml` into a new environment, export it with `conda list --explicit --sha256`, review the differences and repeat all tests.

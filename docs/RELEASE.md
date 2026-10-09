# Release qualification

Current release: **v2.1.1** (2026-10-08). It runs on **Linux x86_64** and covers **bulk RNA-seq of annotated organisms, for libraries without UMIs**.

| Version | Date | What it is |
|---|---|---|
| 2.1.1 | 2026-10-08 | Current release. Same pipeline code as 2.1.0, without the n8n integration |
| 2.1.0 | 2026-10-08 | Added STAR, edgeR and DESeq2, read screening and the larger intake workbook. Superseded by 2.1.1 |
| 2.0.0 | 2026-09-23 | First release: Salmon with limma-voom |

A version is released only after every required check passes with none skipped.

## What was tested

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

The last four rows ran on the 2.1.0 source and were not repeated for 2.1.1, which differs only by the removed n8n files, the test runner and documentation.

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

Main tool versions: Python 3.12.14, Snakemake 9.19.0, R 4.5.2, Salmon 1.10.3, fastp 0.23.4, MultiQC 1.35, tximport 1.38.2, limma 3.66.0, edgeR 4.8.2, DESeq2 1.50.2, variancePartition 1.40.2, clusterProfiler 4.18.4, WGCNA 1.74, STAR 2.7.11b, featureCounts 2.1.1, samtools 1.24, FastQ Screen 0.16.0 and Bowtie2 2.5.4.

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
| Assays | Bulk RNA-seq without UMIs, one library per sample. Other assay types are listed under [Future scope](../README.md#future-scope) |
| Evidence | The checks use synthetic reads and one small public study. They do not prove correctness for every organism, design or contrast |
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

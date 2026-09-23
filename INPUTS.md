# Project inputs

Use the environment and launch instructions in [README.md](README.md). A project lives outside the source checkout and contains three files:

| File | Purpose |
|---|---|
| `config/config.yaml` | Organism, versioned references, assay/layout, model, contrasts and requested analyses |
| `config/samples.tsv` | One row per biological measurement, with its read file or pair and model covariates |
| `config/thresholds.yaml` | Explicit QC and analysis thresholds; project overrides are merged with defaults |

Copy the templates from `config/` into your project, then edit them. Do not run the shipped grapevine/model examples unchanged on a different experiment.

## Sample sheet

A paired-end example is:

```tsv
sample_id	fastq_url	fastq_url_r2	treatment
control_1	/shared/reads/control_1_R1.fastq.gz	/shared/reads/control_1_R2.fastq.gz	control
control_2	/shared/reads/control_2_R1.fastq.gz	/shared/reads/control_2_R2.fastq.gz	control
control_3	/shared/reads/control_3_R1.fastq.gz	/shared/reads/control_3_R2.fastq.gz	control
treated_1	/shared/reads/treated_1_R1.fastq.gz	/shared/reads/treated_1_R2.fastq.gz	treated
treated_2	/shared/reads/treated_2_R1.fastq.gz	/shared/reads/treated_2_R2.fastq.gz	treated
treated_3	/shared/reads/treated_3_R1.fastq.gz	/shared/reads/treated_3_R2.fastq.gz	treated
```

Use unique sample IDs, both mates for every paired library, and columns for every model variable. For single-end input, omit or leave the second-mate column blank and choose `rnaseq_single`. Repeated measurements require their biological-unit column and an explicit random-effect term; repeated rows from one unit are not independent biological replicates.

Local input files must be readable on every worker. HTTPS inputs are downloaded per sample. S3 requires additional AWS tooling and credentials and is outside the locked release's tested input routes. Local source FASTQs are read in place and never deleted by workflow cleanup.

Legacy sample sheets accept one read file/pair per row. For multiple sequencing runs or lanes, set `samples.metadata_dir: metadata` and `samples.sheet: metadata/samples.tsv`, or select workbook tables in the Excel intake. Canonical bulk execution accepts one non-UMI library per sample, with multiple run/lane units. Libraries declare `assay_family=bulk`, layout, strandedness and `umi=no`; each paired unit needs R1 and R2. Strandedness is forward, reverse, unstranded or unknown. Relative canonical TSV URIs resolve beside the metadata tables; relative workbook URIs resolve beside the workbook.

If any `metadata/samples.tsv`, `metadata/libraries.tsv` or `metadata/reads.tsv` file is present, preflight requires all three populated tables. Required cells and parent relationships must be complete; every sample needs a library and every library needs reads. Reused lane numbers in different sequencing runs need distinct `run` values. Canonical execution checks gzip FASTQ structure, mate IDs/counts, declared checksums/sizes and repeated source assignments before merging within the library. Each library retains `read_preparation.json` with source hashes and per-unit fragment counts. Original FASTQs are preserved; successful cleanup removes only owned merged/trimmed copies. Multiple prepared libraries per sample and UMI/barcode processing remain unavailable.

Sample identifiers are text, including leading zeros and a literal `NA`. When preparing metadata in Excel, format the identifier column as **Text before entering IDs**: zeros already removed by Excel cannot be recovered from a numeric value or its display format.

## References and analysis

Provide a transcriptome FASTA and a matching GTF from a recorded reference release. Every quantifiable transcript must map unambiguously to one annotated gene. The pipeline stops on unmatched or duplicated transcript identifiers. Record accession and organism identity; use immutable source URLs. The current default is a transcriptome-only Salmon index. Do not describe it as decoy-aware unless a compatible decoy reference and identifier list were explicitly supplied.

Set the fixed-effects formula, any random-effects term, primary factor and contrasts to match the experiment. Preflight checks estimability and supported capabilities. `analysis.backend: auto` follows the dependence structure; unavailable backend choices block. Optional GO/KEGG/WGCNA settings must match the requested objectives and available annotation.

For GO, provide a compatible OrgDb package and the correct `orgdb.key_type` (`ENTREZID` for suitable public packages, or `GID` for a compatible custom package). Non-Entrez input gene identifiers require a valid annotation mapping. A numeric-looking identifier alone is not biological proof of the correct namespace; review mapping and reference provenance.

## Optional setup helpers

`python3 scripts/setup.py --project-dir PROJECT INPUTS.yaml` can generate the project files from a spreadsheet and read-source directory/list. Start from `scripts/setup_inputs.template.yaml` and review every generated reference URL, sample match and model field. Automatic lookups can change as external services change; preserve the resolved configuration.

The interactive `scripts/new_project.sh` helper uses the same locked runtime and reads an inbox at `~/new_project_inbox/PROJECT_NAME`. `scripts/start_new.sh` is a convenience uploader/wizard launcher whose transfer and interactive path is not covered by the release's end-to-end qualification. Explicit TSV/YAML input is the tested release interface.

Storage is planned from the dataset and filesystem capacity. New setup files have `storage_budget_gb: null`; there is no 20 GB platform cap. A supplied legacy value is advisory quota information. Low or unknown capacity warns without blocking launch or reducing concurrency; real failed writes remain task failures.

## Direct Excel intake: initial bulk support

The existing `intake_template.xlsx` can now drive setup directly. Fill the **Bulk RNA-seq** tab. Choose external files for the existing metadata-plus-folder workflow, or workbook tables to fill **Samples**, **Libraries** and **Reads** starting at row 5. In table mode, clear the external metadata and FASTQ folder answers. Add sample covariate columns as needed, matching the declared model fields. Then run in the prepared analysis environment:

```bash
python3 scripts/setup.py --intake /path/completed_intake.xlsx --project-dir /path/new_project
bash submit.sh -d /path/new_project
```

Use `--intake-sheet 'Bulk RNA-seq'` to select a sheet explicitly. Relative paths resolve beside the workbook. Setup saves the original workbook, checksum, field-to-cell map and resolved setup inputs inside the project's `intake/` directory. It will not overwrite an existing configuration; resume existing projects with `submit.sh`.

This first adapter supports independent bulk designs with an optional declared batch. It preserves enrichment/WGCNA choices and strand expectations, checks the declared smallest group size, and rejects formulas and ambiguous answers. Known organism presets accept common/scientific names; otherwise enter an NCBI taxID. Platform, preparation, tissue and contact email are recorded context only and are explicitly reported as such. They do not yet configure trimming or scheduler mail.

Schema 3 adds source ownership and the canonical tables to schema 2's stable field IDs, explicit raw-FASTQ stage and non-UMI eligibility. Legacy schema 1/2 questionnaires remain readable. Unknown UMI status must be resolved; UMI input and processed-object/instrument input are not supported. Kit and strandedness defaults do not invent a protocol or orientation. Repeated-measures fields and non-bulk routes remain pending. TAG-seq and small-RNA intake are rejected explicitly; use the provisional YAML path for explicit repeated-measures models. See the [intake coverage audit](docs/INTAKE_TEMPLATE_AUDIT.md) and [v2 scope](docs/RELEASE.md).

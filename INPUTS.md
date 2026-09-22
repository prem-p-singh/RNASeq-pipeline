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

The executable workflow currently accepts one read file/pair per sample-sheet row. Multiple lanes need an external, documented merge. The metadata templates for samples/libraries/reads support validation work for the broader roadmap; they are not yet an alternative executable intake contract.

## References and analysis

Provide a transcriptome FASTA and a matching GTF from a recorded reference release. Every quantifiable transcript must map unambiguously to one annotated gene. The pipeline stops on unmatched or duplicated transcript identifiers. Record accession and organism identity; use immutable source URLs. The current default is a transcriptome-only Salmon index. Do not describe it as decoy-aware unless a compatible decoy reference and identifier list were explicitly supplied.

Set the fixed-effects formula, any random-effects term, primary factor and contrasts to match the experiment. Preflight checks estimability and supported capabilities. `analysis.backend: auto` follows the dependence structure; unavailable backend choices block. Optional GO/KEGG/WGCNA settings must match the requested objectives and available annotation.

For GO, provide a compatible OrgDb package and the correct `orgdb.key_type` (`ENTREZID` for suitable public packages, or `GID` for a compatible custom package). Non-Entrez input gene identifiers require a valid annotation mapping. A numeric-looking identifier alone is not biological proof of the correct namespace; review mapping and reference provenance.

## Optional setup helpers

`python3 scripts/setup.py --project-dir PROJECT INPUTS.yaml` can generate the project files from a spreadsheet and read-source directory/list. Start from `scripts/setup_inputs.template.yaml` and review every generated reference URL, sample match and model field. Automatic lookups can change as external services change; preserve the resolved configuration.

The interactive `scripts/new_project.sh` helper uses the same locked runtime and reads an inbox at `~/new_project_inbox/PROJECT_NAME`. `scripts/start_new.sh` is a convenience uploader/wizard launcher whose transfer and interactive path is not covered by the release's end-to-end qualification. Explicit TSV/YAML input is the tested release interface.

Storage is planned from the dataset and filesystem capacity. New setup files have `storage_budget_gb: null`; there is no 20 GB platform cap. A manually supplied legacy value represents an explicit site constraint and produces a migration notice.

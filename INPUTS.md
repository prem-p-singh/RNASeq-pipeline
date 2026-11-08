# Project inputs

Install and launch as described in the [README](README.md). This page explains what you provide.

There are two ways to describe a study:

- **The Excel workbook** ([intake_template.xlsx](intake_template.xlsx)). Recommended.
- **Text files**: a YAML configuration and tab-separated tables. Use these for scripted setups or anything the workbook does not cover.

Both produce the same project folder, which lives outside the source checkout:

| File | Purpose |
|---|---|
| `config/config.yaml` | Organism, references, read layout, model, contrasts and requested analyses |
| Sample records | `config/samples.tsv`, or the three tables under `metadata/` |
| `config/thresholds.yaml` | QC and analysis thresholds; your values override the defaults |

## Excel workbook

Fill the tab for your assay (**Bulk RNA-seq**, **sRNA-seq** or **TAGseq**) with plain answers, not spreadsheet formulas. Then either let the launcher set up and run the project, or set it up first and review it:

```bash
bash submit.sh --intake /path/study.xlsx -d /path/new_project --executor slurm
# or
python3 scripts/setup.py --intake /path/study.xlsx --project-dir /path/new_project
bash submit.sh -d /path/new_project
```

Setup keeps a copy of the workbook, its checksum and a record of which cell each setting came from, under `intake/` in the project. It will not overwrite an existing configuration. If you change the workbook later, the launcher lists the changed settings and the stages they affect, then stops; add `--accept-revision` to apply the change. A changed project name needs a new project folder.

### Where the sample records come from

| Choice | What to fill |
|---|---|
| **Workbook tables** | The Samples, Libraries and Reads tabs, from row 5. Leave the external metadata file and FASTQ folder answers empty |
| **External files** | A metadata file (CSV, TSV or XLSX) and a FASTQ folder, one file or pair per sample. Leave the three tabs empty |

Use workbook tables when a sample was sequenced on more than one run or lane. Relative paths are read relative to the workbook.

### Questions that shape the analysis

- **Analysis goal:** differential expression, expression only (counts, optional WGCNA) or QC only. Enrichment needs the differential expression goal.
- **Quantifier:** `auto`, `salmon` or `star`. STAR needs a declared strand.
- **DE method:** `auto`, `limma_voom`, `edger_ql` or `deseq2`. With `auto`, a new project gets DESeq2 for 3 to 12 biological units in the smallest group and limma-voom above 12. Two units per group needs an explicit method. The method never changes the model.
- **Paired or repeated measures:** answer yes and name the subject column. `fixed block` adds the subject as a fixed term, for paired designs with any DE method. `random` fits a random subject effect with dream. Repeated rows of one subject count once toward the smallest group.
- **Fixed-effects formula:** optional, for example `~ batch + genotype * treatment + age`. It must include the primary factor, a declared batch and a fixed subject block.
- **Contrasts tab:** optional, one row per comparison. A pairwise row gives `factor`, optional `by`, and `numerator` and `denominator` levels; blank levels compare all pairs. A linear row gives `weights` as `coefficient=number` pairs separated by `;`. An empty tab keeps the automatic contrasts.
- **Library protocol:** a named protocol (TruSeq, NEBNext) sets adapters and strand from [config/kit_profiles.yaml](config/kit_profiles.yaml). An answer that conflicts with it is an error. You can instead give R1 and R2 adapters and a minimum read length yourself.
- **Platform:** sets poly-G trimming (on for NovaSeq and NextSeq, off for MiSeq) unless you choose it explicitly.
- **Screening:** yes builds a contamination panel from the project reference and PhiX; add more references as `name=FASTA` pairs.
- **Organism:** a common or scientific name for known organisms, otherwise an NCBI taxID.
- **Custom references:** give a genome and GTF for STAR; a transcriptome and GTF for Salmon, plus the genome for decoys. A partial set is rejected, so two annotation releases cannot be mixed by accident. Leave all of them blank for automatic selection.

Preparation, tissue and contact email are recorded for reference and change nothing. The Bulk RNA-seq, sRNA-seq and TAGseq tabs can be run in this release. Libraries with UMIs are supported only for the QuantSeq FWD-UMI kit. Workbooks made for earlier releases still load.

## Text files

Copy the templates from `config/` into your project and edit them. They are examples: replace the organism, references, samples, model and contrasts with your own.

### Sample sheet

A paired-end example:

```tsv
sample_id	fastq_url	fastq_url_r2	treatment
control_1	/shared/reads/control_1_R1.fastq.gz	/shared/reads/control_1_R2.fastq.gz	control
control_2	/shared/reads/control_2_R1.fastq.gz	/shared/reads/control_2_R2.fastq.gz	control
control_3	/shared/reads/control_3_R1.fastq.gz	/shared/reads/control_3_R2.fastq.gz	control
treated_1	/shared/reads/treated_1_R1.fastq.gz	/shared/reads/treated_1_R2.fastq.gz	treated
treated_2	/shared/reads/treated_2_R1.fastq.gz	/shared/reads/treated_2_R2.fastq.gz	treated
treated_3	/shared/reads/treated_3_R1.fastq.gz	/shared/reads/treated_3_R2.fastq.gz	treated
```

Use unique sample IDs, both mates for every paired library, and a column for every model variable. For single-end input, leave the second-mate column blank and choose `rnaseq_single`. Repeated measurements of one subject are not independent replicates: give the subject column and model it as a fixed block or a random effect.

A sample sheet takes one file or pair per row. For samples sequenced on several runs or lanes, use the three tables instead by setting `samples.metadata_dir: metadata` and `samples.sheet: metadata/samples.tsv`.

### Samples, Libraries and Reads tables

These are `metadata/samples.tsv`, `metadata/libraries.tsv` and `metadata/reads.tsv`; templates are in [templates/](templates/).

- If any of the three is present, all three must be filled in. Every sample needs a library and every library needs reads.
- One library per sample, with `assay_family=bulk`, the layout, the strandedness (forward, reverse, unstranded or unknown) and `umi=no`.
- Each paired read unit needs R1 and R2. Lane numbers reused in different sequencing runs need distinct `run` values.
- Before merging lanes, the workflow checks gzip FASTQ structure, mate IDs and counts, and any checksums or sizes you declared. Each library keeps a `read_preparation.json` with source hashes and fragment counts.
- Relative read paths are read relative to the tables.

### Read files

Local files must be readable on every worker and are read in place; the workflow never deletes them. HTTPS inputs are downloaded per sample. S3 needs extra AWS tooling and credentials and is not a tested input.

### Sample IDs

Sample IDs are text, including leading zeros and a literal `NA`. In Excel, format the ID column as **Text before typing the IDs**: zeros that Excel has already removed cannot be recovered.

## References

Provide a transcriptome FASTA and a matching GTF from one reference release; STAR also needs the genome. Every transcript must map to exactly one annotated gene, and the workflow stops on unmatched or duplicated transcript IDs. Use versioned URLs, because a file that changes at the same URL is not detected.

New projects build the Salmon index with genome decoys. Set `reference.decoys: null` for a transcriptome-only index.

## Model, contrasts and optional analyses

Set the fixed-effects formula, any random-effects term, the primary factor and the contrasts to match the experiment. The pre-run check confirms the design can be estimated and blocks unsupported choices. In an existing configuration `analysis.backend: auto` keeps limma-voom, and a random-effects term selects dream. Details are in [docs/OPTIONS.md](docs/OPTIONS.md).

For GO, provide a compatible OrgDb package and the right `orgdb.key_type` (`ENTREZID` for suitable public packages, `GID` for a compatible custom one). Gene IDs that are not Entrez IDs need a valid mapping. An ID that merely looks numeric is not proof of the right namespace, so check the mapping.

## Storage

Storage is estimated from the dataset and the filesystem. There is no fixed size limit. Low or unknown capacity gives a warning and does not block the launch or reduce parallel jobs; a real failed write still stops the affected step.

## Small RNA (workbook sRNA-seq tab)

- **small_rna_profile:** `animal` (Bowtie 1 to the species' mature miRNAs from miRBase) or `plant` (ShortStack loci on the genome). Both passed their public-data checks with the TruSeq Small RNA kit (FARM job 39055140) and are validated for that kit.
- **adapter_kit:** a named non-UMI small-RNA kit from [config/small_rna_kits.yaml](config/small_rna_kits.yaml); UMI small-RNA kits are rejected, since UMI-aware counting is not implemented.
- **mirbase_prefix:** required for `animal`; optional for `plant`, where it only annotates ShortStack loci with known miRNAs and does not change which loci are called.
- **reference_genome:** required for `plant`; not used by `animal`.
- **Sample and read records:** external files only. Workbook tables are refused on this tab.
- Counts are **reads, not molecules**: no UMI deduplication, no isomiR resolution, and the plant route does not call phased siRNA clusters.

## 3′ end-tag libraries (workbook TAGseq tab)

- **kit:** `QuantSeq FWD`, `QuantSeq REV` or `QuantSeq FWD-UMI` from [config/tagseq_kits.yaml](config/tagseq_kits.yaml). In a project configuration the same name goes in `analysis.kit` with `samples.seq_type: tagseq`. Any other name is refused.
- **reference_genome** and **reference_gtf:** both required; the named-kit route aligns to the genome with STAR and counts with featureCounts.
- **Strand comes from the kit:** forward for FWD and FWD-UMI, reverse for REV. An `expected_libtype` that disagrees is an error; blank defers to the kit.
- **Sample and read records:** external files only. Workbook tables are refused on this tab.
- **Trimming comes from the kit:** a preprocessing setting that contradicts the kit (for example another `trim_front`) is an error. Reads are single-end Read 1 only.
- **UMI (FWD-UMI only):** fastp (`--umi --umi_len 6 --umi_skip 4`) moves the 6 nt UMI into the read name and drops the TATA spacer, then applies the 12 nt front trim. `star/gene_counts.tsv` counts reads; `umi/gene_counts.tsv` and the cohort `counts.tsv` count molecules.
- All three kits passed their public-data checks (FARM jobs 39025685, 39023623, 39023624) and are listed in `tagseq.validated_kits` in `config/spec.yaml`. A `tagseq` project without a kit stays provisional.
- A `tagseq` project without a kit keeps the earlier provisional route unchanged.

## Host and pathogen in the same libraries

Status: `validated` on both the STAR and the Salmon route (public checks on FARM: STAR jobs 39420845 and 39505358, Salmon jobs 39586469 and 39586470). Bulk RNA-seq only (single-end or paired-end); not for TAG-seq or small RNA.

In the workbook (schema 8), the Bulk RNA-seq tab has these answers:

- **dual_organism:** `yes` assigns each read to the host or the pathogen and analyses each organism on its own. The host is the organism of the project; give its genome and GTF as the custom reference, plus its transcripts for Salmon.
- **pathogen_label:** a short name of letters, digits, `_`, `.` and `-`, different from the host label.
- **pathogen_genome** and **pathogen_gtf:** both required. The GTF needs `exon` rows with `gene_id` and `transcript_id`.
- **pathogen_transcriptome:** required unless the quantification method is `star`.
- **pathogen_objective:** `counts only` (default) or `differential expression`. Pathogen enrichment and co-expression are not available.
- **pathogen_min_reads:** fewest pathogen reads for a sample to be analysed for the pathogen (default 500).

In a project configuration:

```yaml
dual_organism:
  enabled: true
  multimap_limit: 20              # assumption
  ambiguous_warn_fraction: 0.05   # assumption
  host:
    label: grapevine
    genome_fasta_url: ...
    gtf_url: ...                  # required
    transcriptome_fasta_url: ...  # required on the Salmon route
    objectives: [gene_expression, differential_expression]  # default: analysis.objectives
  pathogen:
    label: GRBV
    genome_fasta_url: ...
    gtf_url: ...                  # required
    transcriptome_fasta_url: ...  # required on the Salmon route
    objectives: [gene_expression] # default; or add differential_expression
    min_assigned_reads: 500       # assumption
```

- `host.objectives` defaults to the project objectives (`analysis.objectives`); `pathogen.objectives` defaults to counts only (`[gene_expression]`). A project with only the `qc` objective is refused.
- The quantifier, read layout, preprocessing, model and contrasts at the top level apply to both organisms. A STAR project must declare `samples.expected_libtype`.
- `analysis.salmon_position_bias` (default `true`) keeps Salmon position-bias correction on; set it to `false` for a very small transcript set, because Salmon 1.10.3 can crash in that correction when a transcript length class has no reads. The pathogen project always runs with it off and the host project follows this setting.
- Each read gets one label: `host`, `pathogen`, `ambiguous` (best alignments in both organisms) or `unassigned` (no alignment; this includes a pair with one mate in each organism, which STAR does not align). Ambiguous and unassigned reads are counted in `results/assignment.tsv` and enter neither count table.
- Outputs: `results/assignment.tsv`, `results/host/` and `results/pathogen/` (each laid out like a single-organism project), `results/pathogen/samples_not_quantified.tsv` and `results/pathogen/status.json`.
- Samples below `min_assigned_reads` are not quantified for the pathogen and are listed in `samples_not_quantified.tsv`. Pathogen differential expression is reported as not estimable when a level of the primary factor keeps fewer than 2 quantified samples.
- `dual_organism.enabled: false` (the default) leaves every other route unchanged.
- Out of scope: more than two organisms; dual-organism TAG-seq or small RNA; virus-derived small RNA; cross-organism network analysis; pathogen strain or variant analysis. See [docs/OPTIONS.md](docs/OPTIONS.md#host-and-pathogen-in-the-same-libraries) for the known limits.

## Other setup helpers

`python3 scripts/setup.py --project-dir PROJECT INPUTS.yaml` builds the project files from a YAML description; start from `scripts/setup_inputs.template.yaml` and review every generated reference URL, sample match and model field. `scripts/new_project.sh` and `scripts/start_new.sh` are interactive conveniences and are not covered by the release tests.

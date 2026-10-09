# Options in detail

This page explains the choices listed under [Options](../README.md#options) in the README.

## Raw read QC without analysis

An observation-only QC command is available for bulk projects that use the Samples, Libraries and Reads tables. Activate the locked Linux environment, then run:

```bash
python3 "$REPO/scripts/raw_qc.py" -d "$PROJECT" \
  --out "$PROJECT/raw_qc_run1" --threads 2
```

This command checks lane/mate integrity and writes per-sample fastp JSON/HTML, read provenance, `manifest.json` and a cohort `index.html`. It needs no reference, Salmon index, annotation or DE model. It disables trimming/filtering and preserves original reads; temporary merged reads are removed after inspection. Use a new output directory for each run. A failure stops the command with a failed status while previously completed reports remain available.


## Analysis goals and preprocessing

Preprocessing runs separately from quantification and keeps fastp reports plus a read-loss ledger. `analysis.objectives: [qc]` builds the preprocessing report without a reference or DE model; `[gene_expression]` also produces counts without running inference. The standalone command above remains observation-only; the DAG QC objective performs the configured preprocessing. Named kit profiles in [config/kit_profiles.yaml](../config/kit_profiles.yaml) set TruSeq/NEBNext adapters and strand; they are not UMI, end-tag or small-RNA profiles. Fastp overrepresentation is a diagnostic, not a contamination verdict.

## Workbook

The workbook exposes the analysis goal (differential expression, expression only or QC only), Salmon/STAR selection, DE method, subject column for repeated measures (random subject with dream, or a fixed subject block for pairing), an optional fixed-effects formula (interactions, continuous covariates), an optional Contrasts tab, complete custom reference sets, genomic-decoy choice, a named library protocol, platform-based poly-G handling, adapters, minimum retained read length and contamination screening. Custom reference paths resolve beside the workbook; supplying a custom set bypasses automatic reference selection. Workbooks from earlier releases remain readable. Only the Bulk RNA-seq tab can be run.

The main launcher accepts the bulk workbook directly:

```bash
bash "$REPO/submit.sh" --intake /path/study.xlsx -d "$PROJECT" --plan-only
bash "$REPO/submit.sh" --intake /path/study.xlsx -d "$PROJECT" --executor local --cores 4
# Or use --executor slurm with the configured site profile.
```

Planning needs Python with pandas/openpyxl/PyYAML and R with jsonlite/emmeans; it does not build the scientific runtime or download annotation tables. It may query reference metadata services. A normal launch prepares the locked runtime automatically. Resume accepts the current workbook checksum. A changed workbook prints which settings changed and which stages they reach, then stops; add `--accept-revision` to apply it (the previous configuration is kept under `intake/<old checksum>/`, and Snakemake reruns only affected work). A changed project name still needs a new project directory.

## Differential expression

For explicit fixed-effect inference, the configuration accepts `analysis.backend: limma_voom`, `edger_ql`, or `deseq2`; `auto` in an existing configuration retains limma-voom. New projects created by setup record the handbook default instead: DESeq2 for 3 to 12 biological units in the smallest group, limma-voom above 12. Two units per group is exploratory and requires an explicit method choice. Random-effects models select dream. These methods share count provenance and declared contrasts. DESeq2 uses rounded lengthScaledTPM counts from Salmon or raw integer gene counts from STAR, its own size-factor estimation and Wald statistics; edgeR uses TMM and quasi-likelihood F tests. Neither adds another transcript-length offset.

## Contrasts and Salmon decoys

Contrasts may select `numerator` and `denominator` levels instead of `reverse`, or use `type: linear` with named design-coefficient `weights`. Preflight validates these against the actual model. New setup projects default to genomic decoys; existing configurations keep their setting. To enable genome decoys on an existing project, set `reference.decoys: genome` with a compatible `genome_fasta_url`. The builder combines transcriptome and genome, checks identifiers, records hashes and verifies cache reuse. Its default worker request is 64 GB RAM; override rule resources for your reference where appropriate. `decoys: null` explicitly retains the legacy transcriptome-only index. See [Salmon's construction method](https://salmon.readthedocs.io/en/latest/salmon.html).

## STAR and featureCounts

For genomic alignments and exon-level gene counting, set `analysis.quantifier: star` (or request `analysis.objectives: [alignment]` with `quantifier: auto`), supply compatible `reference.genome_fasta_url` and `gtf_url`, and declare strandedness. Legacy paired libraries use `IU`, `ISF` or `ISR`; single-end libraries use `U`, `SF` or `SR`. Canonical libraries use their own declared strand. Unknown strandedness blocks this route. Sample number does not choose the quantifier.

STAR 2.7.11b and featureCounts 2.1.1 produce sorted/indexed BAMs and raw exon-union gene counts. Only uniquely aligned, unambiguously assigned reads/fragments count; paired mates count as one fragment. Counts have **no transcript-length scaling**. Gene lengths are exon-union lengths. The launcher installs and verifies a separate [locked tool module](../environments/star-linux-64.explicit.txt), retaining the core R runtime. Index and alignment workers request 64 GB and 40 GB respectively; these are worker memory requests, not storage limits. Missing declared BAM indexes and gene-length files trigger recovery.

## Read screening

Optional post-preprocessing screening uses FastQ Screen 0.16.0 with Bowtie2 2.5.4. Declare the expected organism and any suspected contaminants as prebuilt Bowtie2 indexes (`index`) or as FASTA sources (`fasta`: a path, URL, `transcriptome` or `genome`) that the workflow indexes under `screening_panel/`:

```yaml
screening:
  enabled: true
  fragments: 100000
  seed: 1
  references:
    - {name: Host, role: expected, index: /references/host_transcriptome}
    - {name: PhiX, role: possible_contaminant, fasta: https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_001422.1&rettype=fasta&retmode=text}
```

Relative paths resolve against the project directory. FASTA members are downloaded or copied, hashed and indexed with the screening module's bowtie2-build; the build version and FASTA checksum are recorded in `screening_panel/<name>/built.txt`. The workbook's screening answer builds the panel automatically: the project transcriptome (genome for STAR) as expected, PhiX, and any extra `name=FASTA` entries as possible contaminants. Include an intended pathogen as `expected`, never as an assumed contaminant.

The launcher prepares the separate [screening runtime](../environments/screen-linux-64.explicit.txt) only when enabled. Uniform seeded sampling retains paired fragments together, then screens mates separately: reported counts are **reads, not fragments**. Per-sample `screening/screening.json`, `.tsv` and `.html` record sampling, reference-file checksums, exclusive/shared matches and zero reads removed. Disabled runs explicitly say `not_performed`. Screening neither excludes samples nor removes reads. Shared matches are ambiguous; matches are not organism abundance or proof of contamination, and absence of a match cannot rule out organisms omitted from the panel. Bowtie2 is not splice-aware: prefer a compatible host transcriptome for RNA screening and interpret genomic screens accordingly. Panel sensitivity and study-specific thresholds still need qualification.

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

## Small RNA

Two small-RNA routes read `samples.seq_type: srnaseq` with `analysis.small_rna_profile: animal` or `plant`, plus a named kit from [config/small_rna_kits.yaml](../config/small_rna_kits.yaml) (currently TruSeq Small RNA; UMI small-RNA kits are rejected, not silently counted without deduplication). The animal route aligns trimmed reads with Bowtie 1 (`-v 0 -a --best --strata`) against the species' mature miRNAs from miRBase, keyed by `organism.mirbase_prefix`, and counts uniquely aligned reads per miRNA; a read aligning to more than one miRNA is reported as multi-mapped and excluded, never split. The plant route runs ShortStack 4.1.0 once for the whole cohort against `reference.genome_fasta_url` (optionally with known miRNAs for annotation), producing one joint locus/count matrix with each locus's `DicerCall` size class and `MIRNA` call and overlapping known miRNAs (`known_miRNAs`). Both routes report **reads, not molecules**: there is no UMI deduplication, no isomiR resolution, and the plant route does not call phased siRNA clusters. The workbook's `sRNA-seq` tab launches either profile directly. The TruSeq Small RNA kit is validated for both profiles by their public-data checks (`srnaseq.validated_kits` in `config/spec.yaml`). See [tests/check_small_rna_public.sh](../tests/check_small_rna_public.sh) for the animal/plant public-data qualification checks (technical controls only; run separately from `tests/run_all.sh` because they need downloaded public data).

## 3′ end-tag libraries (QuantSeq)

The 3′ end-tag route covers three Lexogen QuantSeq kits named in [config/tagseq_kits.yaml](../config/tagseq_kits.yaml): `QuantSeq FWD`, `QuantSeq REV` and `QuantSeq FWD-UMI`. Choose the kit in the workbook's `TAGseq` tab or set `analysis.kit` with `samples.seq_type: tagseq`, a genome (`reference.genome_fasta_url`) and a GTF. The kit sets trimming (12 nt front trim for FWD and FWD-UMI, poly(A) and the Read 1 adapter) and the strand: STAR/featureCounts count sense reads only (`-s 1` for FWD and FWD-UMI, `-s 2` for REV), and a declared `expected_libtype` that disagrees with the kit is an error. For FWD-UMI, fastp moves the 6 nt UMI into the read name and skips the TATA spacer; preprocessing first checks that at least 80% of reads carry TATA after the UMI. `results/quant/<sample>/star/gene_counts.tsv` then holds **reads** and `results/quant/<sample>/umi/gene_counts.tsv` holds **molecules** (umi_tools per-gene deduplication); the cohort `counts.tsv` uses molecules and `counts_provenance.json` records `count_unit: molecules`. No count is length-corrected. All three kits are validated by their public-data checks (`tagseq.validated_kits` in `config/spec.yaml`). The public-data checks are `tests/check_tagseq_public.sh --fwd|--rev|--umi`, on data from `tests/fixtures/prepare_tagseq_public.py`. `tests/check_tagseq_workflow.sh` covers strand and UMI counting on synthetic reads. A `tagseq` project without `analysis.kit` keeps the earlier provisional Salmon behaviour unchanged.

## Host and pathogen in the same libraries

A dual-organism route handles libraries that hold RNA from a host and a pathogen, for example grapevine with a virus or a plant with a fungus. Set `dual_organism.enabled: true` with a `host` and a `pathogen` block, or answer yes to "Host and pathogen in the same libraries" in the workbook's Bulk RNA-seq tab; [INPUTS.md](../INPUTS.md#host-and-pathogen-in-the-same-libraries) lists the settings. It works with single-end and paired-end bulk RNA-seq on the Salmon or the STAR route, not with TAG-seq or small RNA. The route is `validated`: both public studies passed on the STAR route and on the Salmon route (job IDs in [RELEASE.md](RELEASE.md)).

Reads are preprocessed once. STAR then aligns each read or pair to one combined genome in which every sequence name starts with its organism label, and the best alignments decide one of four labels:

| Label | Best alignments are in | What happens to the read |
|---|---|---|
| `host` | the host only | goes to the host analysis |
| `pathogen` | the pathogen only | goes to the pathogen analysis |
| `ambiguous` | both organisms | counted, in neither count table |
| `unassigned` | no alignment, more best hits than `multimap_limit`, or a pair with one mate in each organism (STAR does not align such a pair) | counted, in neither count table |

`results/assignment.tsv` has one row per sample: `reads`, the four label counts, `too_many_loci` (the unassigned reads that had more best hits than the limit), the four fractions, and a `warnings` column that is filled when the ambiguous fraction is above `ambiguous_warn_fraction`. The four label counts must add up to `reads`; a sample where they do not fails. The assignment table is this file; the HTML reports do not include it.

The launcher then runs the existing route once per organism, on that organism's reads against that organism's own reference, as ordinary projects under `dual/host` and `dual/pathogen`. `results/host/` and `results/pathogen/` are each laid out like the results of a single-organism project, with their own counts, normalization, QC and report. The sample sheet, model and contrasts are declared once and used for both. `submit.sh` stops at the first project that fails, and arguments after `--` go to the assignment run and to both organism runs. `--dry-run` shows the parent run only. The split read files stay under `<output>/quant/<sample>/dual/` after the per-organism runs, because a relaunch needs them; they add about one more copy of the trimmed reads.

A sample with fewer than `min_assigned_reads` pathogen reads is not quantified for the pathogen. It is listed in `results/pathogen/samples_not_quantified.tsv` with its pathogen read count and the reason, and it stays in `results/assignment.tsv` and in the host analysis. `results/pathogen/status.json` records which samples were quantified and the state of pathogen differential expression: `not_requested`, `requested`, `not_estimable` (a level of the primary factor has fewer than 2 quantified samples, which is the usual case when the controls carry no pathogen reads) or `skipped` (no sample reached the floor, so no pathogen project is run). The pathogen analysis is counts only by default; differential expression can be requested. Pathogen enrichment and co-expression are not available. When an older `results/pathogen` folder remains from an earlier run, `status.json` is the authority on what the current run produced.

Three defaults are assumptions, not calibrated values: `multimap_limit: 20`, `ambiguous_warn_fraction: 0.05` and `pathogen.min_assigned_reads: 500`. The last was lowered from 10000 after a public grapevine red blotch virus run, where infected berry libraries of 11 to 20 million reads held 906 to 2415 virus reads and healthy ones held 1 to 4.

Known limits of the dual-organism route:

- Out of scope: more than two organisms; dual-organism TAG-seq or small RNA; virus-derived small RNA; cross-organism network analysis; pathogen strain or variant analysis.
- Salmon 1.10.3 with position-bias correction (`--posBias`) can crash with a segmentation fault when a transcript length class has no reads. The crash is intermittent. Measured on a synthetic set of one expressed 5000 nt transcript and two unexpressed 2000 nt transcripts: 45 of 60 runs crashed at 300 reads, 37 of 60 at 3000 reads and 26 of 60 at 20000 reads. The same set with `--gcBias --seqBias` only crashed in 0 of 40 runs, and sets of 2, 5, 20 or 60 expressed transcripts of 2000 nt crashed in 0 of 40 runs with all three corrections. The pathogen project therefore always runs Salmon without position-bias correction. The host project keeps it unless `analysis.salmon_position_bias: false` is set; set that for any project with a very small transcript set. The STAR route is not affected.
- STAR 2.7.11b index builds with an annotation returned 0 for random genomes of 3206 bp, 10 kb and 30 kb. An earlier 10 kb fixture aborted, so that abort is not explained by size alone. On the STAR route the pathogen is indexed alone; the 3206 bp Grapevine red blotch virus genome was indexed and quantified this way in the public check.
- Both organisms need a GTF with `exon` rows that carry `gene_id` and `transcript_id`. The NCBI RefSeq GTF of Grapevine red blotch virus (NC_022002.1) has gene, CDS, start_codon and stop_codon rows and no exon rows, so it needs exon rows added before use; `tests/fixtures/prepare_dual_public.py` shows one way.
- A sample with no host-assigned reads (for example a pathogen-only control) is not supported.
- On the STAR route a read that overlaps two annotated genes is counted for neither (featureCounts is run without `-O`). Compact virus genomes have overlapping genes.

`tests/check_dual_workflow.sh` checks exact labels and per-organism counts on synthetic reads with both routes. `tests/check_dual_public.sh --virus|--cellular` runs pinned public runs (`tests/fixtures/dual_public.json`, prepared by `prepare_dual_public.py`): grapevine berries with and without Grapevine red blotch virus (PRJNA339456) and Arabidopsis leaves with and without *Botrytis cinerea* (PRJNA1270936). The five acceptance criteria are written in the fixture. Both public checks passed on the STAR route and on the Salmon route; job IDs and read counts are in [RELEASE.md](RELEASE.md). The checks do not compare the two routes with each other. For the virus, whose genes overlap and whose libraries are unstranded, per-gene counts differed between them: STAR with featureCounts leaves out reads in overlapping genes, and Salmon divides them between the overlapping transcripts.

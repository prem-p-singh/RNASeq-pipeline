# Public-data qualification

Public examples are authorized for development and validation. A successful bulk test does not qualify another assay or a different chemistry. The release capability matrix must cite completed tests for the exact producer, input stage and protocol it advertises.

## Executable bulk checks

`tests/check_public.sh` uses six distinct biological samples, SRR6357070–SRR6357075 from GSE110004. The [pinned nf-core source](https://raw.githubusercontent.com/nf-core/test-datasets/626c8fab639062eade4b10747e919341cbf9b41a/README.md) provides chromosome-I-filtered, downsampled yeast reads. These are operational smoke tests, not reproduction of the study's genome-wide conclusions.

- `bash tests/check_public.sh`: Salmon, gene aggregation, differential expression and reports.
- `bash tests/check_public.sh --star`: STAR/featureCounts, indexed alignments, raw gene counts, edgeR QL and reports. Requires the STAR module described in README.
- Add `--screen` to exercise declared expected-yeast screening on 10,000 paired fragments per sample, with separate read counts for each mate. This is an operational public check, not a known-contaminant sensitivity study. The independent synthetic mixed-panel check supplies known exclusive/shared/unmapped expectations.
- Set `RNASEQ_PUBLIC_PROJECT` to a new directory to retain downloaded inputs, provenance and outputs.

The STAR fixture takes reverse strandedness from the pinned upstream samplesheet. Biological replicate identities come from the accession-level source README; the upstream workflow-test grouping is not used as biological metadata. Both routes record reference/input identities. Synthetic exon-junction and strand controls establish known count expectations independently of this public smoke test.

## Next representative profiles

These are source candidates, **not implemented or qualified capabilities**. Download manifests, checksums and concrete acceptance results must be added when each route is built.

| Profile | Public source | Required comparison and boundary |
|---|---|---|
| Chromium 3′ v3 PBMC | [10x 1k PBMCs](https://www.10xgenomics.com/datasets/1-k-pbm-cs-from-a-healthy-donor-v-3-chemistry-3-standard-3-0-0) | Barcode/UMI extraction, raw and filtered sparse matrices, cell calling and QC against the published reference outputs. A single donor cannot qualify donor-level DE. |
| Visium fresh-frozen mouse brain | [10x aggregate of mouse-brain sections](https://www.10xgenomics.com/datasets/aggregate-of-mouse-brain-sections-visium-fresh-frozen-whole-transcriptome-1-standard) | Select individual source sections; preserve barcodes, coordinates, image scale and tissue masks. Validate registration before aggregation. Sections are not automatically independent animals. |
| PacBio Iso-Seq | [Official IsoSeq project](https://github.com/PacificBiosciences/IsoSeq) | Select a versioned example with a documented HiFi/FLNC input stage; compare isoform structures and classifications. Starting at FLNC cannot qualify the preceding CCS/demultiplexing stages. |
| ONT transcriptomes | [Official workflow demo](https://epi2me.nanoporetech.com/workflows/wf-transcriptomes/) | Freeze the demo/version and kit; compare alignment, transcript structures and counts. Basecalled reads do not qualify raw-signal modification analysis. |
| QuantSeq FWD, REV and UMI | [Manufacturer protocol documentation](https://www.lexogen.com/docs/quantseq/) | Choose separate public accessions with explicit kit versions, read orientation and UMI structure. No accession is qualified yet; kit variants cannot inherit one another's evidence. |
| Animal small RNA | [GSE94585 synthetic ratiometric pools](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE94585) | Select matching TruSeq technical libraries and published input ratios to measure count recovery and ratio bias. These technical pools cannot qualify biological replication. |
| Plant small RNA, dual RNA and specialized assays | Accession selection remains open | Match the handbook profile and required controls before implementing a producer. Include organism assignment/ambiguity for dual RNA, length/adapter/UMI controls for small RNA, and protocol-specific controls for specialized assays. |

### Frozen biological expectations for the first cell example

The official 10x PBMC v3 page specifies one donor, 1,222 called cells, a 28-base R1 (16-base cell barcode plus 12-base UMI), 91-base transcript R2 and an 8-base I7 sample index; its published outputs used Cell Ranger 3.0.0 with `--expect-cells=1000`. These are comparison metadata for this particular example, not universal 10x defaults. Download URLs, input hashes and version-specific comparison tolerances must still be frozen before execution. [Dataset source](https://www.10xgenomics.com/datasets/1-k-pbm-cs-from-a-healthy-donor-v-3-chemistry-3-standard-3-0-0).

### Selected TruSeq technical controls

The [pinned fixture manifest](../tests/fixtures/small_rna_public.json) selects GSE94585 Lab1 SynthA/SynthB technical replicate 1: **SRR5234383 / GSM2478899** and **SRR5234463 / GSM2478919**. GEO metadata, ENA run aliases/file sizes/MD5s, the publication's supplementary workbook and published exceRpt count tables establish their identities. The reference source is Table S2 of [Giraldez et al.](https://doi.org/10.1038/nbt.4183); the adapter is taken from [Illumina's TruSeq Small RNA documentation](https://support-docs.illumina.com/SHARE/AdapterSequences/Content/TruSeq-SmallRNA.htm).

`python3 tests/fixtures/prepare_small_rna_public.py NEW_DIRECTORY` downloads the pinned files and verifies identities. It derives FASTA and expected pool ratios while recording three duplicated sequence IDs with conflicting alias/biotype labels. The 334 table rows contain 331 distinct IDs/sequences; duplicates retain their original labels in the derivation record. Linux fixture preparation passed in job 38994949, including both FASTQ size/MD5 checks and the pinned source/count-table SHA256 checks. This preparation step is **not a small-RNA analysis route or qualification pass**. The fixture includes 15–90 nt synthetic targets, so a universal 18–30 nt filter would discard intentional controls. These two technical libraries cannot qualify biological differential expression.

## Evidence recorded for each profile

1. Accession, source URL/revision, retrieval date, file size and SHA256; document any deterministic subsetting or reference derivation.
2. Organism, kit/chemistry, input stage, layout, strand, barcode/UMI structure and biological-unit mapping.
3. Exact source revision, environment locks, reference identity, commands and resources.
4. Known synthetic controls plus an independent public output or direct native-tool comparison; state the tolerance before evaluating it.
5. Fresh execution, ordinary resume, missing-output recovery and changed-input invalidation.
6. A capability-specific result: passed, failed, unavailable or untested. Installation, download and job submission are not scientific passes.

All required checks for an advertised release route must run without skips. A reduced public dataset qualifies its tested behavior; larger-cohort performance, biological inference and other protocols need their own evidence.

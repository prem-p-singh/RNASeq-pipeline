# Intake template audit and v2 release gate

Updated: 2026-09-23. Scope: **v2 bulk intake is executable; full multi-assay intake remains incomplete.**

The owner has authorized v2 publication for the current bulk workflow, superseding the earlier full-roadmap publication gate. Template schema version 3 identifies the workbook format; software release v2.0.0 identifies the pipeline. Future assay families require their own executable intake and scientific qualification before being advertised as supported.

## Canonical intake implementation update

Schema 3 adds Samples, Libraries and Reads sheets and an explicit source selector. External-file intake remains compatible. Workbook-table intake records sample/library ownership, text IDs, covariates, run/lane/role, URI and optional byte size/SHA256. Setup writes the three canonical TSVs, and the workflow consumes them directly. Multiple lanes/runs within one non-UMI bulk library are validated and merged with matching mate IDs/counts; source files remain untouched. Per-library read provenance is a required workflow output. Storage estimates include all lanes and the merged copy and remain advisory.

This closes the previously missing canonical bulk read connection. It does not enable multiple prepared libraries per specimen, repeated-measures Excel intake, UMI processing or other assay families. Record-sheet rows start at 5 and may be extended beyond the preformatted rows. Scientific covariates may be added to Samples; input identifiers must be literal text. Conflicting external paths/table records, malformed table relationships and unresolved library eligibility are rejected.

Validation: the schema/parser/generator field checks and both external/inline setup tests passed. The three new sheets and changed questionnaire views were rendered and visually reviewed. Linux job 38577940 passed all 22 strict checks with no failures/skips, including real synthetic FASTQ-to-report processing, lane-source changes and provenance recovery. The final resource delta passed a separate Linux check. These establish the tested bulk scope, not the remaining family-specific roadmap.

## What was checked

The audit compares `intake_template.xlsx`, `config/intake_questions.yaml`, the workbook generator, `scripts/intake.py`, `scripts/setup.py`, metadata validation, configuration defaults and the executable workflow. The broader handbook-driven roadmap describes intended coverage; a questionnaire or validation helper alone is not evidence that an assay runs. Operational instructions are in [the input guide](../INPUTS.md) and [README](../README.md).

## Findings repaired in this pass

| Finding | Repair | Boundary |
|---|---|---|
| README incorrectly said nothing parses the workbook | Replaced with direct `setup.py --intake` instructions, metadata/runtime prerequisites and review/launch steps | Direct Excel setup remains independent bulk only |
| No visible schema version or populated machine field IDs | Added schema marker at `README!B2` and field IDs in column E; generator and parser agree | Schema 3 is a development input format, not a release declaration |
| Assay selection depended on the displayed Project name label | Selection accepts the stable `project_name` ID | Complete-row sorting and label changes preserve interpretation |
| No way to detect an unsupported UMI library at intake | Added required UMI declaration; yes/unknown blocks this bulk adapter | No UMI producer or deduplication was implemented |
| Input processing stage was implicit | Added required stage selection; only raw FASTQ is accepted | Count/object and instrument-data imports remain unavailable |
| Reverse strandedness was preselected without user evidence | Default changed to don't know; no orientation is invented | Library-type detection and disagreement reporting remain the current behavior |
| TAGseq kit/orientation defaults could contradict each other | Removed automatic kit/orientation choices and marked the tab planning-only | Exact kit profiles still require implementation and qualification |
| Help claimed email, tissue, platform, sample N and kit fields controlled unimplemented behavior | Corrected help and made contact optional | Context-only answers stay in provenance but do not configure processing |
| OrgDb help promised general coverage and a missing package field | Help now states limits and the YAML requirement for use_existing | Existing-package selection still cannot be completed entirely in Excel |
| Kit identity could not be recorded on the bulk tab | Added optional exact kit/version text | It is context only, not a qualified kit profile |

Legacy workbooks without a version marker remain schema 1. They lack the new stage/UMI declarations; setup explicitly warns about the legacy raw non-UMI assumption. Current templates require explicit answers. No formula evaluation, macro execution or automatic repair of date-like IDs is introduced.

## Current completeness

| Input area | Current status | Still required for the planned intake |
|---|---|---|
| Project and biological question | Present | Shared study record instead of repeating answers across assay tabs |
| Assay and objective | Bulk selected by worksheet; gene expression/DE are implicit | Explicit assay/subtype, input stage, objective list, per-analysis identity and compatible tool preference |
| Organism | One taxID or recognized local preset name | Multiple taxa, organism roles and species-specific reference/annotation identities |
| Biological sample records | External files or inline Samples with text IDs and covariates | Qualified donor/specimen/repeated-unit semantics for later routes |
| Libraries | Linked Libraries records with layout, strand, UMI declaration and protocol context | Exact processing profiles, barcode/UMI geometry and multiple libraries per sample |
| Reads and other inputs | Explicit library/run/lane/role/URI/checksum/size, or legacy folder scan | Object/image/signal inputs and stage-aware producers for relevant routes |
| Read relationships | Canonical records drive bulk DAG; paired order/count checks and within-library lane merging | Qualified non-bulk geometry and multiple-library count aggregation |
| Independent-group design | Primary factor, optional batch and smallest-group check | Explicit independent unit, model terms, selected levels/subsets, contrast direction and reference levels |
| Paired/repeated designs | Workbook detects and rejects; YAML path exists | Subject ID and explicit fixed/random terms, dependence-preserving validation and route-specific qualification |
| Complex analyses | Configuration supports limited pairwise contrasts | Analyses table for interactions, continuous/time-course effects and objective-specific responses; validated adapters |
| References | Setup resolves online and asks for generated-config review | User-selected accession/version, transcriptome/genome/GTF paths, matching identifiers, decoys and annotation/control resources |
| Enrichment | Combined GO/KEGG toggle and partial OrgDb policy | Package identity, key type, gene-ID mapping, gene sets/background, separate methods, optional module outcomes |
| Coexpression | WGCNA toggle; current gate/defaults apply | Eligible cohort/object selection independent of DE, covariates/traits and exposed approved parameters |
| QC and contamination | Workflow defaults, not full intake control | Expected organisms, contaminants/spike-ins, screen/remove policy, requested QC, inclusion/review policy and applicable controls |
| Execution | Project directory and launcher settings outside workbook | Local/SLURM choice, site profile, cores/RAM/jobs, runtime/reference/cache/scratch paths, retention and offline policy |
| Storage | Advisory estimates implemented | Dataset/read-unit-aware estimates for every added route; no fixed platform storage cap |
| Reproducibility | Workbook checksum/source map and resolved setup inputs | Complete external table/input hashes, reference/tool identities, per-analysis decisions and configuration revision handling |

The workbook is usable for the narrow supported bulk adapter, with external or inline sample metadata and configuration review. It does **not** contain everything needed to drive the full roadmap, and it does not expose every supported advanced configuration option.

## Assay-specific fields and execution still missing

The following are planned intake contracts, not newly implemented scientific support. Requirements depend on the selected objective and exact protocol; unrelated fields should not become mandatory for every study.

| Planned family | Additional intake needed | Execution gate |
|---|---|---|
| Bulk alternative measurements | Gene/transcript/alignment objectives; reference and aligner preference; strand/counting/retention policy | Qualified alignment/count adapters, requested outputs and alternative statistical backends |
| Bulk end-tag/TAG-seq | Exact kit/version and FWD/REV geometry; UMI/spacer fields when applicable; end-site objective and internal-priming policy | Kit-specific preprocessing/counting with correct count semantics |
| Small RNA: plant/animal | Kit adapters, random bases/UMIs, size classes, organism class, target annotations and selected miRNA/siRNA/isomiR/phasing objectives | Distinct qualified producers and objective-specific outputs |
| Viral small RNA | Host/virus references, discovery versus known targets, controls and expected organisms | Assignment/discovery and read-back evidence; no silent host/pathogen filtering |
| Single-cell/single-nucleus | Chemistry/version, barcode whitelist and read geometry, UMI layout, donor/specimen, input stage, cell-calling/intron policy and requested analyses | Barcode-aware producer, cell QC and donor-aware inference |
| Spatial sequencing | Platform/chemistry, specimen/slide/section IDs, coordinates/images/registration, spot/bin geometry and reference | Native producer, matched matrix/spatial objects and specimen-aware analysis |
| Spatial imaging | Vendor/export format, panel/features, segmentation, coordinate system and specimen links | Validated imaging-object import and panel-aware spatial analysis |
| PacBio RNA | Input stage, chemistry/primer information, sample demultiplex metadata, reference and isoform objectives | Stage-aware processing, transcript classification and quantification |
| ONT cDNA/direct RNA | Library type/chemistry, raw signal versus reads, basecaller/model provenance, reference; signal links for tail/modification objectives | Compatible long-read producer and separately qualified signal-dependent modules |
| Dual organism | Species roles, reference namespaces, intended mixture, joint/sequential assignment choice and per-species models | Ambiguity accounting and species-separated measurements/inference |
| Ribo-seq | Footprint/adapter protocol, size classes, depletion references, P-site calibration and matched RNA samples | Periodicity/ORF diagnostics and validated translation-efficiency model |
| CLIP | Exact protocol, UMIs/barcodes, crosslink/read geometry, input/control sample links and replicate groups | Protocol-specific sites/peaks and control/replicate qualification |
| GRO/PRO | Protocol, strand/adapter geometry, spike-ins/controls and transcription/pausing objectives | Nascent-signal processing and objective-specific models |
| TT/SLAM/BRU | Label/enrichment protocol, timing/dose, controls, conversion expectations and kinetic objective | Compatible producer, efficiency diagnostics and response model |
| CAGE | End/read geometry, strand, promoter/TSS objectives and reference | End-aware mapping, TSS clustering and differential usage |
| PAS/3P | End-capture protocol, UMI geometry, internal-priming policy and site/usage objectives | Protocol-qualified site calling/counting and usage analysis |
| Reference-poor/de novo | Reference availability, assembly/annotation objectives, expected taxa and contamination controls | Qualified assembly/import, feature assessment and uncertainty reporting |
| Counts/producer-object import | Producer/version, feature/sample IDs, count units, normalization/offsets, layers and previous processing | Semantic validation that prevents double normalization or invalid downstream models |

## Required next implementation order

1. Establish the shared workbook tables: Study, Samples, Libraries, Reads/Inputs, Analyses, References/Controls and Execution. Preserve the recognizable assay questions as conditional details. Add fields only alongside a defined consumer and validation policy.
2. Connect these records through normalization, recommendation and actual DAG inputs. Reject conflicting external/inline sources, invalid joins, duplicate files and unsupported read structures. Do not count lanes, cells or technical libraries as biological replicates.
3. Complete handbook-aligned raw preprocessing/QC and contamination accounting for the first supported routes. Platform and exact kit selections must actually control their documented operations.
4. Implement and qualify every family in the table, including its requested measurement and downstream response. Expose executable choices only when the route and its input contract exist.
5. Extend model/contrast, reference, annotation and execution controls; keep optional/unavailable modules and actual failures distinguishable.
6. Run per-family fixtures and independent scientific comparisons, then test runtime installation, restart, lineage and generic Linux/SLURM deployment. Reconcile workbook, schema, parser, generator, docs and recommendation output.

## Publication gate

The following checklist governs future full multi-assay coverage. It is not a claim that v2 implements those families. Each additional family needs executable intake-to-result paths and recorded qualification:

- Matching generated workbook/schema/parser with stable IDs and typed literal input validation.
- Complete read/sample/analysis relationships and source provenance.
- Per-route declared support and limitations, with positive and negative fixtures.
- Verified environment/reference identities and reviewed configuration defaults.
- Passing required Linux checks with no skipped release gates, plus the relevant scientific benchmarks.
- A source freeze and release notes describing the exact tested commit; only then create the v2 tag and GitHub release.

Targeted validation passed: `tests/check_intake.py`, `tests/check_setup_helpers.py`, Python syntax checks and `git diff --check`. The intake check verifies schema/workbook field agreement, row sorting, relabeled project selection, literal sample IDs, formula/dropdown rejection, unknown schema rejection, UMI/stage eligibility, generated configuration and preservation of the source workbook. Setup's reference services are mocked in this check. The saved workbook was inspected and changed views were rendered for readability.

The later canonical execution continuation passed the 22-check strict Linux suite, as recorded above. Template checks alone do not establish biological qualification. All planned assay routes must pass their own qualification before they can be added to the release's supported capability table.

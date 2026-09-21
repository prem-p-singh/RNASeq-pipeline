# Project setup — the easy path

Three ways to start a new project, listed from **least effort** to most.

## 🚀 Easiest: drag-and-drop via `start_new.sh`

From your Mac, one command uploads your metadata (and optionally FASTQs),
SSHes to FARM, and runs an interactive wizard. No YAML editing, no taxID
lookups, no remembering column names.

```bash
# Just metadata — FASTQs already on FARM or pulled from a URL you'll paste
scripts/start_new.sh my_project ~/Downloads/my_metadata.xlsx

# With local FASTQs to upload
scripts/start_new.sh my_project ~/data/meta.csv ~/data/fastq/
```

On FARM, you'll be prompted for (all with smart defaults):

| Prompt | What it does | Typical keystrokes |
|---|---|---|
| Organism | Menu: grape / human / mouse / rice / arabidopsis / other | `1` + Enter |
| Sample-ID column | Auto-guessed from your metadata | Enter |
| Primary factor | Auto-listed + guessed | `1` + Enter |
| Sequencing type | Menu: tagseq / rnaseq_single / rnaseq_paired | `1` + Enter |
| Proceed? | Confirms, runs everything | `y` |

That's it. Total interaction: ~6 keystrokes + one `y`.

## ⚙️ Medium: manual YAML (the old path)

Fill in `setup_inputs.yaml` by hand if you want full control:

```bash
cp scripts/setup_inputs.template.yaml setup_inputs.yaml
$EDITOR setup_inputs.yaml
python scripts/setup.py setup_inputs.yaml
./submit.sh
```

## 🔧 Full manual

Write `config/config.yaml` and `config/samples.tsv` by hand, skip the wizard
entirely.

---

# Details: what the wizard generates (all three paths converge here)

---

## What you need to provide (7 fields)

Copy `scripts/setup_inputs.template.yaml` to `setup_inputs.yaml` in the repo root and fill in:

| # | Field | Example |
|---|---|---|
| 1 | `project_name` | `"rice_heat_stress"` |
| 2 | `tax_id` | `4530` (Oryza sativa) |
| 3 | `seq_type` | `"rnaseq_paired"` |
| 4 | `fastq_source` | `"/scratch/fastq/"` or `"s3://bucket/project/fastq/"` |
| 5 | `metadata_file` | `"~/metadata.xlsx"` (or `.csv`) |
| 6 | `sample_id_column` | `"Sample"` |
| 7 | `primary_factor` | `"treatment"` |

Plus two optional fields: `model_fixed_effects` (default `"~ treatment"`) and `model_random_effects` (default `null`).

---

## What you get out

Running `python scripts/setup.py setup_inputs.yaml` produces:

```
config/
  config.yaml               ← fully populated, ready for ./submit.sh
  samples.tsv               ← one row per matched sample
reference/
  annotation_info.tsv       ← pulled from NCBI gene_info
  NCBI_to_ensembl.txt       ← pulled from Ensembl BioMart (optional)
  org.XXX.eg.db/            ← built for non-model organisms via AnnotationForge
```

And appends a line to `gates/decisions.log` for every step.

---

## What the wizard does behind the scenes

| Stage | Action | Source |
|---|---|---|
| Resolve organism | taxID → scientific name → OrgDb package name → KEGG code | NCBI Taxonomy eutils + KEGG REST |
| Resolve reference | taxID → latest RefSeq assembly accession → FASTA/GTF/gene_info URLs | NCBI Datasets API |
| Fetch annotation | Download `*_gene_info.gz`, project to our 4-column schema | NCBI FTP |
| Fetch Ensembl map | Query BioMart for `entrezgene_id ↔ ensembl_gene_id` | Ensembl BioMart |
| Match FASTQs to metadata | Enumerate `fastq_source`, substring-match to `sample_id_column` | local / S3 / HTTP listing |
| Install OrgDb | If Bioconductor has it → print `BiocManager::install(...)` command. Else → invoke `build_orgdb.R` → AnnotationForge builds from NCBI | Bioconductor / AnnotationForge |
| Render config | Merge inputs into `config.template.yaml` | — |

---

## Two ways to run it

**Declarative (recommended):**
```bash
cp scripts/setup_inputs.template.yaml setup_inputs.yaml
$EDITOR setup_inputs.yaml                 # fill in the 7 fields
python scripts/setup.py setup_inputs.yaml
```

**Interactive (no YAML needed):**
```bash
python scripts/setup.py --interactive
# Prompts you for each field at the terminal
```

---

## What happens when something fails

Network calls (NCBI Datasets, Ensembl BioMart, KEGG) are best-effort. If one fails, the wizard:

1. Logs a **WARN** line to stdout and `gates/decisions.log`
2. Writes the output file with placeholders (e.g. `"TODO"` URLs, empty annotation)
3. Continues — downstream stages degrade gracefully:
   - Missing `annotation_info.tsv` → DE tables have no `Description` column
   - Missing `NCBI_to_ensembl.txt` → no Ensembl IDs in outputs
   - Missing `kegg_code` → KEGG stage auto-skips
   - OrgDb not built → GO enrichment errors (only real hard-fail)

You can re-run the wizard any time; it overwrites.

---

## Skipping the wizard

If you already have `config/config.yaml` + `config/samples.tsv` written by hand, you don't need `setup.py` at all. The wizard is strictly a convenience.

---

## Worked example (grape)

The grape GRBV project — preserved in `examples/grape/` — would be regenerated from:

```yaml
project_name: "grape_grbv"
tax_id: 29760
seq_type: "tagseq"
fastq_source: "s3://ucd-core/singh_grape/fastq/"
metadata_file: "~/Correct.Sample_Metadata.xlsx"
sample_id_column: "Sample"
fastq_pattern: "{sample_id}.fastq.gz"
primary_factor: "treatment"
model_fixed_effects: "~ group * stage"
model_random_effects: "(1|vine)"
storage_budget_gb: 20
```

`python scripts/setup.py` would produce the same `config/config.yaml` + `config/samples.tsv` the hand-written grape pipeline used.

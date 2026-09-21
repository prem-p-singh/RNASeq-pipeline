#!/usr/bin/env python3
"""
setup.py — Central project setup wizard.

Takes 7 minimal inputs (see scripts/setup_inputs.template.yaml) and produces
every external file the pipeline needs to run:

    config/config.yaml
    config/samples.tsv
    reference/annotation_info.tsv       (from NCBI)
    reference/<org.XXX.eg.db>/          (built via AnnotationForge if needed)

Usage:
    python scripts/setup.py setup_inputs.yaml
    python scripts/setup.py --interactive
    python scripts/setup.py                     # uses setup_inputs.yaml if present

Design notes
------------
Network calls are best-effort. If an NCBI/Ensembl endpoint fails, the wizard
logs the failure, emits an empty placeholder, and prints a fallback command
the user can run manually. The pipeline does not hard-fail on missing
annotation/ensembl files — only the gene-description merge is skipped.
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import re
import shutil
import subprocess
import sys
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime
from pathlib import Path
from typing import Any

import yaml

# --------------------------------------------------------------------- paths
ROOT = Path(__file__).resolve().parent.parent
CONFIG_DIR = ROOT / "config"
PRESETS_DIR = ROOT / "scripts" / "presets"
REFERENCE_DIR = ROOT / "reference"
TEMPLATE_CFG = CONFIG_DIR / "config.template.yaml"
DECISIONS_LOG = ROOT / "gates" / "decisions.log"


# ======================================================================= I/O
def log(stage: str, msg: str, level: str = "INFO"):
    ts = datetime.now().isoformat(timespec="seconds")
    print(f"[{ts}] {level:>5}  {stage}: {msg}", flush=True)
    DECISIONS_LOG.parent.mkdir(parents=True, exist_ok=True)
    with open(DECISIONS_LOG, "a") as f:
        f.write(f"[{ts}] SETUP/{stage}: {msg}\n")


def load_yaml(path: Path) -> dict:
    with open(path) as f:
        return yaml.safe_load(f)


def dump_yaml(obj: dict, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        yaml.safe_dump(obj, f, sort_keys=False, default_flow_style=False)


# ============================================================ input gathering
REQUIRED_FIELDS = [
    "project_name", "tax_id", "seq_type", "fastq_source",
    "metadata_file", "sample_id_column", "primary_factor",
]

def gather_inputs(args) -> dict:
    if args.inputs_file and Path(args.inputs_file).exists():
        inputs = load_yaml(Path(args.inputs_file))
        log("INPUTS", f"loaded from {args.inputs_file}")
    elif (ROOT / "setup_inputs.yaml").exists() and not args.interactive:
        inputs = load_yaml(ROOT / "setup_inputs.yaml")
        log("INPUTS", "loaded from setup_inputs.yaml")
    else:
        inputs = prompt_inputs()

    missing = [f for f in REQUIRED_FIELDS if not inputs.get(f)]
    if missing:
        sys.exit(f"Missing required fields: {', '.join(missing)}")
    return inputs


def prompt_inputs() -> dict:
    def ask(q, default=None):
        suffix = f" [{default}]" if default is not None else ""
        ans = input(f"{q}{suffix}: ").strip()
        return ans or default

    print("--- Interactive setup ---")
    return {
        "project_name":          ask("Project name"),
        "tax_id":                int(ask("NCBI taxID")),
        "seq_type":              ask("Seq type (tagseq/rnaseq_single/rnaseq_paired)", "tagseq"),
        "fastq_source":          ask("FASTQ source (dir / s3:// / http listing)"),
        "metadata_file":         ask("Metadata spreadsheet path (.csv/.xlsx)"),
        "sample_id_column":      ask("Sample-ID column in metadata", "Sample"),
        "fastq_pattern":         ask("FASTQ filename pattern (use {sample_id}; blank = auto)") or None,
        "primary_factor":        ask("Primary factor (column used for DE contrasts)"),
        "model_fixed_effects":   ask("Fixed-effects formula", "~ treatment"),
        "model_random_effects":  ask("Random-effects term (e.g. (1|donor); blank = none)") or None,
        "storage_budget_gb":     int(ask("HPC storage budget GB", "20")),
    }


# ======================================================== organism resolution
def preset_for_taxid(tax_id: int) -> dict:
    """Curated organism facts from scripts/presets/, keyed by NCBI taxID.

    Returns {} when no preset matches, which the caller treats as "guess and
    warn" rather than "assume correct".
    """
    if not PRESETS_DIR.is_dir():
        return {}
    for path in sorted(PRESETS_DIR.glob("*.yaml")):
        try:
            d = load_yaml(path) or {}
            if int(d.get("tax_id", -1)) == int(tax_id):
                return d
        except (ValueError, TypeError, OSError):
            continue
    return {}


def resolve_organism(tax_id: int) -> dict:
    """
    Look up genus, species, KEGG code, and Bioconductor OrgDb name.

    Strategy:
      1. Hit NCBI Taxonomy efetch (XML) for scientific_name + lineage.
      2. Derive OrgDb package name by Bioconductor convention:
         org.<GenusFirstLetter><species>.eg.db   (e.g. Vvinifera -> org.Vvinifera.eg.db)
         For common organisms this matches an existing Bioconductor package.
      3. Query KEGG REST for the org code (ttp://rest.kegg.jp/list/organism).

    Failures return best-effort defaults rather than crashing.
    """
    log("ORGANISM", f"resolving taxID {tax_id}")
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
        f"db=taxonomy&id={tax_id}&retmode=xml"
    )
    genus = species = None
    try:
        with urllib.request.urlopen(url, timeout=20) as r:
            xml = r.read().decode()
        m = re.search(r"<ScientificName>([^<]+)</ScientificName>", xml)
        if m:
            parts = m.group(1).split()
            if len(parts) >= 2:
                genus, species = parts[0], parts[1]
    except Exception as e:
        log("ORGANISM", f"NCBI taxonomy lookup failed ({e}); using placeholders", "WARN")

    genus = genus or "Genus"
    species = species or "species"

    # Bioconductor OrgDb names are not derivable from the species name:
    # Homo sapiens is org.Hs.eg.db, not org.Hsapiens.eg.db. The curated
    # presets are the only reliable mapping, so consult them by taxID first
    # and only fall back to a constructed name (flagged as a guess).
    preset = preset_for_taxid(tax_id)
    if preset:
        orgdb_pkg = preset.get("orgdb_package") or ""
        log("ORGANISM", f"matched preset '{preset.get('common_name', '?')}' "
                        f"-> {orgdb_pkg}")
    else:
        orgdb_pkg = f"org.{genus[0]}{species}.eg.db"
        log("ORGANISM",
            f"no preset for taxID {tax_id}; guessed OrgDb name '{orgdb_pkg}' "
            f"which may not match Bioconductor. Add a file to scripts/presets/ "
            f"to pin it.", "WARN")

    # KEGG code lookup (tiny text file, ~1 MB)
    kegg_code = preset.get("kegg_code") if preset else None
    if kegg_code is None:
        try:
            with urllib.request.urlopen("https://rest.kegg.jp/list/organism",
                                        timeout=20) as r:
                for line in r.read().decode().splitlines():
                    cols = line.split("\t")
                    if len(cols) >= 3 and f"{genus} {species}" in cols[2]:
                        kegg_code = cols[1]
                        break
        except Exception as e:
            log("ORGANISM", f"KEGG lookup failed ({e}); will skip KEGG", "WARN")

    log("ORGANISM",
        f"genus={genus} species={species} orgdb={orgdb_pkg} kegg={kegg_code or 'None'}")
    return {
        "tax_id": tax_id, "genus": genus, "species": species,
        "orgdb_package": orgdb_pkg, "kegg_code": kegg_code,
    }


# ================================================= NCBI reference URL lookup
def resolve_reference_urls(tax_id: int) -> dict:
    """Query NCBI Datasets API for the latest RefSeq assembly + transcriptome + GTF."""
    log("REFERENCE", f"querying NCBI Datasets for taxID {tax_id}")
    api = (
        f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/taxon/{tax_id}/"
        "dataset_report?filters.assembly_source=refseq&filters.reference_only=true"
        "&page_size=1"
    )
    try:
        with urllib.request.urlopen(api, timeout=30) as r:
            data = json.loads(r.read().decode())
        rpt = data["reports"][0]
        acc = rpt["accession"]
        # NCBI's FTP directory is <accession>_<assembly_name> and every file in
        # it carries that same stem. assembly_name comes back in this very
        # report, so exact filenames need no directory listing and no globbing.
        name = (rpt.get("assembly_info") or {}).get("assembly_name")
        if not name:
            raise KeyError("assembly_info.assembly_name missing from NCBI report")
        stem = f"{acc}_{name}"
        digits = acc.split("_")[1].split(".")[0]        # GCF_030704535.1 -> 030704535
        mid = "/".join(digits[i:i + 3] for i in range(0, 9, 3))
        base = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{acc[:3]}/{mid}/{stem}"
        log("REFERENCE", f"assembly {acc} ({name})")
        return {
            "accession": acc,
            "assembly_name": name,
            "genome_fasta_url":        f"{base}/{stem}_genomic.fna.gz",
            "transcriptome_fasta_url": f"{base}/{stem}_rna.fna.gz",
            "gtf_url":                 f"{base}/{stem}_genomic.gtf.gz",
            "gene_info_url":           f"{base}/{stem}_gene_info.gz",
        }
    except Exception as e:
        log("REFERENCE",
            f"NCBI Datasets lookup failed ({e}); fill reference.*_url in config.yaml manually",
            "WARN")
        return {"accession": None, "assembly_name": None,
                "genome_fasta_url": "TODO",
                "transcriptome_fasta_url": "TODO", "gtf_url": "TODO",
                "gene_info_url": None}


def unresolved_reference_urls(ref: dict) -> list:
    """Problems that would make a fetch job fail. Empty list means good to go.

    curl does not expand '*', so a wildcard URL is a guaranteed Stage 0 failure;
    catch it here instead of at job time.
    """
    problems = []
    for key in ("genome_fasta_url", "transcriptome_fasta_url", "gtf_url"):
        url = ref.get(key) or ""
        if not url or url == "TODO":
            problems.append(f"{key} is unresolved")
        elif "*" in url:
            problems.append(f"{key} still contains a wildcard: {url}")
    return problems


# ============================================= annotation_info.tsv generation
def fetch_annotation_info(gene_info_url: str | None, out_path: Path) -> bool:
    """Download NCBI gene_info, project to our 4-column schema."""
    if not gene_info_url or gene_info_url == "TODO":
        log("ANNOTATION", "no gene_info URL available; skipping", "WARN")
        return False
    log("ANNOTATION", f"fetching {gene_info_url}")
    try:
        tmp = out_path.with_suffix(".tsv.gz")
        urllib.request.urlretrieve(gene_info_url, tmp)
        # NCBI gene_info header: #tax_id GeneID Symbol LocusTag Synonyms dbXrefs ...
        # We produce:  Gene.stable.ID  DB.xref  Description  Molecule.Type
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(tmp, "rt") as fin, open(out_path, "w") as fout:
            fout.write("Gene.stable.ID\tDB.xref\tDescription\tMolecule.Type\n")
            for line in fin:
                if line.startswith("#"):
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 10:
                    continue
                gene_id, symbol, description, mol_type = cols[1], cols[2], cols[8], cols[9]
                fout.write(f"{symbol}\tGeneID:{gene_id}\t{description}\t{mol_type}\n")
        tmp.unlink()
        log("ANNOTATION", f"wrote {out_path}")
        return True
    except Exception as e:
        log("ANNOTATION", f"failed ({e}); skipping", "WARN")
        return False


# ====================================================== OrgDb: report only
# Names Bioconductor actually ships. Kept explicit because these cannot be
# derived from a species name.
BIOCONDUCTOR_ORGDBS = {
    "org.Hs.eg.db", "org.Mm.eg.db", "org.Rn.eg.db", "org.Dm.eg.db",
    "org.Ce.eg.db", "org.Sc.sgd.db", "org.At.tair.db", "org.Dr.eg.db",
    "org.Gg.eg.db", "org.Ss.eg.db", "org.Bt.eg.db", "org.Cf.eg.db",
}


def setup_orgdb(org: dict) -> str:
    """Report how the OrgDb will be obtained. Builds nothing.

    The build belongs to the `build_orgdb` Snakemake rule, which owns the
    tiered fallback and passes the arguments scripts/build_orgdb.R actually
    accepts. Having a second caller here meant two code paths to keep in
    sync, and this one drifted: it passed --out_dir, which the script has
    never accepted, so the "non-model organism" route always failed.
    """
    pkg = org["orgdb_package"]
    if pkg in BIOCONDUCTOR_ORGDBS:
        log("ORGDB", f"{pkg} is a Bioconductor package")
        return "bioconductor"
    log("ORGDB", f"{pkg} is not a stock Bioconductor package; the build_orgdb "
                 f"rule will construct it at run time")
    return "will_build"


# =========================================== samples.tsv from FASTQ + metadata
def list_fastqs(source: str) -> list[str]:
    if source.startswith("s3://"):
        cp = subprocess.run(
            ["aws", "s3", "ls", source], capture_output=True, text=True, check=True
        )
        files = [line.split()[-1] for line in cp.stdout.splitlines() if line.strip()]
        return [source.rstrip("/") + "/" + f for f in files if f.endswith((".fastq.gz", ".fq.gz"))]
    if source.startswith(("http://", "https://")):
        with urllib.request.urlopen(source, timeout=30) as r:
            html = r.read().decode()
        return [
            urllib.parse.urljoin(source, m.group(1))
            for m in re.finditer(r'href="([^"]+\.f(?:ast)?q\.gz)"', html)
        ]
    # local
    src = Path(source).expanduser().resolve()
    return [str(p) for p in src.rglob("*.f*q.gz")]


# Mate markers, most explicit first. The bare _1/_2 forms are only honoured
# immediately before the extension, which is where that convention actually
# appears; matching them anywhere would misread "Lane_1_S3_R1_001.fastq.gz".
_MATE_PATTERNS = [
    (2, re.compile(r"[._]R2[._]")),
    (1, re.compile(r"[._]R1[._]")),
    (2, re.compile(r"[._]2\.f")),
    (1, re.compile(r"[._]1\.f")),
]


def mate_of(basename: str):
    """Return 1, 2, or None for a FASTQ filename."""
    for mate, pat in _MATE_PATTERNS:
        if pat.search(basename):
            return mate
    return None


def files_for_sample(sid: str, fastqs: list) -> list:
    """FASTQs whose basename contains `sid` as a delimited token.

    Requiring a non-alphanumeric boundary (or start/end of name) is what stops
    sample "S1" from also claiming "S10_R1.fastq.gz".
    """
    pat = re.compile(rf"(^|[^A-Za-z0-9]){re.escape(sid)}([^A-Za-z0-9]|$)")
    return [f for f in fastqs if pat.search(os.path.basename(f))]


def build_samples_tsv(inputs: dict, out_path: Path):
    log("SAMPLES", f"scanning {inputs['fastq_source']}")
    fastqs = list_fastqs(inputs["fastq_source"])
    log("SAMPLES", f"found {len(fastqs)} FASTQ files")

    import pandas as pd

    # Load metadata. TSV must not be parsed as CSV or every column collapses
    # into one and the sample-ID lookup fails with a confusing error.
    mpath = Path(inputs["metadata_file"]).expanduser()
    suffix = mpath.suffix.lower()
    if suffix in {".xlsx", ".xls"}:
        try:
            import openpyxl  # noqa: F401
        except ImportError:
            sys.exit("openpyxl + pandas required to read .xlsx")
        mdf = pd.read_excel(mpath)
    elif suffix in {".tsv", ".tab", ".txt"}:
        mdf = pd.read_csv(mpath, sep="\t")
    else:
        mdf = pd.read_csv(mpath)

    sid_col = inputs["sample_id_column"]
    if sid_col not in mdf.columns:
        sys.exit(f"Column '{sid_col}' not found in {mpath}. "
                 f"Available: {list(mdf.columns)}")

    mdf[sid_col] = mdf[sid_col].astype(str)
    dups = mdf[sid_col][mdf[sid_col].duplicated()].unique().tolist()
    if dups:
        sys.exit(f"duplicate sample IDs in {mpath}: {', '.join(dups)}")

    paired = inputs["seq_type"] == "rnaseq_paired"
    pattern = inputs.get("fastq_pattern")

    # Resolve every sample. Problems are collected and reported together: a
    # metadata row is never silently dropped, because a sample vanishing from
    # the sheet is indistinguishable from one that was never sequenced.
    assignments: dict = {}
    problems: list = []

    for sid in mdf[sid_col]:
        hits = files_for_sample(sid, fastqs)
        if pattern:
            expected = pattern.format(sample_id=sid)
            narrowed = [f for f in hits if f.endswith(expected)]
            if narrowed:
                hits = narrowed
        if not hits:
            problems.append(f"{sid}: no FASTQ matched")
            continue

        m1 = [f for f in hits if mate_of(os.path.basename(f)) == 1]
        m2 = [f for f in hits if mate_of(os.path.basename(f)) == 2]
        unmated = [f for f in hits if mate_of(os.path.basename(f)) is None]

        if paired:
            if len(m1) != 1 or len(m2) != 1:
                extra = ""
                if len(m1) > 1 or len(m2) > 1:
                    # ponytail: one read pair per sample. Multi-lane merging is
                    # a sample-sheet schema change (one row per read unit);
                    # until then, merge lanes before running.
                    extra = (" — multiple lanes are not supported; "
                             "concatenate them into one pair first")
                problems.append(
                    f"{sid}: expected exactly one R1 and one R2, found "
                    f"{len(m1)} R1 and {len(m2)} R2{extra}")
                continue
            assignments[sid] = (m1[0], m2[0])
        else:
            if m2:
                problems.append(
                    f"{sid}: found an R2 file but seq_type is "
                    f"'{inputs['seq_type']}'")
                continue
            single = m1 + unmated
            if len(single) != 1:
                problems.append(
                    f"{sid}: expected exactly one FASTQ, found {len(single)}")
                continue
            assignments[sid] = (single[0], "")

    # One file must not serve two samples; that silently duplicates data.
    claimed: dict = {}
    for sid, mates in assignments.items():
        for f in mates:
            if not f:
                continue
            if f in claimed:
                problems.append(f"{f} matched both '{claimed[f]}' and '{sid}'")
            claimed[f] = sid

    if problems:
        sys.exit("Sample sheet could not be built:\n  - " +
                 "\n  - ".join(problems))

    mdf = mdf.rename(columns={sid_col: "sample_id"})
    mdf["fastq_url"] = mdf["sample_id"].map(lambda s: assignments[s][0])
    mdf["fastq_url_r2"] = mdf["sample_id"].map(lambda s: assignments[s][1])

    lead = ["sample_id", "fastq_url", "fastq_url_r2"]
    mdf = mdf[lead + [c for c in mdf.columns if c not in lead]]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    mdf.to_csv(out_path, sep="\t", index=False)
    log("SAMPLES", f"wrote {out_path} ({len(mdf)} samples, "
                   f"{'paired' if paired else 'single'}-end)")


# ========================================== design-driven contrast generation
def formula_vars(formula: str) -> list:
    """Variable names on the right-hand side of an R formula.

    Handles the shapes this wizard emits ("~ a", "~ a * b", "~ a + b:c").
    """
    body = formula.split("~", 1)[-1]
    return [t for t in re.split(r"[+*:()\s~|]+", body) if t and not t.isdigit()]


def build_contrast_specs(inputs: dict, sheet_path: Path) -> list:
    """Generate declarative contrasts from the actual design.

    Replaces the template's contrasts outright so a new project can never
    inherit variables from a different study. Every model variable is checked
    against the generated sample sheet, and anything missing stops setup here
    rather than failing inside the DE job.
    """
    import pandas as pd

    sheet = pd.read_csv(sheet_path, sep="\t", comment="#")
    cols = set(sheet.columns)
    primary = inputs["primary_factor"]
    fixed = inputs["model_fixed_effects"]

    missing = [v for v in formula_vars(fixed) if v not in cols]
    if missing:
        sys.exit(f"model_fixed_effects references column(s) absent from "
                 f"{sheet_path}: {', '.join(missing)} "
                 f"(available: {', '.join(sorted(cols))})")
    rand = inputs.get("model_random_effects")
    if rand:
        for v in re.findall(r"\|\s*([A-Za-z_.][A-Za-z0-9_.]*)", str(rand)):
            if v not in cols:
                sys.exit(f"model_random_effects references '{v}', which is not "
                         f"a column in {sheet_path}")
    if primary not in cols:
        sys.exit(f"primary_factor '{primary}' is not a column in {sheet_path}")
    if sheet[primary].nunique() < 2:
        sys.exit(f"primary_factor '{primary}' has fewer than 2 levels; "
                 f"there is nothing to compare")

    specs = [{"id": f"{primary}_pairwise", "type": "pairwise",
              "factor": primary, "by": None, "reverse": True}]

    # Interaction design: also compare the primary factor within each level of
    # the factor it interacts with.
    if "*" in fixed or ":" in fixed:
        others = [v for v in formula_vars(fixed)
                  if v != primary and v in cols and sheet[v].nunique() >= 2]
        if others:
            by = others[0]
            specs.append({"id": f"{primary}_within_{by}", "type": "pairwise",
                          "factor": primary, "by": by, "reverse": True})
    return specs


# ==================================================== final config.yaml render
def render_config(inputs: dict, org: dict, ref: dict, out_path: Path):
    cfg = load_yaml(TEMPLATE_CFG)

    cfg["project"]["name"] = inputs["project_name"]
    cfg["project"]["output_dir"] = f"results/{inputs['project_name']}/"

    cfg["organism"].update({
        "scientific_name": f"{org['genus']} {org['species']}",
        "tax_id": org["tax_id"],
        "kegg_code": org["kegg_code"],
        "orgdb_package": org["orgdb_package"],
    })

    cfg["reference"].update({
        # accession/assembly_name are what make the OrgDb cache key specific to
        # this assembly; without them every organism collides on "unknown".
        "accession": ref.get("accession"),
        "assembly_name": ref.get("assembly_name"),
        "genome_fasta_url": ref["genome_fasta_url"],
        "transcriptome_fasta_url": ref["transcriptome_fasta_url"],
        "gtf_url": ref["gtf_url"],
    })
    cfg["reference"]["annotation_tsv"]["path"] = "reference/annotation_info.tsv"
    cfg["reference"].pop("ncbi_to_ensembl", None)

    cfg["samples"]["seq_type"] = inputs["seq_type"]

    cfg["model"]["fixed_effects"] = inputs["model_fixed_effects"]
    cfg["model"]["random_effects"] = inputs["model_random_effects"]
    cfg["model"]["primary_factor"] = inputs["primary_factor"]

    # Replace the template's contrasts outright. Inheriting them is how a
    # treatment-only project ended up testing grape's group/stage variables.
    cfg["contrasts"] = build_contrast_specs(inputs, CONFIG_DIR / "samples.tsv")

    cfg["hpc"]["storage_budget_gb"] = inputs["storage_budget_gb"]

    dump_yaml(cfg, out_path)
    log("CONFIG", f"wrote {out_path}")


# =========================================================================== main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inputs_file", nargs="?", default=None)
    ap.add_argument("--interactive", action="store_true")
    args = ap.parse_args()

    inputs = gather_inputs(args)
    org = resolve_organism(int(inputs["tax_id"]))
    ref = resolve_reference_urls(int(inputs["tax_id"]))

    REFERENCE_DIR.mkdir(parents=True, exist_ok=True)
    fetch_annotation_info(ref.get("gene_info_url"), REFERENCE_DIR / "annotation_info.tsv")

    build_samples_tsv(inputs, CONFIG_DIR / "samples.tsv")
    render_config(inputs, org, ref, CONFIG_DIR / "config.yaml")
    orgdb_status = setup_orgdb(org)

    # -------------------------------------------------- final summary
    print("\n" + "=" * 60)
    print(" Setup complete.")
    print("=" * 60)
    print(f"  config/config.yaml          ✓")
    print(f"  config/samples.tsv          ✓")
    print(f"  reference/annotation_info.tsv   "
          f"{'✓' if (REFERENCE_DIR/'annotation_info.tsv').exists() else '⚠ skipped'}")
    if orgdb_status == "bioconductor":
        print(f"  OrgDb:   {org['orgdb_package']} (Bioconductor)")
    else:
        print(f"  OrgDb:   {org['orgdb_package']} — will be built by the "
              f"build_orgdb rule on first run")

    # Stop before anything can be submitted if the reference never resolved.
    # The config is still written so it can be corrected by hand.
    problems = unresolved_reference_urls(ref)
    if problems:
        print("\n" + "=" * 60)
        print(" SETUP INCOMPLETE — reference URLs did not resolve")
        print("=" * 60)
        for p in problems:
            print(f"  - {p}")
        print("\nFix reference.*_url in config/config.yaml, then run ./submit.sh.")
        sys.exit(1)

    print("\nNext:")
    print("  1. Review config/config.yaml and config/samples.tsv")
    print("  2. Launch:  ./submit.sh")


if __name__ == "__main__":
    main()

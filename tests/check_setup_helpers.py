#!/usr/bin/env python3
"""Self-check for the pandas-free helpers in scripts/setup.py.

The things that must not break:
  - a wildcard or TODO reference URL must never be reported as resolved,
    because curl cannot expand '*' and the fetch job would fail
  - formula variable extraction, since contrast generation validates against it

Run:  python3 tests/check_setup_helpers.py
"""
import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
spec = importlib.util.spec_from_file_location("setup_mod", ROOT / "scripts" / "setup.py")
setup = importlib.util.module_from_spec(spec)
spec.loader.exec_module(setup)

# --- formula_vars ------------------------------------------------------
assert setup.formula_vars("~ treatment") == ["treatment"]
assert setup.formula_vars("~ group * stage") == ["group", "stage"]
assert setup.formula_vars("~ batch + group:stage") == ["batch", "group", "stage"]

# --- unresolved_reference_urls: the wildcard guard ---------------------
good = {
    "genome_fasta_url": "https://x/GCF_1.1_ASM1/GCF_1.1_ASM1_genomic.fna.gz",
    "transcriptome_fasta_url": "https://x/GCF_1.1_ASM1/GCF_1.1_ASM1_rna.fna.gz",
    "gtf_url": "https://x/GCF_1.1_ASM1/GCF_1.1_ASM1_genomic.gtf.gz",
}
assert setup.unresolved_reference_urls(good) == [], "concrete URLs must pass"

wild = dict(good, gtf_url="https://x/GCF_1.1_*/GCF_1.1_*_genomic.gtf.gz")
problems = setup.unresolved_reference_urls(wild)
assert len(problems) == 1 and "wildcard" in problems[0], problems

todo = dict(good, transcriptome_fasta_url="TODO")
assert setup.unresolved_reference_urls(todo), "TODO must be reported unresolved"

missing = dict(good)
del missing["genome_fasta_url"]
assert setup.unresolved_reference_urls(missing), "absent URL must be reported"

# --- mate_of: which read of a pair is this file? -----------------------
assert setup.mate_of("S01_R1_001.fastq.gz") == 1
assert setup.mate_of("S01_R2_001.fastq.gz") == 2
assert setup.mate_of("S01_1.fastq.gz") == 1
assert setup.mate_of("S01_2.fastq.gz") == 2
assert setup.mate_of("S01.fastq.gz") is None, "single-end file has no mate"
# a real delivery filename
assert setup.mate_of("60207403414-0ha-NCGM-4408_L001_R1_001.fastq.gz") == 1
# an embedded lane number must not be mistaken for the mate number
assert setup.mate_of("Lane_2_S3_R1_001.fastq.gz") == 1

# --- files_for_sample: S1 must not swallow S10 -------------------------
pool = [
    "/d/S1_R1.fastq.gz", "/d/S1_R2.fastq.gz",
    "/d/S10_R1.fastq.gz", "/d/S10_R2.fastq.gz",
]
s1 = setup.files_for_sample("S1", pool)
assert sorted(s1) == ["/d/S1_R1.fastq.gz", "/d/S1_R2.fastq.gz"], s1
s10 = setup.files_for_sample("S10", pool)
assert sorted(s10) == ["/d/S10_R1.fastq.gz", "/d/S10_R2.fastq.gz"], s10
assert setup.files_for_sample("S2", pool) == [], "no match must be empty, not a guess"

# --- preset_for_taxid: canonical OrgDb names, not derived ones ---------
# The old code built "org.<Initial><species>.eg.db", which gives
# org.Hsapiens.eg.db for human -- a package that does not exist.
assert setup.preset_for_taxid(9606)["orgdb_package"] == "org.Hs.eg.db"
assert setup.preset_for_taxid(10090)["orgdb_package"] == "org.Mm.eg.db"
assert setup.preset_for_taxid(3702)["orgdb_package"] == "org.At.tair.db"
assert setup.preset_for_taxid(29760)["orgdb_package"] == "org.Vvinifera.eg.db"
# an unknown organism returns {} so the caller guesses *and warns*
assert setup.preset_for_taxid(99999999) == {}

# --- the removed cross-species BioMart fetcher must stay removed -------
assert not hasattr(setup, "fetch_ncbi_to_ensembl"), (
    "fetch_ncbi_to_ensembl queried the human/Arabidopsis marts for every "
    "organism; it must not come back without being made species-aware")

print("check_setup_helpers.py: all assertions passed")

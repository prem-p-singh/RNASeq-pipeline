#!/usr/bin/env python3
"""Self-check for scripts/artifacts.py (P0 / R09, R11; acceptance V18).

Master plan 11: "A done flag alone never proves a result exists. Manifests must
be checked on reuse, or outputs explicitly declared in the DAG; merely listing a
missing file in JSON is insufficient." Its change table requires a deleted DE
table with an old flag present to be detected and the affected result rebuilt.

The defect: Snakemake targets de_done.flag, enrichment_done.flag and
wgcna_done.flag, not the result tables. Deleting a DE table left the flag in
place, so the target looked satisfied and the table was never rebuilt.

What must not break:
  - an intact stage is left alone
  - a deleted or emptied result invalidates its stage, clearing flag and
    manifest so Snakemake reschedules it
  - a flag with no manifest cannot stand on its own
  - only the affected stage is invalidated
  - statuses that do not promise a file are not treated as missing artifacts
  - NA and blank paths are not mistaken for missing files
  - the Snakefile performs this before the DAG is built

Run:  python3 tests/check_artifacts.py
"""
import importlib.util
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
spec = importlib.util.spec_from_file_location(
    "artifacts", ROOT / "scripts" / "artifacts.py")
art = importlib.util.module_from_spec(spec)
spec.loader.exec_module(art)


def project(*, de_tables=("c1", "c2"), enrich=True):
    """An output directory that looks like a completed run."""
    d = Path(tempfile.mkdtemp())
    (d / "DE_Results").mkdir()
    rows = ["schema_version\tcontrast_name\tstatus\tanalysis_table\tmeaningful_table"]
    for c in de_tables:
        (d / "DE_Results" / f"{c}_DE_analysis.tsv").write_text("gene_id\tlogFC\ng1\t1.0\n")
        (d / "DE_Results" / f"{c}_DE_meaningful.tsv").write_text("gene_id\tlogFC\ng1\t1.0\n")
        rows.append(f"1\t{c}\tsucceeded\tDE_Results/{c}_DE_analysis.tsv\t"
                    f"DE_Results/{c}_DE_meaningful.tsv")
    (d / "de_manifest.tsv").write_text("\n".join(rows) + "\n")
    (d / "de_done.flag").touch()

    if enrich:
        (d / "Enrichment" / "GO").mkdir(parents=True)
        (d / "Enrichment" / "GO" / "c1_GO.tsv").write_text("ID\tDescription\nGO:1\tx\n")
        (d / "enrichment_status.tsv").write_text(
            "contrast\tmethod\tstatus\treason\tn_terms\toutput_table\tschema_version\n"
            "c1\tGO\tsucceeded\t\t1\tEnrichment/GO/c1_GO.tsv\t1\n"
            "c1\tKEGG\tsucceeded_empty\t\t0\tNA\t1\n"
            "c2\tGO\tskipped_by_policy\trun_go is false\tNA\tNA\t1\n"
            "c2\tKEGG\tunavailable\tno kegg_code\tNA\tNA\t1\n")
        (d / "enrichment_done.flag").touch()
    return d


# --- 1. an intact run is untouched ------------------------------------
d = project()
assert art.invalidate_stale(d) == [], "intact run was invalidated"
assert (d / "de_done.flag").is_file() and (d / "enrichment_done.flag").is_file()

# --- 2. the headline scenario: deleted table, flag still present ------
d = project()
(d / "DE_Results" / "c2_DE_analysis.tsv").unlink()
rep = art.invalidate_stale(d)
stages = {r["stage"] for r in rep}
assert stages == {"differential_expression"}, stages
assert not (d / "de_done.flag").exists(), "flag outlived a deleted result"
assert not (d / "de_manifest.tsv").exists()
assert "DE_Results/c2_DE_analysis.tsv" in rep[0]["missing"], rep[0]
# enrichment was not touched: only the affected stage is invalidated
assert (d / "enrichment_done.flag").is_file()

# --- 3. a zero-byte table is not a result -----------------------------
d = project()
(d / "DE_Results" / "c1_DE_analysis.tsv").write_text("")
rep = art.invalidate_stale(d)
assert rep and rep[0]["stage"] == "differential_expression", rep
assert not (d / "de_done.flag").exists()

# --- 4. a flag with no manifest cannot stand alone --------------------
d = project()
(d / "de_manifest.tsv").unlink()
rep = art.invalidate_stale(d)
assert any(r["reason"] == "flag without manifest" for r in rep), rep
assert not (d / "de_done.flag").exists()

# --- 5. enrichment: only statuses that promise a file are checked -----
d = project()
(d / "Enrichment" / "GO" / "c1_GO.tsv").unlink()
rep = art.invalidate_stale(d)
assert {r["stage"] for r in rep} == {"enrichment"}, rep
assert not (d / "enrichment_done.flag").exists()
assert (d / "de_done.flag").is_file(), "DE was invalidated by an enrichment gap"

# succeeded_empty, skipped_by_policy and unavailable rows carry NA paths and
# must not be read as missing artifacts
d = project()
rep = art.invalidate_stale(d)
assert rep == [], "a non-result status was treated as a missing artifact"

# --- 6. both stages can be invalidated independently ------------------
d = project()
(d / "DE_Results" / "c1_DE_analysis.tsv").unlink()
(d / "Enrichment" / "GO" / "c1_GO.tsv").unlink()
rep = art.invalidate_stale(d)
assert {r["stage"] for r in rep} == {"differential_expression", "enrichment"}, rep
assert not (d / "de_done.flag").exists()
assert not (d / "enrichment_done.flag").exists()

# --- 7. a run that never happened is not an error ---------------------
empty = Path(tempfile.mkdtemp())
assert art.invalidate_stale(empty) == []

# --- 8. the Snakefile runs this before the DAG is built ---------------
snake = (ROOT / "Snakefile").read_text()
assert "invalidate_stale" in snake, "Snakefile does not invalidate stale results"
pos_inv = snake.index("invalidate_stale")
pos_include = snake.index('include: "workflow/rules/')
pos_all = snake.index("rule all:")
assert pos_inv < pos_include < pos_all, (
    "invalidation must run before the rule includes and rule all, or the DAG "
    "is resolved against a stale flag")

# --- 9. the emitted statuses are the master plan 11 taxonomy ----------
enr = (ROOT / "workflow" / "scripts" / "04_enrichment.R").read_text()
allowed = {"succeeded", "succeeded_empty", "skipped_by_policy", "unavailable",
           "failed", "planned", "running", "blocked"}
import re
emitted = set(re.findall(r'record\([^)]*?"(succeeded|succeeded_empty|failed|'
                         r'skipped_by_policy|unavailable|empty|skipped)"', enr, re.S))
stale_names = emitted - allowed
assert not stale_names, f"pre-taxonomy status still emitted: {stale_names}"

print("check_artifacts.py: all assertions passed")

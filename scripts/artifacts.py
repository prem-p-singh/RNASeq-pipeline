#!/usr/bin/env python3
"""artifacts.py — a completion flag never proves a result exists (R09, R11).

Master plan section 11: "A done flag alone never proves a result exists.
Manifests must be checked on reuse, or outputs explicitly declared in the DAG;
merely listing a missing file in JSON is insufficient." Its change table
requires that a deleted DE table with an old flag present is detected and the
affected result rebuilt (acceptance V18).

Why manifests rather than declared outputs. Result table names embed the
resolved contrast, for example `group_within_stage__treated_vs_control__stage.veraison`.
The levels come from the sample sheet through emmeans, so the names are not
known when Snakemake builds the DAG. Declaring them would mean reimplementing
level resolution in the Snakefile, a second copy of logic that must agree with
03_de.R. R09's closure permits either route; this takes the manifest route.

How a missing artifact causes a rebuild rather than a warning. Snakemake reruns
a rule when one of its declared outputs is absent, and it decides that during
DAG construction. So this module runs BEFORE the DAG is built, at Snakefile
parse time: it reads each stage manifest, checks every artifact it lists, and
when one is missing it removes that stage's flag and manifest. Snakemake then
resolves the DAG against a genuinely incomplete stage and reschedules it,
together with anything downstream, in the same invocation.

Pure filesystem work, no Snakemake import, so it is testable on its own.
"""
from __future__ import annotations

import csv
import os
from pathlib import Path

# Per stage: the manifest, the flag that must not outlive it, and which
# manifest columns hold paths to artifacts that must exist.
# Paths in a manifest are relative to the output directory that holds it.
STAGE_MANIFESTS = {
    "differential_expression": {
        "manifest": "de_manifest.tsv",
        "flag": "de_done.flag",
        "path_columns": ("analysis_table", "meaningful_table"),
    },
    "enrichment": {
        "manifest": "enrichment_status.tsv",
        "flag": "enrichment_done.flag",
        "path_columns": ("output_table",),
        # Only rows that claim a result are checked; see status_is_result.
        "status_column": "status",
    },
    "wgcna": {
        "manifest": "wgcna_manifest.tsv",
        "flag": "wgcna_done.flag",
        "path_columns": ("output_table",),
        "status_column": "status",
    },
}

# Statuses that assert a file was produced. Master plan section 11 taxonomy:
# planned, running, succeeded, succeeded_empty, skipped_by_policy, unavailable,
# blocked, failed. Only `succeeded` promises an artifact on disk; an empty
# result is a valid outcome that need not have written a table.
RESULT_STATUSES = frozenset({"succeeded"})


def read_rows(path: Path) -> list[dict]:
    """Rows of a TSV manifest, or [] when it is absent or has no data rows."""
    if not path.is_file():
        return []
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def missing_artifacts(manifest_path: Path, out_dir: Path, spec: dict) -> list[str]:
    """Paths a manifest claims exist but which are absent or empty.

    An empty file counts as missing: a zero-byte result table is not a result,
    and it is the shape a crashed or interrupted write leaves behind.
    """
    status_col = spec.get("status_column")
    missing = []
    rows = read_rows(manifest_path)
    required = set(spec["path_columns"])
    if status_col:
        required.add(status_col)
    if not rows or not required.issubset(rows[0]):
        return [f"invalid or empty manifest: {manifest_path.name}"]
    for row in rows:
        if status_col and row.get(status_col, "") not in RESULT_STATUSES:
            if row.get(status_col) not in {"succeeded_empty", "skipped_by_policy", "unavailable"}:
                missing.append(f"nonterminal or invalid completion status: {row.get(status_col)!r}")
            continue
        for col in spec["path_columns"]:
            rel = (row.get(col) or "").strip()
            if not rel or rel.upper() in {"NA", "NONE", ""}:
                if col != "meaningful_table":
                    missing.append(f"missing required path: {col}")
                continue
            target = out_dir / rel
            if not target.is_file() or target.stat().st_size == 0:
                missing.append(rel)
    return missing


def invalidate_stale(out_dir) -> list[dict]:
    """Drop flags whose manifests point at artifacts that are gone.

    Returns one record per invalidated stage. Removing the flag and manifest is
    what makes Snakemake rebuild: with them absent, the stage's declared outputs
    are incomplete and it is rescheduled. Nothing else is deleted, and a stage
    whose artifacts are all present is left untouched.
    """
    out_dir = Path(out_dir)
    report = []
    for stage, spec in STAGE_MANIFESTS.items():
        manifest = out_dir / spec["manifest"]
        flag = out_dir / spec["flag"]
        if not manifest.is_file():
            # No manifest to contradict. A flag with no manifest is itself
            # unprovable, so it is cleared too.
            if flag.is_file():
                flag.unlink()
                report.append({"stage": stage, "reason": "flag without manifest",
                               "missing": [], "removed": [spec["flag"]]})
            continue

        missing = missing_artifacts(manifest, out_dir, spec)
        if not missing:
            continue

        removed = []
        for p in (flag, manifest):
            if p.is_file():
                p.unlink()
                removed.append(p.name)
        report.append({"stage": stage,
                       "reason": f"{len(missing)} artifact(s) listed but absent",
                       "missing": missing, "removed": removed})
    return report


def format_report(report: list[dict]) -> str:
    lines = []
    for r in report:
        lines.append(f"  {r['stage']}: {r['reason']}")
        for m in r["missing"][:10]:
            lines.append(f"      missing: {m}")
        if len(r["missing"]) > 10:
            lines.append(f"      ... and {len(r['missing']) - 10} more")
        lines.append(f"      cleared: {', '.join(r['removed'])} -> stage will rerun")
    return "\n".join(lines)


def demo():
    """Self-check: a deleted table with its flag intact must invalidate."""
    import tempfile
    d = Path(tempfile.mkdtemp())
    (d / "DE_Results").mkdir()
    (d / "DE_Results" / "c1_DE_analysis.tsv").write_text("gene\tlogFC\ng1\t1\n")
    (d / "de_manifest.tsv").write_text(
        "contrast_name\tanalysis_table\tmeaningful_table\n"
        "c1\tDE_Results/c1_DE_analysis.tsv\tNA\n")
    (d / "de_done.flag").touch()

    assert invalidate_stale(d) == [], "intact stage must not be invalidated"
    assert (d / "de_done.flag").is_file()

    (d / "DE_Results" / "c1_DE_analysis.tsv").unlink()
    rep = invalidate_stale(d)
    assert len(rep) == 1 and rep[0]["stage"] == "differential_expression", rep
    assert not (d / "de_done.flag").exists(), "flag survived a missing artifact"
    assert not (d / "de_manifest.tsv").exists()
    print("artifacts.demo: ok")


if __name__ == "__main__":
    demo()

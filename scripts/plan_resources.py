#!/usr/bin/env python3
"""plan_resources.py — build and enforce the storage plan for one run.

Wraps scripts/resources.py for the launcher. Exits non-zero when the run cannot
fit, after first trying a lower concurrency (RS08), and writes
gates/resource_plan.json, which master plan 12.1 lists as a required artifact.

Replaces the old fixed-budget logic: one platform-wide 20 GB cap, per-assay
FASTQ sizes guessed from a label, and a concurrency limit derived from sample
count. Master plan 10.1: "There is no default 20 GB limit and no sample-count
cutoff that rejects an analysis."

RS12 and master plan 10.7 require a legacy storage_budget_gb to produce a
migration warning rather than being silently erased or silently obeyed: only
the user knows whether that number was a real site quota or the old template
default, so it is surfaced and, when kept, applied as an explicit quota.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import metadata as md          # noqa: E402
import resources as rs         # noqa: E402

REPO = Path(__file__).resolve().parent.parent
GB = 1024 ** 3

# seq_type is still the only assay signal a legacy project carries (U01 splits
# it properly in P1); map it onto the resource model's assay families.
ASSAY_OF_SEQ_TYPE = {
    "tagseq": "bulk_end_tag",
    "rnaseq_single": "bulk",
    "rnaseq_paired": "bulk",
}


def input_uris(proj: Path, cfg: dict) -> list[str]:
    """Every read file this cohort refers to, each once.

    Prefers the three-table form, where a library may have several lanes, and
    falls back to the single sheet's fastq_url columns.
    """
    tables = md.load_tables(proj / "metadata")
    if tables.get("reads"):
        return [r.get("uri", "") for r in tables["reads"] if r.get("uri")]

    sheet = proj / (cfg.get("samples", {}) or {}).get("sheet", "config/samples.tsv")
    uris = []
    for row in md.read_tsv(sheet):
        for col in ("fastq_url", "fastq_url_r2"):
            if row.get(col):
                uris.append(row[col])
    return uris


def n_libraries(proj: Path, cfg: dict) -> int:
    tables = md.load_tables(proj / "metadata")
    if tables.get("libraries"):
        return len(tables["libraries"])
    sheet = proj / (cfg.get("samples", {}) or {}).get("sheet", "config/samples.tsv")
    return len(md.read_tsv(sheet))


def legacy_quota(cfg: dict) -> tuple[int | None, list[str]]:
    """Interpret a legacy hpc.storage_budget_gb (RS12, master plan 10.7)."""
    hpc = cfg.get("hpc") or {}
    if "storage_budget_gb" not in hpc or hpc.get("storage_budget_gb") in (None, ""):
        return None, []
    value = hpc["storage_budget_gb"]
    notes = [
        f"hpc.storage_budget_gb={value} is a legacy setting. The platform no "
        f"longer has a default storage ceiling (master plan 10.1), and this "
        f"value is NOT a planning input.",
        "It is applied as an explicit quota on the results filesystem, because "
        "a real site limit must not be silently erased (RS12). If it was only "
        "the old template default of 20, remove the key and let the plan use "
        "measured free space.",
    ]
    try:
        return int(float(value) * GB), notes
    except (TypeError, ValueError):
        return None, notes + [f"  value {value!r} is not a number; ignored"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-d", "--project-dir", required=True)
    ap.add_argument("-c", "--configfile", default="config/config.yaml")
    ap.add_argument("--profile", default="")
    ap.add_argument("--out", required=True)
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args()

    proj = Path(a.project_dir).resolve()
    cfg = yaml.safe_load((proj / a.configfile).read_text()
                         if not Path(a.configfile).is_absolute()
                         else Path(a.configfile).read_text())

    seq_type = (cfg.get("samples", {}) or {}).get("seq_type", "")
    assay = ASSAY_OF_SEQ_TYPE.get(seq_type, "bulk")
    libs = n_libraries(proj, cfg)
    uris = input_uris(proj, cfg)
    inputs = rs.measure_inputs(uris)

    # Concurrency the scheduler would allow, before storage has its say.
    requested = 1
    if a.profile:
        try:
            requested = int(yaml.safe_load(
                (Path(a.profile) / "config.yaml").read_text()).get("jobs", 1))
        except (OSError, ValueError, AttributeError):
            requested = 1

    out_dir = proj / (cfg.get("project", {}) or {}).get("output_dir", "results")
    ref_cache = (cfg.get("reference", {}) or {}).get("cache_dir")
    paths = {
        "results": out_dir,
        "scratch": proj / "tmp",
        "cache": Path(ref_cache).expanduser() if ref_cache else out_dir,
    }
    local_inputs = [u for u in uris if "://" not in u]
    if local_inputs:
        paths["inputs"] = Path(local_inputs[0]).parent

    quota_bytes, quota_notes = legacy_quota(cfg)
    quotas = {"results": quota_bytes} if quota_bytes is not None else {}

    kwargs = dict(repo_root=REPO, inputs=inputs, n_libraries=libs, assay=assay,
                  concurrency=requested, paths=paths, quotas=quotas)
    report = rs.plan_storage(**kwargs)
    report = rs.reduce_concurrency(report, **kwargs)
    report["legacy_notes"] = quota_notes
    report["seq_type"] = seq_type

    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(json.dumps(report, indent=2) + "\n")

    if not a.quiet:
        print(rs.format_report(report))
        for n in quota_notes:
            print(f"  NOTE: {n}")

    # RS09: if even one job or the final retained data cannot fit, stop and
    # report the requirement rather than quietly choosing another route.
    bad = [m for m in report["mounts"] if m["verdict"] != "ok"]
    if bad:
        print("\nStorage plan cannot be satisfied:", file=sys.stderr)
        for m in bad:
            if m["verdict"] == "insufficient":
                print(f"  {', '.join(m['paths'])}: need "
                      f"{m['required_free_gb']:g} GB, have "
                      f"{m['available_gb']:g} GB, short by "
                      f"{m['shortfall_gb']:g} GB", file=sys.stderr)
            else:
                print(f"  {', '.join(m['paths'])}: {m['verdict']}", file=sys.stderr)
        print("  Concurrency was already reduced as far as it goes "
              f"({report['planned_concurrency']}). Free space or choose another "
              f"location; the analysis route is unchanged (RS08, RS09).",
              file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()

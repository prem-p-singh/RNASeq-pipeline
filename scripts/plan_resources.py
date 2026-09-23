#!/usr/bin/env python3
"""plan_resources.py — report advisory storage estimates for one run.

Wraps scripts/resources.py for the launcher. Capacity concerns only warn;
they do not change concurrency or reject execution. Writes
gates/resource_plan.json, which master plan 12.1 lists as a required artifact.

Replaces the old fixed-budget logic: one platform-wide 20 GB cap, per-assay
FASTQ sizes guessed from a label, and a concurrency limit derived from sample
count. Master plan 10.1: "There is no default 20 GB limit and no sample-count
cutoff that rejects an analysis."

RS12 and master plan 10.7 require a legacy storage_budget_gb to produce a
migration warning rather than being silently erased or silently obeyed: only
the user knows whether that number was a real site quota or the old template
default, so it is surfaced as advisory quota information.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path
from urllib.parse import unquote, urlsplit

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config_resolve as cr
import metadata as md          # noqa: E402
import reference_cache as rc
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
    canonical = md.execution_inputs(proj, cfg)
    if canonical is not None:
        return [unit[role]["uri"] for library in canonical.values()
                for unit in library["units"] for role in ("r1", "r2") if unit.get(role)]

    sheet = proj / (cfg.get("samples", {}) or {}).get("sheet", "config/samples.tsv")
    uris = []
    for row in md.read_tsv(sheet):
        for col in ("fastq_url", "fastq_url_r2"):
            if row.get(col):
                uris.append(row[col])
    return uris


def n_libraries(proj: Path, cfg: dict) -> int:
    canonical = md.execution_inputs(proj, cfg)
    if canonical is not None:
        return len(canonical)
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
        "It is reported as advisory quota information on the results filesystem. "
        "It never blocks launch or lowers concurrency. If it was only "
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

    cfg, _, issues = cr.resolve(REPO, [("project", cfg)])
    problems = cr.blocking(issues, REPO)
    if problems:
        raise SystemExit("Invalid planning configuration: " + "; ".join(f"{code}: {detail}" for code, detail in problems))

    seq_type = (cfg.get("samples", {}) or {}).get("seq_type", "")
    assay = ASSAY_OF_SEQ_TYPE.get(seq_type, "bulk")
    libs = n_libraries(proj, cfg)
    uris = input_uris(proj, cfg)
    uris = [str((proj / md.local_path(u)).resolve())
            if md.local_path(u) is not None else u for u in uris]
    inputs = rs.measure_inputs(uris)

    # Concurrency the scheduler would allow, before storage has its say.
    requested = 1
    if a.profile:
        try:
            requested = int(yaml.safe_load(
                (Path(a.profile) / "config.yaml").read_text()).get("jobs", 1))
        except (OSError, ValueError, AttributeError):
            requested = 1

    configured_limit = (cfg.get("hpc") or {}).get("samples_in_flight")
    if configured_limit is not None:
        if isinstance(configured_limit, bool) or int(configured_limit) < 1:
            raise SystemExit("hpc.samples_in_flight must be a positive integer")
        requested = min(requested, int(configured_limit))

    out_dir = proj / (cfg.get("project", {}) or {}).get("output_dir", "results")
    ref_cache = (cfg.get("reference", {}) or {}).get("cache_dir")
    paths = {
        "results": out_dir,
        "scratch": out_dir,  # Stage 1 stores downloads and trims in quant/<sample>
        "cache": Path(ref_cache).expanduser() if ref_cache else out_dir,
    }
    local_inputs = [u for u in uris if "://" not in u]
    if local_inputs:
        paths["inputs"] = Path(local_inputs[0]).parent

    quota_bytes, quota_notes = legacy_quota(cfg)
    quotas = {"results": quota_bytes} if quota_bytes is not None else {}

    ref_cfg = cfg.get("reference") or {}
    entry = (Path(ref_cache).expanduser() / rc.cache_key(ref_cfg, ref_cfg.get("kmer", 31), ref_cfg.get("decoys"))
             if ref_cache else proj / "reference")
    complete_reference = not rc.entry_problems(entry)
    environment = os.environ.get("RNASEQ_ENV_PREFIX") or os.environ.get("CONDA_PREFIX")
    if environment:
        paths["environment"] = Path(environment)
    cache_present = {"reference": complete_reference, "index_build": complete_reference,
                     "environment": bool(environment and (Path(environment) / "conda-meta").is_dir())}
    reference_sources = []
    for field in ("transcriptome_fasta_url", "gtf_url"):
        uri = ref_cfg.get(field) or ""
        if uri.startswith("file://"):
            source = Path(unquote(urlsplit(uri).path))
            if source.is_file() and source.suffix not in {".gz", ".bz2", ".xz"}:
                reference_sources.append(source.stat().st_size)
    reference_source_bytes = sum(reference_sources) if len(reference_sources) == 2 else None
    kwargs = dict(repo_root=REPO, inputs=inputs, n_libraries=libs, assay=assay,
                  concurrency=requested, paths=paths, quotas=quotas, cache_present=cache_present,
                  reference_source_bytes=reference_source_bytes,
                  retain_downloads=not (cfg.get("hpc") or {}).get("delete_fastq_after_quant", True))
    canonical = md.execution_inputs(proj, cfg)
    if canonical is not None:
        kwargs["library_uris"] = [[unit[role]["uri"] for unit in library["units"]
                                   for role in ("r1", "r2") if unit.get(role)]
                                  for library in canonical.values()]
    report = rs.plan_storage(**kwargs)
    report["legacy_notes"] = quota_notes
    report["seq_type"] = seq_type

    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(json.dumps(report, indent=2) + "\n")

    if not a.quiet:
        print(rs.format_report(report))
        for n in quota_notes:
            print(f"  NOTE: {n}")

    # Quiet mode suppresses the estimate, not actionable capacity cautions.
    if a.quiet:
        for warning in report["warnings"]:
            print(f"CAUTION: {warning}", file=sys.stderr)


if __name__ == "__main__":
    main()

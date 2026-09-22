#!/usr/bin/env python3
"""resources.py — per-filesystem storage planning (R15, U07).

Master plan 10.1 removes the fixed ceiling outright: "There is no default 20 GB
limit and no sample-count cutoff that rejects an analysis." 10.7 requires
hpc.storage_budget_gb and the fixed per-assay FASTQ assumptions to be replaced
by measured planning, and any real site quota to be preserved as an explicit
constraint rather than an implicit default.

What was there before: one number for the whole platform, per-assay FASTQ sizes
guessed from a label, and a concurrency cap derived from sample count. None of
it looked at the data or at the filesystem the run would land on.

The accounting follows 10.3, per filesystem m:

    existing_m        bytes already present; reported, not charged again
    retained_new_m    new inputs and outputs still there after the run
    cache_new_m       reference, index and environment missing from m
    temporary_m       live staging and scratch, at peak concurrency
    peak_new_m        conservative upper bound, as 10.3 permits:
                      retained_new + cache_new + max_concurrent_temporary
    required_free_m   peak_new + uncertainty margin + operational reserve

10.3 calls that sum "clearly labeled as an upper bound because these maxima
need not coincide", which is what `basis` records.

Paths are grouped by device id, so two directories on one filesystem do not
each claim its free space (RS07). Available space is the lower of physical free
bytes and a declared quota (10.3); a quota that cannot be queried is recorded
unknown and never reported as unlimited.
"""
from __future__ import annotations

import os
import shutil
from dataclasses import dataclass, field
from pathlib import Path

import yaml

GB = 1024 ** 3


def _band(entry: dict) -> tuple[float, float, str]:
    return float(entry["low"]), float(entry["high"]), entry.get("source", "assumed")


def load_models(repo_root) -> dict:
    return yaml.safe_load((Path(repo_root) / "config" / "resource_models.yaml").read_text())


# ------------------------------------------------------------- measurement
def measure_inputs(uris) -> dict:
    """Actual bytes for every unique input, or an explicit unknown.

    RS01: "Sum actual unique input objects, both mates and every lane." The
    dict is keyed by URI, so one file referenced twice is counted once.
    Remote objects are not fetched to size them; they are recorded unknown and
    the caller falls back to a labelled band (RS02).
    """
    out = {}
    for uri in dict.fromkeys(uris):          # unique, order preserved
        if not uri:
            continue
        if "://" in uri and not uri.startswith("file://"):
            out[uri] = {"bytes": None, "source": "unknown_remote"}
            continue
        p = Path(uri[len("file://"):] if uri.startswith("file://") else uri)
        try:
            out[uri] = {"bytes": p.stat().st_size, "source": "measured"}
        except OSError:
            out[uri] = {"bytes": None, "source": "unknown_missing"}
    return out


def filesystem_id(path) -> str:
    """Identify the filesystem a path lives on.

    RS07: paths that share a filesystem must be accounted once. Walks up to the
    nearest existing ancestor, since a planned directory may not exist yet.
    """
    p = Path(path).resolve()
    while not p.exists() and p != p.parent:
        p = p.parent
    try:
        return str(os.stat(p).st_dev)
    except OSError:
        return "unknown"


def free_bytes(path) -> int | None:
    p = Path(path).resolve()
    while not p.exists() and p != p.parent:
        p = p.parent
    try:
        return shutil.disk_usage(p).free
    except OSError:
        return None


# ------------------------------------------------------------------ plan
@dataclass
class MountPlan:
    mount: str
    paths: list = field(default_factory=list)
    retained_new: float = 0.0
    cache_new: float = 0.0
    max_concurrent_temporary: float = 0.0
    existing: float = 0.0
    free: int | None = None
    quota_free: int | None = None
    notes: list = field(default_factory=list)

    @property
    def peak_new(self) -> float:
        return self.retained_new + self.cache_new + self.max_concurrent_temporary

    def required_free(self, margin: float, reserve_gb: float) -> float:
        return self.peak_new * (1.0 + margin) + reserve_gb * GB

    @property
    def available(self) -> int | None:
        """Lower of physical free and remaining quota (master plan 10.3)."""
        vals = [v for v in (self.free, self.quota_free) if v is not None]
        return min(vals) if vals else None


def plan_storage(*, repo_root, inputs: dict, n_libraries: int, assay: str,
                 concurrency: int, paths: dict, quotas: dict | None = None,
                 cache_present: dict | None = None,
                 retain_alignments: bool = False) -> dict:
    """Build a per-filesystem plan.

    `paths` maps a role (inputs, results, scratch, cache) to a directory.
    `quotas` optionally maps a role to remaining bytes, for a site limit that
    cannot be read from the filesystem. `cache_present` marks caches already on
    disk, which are reported but not charged as new allocation (RS05).
    """
    models = load_models(repo_root)
    band = models["assays"].get(assay) or models["assays"]["bulk"]
    planning = models["planning"]
    margin = float(planning["uncertainty_margin_fraction"])
    reserve = float(planning["operational_reserve_gb"])
    cache_present = cache_present or {}
    quotas = quotas or {}

    # --- inputs: measured where possible, banded where not (RS01/RS02) ---
    measured = [v["bytes"] for v in inputs.values() if v["bytes"] is not None]
    unknown = [k for k, v in inputs.items() if v["bytes"] is None]
    in_low, in_high, in_src = _band(band["input_gb"])
    input_bytes = sum(measured) + len(unknown) * in_high * GB
    uncertainty = ("measured" if not unknown else
                   f"{len(measured)} measured, {len(unknown)} estimated at the "
                   f"upper band ({in_high:g} GB, source={in_src})")

    ret_low, ret_high, _ = _band(band["retained_gb"])
    tmp_low, tmp_high, _ = _band(band["temporary_gb"])
    coh_low, coh_high, _ = _band(models["cohort"]["aggregate_outputs_gb"])

    # RS04: scratch scales with CONCURRENT jobs, not with cohort size. RS13:
    # more libraries at the same concurrency grow retained data, not the peak.
    concurrent = max(1, min(concurrency, n_libraries))
    temporary = concurrent * tmp_high * GB
    retained = n_libraries * ret_high * GB + coh_high * GB
    if retain_alignments:
        # RS06: retained alignments are their own per-library cost and are not
        # predicted by sample count alone.
        retained += n_libraries * in_high * GB

    ref_low, ref_high, _ = _band(models["cache"]["reference_bundle_gb"])
    idx_low, idx_high, _ = _band(models["cache"]["index_build_overhead_gb"])
    env_low, env_high, _ = _band(models["cache"]["conda_environment_gb"])

    # --- group the roles by filesystem (RS07) ---------------------------
    mounts: dict[str, MountPlan] = {}

    def mp(role) -> MountPlan | None:
        target = paths.get(role)
        if not target:
            return None
        fid = filesystem_id(target)
        m = mounts.setdefault(fid, MountPlan(mount=fid))
        if str(target) not in m.paths:
            m.paths.append(str(target))
        if m.free is None:
            m.free = free_bytes(target)
        q = quotas.get(role)
        if q is not None:
            m.quota_free = q if m.quota_free is None else min(m.quota_free, q)
        return m

    m_in, m_res = mp("inputs"), mp("results")
    m_scratch, m_cache = mp("scratch"), mp("cache")

    if m_in is not None:
        # Inputs already on disk are existing bytes, not a new allocation.
        m_in.existing += input_bytes
        m_in.notes.append(f"inputs: {uncertainty}")
    if m_res is not None:
        m_res.retained_new += retained
        m_res.notes.append(
            f"{n_libraries} libraries x {ret_high:g} GB retained, plus "
            f"{coh_high:g} GB cohort outputs"
            + (", including retained alignments" if retain_alignments else ""))
    if m_scratch is not None:
        m_scratch.max_concurrent_temporary += temporary
        m_scratch.notes.append(
            f"{concurrent} concurrent job(s) x {tmp_high:g} GB working set")
    if m_cache is not None:
        for label, hi, key in (("reference bundle", ref_high, "reference"),
                               ("index build workspace", idx_high, "index_build"),
                               ("conda environment", env_high, "environment")):
            if cache_present.get(key):
                m_cache.existing += hi * GB
                m_cache.notes.append(f"{label} already present, not charged (RS05)")
            else:
                m_cache.cache_new += hi * GB
                m_cache.notes.append(f"{label}: {hi:g} GB")

    # --- verdicts --------------------------------------------------------
    report = {
        "schema_version": 1,
        "basis": ("conservative upper bound: retained_new + cache_new + "
                  "max_concurrent_temporary, per master plan 10.3; these "
                  "maxima need not coincide"),
        "uncertainty": uncertainty,
        "n_libraries": n_libraries,
        "assay": assay,
        "requested_concurrency": concurrency,
        "planned_concurrency": concurrent,
        "mounts": [],
        "blocking": [],
    }

    for fid, m in sorted(mounts.items()):
        need = m.required_free(margin, reserve)
        avail = m.available
        entry = {
            "filesystem": fid,
            "paths": m.paths,
            "existing_gb": round(m.existing / GB, 2),
            "retained_new_gb": round(m.retained_new / GB, 2),
            "cache_new_gb": round(m.cache_new / GB, 2),
            "max_concurrent_temporary_gb": round(m.max_concurrent_temporary / GB, 2),
            "peak_new_gb": round(m.peak_new / GB, 2),
            "required_free_gb": round(need / GB, 2),
            "available_gb": round(avail / GB, 2) if avail is not None else None,
            "quota_known": m.quota_free is not None,
            "notes": m.notes,
        }
        if avail is None:
            entry["verdict"] = "unknown_capacity"
            report["blocking"].append(
                f"{fid}: free space could not be determined for "
                f"{', '.join(m.paths)}; capacity is recorded unknown, never "
                f"assumed unlimited")
        elif need > avail:
            entry["verdict"] = "insufficient"
            entry["shortfall_gb"] = round((need - avail) / GB, 2)
        else:
            entry["verdict"] = "ok"
        report["mounts"].append(entry)

    return report


def reduce_concurrency(report: dict, **plan_kwargs) -> dict:
    """RS08: reduce concurrency before abandoning the scientific route.

    Only scratch scales with concurrency, so this is retried downward while a
    scratch shortfall is the reason. RS09: if even one job, or the final
    retained data, cannot fit, the caller must block and report the requirement
    rather than silently choosing a different route.
    """
    concurrency = plan_kwargs.get("concurrency", 1)
    while concurrency > 1:
        bad = [m for m in report["mounts"] if m["verdict"] == "insufficient"]
        if not bad:
            return report
        # only worth retrying if scratch is what is over
        if not any(m["max_concurrent_temporary_gb"] > 0 for m in bad):
            return report
        concurrency -= 1
        plan_kwargs["concurrency"] = concurrency
        report = plan_storage(**plan_kwargs)
    return report


def format_report(report: dict) -> str:
    lines = [f"Storage plan: {report['n_libraries']} librar"
             f"{'y' if report['n_libraries'] == 1 else 'ies'}, assay "
             f"{report['assay']}, concurrency {report['planned_concurrency']}"
             + (f" (reduced from {report['requested_concurrency']})"
                if report["planned_concurrency"] != report["requested_concurrency"] else ""),
             f"  inputs: {report['uncertainty']}"]
    for m in report["mounts"]:
        avail = ("unknown" if m["available_gb"] is None
                 else f"{m['available_gb']:g} GB available")
        lines.append(
            f"  [{m['verdict']}] {', '.join(m['paths'])}\n"
            f"      need {m['required_free_gb']:g} GB, {avail}"
            + (f", short by {m['shortfall_gb']:g} GB" if "shortfall_gb" in m else ""))
        for n in m["notes"]:
            lines.append(f"      - {n}")
    lines.append(f"  basis: {report['basis']}")
    return "\n".join(lines)


def demo():
    import tempfile
    d = Path(tempfile.mkdtemp())
    (d / "in").mkdir(); (d / "res").mkdir()
    f1 = d / "in" / "a.fq.gz"; f1.write_bytes(b"x" * 1024)
    inputs = measure_inputs([str(f1), str(f1), "s3://bucket/b.fq.gz"])
    assert len(inputs) == 2, "a repeated URI must be counted once"
    assert inputs[str(f1)]["bytes"] == 1024
    assert inputs["s3://bucket/b.fq.gz"]["source"] == "unknown_remote"

    rep = plan_storage(repo_root=Path(__file__).resolve().parent.parent,
                       inputs=inputs, n_libraries=2, assay="bulk_end_tag",
                       concurrency=2,
                       paths={"inputs": d / "in", "results": d / "res",
                              "scratch": d / "res", "cache": d / "res"})
    assert rep["mounts"], rep
    # more libraries at the same concurrency: retained grows, scratch does not
    rep4 = plan_storage(repo_root=Path(__file__).resolve().parent.parent,
                        inputs=inputs, n_libraries=8, assay="bulk_end_tag",
                        concurrency=2,
                        paths={"results": d / "res", "scratch": d / "res"})
    m1 = rep["mounts"][0]; m4 = rep4["mounts"][0]
    assert m4["retained_new_gb"] > m1["retained_new_gb"], "RS13 retained"
    assert m4["max_concurrent_temporary_gb"] == m1["max_concurrent_temporary_gb"], \
        "RS13: scratch peak must not scale with cohort size at fixed concurrency"
    print("resources.demo: ok")


if __name__ == "__main__":
    demo()

#!/usr/bin/env python3
"""Self-check for storage planning (P1 / R15, U07; V16, V17).

Master plan 10.1: "There is no default 20 GB limit and no sample-count cutoff
that rejects an analysis." 10.7 requires hpc.storage_budget_gb and the fixed
per-assay FASTQ assumptions to be replaced by measured planning.

What was there: one platform-wide cap, per-assay FASTQ sizes guessed from a
label, and concurrency derived from sample count. None of it looked at the data
or at the filesystem.

What must not break:
  - inputs are measured where they exist and counted once per unique file (RS01)
  - an unmeasurable input is labelled, never silently treated as zero (RS02)
  - more libraries at fixed concurrency grow retained data, not scratch (RS13)
  - more concurrency grows scratch (RS04)
  - a cache already present is reported but not charged again (RS05)
  - paths sharing a filesystem are accounted once (RS07)
  - unknown capacity is never reported as unlimited
  - concurrency is reduced before the run is refused (RS08), and a genuinely
    impossible run is refused with the shortfall (RS09)
  - a legacy storage_budget_gb warns and acts as an explicit quota (RS12)
  - no default cap is reintroduced anywhere

Run:  python3 tests/check_resources.py
"""
import importlib.util
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent
GB = 1024 ** 3

spec = importlib.util.spec_from_file_location("resources", ROOT / "scripts" / "resources.py")
rs = importlib.util.module_from_spec(spec)
# Register before executing: resources.py defines a dataclass, and dataclasses
# resolve their module through sys.modules at class-creation time.
sys.modules["resources"] = rs
spec.loader.exec_module(rs)


def tmpdir():
    return Path(tempfile.mkdtemp())


# --- 1. measurement (RS01, RS02) --------------------------------------
d = tmpdir()
f = d / "a.fq.gz"
f.write_bytes(b"x" * 4096)
inputs = rs.measure_inputs([str(f), str(f), "s3://b/x.fq.gz", str(d / "gone.fq")])
assert len(inputs) == 3, "a repeated URI must be counted once (RS01)"
assert inputs[str(f)] == {"bytes": 4096, "source": "measured"}
assert inputs["s3://b/x.fq.gz"]["bytes"] is None
assert inputs["s3://b/x.fq.gz"]["source"] == "unknown_remote"
assert inputs[str(d / "gone.fq")]["source"] == "unknown_missing"

# an unknown size is labelled in the report, not silently zero
d2 = tmpdir(); (d2 / "res").mkdir()
rep = rs.plan_storage(repo_root=ROOT, inputs=inputs, n_libraries=3,
                      assay="bulk", concurrency=1,
                      paths={"results": d2 / "res", "scratch": d2 / "res"})
assert "estimated at the upper band" in rep["uncertainty"], rep["uncertainty"]
rep_all = rs.plan_storage(repo_root=ROOT, inputs={str(f): {"bytes": 1, "source": "measured"}},
                          n_libraries=1, assay="bulk", concurrency=1,
                          paths={"results": d2 / "res"})
assert rep_all["uncertainty"] == "measured", rep_all["uncertainty"]


def plan(n_libraries, concurrency, *, assay="bulk_end_tag", paths=None,
         quotas=None, cache_present=None, retain_alignments=False):
    base = tmpdir()
    for sub in ("res", "scr", "cache"):
        (base / sub).mkdir()
    return rs.plan_storage(
        repo_root=ROOT, inputs={}, n_libraries=n_libraries, assay=assay,
        concurrency=concurrency,
        paths=paths or {"results": base / "res", "scratch": base / "scr",
                        "cache": base / "cache"},
        quotas=quotas, cache_present=cache_present,
        retain_alignments=retain_alignments)


# --- 2. scaling behaviour (RS03, RS04, RS13) --------------------------
small, big = plan(4, 2), plan(64, 2)
m_s, m_b = small["mounts"][0], big["mounts"][0]
assert m_b["retained_new_gb"] > m_s["retained_new_gb"], "RS03: retained must grow"
assert m_b["max_concurrent_temporary_gb"] == m_s["max_concurrent_temporary_gb"], (
    "RS13: scratch peak must not scale with cohort size at fixed concurrency")
assert m_b["cache_new_gb"] == m_s["cache_new_gb"], "RS03: cache is counted once"

wide = plan(64, 8)
assert wide["mounts"][0]["max_concurrent_temporary_gb"] > m_b["max_concurrent_temporary_gb"], \
    "RS04: more concurrent jobs must raise the scratch reservation"

# concurrency cannot exceed the number of libraries
assert plan(2, 16)["planned_concurrency"] == 2

# retained alignments are their own cost, not predicted by sample count (RS06)
assert plan(4, 2, retain_alignments=True)["mounts"][0]["retained_new_gb"] > \
       plan(4, 2)["mounts"][0]["retained_new_gb"]

# --- 3. an existing cache is reported, not charged (RS05) -------------
fresh = plan(4, 2)
warm = plan(4, 2, cache_present={"reference": True, "index_build": True,
                                 "environment": True})
assert warm["mounts"][0]["cache_new_gb"] == 0, warm["mounts"][0]
assert warm["mounts"][0]["existing_gb"] > 0, "existing cache must still be reported"
assert warm["mounts"][0]["required_free_gb"] < fresh["mounts"][0]["required_free_gb"]

# --- 4. one filesystem is accounted once (RS07) -----------------------
base = tmpdir()
for sub in ("a", "b", "c"):
    (base / sub).mkdir()
one_fs = rs.plan_storage(repo_root=ROOT, inputs={}, n_libraries=4,
                         assay="bulk_end_tag", concurrency=2,
                         paths={"results": base / "a", "scratch": base / "b",
                                "cache": base / "c"})
assert len(one_fs["mounts"]) == 1, "three paths on one filesystem must be one mount"
assert len(one_fs["mounts"][0]["paths"]) == 3, one_fs["mounts"][0]

# --- 5. unknown capacity is never "unlimited" -------------------------
rep = rs.plan_storage(repo_root=ROOT, inputs={}, n_libraries=1,
                      assay="bulk", concurrency=1,
                      paths={"results": Path("/nonexistent-mount-xyz/results")})
m = rep["mounts"][0]
assert m["available_gb"] is None or m["verdict"] in ("ok", "insufficient")
if m["available_gb"] is None:
    assert m["verdict"] == "unknown_capacity", m
    assert rep["blocking"], "unknown capacity must be recorded as blocking"

# --- 6. a quota below free space wins (master plan 10.3) --------------
tight = plan(4, 2, quotas={"results": 1 * GB})
m = tight["mounts"][0]
assert m["available_gb"] == 1.0, m
assert m["verdict"] == "insufficient" and m["shortfall_gb"] > 0, m

# --- 7. RS08 reduces concurrency before refusing ----------------------
base = tmpdir(); (base / "res").mkdir(); (base / "scr").mkdir()
kwargs = dict(repo_root=ROOT, inputs={}, n_libraries=32, assay="bulk",
              concurrency=16,
              paths={"results": base / "res", "scratch": base / "scr"},
              quotas={"results": 40 * GB})
first = rs.plan_storage(**kwargs)
reduced = rs.reduce_concurrency(first, **kwargs)
assert reduced["planned_concurrency"] <= first["planned_concurrency"], reduced
# the route is unchanged: same assay, same library count
assert reduced["assay"] == first["assay"] and reduced["n_libraries"] == first["n_libraries"]

# --- 8. end to end through plan_resources.py, including RS12 ----------
def project(with_legacy_cap):
    p = tmpdir() / "proj"
    (p / "config").mkdir(parents=True)
    (p / "fastq").mkdir()
    for i in range(4):
        (p / "fastq" / f"S{i}.fq.gz").write_bytes(b"x" * (1024 * 1024))
    hpc = ("hpc: {storage_budget_gb: 20, samples_in_flight: null}"
           if with_legacy_cap else "hpc: {samples_in_flight: null}")
    (p / "config" / "config.yaml").write_text(
        'project: {name: t, output_dir: "results/"}\n'
        "organism: {scientific_name: Vitis vinifera, tax_id: 29760}\n"
        "reference: {accession: GCF_1.1}\n"
        "samples: {sheet: config/samples.tsv, seq_type: tagseq}\n"
        'model: {fixed_effects: "~ cond", random_effects: null, primary_factor: cond}\n'
        + hpc + "\n")
    rows = ["sample_id\tfastq_url\tcond"]
    for i in range(4):
        rows.append(f"S{i}\t{p / 'fastq' / f'S{i}.fq.gz'}\t"
                    f"{'ctl' if i < 2 else 'trt'}")
    (p / "config" / "samples.tsv").write_text("\n".join(rows) + "\n")
    return p


p = project(with_legacy_cap=False)
out_json = p / "gates" / "resource_plan.json"
r = subprocess.run([sys.executable, str(ROOT / "scripts" / "plan_resources.py"),
                    "-d", str(p), "--out", str(out_json)],
                   capture_output=True, text=True)
assert r.returncode == 0, r.stdout + r.stderr
rep = json.loads(out_json.read_text())
assert rep["uncertainty"] == "measured", rep["uncertainty"]
assert rep["n_libraries"] == 4 and rep["assay"] == "bulk_end_tag", rep
assert rep["legacy_notes"] == [], rep["legacy_notes"]
assert "conservative upper bound" in rep["basis"]

p = project(with_legacy_cap=True)
out_json = p / "gates" / "resource_plan.json"
r = subprocess.run([sys.executable, str(ROOT / "scripts" / "plan_resources.py"),
                    "-d", str(p), "--out", str(out_json)],
                   capture_output=True, text=True)
rep = json.loads(out_json.read_text())
assert rep["legacy_notes"], "RS12: a legacy cap must produce a migration note"
joined = " ".join(rep["legacy_notes"])
assert "no longer" in joined and "NOT a planning input" in joined, joined
assert "must not be silently erased" in joined, joined
# and it is applied, not ignored: 20 GB cannot hold reference + index + env
assert any(m["quota_known"] for m in rep["mounts"]), rep["mounts"]
assert r.returncode == 1, "an unsatisfiable quota must block (RS09)"
assert "short by" in (r.stdout + r.stderr)

# --- 9. no default cap is reintroduced --------------------------------
tpl = yaml.safe_load((ROOT / "config" / "config.template.yaml").read_text())
assert tpl["hpc"].get("storage_budget_gb") is None, (
    "config.template.yaml reinstated a default storage cap")
models = yaml.safe_load((ROOT / "config" / "resource_models.yaml").read_text())
flat = json.dumps(models)
assert "storage_budget" not in flat, "resource_models.yaml must not define a cap"
submit = (ROOT / "submit.sh").read_text()
# Strip comments: submit.sh names the setting in the comment explaining that it
# no longer reads it.
submit_code = "\n".join(l for l in submit.splitlines()
                        if not l.strip().startswith("#"))
assert "storage_budget_gb" not in submit_code, (
    "submit.sh still reads a storage budget instead of planning")
assert "plan_resources.py" in submit_code, "submit.sh does not run the storage plan"
# sample-count tiers must no longer drive concurrency (master plan 10.7)
assert "fastq_size_estimate_gb" not in submit_code, (
    "submit.sh still uses per-assay FASTQ guesses")
assert "reserved_gb" not in submit_code, "submit.sh still uses reserved_gb"

# every coefficient declares whether it was measured or assumed (RS02)
for group in ("assays", "cache", "cohort"):
    for name, entry in models[group].items():
        vals = entry.values() if "low" not in entry else [entry]
        for v in vals:
            if isinstance(v, dict) and "low" in v:
                assert "source" in v, f"{group}.{name} has no source label"

print("check_resources.py: all assertions passed")

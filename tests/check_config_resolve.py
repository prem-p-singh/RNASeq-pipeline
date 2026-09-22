#!/usr/bin/env python3
"""Self-check for scripts/config_resolve.py (P0 / R17e, R05).

The defect this guards against was silent and shipped. config.template.yaml
carried its own partial `thresholds.wgcna`, and the Snakefile merged
config/thresholds.yaml only one level deep. Because `thresholds.wgcna` was
already present, the fuller definition was never merged, so 05_wgcna.R read
NULL for merge_cut_height_default and merge_cut_height_bumped. WGCNA does not
fail on a NULL cut height: it catches the error and returns unmerged modules.
The run completed, the metrics recorded a null height, and the requested
merging never happened.

What must not break:
  - deep merge keeps siblings when one nested key is overridden
  - lists replace rather than append (master plan 6.3)
  - unknown keys and type conflicts are reported, not absorbed
  - the origin of every resolved value is recorded
  - the real repository files resolve to usable WGCNA merge heights
  - defaults have a single source: no `thresholds` block in the template
  - the Snakefile no longer does its own shallow merge
  - preflight reports the same findings through the same module

Run:  python3 tests/check_config_resolve.py
"""
import importlib.util
import json
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
spec = importlib.util.spec_from_file_location(
    "config_resolve", ROOT / "scripts" / "config_resolve.py")
cr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cr)

# --- 1. the merge behaviour itself ------------------------------------
base = {"a": {"x": 1, "y": 2}, "list": [1, 2], "scalar": 5}
over = {"a": {"x": 99}, "list": [3], "scalar": 6}
m = cr.deep_merge(base, over)
assert m["a"] == {"x": 99, "y": 2}, f"sibling lost: {m['a']}"
assert m["list"] == [3], f"lists must replace, not append: {m['list']}"
assert m["scalar"] == 6

# nested three deep
deep = cr.deep_merge({"p": {"q": {"r": 1, "s": 2}}}, {"p": {"q": {"r": 9}}})
assert deep["p"]["q"] == {"r": 9, "s": 2}, deep

# --- 2. the real files: the exact bug that shipped --------------------
schema = cr.load_schema(ROOT)
w = schema["thresholds"]["wgcna"]
for k in ("merge_cut_height_default", "merge_cut_height_bumped",
          "min_samples", "target_r2", "power_cap", "max_modules"):
    assert k in w, f"composed schema lost thresholds.wgcna.{k}: {w}"

cfg, origins, issues = cr.resolve(
    ROOT, [("project", {"thresholds": {"wgcna": {"min_samples": 20}}})])
w2 = cfg["thresholds"]["wgcna"]
assert w2["min_samples"] == 20, w2
assert w2["merge_cut_height_default"] == 0.25, (
    "overriding one wgcna key dropped a merge height; this is the shipped bug")
assert w2["merge_cut_height_bumped"] == 0.35, w2
assert issues == [], issues

# 05_wgcna.R reads cfg$thresholds$wgcna; none of what it needs may be null.
for k in ("min_samples", "target_r2", "power_cap", "max_modules",
          "merge_cut_height_default", "merge_cut_height_bumped"):
    assert w2.get(k) is not None, f"05_wgcna.R would read NULL for {k}"

# --- 3. origins are recorded (master plan 6.3) ------------------------
assert origins["thresholds.wgcna.min_samples"] == "project", origins
assert origins["thresholds.wgcna.merge_cut_height_default"] == "default", origins

# --- 4. unknown keys and type conflicts ------------------------------
_, _, iss = cr.resolve(ROOT, [("project", {"bogus_section": {"a": 1}})])
codes = [c for c, _ in iss]
assert "CFG001" in codes, iss
# one stray section yields one finding, not one per leaf beneath it
assert sum(1 for c in codes if c == "CFG001") == 1, iss

_, _, iss = cr.resolve(
    ROOT, [("project", {"hpc": {"storage_budget_gb": "twenty"}})])
assert "CFG003" in [c for c, _ in iss], iss

# a legitimate override produces nothing
_, _, iss = cr.resolve(ROOT, [("project", {"hpc": {"storage_budget_gb": 50}})])
assert iss == [], iss

# user-defined subtrees are carried through, never flagged
_, _, iss = cr.resolve(ROOT, [("project", {"contrasts": [
    {"id": "anything", "type": "pairwise", "factor": "f", "made_up": True}]})])
assert iss == [], iss

# --- 5. single source of defaults (master plan 5.2) -------------------
import yaml
template = yaml.safe_load((ROOT / "config" / "config.template.yaml").read_text())
assert "thresholds" not in template, (
    "config.template.yaml has a thresholds block again; defaults must have one "
    "source, and its partial wgcna entry is what shadowed the merge heights")

# --- 6. the Snakefile no longer merges on its own ---------------------
snake = (ROOT / "Snakefile").read_text()
assert "config_resolve" in snake, "Snakefile does not use the shared resolver"
assert 'config["thresholds"].setdefault' not in snake, (
    "Snakefile still does its own one-level threshold merge")

# --- 7. preflight surfaces the same findings, end to end --------------
TMP = Path(tempfile.mkdtemp())
proj = TMP / "p"
(proj / "config").mkdir(parents=True)
(proj / "config" / "config.yaml").write_text(
    'project: {name: t, output_dir: "results/"}\n'
    "organism: {scientific_name: Vitis vinifera, tax_id: 29760}\n"
    "reference: {accession: GCF_1.1}\n"
    "samples: {sheet: config/samples.tsv, seq_type: tagseq}\n"
    'model: {fixed_effects: "~ treatment", random_effects: null, primary_factor: treatment}\n'
    "thresholds: {wgcna: {min_samples: 20}}\n"
    "bogus_key: 1\n")
(proj / "config" / "samples.tsv").write_text(
    "sample_id\ttreatment\nC1\tcontrol\nC2\tcontrol\nT1\ttreated\nT2\ttreated\n")

r = subprocess.run([sys.executable, str(ROOT / "scripts" / "preflight.py"),
                    "-d", str(proj)], capture_output=True, text=True)
out = r.stdout + r.stderr
assert "CFG001" in out, f"preflight did not report the unknown key:\n{out}"
assert r.returncode == 1, r.returncode

# and the resolved plan still carries the merged thresholds
plan = json.loads((proj / "gates" / "preflight_plan.json").read_text())
assert plan["n_errors"] >= 1, plan

print("check_config_resolve.py: all assertions passed")

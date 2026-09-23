#!/usr/bin/env python3
"""Self-check for the samples/libraries/reads contract (P0-P1 / R04, R17b, R05).

Master plan 6.1 defines three records where the pipeline had one sheet. With
one row per sample there is nowhere to put a second lane, which is why R04
rejects multi-lane libraries outright, and nothing can express that two
libraries came from one specimen.

What must not break (rules IN01-IN04, master plan 6.2):
  - duplicate identifiers and foreign keys that point nowhere are named, with
    the exact offending records
  - one file may back only one read unit
  - a paired library missing a mate is caught per lane, never downgraded
  - several lanes of one library are valid and are grouped in mate order
  - several libraries from one specimen are recorded as not independent
  - legacy samples.tsv converts, with the assumptions reported rather than hidden
  - preflight validates the three tables when present and the single sheet when
    not, so existing projects keep running

Run:  python3 tests/check_metadata.py
"""
import importlib.util
import json
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
spec = importlib.util.spec_from_file_location("metadata", ROOT / "scripts" / "metadata.py")
md = importlib.util.module_from_spec(spec)
spec.loader.exec_module(md)


def codes(tables):
    return [c for c, _ in md.validate(tables)]


def base():
    return {
        "samples": [{"sample_id": "S1", "biological_unit_id": "P1", "condition": "ctl"},
                    {"sample_id": "S2", "biological_unit_id": "P2", "condition": "trt"}],
        "libraries": [{"library_id": "L1", "sample_id": "S1", "layout": "paired"},
                      {"library_id": "L2", "sample_id": "S2", "layout": "paired"}],
        "reads": [
            {"read_unit_id": "r1", "library_id": "L1", "role": "R1", "uri": "/d/a_R1.fq", "lane": "1"},
            {"read_unit_id": "r2", "library_id": "L1", "role": "R2", "uri": "/d/a_R2.fq", "lane": "1"},
            {"read_unit_id": "r3", "library_id": "L2", "role": "R1", "uri": "/d/b_R1.fq", "lane": "1"},
            {"read_unit_id": "r4", "library_id": "L2", "role": "R2", "uri": "/d/b_R2.fq", "lane": "1"},
        ],
    }


# --- 1. a correct set is silent ---------------------------------------
assert md.validate(base()) == [], md.validate(base())

# --- 2. the templates shipped with the repo are valid ------------------
tpl = {n: md.read_tsv(ROOT / "templates" / f"{n}.tsv")
       for n in ("samples", "libraries", "reads")}
assert all(tpl.values()), "a template is empty"
assert md.validate(tpl) == [], md.validate(tpl)

# --- 3. duplicate identifiers, named (IN01) ---------------------------
t = base()
t["libraries"].append({"library_id": "L1", "sample_id": "S2", "layout": "single"})
iss = md.validate(t)
assert "MET008" in [c for c, _ in iss], iss
detail = [d for c, d in iss if c == "MET008"][0]
assert "L1" in detail and "row" in detail, detail

# --- 4. foreign keys that point nowhere (IN01) ------------------------
t = base()
t["libraries"][0]["sample_id"] = "S99"
iss = md.validate(t)
assert "MET007" in [c for c, _ in iss], iss
assert "S99" in [d for c, d in iss if c == "MET007"][0]

t = base()
t["reads"][0]["library_id"] = "L99"
assert "MET007" in codes(t)

# --- 5. one file may back only one read unit (IN01) -------------------
t = base()
t["reads"][2]["uri"] = t["reads"][0]["uri"]
iss = md.validate(t)
assert "MET009" in [c for c, _ in iss], iss
d = [x for c, x in iss if c == "MET009"][0]
assert "r1" in d and "r3" in d, d

# --- 6. a paired library missing a mate (IN02) ------------------------
t = base()
t["reads"] = [r for r in t["reads"] if r["read_unit_id"] != "r2"]
iss = md.validate(t)
assert "MET010" in [c for c, _ in iss], iss
assert "R2" in [d for c, d in iss if c == "MET010"][0]

# a paired library must not be expressed as single-end
t = base()
t["reads"][1]["role"] = "single"
assert "MET013" in codes(t)

# and a single-end library must not carry an R2
t = base()
t["libraries"][0]["layout"] = "single"
t["reads"] = [r for r in t["reads"] if r["read_unit_id"] != "r2"]
t["reads"].append({"read_unit_id": "r9", "library_id": "L1", "role": "R2",
                   "uri": "/d/x_R2.fq", "lane": "1"})
assert "MET013" in codes(t)

# --- 7. multi-lane is valid, and mates are checked per lane (IN03) ----
t = base()
t["reads"] += [
    {"read_unit_id": "r5", "library_id": "L1", "role": "R1", "uri": "/d/a2_R1.fq", "lane": "2"},
    {"read_unit_id": "r6", "library_id": "L1", "role": "R2", "uri": "/d/a2_R2.fq", "lane": "2"},
]
assert md.validate(t) == [], md.validate(t)
groups = md.lane_groups(t)
assert [g["lane"] for g in groups["L1"]] == ["1", "2"], groups["L1"]
assert [r["role"] for r in groups["L1"][1]["reads"]] == ["R1", "R2"]

# Lane 1 in distinct runs must survive grouping as separate read pairs.
t_runs = base()
for read in t_runs["reads"]:
    read["run"] = "run_a"
for read in t_runs["reads"][:2]:
    t_runs["reads"].append(dict(read, read_unit_id=read["read_unit_id"] + "b",
                                uri=read["uri"] + ".b", run="run_b"))
assert md.validate(t_runs) == []
groups = md.lane_groups(t_runs)["L1"]
assert [(g["run"], g["lane"]) for g in groups] == [("run_a", "1"), ("run_b", "1")]
assert sum(len(g["reads"]) for g in groups) == 4
# Duplicate roles remain visible in grouping and fail validation.
t_runs["reads"].append(dict(t_runs["reads"][0], read_unit_id="duplicate", uri="/other"))
assert "MET013" in codes(t_runs)
assert len(md.lane_groups(t_runs)["L1"][0]["reads"]) == 3

for name in ("samples", "libraries", "reads"):
    incomplete = base()
    incomplete[name] = []
    assert "MET001" in codes(incomplete), name
for name, field in (("libraries", "sample_id"), ("libraries", "layout"),
                    ("reads", "library_id"), ("reads", "uri"), ("reads", "role")):
    incomplete = base()
    incomplete[name][0][field] = ""
    assert "MET003" in codes(incomplete), (name, field)
t_empty = base()
t_empty["reads"] = t_empty["reads"][2:]
assert "MET013" in codes(t_empty)
t_empty = base()
t_empty["libraries"][0]["layout"] = "single"
t_empty["reads"][1]["role"] = "single"
assert "MET013" in codes(t_empty)

# a lane missing its mate is caught, and the message names that lane
t["reads"] = [r for r in t["reads"] if r["read_unit_id"] != "r6"]
iss = md.validate(t)
det = [d for c, d in iss if c == "MET010"]
assert det and "lane 2" in det[0], det

# --- 8. several libraries from one specimen (IN04) --------------------
t = base()
t["libraries"].append({"library_id": "L1b", "sample_id": "S1", "layout": "paired"})
t["reads"] += [
    {"read_unit_id": "r7", "library_id": "L1b", "role": "R1", "uri": "/d/c_R1.fq", "lane": "1"},
    {"read_unit_id": "r8", "library_id": "L1b", "role": "R2", "uri": "/d/c_R2.fq", "lane": "1"},
]
iss = md.validate(t)
assert "MET012" in [c for c, _ in iss], iss
d = [x for c, x in iss if c == "MET012"][0]
assert "not" in d and "independent" in d, d
# it is a warning, not a blocker: the relationship is preserved, not rejected
spec_yaml = __import__("yaml").safe_load((ROOT / "config" / "spec.yaml").read_text())
assert spec_yaml["issues"]["MET012"]["severity"] == "warning"

# --- 9. unknown role, blank key ---------------------------------------
t = base(); t["reads"][0]["role"] = "mate1"
assert "MET011" in codes(t)
t = base(); t["samples"][0]["sample_id"] = ""
assert "MET003" in codes(t)

# --- 10. legacy conversion, with its assumptions stated ---------------
legacy = [
    {"sample_id": "A", "fastq_url": "/d/A_R1.fq.gz", "fastq_url_r2": "/d/A_R2.fq.gz",
     "treatment": "ctl"},
    {"sample_id": "B", "fastq_url": "/d/B_R1.fq.gz", "fastq_url_r2": "/d/B_R2.fq.gz",
     "treatment": "trt"},
]
tables, notes = md.from_legacy(legacy, "rnaseq_paired")
assert md.validate(tables) == [], md.validate(tables)
assert len(tables["reads"]) == 4 and len(tables["libraries"]) == 2
assert tables["libraries"][0]["layout"] == "paired"
assert tables["samples"][0]["treatment"] == "ctl", "covariate not carried"
assert any("assumption" in n for n in notes), notes

# single-end legacy
tables, notes = md.from_legacy(
    [{"sample_id": "A", "fastq_url": "/d/A.fq.gz", "fastq_url_r2": ""}], "tagseq")
assert md.validate(tables) == [], md.validate(tables)
assert tables["reads"][0]["role"] == "single"
assert tables["libraries"][0]["assay_family"] == "bulk_end_tag"
assert any("kit is unknown" in n for n in notes), notes

# a paired legacy sheet with a blank R2 converts, and the gap is reported
# rather than quietly becoming single-end
tables, notes = md.from_legacy(
    [{"sample_id": "A", "fastq_url": "/d/A_R1.fq.gz", "fastq_url_r2": ""}],
    "rnaseq_paired")
assert any("fastq_url_r2 is empty" in n for n in notes), notes
assert tables["libraries"][0]["layout"] == "paired", "must not downgrade (IN02)"

# --- 11. round trip through files -------------------------------------
tmp = Path(tempfile.mkdtemp())
tables = base()
md.write_tables(tables, tmp / "metadata")
back = md.load_tables(tmp / "metadata")
assert {k: len(v) for k, v in back.items()} == {"samples": 2, "libraries": 2, "reads": 4}
assert md.validate(back) == [], md.validate(back)

# --- 12. preflight uses the tables when present, the sheet when not ----
def project(with_tables):
    p = Path(tempfile.mkdtemp()) / "proj"
    (p / "config").mkdir(parents=True)
    (p / "config" / "config.yaml").write_text(
        'project: {name: t, output_dir: "results/"}\n'
        "organism: {scientific_name: Vitis vinifera, tax_id: 29760}\n"
        "reference: {accession: GCF_1.1}\n"
        "samples: {sheet: config/samples.tsv, seq_type: tagseq}\n"
        'model: {fixed_effects: "~ condition", random_effects: null, primary_factor: condition}\n')
    (p / "config" / "samples.tsv").write_text(
        "sample_id\tcondition\nS1\tctl\nS2\tctl\nS3\ttrt\nS4\ttrt\n")
    if with_tables:
        t = base()
        # break a foreign key so preflight must report it
        t["libraries"][0]["sample_id"] = "MISSING"
        md.write_tables(t, p / "metadata")
    return p


p = project(with_tables=False)
r = subprocess.run([sys.executable, str(ROOT / "scripts" / "preflight.py"), "-d", str(p)],
                   capture_output=True, text=True)
plan = json.loads((p / "gates" / "preflight_plan.json").read_text())
assert plan["metadata_tables"] is None, "single-sheet project must not claim tables"
assert r.returncode == 0, r.stdout + r.stderr

p = project(with_tables=True)
r = subprocess.run([sys.executable, str(ROOT / "scripts" / "preflight.py"), "-d", str(p)],
                   capture_output=True, text=True)
out = r.stdout + r.stderr
assert "MET007" in out, out
assert r.returncode == 1
plan = json.loads((p / "gates" / "preflight_plan.json").read_text())
assert plan["metadata_tables"]["form"] == "three_table", plan["metadata_tables"]
assert plan["metadata_tables"]["counts"]["reads"] == 4, plan["metadata_tables"]

# Empty/header-only canonical files must not silently fall back to legacy.
p = project(with_tables=False)
(p / "metadata").mkdir()
(p / "metadata/samples.tsv").write_text("sample_id\n")
r = subprocess.run([sys.executable, str(ROOT / "scripts/preflight.py"), "-d", str(p)],
                   capture_output=True, text=True)
assert r.returncode == 1 and "MET001" in r.stdout + r.stderr

print("check_metadata.py: all assertions passed")

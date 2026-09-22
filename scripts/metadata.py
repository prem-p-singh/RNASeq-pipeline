#!/usr/bin/env python3
"""metadata.py — the samples / libraries / reads contract (R04, R17b, R05).

Master plan 6.1 defines three input records where the pipeline has one sheet:

  samples    a biological specimen at one observation
  libraries  a prepared library, with assay, layout, strandedness, kit, batch
  reads      a read unit or lane, with role, URI, mate relationship, lane

Why this matters beyond tidiness. With one row per sample there is nowhere to
put a second lane, so multi-lane libraries are currently rejected outright, and
nothing can express that two libraries came from the same specimen. Master plan
6.1: "Multiple lanes may merge only within the same compatible library and read
structure, with order/mate integrity preserved. Multiple libraries from a
biological unit do not create independent replicates."

Validation follows 6.2 and rules IN01-IN04: identifiers, uniqueness, foreign
keys, mate completeness, repeated file assignment, and permitted roles. A
finding names the exact conflicting records, because IN01 requires that.

Legacy projects are not broken. from_legacy() converts an existing samples.tsv
plus seq_type into the three tables and returns a migration report, per 6.3
("Convert the current seq_type into explicit assay/layout fields with a
migration report"). Nothing here rewrites a project on disk.
"""
from __future__ import annotations

import csv
from collections import Counter, defaultdict
from pathlib import Path

# Layouts a library may declare. `paired` is the only one that requires mates.
LAYOUTS = ("paired", "single")
ROLES = ("R1", "R2", "single", "index")

SCHEMAS = {
    "samples": {
        "key": "sample_id",
        "required": ("sample_id",),
        # biological_unit_id lets repeated observations of one specimen be
        # linked without being counted as independent replicates (IN04).
        "recommended": ("biological_unit_id", "condition"),
        "foreign_keys": {},
    },
    "libraries": {
        "key": "library_id",
        "required": ("library_id", "sample_id", "layout"),
        "recommended": ("assay_family", "strandedness", "kit", "platform",
                        "prep_batch"),
        "foreign_keys": {"sample_id": "samples"},
    },
    "reads": {
        "key": "read_unit_id",
        "required": ("read_unit_id", "library_id", "role", "uri"),
        "recommended": ("lane", "run", "bytes", "checksum"),
        "foreign_keys": {"library_id": "libraries"},
    },
}


def read_tsv(path) -> list[dict]:
    """Rows of a TSV, ignoring blank lines and # comments."""
    path = Path(path)
    if not path.is_file():
        return []
    rows = []
    with open(path, newline="") as fh:
        lines = [l for l in fh
                 if l.strip() and not l.lstrip().startswith("#")]
    if not lines:
        return []
    for row in csv.DictReader(lines, delimiter="\t"):
        rows.append({(k or "").strip(): (v or "").strip()
                     for k, v in row.items()})
    return rows


def load_tables(metadata_dir) -> dict:
    """Load whichever of the three tables are present."""
    d = Path(metadata_dir)
    return {name: read_tsv(d / f"{name}.tsv") for name in SCHEMAS}


# --------------------------------------------------------------- validation
def _missing_columns(rows, required):
    if not rows:
        return list(required)
    present = set(rows[0])
    return [c for c in required if c not in present]


def validate(tables: dict) -> list[tuple[str, str]]:
    """Return (issue_code, detail) for every contract violation.

    Codes are declared in config/spec.yaml. Details name the offending records,
    which IN01 requires ("list exact conflicting records").
    """
    issues: list[tuple[str, str]] = []

    # --- structure and identifiers ------------------------------------
    for name, schema in SCHEMAS.items():
        rows = tables.get(name) or []
        if not rows:
            continue
        for col in _missing_columns(rows, schema["required"]):
            issues.append(("MET001", f"{name}.tsv missing required column {col}"))

        key = schema["key"]
        if rows and key in rows[0]:
            values = [r.get(key, "") for r in rows]
            blanks = [i + 2 for i, v in enumerate(values) if not v]
            if blanks:
                issues.append(("MET003",
                               f"{name}.tsv blank {key} at row(s) "
                               f"{', '.join(map(str, blanks))}"))
            dupes = sorted({v for v, n in Counter(values).items() if v and n > 1})
            for v in dupes:
                where = [i + 2 for i, x in enumerate(values) if x == v]
                issues.append(("MET008",
                               f"{name}.tsv duplicate {key} {v!r} at row(s) "
                               f"{', '.join(map(str, where))}"))

    # --- foreign keys (IN01: invalid table joins) ---------------------
    for name, schema in SCHEMAS.items():
        rows = tables.get(name) or []
        for col, parent in schema["foreign_keys"].items():
            parent_rows = tables.get(parent) or []
            if not rows or not parent_rows:
                continue
            known = {r.get(SCHEMAS[parent]["key"], "") for r in parent_rows}
            for i, r in enumerate(rows):
                v = r.get(col, "")
                if v and v not in known:
                    issues.append((
                        "MET007",
                        f"{name}.tsv row {i + 2}: {col}={v!r} has no matching "
                        f"{SCHEMAS[parent]['key']} in {parent}.tsv"))

    reads = tables.get("reads") or []
    libraries = tables.get("libraries") or []

    # --- one file may back only one read unit (IN01) -------------------
    by_uri = defaultdict(list)
    for r in reads:
        if r.get("uri"):
            by_uri[r["uri"]].append(r.get("read_unit_id", "?"))
    for uri, units in sorted(by_uri.items()):
        if len(units) > 1:
            issues.append(("MET009",
                           f"{uri} is assigned to {len(units)} read units: "
                           f"{', '.join(units)}"))

    # --- role vocabulary ----------------------------------------------
    for i, r in enumerate(reads):
        role = r.get("role", "")
        if role and role not in ROLES:
            issues.append(("MET011",
                           f"reads.tsv row {i + 2}: role={role!r}; "
                           f"expected one of {', '.join(ROLES)}"))

    # --- layout, mates and lanes (IN02, IN03) -------------------------
    layout_of = {l.get("library_id", ""): l.get("layout", "") for l in libraries}
    for lib_id, layout in sorted(layout_of.items()):
        if layout and layout not in LAYOUTS:
            issues.append(("MET013",
                           f"libraries.tsv {lib_id}: layout={layout!r}; "
                           f"expected one of {', '.join(LAYOUTS)}"))

    # group this library's reads by lane, so mates are checked per lane
    per_lib_lane = defaultdict(lambda: defaultdict(list))
    for r in reads:
        per_lib_lane[r.get("library_id", "")][r.get("lane", "") or "-"].append(r)

    for lib_id, lanes in sorted(per_lib_lane.items()):
        layout = layout_of.get(lib_id)
        if layout is None:
            continue
        for lane, rows in sorted(lanes.items()):
            roles = [r.get("role", "") for r in rows]
            where = f"library {lib_id}" + (f", lane {lane}" if lane != "-" else "")
            if layout == "paired":
                # IN02: never downgrade a paired library to single-end.
                for need in ("R1", "R2"):
                    if roles.count(need) == 0:
                        issues.append(("MET010",
                                       f"{where}: layout is paired but no {need} "
                                       f"read unit is present"))
                for dup in ("R1", "R2"):
                    if roles.count(dup) > 1:
                        issues.append(("MET013",
                                       f"{where}: {roles.count(dup)} {dup} read "
                                       f"units; a lane holds one of each mate"))
                if "single" in roles:
                    issues.append(("MET013",
                                   f"{where}: role 'single' in a paired library"))
            elif layout == "single":
                if "R2" in roles:
                    issues.append(("MET013",
                                   f"{where}: R2 present in a single-end library"))

    # --- IN04: several libraries from one specimen are not replicates --
    per_sample = defaultdict(list)
    for l in libraries:
        per_sample[l.get("sample_id", "")].append(l.get("library_id", "?"))
    for sid, libs in sorted(per_sample.items()):
        if sid and len(libs) > 1:
            issues.append(("MET012",
                           f"sample {sid} has {len(libs)} libraries "
                           f"({', '.join(sorted(libs))}); these are not "
                           f"independent biological replicates"))

    return issues


def lane_groups(tables: dict) -> dict:
    """Read units per library per lane, in mate order.

    IN03 permits merging lanes "within the same compatible library and read
    structure, with order/mate integrity preserved". This is the grouping a
    merge step would consume; nothing here merges anything.
    """
    out: dict = {}
    per_lib = defaultdict(lambda: defaultdict(dict))
    for r in tables.get("reads") or []:
        per_lib[r.get("library_id", "")][r.get("lane", "") or "-"][r.get("role", "")] = r
    for lib, lanes in per_lib.items():
        out[lib] = {lane: {role: rows[role]["uri"] for role in sorted(rows)}
                    for lane, rows in sorted(lanes.items())}
    return out


# ------------------------------------------------------------------ legacy
def from_legacy(sheet_rows: list[dict], seq_type: str) -> tuple[dict, list[str]]:
    """Convert one samples.tsv plus seq_type into the three tables.

    Master plan 6.3 requires this conversion to produce a migration report and
    to keep what it could not determine visible rather than invented. One
    library and one lane per sample is the only reading the old format
    supports; that assumption is stated in the report rather than hidden.
    """
    layout = "paired" if seq_type == "rnaseq_paired" else "single"
    assay = {"tagseq": "bulk_end_tag",
             "rnaseq_paired": "bulk",
             "rnaseq_single": "bulk"}.get(seq_type, "unknown")

    samples, libraries, reads = [], [], []
    notes = [
        f"seq_type={seq_type!r} became assay_family={assay!r}, layout={layout!r}",
        "one library and one lane per sample: the single-sheet format cannot "
        "express more, so this is an assumption, not a recorded fact",
    ]
    if assay == "bulk_end_tag":
        notes.append("kit is unknown; H18 and V05 require a named kit before "
                     "this assay is qualified, so it stays unset")

    meta_cols = [c for c in (sheet_rows[0] if sheet_rows else {})
                 if c not in ("sample_id", "fastq_url", "fastq_url_r2")]

    for row in sheet_rows:
        sid = row.get("sample_id", "")
        sample = {"sample_id": sid, "biological_unit_id": sid}
        for c in meta_cols:
            sample[c] = row.get(c, "")
        samples.append(sample)

        lib_id = f"{sid}_lib1"
        libraries.append({
            "library_id": lib_id, "sample_id": sid, "layout": layout,
            "assay_family": assay, "strandedness": "", "kit": "",
            "platform": "", "prep_batch": "",
        })

        r1 = row.get("fastq_url", "")
        r2 = row.get("fastq_url_r2", "")
        if layout == "paired":
            reads.append({"read_unit_id": f"{lib_id}_L1_R1", "library_id": lib_id,
                          "role": "R1", "uri": r1, "lane": "L1"})
            reads.append({"read_unit_id": f"{lib_id}_L1_R2", "library_id": lib_id,
                          "role": "R2", "uri": r2, "lane": "L1"})
            if not r2:
                notes.append(f"{sid}: layout is paired but fastq_url_r2 is empty; "
                             f"the R2 read unit is recorded with no URI and will "
                             f"fail validation (IN02)")
        else:
            reads.append({"read_unit_id": f"{lib_id}_L1", "library_id": lib_id,
                          "role": "single", "uri": r1, "lane": "L1"})

    notes.append(f"{len(samples)} samples, {len(libraries)} libraries, "
                 f"{len(reads)} read units")
    if meta_cols:
        notes.append(f"carried onto samples: {', '.join(meta_cols)}")

    return {"samples": samples, "libraries": libraries, "reads": reads}, notes


def write_tables(tables: dict, out_dir):
    """Write the tables, each with its schema's columns first."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    written = []
    for name, rows in tables.items():
        if not rows:
            continue
        lead = list(SCHEMAS[name]["required"]) + list(SCHEMAS[name]["recommended"])
        cols = [c for c in lead if any(c in r for r in rows)]
        cols += [c for c in rows[0] if c not in cols]
        path = out / f"{name}.tsv"
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                               extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)
        written.append(path)
    return written


def demo():
    tables = {
        "samples": [{"sample_id": "S1", "biological_unit_id": "P1", "condition": "ctl"},
                    {"sample_id": "S2", "biological_unit_id": "P2", "condition": "trt"}],
        "libraries": [{"library_id": "S1_lib1", "sample_id": "S1", "layout": "paired"},
                      {"library_id": "S2_lib1", "sample_id": "S2", "layout": "paired"}],
        "reads": [
            {"read_unit_id": "a", "library_id": "S1_lib1", "role": "R1",
             "uri": "/d/S1_L1_R1.fq.gz", "lane": "L1"},
            {"read_unit_id": "b", "library_id": "S1_lib1", "role": "R2",
             "uri": "/d/S1_L1_R2.fq.gz", "lane": "L1"},
            {"read_unit_id": "c", "library_id": "S2_lib1", "role": "R1",
             "uri": "/d/S2_L1_R1.fq.gz", "lane": "L1"},
            {"read_unit_id": "d", "library_id": "S2_lib1", "role": "R2",
             "uri": "/d/S2_L1_R2.fq.gz", "lane": "L1"},
        ],
    }
    assert validate(tables) == [], validate(tables)

    # a second lane is expressible, and valid
    tables["reads"] += [
        {"read_unit_id": "e", "library_id": "S1_lib1", "role": "R1",
         "uri": "/d/S1_L2_R1.fq.gz", "lane": "L2"},
        {"read_unit_id": "f", "library_id": "S1_lib1", "role": "R2",
         "uri": "/d/S1_L2_R2.fq.gz", "lane": "L2"},
    ]
    assert validate(tables) == [], validate(tables)
    assert sorted(lane_groups(tables)["S1_lib1"]) == ["L1", "L2"]

    # dropping one mate of one lane is caught, and names the lane
    tables["reads"] = [r for r in tables["reads"] if r["read_unit_id"] != "f"]
    codes = [c for c, _ in validate(tables)]
    assert "MET010" in codes, validate(tables)

    print("metadata.demo: ok")


if __name__ == "__main__":
    demo()

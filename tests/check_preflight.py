#!/usr/bin/env python3
"""Self-check for scripts/preflight.py (upgrade plan Phase 1).

The point of preflight is that an unrunnable project is rejected by name in
seconds, instead of after every FASTQ has been quantified. So these cases assert
on the issue CODE, not just on a non-zero exit:

  - unknown / mistyped config keys
  - an assay or capability outside config/spec.yaml
  - sample-sheet violations (duplicate ids, missing model columns, one level)
  - designs that cannot be fitted: rank-deficient, no residual df, n < 2
  - the design family and backend chosen from the config
  - the emitted plan naming which stages will and will not run, and why

Acceptance scenario 5 of the upgrade plan ("an unestimable interaction fails
during preflight, before expensive analysis") is what the DSN cases cover.

Run:  python3 tests/check_preflight.py
"""
import json
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PREFLIGHT = ROOT / "scripts" / "preflight.py"
TMP = Path(tempfile.mkdtemp())

BASE_CONFIG = """\
project: {{name: t, description: d, output_dir: "results/"}}
organism: {{common_name: grape, scientific_name: Vitis vinifera, tax_id: 29760,
  kegg_code: vvi, orgdb_package: org.Vvinifera.eg.db}}
reference: {{accession: GCF_1.1, assembly_name: ASM1}}
samples: {{sheet: config/samples.tsv, seq_type: {seq_type}}}
model: {{fixed_effects: "{fixed}", random_effects: {random}, primary_factor: {primary}}}
downstream: {{run_go: true, run_kegg: true, run_wgcna: false}}
orgdb: {{strategy: auto}}
hpc: {{storage_budget_gb: 20, delete_fastq_after_quant: true, samples_in_flight: null}}
{extra}"""


def project(name, sheet, *, seq_type="tagseq", fixed="~ treatment",
            random="null", primary="treatment", extra=""):
    p = TMP / name
    (p / "config").mkdir(parents=True, exist_ok=True)
    (p / "config" / "config.yaml").write_text(BASE_CONFIG.format(
        seq_type=seq_type, fixed=fixed, random=random, primary=primary,
        extra=extra))
    (p / "config" / "samples.tsv").write_text(sheet)
    return p


def run(p):
    r = subprocess.run([sys.executable, str(PREFLIGHT), "-d", str(p)],
                       capture_output=True, text=True)
    out = r.stdout + r.stderr
    codes = set()
    for line in out.splitlines():
        for tok in line.split():
            if len(tok) == 6 and tok[:3].isalpha() and tok[3:].isdigit():
                codes.add(tok)
    plan_path = p / "gates" / "preflight_plan.json"
    plan = json.loads(plan_path.read_text()) if plan_path.is_file() else None
    return r.returncode, codes, plan, out


def expect(label, p, *, rc, has=(), lacks=()):
    got_rc, codes, plan, out = run(p)
    if got_rc != rc:
        sys.exit(f"FAIL {label}: expected exit {rc}, got {got_rc}\n{out}")
    for c in has:
        if c not in codes:
            sys.exit(f"FAIL {label}: expected issue {c}, got {sorted(codes)}\n{out}")
    for c in lacks:
        if c in codes:
            sys.exit(f"FAIL {label}: did not expect {c}\n{out}")
    return plan, codes, out


def sheet(rows, header="sample_id\ttreatment"):
    return header + "\n" + "".join(r + "\n" for r in rows)


BALANCED = sheet(["C1\tcontrol", "C2\tcontrol", "C3\tcontrol",
                  "T1\ttreated", "T2\ttreated", "T3\ttreated"])

# --- 1. a healthy project passes, with the design actually checked -----
plan, codes, _ = expect("healthy", project("ok", BALANCED), rc=0)
assert plan["design"]["checked"] and plan["design"]["estimable"], plan["design"]
assert plan["design"]["n_bio_replicates"] == 3, plan["design"]
assert plan["design"]["residual_df"] == 4, plan["design"]
assert plan["resolved"]["design"] == "independent_groups"
assert plan["resolved"]["backend"] == "limma_voom"
# TAGseq must not be length-corrected
assert plan["resolved"]["count_treatment"] == "no", plan["resolved"]
planned = {s["stage"] for s in plan["stages"] if s["planned"]}
assert "differential_expression" in planned and "enrichment" in planned
notplanned = {s["stage"]: s["reason"] for s in plan["stages"] if not s["planned"]}
assert notplanned == {"wgcna": "downstream.run_wgcna is false"}, notplanned

# --- 2. unknown and mistyped config keys ------------------------------
expect("unknown key", project("badkey", BALANCED, extra="bogus_section: {a: 1}\n"),
       rc=1, has=["CFG001"])
# storage_budget_gb now defaults to null (no default cap), and a null default
# constrains no type, so the type case uses a key that still commits to one.
expect("wrong type",
       project("badtype", BALANCED,
               extra="samples: {sheet: config/samples.tsv, seq_type: tagseq, "
                     "min_reads_on_genes: \"many\"}\n"),
       rc=1, has=["CFG003"])

# --- 3. scope: an assay the pipeline does not implement ---------------
expect("srnaseq", project("srna", BALANCED, seq_type="srnaseq"),
       rc=1, has=["SCP001"])
expect("nonsense assay", project("nonsense", BALANCED, seq_type="nanopore_direct"),
       rc=1, has=["SCP001"])
# bulk paired is in scope and carries no kit warning
_, codes, _ = expect("bulk paired", project("bulk", BALANCED, seq_type="rnaseq_paired"),
                     rc=0, lacks=["SCP004"])

# --- 4. sample-sheet contract ----------------------------------------
expect("duplicate id",
       project("dupe", sheet(["C1\tcontrol", "C1\tcontrol", "C2\tcontrol",
                              "T1\ttreated", "T2\ttreated", "T3\ttreated"])),
       rc=1, has=["MET002"])
expect("model column absent",
       project("nocol", BALANCED, fixed="~ genotype", primary="genotype"),
       rc=1, has=["MET004"])
expect("single level",
       project("onelevel", sheet(["C1\tcontrol", "C2\tcontrol", "C3\tcontrol"])),
       rc=1, has=["MET005"])
expect("random unit absent",
       project("norandom", BALANCED, random='"(1|vine)"'),
       rc=1, has=["DSN005"])

# --- 5. designs that cannot be fitted (acceptance scenario 5) ---------
# treatment and batch move together, so their effects cannot be separated.
CONFOUNDED = sheet(
    ["C1\tcontrol\tb1", "C2\tcontrol\tb1", "C3\tcontrol\tb1",
     "T1\ttreated\tb2", "T2\ttreated\tb2", "T3\ttreated\tb2"],
    header="sample_id\ttreatment\tbatch")
expect("confounded batch",
       project("confounded", CONFOUNDED, fixed="~ treatment + batch"),
       rc=1, has=["DSN001"])

# validate_design checks in order: rank, then residual df, then replication.
# Each case below is built to reach a different one of those three.

# 2 samples, 2 coefficients: residual df is 0, so this stops at DSN002 before
# the replicate count is ever consulted.
expect("no residual df",
       project("norep", sheet(["C1\tcontrol", "T1\ttreated"])),
       rc=1, has=["DSN002"])

# 3 samples, 2 coefficients: 1 residual df, so it gets past DSN002, but the
# smallest group has a single sample.
expect("one sample in a group",
       project("onerep", sheet(["C1\tcontrol", "C2\tcontrol", "T1\ttreated"])),
       rc=1, has=["DSN003"])

# 3 samples, 3 levels: saturated, no residual df.
expect("saturated",
       project("saturated", sheet(["S1\ta", "S2\tb", "S3\tc"])),
       rc=1, has=["DSN002"])

# an unestimable interaction: no cell has both levels of each factor
CELLS = sheet(
    ["S1\tctl\tearly", "S2\tctl\tearly", "S3\ttrt\tlate", "S4\ttrt\tlate"],
    header="sample_id\ttreatment\tstage")
expect("unestimable interaction",
       project("interaction", CELLS, fixed="~ treatment * stage"),
       rc=1, has=["DSN001"])

# --- 6. repeated measures picks dream, and counts units not rows ------
PAIRED = sheet(
    ["S1\tctl\tv1", "S2\tctl\tv2", "S3\tctl\tv3",
     "S4\ttrt\tv1", "S5\ttrt\tv2", "S6\ttrt\tv3"],
    header="sample_id\ttreatment\tvine")
plan, _, _ = expect("repeated measures",
                    project("repeated", PAIRED, random='"(1|vine)"'), rc=0)
assert plan["resolved"]["design"] == "repeated_measures", plan["resolved"]
assert plan["resolved"]["backend"] == "dream", plan["resolved"]
assert plan["design"]["biological_unit"] == "vine", plan["design"]

# --- 7. factorial is recognised from the formula ----------------------
FACTORIAL = sheet(
    ["S1\tctl\tearly", "S2\tctl\tearly", "S3\tctl\tlate", "S4\tctl\tlate",
     "S5\ttrt\tearly", "S6\ttrt\tearly", "S7\ttrt\tlate", "S8\ttrt\tlate"],
    header="sample_id\ttreatment\tstage")
plan, _, _ = expect("factorial",
                    project("factorial", FACTORIAL, fixed="~ treatment * stage"),
                    rc=0)
assert plan["resolved"]["design"] == "factorial", plan["resolved"]
assert plan["design"]["estimable"], plan["design"]

# --- 8. a stage that cannot run says why ------------------------------
plan, _, _ = expect("enrichment off",
                    project("noenrich", BALANCED,
                            extra="downstream: {run_go: false, run_kegg: false, "
                                  "run_wgcna: false}\n"),
                    rc=0)
reasons = {s["stage"]: s["reason"] for s in plan["stages"] if not s["planned"]}
assert "enrichment" in reasons and "run_go" in reasons["enrichment"], reasons

# --- 9. the issue table is written and matches the plan counts --------
p = project("table", BALANCED, extra="bogus: 1\n")
rc, codes, plan, _ = run(p)
rows = (p / "gates" / "preflight_issues.tsv").read_text().splitlines()
assert rows[0].split("\t") == ["code", "severity", "scope", "message",
                               "detail", "remedy"], rows[0]
assert len(rows) - 1 == plan["n_errors"] + plan["n_warnings"], (len(rows), plan)

# --- 10. every code preflight can emit is declared in spec.yaml -------
import yaml
spec = yaml.safe_load((ROOT / "config" / "spec.yaml").read_text())
src = PREFLIGHT.read_text()
import re
emitted = set(re.findall(r'iss\.add\("([A-Z]{3}\d{3})"', src))
undeclared = emitted - set(spec["issues"])
assert not undeclared, f"codes emitted but not declared in spec.yaml: {undeclared}"

print(f"check_preflight.py: all assertions passed ({len(emitted)} issue codes exercised or declared)")

#!/usr/bin/env python3
"""preflight.py — validate a project before anything expensive runs.

Upgrade plan Phase 1. Checks the configuration, the sample sheet, the requested
assay/design/backend combination against config/spec.yaml, and the estimability
of the model, then writes an issue table and a machine-readable plan.

Exits 1 if any error-severity issue is found. Warnings do not block.

The design check matters most: without it, an unestimable model is only caught by
03_de.R, which runs after every FASTQ has been quantified. Here it costs seconds.

Usage:
    python3 scripts/preflight.py --project-dir ~/rnaseq_projects/my_study
    python3 scripts/preflight.py -d <project> --configfile config/config.yaml

Outputs, inside the project:
    gates/preflight_issues.tsv    code, severity, scope, message, detail, remedy
    gates/preflight_plan.json     counts, design facts, stages, resolved settings
"""
from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import yaml

import config_resolve
import metadata as metadata_tables
from recommend import recommend

ROOT = Path(__file__).resolve().parent.parent
SPEC = ROOT / "config" / "spec.yaml"
DESIGN_LIB = ROOT / "workflow" / "scripts" / "_design.R"


# ----------------------------------------------------------------- issue table
class Issues:
    """Accumulates findings. Codes and severities come from spec.yaml."""

    def __init__(self, catalog: dict):
        self.catalog = catalog
        self.rows: list[dict] = []

    def add(self, code: str, scope: str, detail: str = ""):
        entry = self.catalog.get(code)
        if entry is None:
            # An unknown code is itself a bug; surface it rather than dropping it.
            entry = {"severity": "error", "message": f"Undeclared issue code {code}"}
        self.rows.append({
            "code": code,
            "severity": entry["severity"],
            "scope": scope,
            "message": entry["message"],
            "detail": detail,
            "remedy": entry.get("remedy", ""),
        })

    @property
    def errors(self):
        return [r for r in self.rows if r["severity"] == "error"]

    @property
    def warnings(self):
        return [r for r in self.rows if r["severity"] == "warning"]

    def write(self, path: Path):
        cols = ["code", "severity", "scope", "message", "detail", "remedy"]
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as fh:
            fh.write("\t".join(cols) + "\n")
            for r in self.rows:
                fh.write("\t".join(str(r[c]).replace("\t", " ") for c in cols) + "\n")


# ------------------------------------------------------------------- config
def check_config(cfg: dict, repo_root, iss: Issues) -> dict:
    """Resolve the project config and report what resolution found.

    Shares scripts/config_resolve.py with the Snakefile and setup, so the three
    cannot disagree about defaults, unknown keys or types (R17e, R05; master
    plan 6.3). Returns the resolved configuration, which every later check uses
    in place of the raw file: validating the raw config would miss a default
    that only appears after merging.
    """
    resolved, origins, issues = config_resolve.resolve(
        repo_root, [("project", cfg)])
    for code, detail in issues:
        iss.add(code, "config", detail)

    for section in ("project", "samples", "model", "reference", "organism"):
        if section not in resolved:
            iss.add("CFG002", "config", section)

    return resolved


# -------------------------------------------------------------------- scope
def formula_vars(formula: str) -> list[str]:
    """Variable names in an R-style formula, right-hand side only."""
    rhs = formula.split("~", 1)[-1]
    seen, out = set(), []
    for tok in re.findall(r"[A-Za-z._][A-Za-z0-9._]*", rhs):
        if tok not in seen:
            seen.add(tok)
            out.append(tok)
    return out


def random_effect_unit(term: str | None) -> str | None:
    """The grouping column in a term like '(1|vine)'."""
    if not term:
        return None
    m = re.search(r"\|\s*([A-Za-z._][A-Za-z0-9._]*)", term)
    return m.group(1) if m else None


def classify_design(cfg: dict) -> str:
    """Which design family this config asks for."""
    if random_effect_unit(cfg.get("model", {}).get("random_effects")):
        return "repeated_measures"
    fixed = cfg.get("model", {}).get("fixed_effects", "")
    if "*" in fixed or ":" in fixed:
        return "factorial"
    return "independent_groups"


def check_scope(cfg: dict, spec: dict, iss: Issues) -> dict:
    """Assay, design and backend against the supported matrix."""
    scope = spec["scope"]
    assay = cfg.get("samples", {}).get("seq_type")
    design = classify_design(cfg)

    a = scope["assays"].get(assay)
    if a is None:
        iss.add("SCP001", "scope",
                f"{assay!r}; supported: {', '.join(sorted(scope['assays']))}")
    elif not a.get("supported"):
        iss.add("SCP001", "scope", f"{assay!r}: {a.get('reason', 'not supported')}")
    elif assay == "tagseq" and not a.get("validated_kits"):
        iss.add("SCP004", "scope",
                "tagseq: scope.assays.tagseq.validated_kits is empty, so the "
                "count treatment is unverified for your kit")

    d = scope["designs"].get(design)
    backend = None
    if d is None:
        iss.add("SCP003", "scope", f"design {design!r} is not in scope.designs")
    else:
        backend = d["backends"][0]
        needs_random = design == "repeated_measures"
        if scope["backends"][backend]["handles_random_effects"] != needs_random:
            iss.add("SCP002", "scope",
                    f"design {design} with backend {backend}")

    return {"assay": assay, "design": design, "backend": backend,
            "count_treatment": (a or {}).get("count_treatment")}


# ----------------------------------------------------------------- metadata
def read_sheet(path: Path) -> tuple[list[str], list[list[str]]]:
    header, rows = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = line.split("\t")
            if header is None:
                header = fields
            else:
                rows.append(fields)
    return header or [], rows


def check_metadata(cfg: dict, spec: dict, sheet_path: Path, iss: Issues) -> dict:
    contract = spec["metadata"]["samples_tsv"]
    header, rows = read_sheet(sheet_path)

    for col in contract["required_columns"]:
        if cfg.get("samples", {}).get("metadata_dir") and col in ("fastq_url", "fastq_url_r2"):
            continue
        if col not in header:
            iss.add("MET001", "samples.tsv", col)

    def column(name):
        if name not in header:
            return []
        i = header.index(name)
        return [r[i] if i < len(r) else "" for r in rows]

    for col in contract["unique_columns"]:
        vals = column(col)
        dupes = sorted({v for v in vals if vals.count(v) > 1})
        if dupes:
            iss.add("MET002", "samples.tsv", f"{col}: {', '.join(dupes)}")

    for col in contract["forbid_blank_in"]:
        blanks = [i + 2 for i, v in enumerate(column(col)) if not v.strip()]
        if blanks:
            iss.add("MET003", "samples.tsv",
                    f"{col}: row(s) {', '.join(map(str, blanks))}")

    model = cfg.get("model", {})
    wanted = formula_vars(model.get("fixed_effects", ""))
    unit = random_effect_unit(model.get("random_effects"))
    primary = model.get("primary_factor")
    if primary and primary not in wanted:
        wanted.append(primary)

    for var in wanted:
        if var not in header:
            iss.add("MET004", "model", f"{var!r} (fixed_effects / primary_factor)")
        elif len({v for v in column(var) if v.strip()}) < 2:
            iss.add("MET005", "model", f"{var!r} has one level in the sheet")

    if unit and unit not in header:
        iss.add("DSN005", "model", f"{unit!r} (random_effects)")

    declared = set(wanted) | ({unit} if unit else set())
    unused = [c for c in header
              if c not in declared and c not in ("sample_id", "fastq", "fastq_r2")]
    if unused:
        iss.add("MET006", "samples.tsv", ", ".join(unused))

    return {"n_rows": len(rows), "columns": header,
            "model_vars": wanted, "random_unit": unit}


# ------------------------------------------------------------------- design
R_DESIGN = r'''
args <- commandArgs(trailingOnly = TRUE)
source(args[1])
sheet <- read_sample_sheet(args[2])
fixed <- args[3]; primary <- args[4]
random <- if (length(args) >= 5 && nzchar(args[5])) args[5] else NULL

res <- tryCatch({
  reps <- biological_replicates(sheet, primary, random)
  mm <- model.matrix(as.formula(fixed), data = as.data.frame(sheet))
  d <- validate_design(mm, reps$n, primary, reps$unit)
  list(ok = TRUE, n_samples = nrow(mm), rank = d$rank,
       residual_df = d$residual_df, n_bio_replicates = reps$n,
       biological_unit = reps$unit)
}, error = function(e) list(ok = FALSE, error = conditionMessage(e)))

cat(jsonlite::toJSON(res, auto_unbox = TRUE))
'''


def check_metadata_tables(proj: Path, iss: Issues, cfg=None) -> dict | None:
    """Validate the three-table contract when a project uses it.

    Master plan 6.1 defines samples / libraries / reads. A project that still
    has only samples.tsv keeps working: this returns None and the caller falls
    back to the single-sheet checks, which the legacy form is all that supports.
    """
    mdir = proj / ((cfg or {}).get("samples", {}).get("metadata_dir") or "metadata")
    try:
        tables = metadata_tables.load_tables(mdir)
    except (OSError, ValueError) as exc:
        iss.add("MET001", "metadata", str(exc))
        return None
    if not any((mdir / f"{name}.tsv").exists() for name in metadata_tables.SCHEMAS):
        return None

    present = [n for n, rows in tables.items() if rows]
    for code, detail in metadata_tables.validate(tables):
        iss.add(code, "metadata", detail)

    lanes = metadata_tables.lane_groups(tables)
    multi_lane = {lib: [{"run": g["run"], "lane": g["lane"]} for g in groups]
                  for lib, groups in lanes.items() if len(groups) > 1}
    return {
        "form": "three_table",
        "tables_present": present,
        "counts": {n: len(rows) for n, rows in tables.items()},
        "multi_lane_libraries": multi_lane,
    }


def check_design(cfg: dict, sheet_path: Path, iss: Issues) -> dict:
    """Estimability, via the same _design.R the DE stage uses.

    Calling into R rather than reimplementing model.matrix and the QR rank check
    in Python: two implementations would eventually disagree, and disagreeing
    with the stage that actually fits the model is the one outcome worse than
    having no preflight check.
    """
    if shutil.which("Rscript") is None:
        iss.add("ENV001", "design", "Rscript (needed for the estimability check)")
        return {"checked": False, "reason": "Rscript not found"}

    model = cfg.get("model", {})
    proc = subprocess.run(
        # No "--args": with -e, Rscript passes it through as a trailing
        # argument, so it would land in args[1] where the library path belongs.
        ["Rscript", "-e", R_DESIGN, str(DESIGN_LIB), str(sheet_path),
         model.get("fixed_effects", "~ 1"), model.get("primary_factor", ""),
         model.get("random_effects") or ""],
        capture_output=True, text=True)

    payload = proc.stdout.strip()
    start = payload.find("{")
    if start < 0:
        iss.add("ENV001", "design",
                f"estimability check produced no result: "
                f"{(proc.stderr or payload).strip()[:300]}")
        return {"checked": False, "reason": "no JSON from Rscript"}

    res = json.loads(payload[start:])
    if not res.get("ok"):
        msg = res.get("error", "unknown error")
        low = msg.lower()
        code = ("DSN001" if "rank-deficient" in low
                else "DSN002" if "residual degrees" in low
                else "DSN003" if "biological replicate" in low
                else "DSN004")
        iss.add(code, "design", msg)
        return {"checked": True, "estimable": False, "error": msg}

    res["checked"] = True
    res["estimable"] = True
    return res


# --------------------------------------------------------------------- main
def stage_plan(cfg: dict, resolved: dict, design: dict) -> list[dict]:
    """Which stages will run, and why any will not."""
    down = cfg.get("downstream", {}) or {}
    orgdb_strategy = (cfg.get("orgdb", {}) or {}).get("strategy", "auto")
    estimable = design.get("estimable", False)

    def stage(name, run, reason=""):
        return {"stage": name, "planned": bool(run), "reason": reason}

    def enrichment_plan(down, orgdb_strategy, estimable):
        """(planned, reason). Enrichment consumes DE output, so DE blocks it."""
        if not (down.get("run_go") or down.get("run_kegg")):
            return False, "downstream.run_go and run_kegg are both false"
        if orgdb_strategy == "skip":
            return False, "orgdb.strategy=skip"
        if not estimable:
            return False, "depends on differential_expression, which is not planned"
        return True, ""

    return [
        stage("reference", True),
        stage("quantification", True),
        stage("aggregation", True),
        stage("differential_expression", estimable,
              "" if estimable else "design is not estimable"),
        stage("enrichment", *enrichment_plan(down, orgdb_strategy, estimable)),
        stage("wgcna", bool(down.get("run_wgcna")),
              "" if down.get("run_wgcna") else "downstream.run_wgcna is false"),
        stage("report", True),
    ]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-d", "--project-dir", default=".",
                    help="project directory (default: current directory)")
    ap.add_argument("-c", "--configfile", default="config/config.yaml",
                    help="config path, relative to the project directory")
    ap.add_argument("--quiet", action="store_true",
                    help="only print the summary line and any errors")
    args = ap.parse_args()

    proj = Path(args.project_dir).resolve()
    if not proj.is_dir():
        sys.exit(f"preflight: project directory not found: {proj}")
    cfg_path = proj / args.configfile
    if not cfg_path.is_file():
        sys.exit(f"preflight: config not found: {cfg_path}")

    spec = yaml.safe_load(open(SPEC))
    cfg = yaml.safe_load(open(cfg_path))
    iss = Issues(spec["issues"])

    # Every later check runs against the RESOLVED config, not the raw file:
    # a default that only exists after merging would otherwise look missing.
    cfg = check_config(cfg, ROOT, iss)
    resolved = check_scope(cfg, spec, iss)

    sheet_rel = (cfg.get("samples", {}) or {}).get("sheet", "config/samples.tsv")
    sheet_path = proj / sheet_rel
    meta: dict = {}
    design: dict = {"checked": False, "reason": "sample sheet unreadable"}
    # Three-table form takes precedence when present (master plan 6.1); the
    # single sheet remains supported so existing projects still run.
    tables_meta = check_metadata_tables(proj, iss, cfg)
    try:
        metadata_tables.execution_inputs(proj, cfg)
    except (OSError, ValueError) as exc:
        iss.add("MET013", "execution inputs", str(exc))

    if not sheet_path.is_file():
        iss.add("MET001", "samples.tsv", f"file not found: {sheet_path}")
    else:
        meta = check_metadata(cfg, spec, sheet_path, iss)
        # Only worth asking R about estimability once the columns exist.
        if not [r for r in iss.rows
                if r["code"] in ("MET004", "MET005", "DSN005")]:
            design = check_design(cfg, sheet_path, iss)
        else:
            design = {"checked": False,
                      "reason": "model variables missing or single-level"}

    recommendation = recommend(cfg, meta.get("n_rows"))
    for reason in recommendation["issues"]:
        iss.add("SCP003", "recommendation", reason)

    plan = {
        "spec_version": spec["spec_version"],
        "project_dir": str(proj),
        "config": str(cfg_path),
        "resolved": resolved,
        "recommendation": recommendation,
        "metadata": meta,
        "metadata_tables": tables_meta,
        "design": design,
        "stages": stage_plan(cfg, resolved, design),
        "n_errors": len(iss.errors),
        "n_warnings": len(iss.warnings),
    }

    issues_path = proj / "gates" / "preflight_issues.tsv"
    plan_path = proj / "gates" / "preflight_plan.json"
    iss.write(issues_path)
    plan_path.parent.mkdir(parents=True, exist_ok=True)
    (plan_path.parent / "recommendation.json").write_text(json.dumps(recommendation, indent=2) + "\n")
    with open(plan_path, "w") as fh:
        json.dump(plan, fh, indent=2)
        fh.write("\n")

    if not args.quiet:
        print(f"preflight: {proj}")
        print(f"  assay={resolved['assay']}  design={resolved['design']}  "
              f"backend={resolved['backend']}  "
              f"count_treatment={resolved['count_treatment']}")
        if meta:
            print(f"  samples={meta['n_rows']}  model_vars="
                  f"{','.join(meta['model_vars']) or '-'}  "
                  f"unit={meta['random_unit'] or 'sample'}")
        if design.get("checked") and design.get("estimable"):
            print(f"  design: {design['rank']} coefficients, "
                  f"{design['residual_df']} residual df, "
                  f"{design['n_bio_replicates']} biological replicates "
                  f"per {design['biological_unit']}")
        elif not design.get("checked"):
            print(f"  design: not checked ({design.get('reason')})")
        print("  stages: " + ", ".join(
            s["stage"] for s in plan["stages"] if s["planned"]))
        skipped = [s for s in plan["stages"] if not s["planned"]]
        for s in skipped:
            print(f"    not planned: {s['stage']} ({s['reason']})")

    for r in iss.rows:
        stream = sys.stderr if r["severity"] == "error" else sys.stdout
        line = f"  {r['severity'].upper():7} {r['code']}  {r['scope']}: {r['message']}"
        if r["detail"]:
            line += f" -> {r['detail']}"
        print(line, file=stream)
        if r["remedy"] and r["severity"] == "error":
            print(f"          {r['remedy']}", file=stream)

    print(f"preflight: {len(iss.errors)} error(s), {len(iss.warnings)} warning(s); "
          f"details in {issues_path.relative_to(proj)}")
    sys.exit(1 if iss.errors else 0)


if __name__ == "__main__":
    main()

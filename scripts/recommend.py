#!/usr/bin/env python3
"""Explain the current supported route; never substitute bulk for another assay."""
import hashlib
from pathlib import Path
import yaml

POLICY = Path(__file__).resolve().parent.parent / "config/recommendation_rules.yaml"


def recommend(cfg, n_samples=None):
    policy = yaml.safe_load(POLICY.read_text())
    seq = cfg.get("samples", {}).get("seq_type")
    analysis = cfg.get("analysis") or {}
    objectives = analysis.get("objectives", ["gene_expression", "differential_expression"])
    issues = []
    route = policy["routes"].get(seq)
    if route is None:
        issues.append(f"No implemented route for assay/layout {seq!r}")
    if not isinstance(objectives, list) or not objectives or not all(isinstance(x, str) for x in objectives):
        issues.append("analysis.objectives must be a nonempty list of objective names")
        objectives = []
    downstream = cfg.get("downstream") or {}
    if "enrichment" in objectives and not (downstream.get("run_go") or downstream.get("run_kegg")):
        issues.append("Enrichment objective requires downstream.run_go or downstream.run_kegg")
    if "coexpression" in objectives and not downstream.get("run_wgcna"):
        issues.append("Coexpression objective requires downstream.run_wgcna")
    unknown = sorted(set(objectives) - set(policy["objectives"]))
    if unknown:
        issues.append(f"Requested objectives have no implemented route: {', '.join(unknown)}")
    if analysis.get("umi"):
        issues.append("UMI extraction/deduplication has no qualified implementation")
    backend = policy["backends"]["random" if cfg.get("model", {}).get("random_effects") else "fixed"]
    requested = analysis.get("backend", "auto")
    if requested not in ("auto", backend["name"]):
        issues.append(f"Backend {requested!r} cannot execute the declared design in this release; requires {backend['name']}")
    return {
        "policy_version": policy["policy_version"],
        "policy_sha256": hashlib.sha256(POLICY.read_bytes()).hexdigest(),
        "status": "blocked" if issues else "selected",
        "route": route["route"] if route else None,
        "backend": backend["name"],
        "count_treatment": route["count_treatment"] if route else None,
        "rule_ids": (route["rule_ids"] if route else []) + backend["rule_ids"] + ["ST11"],
        "handbook_ids": route["handbook_ids"] if route else [],
        "inputs": {"seq_type": seq, "n_samples": n_samples, "objectives": objectives,
                   "fixed_effects": cfg.get("model", {}).get("fixed_effects"),
                   "random_effects": cfg.get("model", {}).get("random_effects")},
        "reasons": ["Annotated gene-expression objectives use the implemented Salmon route.",
                    "The declared dependence structure selects the statistical backend.",
                    "Sample count informs replication/QC checks; it does not silently change the model."],
        "qualification": "provisional; no named kit qualified" if seq == "tagseq" else "see release validation evidence",
        "alternatives": policy["alternatives"], "issues": issues,
    }

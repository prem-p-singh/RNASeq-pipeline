#!/usr/bin/env python3
"""Explain the current supported route; never substitute bulk for another assay."""
import hashlib
from pathlib import Path
import yaml
import dual_organism
from small_rna import load_kit
import tagseq

POLICY = Path(__file__).resolve().parent.parent / "config/recommendation_rules.yaml"


def requires_inference(cfg):
    objectives = (cfg.get("analysis") or {}).get("objectives", ["gene_expression", "differential_expression"])
    return isinstance(objectives, list) and bool(set(x for x in objectives if isinstance(x, str)) & {"differential_expression", "enrichment"})


def qualification(seq, analysis):
    """Kit status from its table: validated only after the kit's own public-data check passed."""
    if seq not in ("tagseq", "srnaseq"):
        return "see release validation evidence"
    kit = analysis.get("kit" if seq == "tagseq" else "small_rna_kit")
    if not kit:
        return "provisional; no named kit qualified"
    try:
        status = (tagseq.load_kit(kit) if seq == "tagseq" else load_kit(kit))["status"]
    except ValueError:
        return "blocked: unknown or unsupported kit"
    return "validated; public-data check passed" if status == "validated" else "experimental; no public-data check passed"


def handbook_backend(min_units):
    """Handbook default for a new independent fixed-effect design (ST01/ST02).

    Applied once at project setup and written into the config, so existing
    projects keep their declared backend. Two units per group is exploratory
    and must be chosen explicitly; fewer cannot be tested.
    """
    if min_units is None or min_units < 3:
        return None
    return "deseq2" if min_units <= 12 else "limma_voom"


def recommend(cfg, n_samples=None):
    policy = yaml.safe_load(POLICY.read_text())
    seq = cfg.get("samples", {}).get("seq_type")
    analysis = cfg.get("analysis") or {}
    objectives = analysis.get("objectives", ["gene_expression", "differential_expression"])
    issues = []
    if (cfg.get("dual_organism") or {}).get("enabled"):
        # Reads are assigned to an organism first; each organism then gets its own route.
        try:
            dual_organism.validate(cfg)
        except (ValueError, TypeError) as exc:
            issues.append(str(exc))
        # The shared checks that still apply before any per-organism route is chosen.
        if seq not in ("rnaseq_single", "rnaseq_paired", "tagseq", "srnaseq"):
            issues.append(f"No implemented route for assay/layout {seq!r}")
        if analysis.get("umi"):
            issues.append("UMI extraction/deduplication has no qualified implementation")
        if analysis.get("kit"):
            issues.append("analysis.kit applies to seq_type tagseq only")
        if not isinstance(objectives, list) or not objectives or not all(isinstance(x, str) for x in objectives):
            issues.append("analysis.objectives must be a nonempty list of objective names")
        elif set(objectives) - set(policy["objectives"]):
            issues.append("Requested objectives have no implemented route: "
                          + ", ".join(sorted(set(objectives) - set(policy["objectives"]))))
        dual = policy["dual_route"]
        return {
            "policy_version": policy["policy_version"],
            "policy_sha256": hashlib.sha256(POLICY.read_bytes()).hexdigest(),
            "status": "blocked" if issues else "selected",
            "route": dual["route"], "backend": None, "count_treatment": None,
            "rule_ids": dual["rule_ids"], "handbook_ids": dual["handbook_ids"],
            "inputs": {"seq_type": seq, "n_samples": n_samples, "objectives": objectives,
                       "fixed_effects": cfg.get("model", {}).get("fixed_effects"),
                       "random_effects": cfg.get("model", {}).get("random_effects")},
            "reasons": ["dual_organism.enabled: STAR assigns each read to host or pathogen on the combined genome."],
            "qualification": "experimental; no public-data check passed",
            "alternatives": policy["alternatives"], "issues": issues,
        }
    route = policy["routes"].get(seq)
    if seq == "srnaseq":
        route = policy["small_rna_routes"].get(analysis.get("small_rna_profile"))
        if route is None:
            issues.append(f"analysis.small_rna_profile must be one of {sorted(policy['small_rna_routes'])}")
        ref, org = cfg.get("reference") or {}, cfg.get("organism") or {}
        known = ((ref.get("mature_mirna_url"), "reference.mature_mirna_url"),
                 (org.get("mirbase_prefix"), "organism.mirbase_prefix"))
        # Plant loci come from the genome; known miRNAs only annotate them, so they are optional as a pair.
        if analysis.get("small_rna_profile") == "plant":
            known = ((ref.get("genome_fasta_url"), "reference.genome_fasta_url"),) + (known if any(v for v, _ in known) else ())
        for value, name in ((analysis.get("small_rna_kit"), "analysis.small_rna_kit"),) + known:
            if not value:
                issues.append(f"Small-RNA routes require {name}")
        if analysis.get("small_rna_kit"):
            try:
                load_kit(analysis["small_rna_kit"])
            except ValueError as exc:
                issues.append(str(exc))
        if isinstance(objectives, list) and "enrichment" in objectives:
            issues.append("GO/KEGG enrichment needs gene identifiers; small-RNA counts are miRNA/locus identifiers")
    elif route is None:
        issues.append(f"No implemented route for assay/layout {seq!r}")
    if not isinstance(objectives, list) or not objectives or not all(isinstance(x, str) for x in objectives):
        issues.append("analysis.objectives must be a nonempty list of objective names")
        objectives = []
    qc_only = bool(objectives) and set(objectives) == {"qc"}
    inference = requires_inference(cfg)
    quantifier = analysis.get("quantifier", "auto")
    if quantifier == "auto":
        quantifier = "star" if "alignment" in objectives else "salmon"
    kit = None
    if analysis.get("kit"):
        try:
            kit = tagseq.load_kit(analysis["kit"])
        except ValueError as exc:
            issues.append(str(exc))
        if seq != "tagseq":
            issues.append("analysis.kit applies to seq_type tagseq only")
        if not all((cfg.get("reference") or {}).get(k) for k in ("genome_fasta_url", "gtf_url")):
            issues.append("TAG-seq kit route requires reference.genome_fasta_url and gtf_url")
        if quantifier != "star" and analysis.get("quantifier", "auto") != "auto":
            issues.append("Named TAG-seq kits use the STAR route; set analysis.quantifier to auto or star")
    tagseq_route = seq == "tagseq" and kit is not None and not qc_only
    if tagseq_route:
        quantifier = "star"  # strand comes from the kit, not expected_libtype
    if quantifier not in ("salmon", "star"):
        issues.append("analysis.quantifier must be auto, salmon or star")
    if "alignment" in objectives and quantifier != "star":
        issues.append("Alignment output requires the STAR route")
    if quantifier == "star" and not qc_only and not tagseq_route:
        if seq not in ("rnaseq_single", "rnaseq_paired"):
            issues.append("STAR counting is implemented for non-UMI bulk libraries only")
        if not (cfg.get("reference") or {}).get("genome_fasta_url"):
            issues.append("STAR requires reference.genome_fasta_url")
        strand = (cfg.get("samples") or {}).get("expected_libtype")
        allowed_strands = ("IU", "ISF", "ISR") if seq == "rnaseq_paired" else ("U", "SF", "SR")
        if not (cfg.get("samples") or {}).get("metadata_dir") and strand not in allowed_strands:
            issues.append("STAR requires a declared expected_libtype matching the read layout; strand cannot be guessed")
    downstream = cfg.get("downstream") or {}
    if "enrichment" in objectives and not (downstream.get("run_go") or downstream.get("run_kegg")):
        issues.append("Enrichment objective requires downstream.run_go or downstream.run_kegg")
    if "coexpression" in objectives and not downstream.get("run_wgcna"):
        issues.append("Coexpression objective requires downstream.run_wgcna")
    unknown = sorted(set(objectives) - set(policy["objectives"]))
    if unknown:
        issues.append(f"Requested objectives have no implemented route: {', '.join(unknown)}")
    if analysis.get("umi") and not (kit and kit["umi_len"]):
        issues.append("UMI extraction/deduplication has no qualified implementation")
    backend = policy["backends"]["random" if cfg.get("model", {}).get("random_effects") else "fixed"]
    requested = analysis.get("backend", "auto")
    if not cfg.get("model", {}).get("random_effects") and requested in policy["fixed_backends"]:
        backend = {"name": requested, "rule_ids": ["ST01", "ST02"]}
    if inference and requested not in ("auto", backend["name"]):
        issues.append(f"Backend {requested!r} cannot execute the declared design in this release; requires {backend['name']}")
    return {
        "policy_version": policy["policy_version"],
        "policy_sha256": hashlib.sha256(POLICY.read_bytes()).hexdigest(),
        "status": "blocked" if issues else "selected",
        "route": ("fastp_qc" if qc_only else "star_tagseq" if tagseq_route
                  else "star_counts" if quantifier == "star" else route["route"]) if route else None,
        "backend": backend["name"] if inference else None,
        "count_treatment": (("umi_gene_counts" if kit["umi_len"] else "raw_gene_counts") if tagseq_route
                            else "raw_gene_counts" if quantifier == "star" else route["count_treatment"]) if route and not qc_only else None,
        "rule_ids": (route["rule_ids"] if route else []) + backend["rule_ids"] + ["ST11"],
        "handbook_ids": route["handbook_ids"] if route else [],
        "inputs": {"seq_type": seq, "n_samples": n_samples, "objectives": objectives,
                   "fixed_effects": cfg.get("model", {}).get("fixed_effects"),
                   "random_effects": cfg.get("model", {}).get("random_effects")},
        "reasons": ["QC-only preprocessing needs no quantifier, reference or statistical model." if qc_only else
                    f"Selected {('ShortStack' if analysis.get('small_rna_profile') == 'plant' else 'bowtie') if seq == 'srnaseq' else quantifier} for the declared expression/alignment objectives.",
                    "The declared dependence structure selects the statistical backend." if inference else
                    "No differential-inference objective was requested; no DE model will be fitted.",
                    "Sample count informs replication/QC checks; it does not silently change the model."],
        "qualification": qualification(seq, analysis),
        "alternatives": policy["alternatives"], "issues": issues,
    }

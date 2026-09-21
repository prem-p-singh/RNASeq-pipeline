# =============================================================================
# Stages 2-5 — runs ONCE after all per-sample quants finish.
# Memory-bound not storage-bound. Safe to run on login node or small partition.
#
# Scripts referenced here are parameterized refactors of the project-specific
# .Rmd files (analysis_20{21,22}.Rmd, GO_*.Rmd, KEGG_*.Rmd, WGCNA.Rmd).
# They are written next — these rules currently wire up the DAG.
# =============================================================================

# --- Stage 2: aggregate quants to a count matrix -----------------------
rule aggregate_counts:
    input:
        quants = expand(str(QUANT / "{sample}/quant.sf"), sample=SAMPLES),
        # Per-sample QC drives the inclusion decision, so declare it: a changed
        # metrics.json must re-run aggregation, not be silently ignored.
        sample_metrics = expand(str(QUANT / "{sample}/metrics.json"), sample=SAMPLES),
        gtf = REF / "annotation.gtf",
        sheet = config["samples"]["sheet"],
    output:
        counts = OUT / "counts.tsv",
        metrics = METRICS / "aggregate.json",
        disposition = OUT / "sample_disposition.tsv",
    resources:
        mem_mb = 8000,
        runtime = 30,
    script:
        "../scripts/02_aggregate.R"


# --- Stage 3: differential expression (adaptive backend) --------------
rule differential_expression:
    input:
        counts = OUT / "counts.tsv",
        metrics = METRICS / "aggregate.json",
        sheet = config["samples"]["sheet"],
    output:
        done = touch(OUT / "de_done.flag"),
        metrics = METRICS / "de.json",
        manifest = OUT / "de_manifest.tsv",
    params:
        model_fixed = config["model"]["fixed_effects"],
        model_random = config["model"].get("random_effects"),
        primary = config["model"]["primary_factor"],
    resources:
        mem_mb = 16000,
        runtime = 120,
    script:
        "../scripts/03_de.R"


# --- Stage 4: enrichment (GO + KEGG, with simplify) -------------------
# Conditional dependency: only wait for the OrgDb sentinel when enrichment
# is actually requested (strategy != "skip"). If skipped, the rule still
# runs but 04_enrichment.R will see orgdb missing and emit skip metrics.
def enrichment_inputs(wildcards):
    # The manifest, not a directory glob, defines which contrasts exist.
    deps = {
        "de_flag": str(OUT / "de_done.flag"),
        "manifest": str(OUT / "de_manifest.tsv"),
    }
    if orgdb_required():
        deps["orgdb_sentinel"] = orgdb_sentinel()
    return deps


rule enrichment:
    input:
        unpack(enrichment_inputs),
    output:
        done = touch(OUT / "enrichment_done.flag"),
        metrics = METRICS / "enrichment.json",
        status = OUT / "enrichment_status.tsv",
    params:
        orgdb = config["organism"]["orgdb_package"],
        kegg_code = config["organism"].get("kegg_code"),
        run_go = config["downstream"]["run_go"],
        run_kegg = config["downstream"]["run_kegg"],
    resources:
        mem_mb = 8000,
        runtime = 60,
    script:
        "../scripts/04_enrichment.R"


# --- Stage 5: WGCNA (auto-skip if n_samples too small) ----------------
rule wgcna:
    input:
        counts = OUT / "counts.tsv",
        de_flag = OUT / "de_done.flag",
    output:
        done = touch(OUT / "wgcna_done.flag"),
        metrics = METRICS / "wgcna.json",
    params:
        run_wgcna = config["downstream"]["run_wgcna"],
        min_samples = config["thresholds"]["wgcna"]["min_samples"],
    # The script reads both: nThreads is passed to WGCNA instead of a hardcoded
    # 2 (which asked SLURM for 1 CPU and then used two), and mem_mb sizes
    # maxBlockSize via WGCNA::blockSize.
    threads: 4
    resources:
        mem_mb = 16000,
        runtime = 180,
    script:
        "../scripts/05_wgcna.R"

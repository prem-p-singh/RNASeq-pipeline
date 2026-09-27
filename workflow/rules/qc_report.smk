# =============================================================================
# QC report — auto-generated comparative charts + MultiQC.
# Runs once after all samples are quantified. Produces:
#   results/qc_report/comparative_charts.png  (before/after clean + alignment fate)
#   results/qc_report/qc_charts.html          (self-contained report)
#   results/qc_report/alignment_summary.tsv   (per-sample table)
#   results/qc_report/multiqc/multiqc_report.html
# Depends only on per-sample outputs, so it runs right after Stage 1 (no DE needed).
# =============================================================================

rule preprocessing_report:
    input:
        reports = expand(str(QUANT / "{sample}/fastp.json"), sample=SAMPLES),
        ledgers = expand(str(QUANT / "{sample}/preprocessing.json"), sample=SAMPLES),
    output:
        html = OUT / "preprocessing_report/multiqc_report.html",
    threads: 1
    resources:
        mem_mb = 2000,
        runtime = 30,
    log:
        "logs/preprocessing_report.log",
    shell:
        "multiqc {input.reports:q} -o {OUT:q}/preprocessing_report -n multiqc_report -f > {log:q} 2>&1"

rule qc_report:
    input:
        quant = expand(str(QUANT / "{sample}" / QUANT_PRODUCT), sample=SAMPLES),
        star_reports = expand(str(QUANT / "{sample}/star/{report}"), sample=SAMPLES,
                              report=["Log.final.out", "featureCounts.tsv.summary"]) if STAR_ROUTE else [],
        sheet = config["samples"]["sheet"],
        # R17c: the report states reference construction from this record only.
        reference_lock = REF / ("star_idx/reference.lock.json" if STAR_ROUTE else "reference.lock.json"),
        report_script = REPO_DIR / "workflow/scripts/qc_report.py",
    output:
        done = OUT / "qc_report/qc_report.done",
        tsv  = OUT / "qc_report/alignment_summary.tsv",
        png  = OUT / "qc_report/comparative_charts.png",
        html = OUT / "qc_report/qc_charts.html",
        mqc  = OUT / "qc_report/multiqc/multiqc_report.html",
    params:
        runtime_identity = RUNTIME_ID,
        quantifier = "star" if STAR_ROUTE else "salmon",
        reports = [str(QUANT / sample / filename) for sample in SAMPLES
                   for filename in (["fastp.json", "star/Log.final.out", "star/featureCounts.tsv.summary"]
                                    if STAR_ROUTE else ["fastp.json", "aux_info/meta_info.json"])],
        quantdir = lambda wc: str(QUANT),
        outdir   = lambda wc: str(OUT / "qc_report"),
        # Report identity comes from this run's config, so the output describes
        # the run that produced it instead of carrying a baked-in project name.
        project  = lambda wc: config["project"].get("name") or "",
        species  = lambda wc: (config.get("organism", {}).get("scientific_name")
                               or config.get("organism", {}).get("common_name") or ""),

    threads: 2
    resources:
        mem_mb = 4000,
        runtime = 30,
    log:
        "logs/qc_report.log",
    shell:
        """
        set -euo pipefail
        mkdir -p {params.outdir:q}
        echo "[qc_report] MultiQC over {params.quantdir:q}" > {log:q}
        multiqc {params.reports:q} -o {params.outdir:q}/multiqc -n multiqc_report -f >> {log:q} 2>&1
        echo "[qc_report] comparative charts" >> {log:q}
        python {REPO_DIR:q}/workflow/scripts/qc_report.py \
            --quant {params.quantdir:q} \
            --quantifier {params.quantifier:q} \
            --out {params.outdir:q} \
            --samples {input.sheet:q} \
            --project={params.project:q} \
            --species={params.species:q} \
            --reference-lock {input.reference_lock:q} >> {log:q} 2>&1
        touch {output.done:q}
        """

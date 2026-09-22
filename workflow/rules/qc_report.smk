# =============================================================================
# QC report — auto-generated comparative charts + MultiQC.
# Runs once after all samples are quantified. Produces:
#   results/qc_report/comparative_charts.png  (before/after clean + alignment fate)
#   results/qc_report/qc_charts.html          (self-contained report)
#   results/qc_report/alignment_summary.tsv   (per-sample table)
#   results/qc_report/multiqc/multiqc_report.html
# Depends only on per-sample outputs, so it runs right after Stage 1 (no DE needed).
# =============================================================================

rule qc_report:
    input:
        quant = expand(str(QUANT / "{sample}/quant.sf"), sample=SAMPLES),
        sheet = config["samples"]["sheet"],
        # R17c: the report states reference construction from this record only.
        reference_lock = REF / "reference.lock.json",
        report_script = REPO_DIR / "workflow/scripts/qc_report.py",
    output:
        done = OUT / "qc_report/qc_report.done",
        tsv  = OUT / "qc_report/alignment_summary.tsv",
        png  = OUT / "qc_report/comparative_charts.png",
        html = OUT / "qc_report/qc_charts.html",
        mqc  = OUT / "qc_report/multiqc/multiqc_report.html",
    params:
        runtime_identity = RUNTIME_ID,
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
        multiqc {params.quantdir:q} -o {params.outdir:q}/multiqc -n multiqc_report -f >> {log:q} 2>&1
        echo "[qc_report] comparative charts" >> {log:q}
        python {REPO_DIR:q}/workflow/scripts/qc_report.py \
            --quant {params.quantdir:q} \
            --out {params.outdir:q} \
            --samples {input.sheet:q} \
            --project={params.project:q} \
            --species={params.species:q} \
            --reference-lock {input.reference_lock:q} >> {log:q} 2>&1
        touch {output.done:q}
        """

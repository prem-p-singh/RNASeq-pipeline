# =============================================================================
# Stage 1 — Per-sample QC + trim + quant.
# One SLURM task per sample. fastp clean -> Salmon (decoy-aware OK).
# Supports single-end (tagseq, rnaseq_single) and paired-end (rnaseq_paired).
# =============================================================================

rule qc_quant_sample:
    """Locate reads -> fastp -> salmon -> validate -> emit metrics.json."""
    input:
        idx = REF / "salmon_idx",
    output:
        quant = QUANT / "{sample}/quant.sf",
        metrics = QUANT / "{sample}/metrics.json",
    params:
        fastq_url = get_fastq_url,
        fastq_url_r2 = get_fastq_url_r2,
        seq_type = config["samples"]["seq_type"],
        threads = 4,
        sample_dir = lambda wc: str(QUANT / wc.sample),
        min_map_rate = lambda wc: config["thresholds"]["sample_qc"]["mapping_rate_min"],
        expected_libtype = lambda wc: config["samples"].get("expected_libtype", "") or "",
        # Owned intermediates (downloaded reads, trimmed reads) are removed once
        # quant.sf is verified. User-supplied FASTQs read in place are untouched.
        delete_intermediates = lambda wc: str(
            config["hpc"].get("delete_fastq_after_quant", True)
        ).lower(),
    threads: 4
    resources:
        mem_mb = 16000,    # decoy-aware salmon + --gcBias needs headroom; 6 GB stalls
        runtime = 90,      # minutes
        tmpdir = "tmp",
    log:
        "logs/qc_quant/{sample}.log",
    shell:
        """
        set -euo pipefail
        bash workflow/scripts/01_qc_quant.sh \
            --sample {wildcards.sample} \
            --url "{params.fastq_url}" \
            --url2 "{params.fastq_url_r2}" \
            --seq-type {params.seq_type} \
            --index {input.idx} \
            --outdir {params.sample_dir} \
            --threads {params.threads} \
            --min-map-rate {params.min_map_rate} \
            --expected-libtype "{params.expected_libtype}" \
            --delete-intermediates {params.delete_intermediates} \
            > {log} 2>&1
        """

# =============================================================================
# Stage 1 — Per-sample QC + trim + quant.
# One SLURM task per sample. fastp clean -> Salmon (decoy-aware OK).
# Supports single-end (tagseq, rnaseq_single) and paired-end (rnaseq_paired).
# =============================================================================

rule qc_quant_sample:
    """Locate reads -> fastp -> salmon -> validate -> emit metrics.json."""
    input:
        idx = REF / "salmon_idx",
        reads = local_read_inputs,
        stage_script = REPO_DIR / "workflow/scripts/01_qc_quant.sh",
        tools_helper = REPO_DIR / "workflow/scripts/_tools.sh",
        read_helper = REPO_DIR / "scripts/prepare_reads.py",
    output:
        quant = QUANT / "{sample}/quant.sf",
        metrics = QUANT / "{sample}/metrics.json",
        provenance = [str(QUANT / "{sample}/read_preparation.json")] if CANONICAL_INPUTS is not None else [],
    params:
        runtime_identity = RUNTIME_ID,
        fastq_url = lambda wc: shlex.quote(str(get_fastq_url(wc))),
        reads_json = lambda wc: shlex.quote(json.dumps(CANONICAL_INPUTS[wc.sample]) if CANONICAL_INPUTS is not None else ""),
        fastq_url_r2 = lambda wc: shlex.quote(str(get_fastq_url_r2(wc))),
        seq_type = config["samples"]["seq_type"],
        sample_dir = lambda wc: str(QUANT / wc.sample),
        min_map_rate = lambda wc: config["thresholds"]["sample_qc"]["mapping_rate_min"],
        expected_libtype = lambda wc: shlex.quote(CANONICAL_INPUTS[wc.sample]["expected_libtype"]
            if CANONICAL_INPUTS is not None else config["samples"].get("expected_libtype", "") or ""),
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
        bash {REPO_DIR:q}/workflow/scripts/01_qc_quant.sh \
            --sample {wildcards.sample:q} \
            --url {params.fastq_url} \
            --url2 {params.fastq_url_r2} \
            --reads-json {params.reads_json} \
            --seq-type {params.seq_type} \
            --index {input.idx:q} \
            --outdir {params.sample_dir:q} \
            --threads {threads} \
            --min-map-rate {params.min_map_rate} \
            --expected-libtype {params.expected_libtype} \
            --delete-intermediates {params.delete_intermediates} \
            > {log:q} 2>&1
        """

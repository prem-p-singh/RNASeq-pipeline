# =============================================================================
# Stage 1 — Per-sample QC + trim + quant.
# One SLURM task per sample. fastp clean -> Salmon (decoy-aware OK).
# Supports single-end (tagseq, rnaseq_single) and paired-end (rnaseq_paired).
# =============================================================================

# Reports are durable; only trimmed FASTQs may be removed after all consumers.
_trim_patterns = ([str(QUANT / "{sample}/{sample}_R1.trim.fastq.gz"),
                   str(QUANT / "{sample}/{sample}_R2.trim.fastq.gz")]
                  if config["samples"]["seq_type"] == "rnaseq_paired"
                  else [str(QUANT / "{sample}/{sample}.trim.fastq.gz")])

rule preprocess_sample:
    input:
        reads = local_read_inputs,
        stage_script = REPO_DIR / "workflow/scripts/01_qc_quant.sh",
        processing_script = REPO_DIR / "scripts/preprocessing.py",
        read_helper = REPO_DIR / "scripts/prepare_reads.py",
        metadata_helper = REPO_DIR / "scripts/metadata.py",
        tools_helper = REPO_DIR / "workflow/scripts/_tools.sh",
    output:
        trims = temp(_trim_patterns) if config["hpc"].get("delete_fastq_after_quant", True) else _trim_patterns,
        json = QUANT / "{sample}/fastp.json",
        html = QUANT / "{sample}/fastp.html",
        ledger = QUANT / "{sample}/preprocessing.json",
        provenance = [str(QUANT / "{sample}/read_preparation.json")] if CANONICAL_INPUTS is not None else [],
    params:
        runtime_identity = RUNTIME_ID,
        url = lambda wc: shlex.quote(str(get_fastq_url(wc))),
        url2 = lambda wc: shlex.quote(str(get_fastq_url_r2(wc))),
        reads_json = lambda wc: shlex.quote(json.dumps(CANONICAL_INPUTS[wc.sample]) if CANONICAL_INPUTS is not None else ""),
        delete_intermediates = str(config["hpc"].get("delete_fastq_after_quant", True)).lower(),
        policy = lambda wc: shlex.quote(json.dumps(config["preprocessing"], sort_keys=True)),
        seq_type = config["samples"]["seq_type"],
        outdir = lambda wc: str(QUANT / wc.sample),
    threads: 4
    resources:
        mem_mb = 4000,
        runtime = 90,
    log:
        "logs/preprocess/{sample}.log",
    shell:
        """
        bash {REPO_DIR:q}/workflow/scripts/01_qc_quant.sh --stage preprocess \
            --sample {wildcards.sample:q} --url {params.url} --url2 {params.url2} \
            --reads-json {params.reads_json} --processing-json {params.policy} \
            --delete-intermediates {params.delete_intermediates} \
            --seq-type {params.seq_type} --outdir {params.outdir:q} \
            --threads {threads} > {log:q} 2>&1
        """

rule qc_quant_sample:
    input:
        idx = REF / "salmon_idx",
        trims = lambda wc: [p.format(sample=wc.sample) for p in _trim_patterns],
        preprocessing = QUANT / "{sample}/preprocessing.json",
        fastp_json = QUANT / "{sample}/fastp.json",
        fastp_html = QUANT / "{sample}/fastp.html",
        stage_script = REPO_DIR / "workflow/scripts/01_qc_quant.sh",
        tools_helper = REPO_DIR / "workflow/scripts/_tools.sh",
    output:
        quant = QUANT / "{sample}/quant.sf",
        metrics = QUANT / "{sample}/metrics.json",
    params:
        runtime_identity = RUNTIME_ID,
        seq_type = config["samples"]["seq_type"],
        sample_dir = lambda wc: str(QUANT / wc.sample),
        min_map_rate = config["thresholds"]["sample_qc"]["mapping_rate_min"],
        expected_libtype = lambda wc: shlex.quote(CANONICAL_INPUTS[wc.sample]["expected_libtype"]
            if CANONICAL_INPUTS is not None else config["samples"].get("expected_libtype", "") or ""),
    threads: 4
    resources:
        mem_mb = 16000,
        runtime = 90,
        tmpdir = "tmp",
    log:
        "logs/qc_quant/{sample}.log",
    shell:
        """
        bash {REPO_DIR:q}/workflow/scripts/01_qc_quant.sh --stage quant \
            --sample {wildcards.sample:q} --seq-type {params.seq_type} \
            --index {input.idx:q} --outdir {params.sample_dir:q} \
            --threads {threads} --min-map-rate {params.min_map_rate} \
            --expected-libtype {params.expected_libtype} \
            --delete-intermediates false > {log:q} 2>&1
        """

rule fetch_screen_fasta:
    output: fa = "screening_panel/{name}/source.fa"
    params:
        url = lambda wc: SCREEN_BUILT[wc.name] if "://" in SCREEN_BUILT[wc.name]
              else Path(SCREEN_BUILT[wc.name]).expanduser().resolve().as_uri(),
    shell:
        "bash {REPO_DIR:q}/workflow/scripts/fetch_reference_file.sh --url {params.url:q} --out {output.fa:q} --kind fasta"

rule build_screen_index:
    input:
        fa = lambda wc: {"transcriptome": REF / "transcriptome.fa", "genome": REF / "genome.fa"}.get(
            SCREEN_BUILT[wc.name], f"screening_panel/{wc.name}/source.fa"),
    output: done = "screening_panel/{name}/built.txt"
    params: base = "screening_panel/{name}/{name}"
    threads: 4
    resources:
        mem_mb = 16000,
        runtime = 240,
    log: "logs/screening/index_{name}.log"
    shell:
        """
        set -euo pipefail
        bowtie2-build --threads {threads} {input.fa:q} {params.base:q} > {log:q} 2>&1
        {{ bowtie2-build --version | head -1; sha256sum {input.fa:q}; }} > {output.done:q}
        """

# A second consumer keeps temporary trimmed reads alive until screening finishes.
# Disabled screening still publishes an explicit status; it requires no tools.
rule screen_sample:
    input:
        trims = lambda wc: [p.format(sample=wc.sample) for p in _trim_patterns] if SCREEN_POLICY["enabled"] else [],
        preprocessing = QUANT / "{sample}/preprocessing.json",
        indexes = SCREEN_INDEXES,
        script = REPO_DIR / "scripts/screen_reads.py",
        reader = REPO_DIR / "scripts/prepare_reads.py",
    output:
        json = QUANT / "{sample}/screening/screening.json",
        tsv = QUANT / "{sample}/screening/screening.tsv",
        html = QUANT / "{sample}/screening/screening.html",
    params:
        runtime_identity = RUNTIME_ID,
        policy = json.dumps(SCREEN_RUN_POLICY, sort_keys=True),
        # Disabled mode does not read these paths (trimmed reads may be removed).
        reads = lambda wc: [p.format(sample=wc.sample) for p in _trim_patterns],
        outdir = lambda wc: str(QUANT / wc.sample / "screening"),
    threads: 2
    resources:
        mem_mb = 8000,
        runtime = 90,
    log: "logs/screening/{sample}.log"
    shell:
        "python {input.script:q} --reads {params.reads:q} --sample {wildcards.sample:q} "
        "--outdir {params.outdir:q} --policy-json {params.policy:q} --threads {threads} > {log:q} 2>&1"

# =============================================================================
# Stage 0 — Retrieve reference genome, annotation, and build Salmon index.
# Runs once per project. Outputs persist in reference/ (~1 GB total).
# =============================================================================

rule fetch_transcriptome:
    """Download and decompress the transcriptome FASTA."""
    output:
        fa = REF / "transcriptome.fa",
    params:
        url = config["reference"]["transcriptome_fasta_url"],
    shell:
        """
        set -euo pipefail
        mkdir -p {REF}
        curl -sSL "{params.url}" | gunzip -c > {output.fa}
        """

rule fetch_gtf:
    output:
        gtf = REF / "annotation.gtf",
    params:
        url = config["reference"]["gtf_url"],
    shell:
        """
        set -euo pipefail
        curl -sSL "{params.url}" | gunzip -c > {output.gtf}
        """

rule salmon_index:
    """Build a Salmon index.

    This builds a TRANSCRIPTOME-ONLY index: no `-d` decoy file is supplied.
    Nothing downstream may describe it as decoy-aware (RF02, D11). The
    reference_lock rule below records that fact, and the QC report prints the
    construction from that record rather than from a fixed string.
    """
    input:
        fa = REF / "transcriptome.fa",
    output:
        idx = directory(REF / "salmon_idx"),
        sentinel = REF / "salmon_idx/info.json",
    threads: 4
    resources:
        mem_mb = 8000,
        runtime = 60,
    shell:
        """
        set -euo pipefail
        source {REPO_DIR}/workflow/scripts/_tools.sh
        ensure_tools salmon
        salmon --version
        salmon index -t {input.fa} -i {output.idx} -k 31 --threads {threads}
        """


rule reference_lock:
    """Record the reference bundle's identity (RF01, master plan 12.1).

    Hashes the bundle files and captures the index parameters actually used,
    so every downstream claim about the reference is traceable to one record
    instead of being asserted in prose.
    """
    input:
        fa = REF / "transcriptome.fa",
        gtf = REF / "annotation.gtf",
        sentinel = REF / "salmon_idx/info.json",
    output:
        lock = REF / "reference.lock.json",
    params:
        idx = lambda wc: str(REF / "salmon_idx"),
        organism = lambda wc: (config.get("organism", {}) or {}).get("scientific_name") or "",
        tax_id = lambda wc: (config.get("organism", {}) or {}).get("tax_id") or 0,
        accession = lambda wc: (config.get("reference", {}) or {}).get("accession") or "",
        assembly = lambda wc: (config.get("reference", {}) or {}).get("assembly_name") or "",
        tx_url = lambda wc: (config.get("reference", {}) or {}).get("transcriptome_fasta_url") or "",
        gtf_url = lambda wc: (config.get("reference", {}) or {}).get("gtf_url") or "",
        genome_url = lambda wc: (config.get("reference", {}) or {}).get("genome_fasta_url") or "",
    resources:
        mem_mb = 2000,
        runtime = 30,
    shell:
        """
        set -euo pipefail
        source {REPO_DIR}/workflow/scripts/_tools.sh
        ensure_tools salmon
        python3 {REPO_DIR}/workflow/scripts/reference_lock.py \
            --out {output.lock} \
            --transcriptome {input.fa} \
            --gtf {input.gtf} \
            --index {params.idx} \
            --organism "{params.organism}" \
            --tax-id {params.tax_id} \
            --accession "{params.accession}" \
            --assembly-name "{params.assembly}" \
            --transcriptome-url "{params.tx_url}" \
            --gtf-url "{params.gtf_url}" \
            --genome-url "{params.genome_url}" \
            --kmer 31
        """

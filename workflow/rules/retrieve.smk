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
    """Build a Salmon index. Uses --gencode off by default — set in config if needed."""
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

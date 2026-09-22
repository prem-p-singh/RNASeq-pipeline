import json

# =============================================================================
# Stage 0 — Retrieve reference genome, annotation, and build Salmon index.
# Runs once per project. Outputs persist in reference/ (~1 GB total).
# =============================================================================

# R03 requires atomic transfer and format checks. `curl | gunzip > out`
# streams straight onto the final path, so an interrupted transfer leaves a
# truncated file there, and nothing verifies the payload is the format it
# claims. fetch_reference_file downloads beside the target, checks the archive,
# then renames: a reader sees the final name only once the content is whole
# (master plan 5.3, "temporary build directory, validation, atomic publication").
rule fetch_transcriptome:
    """Download, verify and decompress the transcriptome FASTA."""
    output:
        fa = REF / "transcriptome.fa",
    params:
        url = config["reference"]["transcriptome_fasta_url"],
        kind = "fasta",
    shell:
        """
        set -euo pipefail
        bash {REPO_DIR:q}/workflow/scripts/fetch_reference_file.sh \
            --url {params.url:q} --out {output.fa:q} --kind {params.kind}
        """

rule fetch_gtf:
    """Download, verify and decompress the annotation."""
    output:
        gtf = REF / "annotation.gtf",
    params:
        url = config["reference"]["gtf_url"],
        kind = "gtf",
    shell:
        """
        set -euo pipefail
        bash {REPO_DIR:q}/workflow/scripts/fetch_reference_file.sh \
            --url {params.url:q} --out {output.gtf:q} --kind {params.kind}
        """

rule salmon_index:
    """Build a Salmon index, under a lock, published atomically.

    This builds a TRANSCRIPTOME-ONLY index unless reference.decoys is set: no
    `-d` file is supplied by default. Nothing downstream may describe it as
    decoy-aware (RF02, D11).

    The index is built into a staging directory and renamed into place, so a
    concurrent reader never sees a half-written index, and a build lock stops
    two runs indexing the same cache key at once (master plan 5.3, R17d).
    reference.lock.json is written here rather than afterwards, because this is
    the only point where the Salmon version that actually did the build is
    known (WORKING_PLAN 3.1).
    """
    input:
        fa = REF / "transcriptome.fa",
        gtf = REF / "annotation.gtf",
    output:
        idx = directory(REF / "salmon_idx"),
        sentinel = REF / "salmon_idx/info.json",
        lock = REF / "reference.lock.json",
    params:
        staging = lambda wc: str(REF / ".salmon_idx.building"),
        lockdir = lambda wc: str(REF / "salmon_idx"),
        kmer = lambda wc: SALMON_KMER,
        decoys = lambda wc: SALMON_DECOYS or "",
        key_inputs = lambda wc: json.dumps(
            _refcache.key_inputs(config.get("reference") or {},
                                 SALMON_KMER, SALMON_DECOYS)),
        organism = lambda wc: (config.get("organism", {}) or {}).get("scientific_name") or "",
        tax_id = lambda wc: (config.get("organism", {}) or {}).get("tax_id") or 0,
        accession = lambda wc: (config.get("reference", {}) or {}).get("accession") or "",
        assembly = lambda wc: (config.get("reference", {}) or {}).get("assembly_name") or "",
        tx_url = lambda wc: (config.get("reference", {}) or {}).get("transcriptome_fasta_url") or "",
        gtf_url = lambda wc: (config.get("reference", {}) or {}).get("gtf_url") or "",
        genome_url = lambda wc: (config.get("reference", {}) or {}).get("genome_fasta_url") or "",
    threads: 4
    resources:
        mem_mb = 8000,
        runtime = 60,
    shell:
        """
        set -euo pipefail
        source {REPO_DIR:q}/workflow/scripts/_tools.sh
        ensure_tools salmon

        python3 {REPO_DIR:q}/workflow/scripts/build_salmon_index.py \
            --transcriptome {input.fa:q} \
            --gtf {input.gtf:q} \
            --final-index {output.idx:q} \
            --staging {params.staging:q} \
            --lock-out {output.lock:q} \
            --kmer {params.kmer} \
            --decoys={params.decoys:q} \
            --threads {threads} \
            --key-inputs {params.key_inputs:q} \
            --organism={params.organism:q} \
            --tax-id {params.tax_id} \
            --accession={params.accession:q} \
            --assembly-name={params.assembly:q} \
            --transcriptome-url={params.tx_url:q} \
            --gtf-url={params.gtf_url:q} \
            --genome-url={params.genome_url:q} \
            --helper {REPO_DIR:q}/workflow/scripts/reference_lock.py
        """

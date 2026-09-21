# =============================================================================
# build_orgdb.smk
# -----------------------------------------------------------------------------
# Builds an OrgDb package on demand using the tiered fallback in
# scripts/build_orgdb.R. The OrgDb is needed for Stage 4 (GO enrichment).
#
# Sentinel output lives at:
#   <cache_dir>/<taxid>_<assembly>/.orgdb_built
# so the same organism is built once and reused across projects.
# =============================================================================

from pathlib import Path


# --- helper: where to put the sentinel file for THIS project's organism ----
def orgdb_cache_dir():
    """Return the cache dir path from config, with ~ expanded."""
    return Path(config["orgdb"]["cache_dir"]).expanduser()


def orgdb_sentinel():
    """Return the per-organism sentinel path. The build rule writes this on success."""
    tax_id = config["organism"]["tax_id"]
    accession = config["reference"].get("accession", "unknown")
    return str(orgdb_cache_dir() / f"{tax_id}_{accession}" / ".orgdb_built")


def orgdb_required() -> bool:
    """Is an OrgDb actually needed for this run?

    Only GO enrichment consumes it. Asking for the build when GO is switched
    off costs a download (or a multi-hour eggNOG run) that nothing reads.
    """
    if config.get("orgdb", {}).get("strategy", "auto") == "skip":
        return False
    return bool(config.get("downstream", {}).get("run_go", False))


# --- rule -------------------------------------------------------------------
rule build_orgdb:
    """Ensure the OrgDb is installed in R; build it if needed."""
    input:
        # If eggNOG path is needed (Tier 3b), we need the transcriptome
        # already fetched. retrieve.smk handles that.
        transcriptome = REF / "transcriptome.fa",
    output:
        sentinel = orgdb_sentinel(),
    params:
        tax_id        = config["organism"]["tax_id"],
        # Best-effort split of "Genus species" → genus + species
        genus         = lambda wc: config["organism"]["scientific_name"].split()[0],
        species       = lambda wc: (config["organism"]["scientific_name"].split() + [""])[1],
        assembly      = lambda wc: config["reference"].get("accession", "unknown"),
        orgdb_package = config["organism"].get("orgdb_package", ""),
        cache_dir     = str(orgdb_cache_dir()),
        strategy      = config["orgdb"]["strategy"],
        # Tier 3b inputs — currently unused (eggNOG not installed) but wired up
        proteome_fa   = lambda wc: str(REF / "proteome.fa"),
        eggnog_db     = lambda wc: str(config["orgdb"].get("eggnog_db") or ""),
    resources:
        mem_mb = 8000,
        runtime = 30,      # 30 min covers Tier 3a; bump to 240 if Tier 3b is in use
    log:
        "logs/build_orgdb.log",
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.sentinel})

        Rscript {REPO_DIR}/scripts/build_orgdb.R \
            --tax_id {params.tax_id} \
            --genus "{params.genus}" \
            --species "{params.species}" \
            --assembly "{params.assembly}" \
            --orgdb_package "{params.orgdb_package}" \
            --cache_dir "{params.cache_dir}" \
            --strategy "{params.strategy}" \
            --proteome_fa "{params.proteome_fa}" \
            --eggnog_db "{params.eggnog_db}" \
            --sentinel "{output.sentinel}" \
            > {log} 2>&1
        """

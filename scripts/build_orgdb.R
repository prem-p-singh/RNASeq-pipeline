#!/usr/bin/env Rscript
# ===========================================================================
# build_orgdb.R — tiered OrgDb builder
# ---------------------------------------------------------------------------
# This script makes sure an OrgDb package (org.Xxx.eg.db) is available for
# the GO enrichment stage. Annotation availability varies by organism:
#
#   Tier 1  Bioconductor package already installed in R          instant
#   Tier 2  We built it before and cached it                      instant
#   Tier 3a NCBI Gene records (makeOrgPackageFromNCBI)            ~1 min
#   Tier 3b eggNOG-mapper on the proteome (makeOrgPackage)        1-3 hours
#
# The script tries each tier in the order above and stops at the first one
# that works. Results from Tier 3a/3b get cached to disk so the next project
# on the same organism is instant (Tier 2 hit).
#
# Designed to be readable — open it, read top to bottom, tweak freely.
# ===========================================================================

suppressPackageStartupMessages({
  library(optparse)
})

# ---------------------------------------------------------------------------
# Parse command-line arguments
# ---------------------------------------------------------------------------
option_list <- list(
  make_option("--tax_id",         type = "integer",   help = "NCBI taxonomy ID"),
  make_option("--genus",          type = "character", help = "Genus (e.g. Homo)"),
  make_option("--species",        type = "character", help = "Species epithet (e.g. sapiens)"),
  make_option("--assembly",       type = "character", default = "unknown",
              help = "Assembly accession (e.g. GCF_004011695.2); used in cache key"),
  make_option("--orgdb_package",  type = "character", default = "",
              help = "Bioconductor OrgDb name to check first (e.g. org.Hs.eg.db)"),
  make_option("--cache_dir",      type = "character", default = "~/orgdb_cache",
              help = "Directory where built OrgDb packages live"),
  make_option("--strategy",       type = "character", default = "auto",
              help = "auto | force_build | use_existing | skip"),
  make_option("--proteome_fa",    type = "character", default = "",
              help = "Path to translated proteome FASTA (Tier 3b input)"),
  make_option("--eggnog_db",      type = "character", default = "",
              help = "Path to eggNOG-mapper database (Tier 3b)"),
  make_option("--sentinel",       type = "character", default = "",
              help = "Path to write a sentinel file on success")
)
opt <- parse_args(OptionParser(option_list = option_list))

# Make every path absolute before anything calls setwd(). The build tiers
# change directory into the cache, after which a relative --proteome_fa or
# --sentinel would resolve somewhere else entirely.
absolutise <- function(p) {
  if (is.null(p) || !nzchar(p)) return(p)
  p <- path.expand(p)
  if (substr(p, 1, 1) == "/") p else file.path(getwd(), p)
}
opt$cache_dir   <- absolutise(opt$cache_dir)
opt$proteome_fa <- absolutise(opt$proteome_fa)
opt$eggnog_db   <- absolutise(opt$eggnog_db)
opt$sentinel    <- absolutise(opt$sentinel)

# This script's own directory, so it can find its siblings (parse_eggnog.R)
# no matter what the working directory has been changed to.
SCRIPT_DIR <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f) >= 1) {
    dirname(normalizePath(sub("^--file=", "", f[1])))
  } else {
    getwd()
  }
})


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

log_msg <- function(stage, msg) {
  cat(sprintf("[%s] [orgdb/%s] %s\n",
              format(Sys.time(), "%FT%T"), stage, msg))
}

# Where this organism's cache lives (e.g. ~/orgdb_cache/3827_GCF_004011695.2/)
cache_path_for_organism <- function() {
  file.path(opt$cache_dir, paste0(opt$tax_id, "_", opt$assembly))
}

# Annotation packages are data products, installed outside the locked runtime.
orgdb_lib <- file.path(cache_path_for_organism(), "R-library")
dir.create(orgdb_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(orgdb_lib, .libPaths()))

# A successful install() call is not proof the package is usable. Check it
# loads before any sentinel claims the build succeeded, otherwise a cached
# marker will make every future run skip a package that was never there.
verify_installed <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    log_msg("verify", paste0("install reported success but '", pkg,
                             "' is not loadable; refusing to write a sentinel"))
    quit(save = "no", status = 1)
  }
  log_msg("verify", paste0(pkg, " loads"))
}

# Write a small sentinel file marking which tier we ended on
write_sentinel <- function(tier, source_note = "") {
  if (nchar(opt$sentinel) == 0) return(invisible(NULL))
  dir.create(dirname(opt$sentinel), recursive = TRUE, showWarnings = FALSE)
  msg <- sprintf("tier=%s\norgdb_package=%s\ntax_id=%d\nassembly=%s\nbuilt_on=%s\nsource=%s\n",
                 tier, opt$orgdb_package, opt$tax_id, opt$assembly,
                 format(Sys.time(), "%FT%T"), source_note)
  writeLines(msg, opt$sentinel)
  log_msg("done", paste0("wrote sentinel: ", opt$sentinel))
}


# ---------------------------------------------------------------------------
# Strategy: "skip"  →  do nothing, write a sentinel saying we skipped
# ---------------------------------------------------------------------------
if (opt$strategy == "skip") {
  log_msg("skip", "strategy=skip; not building any OrgDb (enrichment will be skipped)")
  write_sentinel("skip", "strategy=skip")
  quit(save = "no", status = 0)
}


# ---------------------------------------------------------------------------
# TIER 1 — Bioconductor package already installed
# ---------------------------------------------------------------------------
# If the user named a Bioconductor package and it's installed, just use it.
# Skipped when strategy = "force_build" (we want a fresh build from reference).
# ---------------------------------------------------------------------------
if (opt$strategy %in% c("auto", "use_existing") &&
    nchar(opt$orgdb_package) > 0) {

  log_msg("tier1", paste0("checking Bioconductor for: ", opt$orgdb_package))

  if (requireNamespace(opt$orgdb_package, quietly = TRUE)) {
    log_msg("tier1", paste0("found ", opt$orgdb_package, " — using it"))
    write_sentinel("1_bioconductor", opt$orgdb_package)
    quit(save = "no", status = 0)
  } else {
    log_msg("tier1", paste0(opt$orgdb_package, " not installed in R"))

    if (opt$strategy == "use_existing") {
      log_msg("tier1",
              paste0("strategy=use_existing but ", opt$orgdb_package,
                     " is not installed. Install it via BiocManager or change strategy."))
      quit(save = "no", status = 1)
    }
  }
}


# ---------------------------------------------------------------------------
# TIER 2 — cache hit
# ---------------------------------------------------------------------------
# Did we build this OrgDb before? Each successful build leaves the package
# source dir + a .built marker in <cache_dir>/<tax_id>_<assembly>/.
# ---------------------------------------------------------------------------
cache_dir_for_org <- cache_path_for_organism()
cache_marker <- file.path(cache_dir_for_org, ".built")

if (opt$strategy != "force_build" && file.exists(cache_marker)) {
  log_msg("tier2", paste0("cache hit at ", cache_dir_for_org))

  # Cache contains the OrgDb package as an installable source dir
  pkg_dirs <- list.dirs(cache_dir_for_org, recursive = FALSE)
  pkg_dirs <- pkg_dirs[grepl("^org\\.", basename(pkg_dirs))]

  if (length(pkg_dirs) >= 1) {
    pkg_dir <- pkg_dirs[1]
    pkg_name <- basename(pkg_dir)

    if (!requireNamespace(pkg_name, quietly = TRUE)) {
      log_msg("tier2", paste0("installing cached package: ", pkg_name))
      install.packages(pkg_dir, repos = NULL, type = "source", lib = orgdb_lib)
    }
    verify_installed(pkg_name)
    log_msg("tier2", paste0("using cached ", pkg_name))
    write_sentinel("2_cache", pkg_dir)
    quit(save = "no", status = 0)
  }
  # If cache exists but no usable package — fall through and rebuild
  log_msg("tier2", "cache marker exists but no installable package; rebuilding")
}


# ---------------------------------------------------------------------------
# TIER 3a — build from NCBI Gene records
# ---------------------------------------------------------------------------
# AnnotationForge::makeOrgPackageFromNCBI() pulls gene + GO + KEGG data from
# NCBI Gene for the given tax_id, builds an org.<Xx>.eg.db source tree, and
# we install it from source.
# ---------------------------------------------------------------------------
if (opt$strategy %in% c("auto", "force_build")) {

  log_msg("tier3a", paste0("attempting NCBI Gene build for tax_id=",
                            opt$tax_id, " (", opt$genus, " ", opt$species, ")"))

  if (!requireNamespace("AnnotationForge", quietly = TRUE)) {
    stop("AnnotationForge is missing; repair the release runtime with scripts/bootstrap.sh")
  }
  library(AnnotationForge)

  # Build in the cache directory so the source tree is preserved + reusable
  dir.create(cache_dir_for_org, recursive = TRUE, showWarnings = FALSE)
  old_wd <- setwd(cache_dir_for_org)
  on.exit(setwd(old_wd), add = TRUE)

  pkg_name <- tryCatch(
    makeOrgPackageFromNCBI(
      version    = "0.1",
      author     = "Auto-generated <noreply@example.org>",
      maintainer = "Auto-generated <noreply@example.org>",
      outputDir  = ".",
      tax_id     = as.character(opt$tax_id),
      genus      = opt$genus,
      species    = opt$species
    ),
    error = function(e) {
      log_msg("tier3a", paste0("NCBI build failed: ", e$message))
      return(NULL)
    }
  )

  if (!is.null(pkg_name)) {
    log_msg("tier3a", paste0("built ", pkg_name, " — installing from source"))
    install.packages(pkg_name, repos = NULL, type = "source", lib = orgdb_lib)
    verify_installed(pkg_name)

    # Drop a marker so Tier 2 finds it next time
    file.create(cache_marker)
    log_msg("tier3a", paste0("cached at ", cache_dir_for_org))
    write_sentinel("3a_ncbi", file.path(cache_dir_for_org, pkg_name))
    quit(save = "no", status = 0)
  }
  # Fall through to Tier 3b if NCBI failed
}


# ---------------------------------------------------------------------------
# TIER 3b — build from eggNOG-mapper output
# ---------------------------------------------------------------------------
# This is the universal fallback. It works for ANY genome, including custom
# assemblies and non-model organisms. Cost: ~1-3 hours of eggNOG-mapper +
# ~50 GB database (downloaded once per machine).
#
# This tier is NOT yet runnable because eggnog-mapper is not installed in
# the analysis Conda environment and the eggNOG database has not been downloaded. The
# code below is here so it's wired up and ready when those two prereqs
# arrive. Until then, it prints a helpful error and exits.
# ---------------------------------------------------------------------------
if (opt$strategy %in% c("auto", "force_build")) {

  log_msg("tier3b", "NCBI build failed; would fall back to eggNOG-mapper here")

  eggnog_runner <- Sys.which("emapper.py")

  if (nchar(eggnog_runner) == 0) {
    log_msg("tier3b",
            "eggnog-mapper not installed. To enable Tier 3b:")
    log_msg("tier3b",
            "  conda install -n rnaseq-pipeline -c bioconda eggnog-mapper")
    log_msg("tier3b",
            "  download_eggnog_data.py --data_dir /path/to/db    # ~50 GB")
    log_msg("tier3b",
            "  then set orgdb.eggnog_db in config.yaml")
    quit(save = "no", status = 1)
  }

  if (nchar(opt$eggnog_db) == 0) {
    log_msg("tier3b", "no --eggnog_db path provided; cannot run")
    quit(save = "no", status = 1)
  }

  if (nchar(opt$proteome_fa) == 0 || !file.exists(opt$proteome_fa)) {
    log_msg("tier3b", paste0("--proteome_fa not found: ", opt$proteome_fa))
    quit(save = "no", status = 1)
  }

  # === Tier 3b implementation (ready to use when prereqs land) ===
  # 1. Run eggnog-mapper on the proteome
  log_msg("tier3b", paste0("running emapper.py on ", opt$proteome_fa))
  emapper_out <- file.path(cache_dir_for_org, "eggnog")
  dir.create(emapper_out, recursive = TRUE, showWarnings = FALSE)
  emapper_cmd <- sprintf(
    "%s -i %s --output eggnog --output_dir %s --data_dir %s --cpu 4 -m diamond",
    eggnog_runner, opt$proteome_fa, emapper_out, opt$eggnog_db
  )
  rc <- system(emapper_cmd)
  if (rc != 0) {
    log_msg("tier3b", "emapper.py failed")
    quit(save = "no", status = 1)
  }

  # 2. Parse the .emapper.annotations TSV → tables for AnnotationForge
  annotations_file <- file.path(emapper_out, "eggnog.emapper.annotations")
  log_msg("tier3b", paste0("parsing ", annotations_file))
  source(file.path(SCRIPT_DIR, "parse_eggnog.R"))
  tables <- parse_eggnog_to_tables(annotations_file)

  # 3. Build OrgDb from tables
  log_msg("tier3b", "building OrgDb from eggNOG annotations")
  library(AnnotationForge)
  old_wd <- setwd(cache_dir_for_org)
  on.exit(setwd(old_wd), add = TRUE)

  pkg_name <- makeOrgPackage(
    gene_info = tables$gene_info,
    go        = tables$go,
    version   = "0.1",
    author    = "Auto-generated <noreply@example.org>",
    maintainer= "Auto-generated <noreply@example.org>",
    outputDir = ".",
    tax_id    = as.character(opt$tax_id),
    genus     = opt$genus,
    species   = opt$species,
    goTable   = "go"
  )

  install.packages(pkg_name, repos = NULL, type = "source", lib = orgdb_lib)
  verify_installed(pkg_name)
  file.create(cache_marker)
  log_msg("tier3b", paste0("cached at ", cache_dir_for_org))
  write_sentinel("3b_eggnog", file.path(cache_dir_for_org, pkg_name))
  quit(save = "no", status = 0)
}


# ---------------------------------------------------------------------------
# Fell through every tier — nothing worked
# ---------------------------------------------------------------------------
log_msg("fail", "no tier succeeded; OrgDb not available")
quit(save = "no", status = 1)

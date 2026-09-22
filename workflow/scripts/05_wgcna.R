#!/usr/bin/env Rscript
# 05_wgcna.R — Adaptive WGCNA co-expression network.
#
# Gates:
#   * Skipped if run_wgcna = false or n_samples < min_samples.
#   * Soft-threshold power: lowest with scale-free R² >= target_r2, capped at 30.
#   * If >max_modules detected, re-run with bumped mergeCutHeight.
#
# Inputs (named):  counts, de_flag
# Outputs (named): done, metrics
# Params:          run_wgcna, min_samples

suppressPackageStartupMessages({
  library(edgeR)
  library(dplyr)
  library(readr)
  library(jsonlite)
  library(tibble)
})

cfg         <- snakemake@config
thr         <- cfg$thresholds$wgcna
counts_path <- snakemake@input$counts
out_done    <- snakemake@output$done
out_metrics <- snakemake@output$metrics
run_wgcna   <- isTRUE(snakemake@params$run_wgcna)
min_samples <- snakemake@params$min_samples

out_dir     <- file.path(cfg$project$output_dir, "WGCNA")

write_skip <- function(reason) {
  dir.create(dirname(out_metrics), recursive = TRUE, showWarnings = FALSE)
  write_json(list(skipped = TRUE, reason = reason),
             out_metrics, pretty = TRUE, auto_unbox = TRUE)
  cat(sprintf("[%s] WGCNA_SKIP: %s\n",
              format(Sys.time(), "%FT%T"), reason),
      file = "gates/decisions.log", append = TRUE)
  file.create(out_done)
}

# --- Gate: user opt-out -----------------------------------------------
if (!run_wgcna) { write_skip("run_wgcna = false"); quit(save = "no") }

# --- Gate: sample count -----------------------------------------------
counts <- read_tsv(counts_path, show_col_types = FALSE) %>%
  as.data.frame() %>% column_to_rownames("gene_id") %>% as.matrix()
n <- ncol(counts)
if (n < min_samples) {
  write_skip(sprintf("n_samples=%d < min_samples=%d", n, min_samples))
  quit(save = "no")
}

suppressPackageStartupMessages({
  library(WGCNA)
})
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- Prepare expression matrix (log2 CPM, filtered) -------------------
d <- DGEList(counts) %>% calcNormFactors()
# Use the same filter the DE stage did, approximated here
keep <- filterByExpr(d, design = matrix(1, ncol(counts), 1))
d <- d[keep, ]
datExpr <- t(cpm(d, log = TRUE, prior.count = 0.5))
message("datExpr: ", nrow(datExpr), " samples x ", ncol(datExpr), " genes")

# --- Validate the expression matrix -----------------------------------
# Must run before pickSoftThreshold: bicor on a zero-variance gene returns NA
# and poisons the scale-free fit. goodSamplesGenes is WGCNA's own check for
# genes/samples with too many missing or zero-variance entries.
gsg <- goodSamplesGenes(datExpr, verbose = 0)
n_genes_dropped   <- sum(!gsg$goodGenes)
n_samples_dropped <- sum(!gsg$goodSamples)
if (!gsg$allOK) {
  message("Dropping ", n_genes_dropped, " genes and ", n_samples_dropped,
          " samples that failed goodSamplesGenes")
  cat(sprintf("[%s] WGCNA_FILTER: dropped %d genes, %d samples (goodSamplesGenes)\n",
              format(Sys.time(), "%FT%T"), n_genes_dropped, n_samples_dropped),
      file = "gates/decisions.log", append = TRUE)
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
}
# Re-gate: the filter can push the cohort back under the minimum.
n <- nrow(datExpr)
if (n < min_samples) {
  write_skip(sprintf("n_samples=%d < min_samples=%d after goodSamplesGenes filter",
                     n, min_samples))
  quit(save = "no")
}

# --- Size the computation to the allocation ----------------------------
# WGCNA holds an n_genes^2 matrix of doubles per block, several times over.
# blockSize() derives a safe per-block gene count from the memory this job
# actually requested; blockwiseModules then pre-clusters genes into blocks.
mem_mb    <- snakemake@resources[["mem_mb"]]
if (is.null(mem_mb)) mem_mb <- 8000   # standalone run outside SLURM
n_threads_requested <- max(1L, as.integer(snakemake@threads))
n_threads <- n_threads_requested
max_block <- WGCNA::blockSize(ncol(datExpr), rectangularBlocks = TRUE,
                              maxMemoryAllocation = mem_mb * 1024^2,
                              overheadFactor = 3)
message("maxBlockSize = ", max_block, " genes (from mem_mb=", mem_mb,
        "), nThreads = ", n_threads)

# --- Threading, guarded against an unknown core count -----------------
# WGCNA::enableWGCNAThreads() calls parallel::detectCores() unconditionally,
# before it looks at nThreads, and then evaluates `nThreads > nCores`. Where
# detectCores() returns NA (restricted shells, some containers and schedulers)
# that comparison is NA, and the call aborts with "missing value where
# TRUE/FALSE needed" even though nThreads was supplied explicitly.
#
# Verified against WGCNA 1.73 by forcing detectCores() to NA: only
# enableWGCNAThreads fails; disableWGCNAThreads() and
# blockwiseModules(nThreads = 1) are unaffected. So an unusable core count
# costs parallelism, not the analysis.
n_cores <- suppressWarnings(parallel::detectCores())
if (!is.na(n_cores) && n_cores >= 2L && n_threads >= 2L) {
  n_threads <- min(n_threads, n_cores)
  enableWGCNAThreads(nThreads = n_threads)
} else {
  disableWGCNAThreads()
  message("WGCNA threading disabled (detectCores=", n_cores,
          ", requested threads=", n_threads, "); running single-threaded")
  n_threads <- 1L
}

# --- Adaptive soft-threshold power ------------------------------------
pow_scan <- pickSoftThreshold(datExpr,
                              networkType = "signed",
                              powerVector = 1:60,
                              corFnc = "bicor",
                              RsquaredCut = thr$target_r2,
                              verbose = 0)
power_estimate <- pow_scan$powerEstimate
power_use <- if (is.na(power_estimate) || power_estimate > thr$power_cap)
               thr$power_cap else power_estimate
r2_at_used <- pow_scan$fitIndices[match(power_use, pow_scan$fitIndices$Power), "SFT.R.sq"]

if (is.na(power_estimate)) {
  cat(sprintf("[%s] WGCNA_POWER: no power hit R^2 >= %.2f; using cap=%d (R^2=%.2f)\n",
              format(Sys.time(), "%FT%T"), thr$target_r2, thr$power_cap, r2_at_used),
      file = "gates/decisions.log", append = TRUE)
}
message("Using soft-threshold power = ", power_use, " (R^2 = ",
        round(r2_at_used, 3), ")")

# --- Module detection (with auto-merge gate) --------------------------
detect_modules <- function(mergeCutHeight) {
  blockwiseModules(datExpr, power = power_use,
                   corType = "bicor", networkType = "signed",
                   mergeCutHeight = mergeCutHeight,
                   maxBlockSize = max_block,
                   nThreads = n_threads,
                   verbose = 0)
}
# AN06: "Preserve declared parameter policy; optional sensitivity analysis is a
# new recorded analysis, not hidden retuning."
#
# This previously re-ran with a bumped merge height whenever the module count
# exceeded a cap, and the second network replaced the first. The published
# result was therefore whichever height happened to produce a more attractive
# module count, which is tuning an analysis to its own output. The declared
# height is now the reported result, full stop.
merge_used <- thr$merge_cut_height
if (is.null(merge_used)) {
  stop("thresholds.wgcna.merge_cut_height is not set. It is the declared ",
       "parameter policy and has no automatic fallback.")
}
net <- detect_modules(merge_used)
n_modules <- length(unique(net$colors))

# A high module count is a diagnostic about the data, not a reason to change
# the analysis. Recorded and reported; nothing is re-run because of it.
modules_exceed_diagnostic <- !is.null(thr$max_modules_diagnostic) &&
  n_modules > thr$max_modules_diagnostic
if (modules_exceed_diagnostic) {
  message("Module count ", n_modules, " exceeds the diagnostic threshold ",
          thr$max_modules_diagnostic,
          ". Reported as-is; the merge height was NOT changed (AN06).")
  cat(sprintf(paste0("[%s] WGCNA_DIAGNOSTIC: %d modules at declared ",
                     "mergeCutHeight %.2f exceeds max_modules_diagnostic %d; ",
                     "no retuning performed\n"),
              format(Sys.time(), "%FT%T"), n_modules, merge_used,
              thr$max_modules_diagnostic),
      file = "gates/decisions.log", append = TRUE)
}

# Optional sensitivity analysis. Each requested height is an ADDITIONAL
# recorded result; none of them replaces the declared one.
sensitivity <- list()
heights <- thr$sensitivity_merge_heights
if (!is.null(heights) && length(heights) > 0) {
  for (h in heights) {
    alt <- detect_modules(h)
    k <- length(unique(alt$colors))
    sensitivity[[length(sensitivity) + 1L]] <- list(
      merge_cut_height = h, n_modules = k)
    message("Sensitivity: mergeCutHeight ", h, " gives ", k, " modules ",
            "(recorded alongside the declared result, not in place of it)")
    write_tsv(data.frame(gene_id = colnames(datExpr), module = alt$colors),
              file.path(out_dir, sprintf("module_assignments_h%s.tsv", h)))
  }
}

# --- Outputs -----------------------------------------------------------
MEs <- net$MEs
MEs$sample_id <- rownames(MEs)
write_tsv(MEs[, c("sample_id", setdiff(names(MEs), "sample_id"))],
          file.path(out_dir, "MEs.tsv"))

module_map <- data.frame(
  gene_id = colnames(datExpr),
  module  = net$colors
)
write_tsv(module_map, file.path(out_dir, "module_assignments.tsv"))

# Save the full network object in case user wants to explore further
saveRDS(net, file.path(out_dir, "network.rds"))

metrics <- list(
  skipped            = FALSE,
  n_samples          = n,
  n_genes            = ncol(datExpr),
  n_genes_dropped    = n_genes_dropped,
  n_samples_dropped  = n_samples_dropped,
  power_estimate     = power_estimate,
  power_used         = power_use,
  r2_at_power_used   = r2_at_used,
  n_modules          = n_modules,
  merge_cut_height          = merge_used,
  merge_cut_height_policy   = "declared; never retuned to reach a module count (AN06)",
  max_modules_diagnostic    = thr$max_modules_diagnostic,
  modules_exceed_diagnostic = modules_exceed_diagnostic,
  sensitivity_analyses      = sensitivity,
  max_block_size     = max_block,
  n_blocks           = length(unique(net$blocks)),
  mem_mb             = mem_mb,
  n_threads_requested = n_threads_requested,
  n_threads          = n_threads,
  detected_cores     = if (is.na(n_cores)) NA_integer_ else as.integer(n_cores)
)
write_json(metrics, out_metrics, pretty = TRUE, auto_unbox = TRUE)
file.create(out_done)
message("WGCNA done. Modules: ", n_modules)

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
n_threads <- max(1L, as.integer(snakemake@threads))
max_block <- WGCNA::blockSize(ncol(datExpr), rectangularBlocks = TRUE,
                              maxMemoryAllocation = mem_mb * 1024^2,
                              overheadFactor = 3)
message("maxBlockSize = ", max_block, " genes (from mem_mb=", mem_mb,
        "), nThreads = ", n_threads)

# --- Adaptive soft-threshold power ------------------------------------
enableWGCNAThreads(nThreads = n_threads)
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
# merge_used tracks the height actually passed to the run that produced `net`.
# Deriving it afterwards from n_modules was wrong: a successful bump lowers
# n_modules below max_modules, so the metric reported the default height.
merge_used <- thr$merge_cut_height_default
net <- detect_modules(merge_used)
n_modules <- length(unique(net$colors))
if (n_modules > thr$max_modules) {
  message("Too many modules (", n_modules, "); re-running with bumped mergeCutHeight")
  cat(sprintf("[%s] WGCNA_MERGE: bumped mergeCutHeight to %.2f (was %d modules)\n",
              format(Sys.time(), "%FT%T"),
              thr$merge_cut_height_bumped, n_modules),
      file = "gates/decisions.log", append = TRUE)
  merge_used <- thr$merge_cut_height_bumped
  net <- detect_modules(merge_used)
  n_modules <- length(unique(net$colors))
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
  merge_cut_height   = merge_used,
  max_block_size     = max_block,
  n_blocks           = length(unique(net$blocks)),
  mem_mb             = mem_mb,
  n_threads          = n_threads
)
write_json(metrics, out_metrics, pretty = TRUE, auto_unbox = TRUE)
file.create(out_done)
message("WGCNA done. Modules: ", n_modules)

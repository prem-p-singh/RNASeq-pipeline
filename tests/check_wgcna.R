#!/usr/bin/env Rscript
# Self-check for workflow/scripts/05_wgcna.R.
#
# Runs the real script against a mock `snakemake` object and synthetic counts,
# then asserts on the metrics it writes. What must not break:
#   - merge_cut_height records the height actually used, not one re-derived
#     from the post-bump module count
#   - maxBlockSize follows the memory allocation, so a small allocation splits
#     the matrix into several blocks instead of forcing one
#   - goodSamplesGenes runs before pickSoftThreshold and its drops are recorded
#   - an unusable core count costs parallelism, not the analysis: WGCNA's
#     enableWGCNAThreads calls detectCores() unconditionally and then evaluates
#     `nThreads > nCores`, which aborts on NA. Guarded in 05_wgcna.R; run this
#     file with parallel::detectCores() forced to NA to exercise that path.
#
# Run:  Rscript tests/check_wgcna.R
suppressPackageStartupMessages({
  library(jsonlite)
})

ROOT   <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))), ".."))
SCRIPT <- file.path(ROOT, "workflow", "scripts", "05_wgcna.R")
stopifnot(file.exists(SCRIPT))

setClass("MockSnakemake",
         representation(config = "list", input = "list", output = "list",
                        params = "list", resources = "list", threads = "numeric"))

# --- Synthetic counts: 3 planted co-expression groups ------------------
make_counts <- function(n_samples, n_genes, seed = 1, constant_genes = 0) {
  set.seed(seed)
  lat <- matrix(rnorm(n_samples * 3), n_samples, 3)
  mu  <- lat[, rep(1:3, length.out = n_genes)] * 0.8
  m   <- matrix(rnbinom(n_samples * n_genes, mu = 200 * exp(mu), size = 8),
                n_samples, n_genes)
  if (constant_genes > 0) {
    # Constant RAW counts. Note these are not constant after normalisation:
    # each sample has its own library size and norm factor, so log-CPM varies.
    # See case 4 on why no positive zero-variance case is asserted.
    m[, seq_len(constant_genes)] <- 500L
  }
  out <- as.data.frame(t(m))
  colnames(out) <- paste0("S", seq_len(n_samples))
  cbind(gene_id = paste0("g", seq_len(n_genes)), out)
}

# Runs the script in a throwaway directory and returns the parsed metrics.
run_wgcna <- function(counts_df, mem_mb, max_modules, threads = 2,
                      sensitivity = list()) {
  wd <- file.path(tempdir(), paste0("wgcna_", as.integer(runif(1, 1, 1e9))))
  dir.create(file.path(wd, "gates"), recursive = TRUE)
  old <- setwd(wd); on.exit(setwd(old), add = TRUE)

  write.table(counts_df, "counts.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

  snakemake <- new("MockSnakemake",
    config = list(
      project = list(output_dir = "results"),
      thresholds = list(wgcna = list(
        min_samples = 15, target_r2 = 0.85, power_cap = 30,
        max_modules_diagnostic = max_modules,
        merge_cut_height = 0.25,
        sensitivity_merge_heights = sensitivity))),
    input     = list(counts = "counts.tsv"),
    output    = list(done = "wgcna.done", metrics = "wgcna.json"),
    params    = list(run_wgcna = TRUE, min_samples = 15),
    resources = list(mem_mb = mem_mb),
    threads   = threads)

  env <- new.env(parent = globalenv())
  assign("snakemake", snakemake, envir = env)
  sys.source(SCRIPT, envir = env)
  fromJSON("wgcna.json")
}

cts <- make_counts(20, 900)

# --- 1. the declared height is the result, whatever the module count ---
# AN06 forbids retuning to reach an attractive module count. Previously a
# module count above the cap triggered a second run at a bumped height whose
# network replaced the first, so the published result depended on its own
# output. max_modules_diagnostic = 0 forces the diagnostic to trip.
m <- run_wgcna(cts, mem_mb = 8000, max_modules = 0)
stopifnot(identical(m$skipped, FALSE))
if (!isTRUE(all.equal(m$merge_cut_height, 0.25)))
  stop("declared merge_cut_height was not used: ", m$merge_cut_height)
stopifnot(isTRUE(m$modules_exceed_diagnostic))
stopifnot(length(m$sensitivity_analyses) == 0)

# --- 2. a module count within the threshold behaves identically --------
m2 <- run_wgcna(cts, mem_mb = 8000, max_modules = 999)
if (!isTRUE(all.equal(m2$merge_cut_height, 0.25)))
  stop("declared merge_cut_height was not used: ", m2$merge_cut_height)
stopifnot(!isTRUE(m2$modules_exceed_diagnostic))
# the network itself must not depend on whether the diagnostic tripped
stopifnot(identical(m$n_modules, m2$n_modules))

m <- m2
if (m$n_blocks != 1)
  stop("900 genes at 8 GB should fit one block, got ", m$n_blocks)
if (m$max_block_size < m$n_genes)
  stop("maxBlockSize ", m$max_block_size, " should cover all ", m$n_genes, " genes at 8 GB")
# Threads: the request comes from snakemake@threads; what was USED may be
# lower, because enableWGCNAThreads is skipped when parallel::detectCores()
# cannot report a usable core count. Both are recorded so the metric describes
# the run rather than the intent (R16).
if (m$n_threads_requested != 2)
  stop("n_threads_requested not taken from snakemake@threads: ", m$n_threads_requested)
if (!(m$n_threads >= 1 && m$n_threads <= m$n_threads_requested))
  stop("n_threads used (", m$n_threads, ") outside 1..", m$n_threads_requested)

# --- 2b. sensitivity is additional, never a replacement ----------------
ms <- run_wgcna(cts, mem_mb = 8000, max_modules = 999,
                sensitivity = list(0.4, 0.9))
stopifnot(isTRUE(all.equal(ms$merge_cut_height, 0.25)))
stopifnot(nrow(ms$sensitivity_analyses) == 2 ||
          length(ms$sensitivity_analyses) == 2)
# the declared result is unchanged by running a sensitivity analysis
stopifnot(identical(ms$n_modules, m2$n_modules))

# --- 3. a small allocation must split into blocks, not force one ------
# 5 MB admits only a few hundred genes, so the same matrix has to be blocked.
m <- run_wgcna(cts, mem_mb = 5, max_modules = 999)
if (m$max_block_size >= m$n_genes)
  stop("maxBlockSize ", m$max_block_size, " ignored the 5 MB allocation")
if (m$n_blocks < 2)
  stop("expected >1 block at 5 MB, got ", m$n_blocks)

# --- 4. the goodSamplesGenes guard ran and reported ------------------
# No positive-drop case is asserted because none is reachable through the
# upstream filters: calcNormFactors errors on NA counts, filterByExpr removes
# all-zero genes, and a constant raw count is NOT constant in log-CPM (each
# sample gets its own normalisation factor). The guard is kept because it is
# WGCNA's documented precondition for network construction, and asserting the
# counts are reported means deleting the block fails this test.
m <- run_wgcna(make_counts(20, 900, constant_genes = 5), mem_mb = 8000, max_modules = 999)
stopifnot(is.numeric(m$n_genes_dropped), is.numeric(m$n_samples_dropped))
if (m$n_samples_dropped != 0)
  stop("clean input dropped ", m$n_samples_dropped, " samples; filter is too aggressive")

cat("check_wgcna.R: all assertions passed\n")

#!/usr/bin/env Rscript
# =============================================================================
# _design.R — study-design checks, shared by preflight and the DE stage.
#
# These two functions used to live in 02_aggregate.R (biological_replicates) and
# 03_de.R (validate_design). They are here because three callers need them and
# two copies would drift:
#
#   scripts/preflight_design.R   before anything expensive runs
#   02_aggregate.R               records the replicate count in its metrics
#   03_de.R                      refuses to fit an unsupported design
#
# Preflight and Stage 3 reaching different verdicts on the same sheet would be
# worse than having no preflight check at all, so they call the same code.
#
# Sourced, not run. No snakemake@ references, no library() calls: base R only,
# so preflight can use it without the Bioconductor stack installed.
# =============================================================================


# Keep identifier spelling while retaining missing-value and numeric semantics
# for model covariates. Shared by preflight, aggregation and differential testing.
read_sample_sheet <- function(path) {
  sheet <- read.delim(path, comment.char = "#", colClasses = "character",
                     na.strings = NULL, check.names = FALSE)
  for (name in setdiff(names(sheet), "sample_id")) {
    sheet[[name]] <- type.convert(sheet[[name]], as.is = TRUE, na.strings = c("", "NA"))
  }
  sheet
}

#' Independent biological replicates in the smallest level of `primary`.
#'
#' Rows are not replicates. When a random-effects term names a grouping unit
#' (e.g. "(1|vine)"), repeated measurements of one unit count once, because
#' technical or longitudinal repeats do not add biological replication.
#'
#' @param sheet data.frame of the sample sheet.
#' @param primary name of the primary factor column.
#' @param random_effects random-effects term as a string, or NULL.
#' @return list(n, unit, rows_per_group). `n` is NA when `primary` is absent.
biological_replicates <- function(sheet, primary, random_effects) {
  unit <- NA_character_
  if (!is.null(random_effects) && nzchar(random_effects)) {
    m <- regmatches(random_effects,
                    regexpr("\\|[[:space:]]*[A-Za-z._][A-Za-z0-9._]*",
                            random_effects))
    if (length(m) == 1) unit <- trimws(sub("\\|", "", m))
  }
  rows_per_group <- if (primary %in% names(sheet)) {
    min(table(sheet[[primary]]))
  } else NA_integer_
  n <- if (!is.na(unit) && unit %in% names(sheet) && primary %in% names(sheet)) {
    min(tapply(sheet[[unit]], sheet[[primary]],
               function(x) length(unique(x))))
  } else {
    rows_per_group
  }
  list(n = n,
       unit = if (is.na(unit)) "sample" else unit,
       rows_per_group = rows_per_group)
}


#' Is the requested model estimable on this design matrix?
#'
#' Stops rather than substituting. A design that cannot support the requested
#' model is a study-design problem; quietly falling back to a simpler test would
#' answer a different question and drop the interactions, blocking and random
#' effects the plan asked for.
#'
#' @return list(rank, residual_df) on success; stops with a specific reason
#'   otherwise.
validate_design <- function(mm, n_bio_rep, primary, unit = "sample") {
  mm_qr <- qr(mm)
  if (mm_qr$rank < ncol(mm)) {
    aliased <- colnames(mm)[mm_qr$pivot[(mm_qr$rank + 1L):ncol(mm)]]
    stop("Design is rank-deficient: ", ncol(mm), " coefficients but rank ",
         mm_qr$rank, ". Confounded/aliased term(s): ",
         paste(aliased, collapse = ", "),
         ". Drop the confounded term or supply a design that can estimate it.")
  }
  residual_df <- nrow(mm) - mm_qr$rank
  if (residual_df < 1L) {
    stop("No residual degrees of freedom (", nrow(mm), " samples, ",
         mm_qr$rank, " coefficients): variance cannot be estimated. ",
         "Reduce model terms or add replicates.")
  }
  if (!is.null(n_bio_rep) && !is.na(n_bio_rep) && n_bio_rep < 2L) {
    stop("Fewer than 2 independent biological replicates in the smallest '",
         primary, "' group (n = ", n_bio_rep, ", unit = ", unit,
         "). Differential expression is not supported for this design.")
  }
  list(rank = mm_qr$rank, residual_df = residual_df)
}

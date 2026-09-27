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


build_contrast_matrix <- function(dummy_fit, specs, sheet, model_vars) {
  if (is.null(specs) || length(specs) == 0) {
    stop("config$contrasts is empty - there is nothing to test.")
  }
  # Reject the removed R-expression format loudly rather than guessing.
  looks_legacy <- vapply(specs,
                         function(s) is.character(s) && length(s) == 1L,
                         logical(1))
  if (any(looks_legacy)) {
    stop("config$contrasts uses the removed R-expression format, e.g. ",
         "\"pairs(emmeans(fit, ~group | stage))\". Rewrite it as declarative ",
         "records - see the `contrasts:` block in config/config.template.yaml.")
  }

  mats <- list()
  dirs <- list()

  for (spec in specs) {
    id  <- spec$id
    typ <- if (is.null(spec$type)) "pairwise" else spec$type
    fac <- spec$factor
    by  <- spec$by
    rev <- isTRUE(spec$reverse)
    if (!is.null(by) && !nzchar(by)) by <- NULL

    if (is.null(id)  || !nzchar(id))  stop("every contrast needs an 'id'")
    if (!grepl("^[A-Za-z0-9_][A-Za-z0-9_.-]*$", id) || id %in% names(mats)) {
      stop("contrast IDs must be unique safe identifiers")
    }
    allowed <- if (identical(typ, "linear")) c("id", "type", "weights") else
      c("id", "type", "factor", "by", "reverse", "numerator", "denominator")
    if (length(setdiff(names(spec), allowed))) stop("Unknown contrast fields: ", paste(setdiff(names(spec), allowed), collapse = ", "))
    if (identical(typ, "linear")) {
      w <- unlist(spec$weights)
      coefficients <- names(coef(dummy_fit))
      if (!is.numeric(w) || !length(w) || is.null(names(w)) ||
          anyDuplicated(names(w)) || any(!is.finite(w)) || all(w == 0) ||
          any(!names(w) %in% coefficients)) {
        stop("linear contrast '", id, "' needs finite, nonzero named weights matching design coefficients")
      }
      cm <- matrix(0, nrow = 1, ncol = length(coefficients),
                   dimnames = list(id, coefficients))
      cm[1, names(w)] <- w
      mats[[id]] <- cm
      dirs[[id]] <- data.frame(contrast_name = id, spec_id = id,
        factor = NA_character_, by_variable = NA_character_, by_level = NA_character_,
        numerator_level = NA_character_, denominator_level = NA_character_,
        interpretation = paste(paste(names(w), w, sep = " * "), collapse = " + "),
        stringsAsFactors = FALSE)
      next
    }
    if (is.null(fac) || !nzchar(fac)) stop("contrast '", id, "' needs a 'factor'")
    if (!identical(typ, "pairwise")) {
      stop("contrast '", id, "': unsupported type '", typ,
           "' (supported: pairwise, linear)")
    }

    # Validate every referenced variable against BOTH the sample sheet and the
    # fitted model, so an inherited or stale contrast cannot reach the fit.
    for (v in c(fac, by)) {
      if (!v %in% names(sheet)) {
        stop("contrast '", id, "': '", v, "' is not a column in the sample ",
             "sheet (have: ", paste(names(sheet), collapse = ", "), ")")
      }
      if (!v %in% model_vars) {
        stop("contrast '", id, "': '", v, "' is not a term in the model ",
             "formula (have: ", paste(model_vars, collapse = ", "), ")")
      }
      if (nlevels(factor(sheet[[v]])) < 2L) {
        stop("contrast '", id, "': '", v, "' has fewer than 2 levels")
      }
    }

    em <- if (is.null(by)) {
      emmeans::emmeans(dummy_fit, specs = fac)
    } else {
      emmeans::emmeans(dummy_fit, specs = fac, by = by)
    }
    pr <- pairs(em, reverse = rev)   # emmeans S3 method; reverse sets direction

    cm   <- pr@linfct
    grid <- as.data.frame(pr@grid)
    labs <- as.character(grid$contrast)

    # emmeans labels pairwise contrasts "<numerator> - <denominator>".
    halves <- strsplit(labs, " - ", fixed = TRUE)
    bad <- lengths(halves) != 2L
    if (any(bad)) {
      stop("contrast '", id, "': cannot parse direction from emmeans label(s): ",
           paste(labs[bad], collapse = "; "),
           " (a factor level probably contains ' - ')")
    }
    numer <- vapply(halves, `[`, character(1), 1L)
    denom <- vapply(halves, `[`, character(1), 2L)

    if (!is.null(spec$numerator) || !is.null(spec$denominator)) {
      n <- spec$numerator; d <- spec$denominator
      if (is.null(n) || is.null(d) || length(n) != 1L || length(d) != 1L ||
          identical(n, d) || !all(c(n, d) %in% as.character(sheet[[fac]])))
        stop("Explicit pair needs two distinct observed numerator/denominator levels")
      if (!is.null(spec$reverse)) stop("Use numerator/denominator or reverse, not both")
      selected <- (numer == n & denom == d) | (numer == d & denom == n)
      if (!any(selected)) stop("Requested pair is not estimable")
      sign <- ifelse(numer[selected] == n, 1, -1)
      cm <- cm[selected, , drop = FALSE] * sign
      grid <- grid[selected, , drop = FALSE]
      numer <- rep(n, nrow(cm)); denom <- rep(d, nrow(cm))
    }

    safe <- function(x) gsub("[^A-Za-z0-9._]+", ".", x)
    nm <- paste0(id, "__", safe(numer), "_vs_", safe(denom))
    if (!is.null(by)) {
      nm <- paste0(nm, "__", safe(by), ".", safe(as.character(grid[[by]])))
    }
    rownames(cm) <- nm

    mats[[id]] <- cm
    dirs[[id]] <- data.frame(
      contrast_name     = nm,
      spec_id           = id,
      factor            = fac,
      by_variable       = if (is.null(by)) NA_character_ else by,
      by_level          = if (is.null(by)) NA_character_ else as.character(grid[[by]]),
      numerator_level   = numer,
      denominator_level = denom,
      interpretation    = paste0("positive logFC = higher in '", numer,
                                 "' than in '", denom, "'"),
      stringsAsFactors  = FALSE
    )
  }

  cm_all <- do.call(rbind, mats)
  dup <- duplicated(rownames(cm_all))
  if (any(dup)) {
    stop("duplicate contrast names generated: ",
         paste(unique(rownames(cm_all)[dup]), collapse = ", "))
  }
  list(matrix = cm_all, directions = do.call(rbind, dirs))
}

#!/usr/bin/env Rscript
# 03_de.R — Adaptive differential expression.
#
# Reads aggregate.json to decide:
#   - which backend (limma-voom, or dream when the model has random effects)
#   - whether to auto-add a batch covariate (PCA-based)
# Validates that the design is estimable and stops if it is not; there is no
# automatic fallback to a simpler test.
# Then fits the model from config$model, generates contrasts via emmeans
# from config$contrasts, and writes one DE table per contrast.
#
# Inputs (named):  counts, metrics, sheet
# Outputs (named): done (flag), metrics
# Params:          model_fixed, model_random, primary

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(variancePartition)
  library(BiocParallel)
  library(emmeans)
  library(dplyr)
  library(readr)
  library(jsonlite)
  library(tibble)
})

counts_path  <- snakemake@input$counts
metrics_in   <- snakemake@input$metrics
sheet_path   <- snakemake@input$sheet
out_done     <- snakemake@output$done
out_metrics  <- snakemake@output$metrics
out_manifest <- snakemake@output$manifest
cfg          <- snakemake@config
thr          <- cfg$thresholds

model_fixed   <- snakemake@params$model_fixed
model_random  <- snakemake@params$model_random
primary       <- snakemake@params$primary

agg <- fromJSON(metrics_in)
de_backend <- agg$de_backend
message("DE backend: ", de_backend)

# --- Load data ---------------------------------------------------------
counts <- read_tsv(counts_path, show_col_types = FALSE) %>%
  as.data.frame() %>% tibble::column_to_rownames("gene_id") %>% as.matrix()
sheet  <- read_tsv(sheet_path, comment = "#", show_col_types = FALSE)
sheet  <- sheet[match(colnames(counts), sheet$sample_id), ]
stopifnot(all(sheet$sample_id == colnames(counts)))

# Convert categorical covariates
for (nm in setdiff(names(sheet), c("sample_id", "fastq_url"))) {
  if (is.character(sheet[[nm]])) sheet[[nm]] <- factor(sheet[[nm]])
}

# Normalize + filter
d0 <- DGEList(counts)
d0 <- calcNormFactors(d0)
design_vars <- all.vars(as.formula(model_fixed))
mm_for_filter <- model.matrix(as.formula(model_fixed), data = as.data.frame(sheet))
keep <- filterByExpr(d0, mm_for_filter)
d <- d0[keep, ]
message("Kept ", nrow(d), " / ", nrow(d0), " genes after filterByExpr")

# --- Batch detection gate ---------------------------------------------
batch_col <- "batch"
auto_batch_added <- FALSE
if (batch_col %in% names(sheet) && length(unique(sheet[[batch_col]])) > 1) {
  logcpm <- cpm(d, log = TRUE, prior.count = 2)
  pcs <- prcomp(t(logcpm), scale. = TRUE)
  var_by_batch <- sapply(1:2, function(i) {
    summary(lm(pcs$x[, i] ~ sheet[[batch_col]]))$r.squared
  })
  message(sprintf("PC1 var explained by batch: %.2f  PC2: %.2f",
                  var_by_batch[1], var_by_batch[2]))
  if (max(var_by_batch) > thr$batch_correction$pc_var_fraction_trigger) {
    # Not confounded with primary factor?
    tab <- table(sheet[[batch_col]], sheet[[primary]])
    if (min(rowSums(tab > 0)) > 1) {
      model_fixed <- gsub("^~\\s*", "~ batch + ", model_fixed)
      auto_batch_added <- TRUE
      cat(sprintf("[%s] BATCH_ADD: auto-added 'batch' covariate (max R^2=%.2f)\n",
                  format(Sys.time(), "%FT%T"), max(var_by_batch)),
          file = "gates/decisions.log", append = TRUE)
      message("Auto-added batch covariate. New formula: ", model_fixed)
    } else {
      message("Batch is confounded with ", primary, " — not auto-adding")
    }
  }
}

# --- Fit + contrasts --------------------------------------------------
de_dir <- file.path(cfg$project$output_dir, "DE_Results")
dir.create(de_dir, recursive = TRUE, showWarnings = FALSE)

# Clear stale result tables before writing. A rerun with renamed contrasts
# would otherwise leave the previous run's tables in place, and anything that
# discovers results by scanning the directory would then mix obsolete
# contrasts in with current ones.
stale <- list.files(de_dir, pattern = "_DE_(analysis|meaningful)\\.tsv$",
                    full.names = TRUE)
if (length(stale) > 0) {
  file.remove(stale)
  message("Removed ", length(stale), " stale DE table(s) from a previous run")
}

form_fixed <- as.formula(model_fixed)
mm <- model.matrix(form_fixed, data = as.data.frame(sheet))

# --- Design validation -------------------------------------------------
# Stop rather than substitute. A design that cannot support the requested
# model is a study-design problem; quietly falling back to a simpler test
# would answer a different question and drop the interactions, blocking and
# random effects the plan asked for.
# CT02: "Record scaling; do not add a second transcript-length correction."
# The count matrix arrives already carrying whatever correction 02_aggregate.R
# applied. Reading that record, rather than re-deriving it from seq_type, is
# what stops a backend from adding offsets on top of abundance-derived counts.
assert_no_double_correction <- function(prov) {
  if (is.null(prov)) {
    return(list(ok = TRUE, applied = NA, reason = "no provenance record found"))
  }
  if (isTRUE(prov$further_length_correction_permitted)) {
    return(list(ok = TRUE, applied = isTRUE(prov$length_correction_applied),
                reason = "provenance permits a further correction"))
  }
  list(ok = FALSE,
       applied = isTRUE(prov$length_correction_applied),
       reason = paste0("counts were produced with countsFromAbundance='",
                       prov$counts_from_abundance,
                       "'; a further transcript-length correction is refused: ",
                       prov$rationale))
}

count_prov <- NULL
prov_file <- file.path(dirname(counts_path), "counts_provenance.json")
if (file.exists(prov_file)) {
  count_prov <- fromJSON(prov_file)
  guard <- assert_no_double_correction(count_prov)
  message("Count provenance: countsFromAbundance=",
          count_prov$counts_from_abundance,
          ", length correction already applied=", guard$applied,
          ", further correction permitted=",
          isTRUE(count_prov$further_length_correction_permitted))
} else {
  message("No counts_provenance.json beside ", counts_path,
          "; count semantics are unrecorded and no offsets will be added")
}

source(snakemake@params$design_lib)   # validate_design()
design <- validate_design(
  mm, agg$n_bio_replicates_min, primary,
  if (is.null(agg$biological_unit)) "sample" else agg$biological_unit)
message("Design OK: ", nrow(mm), " samples, ", design$rank, " coefficients, ",
        design$residual_df, " residual df, ", agg$n_bio_replicates_min,
        " biological replicates in the smallest group")

# --- Declarative contrast construction --------------------------------
# Contrasts arrive from config as structured records, never as R source text.
# Everything is built with real function calls, so a stale or malformed spec
# produces a clear validation error in preflight instead of an eval() failure
# or a silently wrong comparison.

build_contrasts <- function(dummy_fit, specs, sheet, model_vars) {
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
    if (is.null(fac) || !nzchar(fac)) stop("contrast '", id, "' needs a 'factor'")
    if (!identical(typ, "pairwise")) {
      stop("contrast '", id, "': unsupported type '", typ,
           "' (supported: pairwise)")
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

fit <- NULL
contrast_sig_counts <- list()

# voom / dream. There is no exactTest arm: a design that cannot support the
# requested model stops in the design validation above rather than being
# rerouted to a test that answers a different question.
form_full <- if (de_backend == "dream" && !is.null(model_random) &&
                 nzchar(model_random)) {
  as.formula(paste(model_fixed, "+", model_random))
} else {
  form_fixed
}

if (de_backend == "dream") {
  message("Running variancePartition::dream")
  param <- SnowParam(2, "SOCK", progressbar = FALSE)
  vobj <- voomWithDreamWeights(d, form_full, as.data.frame(sheet),
                               BPPARAM = param)
} else {
  message("Running limma-voom")
  vobj <- voom(d, mm)
}

# Build contrasts from the declarative specs. The dummy lm supplies the
# coefficient structure emmeans needs; its response is irrelevant, and its
# column order matches `mm` because both use the same fixed-effects formula.
z <- rnorm(nrow(mm))
dummy <- lm(as.formula(paste0("z ", model_fixed)), data = as.data.frame(sheet))
built <- build_contrasts(dummy, cfg$contrasts, sheet, all.vars(form_fixed))
contrast_matrix <- built$matrix

# Record the resolved direction of every contrast alongside the results, so
# the sign of each logFC is interpretable without re-deriving it.
write_tsv(built$directions,
          file.path(cfg$project$output_dir, "DE_Results",
                    "contrast_directions.tsv"))
message("Built ", nrow(contrast_matrix), " contrasts from ",
        length(cfg$contrasts), " spec(s); directions written to ",
        "DE_Results/contrast_directions.tsv")

if (de_backend == "dream") {
  fit <- dream(vobj, form_full, as.data.frame(sheet), t(contrast_matrix))
  fit <- eBayes(fit)
} else {
  fit <- lmFit(vobj, mm)
  fit2 <- contrasts.fit(fit, t(contrast_matrix))
  fit <- eBayes(fit2)
}

meaningful_thr <- cfg$thresholds$de_meaningful
manifest_rows <- list()
for (nm in rownames(contrast_matrix)) {
  tt <- topTable(fit, coef = nm, n = Inf, sort.by = "P")
  tt$gene_id <- rownames(tt)
  write_tsv(tt, file.path(de_dir, paste0(nm, "_DE_analysis.tsv")))
  contrast_sig_counts[[nm]] <- sum(tt$adj.P.Val < 0.05, na.rm = TRUE)

  # The manifest, not the directory listing, is what downstream stages read.
  manifest_rows[[nm]] <- data.frame(
    contrast_name    = nm,
    analysis_table   = file.path("DE_Results", paste0(nm, "_DE_analysis.tsv")),
    meaningful_table = if (is.null(meaningful_thr)) NA_character_ else
                       file.path("DE_Results", paste0(nm, "_DE_meaningful.tsv")),
    statistic        = "t",
    n_genes_tested   = nrow(tt),
    n_significant    = contrast_sig_counts[[nm]],
    stringsAsFactors = FALSE
  )

  # Handbook section 10.2.3: filter on BOTH FDR and effect size for the
  # "biologically meaningful" gene list. A 1.05x change with q=0.001 is
  # statistically real but biologically uninteresting.
  if (!is.null(meaningful_thr)) {
    meaningful <- tt[
      !is.na(tt$adj.P.Val) &
      tt$adj.P.Val < meaningful_thr$adjp_max &
      abs(tt$logFC) >= meaningful_thr$log2fc_min,
    ]
    write_tsv(meaningful, file.path(de_dir, paste0(nm, "_DE_meaningful.tsv")))
  }
}

# --- Result manifest ---------------------------------------------------
# Downstream stages read this, never the directory listing, so a renamed or
# removed contrast cannot resurface as an obsolete result table.
write_tsv(do.call(rbind, manifest_rows), out_manifest)
message("Wrote ", out_manifest, " (", length(manifest_rows), " contrasts)")

# --- Sanity check contrasts -------------------------------------------
too_many <- names(contrast_sig_counts)[
  sapply(contrast_sig_counts, function(x) x > thr$de_sanity$sig_genes_too_many)]
zero_sig <- names(contrast_sig_counts)[sapply(contrast_sig_counts, `==`, 0)]

if (length(too_many) > 0) {
  cat(sprintf("[%s] DE_WARN: contrasts with > %d sig genes: %s — re-check batch\n",
              format(Sys.time(), "%FT%T"),
              thr$de_sanity$sig_genes_too_many,
              paste(too_many, collapse = ",")),
      file = "gates/decisions.log", append = TRUE)
}

metrics <- list(
  backend            = de_backend,
  auto_batch_added   = auto_batch_added,
  contrast_sig_counts = contrast_sig_counts,
  flagged_too_many    = as.list(too_many),
  flagged_zero_sig    = as.list(zero_sig)
)
write_json(metrics, out_metrics, pretty = TRUE, auto_unbox = TRUE)
file.create(out_done)
message("DE done.")

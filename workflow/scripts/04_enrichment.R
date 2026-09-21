#!/usr/bin/env Rscript
# 04_enrichment.R — GO + KEGG GSEA, with adaptive simplify().
#
# Per-contrast: reads DE_Results/*.tsv, ranks genes, runs gseGO and (if kegg
# code present) gseKEGG. If >N significant terms, auto-applies simplify()
# to collapse redundant GO terms via GOSemSim similarity.
#
# Inputs (named):  de_flag
# Outputs (named): done, metrics
# Params:          orgdb, kegg_code, run_go, run_kegg

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(dplyr)
  library(readr)
  library(jsonlite)
})

cfg     <- snakemake@config
thr     <- cfg$thresholds$enrichment
out_done    <- snakemake@output$done
out_metrics <- snakemake@output$metrics
out_status  <- snakemake@output$status

orgdb     <- snakemake@params$orgdb
kegg_code <- snakemake@params$kegg_code
run_go    <- isTRUE(snakemake@params$run_go)
run_kegg  <- isTRUE(snakemake@params$run_kegg)

out_dir   <- cfg$project$output_dir
go_dir    <- file.path(out_dir, "GO_results")
kegg_dir  <- file.path(out_dir, "KEGG_results")
dir.create(go_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(kegg_dir, recursive = TRUE, showWarnings = FALSE)

# --- Optional annotation for NCBI ID mapping ---------------------------
anno <- NULL
anno_path <- cfg$reference$annotation_tsv$path
if (!is.null(anno_path) && file.exists(anno_path)) {
  anno <- read_tsv(anno_path, show_col_types = FALSE)
  anno$NCBI <- gsub("GeneID:", "", anno$DB.xref, fixed = TRUE)
}

# --- Skip semantics ----------------------------------------------------
# "skipped" (switched off, or a prerequisite is absent), "failed" (errored),
# and "empty" (ran fine, found nothing) are three different outcomes. The
# status table keeps them apart instead of collapsing all three into a flag.
orgdb_strategy <- if (is.null(cfg$orgdb$strategy)) "auto" else cfg$orgdb$strategy

# Empty string means "do not skip". R's lazy argument evaluation means
# orgdb_installed is never forced when GO is already switched off, so this
# does not load the namespace just to decide not to use it.
go_skip_because <- function(run_go, orgdb_strategy, orgdb_installed, orgdb_name) {
  if (!isTRUE(run_go))                   return("downstream.run_go is false")
  if (identical(orgdb_strategy, "skip")) return("orgdb.strategy is 'skip'")
  if (!isTRUE(orgdb_installed))          return(paste0("OrgDb '", orgdb_name,
                                                       "' is not installed"))
  ""
}

kegg_skip_because <- function(run_kegg, kegg_code) {
  if (!isTRUE(run_kegg))                        return("downstream.run_kegg is false")
  if (is.null(kegg_code) || !nzchar(kegg_code)) return("organism.kegg_code is not set")
  ""
}

go_skip_reason <- go_skip_because(run_go, orgdb_strategy,
                                  requireNamespace(orgdb, quietly = TRUE), orgdb)
run_go <- run_go && !nzchar(go_skip_reason)

kegg_skip_reason <- kegg_skip_because(run_kegg, kegg_code)
run_kegg <- run_kegg && !nzchar(kegg_skip_reason)

orgdb_loaded <- FALSE
if (run_go) {
  library(orgdb, character.only = TRUE)
  orgdb_loaded <- TRUE
}
if (nzchar(go_skip_reason))   message("GO skipped: ", go_skip_reason)
if (nzchar(kegg_skip_reason)) message("KEGG skipped: ", kegg_skip_reason)

# --- Contrasts come from the DE manifest, never a directory listing ----
manifest <- read_tsv(snakemake@input$manifest, show_col_types = FALSE)
if (nrow(manifest) == 0) {
  stop("DE manifest lists no contrasts: ", snakemake@input$manifest)
}

# Identifier namespace is a property of the run, not of a contrast: either an
# annotation table maps gene_id -> Entrez, or gene_id must already be Entrez.
id_namespace <- if (is.null(anno)) {
  "entrez_assumed_from_gene_id"
} else {
  "entrez_via_annotation_tsv"
}
min_id_mapping_rate <- if (is.null(thr$min_id_mapping_rate)) 0.30 else
                       thr$min_id_mapping_rate
message("Identifier namespace: ", id_namespace,
        " (minimum usable mapping rate ", 100 * min_id_mapping_rate, "%)")

status_rows <- list()
record <- function(contrast, method, status, reason = "", n_terms = NA_integer_) {
  status_rows[[length(status_rows) + 1L]] <<- data.frame(
    contrast = contrast, method = method, status = status,
    reason = reason, n_terms = n_terms, stringsAsFactors = FALSE)
}

for (i in seq_len(nrow(manifest))) {
  contrast <- manifest$contrast_name[i]
  path     <- file.path(out_dir, manifest$analysis_table[i])
  rank_col <- manifest$statistic[i]
  message("Enrichment for: ", contrast)

  if (!file.exists(path)) {
    msg <- paste0("DE table listed in the manifest is missing: ", path)
    record(contrast, "GO", "failed", msg)
    record(contrast, "KEGG", "failed", msg)
    next
  }
  de <- read_tsv(path, show_col_types = FALSE)

  # The ranking statistic must exist. A table without it used to be skipped
  # silently, which looks identical to "analysed and found nothing".
  if (!rank_col %in% names(de)) {
    msg <- paste0("ranking column '", rank_col, "' absent from ", basename(path))
    record(contrast, "GO", "failed", msg)
    record(contrast, "KEGG", "failed", msg)
    next
  }

  n_input <- nrow(de)
  if (!is.null(anno)) {
    de <- de %>% left_join(anno %>% select(Gene.stable.ID, NCBI),
                           by = c("gene_id" = "Gene.stable.ID"))
  } else {
    # Without an annotation table the only way gene_id is usable directly is
    # if it already IS an Entrez ID. Test that instead of assuming it: the
    # previous behaviour silently relabelled gene symbols as Entrez IDs, which
    # enriches against whatever those numbers happen to mean.
    de$NCBI <- ifelse(grepl("^[0-9]+$", de$gene_id), de$gene_id, NA_character_)
  }
  de <- de %>% filter(!is.na(NCBI), !duplicated(NCBI))

  coverage <- if (n_input > 0) nrow(de) / n_input else 0
  if (coverage < min_id_mapping_rate) {
    msg <- sprintf(
      "only %.1f%% of %d genes mapped to usable IDs (minimum %.0f%%, namespace: %s)",
      100 * coverage, n_input, 100 * min_id_mapping_rate, id_namespace)
    record(contrast, "GO", "failed", msg)
    record(contrast, "KEGG", "failed", msg)
    next
  }
  message(sprintf("  %d/%d genes mapped (%.1f%%)", nrow(de), n_input,
                  100 * coverage))

  # --- GO GSEA (non-directional, ranked by |statistic|) ---------------
  if (!run_go) {
    record(contrast, "GO", "skipped", go_skip_reason)
  } else {
    gl <- abs(de[[rank_col]]); names(gl) <- de$NCBI
    gl <- sort(gl, decreasing = TRUE)
    res <- tryCatch(
      gseGO(geneList = gl, ont = "all", OrgDb = orgdb,
            pvalueCutoff = 1, scoreType = "pos", verbose = FALSE),
      error = function(e) e)
    if (inherits(res, "error")) {
      record(contrast, "GO", "failed", conditionMessage(res))
    } else {
      tab <- as.data.frame(res)
      n_sig <- sum(tab$p.adjust < 0.05, na.rm = TRUE)
      if (n_sig > thr$sig_terms_simplify_trigger) {
        res <- simplify(res, cutoff = thr$go_simplify_cutoff,
                        by = thr$go_simplify_by,
                        select_fun = get(thr$go_simplify_select_fun))
        tab <- as.data.frame(res)
        cat(sprintf("[%s] GO_SIMPLIFY: %s had %d sig terms, collapsed\n",
                    format(Sys.time(), "%FT%T"), contrast, n_sig),
            file = "gates/decisions.log", append = TRUE)
      }
      write_tsv(tab, file.path(go_dir, paste0(contrast, "_GO.tsv")))
      record(contrast, "GO",
             if (nrow(tab) == 0) "empty" else "succeeded", "", nrow(tab))
    }
  }

  # --- KEGG GSEA (directional, signed statistic) ----------------------
  if (!run_kegg) {
    record(contrast, "KEGG", "skipped", kegg_skip_reason)
  } else {
    gl <- de[[rank_col]]; names(gl) <- de$NCBI
    gl <- sort(gl, decreasing = TRUE)
    res <- tryCatch(
      gseKEGG(geneList = gl, organism = kegg_code,
              pvalueCutoff = 1, verbose = FALSE),
      error = function(e) e)
    if (inherits(res, "error")) {
      record(contrast, "KEGG", "failed", conditionMessage(res))
    } else {
      tab <- as.data.frame(res)
      write_tsv(tab, file.path(kegg_dir, paste0(contrast, "_KEGG.tsv")))
      record(contrast, "KEGG",
             if (nrow(tab) == 0) "empty" else "succeeded", "", nrow(tab))
    }
  }
}

status <- do.call(rbind, status_rows)
write_tsv(status, out_status)

n_of <- function(s) sum(status$status == s)
metrics <- list(
  orgdb_loaded     = orgdb_loaded,
  id_namespace     = id_namespace,
  min_id_mapping_rate = min_id_mapping_rate,
  go_run           = run_go,
  kegg_run         = run_kegg,
  go_skip_reason   = go_skip_reason,
  kegg_skip_reason = kegg_skip_reason,
  n_succeeded      = n_of("succeeded"),
  n_empty          = n_of("empty"),
  n_failed         = n_of("failed"),
  n_skipped        = n_of("skipped")
)
write_json(metrics, out_metrics, pretty = TRUE, auto_unbox = TRUE)

failed <- status[status$status == "failed", , drop = FALSE]
for (i in seq_len(nrow(failed))) {
  cat(sprintf("[%s] ENRICHMENT_FAILED: %s/%s - %s\n",
              format(Sys.time(), "%FT%T"),
              failed$contrast[i], failed$method[i], failed$reason[i]),
      file = "gates/decisions.log", append = TRUE)
}

# A run where every attempt errored is a failed stage, not a success with a
# done flag. Skipped-only runs are fine: nothing was attempted.
attempted <- status[status$status != "skipped", , drop = FALSE]
if (nrow(attempted) > 0 && all(attempted$status == "failed")) {
  stop("every attempted enrichment failed; see ", out_status)
}

file.create(out_done)
message("Enrichment: ", n_of("succeeded"), " succeeded, ", n_of("empty"),
        " empty, ", n_of("failed"), " failed, ", n_of("skipped"), " skipped")

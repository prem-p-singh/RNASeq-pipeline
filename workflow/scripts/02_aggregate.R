#!/usr/bin/env Rscript
# 02_aggregate.R — tximport: Salmon quants -> gene-level count matrix.
# Invoked by Snakemake. Access via the `snakemake` object.
#
# Inputs (named):  quants (list), gtf, sheet
# Outputs (named): counts, metrics
# Config paths:    reference.annotation_tsv.path (optional gene-info merge)

suppressPackageStartupMessages({
  library(tximport)
  library(GenomicFeatures)
  library(dplyr)
  library(readr)
  library(jsonlite)
})

quant_files <- unlist(snakemake@input$quants)
gtf_path    <- snakemake@input$gtf
sheet_path  <- snakemake@input$sheet
out_counts  <- snakemake@output$counts
out_metrics <- snakemake@output$metrics
out_disposition <- snakemake@output$disposition

cfg <- snakemake@config

# --- Build tx2gene from GTF --------------------------------------------
# makeTxDbFromGFF moved from GenomicFeatures to txdbmaker in Bioconductor 3.19.
# Pick whichever the installed release provides so this works on both.
make_txdb <- if (requireNamespace("txdbmaker", quietly = TRUE)) {
  txdbmaker::makeTxDbFromGFF
} else {
  GenomicFeatures::makeTxDbFromGFF
}

message("Building tx2gene from ", gtf_path)
txdb <- make_txdb(gtf_path, format = "gtf")

# Qualify these explicitly: dplyr is attached after GenomicFeatures and its
# select()/keys() would otherwise mask the AnnotationDbi methods, which have
# no TxDb support and fail at runtime.
k <- AnnotationDbi::keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, keys = k,
                                 columns = "GENEID", keytype = "TXNAME")

# --- Read sample sheet -------------------------------------------------
sheet <- read_tsv(sheet_path, comment = "#", show_col_types = FALSE)
stopifnot("sample_id" %in% names(sheet))

# --- Match quant files to samples --------------------------------------
# quant_files come as "results/.../quant/{sample}/quant.sf"
sample_ids <- basename(dirname(quant_files))
names(quant_files) <- sample_ids
quant_files <- quant_files[match(sheet$sample_id, names(quant_files))]
quant_files <- quant_files[!is.na(quant_files)]
if (length(quant_files) == 0) stop("No quant files matched sample sheet")

# Assay-specific count handling. In full-length RNA-seq a shift in isoform
# usage changes the effective length of a gene, which confounds gene-level DE.
# tximport's "lengthScaledTPM" returns bias-corrected counts that limma-voom
# and edgeR can consume directly, with no separate offset matrix to carry
# across the script boundary.
# TAGseq is the documented exception: 3'-tag counts are one read per molecule
# and do not scale with transcript length, so applying the correction would
# distort them.
seq_type <- cfg$samples$seq_type
counts_from_abundance <- if (identical(seq_type, "tagseq")) "no" else "lengthScaledTPM"

message("Importing ", length(quant_files), " quant.sf files ",
        "(seq_type=", seq_type, ", countsFromAbundance=",
        counts_from_abundance, ")")
txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene,
                countsFromAbundance = counts_from_abundance,
                ignoreAfterBar = TRUE)

counts <- round(txi$counts)

# CT01/CT02: record the count semantics with the counts, so a downstream stage
# can tell what was already done rather than inferring it from seq_type. CT02:
# "Record scaling; do not add a second transcript-length correction." The
# `length_correction_applied` flag is what 03_de.R checks before it would add
# offsets of its own.
count_provenance <- list(
  schema_version            = 1L,
  producer                  = "02_aggregate.R",
  seq_type                  = seq_type,
  counts_from_abundance     = counts_from_abundance,
  length_correction_applied = !identical(counts_from_abundance, "no"),
  # CT03: 3-prime counts are molecule-like and must never receive the
  # conventional full-length correction, here or later.
  further_length_correction_permitted = FALSE,
  rationale = if (identical(counts_from_abundance, "no"))
      paste("3-prime tag counts do not scale with transcript length;",
            "no length correction was applied and none may be added (CT03)")
    else
      paste("abundance-derived counts already carry the transcript-length",
            "correction (lengthScaledTPM); adding another would double it (CT02)")
)

# --- Optional gene-info merge -----------------------------------------
anno_path <- cfg$reference$annotation_tsv$path
if (!is.null(anno_path) && file.exists(anno_path)) {
  anno <- read_tsv(anno_path, show_col_types = FALSE)
  message("Merging annotation from ", anno_path, " (", nrow(anno), " rows)")
  # Keep just the count matrix; annotation applied downstream in 03_de.R
}

# --- Sample disposition ------------------------------------------------
# Stage 1 records per-sample QC but has no authority to act on it. The decision
# is made here, written down with a reason, and enforced by dropping excluded
# samples from the count matrix, so DE, enrichment and WGCNA all see the same
# cohort and the design is re-validated against it in 03_de.R.
#
# Exclusion is an explicit policy: when sample_policy is absent from the config
# the safe default is to keep every sample and only record the flag.
sample_disposition <- function(qc, min_map, min_reads, policy_on) {
  qc$reads_mapped <- qc$reads_processed * qc$mapping_rate

  # QC08: "Required QC metric unavailable | Mark missing and block dependent
  # automatic decision if essential; never substitute zero."
  #
  # This previously read `!is.na(x) & x < threshold`, so a sample whose metric
  # was missing produced FALSE, an empty flag, and include = TRUE. Absent was
  # silently treated as passed, which is the one reading QC08 forbids. Missing
  # is now its own outcome, and the automatic decision is withheld rather than
  # guessed in either direction.
  unavailable <- is.na(qc$mapping_rate) | is.na(qc$reads_processed)

  low_map  <- !unavailable & qc$mapping_rate < min_map
  low_read <- !unavailable & qc$reads_mapped < min_reads

  flag <- rep("", nrow(qc))
  flag[low_map] <- sprintf("mapping rate %.3f < %.2f",
                           qc$mapping_rate[low_map], min_map)
  only_read <- low_read & !low_map
  flag[only_read] <- sprintf("mapped reads %.0f < %.0f",
                             qc$reads_mapped[only_read], min_reads)
  flag[unavailable] <- vapply(which(unavailable), function(i) {
    absent <- c("mapping_rate", "num_reads_processed")[
      c(is.na(qc$mapping_rate[i]), is.na(qc$reads_processed[i]))]
    paste0("QC metric unavailable: ", paste(absent, collapse = ", "))
  }, character(1))

  qc$flag <- flag
  qc$status <- ifelse(unavailable, "unavailable",
                      ifelse(flag == "", "ok", "flagged"))
  # include is NA where no automatic decision may be made. Downstream must
  # resolve that explicitly; it is not a silent TRUE.
  qc$include <- if (isTRUE(policy_on)) {
    ifelse(unavailable, NA, flag == "")
  } else {
    ifelse(unavailable, NA, TRUE)
  }
  qc
}

qc_rows <- lapply(quant_files, function(q) {
  p <- file.path(dirname(q), "metrics.json")
  if (!file.exists(p)) return(NULL)
  m <- fromJSON(p)
  data.frame(sample_id       = m$sample,
             mapping_rate    = as.numeric(m$mapping_rate),
             reads_processed = as.numeric(m$num_reads_processed),
             stringsAsFactors = FALSE)
})
qc <- do.call(rbind, qc_rows[!vapply(qc_rows, is.null, logical(1))])
if (is.null(qc)) stop("no per-sample metrics.json found beside the quant files")

min_reads <- if (identical(seq_type, "tagseq")) {
  cfg$thresholds$sample_qc$min_reads_on_genes_tagseq
} else {
  cfg$thresholds$sample_qc$min_reads_on_genes_rnaseq
}
# Reconcile the cohort against its QC records before any decision is made.
# qc_rows drops a sample whose metrics.json is absent, so without this a sample
# could be analysed with no QC record at all: not flagged, not excluded, simply
# unexamined. That is the same failure as QC08 in a different place.
no_metrics <- setdiff(colnames(counts), qc$sample_id)
if (length(no_metrics) > 0) {
  qc <- rbind(qc, data.frame(
    sample_id       = no_metrics,
    mapping_rate    = NA_real_,
    reads_processed = NA_real_,
    stringsAsFactors = FALSE))
  message("No metrics.json for ", length(no_metrics), " sample(s): ",
          paste(no_metrics, collapse = ", "),
          " - recorded as unavailable, not as passing")
}

policy_on <- isTRUE(cfg$sample_policy$exclude_failing_qc)
disposition <- sample_disposition(qc,
                                  cfg$thresholds$sample_qc$mapping_rate_min,
                                  min_reads,
                                  policy_on)
write_tsv(disposition, out_disposition)

# QC08: block the dependent automatic decision rather than guessing it. With
# the exclusion policy on, inclusion depends on metrics that are missing, so
# the cohort cannot be settled automatically.
undecided <- disposition$sample_id[is.na(disposition$include)]
if (policy_on && length(undecided) > 0) {
  stop("sample_policy.exclude_failing_qc is on, but QC metrics are unavailable ",
       "for ", length(undecided), " sample(s): ",
       paste(undecided, collapse = ", "),
       ". Inclusion cannot be decided automatically. Re-run Stage 1 for these ",
       "samples, or set sample_policy.exclude_failing_qc to false to analyse ",
       "the full cohort with the gap recorded in ", out_disposition, ".")
}

excluded <- disposition$sample_id[!is.na(disposition$include) &
                                  !disposition$include]
if (length(excluded) > 0) {
  counts <- counts[, !(colnames(counts) %in% excluded), drop = FALSE]
  for (i in which(!disposition$include)) {
    cat(sprintf("[%s] SAMPLE_EXCLUDED: %s (%s)\n",
                format(Sys.time(), "%FT%T"),
                disposition$sample_id[i], disposition$flag[i]),
        file = "gates/decisions.log", append = TRUE)
  }
}
if (ncol(counts) == 0) stop("every sample was excluded by the QC policy")
message("Cohort: ", ncol(counts), " samples retained, ",
        length(excluded), " excluded (see ", out_disposition, ")")

# --- Count provenance and gene lengths (CT01) --------------------------
prov_path <- file.path(dirname(out_counts), "counts_provenance.json")
write_json(count_provenance, prov_path, pretty = TRUE, auto_unbox = TRUE)
message("Wrote ", prov_path, " (countsFromAbundance=", counts_from_abundance, ")")

# CT01 asks that length information be preserved, not just counts. tximport
# returns it per gene; keeping it means a backend that wants offsets has them
# without re-importing.
if (!is.null(txi$length)) {
  len_df <- as.data.frame(txi$length[, colnames(counts), drop = FALSE])
  len_df$gene_id <- rownames(txi$length)
  write_tsv(len_df[, c("gene_id", setdiff(names(len_df), "gene_id"))],
            file.path(dirname(out_counts), "gene_lengths.tsv"))
}

# Preserve the gene-level import object for downstream adapters and provenance.
for (key in c("counts", "abundance", "length")) {
  txi[[key]] <- txi[[key]][, colnames(counts), drop = FALSE]
}
saveRDS(list(import = txi, provenance = count_provenance),
        file.path(dirname(out_counts), "gene_import.rds"))

# --- Write count matrix ------------------------------------------------
out_df <- as.data.frame(counts)
out_df$gene_id <- rownames(out_df)
out_df <- out_df[, c("gene_id", setdiff(names(out_df), "gene_id"))]
write_tsv(out_df, out_counts)
message("Wrote ", out_counts, "  (", nrow(out_df), " genes x ",
        ncol(out_df) - 1, " samples)")

# --- Metrics + decision gates -----------------------------------------
libs <- colSums(counts)
med <- median(libs)
cv  <- sd(libs) / mean(libs)

thr <- cfg$thresholds$aggregate
cv_warn       <- isTRUE(cv > thr$library_size_cv_warn)
small_samples <- names(libs)[libs < med * thr$sample_min_fraction_of_median]

# Pick the DE backend from the declared design only. There is deliberately no
# low-replication rescue: swapping in a simpler test answers a different
# question than the analysis plan asked, and silently drops interactions,
# blocking and random effects. 03_de.R validates that the design is estimable
# and stops when it is not.
primary <- cfg$model$primary_factor
kept    <- sheet[sheet$sample_id %in% colnames(counts), , drop = FALSE]
has_random <- !is.null(cfg$model$random_effects) &&
              nzchar(cfg$model$random_effects)

# Independent biological replicates, not rows. With a random-effects grouping
# variable (e.g. "(1|vine)") one subject contributes several rows, and those
# rows are not independent evidence about the primary factor.
source(snakemake@params$design_lib)   # biological_replicates()
reps <- biological_replicates(kept, primary, cfg$model$random_effects)
de_backend <- if (has_random) "dream" else "limma_voom"

metrics <- list(
  n_samples         = ncol(counts),
  n_genes           = nrow(counts),
  n_excluded        = length(excluded),
  excluded_samples  = as.list(excluded),
  counts_from_abundance = counts_from_abundance,
  library_size      = list(min = min(libs), median = med, max = max(libs)),
  library_size_cv   = cv,
  cv_warn           = cv_warn,
  flagged_samples   = as.list(small_samples),
  n_per_group_min   = reps$rows_per_group,
  n_bio_replicates_min = reps$n,
  biological_unit   = reps$unit,
  has_random_effect = has_random,
  de_backend        = de_backend
)
write_json(metrics, out_metrics, pretty = TRUE, auto_unbox = TRUE)
message("Wrote ", out_metrics)

# --- Log decision -----------------------------------------------------
cat(sprintf("[%s] DE_BACKEND: %s (rows/group=%s, biological replicates/group=%s, unit=%s)\n",
            format(Sys.time(), "%FT%T"), de_backend,
            reps$rows_per_group, reps$n, reps$unit),
    file = "gates/decisions.log", append = TRUE)
if (length(small_samples) > 0) {
  cat(sprintf("[%s] SAMPLE_FLAG: %s have < median/10 library size\n",
              format(Sys.time(), "%FT%T"),
              paste(small_samples, collapse = ",")),
      file = "gates/decisions.log", append = TRUE)
}

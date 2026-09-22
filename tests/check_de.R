#!/usr/bin/env Rscript
# Self-checks for the non-trivial logic in the Stage 2/3 scripts:
#   03_de.R       build_contrasts()
#   _design.R     validate_design(), biological_replicates()
#
# What must not break:
#   - a contrast named "treated_vs_control" must actually mean treated minus
#     control; a sign flip silently inverts every biological conclusion
#   - an unsupported design must STOP, not fall through to a different test
#   - repeated measurements of one subject must not count as replication
#
# Run:  Rscript tests/check_de.R

suppressPackageStartupMessages(library(emmeans))

# 03_de.R is a Snakemake script (top-level snakemake@ references) so it cannot
# be sourced whole. Pull out just the function definitions we are checking.
load_fns <- function(path, names) {
  found <- character(0)
  for (e in parse(file = path)) {
    # is.name() matters: a top-level `out_df$gene_id <- ...` has a call on the
    # left, and as.character() on that returns length 3.
    if (is.call(e) && identical(as.character(e[[1]]), "<-") &&
        is.name(e[[2]]) && as.character(e[[2]]) %in% names) {
      eval(e, envir = globalenv())
      found <- c(found, as.character(e[[2]]))
    }
  }
  missing <- setdiff(names, found)
  if (length(missing)) {
    stop("not found in ", path, ": ", paste(missing, collapse = ", "))
  }
}
# _design.R is a plain sourceable file, so it loads directly.
source(file.path("workflow", "scripts", "_design.R"))

load_fns(file.path("workflow", "scripts", "03_de.R"),
         c("build_contrasts", "assert_no_double_correction"))
load_fns(file.path("workflow", "scripts", "02_aggregate.R"),
         c("sample_disposition"))
load_fns(file.path("workflow", "scripts", "04_enrichment.R"),
         c("go_skip_because", "kegg_skip_because"))

# ======================================================================
# build_contrasts(): direction is verified numerically, not by label
# ======================================================================
sheet <- data.frame(
  sample_id = paste0("S", 1:8),
  treatment = factor(rep(c("control", "treated"), each = 4)),
  stage     = factor(rep(c("early", "late"), times = 4)),
  stringsAsFactors = FALSE
)
model_vars <- c("treatment", "stage")
set.seed(1)
sheet$z <- rnorm(nrow(sheet))
dummy <- lm(z ~ treatment * stage, data = sheet)

# A response where treated is higher by exactly 10, so applying the contrast
# vector to these coefficients must give ~ +10 for "treated - control".
sheet$y <- ifelse(sheet$treatment == "treated", 10, 0) + rnorm(nrow(sheet), sd = 0.01)
real <- lm(y ~ treatment * stage, data = sheet)
estimate <- function(cm) as.vector(cm %*% coef(real))

b <- build_contrasts(dummy,
                     list(list(id = "trt", type = "pairwise",
                               factor = "treatment", by = NULL, reverse = TRUE)),
                     sheet, model_vars)
stopifnot(nrow(b$matrix) == 1)
stopifnot(grepl("treated_vs_control", rownames(b$matrix)[1]))
stopifnot(b$directions$numerator_level == "treated")
stopifnot(b$directions$denominator_level == "control")
stopifnot(abs(estimate(b$matrix) - 10) < 0.5)

# reverse = FALSE must flip the sign
bf <- build_contrasts(dummy,
                      list(list(id = "trt", type = "pairwise",
                                factor = "treatment", by = NULL, reverse = FALSE)),
                      sheet, model_vars)
stopifnot(abs(estimate(bf$matrix) + 10) < 0.5)

# within-stage keeps direction in every stratum
bw <- build_contrasts(dummy,
                      list(list(id = "trt_by_stage", type = "pairwise",
                                factor = "treatment", by = "stage", reverse = TRUE)),
                      sheet, model_vars)
stopifnot(nrow(bw$matrix) == 2)
stopifnot(all(abs(estimate(bw$matrix) - 10) < 0.5))
stopifnot(all(bw$directions$by_variable == "stage"))

# a contrast naming a column that does not exist must stop
bad <- tryCatch({
  build_contrasts(dummy, list(list(id = "stale", type = "pairwise",
                                   factor = "group", by = NULL, reverse = TRUE)),
                  sheet, model_vars); "no error"
}, error = function(e) conditionMessage(e))
stopifnot(grepl("not a column in the sample sheet", bad))

# the removed R-expression format must be rejected, not guessed at
legacy <- tryCatch({
  build_contrasts(dummy, list(between_groups = "pairs(emmeans(fit, ~group))"),
                  sheet, model_vars); "no error"
}, error = function(e) conditionMessage(e))
stopifnot(grepl("removed R-expression format", legacy))

# ======================================================================
# validate_design(): unsupported designs stop instead of being rerouted
# ======================================================================
mm_ok <- model.matrix(~ treatment * stage, data = sheet)   # 8 rows, 4 coefs
d <- validate_design(mm_ok, n_bio_rep = 4, primary = "treatment")
stopifnot(d$rank == 4, d$residual_df == 4)

# perfectly confounded covariate -> rank-deficient -> must stop
conf <- sheet
conf$batch <- factor(ifelse(conf$treatment == "treated", "B2", "B1"))
mm_conf <- model.matrix(~ treatment + batch, data = conf)
rank_err <- tryCatch({
  validate_design(mm_conf, n_bio_rep = 4, primary = "treatment"); "no error"
}, error = function(e) conditionMessage(e))
stopifnot(grepl("rank-deficient", rank_err))

# saturated model (no residual df) -> must stop
small <- sheet[c(1, 5), ]
mm_sat <- model.matrix(~ treatment, data = small)          # 2 rows, 2 coefs
df_err <- tryCatch({
  validate_design(mm_sat, n_bio_rep = 1, primary = "treatment"); "no error"
}, error = function(e) conditionMessage(e))
stopifnot(grepl("residual degrees of freedom", df_err))

# a single biological replicate -> must stop, and name the unit
rep_err <- tryCatch({
  validate_design(mm_ok, n_bio_rep = 1, primary = "treatment", unit = "vine")
  "no error"
}, error = function(e) conditionMessage(e))
stopifnot(grepl("biological replicates", rep_err))
stopifnot(grepl("vine", rep_err))

# ======================================================================
# biological_replicates(): repeated measures must not be counted as
# independent replication
# ======================================================================
# no random effect -> one row is one independent unit
r <- biological_replicates(sheet, "treatment", NULL)
stopifnot(r$n == 4, r$rows_per_group == 4, r$unit == "sample")

# same vine measured twice per group -> 4 rows but only 2 independent units
rm_sheet <- sheet
rm_sheet$vine <- factor(c("v1", "v1", "v2", "v2", "v3", "v3", "v4", "v4"))
r <- biological_replicates(rm_sheet, "treatment", "(1|vine)")
stopifnot(r$n == 2, r$rows_per_group == 4, r$unit == "vine")

# whitespace in the random-effects term must still parse
r <- biological_replicates(rm_sheet, "treatment", "(1 | vine)")
stopifnot(r$n == 2, r$unit == "vine")

# a grouping variable absent from the sheet falls back to rows, not an error
r <- biological_replicates(sheet, "treatment", "(1|donor)")
stopifnot(r$n == 4, r$unit == "donor")

# ======================================================================
# sample_disposition(): QC findings must actually control inclusion
# ======================================================================
qc <- data.frame(
  sample_id       = c("ok", "lowmap", "lowdepth"),
  mapping_rate    = c(0.90,   0.20,     0.90),
  reads_processed = c(3e7,    3e7,      1e6),
  stringsAsFactors = FALSE
)

# policy on -> the two failing samples are dropped, with reasons
d <- sample_disposition(qc, min_map = 0.60, min_reads = 2e7, policy_on = TRUE)
stopifnot(d$include == c(TRUE, FALSE, FALSE))
stopifnot(d$flag[1] == "")
stopifnot(grepl("mapping rate", d$flag[2]))
stopifnot(grepl("mapped reads", d$flag[3]))

# policy off -> same flags recorded, nothing dropped
d <- sample_disposition(qc, min_map = 0.60, min_reads = 2e7, policy_on = FALSE)
stopifnot(all(d$include))
stopifnot(grepl("mapping rate", d$flag[2]))   # still reported

# a sample failing both criteria reports the mapping-rate reason, not two
both <- data.frame(sample_id = "bad", mapping_rate = 0.1,
                   reads_processed = 1e5, stringsAsFactors = FALSE)
d <- sample_disposition(both, min_map = 0.60, min_reads = 2e7, policy_on = TRUE)
stopifnot(!d$include, grepl("mapping rate", d$flag))

# QC08: an absent metric is not a metric that passed.
#
# This previously read `!is.na(x) & x < threshold`, so a missing mapping rate
# gave FALSE, an empty flag and include = TRUE: the sample entered the analysis
# as though it had cleared QC. The automatic decision is now withheld (NA)
# instead of guessed, and the reason names the metric that is missing.
gap <- data.frame(
  sample_id       = c("ok", "no_map", "no_reads"),
  mapping_rate    = c(0.90, NA,       0.90),
  reads_processed = c(3e7,  3e7,      NA),
  stringsAsFactors = FALSE
)
d <- sample_disposition(gap, min_map = 0.60, min_reads = 2e7, policy_on = TRUE)
stopifnot(identical(d$status, c("ok", "unavailable", "unavailable")))
stopifnot(isTRUE(d$include[1]))
stopifnot(is.na(d$include[2]), is.na(d$include[3]))
stopifnot(!isTRUE(d$include[2]))          # must never read as an implicit pass
stopifnot(grepl("unavailable", d$flag[2]), grepl("mapping_rate", d$flag[2]))
stopifnot(grepl("num_reads_processed", d$flag[3]))

# with the policy off no automatic decision is being made, so the gap is
# recorded but nothing is silently asserted about those samples either
d <- sample_disposition(gap, min_map = 0.60, min_reads = 2e7, policy_on = FALSE)
stopifnot(isTRUE(d$include[1]), is.na(d$include[2]))
stopifnot(identical(d$status[2], "unavailable"))

# status separates the three outcomes that used to share an empty flag
d <- sample_disposition(qc, min_map = 0.60, min_reads = 2e7, policy_on = TRUE)
stopifnot(identical(d$status, c("ok", "flagged", "flagged")))

# ======================================================================
# ST14: a cohort that QC has changed must be revalidated before fitting
# ======================================================================
# validate_design is what 03_de.R applies AFTER exclusion, using the replicate
# count recomputed on the retained samples. A design that was estimable for the
# full cohort can stop being estimable once QC removes samples, and that must
# block rather than proceed.
full <- data.frame(
  sample_id = paste0("S", 1:6),
  treatment = factor(rep(c("control", "treated"), each = 3)),
  stringsAsFactors = FALSE
)
mm_full <- model.matrix(~ treatment, data = full)
reps_full <- biological_replicates(full, "treatment", NULL)
stopifnot(reps_full$n == 3)
ok <- validate_design(mm_full, reps_full$n, "treatment", reps_full$unit)
stopifnot(ok$residual_df == 4)

# QC removes two treated samples: 3 vs 1 leaves one replicate in a group
kept <- full[full$sample_id %in% c("S1", "S2", "S3", "S4"), ]
mm_kept <- model.matrix(~ treatment, data = kept)
reps_kept <- biological_replicates(kept, "treatment", NULL)
stopifnot(reps_kept$n == 1)
err <- tryCatch({
  validate_design(mm_kept, reps_kept$n, "treatment", reps_kept$unit)
  ""
}, error = function(e) conditionMessage(e))
stopifnot(grepl("biological replicate", err))

# QC removes an entire level: the contrast has nothing left to compare
one_level <- full[full$treatment == "control", ]
stopifnot(length(unique(one_level$treatment[drop = TRUE])) == 1)
err <- tryCatch({
  mm <- model.matrix(~ treatment, data = droplevels(one_level))
  validate_design(mm, biological_replicates(one_level, "treatment", NULL)$n,
                  "treatment", "sample")
  ""
}, error = function(e) conditionMessage(e))
stopifnot(nzchar(err))

# ======================================================================
# skip semantics: switched-off and unavailable must be distinguishable,
# and neither may be reported as a completed analysis
# ======================================================================
# GO: everything available -> no skip
stopifnot(go_skip_because(TRUE, "auto", TRUE, "org.Xx.eg.db") == "")
# switched off wins over everything else
stopifnot(grepl("run_go is false", go_skip_because(FALSE, "auto", TRUE, "org.Xx.eg.db")))
# orgdb.strategy=skip must actually disable GO (previously it did not)
stopifnot(grepl("strategy is 'skip'", go_skip_because(TRUE, "skip", TRUE, "org.Xx.eg.db")))
# package missing is its own reason, not a silent pass
stopifnot(grepl("is not installed", go_skip_because(TRUE, "auto", FALSE, "org.Xx.eg.db")))
# laziness: when GO is off the installed-check must never be evaluated
stopifnot(go_skip_because(FALSE, "auto", stop("must not be evaluated"),
                          "org.Xx.eg.db") != "")

# KEGG
stopifnot(kegg_skip_because(TRUE, "vvi") == "")
stopifnot(grepl("run_kegg is false", kegg_skip_because(FALSE, "vvi")))
stopifnot(grepl("kegg_code is not set", kegg_skip_because(TRUE, NULL)))
stopifnot(grepl("kegg_code is not set", kegg_skip_because(TRUE, "")))

# ======================================================================
# CT02: a second transcript-length correction must be refused
# ======================================================================
# The count matrix arrives already carrying whatever correction aggregation
# applied. 03_de.R reads that record rather than re-deriving it from seq_type,
# so a backend cannot add offsets on top of abundance-derived counts.
bulk <- list(counts_from_abundance = "lengthScaledTPM",
             length_correction_applied = TRUE,
             further_length_correction_permitted = FALSE,
             rationale = "abundance-derived counts already carry it")
g <- assert_no_double_correction(bulk)
stopifnot(!g$ok, isTRUE(g$applied))
stopifnot(grepl("refused", g$reason), grepl("lengthScaledTPM", g$reason))

# CT03: 3-prime counts were never corrected and still may not be
tag <- list(counts_from_abundance = "no",
            length_correction_applied = FALSE,
            further_length_correction_permitted = FALSE,
            rationale = "3-prime tag counts do not scale with transcript length")
g <- assert_no_double_correction(tag)
stopifnot(!g$ok, !isTRUE(g$applied))
stopifnot(grepl("do not scale", g$reason))

# an absent record is reported as absent, not assumed safe either way
g <- assert_no_double_correction(NULL)
stopifnot(g$ok, is.na(g$applied), grepl("no provenance", g$reason))

# a record that explicitly permits one is honoured
ok <- list(counts_from_abundance = "no", length_correction_applied = FALSE,
           further_length_correction_permitted = TRUE, rationale = "raw counts")
stopifnot(assert_no_double_correction(ok)$ok)

cat("check_de.R: all assertions passed\n")

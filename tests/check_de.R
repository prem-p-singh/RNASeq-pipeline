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
         c("build_contrasts"))
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

cat("check_de.R: all assertions passed\n")

#!/usr/bin/env Rscript
# Real gseGO against a tiny, locally built OrgDb. No network or species database
# download. Verifies signs, recorded output paths, empty results and failures.
args <- commandArgs(FALSE)
repo <- normalizePath(file.path(dirname(sub("^--file=", "", args[grepl("^--file=", args)])), ".."))
packages <- c("AnnotationForge", "clusterProfiler", "GO.db", "jsonlite", "readr", "dplyr")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) { cat("needs", paste(missing, collapse=", "), "\n"); quit(status=77) }
work <- tempfile("enrichment-")
dir.create(work)
old <- setwd(work)
dir.create("lib")
lib <- normalizePath("lib")
.libPaths(c(lib, .libPaths()))
gids <- as.character(1:200)
gene_info <- data.frame(GID=gids, SYMBOL=paste0("G", gids), GENENAME=paste0("Gene ", gids))
go <- data.frame(GID=gids, GO=rep(c("GO:0006281", "GO:0006412", "GO:0006955", "GO:0006810"), each=50), EVIDENCE="IDA")
AnnotationForge::makeOrgPackage(gene_info=gene_info, go=go, version="0.0.1",
  maintainer="Test Author <test@example.org>", author="Test Author", outputDir=work,
  tax_id="9606", genus="Fixture", species="test", goTable="go")
pkg_path <- list.dirs(work, recursive=FALSE, full.names=TRUE)
pkg_path <- pkg_path[grepl("org\\..*\\.eg\\.db$", pkg_path)]
stopifnot(length(pkg_path) == 1)
install.packages(pkg_path, repos=NULL, type="source", lib=lib, quiet=TRUE)
pkg <- basename(pkg_path)
setClass("EnrichmentFixture", representation(config="list", input="list", output="list", params="list"))

run_case <- function(name, n=200, missing_code=FALSE, missing_rank=FALSE) {
  dir.create(name); setwd(file.path(work, name))
  on.exit(setwd(work), add=TRUE)
  dir.create("results"); dir.create("metrics"); dir.create("gates")
  set.seed(7)
  ranks <- seq(8, -8, length.out=n) + rnorm(n, sd=0.01)
  table <- data.frame(gene_id=gids[seq_len(n)], t=ranks)
  if (missing_rank) names(table)[2] <- "wrong"
  readr::write_tsv(table, "results/de.tsv")
  readr::write_tsv(data.frame(contrast_name="treated_vs_control", analysis_table="de.tsv", statistic="t"), "results/de_manifest.tsv")
  smk <- new("EnrichmentFixture",
    config=list(project=list(output_dir="results"), reference=list(annotation_tsv=list(path=NULL)),
      orgdb=list(strategy="use_existing", key_type="GID"), thresholds=list(enrichment=list(min_id_mapping_rate=0.3,
        sig_terms_simplify_trigger=100000, go_simplify_cutoff=0.7, go_simplify_by="p.adjust", go_simplify_select_fun="min"))),
    input=list(manifest="results/de_manifest.tsv"),
    output=list(done="results/enrichment_done.flag", metrics="metrics/enrichment.json", status="results/enrichment_status.tsv"),
    params=list(orgdb=pkg, kegg_code=NULL, run_go=TRUE, run_kegg=missing_code))
  env <- new.env(parent=globalenv()); env$snakemake <- smk
  failed <- tryCatch({sys.source(file.path(repo, "workflow/scripts/04_enrichment.R"), env); FALSE}, error=function(e) { message("Stage error: ", conditionMessage(e)); TRUE })
  stopifnot(identical(failed, missing_rank))
  status <- readr::read_tsv("results/enrichment_status.tsv", show_col_types=FALSE)
  if (missing_rank) {
    stopifnot(!file.exists("results/enrichment_done.flag"), "failed" %in% status$status)
    return(invisible(NULL))
  }
  if (missing_code) stopifnot(status$status[status$method=="KEGG"] == "unavailable")
  row <- status[status$method=="GO", ]
  if (n < 10) {
    stopifnot(row$status == "succeeded_empty")
  } else {
    stopifnot(row$status == "succeeded", file.exists(file.path("results", row$output_table)))
    tab <- readr::read_tsv(file.path("results", row$output_table), show_col_types=FALSE)
    stopifnot(any(tab$NES > 0), any(tab$NES < 0))
  }
}
run_case("success", missing_code=TRUE)
run_case("empty", n=5)
run_case("failure", missing_rank=TRUE)
setwd(old)
unlink(work, recursive=TRUE)
cat("check_enrichment.R: signed GO, output paths, empty, unavailable and failed passed\n")

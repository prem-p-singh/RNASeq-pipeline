#!/usr/bin/env Rscript
# ===========================================================================
# parse_eggnog.R — turn eggNOG-mapper output into AnnotationForge tables
# ---------------------------------------------------------------------------
# eggNOG-mapper writes a TSV file called <prefix>.emapper.annotations with
# one row per query gene. We need two tables for AnnotationForge::makeOrgPackage:
#
#   gene_info  data.frame:   GID, SYMBOL, GENENAME
#   go         data.frame:   GID, GO, EVIDENCE
#
# This file is loaded (sourced) by build_orgdb.R when Tier 3b runs. It is not
# meant to be run on its own.
# ===========================================================================

parse_eggnog_to_tables <- function(annotations_file) {

  # eggNOG-mapper annotations file is TSV with a "#query" header line.
  # Comment lines start with "##"; we read the data line by line.
  lines <- readLines(annotations_file)
  data_lines <- lines[!grepl("^##", lines)]

  # First non-## line is the header, starts with "#query"
  header_line <- data_lines[1]
  header <- strsplit(sub("^#", "", header_line), "\t", fixed = TRUE)[[1]]
  body <- data_lines[-1]

  # Parse the body into a data.frame
  df <- read.table(
    text = body,
    sep = "\t",
    header = FALSE,
    col.names = header,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    fill = TRUE
  )

  # --- gene_info table ----------------------------------------------------
  # Columns we need:
  #   GID       gene ID (the query)
  #   SYMBOL    preferred name (Preferred_name column)
  #   GENENAME  free-text description (Description column)
  gene_info <- data.frame(
    GID      = df$query,
    SYMBOL   = ifelse(df$Preferred_name == "-" | df$Preferred_name == "",
                      df$query, df$Preferred_name),
    GENENAME = ifelse(df$Description == "-" | df$Description == "",
                      "uncharacterized", df$Description),
    stringsAsFactors = FALSE
  )
  # Drop duplicates by GID
  gene_info <- gene_info[!duplicated(gene_info$GID), ]

  # --- go table ----------------------------------------------------------
  # eggNOG's GOs column is comma-separated (e.g. "GO:0003674,GO:0008152").
  # AnnotationForge wants one row per (GID, GO) pair plus an EVIDENCE code.
  go_rows <- list()
  for (i in seq_len(nrow(df))) {
    gid <- df$query[i]
    gos <- df$GOs[i]
    if (is.na(gos) || gos == "-" || gos == "") next
    for (g in strsplit(gos, ",", fixed = TRUE)[[1]]) {
      go_rows[[length(go_rows) + 1]] <- data.frame(
        GID      = gid,
        GO       = trimws(g),
        EVIDENCE = "IEA",   # IEA = Inferred from Electronic Annotation
        stringsAsFactors = FALSE
      )
    }
  }
  go <- if (length(go_rows) > 0) do.call(rbind, go_rows) else
        data.frame(GID = character(), GO = character(), EVIDENCE = character())

  list(gene_info = gene_info, go = go)
}

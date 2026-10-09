#!/usr/bin/env bash
# Interaction (DESeq2), subject-blocked pairing (edgeR QL) and a continuous
# covariate (limma-voom) through the full DE stage, each against a direct fit.
set -euo pipefail
REPO=$(cd "$(dirname "$0")/.." && pwd)
. "$REPO/tests/_gate.sh"
gate_runtime r
TASK_TMP=$(mktemp -d)
trap 'rm -rf "$TASK_TMP"' EXIT

# Kept samples after fixture QC: control S01 S02 S04, treated S05 S06 S08.
for case in interaction paired continuous; do
    proj="$TASK_TMP/$case"
    python3 "$REPO/tests/fixtures/tiny_project.py" build "$proj" "$REPO" rnaseq
    python3 - "$proj" "$case" <<'PY'
import csv, sys, yaml
from pathlib import Path
p, case = Path(sys.argv[1]), sys.argv[2]
cfg_path = p / 'config/config.yaml'; c = yaml.safe_load(cfg_path.read_text())
sheet = p / 'config/samples.tsv'
with sheet.open() as f: rows = list(csv.DictReader(f, delimiter='\t'))
extra = {'S01': ('F', 'D1', 31), 'S02': ('M', 'D2', 44), 'S03': ('F', 'D4', 52), 'S04': ('F', 'D3', 29),
         'S05': ('F', 'D1', 37), 'S06': ('M', 'D2', 58), 'S07': ('M', 'D4', 41), 'S08': ('F', 'D3', 49)}
for row in rows: row['sex'], row['donor'], row['age'] = extra[row['sample_id']]
with sheet.open('w') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t'); w.writeheader(); w.writerows(rows)
backend, formula, contrasts = {
    'interaction': ('deseq2', '~ sex * treatment',
                    [{'id': 'trt', 'type': 'pairwise', 'factor': 'treatment', 'by': 'sex',
                      'numerator': 'treated', 'denominator': 'control'},
                     {'id': 'inter', 'type': 'linear', 'weights': {'sexM:treatmenttreated': 1}}]),
    'paired': ('edger_ql', '~ donor + treatment',
               [{'id': 'trt', 'type': 'pairwise', 'factor': 'treatment', 'numerator': 'treated', 'denominator': 'control'}]),
    'continuous': ('limma_voom', '~ treatment + age', [{'id': 'age', 'type': 'linear', 'weights': {'age': 1}}]),
}[case]
c['analysis'] = {'objectives': ['differential_expression'], 'backend': backend, 'umi': False}
c['model']['fixed_effects'] = formula
c['contrasts'] = contrasts
cfg_path.write_text(yaml.safe_dump(c))
PY
    if ! run_workflow "$proj" "$REPO" "$TASK_TMP/$case.log" results/de_done.flag; then
        tail -70 "$TASK_TMP/$case.log"; exit 1
    fi
done

Rscript - "$TASK_TMP" <<'R'
suppressPackageStartupMessages({library(edgeR); library(limma); library(DESeq2); library(jsonlite)})
root <- commandArgs(TRUE)[1]
load_case <- function(case, formula) {
  p <- file.path(root, case)
  counts <- as.matrix(read.delim(file.path(p, "results/counts.tsv"), row.names = 1, check.names = FALSE))
  sheet <- read.delim(file.path(p, "config/samples.tsv"))
  sheet <- sheet[match(colnames(counts), sheet$sample_id), ]
  mm <- model.matrix(formula, sheet)
  d <- DGEList(counts); d <- calcNormFactors(d)
  list(p = p, sheet = sheet, mm = mm, d = d, keep = filterByExpr(d, mm))
}
table_for <- function(p, name) read.delim(file.path(p, "results/DE_Results", paste0(name, "_DE_analysis.tsv")), check.names = FALSE)
agree <- function(actual, logfc, pvalue, padj, ids, label) {
  i <- match(actual$gene_id, ids)
  stopifnot(!anyNA(i), max(abs(actual$logFC - logfc[i])) < 1e-8,
            max(abs(actual$P.Value - pvalue[i])) < 1e-8, max(abs(actual$adj.P.Val - padj[i])) < 1e-8)
  cat(label, "agrees with the direct fit\n")
}

# Interaction: treatment effect within each sex, and their difference.
x <- load_case("interaction", ~ sex * treatment)
stopifnot(identical(colnames(x$mm), c("(Intercept)", "sexM", "treatmenttreated", "sexM:treatmenttreated")))
coldata <- x$sheet; rownames(coldata) <- colnames(x$d)
dds <- DESeqDataSetFromMatrix(round(x$d$counts[x$keep, ]), coldata, design = x$mm)
dds <- DESeq(dds, betaPrior = FALSE, minReplicatesForReplace = Inf, parallel = FALSE, quiet = TRUE)
dirs <- read.delim(file.path(x$p, "results/DE_Results/contrast_directions.tsv"))
stopifnot(nrow(dirs) == 3, all(dirs$numerator_level[dirs$spec_id == "trt"] == "treated"))
vectors <- list(F = c(0, 0, 1, 0), M = c(0, 0, 1, 1))
for (i in which(dirs$spec_id == "trt")) {
  r <- results(dds, contrast = vectors[[dirs$by_level[i]]], independentFiltering = FALSE)
  agree(table_for(x$p, dirs$contrast_name[i]), r$log2FoldChange, r$pvalue, r$padj, rownames(r),
        paste("DESeq2 treatment within sex", dirs$by_level[i]))
}
r <- results(dds, contrast = c(0, 0, 0, 1), independentFiltering = FALSE)
agree(table_for(x$p, "inter"), r$log2FoldChange, r$pvalue, r$padj, rownames(r), "DESeq2 interaction coefficient")

# Paired donors as a fixed block.
x <- load_case("paired", ~ donor + treatment)
d <- x$d[x$keep, , keep.lib.sizes = FALSE]; d <- calcNormFactors(d)
d <- estimateDisp(d, x$mm, robust = TRUE)
tt <- topTags(glmQLFTest(glmQLFit(d, x$mm, robust = TRUE), coef = "treatmenttreated"), n = Inf)$table
name <- read.delim(file.path(x$p, "results/DE_Results/contrast_directions.tsv"))$contrast_name
agree(table_for(x$p, name), tt$logFC, tt$PValue, tt$FDR, rownames(tt), "edgeR QL donor-blocked pairing")
stopifnot(fromJSON(file.path(x$p, "metrics/de.json"))$backend == "edger_ql")

# Continuous covariate slope.
x <- load_case("continuous", ~ treatment + age)
fit <- eBayes(lmFit(voom(x$d[x$keep, ], x$mm), x$mm))
tt <- topTable(fit, coef = "age", n = Inf, sort.by = "none")
agree(table_for(x$p, "age"), tt$logFC, tt$P.Value, tt$adj.P.Val, rownames(tt), "limma-voom continuous age slope")
R

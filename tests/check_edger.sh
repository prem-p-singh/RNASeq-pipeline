#!/usr/bin/env bash
set -euo pipefail
REPO=$(cd "$(dirname "$0")/.." && pwd)
. "$REPO/tests/_gate.sh"
gate_runtime r
TASK_TMP=$(mktemp -d)
trap 'rm -rf "$TASK_TMP"' EXIT
proj="$TASK_TMP/project"
python3 "$REPO/tests/fixtures/tiny_project.py" build "$proj" "$REPO" rnaseq
python3 - "$proj" <<'PY'
import sys,yaml
from pathlib import Path
p=Path(sys.argv[1])/'config/config.yaml'; c=yaml.safe_load(p.read_text())
c['analysis']={'objectives':['differential_expression'], 'backend':'edger_ql', 'umi':False}
c['contrasts']=[{'id':'effect','type':'linear','weights':{'treatmenttreated':1}}]
p.write_text(yaml.safe_dump(c))
PY
if ! run_workflow "$proj" "$REPO" "$TASK_TMP/run.log" results/de_done.flag; then
    tail -70 "$TASK_TMP/run.log"; exit 1
fi
Rscript - "$proj" <<'R'
suppressPackageStartupMessages({library(edgeR);library(jsonlite)})
p <- commandArgs(TRUE)[1]
counts <- as.matrix(read.delim(file.path(p,"results/counts.tsv"),row.names=1,check.names=FALSE))
sheet <- read.delim(file.path(p,"config/samples.tsv")); sheet <- sheet[match(colnames(counts),sheet$sample_id),]
mm <- model.matrix(~treatment, sheet)
d <- DGEList(counts); d <- calcNormFactors(d)
d <- d[filterByExpr(d,mm),,keep.lib.sizes=FALSE]; d <- calcNormFactors(d)
d <- estimateDisp(d,mm,robust=TRUE)
expected <- topTags(glmQLFTest(glmQLFit(d,mm,robust=TRUE),coef=2),n=Inf)$table
actual <- read.delim(file.path(p,"results/DE_Results/effect_DE_analysis.tsv"),check.names=FALSE)
expected <- expected[match(actual$gene_id,rownames(expected)),]
stopifnot(max(abs(actual$logFC-expected$logFC)) < 1e-8,
          max(abs(actual$P.Value-expected$PValue)) < 1e-8,
          max(abs(actual$adj.P.Val-expected$FDR)) < 1e-8)
manifest <- read.delim(file.path(p,"results/de_manifest.tsv"))
stopifnot(manifest$statistic == "signed_sqrt_F",
          fromJSON(file.path(p,"metrics/de.json"))$backend == "edger_ql")
cat("edgeR QL full adapter agrees with independent direct GLM fit\n")
R

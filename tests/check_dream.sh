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
c['analysis']={'objectives':['differential_expression'], 'backend':'auto', 'umi':False}
c['model']['random_effects']='(1|donor)'
import csv
sheet=p.parent/'samples.tsv'
with sheet.open() as f: rows=list(csv.DictReader(f,delimiter='\t'))
for i,row in enumerate(rows): row['donor']='D'+str(i%4+1)
with sheet.open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t'); w.writeheader(); w.writerows(rows)
c['contrasts']=[{'id':'effect','type':'linear','weights':{'treatmenttreated':1}}]
p.write_text(yaml.safe_dump(c))
PY
if ! run_workflow "$proj" "$REPO" "$TASK_TMP/run.log" results/de_done.flag; then
    tail -70 "$TASK_TMP/run.log"; exit 1
fi
Rscript - "$proj" <<'R'
suppressPackageStartupMessages({library(edgeR);library(limma);library(variancePartition);library(BiocParallel);library(jsonlite)})
p <- commandArgs(TRUE)[1]
counts <- as.matrix(read.delim(file.path(p,"results/counts.tsv"),row.names=1,check.names=FALSE))
sheet <- read.delim(file.path(p,"config/samples.tsv")); sheet <- sheet[match(colnames(counts),sheet$sample_id),]
mm <- model.matrix(~treatment, sheet)
d <- DGEList(counts); d <- calcNormFactors(d)
d <- d[filterByExpr(d,mm),]
sheet$treatment <- factor(sheet$treatment); sheet$donor <- factor(sheet$donor)
form <- ~ treatment + (1|donor)
v <- voomWithDreamWeights(d,form,sheet,BPPARAM=SerialParam())
L <- matrix(c(0,1),ncol=1,dimnames=list(colnames(mm),'effect'))
fit <- eBayes(dream(v,form,sheet,L,BPPARAM=SerialParam()))
expected <- topTable(fit,coef='effect',n=Inf)
actual <- read.delim(file.path(p,"results/DE_Results/effect_DE_analysis.tsv"),check.names=FALSE)
expected <- expected[match(actual$gene_id,rownames(expected)),]
stopifnot(max(abs(actual$logFC-expected$logFC)) < 1e-6,
          max(abs(actual$P.Value-expected$P.Value)) < 1e-6,
          max(abs(actual$adj.P.Val-expected$adj.P.Val)) < 1e-6)
manifest <- read.delim(file.path(p,"results/de_manifest.tsv"))
stopifnot(manifest$statistic == "t",
          fromJSON(file.path(p,"metrics/de.json"))$backend == "dream")
cat("dream paired-donor adapter agrees with independent direct mixed-model fit\n")
R

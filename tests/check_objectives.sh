#!/usr/bin/env bash
set -euo pipefail
REPO=$(cd "$(dirname "$0")/.." && pwd)
. "$REPO/tests/_gate.sh"
gate_runtime r
command -v fastp >/dev/null || exit 77
TASK_TMP=$(mktemp -d)
trap 'rm -rf "$TASK_TMP"' EXIT
proj="$TASK_TMP/project"
python3 "$REPO/tests/fixtures/fastq_project.py" build "$proj" "$REPO" single
python3 - "$proj" <<'PY'
import sys,yaml
from pathlib import Path
p=Path(sys.argv[1]); f=p/'config/config.yaml'; c=yaml.safe_load(f.read_text())
(p/'original-config.yaml').write_text(f.read_text())
c['analysis']={'objectives':['qc'], 'backend':'auto', 'umi':False}
c['model']={'fixed_effects':'~ absent', 'primary_factor':'absent', 'random_effects':None}
c['reference']['transcriptome_fasta_url']='https://invalid.test/does-not-exist.fa.gz'
c['reference']['gtf_url']='https://invalid.test/does-not-exist.gtf.gz'
c['preprocessing']={'poly_g':'off'}
f.write_text(yaml.safe_dump(c))
PY
python3 "$REPO/scripts/preflight.py" -d "$proj"
if ! run_workflow "$proj" "$REPO" "$TASK_TMP/qc.log"; then tail -60 "$TASK_TMP/qc.log"; exit 1; fi
test -s "$proj/results/preprocessing_report/multiqc_report.html"
test ! -e "$proj/results/counts.tsv"
test ! -e "$proj/reference/salmon_idx"
python3 - "$proj" <<'PY'
import json,sys,yaml
from pathlib import Path
p=Path(sys.argv[1]); f=p/'config/config.yaml'; c=yaml.safe_load(f.read_text())
for ledger in (p/'results/quant').glob('*/preprocessing.json'):
 d=json.loads(ledger.read_text()); assert d['policy']['poly_g']=='off'
 assert d['reads_before'] >= d['reads_after'] > 0
c['analysis']['objectives']=['gene_expression']
c['reference']=yaml.safe_load((p/'original-config.yaml').read_text())['reference']
f.write_text(yaml.safe_dump(c))
PY
python3 "$REPO/scripts/preflight.py" -d "$proj"
if ! run_workflow "$proj" "$REPO" "$TASK_TMP/counts.log"; then tail -60 "$TASK_TMP/counts.log"; exit 1; fi
test -s "$proj/results/counts.tsv"
test ! -e "$proj/results/de_done.flag"
test ! -e "$proj/results/wgcna_done.flag"
echo 'QC-only and count-only workflows passed without a usable DE model; QC required no reference'

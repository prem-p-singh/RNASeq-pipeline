#!/usr/bin/env bash
# Networked release smoke test. Set RNASEQ_PUBLIC_PROJECT to preserve all outputs.
set -euo pipefail
REPO=$(cd "$(dirname "$0")/.." && pwd)
. "$REPO/tests/_gate.sh"
gate_runtime r
proj=${RNASEQ_PUBLIC_PROJECT:-$(mktemp -d)/public}
python3 "$REPO/tests/fixtures/public_project.py" "$proj" "$REPO"
python3 "$REPO/scripts/preflight.py" -d "$proj"
if ! run_workflow "$proj" "$REPO" "$proj/run.log"; then
    tail -80 "$proj/run.log"; exit 1
fi
python3 - "$proj" <<'PY'
import csv, json, sys
from pathlib import Path
p=Path(sys.argv[1])
agg=json.loads((p/'metrics/aggregate.json').read_text())
assert agg['n_samples']==6 and agg['n_genes']>50, agg
manifest=list(csv.DictReader((p/'results/de_manifest.tsv').open(), delimiter='\t'))
assert manifest
for row in manifest:
    assert (p/'results'/row['analysis_table']).stat().st_size > 0
for i in range(6):
    m=json.loads((p/f'results/quant/SRR{6357070+i}/metrics.json').read_text())
    assert m['num_reads_processed']>10000 and m['mapping_rate']>0.2, m
assert (p/'results/qc_report/qc_charts.html').stat().st_size > 0
print('GSE110004 smoke: six biological samples, gene counts, DE tables and QC report passed')
PY
echo "Public test project: $proj"

#!/usr/bin/env bash
set -euo pipefail
REPO=$(cd "$(dirname "$0")/.." && pwd)
. "$REPO/tests/_gate.sh"
gate_runtime r
for tool in salmon fastp multiqc; do
    command -v "$tool" >/dev/null || { echo "needs $tool"; exit 77; }
done
python3 -c 'import matplotlib' 2>/dev/null || { echo "needs matplotlib"; exit 77; }
TASK_TMP=$(mktemp -d)
cleanup() {
    rc=$?
    if [ "$rc" -ne 0 ]; then
        for f in "$TASK_TMP"/*.log; do echo "Log: $f"; tail -60 "$f"; done
    fi
    rm -rf "$TASK_TMP"
}
trap cleanup EXIT
for layout in paired single; do
    proj="$TASK_TMP/$layout project"
    python3 "$REPO/tests/fixtures/fastq_project.py" build "$proj" "$REPO" "$layout"
    if [ "$layout" = paired ]; then
        python3 "$REPO/tests/fixtures/fastq_project.py" canonicalize "$proj" "$REPO"
    fi
    if ! run_workflow "$proj" "$REPO" "$TASK_TMP/$layout.log"; then
        tail -70 "$TASK_TMP/$layout.log"; exit 1
    fi
    python3 "$REPO/tests/fixtures/fastq_project.py" verify "$proj"
    before=$(python3 -c 'import os,sys; print(os.stat(sys.argv[1]).st_mtime_ns)' "$proj/results/counts.tsv")
    run_workflow "$proj" "$REPO" "$TASK_TMP/resume.log"
    after=$(python3 -c 'import os,sys; print(os.stat(sys.argv[1]).st_mtime_ns)' "$proj/results/counts.tsv")
    [ "$before" = "$after" ] || { echo 'Unchanged resume rebuilt counts'; exit 1; }
    rm "$proj/results/qc_report/qc_charts.html"
    run_workflow "$proj" "$REPO" "$TASK_TMP/repair.log"
    test -s "$proj/results/qc_report/qc_charts.html"
    if [ "$layout" = paired ]; then
        rm "$proj/results/quant/S1/read_preparation.json"
        run_workflow "$proj" "$REPO" "$TASK_TMP/provenance-repair.log"
        test -s "$proj/results/quant/S1/read_preparation.json"
        python3 - "$proj" <<'CHANGE'
import gzip, sys, yaml
from pathlib import Path
p=Path(sys.argv[1])
c=p/'config/config.yaml'
cfg=yaml.safe_load(c.read_text()); cfg['contrasts'][0]['reverse']=False
c.write_text(yaml.safe_dump(cfg))
f=p/'inputs/run1_S1_R1.fq.gz'
with gzip.open(f, 'rt') as h: text=h.read()
with gzip.open(f, 'wt') as h: h.write(text.replace('\n+\nII', '\n+\nHI', 1))
CHANGE
        quant_before=$(python3 -c 'import os,sys; print(os.stat(sys.argv[1]).st_mtime_ns)' "$proj/results/quant/S1/quant.sf")
        run_workflow "$proj" "$REPO" "$TASK_TMP/changed.log"
        python3 - "$proj" "$quant_before" <<'VERIFY'
import csv, sys
from pathlib import Path
p=Path(sys.argv[1])
assert (p/'results/quant/S1/quant.sf').stat().st_mtime_ns != int(sys.argv[2]), 'Changed local FASTQ reused stale quantification'
import json
ledger=json.loads((p/'results/quant/S1/read_preparation.json').read_text())
assert len(ledger['units']) == 2 and ledger['fragments'] > 0
assert not (p/'results/quant/S1/prepared_R1.fastq.gz').exists(), 'Owned merged reads were not cleaned'
row=next(csv.DictReader((p/'results/de_manifest.tsv').open(), delimiter='\t'))
rows={x['gene_id']: x for x in csv.DictReader((p/'results'/row['analysis_table']).open(), delimiter='\t')}
assert float(rows['G000']['logFC']) < -1, 'Reversed contrast reused stale direction'
VERIFY
    fi
done
echo 'check_end_to_end.sh: real FASTQ-to-report SE/PE and restart passed'

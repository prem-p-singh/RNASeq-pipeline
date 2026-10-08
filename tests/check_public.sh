#!/usr/bin/env bash
# Networked release smoke test. Set RNASEQ_PUBLIC_PROJECT to preserve all outputs.
set -euo pipefail
REPO=$(cd "$(dirname "$0")/.." && pwd)
. "$REPO/tests/_gate.sh"
gate_runtime r
proj=${RNASEQ_PUBLIC_PROJECT:-$(mktemp -d)/public}
# --verify checks a project that another launcher already ran.
if [ "${1:-}" != --verify ]; then
    python3 "$REPO/tests/fixtures/public_project.py" "$proj" "$REPO" "$@"
    python3 "$REPO/scripts/preflight.py" -d "$proj"
    if ! run_workflow "$proj" "$REPO" "$proj/run.log"; then
        tail -80 "$proj/run.log"; exit 1
    fi
fi
python3 - "$proj" <<'PY'
import csv, json, sys, yaml
from pathlib import Path
p=Path(sys.argv[1])
cfg=yaml.safe_load((p/'config/config.yaml').read_text())
star=cfg.get('analysis',{}).get('quantifier')=='star'
agg=json.loads((p/'metrics/aggregate.json').read_text())
assert agg['n_samples']==6 and agg['n_genes']>50, agg
manifest=list(csv.DictReader((p/'results/de_manifest.tsv').open(), delimiter='\t'))
assert manifest
for row in manifest:
    assert (p/'results'/row['analysis_table']).stat().st_size > 0
for i in range(6):
    root=p/f'results/quant/SRR{6357070+i}'
    screen=json.loads((root/'screening/screening.json').read_text())
    if cfg.get('screening',{}).get('enabled'):
        assert screen['status']=='completed' and screen['reads_removed']==0
        assert screen['fragments_sampled']==10000 and len(screen['results'])==2
        assert all(r['exclusive_to_reference']>1000 for r in screen['results']), screen
    else:
        assert screen['status']=='not_performed'
    if star: root=root/'star'
    m=json.loads((root/'metrics.json').read_text())
    assert m['num_reads_processed']>10000 and m['mapping_rate']>0.2, m
    if star:
        assert m['num_assigned']>1000, m
        assert (root/'Aligned.sortedByCoord.out.bam.bai').is_file()
assert (p/'results/qc_report/qc_charts.html').stat().st_size > 0
print(f'GSE110004 {"STAR" if star else "Salmon"} smoke: six biological samples, gene counts, DE tables and QC report passed')
PY
echo "Public test project: $proj"

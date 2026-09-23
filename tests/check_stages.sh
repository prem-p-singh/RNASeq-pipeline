#!/usr/bin/env bash
# Runs Stages 2 to 5 on a tiny synthetic cohort and checks the answers.
#
# WORKING_PLAN 3.3: "a tiny synthetic fixture through stages 0 to 5, which is
# the first execution of stages 2 to 5 and the evidence R01 requires". Stages 0
# and 1 need salmon and fastp; tests/fixtures/tiny_project.py writes their
# products directly, and everything downstream runs for real through Snakemake.
#
# What it would catch: tx2gene aggregation that takes one transcript instead of
# the sum, a lost gene when a quant.sf name carries a "|", a QC flag that does
# not actually control inclusion, a contrast whose sign is the reverse of its
# name, a count matrix that silently keeps an excluded sample, and a policy skip
# reported as a failure.
#
# Needs snakemake and the R stack (see the gate below). Exits 77 (skipped) without
# them, which is the normal state until the environment of WORKING_PLAN 3.2
# item 3 exists.
#
# Run:  bash tests/check_stages.sh
set -uo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
FIXTURE="$REPO/tests/fixtures/tiny_project.py"
. "$REPO/tests/_gate.sh"
gate_runtime r

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
PROJ="$TMP/proj"

python3 "$FIXTURE" build "$PROJ" "$REPO" || { echo "FAIL: fixture build" >&2; exit 1; }
python3 - "$PROJ" <<'BATCH'
import csv, sys, yaml
from pathlib import Path
p = Path(sys.argv[1])
sheet = p / 'config/samples.tsv'
with sheet.open() as f:
    rows = list(csv.DictReader(f, delimiter='\t'))
for i, row in enumerate(rows):
    row['batch'] = 'a' if i % 2 else 'b'
with sheet.open('w') as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
    writer.writeheader(); writer.writerows(rows)
c = p / 'config/config.yaml'
cfg = yaml.safe_load(c.read_text())
cfg.setdefault('thresholds', {})['batch_correction'] = {'pc_var_fraction_trigger': 0}
c.write_text(yaml.safe_dump(cfg))
BATCH

run_workflow "$PROJ" "$REPO" "$TMP/run.log" \
    results/counts.tsv results/de_done.flag \
    results/enrichment_done.flag results/wgcna_done.flag
rc=$?
if [ $rc -ne 0 ]; then
    echo "FAIL: stages 2-5 did not complete (exit $rc)" >&2
    tail -40 "$TMP/run.log" >&2
    exit 1
fi

python3 "$FIXTURE" verify "$PROJ" || exit 1
python3 - "$PROJ" <<'CHECK_BATCH'
import json, sys, yaml
from pathlib import Path
p = Path(sys.argv[1])
metrics = json.loads((p / 'metrics/de.json').read_text())
cfg = yaml.safe_load((p / 'config/config.yaml').read_text())
assert metrics['model_fixed'] == cfg['model']['fixed_effects'], 'PCA changed the declared model'
assert metrics['auto_batch_added'] is False
assert metrics['batch_diagnostic']['review_suggested'] is True
assert 'BATCH_REVIEW' in (p / 'gates/decisions.log').read_text()
CHECK_BATCH
echo "check_stages.sh: stages 2-5 ran and every expected answer held"

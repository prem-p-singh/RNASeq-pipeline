#!/usr/bin/env bash
set -euo pipefail
REPO=$(cd "$(dirname "$0")/.." && pwd)
. "$REPO/tests/_gate.sh"
gate_runtime
for tool in fastp fastq_screen bowtie2 bowtie2-build; do command -v "$tool" >/dev/null || exit 77; done
TASK_TMP=$(mktemp -d)
trap 'rm -rf "$TASK_TMP"' EXIT
proj="$TASK_TMP/project"
python3 "$REPO/tests/fixtures/fastq_project.py" build "$proj" "$REPO" paired
bowtie2-build "$proj/inputs/transcripts.fa" "$proj/inputs/screen" > "$TASK_TMP/build.log" 2>&1
python3 - "$proj" <<'PY'
import sys,yaml
from pathlib import Path
p=Path(sys.argv[1]);f=p/'config/config.yaml';c=yaml.safe_load(f.read_text())
c['analysis']={'objectives':['qc']}
c['screening']={'enabled':True,'fragments':100,'seed':42,
 'references':[{'name':'Expected','role':'expected','index':'inputs/screen'}]}
f.write_text(yaml.safe_dump(c))
PY
# Relative panel paths resolve against the project even when preflight is run elsewhere.
python3 "$REPO/scripts/preflight.py" -d "$proj"
if ! run_workflow "$proj" "$REPO" "$TASK_TMP/run.log"; then
 tail -70 "$TASK_TMP/run.log"; cat "$proj"/logs/screening/*.log; exit 1
fi
python3 - "$proj" <<'PY'
import json,sys
from pathlib import Path
p=Path(sys.argv[1]);ledgers=list((p/'results/quant').glob('*/screening/screening.json'))
assert len(ledgers)==6,len(ledgers)
for ledger in ledgers:
 d=json.loads(ledger.read_text()); assert d['status']=='completed'
 assert d['fragments_sampled']==100 and d['reads_removed']==0
 assert len(d['results'])==2 and all(r['reads']==100 for r in d['results'])
 assert all(r['exclusive_to_reference']>80 for r in d['results']),d
assert not list((p/'results/quant').rglob('*.trim.fastq.gz'))
assert not (p/'results/counts.tsv').exists()
(p/'times.json').write_text(json.dumps({str(x):x.stat().st_mtime_ns for x in ledgers}))
PY
run_workflow "$proj" "$REPO" "$TASK_TMP/resume.log"
python3 - "$proj" <<'PY'
import json,sys
from pathlib import Path
p=Path(sys.argv[1]);old=json.loads((p/'times.json').read_text())
assert all(Path(f).stat().st_mtime_ns==n for f,n in old.items()),'Screen resume rebuilt outputs'
Path(next(iter(old))).with_name('screening.tsv').unlink()
PY
run_workflow "$proj" "$REPO" "$TASK_TMP/recovery.log"
python3 - "$proj" <<'PY'
import sys
from pathlib import Path
p=Path(sys.argv[1]);assert len(list((p/'results/quant').glob('*/screening/screening.tsv')))==6
assert not list((p/'results/quant').rglob('*.trim.fastq.gz'))
PY
echo 'Screen QC-only DAG, temporary read lifetime, resume and missing-output recovery passed'

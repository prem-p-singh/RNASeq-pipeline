#!/usr/bin/env bash
set -euo pipefail
REPO=$(cd "$(dirname "$0")/.." && pwd)
. "$REPO/tests/_gate.sh"
gate_runtime r
for tool in STAR featureCounts samtools; do command -v "$tool" >/dev/null || exit 77; done
TASK_TMP=$(mktemp -d)
trap 'rm -rf "$TASK_TMP"' EXIT
proj="$TASK_TMP/project"
python3 "$REPO/tests/fixtures/fastq_project.py" build "$proj" "$REPO" paired
python3 - "$proj" <<'PY'
from pathlib import Path
import sys,yaml
p=Path(sys.argv[1]);cfg=p/'config/config.yaml';c=yaml.safe_load(cfg.read_text())
records=(p/'inputs/transcripts.fa').read_text().split('>')[1:]
genome=p/'inputs/genome.fa'
genome.write_text('>chr1\n'+''.join(''.join(r.splitlines()[1:])+'N'*1000 for r in records)+'\n')
c['reference']['genome_fasta_url']=genome.as_uri()
c['samples']['expected_libtype']='IU'
c['analysis']={'objectives':['gene_expression','differential_expression'],'quantifier':'star','backend':'edger_ql'}
cfg.write_text(yaml.safe_dump(c))
PY
python3 "$REPO/scripts/preflight.py" -d "$proj"
if ! run_workflow "$proj" "$REPO" "$TASK_TMP/run.log"; then
 tail -80 "$TASK_TMP/run.log"
 for f in "$proj"/logs/star/*.log "$proj"/logs/star_index.log; do test ! -f "$f" || tail -40 "$f"; done
 exit 1
fi
python3 - "$proj" <<'PY'
from pathlib import Path
import csv,json,sys
p=Path(sys.argv[1]); expected=json.loads((p/'expected.json').read_text())
with (p/'results/counts.tsv').open() as f: rows=list(csv.DictReader(f,delimiter='\t'))
for row in rows:
 for sample in expected:
  assert int(row[sample])==expected[sample][row['gene_id']],(sample,row)
prov=json.loads((p/'results/counts_provenance.json').read_text())
assert prov['quantifier']=='STAR_featureCounts' and not prov['length_correction_applied']
assert prov['count_unit']=='fragments'
assert not (p/'reference/salmon_idx').exists()
assert (p/'results/qc_report/qc_charts.html').is_file()
assert json.loads((p/'metrics/de.json').read_text())['backend']=='edger_ql'
for sample in expected:
 assert (p/f'results/quant/{sample}/star/Aligned.sortedByCoord.out.bam.bai').is_file()
print('STAR raw-read workflow: exact independent fragment counts, no length scaling, DE and QC passed')
PY
before=$(python3 -c 'import os,sys;print(os.stat(sys.argv[1]).st_mtime_ns)' "$proj/results/counts.tsv")
run_workflow "$proj" "$REPO" "$TASK_TMP/resume.log"
after=$(python3 -c 'import os,sys;print(os.stat(sys.argv[1]).st_mtime_ns)' "$proj/results/counts.tsv")
[ "$before" = "$after" ] || { echo 'STAR resume rebuilt counts'; exit 1; }
missing=$(find "$proj/results/quant" -name '*.bam.bai' -print -quit)
rm "$missing"
run_workflow "$proj" "$REPO" "$TASK_TMP/recovery.log"
test -s "$missing"
echo 'STAR missing BAM-index recovery passed'

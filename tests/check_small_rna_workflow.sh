#!/usr/bin/env bash
# Synthetic small-RNA animal route check, run through results/counts.tsv.
# A local mature-miRNA FASTA (miRBase layout, RNA alphabet, prefix "fix") holds:
#   fix-miR-1                 unique sequence
#   fix-miR-3                 unique sequence, wrapped over two lines
#   fix-miR-2a, fix-miR-2b    one identical sequence under two IDs (multi-mapping)
#   oth-miR-1                 another species, removed by the prefix filter
# Reads are TruSeq Small RNA geometry (insert + 3' adapter + filler). Per sample:
#   S1: miR-1 x5, miR-3 x2, miR-2 x3, oth x1, random x1
#   S2: miR-1 x3, miR-3 x4, miR-2 x1, oth x2
#   S3: miR-1 x6,           miR-2 x2,         random x2
# plus one 10 nt insert per sample, below the kit's 18 nt minimum and lost in fastp.
# oth and random reads are unaligned. Expected counts are in expected.json.
set -euo pipefail
REPO=$(cd "$(dirname "$0")/.." && pwd)
. "$REPO/tests/_gate.sh"
gate_runtime r
for tool in fastp bowtie bowtie-build; do
    command -v "$tool" >/dev/null || { echo "needs $tool on PATH"; exit 77; }
done
# The Snakefile hashes the srna module lock into its runtime identity.
test -f "$REPO/environments/srna-linux-64.explicit.txt" \
    || { echo "needs environments/srna-linux-64.explicit.txt (generated on the Linux host)"; exit 77; }
TASK_TMP=$(mktemp -d)
trap 'rm -rf "$TASK_TMP"' EXIT
proj="$TASK_TMP/srna"

python3 "$REPO/tests/fixtures/small_rna_workflow_project.py" "$proj" "$REPO"
python3 "$REPO/scripts/preflight.py" -d "$proj"
if ! run_workflow "$proj" "$REPO" "$proj/run.log" results/counts.tsv; then
    tail -80 "$proj/run.log"
    for f in "$proj"/logs/small_rna/*.log "$proj"/logs/preprocess/*.log "$proj"/logs/mature_mirna_reference.log; do
        test ! -f "$f" || tail -40 "$f"
    done
    exit 1
fi

python3 - "$proj" <<'PY'
import csv, json, sys
from pathlib import Path
proj = Path(sys.argv[1])
expected = json.loads((proj / 'expected.json').read_text())
for s, e in expected.items():
    q = proj / 'results/quant' / s
    with open(q / 'small_rna/counts.tsv') as f:
        counts = {r['gene_id']: int(r['count']) for r in csv.DictReader(f, delimiter='\t')}
    assert counts == e['counts'], (s, counts)
    m = json.loads((q / 'small_rna/metrics.json').read_text())
    assert m['alignment_summary'] == e['alignment'], (s, m['alignment_summary'])
    ledger = json.loads((q / 'preprocessing.json').read_text())
    assert ledger['size_selection']['too_short_reads'] == 1, (s, ledger['size_selection'])
    lost = sum(n for k, n in ledger['filtering_result'].items() if k != 'passed_filter_reads')
    assert (ledger['reads_before'], ledger['reads_after'], lost) == (e['reads_before'], e['reads_before'] - 1, 1), (s, ledger)
    assert m['num_reads_processed'] == ledger['reads_after'], (s, m)
with open(proj / 'results/counts.tsv') as f:
    rows = {r['gene_id']: r for r in csv.DictReader(f, delimiter='\t')}
assert set(rows) == {'fix-miR-1', 'fix-miR-3'}, sorted(rows)
for s, e in expected.items():
    for gene in rows:
        assert int(rows[gene][s]) == e['counts'].get(gene, 0), (s, gene, rows[gene])
prov = json.loads((proj / 'results/counts_provenance.json').read_text())
assert prov['quantifier'] == 'bowtie_mature_mirna', prov
print('small-RNA animal route: unique hits counted, multi-mapped and unaligned reads reported, short read lost in preprocessing')
PY

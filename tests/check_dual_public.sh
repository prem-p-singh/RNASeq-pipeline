#!/usr/bin/env bash
# Public-data check for the dual-organism route: grapevine with GRBV (--virus) or
# Arabidopsis with Botrytis cinerea (--cellular).
#
# Runs on data already prepared by tests/fixtures/prepare_dual_public.py (runs
# downsampled where the study sets a cap, one infected group and one healthy or mock group). The three phases are run
# the way the launcher runs them: the parent workflow (assign and split),
# `dual_organism.py children`, then each child project. Between the first two the
# assignment is run a second time and its counts compared. The workflow fetches the
# pinned references itself. The five acceptance criteria are the ones written in
# tests/fixtures/dual_public.json before any run.
# Skips (exit 77) without prepared data; not in run_all.sh.
#
# Run:
#   python3 tests/fixtures/prepare_dual_public.py virus DATA_DIR    # once
#   bash tests/check_dual_public.sh --virus DATA_DIR
# (--cellular with "cellular"). A third argument, star (default) or salmon, sets
# analysis.quantifier for both organisms. RNASEQ_DUAL_PUBLIC_PROJECT names the
# project folder to create; SMK_EXTRA passes further Snakemake options.
set -euo pipefail
REPO=$(cd "$(dirname "$0")/.." && pwd)
usage="Usage: check_dual_public.sh --virus|--cellular PREPARED_DATA_DIR [star|salmon]"
case "${1:-}" in
    --virus) study=virus ;;
    --cellular) study=cellular ;;
    *) echo "$usage" >&2; exit 2 ;;
esac
data=${2:-}
quantifier=${3:-star}
[ -n "$data" ] || { echo "$usage" >&2; exit 2; }
case "$quantifier" in star|salmon) ;; *) echo "$usage" >&2; exit 2 ;; esac
if [ ! -f "$data/samples.tsv" ]; then
    echo "needs prepared public data at $data (run tests/fixtures/prepare_dual_public.py $study $data first)"
    exit 77
fi

. "$REPO/tests/_gate.sh"
gate_runtime r
# The assignment always uses STAR; the Salmon route also needs salmon.
for tool in STAR featureCounts samtools fastp multiqc $([ "$quantifier" = salmon ] && echo salmon); do
    command -v "$tool" >/dev/null || { echo "needs $tool on PATH"; exit 77; }
done

proj=${RNASEQ_DUAL_PUBLIC_PROJECT:-$(mktemp -d)/dual_public}
python3 - "$study" "$data" "$proj" "$REPO" "$quantifier" <<'PY'
import csv, json, sys
from pathlib import Path
import yaml
study, data, proj, repo, quantifier = sys.argv[1], Path(sys.argv[2]).resolve(), Path(sys.argv[3]), Path(sys.argv[4]), sys.argv[5]
ledger = json.loads((data / 'download_manifest.json').read_text())
if ledger['key'] != study or ledger['status'] != 'downloaded_not_qualified':
    raise SystemExit(f'{data} holds {ledger["key"]!r} with status {ledger["status"]!r}; expected {study!r} downloaded_not_qualified')
profile = json.loads((data / 'profile.json').read_text())
with open(data / 'samples.tsv') as f:
    rows = list(csv.DictReader(f, delimiter='\t'))
groups = sorted({r['condition'] for r in rows})
assert len(rows) == len(profile['runs']) and len(groups) == 2 and 'infected' in groups, (len(rows), groups)
(proj / 'config').mkdir(parents=True)
with open(proj / 'config/samples.tsv', 'w') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['sample_id', 'fastq_url', 'fastq_url_r2', 'treatment'])
    w.writerows([r['sample_id'], r['fastq_path'], r['fastq_path_r2'], r['condition']] for r in rows)


def organism(role):
    o = profile[role]
    urls = {'genome_fasta_url': o['genome_fasta']['url']}
    if 'gtf' in o:
        urls.update(gtf_url=o['gtf']['url'], transcriptome_fasta_url=o['transcriptome_fasta']['url'])
    else:  # annotation written by prepare_dual_public.py from a GTF without exon rows
        urls.update(gtf_url=(data / 'pathogen.exons.gtf').as_uri(), transcriptome_fasta_url=(data / 'pathogen.transcripts.fa').as_uri())
    return dict(label=o['label'], **urls)


host = profile['host']
cfg = {
    'project': {'name': f'dual_public_{profile["study"]}', 'output_dir': 'results/'},
    'organism': {'scientific_name': host['scientific_name'], 'tax_id': host['tax_id'], 'kegg_code': None, 'orgdb_package': None},
    'reference': {'accession': f'dual_public_{study}', 'cache_dir': None, 'annotation_tsv': {'path': None}},
    'samples': {'seq_type': profile['seq_type'], 'sheet': 'config/samples.tsv', 'expected_libtype': profile['expected_libtype']},
    'analysis': {'objectives': ['gene_expression', 'differential_expression'], 'quantifier': quantifier},
    'model': {'fixed_effects': '~ treatment', 'random_effects': None, 'primary_factor': 'treatment'},
    'contrasts': [{'id': 'treatment', 'type': 'pairwise', 'factor': 'treatment', 'reverse': True}],
    'downstream': {'run_go': False, 'run_kegg': False, 'run_wgcna': False},
    'orgdb': {'strategy': 'skip'},
    'hpc': {'delete_fastq_after_quant': True},
    # Runs are downsampled or hold fewer than 20M reads; the 20M production floor would exclude host samples. Same value as
    # tests/check_tagseq_public.sh uses for its downsampled runs. The pathogen floor is min_assigned_reads.
    'thresholds': {'sample_qc': {'min_reads_on_genes_rnaseq': 100000}},
    # min_assigned_reads is left at the route default. Pathogen differential expression is requested so
    # that the run shows whether it completes or is reported as not estimable.
    'dual_organism': {'enabled': True, 'host': organism('host'),
                      'pathogen': dict(organism('pathogen'), objectives=['gene_expression', 'differential_expression'])},
}
(proj / 'config/config.yaml').write_text(yaml.safe_dump(cfg))
(proj / 'config/thresholds.yaml').write_text((repo / 'config/thresholds.yaml').read_text())
PY

fail_log() {  # $1 = what failed, $2 = project dir, $3 = log
    echo "FAIL: $1"
    tail -60 "$3"
    # The failing job's log is among the most recently written.
    for f in $(ls -t "$2"/logs/*.log "$2"/logs/*/*.log 2>/dev/null | head -6); do
        echo "=== $f ==="; tail -40 "$f"
    done
    exit 1
}
run_case() {  # $1 = project dir (parent or child)
    echo "--- preflight and workflow: $1 ---"
    python3 "$REPO/scripts/preflight.py" -d "$1"
    run_workflow "$1" "$REPO" "$1/run.log" || fail_log "workflow in $1" "$1" "$1/run.log"
}

run_case "$proj"
cp "$proj/results/assignment.tsv" "$proj/assignment_first.tsv"
echo "--- repeat of the assignment ---"
touch "$proj/before_repeat"
SMK_EXTRA="${SMK_EXTRA:-} --forcerun dual_assign" run_workflow "$proj" "$REPO" "$proj/repeat.log" results/assignment.tsv \
    || fail_log "repeat of the assignment in $proj" "$proj" "$proj/repeat.log"
# The repeat proves nothing unless dual_assign ran again for every sample.
expected=$(($(wc -l < "$proj/config/samples.tsv") - 1))
rewritten=$(find "$proj/results/quant" -path '*/dual/assignment.json' -newer "$proj/before_repeat" | wc -l)
if [ "$rewritten" -ne "$expected" ]; then
    echo "FAIL: the repeat did not re-run dual_assign for every sample"
    echo "  expected: $expected results/quant/<sample>/dual/assignment.json written after the repeat started"
    echo "  found:    $((rewritten))"
    exit 1
fi

python3 "$REPO/scripts/dual_organism.py" children --project "$proj"
for role in host pathogen; do
    [ -f "$proj/dual/$role/config/config.yaml" ] || { echo "no $role project (see results/pathogen/status.json)"; continue; }
    run_case "$proj/dual/$role"
done

# Every criterion is evaluated and each mismatch is printed as expected against found.
python3 - "$study" "$proj" "$quantifier" <<'PY'
import csv, json, sys
from pathlib import Path
study, p, quantifier = sys.argv[1], Path(sys.argv[2]), sys.argv[3]
failures = []


def check(what, found, expected, ok=None):
    if not (found == expected if ok is None else ok):
        failures.append(what)
        print(f'FAIL {what}\n  expected: {expected}\n  found:    {found}')


def rows(path):
    with open(path, newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def columns(path):
    with open(path, newline='') as f:
        return next(csv.reader(f, delimiter='\t'))[1:]


def de_tables(root):
    manifest = rows(root / 'de_manifest.tsv') if (root / 'de_manifest.tsv').is_file() else []
    return bool(manifest) and all((root / r['analysis_table']).is_file() and (root / r['analysis_table']).stat().st_size > 0 for r in manifest)


COUNTS = ('reads', 'host', 'pathogen', 'ambiguous', 'unassigned', 'too_many_loci')
sheet = rows(p / 'config/samples.tsv')
samples = [r['sample_id'] for r in sheet]
infected = [r['sample_id'] for r in sheet if r['treatment'] == 'infected']
controls = [s for s in samples if s not in infected]
table = {r['sample']: r for r in rows(p / 'results/assignment.tsv')}
check('assignment.tsv samples', sorted(table), sorted(samples))

# 1. Conservation.
for s in samples:
    c = {k: int(table[s][k]) for k in COUNTS}
    check(f'1 conservation {s}: host + pathogen + ambiguous + unassigned', c['host'] + c['pathogen'] + c['ambiguous'] + c['unassigned'],
          f"{c['reads']} (reads)", c['host'] + c['pathogen'] + c['ambiguous'] + c['unassigned'] == c['reads'] and c['reads'] > 0)

# 2. Every control has a lower pathogen-assigned fraction than every infected sample.
empty = [s for s in samples if int(table[s]['reads']) == 0]
check('2 samples with no reads (a pathogen fraction needs reads)', empty, [])
fraction = {s: None if s in empty else int(table[s]['pathogen']) / int(table[s]['reads']) for s in samples}
check('2 pathogen fraction, highest control below lowest infected',
      {'controls': {s: fraction[s] for s in controls}, 'infected': {s: fraction[s] for s in infected}},
      'max(controls) < min(infected)', not empty and max(fraction[s] for s in controls) < min(fraction[s] for s in infected))

# 3. Count tables.
check('3 host counts.tsv columns', columns(p / 'results/host/counts.tsv'), samples)
status = json.loads((p / 'results/pathogen/status.json').read_text())
floor = status['min_assigned_reads']
enough = [s for s in samples if int(table[s]['pathogen']) >= floor]
check(f'3 pathogen status.json samples (at least {floor} pathogen reads)', status['samples'], enough)
check('3 samples_not_quantified.tsv', [r['sample'] for r in rows(p / 'results/pathogen/samples_not_quantified.tsv')],
      [s for s in samples if s not in enough])
pathogen_counts = p / 'results/pathogen/counts.tsv'
if enough:
    check('3 pathogen counts.tsv columns', columns(pathogen_counts) if pathogen_counts.is_file() else 'no counts.tsv', enough)
else:
    check('3 pathogen project when no sample reaches the floor', status['status'], 'skipped')

# 4. Differential expression.
check('4 host differential expression tables (results/host/de_manifest.tsv)', de_tables(p / 'results/host'), True)
if status['de'] == 'requested':
    check('4 pathogen differential expression tables (results/pathogen/de_manifest.tsv)', de_tables(p / 'results/pathogen'), True)
else:
    check('4 pathogen differential expression', (status['de'], status['reason']), 'not_estimable or skipped, with a reason',
          status['de'] in ('not_estimable', 'skipped') and bool(status['reason']))

# 5. The repeat of the assignment gives the same counts.
first = {r['sample']: {k: r[k] for k in COUNTS} for r in rows(p / 'assignment_first.tsv')}
check('5 repeat assignment counts', {s: {k: table[s][k] for k in COUNTS} for s in table}, first)

record = dict(study=study, quantifier=quantifier, min_assigned_reads=floor, pathogen_status=status,
              samples={s: dict(group='infected' if s in infected else 'control', **{k: int(table[s][k]) for k in COUNTS},
                               pathogen_fraction=fraction[s]) for s in samples})
(p / 'dual_public_summary.json').write_text(json.dumps(record, indent=2) + '\n')
for s, row in record['samples'].items():
    print(s, *(f'{k}={v}' for k, v in row.items()))
print(f"pathogen: {status['status']}, {len(status['samples'])} sample(s) quantified, differential expression {status['de']}"
      + (f" ({status['reason']})" if status['reason'] else ''))
# Without one, the run shows nothing about pathogen quantification.
check(f'infected samples with at least {floor} pathogen reads', [s for s in infected if s in enough], 'at least one',
      any(s in enough for s in infected))
if failures:
    sys.exit(f'{len(failures)} criterion check(s) failed for {study} in {p}')
print(f'{study} dual-organism public check ({quantifier}): all five acceptance criteria held')
PY
echo "Public dual-organism project: $proj"

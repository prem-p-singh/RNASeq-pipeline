#!/usr/bin/env bash
# Synthetic host-plus-pathogen check, run the way the launcher runs it: the parent
# workflow (assign and split), `dual_organism.py children`, then each child project.
# Done once with analysis.quantifier star and once with salmon.
#
# Host genome 40 kb and pathogen genome 30 kb, both one contig named chr1, random
# sequence with one identical 2 kb segment in each. All genes are on the plus strand
# and all reads are 100 nt plus-strand copies with no N bases.
#   host      H1 5001-10000, H2 15001-17000 (no reads), HS 30001-32000 (the shared segment)
#   pathogen  P1 4001-8000;  the shared segment is 20001-22000, in no gene
# Six single-end samples, S1-S3 mock and S4-S6 infected. Per sample:
#   300 reads inside H1                -> host
#    60 reads inside the shared 2 kb   -> ambiguous (in neither count table)
#    40 random 100-mers                -> unassigned
#   200 reads inside P1, infected only -> pathogen
# min_assigned_reads is 100, so the pathogen project holds S4-S6 only.
# A further paired-end project has one library: 30 pairs with both mates in H1 (host)
# and 20 pairs with mate 1 in H1 and mate 2 in P1 (ambiguous, design section 4).
#
# Build the fixture alone (python3 and PyYAML only) and print read counts per class:
#   bash tests/check_dual_workflow.sh --fixture DIR
# KEEP_TMP=1 keeps the working folder and prints its path.
set -euo pipefail
REPO=$(cd "$(dirname "$0")/.." && pwd)

build() {  # $1 = folder to create
python3 - "$1" "$REPO" <<'PY'
import csv, gzip, random, sys
from collections import Counter
from pathlib import Path
import yaml
root, repo = Path(sys.argv[1]).resolve(), Path(sys.argv[2])
rng = random.Random(4102026)
READ = 100
rand = lambda n: ''.join(rng.choices('ACGT', k=n))
rc = lambda s: s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
shared = rand(2000)
# STAR 2.7.11b aborts indexing a genome under about 30 kb with a GTF; each child indexes one genome alone.
host = rand(30000) + shared + rand(8000)
pathogen = rand(20000) + shared + rand(8000)
H1, P1 = (5000, 10000), (4000, 8000)  # 0-based half-open


def pick(seq, lo, hi):
    start = rng.randrange(lo, hi - READ + 1)
    return seq[start:start + READ]


inputs = root / 'inputs'
inputs.mkdir(parents=True)
(inputs / 'host.fa').write_text(f'>chr1\n{host}\n')
(inputs / 'pathogen.fa').write_text(f'>chr1\n{pathogen}\n')
gtf = 'chr1\tfixture\texon\t{}\t{}\t.\t+\t.\tgene_id "{}"; transcript_id "{}_tx";\n'
# H2 has no reads, so the host has one expressed gene and two unexpressed genes of another length. That is
# the shape on which Salmon 1.10.3 with --posBias crashes (FARM jobs 39410501 and 39411257).
(inputs / 'host.gtf').write_text(gtf.format(5001, 10000, 'H1', 'H1') + gtf.format(15001, 17000, 'H2', 'H2')
                                 + gtf.format(30001, 32000, 'HS', 'HS'))
(inputs / 'pathogen.gtf').write_text(gtf.format(4001, 8000, 'P1', 'P1'))
(inputs / 'host.transcripts.fa').write_text(f'>H1_tx\n{host[slice(*H1)]}\n>H2_tx\n{host[15000:17000]}\n>HS_tx\n{shared}\n')
(inputs / 'pathogen.transcripts.fa').write_text(f'>P1_tx\n{pathogen[slice(*P1)]}\n')

record = lambda name, seq: f'@{name}\n{seq}\n+\n{"I" * len(seq)}\n'
single = []
for s in range(6):
    sid = f'S{s + 1}'
    reads = [(f'{sid}_h1_{i}', pick(host, *H1)) for i in range(300)]
    reads += [(f'{sid}_sh_{i}', pick(shared, 0, len(shared))) for i in range(60)]
    reads += [(f'{sid}_rnd_{i}', rand(READ)) for i in range(40)]
    if s >= 3:
        reads += [(f'{sid}_p1_{i}', pick(pathogen, *P1)) for i in range(200)]
    rng.shuffle(reads)  # classes interleaved, so the split depends on STAR keeping input order
    with gzip.open(inputs / f'{sid}.fq.gz', 'wt') as f:
        f.writelines(record(name, seq) for name, seq in reads)
    single.append([sid, str(inputs / f'{sid}.fq.gz'), '', 'mock' if s < 3 else 'infected'])
pairs = []
for i in range(30):
    start = rng.randrange(H1[0], H1[1] - 300 + 1)
    pairs.append((f'PE1_hh_{i}', host[start:start + READ], rc(host[start + 200:start + 300])))
pairs += [(f'PE1_hp_{i}', pick(host, *H1), rc(pick(pathogen, *P1))) for i in range(20)]
rng.shuffle(pairs)
with gzip.open(inputs / 'PE1_R1.fq.gz', 'wt') as a, gzip.open(inputs / 'PE1_R2.fq.gz', 'wt') as b:
    for name, left, right in pairs:
        a.write(record(name + '/1', left)); b.write(record(name + '/2', right))
paired = [['PE1', str(inputs / 'PE1_R1.fq.gz'), str(inputs / 'PE1_R2.fq.gz'), 'infected']]

organism = lambda role, label: {
    'label': label, 'genome_fasta_url': (inputs / f'{role}.fa').as_uri(), 'gtf_url': (inputs / f'{role}.gtf').as_uri(),
    'transcriptome_fasta_url': (inputs / f'{role}.transcripts.fa').as_uri()}
for name, quantifier, layout, sheet in (('star', 'star', 'single', single), ('salmon', 'salmon', 'single', single),
                                        ('paired', 'star', 'paired', paired)):
    proj = root / name
    (proj / 'config').mkdir(parents=True)
    with open(proj / 'config/samples.tsv', 'w') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['sample_id', 'fastq_url', 'fastq_url_r2', 'treatment']); w.writerows(sheet)
    cfg = {
        'project': {'name': f'dual_{name}', 'output_dir': 'results/'},
        'organism': {'scientific_name': 'Fixtura synthetica', 'tax_id': 1, 'kegg_code': None,
                     'orgdb_package': 'org.Fixture.eg.db'},
        'reference': {'accession': 'FIXTURE_DUAL_1', 'cache_dir': None, 'annotation_tsv': {'path': None}},
        'samples': {'seq_type': f'rnaseq_{layout}', 'sheet': 'config/samples.tsv',
                    'expected_libtype': 'SF' if layout == 'single' else 'ISF'},
        # Salmon 1.10.3 can crash in position-bias correction when a length class has no reads: H2 and HS have none.
        'analysis': dict({'objectives': ['gene_expression'], 'quantifier': quantifier},
                         **({'salmon_position_bias': False} if quantifier == 'salmon' else {})),
        'model': {'fixed_effects': '~ treatment', 'random_effects': None, 'primary_factor': 'treatment'},
        'contrasts': [{'id': 'treatment', 'type': 'pairwise', 'factor': 'treatment', 'reverse': True}],
        'downstream': {'run_go': False, 'run_kegg': False, 'run_wgcna': False},
        'orgdb': {'strategy': 'skip'},
        'hpc': {'delete_fastq_after_quant': True},
        'thresholds': {'sample_qc': {'min_reads_on_genes_rnaseq': 1}},
        # Three genes, one expressed, cannot carry a differential-expression fit: both organisms are counts only.
        'dual_organism': {'enabled': True,
                          'host': dict(organism('host', 'hostorg'), objectives=['gene_expression']),
                          'pathogen': dict(organism('pathogen', 'pathorg'), min_assigned_reads=100)},
    }
    (proj / 'config/config.yaml').write_text(yaml.safe_dump(cfg))
    (proj / 'config/thresholds.yaml').write_text((repo / 'config/thresholds.yaml').read_text())

print('sample\th1\tsh\trnd\tp1\thh\thp')
for path in sorted(inputs.glob('*.fq.gz')):
    if path.name.endswith('_R2.fq.gz'):
        continue
    with gzip.open(path, 'rt') as f:
        n = Counter(line.split('_')[1] for line in list(f)[0::4])
    print(path.name.split('.')[0].removesuffix('_R1'), *(n[k] for k in ('h1', 'sh', 'rnd', 'p1', 'hh', 'hp')), sep='\t')
PY
}

if [ "${1:-}" = --fixture ]; then
    build "${2:?usage: check_dual_workflow.sh --fixture DIR}"
    exit 0
fi

. "$REPO/tests/_gate.sh"
gate_runtime r
for tool in STAR featureCounts samtools fastp salmon multiqc; do
    command -v "$tool" >/dev/null || { echo "needs $tool on PATH"; exit 77; }
done
python3 -c 'import matplotlib' 2>/dev/null || { echo "needs matplotlib"; exit 77; }
TASK_TMP=$(mktemp -d)
trap 'if [ -n "${KEEP_TMP:-}" ]; then echo "kept $TASK_TMP"; else rm -rf "$TASK_TMP"; fi' EXIT

run_case() {  # $1 = project dir (parent or child), rest = targets
    local proj=$1; shift
    echo "--- preflight and workflow: ${proj#"$TASK_TMP"/} ---"
    python3 "$REPO/scripts/preflight.py" -d "$proj"
    if ! run_workflow "$proj" "$REPO" "$proj/run.log" "$@"; then
        echo "FAIL: workflow in $proj"
        tail -60 "$proj/run.log"
        # The failing job's log is among the most recently written.
        for f in $(ls -t "$proj"/logs/*.log "$proj"/logs/*/*.log 2>/dev/null | head -6); do
            echo "=== $f ==="; tail -40 "$f"
        done
        exit 1
    fi
}

# Parent outputs with size and modification time, to show that a later step left them alone.
stamp() {
    python3 -c 'import glob, os, sys
for p in sorted(glob.glob(sys.argv[1] + "/results/assignment.tsv") + glob.glob(sys.argv[1] + "/results/quant/*/dual/*")):
    print(p, os.stat(p).st_size, os.stat(p).st_mtime_ns)' "$1"
}
unchanged() {  # $1 = stamp before, $2 = stamp after, $3 = what happened in between
    [ "$1" = "$2" ] && return 0
    echo "FAIL: parent outputs changed after $3 (expected: no change)"
    diff <(echo "$1") <(echo "$2") || true
    exit 1
}

variant() {  # $1 = project, $2 = new config file name, $3 = dual_organism.pathogen key, $4 = value (YAML)
python3 - "$@" <<'PY'
import sys, yaml
from pathlib import Path
proj, name, key, value = sys.argv[1:]
cfg = yaml.safe_load((Path(proj) / 'config/config.yaml').read_text())
cfg['dual_organism']['pathogen'][key] = yaml.safe_load(value)
(Path(proj) / 'config' / name).write_text(yaml.safe_dump(cfg))
PY
}

# Every assertion is evaluated and each mismatch is printed as expected against found.
verify() {  # $1 = mode, $2 = project dir, $3 = mode argument
python3 - "$@" <<'PY'
import csv, gzip, json, sys
from pathlib import Path
mode, proj, arg = sys.argv[1], Path(sys.argv[2]), (sys.argv[3:] or [''])[0]
MOCK, INFECTED = ['S1', 'S2', 'S3'], ['S4', 'S5', 'S6']
failures = []


def check(what, found, expected, ok=None):
    if not (found == expected if ok is None else ok):
        failures.append(what)
        print(f'FAIL {what}\n  expected: {expected}\n  found:    {found}')


def names(path):
    with gzip.open(path, 'rt') as f:
        return [line[1:].split()[0].removesuffix('/1').removesuffix('/2') for line in list(f)[0::4]]


def rows(path):
    with open(path, newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def labels(sample, expected, mates=('',)):
    """assignment.json counts and the split files against the classes in the input read names."""
    d = proj / 'results/quant' / sample / 'dual'
    expected = dict(expected, reads=sum(expected.values()), too_many_loci=0)
    record = json.loads((d / 'assignment.json').read_text())
    check(f'{sample} assignment.json counts', record['counts'], expected)
    check(f'{sample} assignment.json labels', record['labels'], {'host': 'hostorg', 'pathogen': 'pathorg'})
    source = names(proj.parent / 'inputs' / f'{sample}{mates[0] or ""}.fq.gz')
    for role, classes in (('host', ('h1', 'hh')), ('pathogen', ('p1',))):
        want = [n for n in source if n.split('_')[1] in classes]
        for mate in mates:
            got = names(d / f'{role}{mate}.fastq.gz')
            check(f'{sample} {role}{mate}.fastq.gz read names, in input order', f'{len(got)} reads, first {got[:3]}',
                  f'{len(want)} reads, first {want[:3]}', got == want)
    return record


if mode == 'parent':
    for s in MOCK + INFECTED:
        labels(s, dict(host=300, pathogen=200 if s in INFECTED else 0, ambiguous=60, unassigned=40))
    table = rows(proj / 'results/assignment.tsv')
    check('assignment.tsv samples (6 rows)', [r['sample'] for r in table], MOCK + INFECTED)
    for r in table:
        want = dict(reads=600, host=300, pathogen=200, ambiguous=60, unassigned=40, too_many_loci=0) if r['sample'] in INFECTED \
            else dict(reads=400, host=300, pathogen=0, ambiguous=60, unassigned=40, too_many_loci=0)
        check(f"assignment.tsv row {r['sample']}", {k: int(r[k]) for k in want}, want)
        # 60 of 400 or 600 reads is above the default 0.05 warning fraction.
        check(f"assignment.tsv warning {r['sample']}", r['warnings'], 'ambiguous fraction ... above 0.05',
              r['warnings'].startswith('ambiguous fraction') and r['warnings'].endswith('above 0.05'))

elif mode == 'paired':
    # Design section 4: mates pointing to different organisms make the pair ambiguous.
    # If the 20 cross-organism pairs come back unassigned, STAR left them unaligned: with default
    # filters one aligned mate of 100 nt is under 66% of the 200 nt pair. That is a pipeline or
    # design question (see the task 5 report), not a fixture error.
    labels('PE1', dict(host=30, pathogen=0, ambiguous=20, unassigned=0), mates=('_R1', '_R2'))
    check('paired assignment.tsv samples', [r['sample'] for r in rows(proj / 'results/assignment.tsv')], ['PE1'])

elif mode == 'status':
    status = json.loads((proj / 'results/pathogen/status.json').read_text())
    missed = rows(proj / 'results/pathogen/samples_not_quantified.tsv')
    child = proj / 'dual/pathogen/config'
    if arg == 'skipped':  # min_assigned_reads 201: no sample qualifies
        check('status.json with min_assigned_reads 201', {k: status[k] for k in ('status', 'de', 'samples', 'not_quantified', 'min_assigned_reads')},
              dict(status='skipped', de='skipped', samples=[], not_quantified=MOCK + INFECTED, min_assigned_reads=201))
        check('samples_not_quantified.tsv with min_assigned_reads 201', [r['sample'] for r in missed], MOCK + INFECTED)
        check('pathogen child config removed when skipped', (child / 'config.yaml').exists(), False)
    else:
        check(f'status.json ({arg})', {k: status[k] for k in ('status', 'de', 'samples', 'not_quantified', 'min_assigned_reads')},
              dict(status='quantified', de=arg, samples=INFECTED, not_quantified=MOCK, min_assigned_reads=100))
        check(f'status.json reason ({arg})', status['reason'], "names level 'mock'" if arg == 'not_estimable' else '',
              'mock' in status['reason'] if arg == 'not_estimable' else status['reason'] == '')
        check('samples_not_quantified.tsv', [(r['sample'], r['pathogen_reads']) for r in missed], [(s, '0') for s in MOCK])
        check('samples_not_quantified.tsv reasons', [r['reason'] for r in missed], ['below min_assigned_reads 100'] * 3)
        check('pathogen child samples', [r['sample_id'] for r in rows(child / 'samples.tsv')], INFECTED)
        check('host child samples', [r['sample_id'] for r in rows(proj / 'dual/host/config/samples.tsv')], MOCK + INFECTED)
        # A design that cannot be estimated leaves the pathogen project as counts only.
        import yaml
        check('pathogen child objectives', yaml.safe_load((child / 'config.yaml').read_text())['analysis']['objectives'], ['gene_expression'])

elif mode == 'children':  # arg = star or salmon
    for role, samples, genes in (('host', MOCK + INFECTED, {'H1': 300, 'H2': 0, 'HS': 0}), ('pathogen', INFECTED, {'P1': 200})):
        out, n = proj / 'results' / role, sum(genes.values())
        for s in samples:
            q = out / 'quant' / s
            ledger = json.loads((q / 'preprocessing.json').read_text())
            # Pass-through preprocessing: fastp must not drop or change a read (reads carry no N bases).
            check(f'{role} {s} preprocessing.json passthrough, reads before, after',
                  (ledger['policy']['passthrough'], ledger['reads_before'], ledger['reads_after']), (True, n, n))
            if arg == 'star':
                with open(q / 'star/gene_counts.tsv') as f:
                    got = {r[0]: int(r[1]) for r in list(csv.reader(f, delimiter='\t'))[1:]}
                check(f'{role} {s} star/gene_counts.tsv', got, genes)
                m = json.loads((q / 'star/metrics.json').read_text())
                check(f'{role} {s} star/metrics.json strandedness, reads, mapped, assigned',
                      (m['strandedness'], m['num_reads_processed'], m['num_mapped'], m['num_assigned']), (1, n, n, n))
            else:
                m = json.loads((q / 'metrics.json').read_text())
                check(f'{role} {s} Salmon metrics.json num_reads_processed', m['num_reads_processed'], n)
                check(f'{role} {s} Salmon bias flags, no --posBias', m['salmon_bias_flags'], '--gcBias --seqBias')
                check(f'{role} {s} quantified by Salmon, not STAR', ((q / 'quant.sf').is_file(), (q / 'star').exists()), (True, False))
        table = rows(out / 'counts.tsv')
        check(f'{role} counts.tsv columns', list(table[0]), ['gene_id'] + samples)
        check(f'{role} counts.tsv genes', sorted(r['gene_id'] for r in table), sorted(genes))
        for r in table:
            want = genes.get(r['gene_id'])
            for s in samples:
                # STAR counts are exact. Salmon counts are estimates: 3% as in tests/fixtures/fastq_project.py; zero stays zero.
                slack = 0 if arg == 'star' else 0.03 * want
                check(f"{role} counts.tsv {r['gene_id']} {s}", r[s], f'{want}' + (f' +/- {slack:g}' if slack else ''),
                      abs(float(r[s]) - want) <= slack)
        prov = json.loads((out / 'counts_provenance.json').read_text())
        if arg == 'star':
            check(f'{role} counts_provenance.json', (prov.get('quantifier'), prov.get('count_unit')), ('STAR_featureCounts', 'reads'))
        else:
            check(f'{role} counts_provenance.json counts_from_abundance', prov['counts_from_abundance'], 'lengthScaledTPM')

if failures:
    sys.exit(f'{len(failures)} assertion(s) failed in {mode} {arg}'.rstrip() + f' for {proj}')
print('ok:', mode, *([arg] if arg else []), proj.name)
PY
}

build "$TASK_TMP"

for quantifier in star salmon; do
    proj="$TASK_TMP/$quantifier"
    run_case "$proj"
    verify parent "$proj"
    before=$(stamp "$proj")

    python3 "$REPO/scripts/dual_organism.py" children --project "$proj"
    verify status "$proj" not_requested
    for role in host pathogen; do
        test -f "$proj/dual/$role/config/config.yaml" || { echo "FAIL: no $role child project in $proj/dual"; exit 1; }
        run_case "$proj/dual/$role"
    done
    verify children "$proj" "$quantifier"
    unchanged "$before" "$(stamp "$proj")" "the child projects ran"

    # An unchanged relaunch of the parent assigns nothing again.
    run_workflow "$proj" "$REPO" "$proj/resume.log" || { echo "FAIL: parent relaunch in $proj"; tail -60 "$proj/resume.log"; exit 1; }
    unchanged "$before" "$(stamp "$proj")" "an unchanged parent relaunch"

    # Pathogen differential expression needs 2 quantified samples per level; mock has none.
    variant "$proj" de.yaml objectives '[gene_expression, differential_expression]'
    python3 "$REPO/scripts/dual_organism.py" children --project "$proj" --configfile config/de.yaml
    verify status "$proj" not_estimable

    # Only the floor changes: the assignment is not redone and no sample reaches 201 pathogen reads.
    variant "$proj" floor.yaml min_assigned_reads 201
    python3 "$REPO/scripts/dual_organism.py" children --project "$proj" --configfile config/floor.yaml
    verify status "$proj" skipped
    unchanged "$before" "$(stamp "$proj")" "children with another min_assigned_reads"
done

run_case "$TASK_TMP/paired" results/assignment.tsv
verify paired "$TASK_TMP/paired"

echo 'dual-organism workflow: exact labels, split reads, per-organism counts with STAR and Salmon, pathogen status and unchanged relaunch passed'

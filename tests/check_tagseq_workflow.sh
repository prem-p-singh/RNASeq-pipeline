#!/usr/bin/env bash
# Synthetic QuantSeq route check: one 10 kb contig, gene A on the plus strand and
# gene B on the minus strand, six single-end samples in kit read geometry.
#   FWD:     12 nt + 50 nt plus-strand insert; 100 reads in A (sense) and 20 in B
#            (antisense to B) -> A=100, B=0, no strand warning
#   REV:     the same reads declared REV -> A=0, B=20 (positive control); 20/120
#            assigned is below the 20% floor, so metrics.json warns of a strand mismatch
#   FWD-UMI: 6 nt UMI + TATA + 12 nt + insert; 20 molecules x 5 identical PCR
#            duplicates -> 100 reads but 20 molecules on A, plus one read too
#            short after UMI, spacer and front trim, lost in the preprocessing ledger
# The B reads are in the FWD/REV libraries only; the UMI library has A reads only.
# Without umi_tools, FWD and REV still run and assert, then the script exits 77
# so run_all.sh reports the UMI case as SKIP (and --strict fails).
set -euo pipefail
REPO=$(cd "$(dirname "$0")/.." && pwd)
. "$REPO/tests/_gate.sh"
gate_runtime r
for tool in STAR featureCounts samtools fastp; do
    command -v "$tool" >/dev/null || { echo "needs $tool on PATH"; exit 77; }
done
TASK_TMP=$(mktemp -d)
trap 'rm -rf "$TASK_TMP"' EXIT

build() {  # $1 = project dir, $2 = kit name
python3 - "$1" "$2" "$REPO" <<'PY'
import csv, gzip, json, random, sys
from pathlib import Path
import yaml
proj, kit, repo = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3])
rng = random.Random(27092026)
genome = ''.join(rng.choices('ACGT', k=30000))  # STAR 2.7.11b aborts indexing a 10 kb genome with a GTF (job 39023968)
(proj / 'inputs').mkdir(parents=True); (proj / 'config').mkdir()
(proj / 'inputs/genome.fa').write_text(f'>chr1\n{genome}\n')
(proj / 'inputs/genes.gtf').write_text(
    'chr1\tfixture\texon\t1001\t3000\t.\t+\t.\tgene_id "A"; transcript_id "A_tx";\n'
    'chr1\tfixture\texon\t5001\t7000\t.\t-\t.\tgene_id "B"; transcript_id "B_tx";\n')
umi = kit == 'QuantSeq FWD-UMI'
umis = []
while len(umis) < 21:  # pairwise >= 3 mismatches, so no directional UMI merging
    u = ''.join(rng.choices('ACGT', k=6))
    if all(sum(a != b for a, b in zip(u, v)) >= 3 for v in umis):
        umis.append(u)
sheet, expected = [], {}
for s in range(6):
    sid = f'S{s+1}'
    reads = []
    for m in range(20 if umi else 100):
        start = rng.randrange(1000, 2950)  # inside A (1-based 1001-3000)
        seq = (umis[m] + 'TATA' if umi else '') + ''.join(rng.choices('ACGT', k=12)) + genome[start:start + 50]
        reads += [(f'{sid}_m{m}_d{d}', seq) for d in range(5 if umi else 1)]
    for b in range(0 if umi else 20):
        start = rng.randrange(5000, 6950)  # inside B (1-based 5001-7000), plus strand = antisense to B
        reads.append((f'{sid}_b{b}', ''.join(rng.choices('ACGT', k=12)) + genome[start:start + 50]))
    if umi:  # 6 + 4 + 12 + 10 nt: 10 nt left after processing, below minimum_length 20
        reads.append((f'{sid}_short', umis[20] + 'TATA' + ''.join(rng.choices('ACGT', k=22))))
        expected[sid] = {seq[:6]: seq[-50:] for _, seq in reads[:-1]}  # UMI -> insert
    path = proj / f'inputs/{sid}.fq.gz'
    with gzip.open(path, 'wt') as f:
        f.writelines(f'@{name}\n{seq}\n+\n{"I" * len(seq)}\n' for name, seq in reads)
    sheet.append([sid, str(path), '', 'control' if s < 3 else 'treated'])
with open(proj / 'config/samples.tsv', 'w') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['sample_id', 'fastq_url', 'fastq_url_r2', 'treatment']); w.writerows(sheet)
cfg = {
    'project': {'name': 'tagseq_' + kit.split()[-1].replace('-', '_').lower(), 'output_dir': 'results/'},
    'organism': {'scientific_name': 'Fixtura synthetica', 'tax_id': 1, 'kegg_code': None,
                 'orgdb_package': 'org.Fixture.eg.db'},
    'reference': {'accession': 'FIXTURE_TAGSEQ_1', 'cache_dir': None, 'transcriptome_fasta_url': None,
                  'genome_fasta_url': (proj / 'inputs/genome.fa').as_uri(),
                  'gtf_url': (proj / 'inputs/genes.gtf').as_uri(), 'annotation_tsv': {'path': None}},
    'samples': {'seq_type': 'tagseq', 'sheet': 'config/samples.tsv'},
    'analysis': {'objectives': ['gene_expression'], 'quantifier': 'auto', 'kit': kit},
    'model': {'fixed_effects': '~ treatment', 'random_effects': None, 'primary_factor': 'treatment'},
    'contrasts': [{'id': 'treatment', 'type': 'pairwise', 'factor': 'treatment', 'reverse': True}],
    'downstream': {'run_go': False, 'run_kegg': False, 'run_wgcna': False},
    'orgdb': {'strategy': 'skip'},
    'hpc': {'delete_fastq_after_quant': False},
    'thresholds': {'sample_qc': {'min_reads_on_genes_tagseq': 1}},
}
(proj / 'config/config.yaml').write_text(yaml.safe_dump(cfg))
(proj / 'config/thresholds.yaml').write_text((repo / 'config/thresholds.yaml').read_text())
(proj / 'expected.json').write_text(json.dumps(expected))
PY
}

run_case() {  # $1 = project dir, rest = targets
    local proj=$1; shift
    python3 "$REPO/scripts/preflight.py" -d "$proj"
    if ! run_workflow "$proj" "$REPO" "$proj/run.log" "$@"; then
        tail -80 "$proj/run.log"
        for f in "$proj"/logs/star/*.log "$proj"/logs/preprocess/*.log "$proj"/logs/star_index.log; do
            test ! -f "$f" || tail -40 "$f"
        done
        exit 1
    fi
}

per_sample() {  # $1 = product under results/quant/{sample}/
    for s in S1 S2 S3 S4 S5 S6; do echo "results/quant/$s/$1"; done
}

for kit in FWD REV; do
    build "$TASK_TMP/$kit" "QuantSeq $kit"
done
run_case "$TASK_TMP/FWD" results/counts.tsv
run_case "$TASK_TMP/REV" $(per_sample star/metrics.json)

python3 - "$TASK_TMP" <<'PY'
import csv, json, sys
from pathlib import Path
root = Path(sys.argv[1])
samples = [f'S{i}' for i in range(1, 7)]

def table(path):
    with open(path) as f:
        return {r[0]: int(r[1]) for r in list(csv.reader(f, delimiter='\t'))[1:]}

def star(proj, s):
    q = proj / 'results/quant' / s / 'star'
    return table(q / 'gene_counts.tsv'), json.loads((q / 'metrics.json').read_text())

for s in samples:
    fwd, m = star(root / 'FWD', s)
    assert fwd == {'A': 100, 'B': 0}, (s, fwd)
    assert m['strandedness'] == 1 and m['mapping_rate'] == 1.0 and m['num_assigned'] == 100, m
    assert m['assignment_summary']['Unassigned_NoFeatures'] == 20, m
    assert m['warnings'] == [], m['warnings']
    rev, m = star(root / 'REV', s)
    assert rev == {'A': 0, 'B': 20}, (s, rev)
    assert m['strandedness'] == 2 and m['mapping_rate'] == 1.0 and m['num_assigned'] == 20, m
    assert m['assignment_summary']['Unassigned_NoFeatures'] == 100, m
    assert len(m['warnings']) == 1 and 'possible strand mismatch' in m['warnings'][0], m['warnings']
print('QuantSeq strand: FWD counts A=100, B=0; REV on the same reads counts A=0, B=20 and warns')

p = root / 'FWD'
with open(p / 'results/counts.tsv') as f:
    rows = {r['gene_id']: r for r in csv.DictReader(f, delimiter='\t')}
assert all(int(rows['A'][s]) == 100 and int(rows['B'][s]) == 0 for s in samples), rows
prov = json.loads((p / 'results/counts_provenance.json').read_text())
assert prov['quantifier'] == 'STAR_featureCounts' and prov['count_unit'] == 'reads', prov
PY

if ! command -v umi_tools >/dev/null; then
    echo "FWD and REV passed; FWD-UMI case needs umi_tools on PATH"
    exit 77
fi
build "$TASK_TMP/FWD-UMI" "QuantSeq FWD-UMI"
run_case "$TASK_TMP/FWD-UMI" results/counts.tsv

python3 - "$TASK_TMP" <<'PY'
import csv, gzip, json, sys
from pathlib import Path
root = Path(sys.argv[1])
samples = [f'S{i}' for i in range(1, 7)]

def table(path):
    with open(path) as f:
        return {r[0]: int(r[1]) for r in list(csv.reader(f, delimiter='\t'))[1:]}

p = root / 'FWD-UMI'
expected = json.loads((p / 'expected.json').read_text())
for s in samples:
    reads = table(p / 'results/quant' / s / 'star/gene_counts.tsv')
    assert reads == {'A': 100, 'B': 0}, (s, reads)
    assert table(p / 'results/quant' / s / 'umi/gene_counts.tsv') == {'A': 20, 'B': 0}, s
    ledger = json.loads((p / 'results/quant' / s / 'preprocessing.json').read_text())
    assert ledger['reads_before'] == 101 and ledger['reads_after'] == 100, ledger
    assert ledger['size_selection']['too_short_reads'] == 1, ledger['size_selection']
    with gzip.open(p / 'results/quant' / s / f'{s}.trim.fastq.gz', 'rt') as f:
        lines = f.read().splitlines()
    assert len(lines) == 400, len(lines)
    for name, seq in zip(lines[0::4], lines[1::4]):
        insert = expected[s][name.split()[0].rsplit(':', 1)[1]]  # UMI moved into the read name
        # UMI, TATA and 12 nt removed; tolerate a few 3' bases clipped as adapter
        assert len(seq) >= 40 and insert.startswith(seq), (name, seq)
with open(p / 'results/counts.tsv') as f:
    rows = {r['gene_id']: r for r in csv.DictReader(f, delimiter='\t')}
assert all(int(rows['A'][s]) == 20 and int(rows['B'][s]) == 0 for s in samples), rows
prov = json.loads((p / 'results/counts_provenance.json').read_text())
assert prov['quantifier'] == 'STAR_featureCounts_umi_tools' and prov['count_unit'] == 'molecules', prov
print('QuantSeq FWD-UMI: 100 reads collapse to 20 molecules; one short read recorded as lost')
PY

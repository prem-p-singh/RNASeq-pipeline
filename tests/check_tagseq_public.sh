#!/usr/bin/env bash
# Public-data check for one QuantSeq kit (FWD, REV or FWD-UMI).
#
# Separate from tests/check_public.sh, which hardcodes the bulk GSE110004 IDs.
# Runs on data already prepared by tests/fixtures/prepare_tagseq_public.py (six
# downsampled runs, two groups of three); the workflow fetches the pinned Ensembl
# genome and GTF itself. Skips (exit 77) without prepared data; not in run_all.sh.
#
# Run:
#   python3 tests/fixtures/prepare_tagseq_public.py "QuantSeq FWD" DATA_DIR    # once
#   bash tests/check_tagseq_public.sh --fwd DATA_DIR
# (--rev with "QuantSeq REV", --umi with "QuantSeq FWD-UMI")
set -euo pipefail
REPO=$(cd "$(dirname "$0")/.." && pwd)
usage="Usage: check_tagseq_public.sh --fwd|--rev|--umi PREPARED_DATA_DIR"
case "${1:-}" in
    --fwd) kit="QuantSeq FWD" ;;
    --rev) kit="QuantSeq REV" ;;
    --umi) kit="QuantSeq FWD-UMI" ;;
    *) echo "$usage" >&2; exit 2 ;;
esac
data=${2:-}
[ -n "$data" ] || { echo "$usage" >&2; exit 2; }
if [ ! -f "$data/samples.tsv" ]; then
    echo "needs prepared public data at $data (run tests/fixtures/prepare_tagseq_public.py \"$kit\" $data first)"
    exit 77
fi

. "$REPO/tests/_gate.sh"
gate_runtime r

proj=${RNASEQ_TAGSEQ_PUBLIC_PROJECT:-$(mktemp -d)/tagseq_public}
python3 - "$kit" "$data" "$proj" "$REPO" <<'PY'
import csv, json, sys
from pathlib import Path
import yaml
kit, data, proj, repo = sys.argv[1], Path(sys.argv[2]).resolve(), Path(sys.argv[3]), Path(sys.argv[4])
ledger = json.loads((data / 'download_manifest.json').read_text())
if ledger['kit'] != kit or ledger['status'] != 'downloaded_not_qualified':
    raise SystemExit(f'{data} holds {ledger["kit"]!r} with status {ledger["status"]!r}; expected {kit!r} downloaded_not_qualified')
profile = json.loads((data / 'profile.json').read_text())
with open(data / 'samples.tsv') as f:
    rows = list(csv.DictReader(f, delimiter='\t'))
groups = sorted({r['condition'] for r in rows})
assert len(rows) == 6 and len(groups) == 2, (len(rows), groups)
(proj / 'config').mkdir(parents=True)
with open(proj / 'config/samples.tsv', 'w') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['sample_id', 'fastq_url', 'fastq_url_r2', 'treatment'])
    w.writerows([r['sample_id'], r['fastq_path'], '', r['condition']] for r in rows)
ref = profile['reference']
cfg = {
    'project': {'name': 'tagseq_public_' + profile['study'], 'output_dir': 'results/'},
    'organism': {'scientific_name': profile['organism'],
                 'tax_id': {'Mus musculus': 10090, 'Homo sapiens': 9606}[profile['organism']],
                 'kegg_code': None, 'orgdb_package': None},
    'reference': {'accession': f'ensembl_{ref["ensembl_release"]}_{profile["organism"].replace(" ", "_")}',
                  'cache_dir': None, 'transcriptome_fasta_url': None,
                  'genome_fasta_url': ref['genome_fasta_url'], 'gtf_url': ref['gtf_url'],
                  'annotation_tsv': {'path': None}},
    'samples': {'seq_type': 'tagseq', 'sheet': 'config/samples.tsv'},
    'analysis': {'objectives': ['gene_expression', 'differential_expression'], 'quantifier': 'auto', 'kit': kit},
    'model': {'fixed_effects': '~ treatment', 'random_effects': None, 'primary_factor': 'treatment'},
    'contrasts': [{'id': 'group', 'type': 'pairwise', 'factor': 'treatment', 'reverse': True}],
    'downstream': {'run_go': False, 'run_kegg': False, 'run_wgcna': False},
    'orgdb': {'strategy': 'skip'},
    'hpc': {'delete_fastq_after_quant': False},
    # Runs are downsampled to at most 2M reads; the 3M production floor would exclude all of them.
    'thresholds': {'sample_qc': {'min_reads_on_genes_tagseq': 100000}},
}
(proj / 'config/config.yaml').write_text(yaml.safe_dump(cfg))
(proj / 'config/thresholds.yaml').write_text((repo / 'config/thresholds.yaml').read_text())
PY
python3 "$REPO/scripts/preflight.py" -d "$proj"
if ! run_workflow "$proj" "$REPO" "$proj/run.log"; then
    tail -80 "$proj/run.log"; exit 1
fi

python3 - "$kit" "$proj" <<'PY'
import csv, json, sys
from pathlib import Path
kit, p = sys.argv[1], Path(sys.argv[2])
umi = kit == 'QuantSeq FWD-UMI'
with open(p / 'config/samples.tsv') as f:
    samples = [r['sample_id'] for r in csv.DictReader(f, delimiter='\t')]
agg = json.loads((p / 'metrics/aggregate.json').read_text())
assert agg['n_samples'] == 6 and agg['n_genes'] > 1000, agg
with open(p / 'results/counts.tsv') as f:
    assert next(csv.reader(f, delimiter='\t'))[1:] == samples, 'samples missing from counts'
prov = json.loads((p / 'results/counts_provenance.json').read_text())
assert prov['count_unit'] == ('molecules' if umi else 'reads'), prov
manifest = list(csv.DictReader((p / 'results/de_manifest.tsv').open(), delimiter='\t'))
assert manifest and all((p / 'results' / r['analysis_table']).stat().st_size > 0 for r in manifest), manifest
assert (p / 'results/qc_report/qc_charts.html').stat().st_size > 0
record = {}
for s in samples:
    m = json.loads((p / 'results/quant' / s / 'star/metrics.json').read_text())
    assert m['strandedness'] == (2 if kit == 'QuantSeq REV' else 1) and m['num_assigned'] > 100000, m
    row = dict(reads=m['num_reads_processed'], mapping_rate=round(m['mapping_rate'], 4),
               assignment_rate=round(m['assignment_rate'], 4))
    if umi:
        with open(p / 'results/quant' / s / 'umi/gene_counts.tsv') as f:
            molecules = sum(int(r['count']) for r in csv.DictReader(f, delimiter='\t'))
        assert 0 < molecules <= m['num_assigned'], (s, molecules, m['num_assigned'])
        row.update(molecules=molecules, duplication_rate=round(1 - molecules / m['num_assigned'], 4))
    record[s] = row
(p / 'tagseq_public_summary.json').write_text(json.dumps(dict(kit=kit, samples=record), indent=2) + '\n')
for s, row in record.items():
    print(s, *(f'{k}={v}' for k, v in row.items()))
print(f'{kit} public check: six samples reached counts, DE tables and the QC report')
PY
echo "Public TAG-seq project: $proj"

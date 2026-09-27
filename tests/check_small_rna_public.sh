#!/usr/bin/env bash
# Public-data qualification for the small-RNA animal and plant routes.
#
# Separate from tests/check_public.sh, which hardcodes the bulk GSE110004
# sample IDs and downloads on every run. This check instead runs against data
# already prepared by tests/fixtures/prepare_small_rna_public.py (animal) or
# tests/fixtures/prepare_plant_small_rna_public.py (plant): no network here.
# It skips (exit 77) when that prepared data directory does not exist yet, so
# it never blocks tests/run_all.sh and is not on that suite's required list.
#
# Run:
#   python3 tests/fixtures/prepare_small_rna_public.py DATA_DIR        # once
#   bash tests/check_small_rna_public.sh --animal DATA_DIR
#
#   python3 tests/fixtures/prepare_plant_small_rna_public.py DATA_DIR  # once
#   bash tests/check_small_rna_public.sh --plant DATA_DIR
set -euo pipefail
REPO=$(cd "$(dirname "$0")/.." && pwd)

mode=${1:-}
data=${2:-}
case "$mode" in
    --animal|--plant) ;;
    *) echo "Usage: check_small_rna_public.sh --animal|--plant PREPARED_DATA_DIR" >&2; exit 2 ;;
esac
[ -n "$data" ] || { echo "Usage: check_small_rna_public.sh --animal|--plant PREPARED_DATA_DIR" >&2; exit 2; }

if [ ! -d "$data" ]; then
    echo "needs prepared public data at $data (run tests/fixtures/prepare_small_rna_public.py or prepare_plant_small_rna_public.py first)"
    exit 77
fi

. "$REPO/tests/_gate.sh"
gate_runtime r

proj=${RNASEQ_SMALL_RNA_PUBLIC_PROJECT:-$(mktemp -d)/small_rna_public}
python3 "$REPO/tests/fixtures/small_rna_public_project.py" "$mode" "$data" "$proj" "$REPO"
python3 "$REPO/scripts/preflight.py" -d "$proj"
if ! run_workflow "$proj" "$REPO" "$proj/run.log"; then
    tail -80 "$proj/run.log"; exit 1
fi

if [ "$mode" = "--plant" ]; then
    # Reproducibility: a forced repeat of the joint ShortStack run must give an
    # identical count matrix (acceptance_plan / acceptance_thresholds.repeat_run_identical_counts).
    cp "$proj/results/small_rna_counts.tsv" "$proj/small_rna_counts.first.tsv"
    SMK_EXTRA="-f" run_workflow "$proj" "$REPO" "$proj/rerun.log" results/small_rna_counts.tsv
    if ! cmp -s "$proj/small_rna_counts.first.tsv" "$proj/results/small_rna_counts.tsv"; then
        echo "Repeat ShortStack run produced a different count matrix" >&2; exit 1
    fi
fi

python3 - "$mode" "$data" "$proj" "$REPO" <<'PY'
import csv
import gzip
import json
import math
import statistics
import sys
from pathlib import Path

mode, data, proj, repo = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3]), Path(sys.argv[4])


def spearman(xs, ys):
    def ranks(values):
        order = sorted(range(len(values)), key=lambda i: values[i])
        r = [0.0] * len(values)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
                j += 1
            avg = (i + j) / 2 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    return statistics.correlation(ranks(xs), ranks(ys))


if mode == "--animal":
    manifest = json.loads((repo / "tests/fixtures/small_rna_public.json").read_text())
    thresholds = manifest["acceptance_thresholds"]
    min_reads = thresholds["min_reads_per_pool_for_ratio_check"]
    expected = {}
    with open(data / "expected_ratios.tsv") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            expected[row["feature_id"]] = float(row["expected_A_over_B"])

    observed, published = {}, {}
    for run in manifest["runs"]:
        sample = run["pool"]
        with open(proj / f"results/quant/{sample}/small_rna/counts.tsv") as f:
            observed[sample] = {r["gene_id"].removeprefix("pub-"): int(r["count"])
                                 for r in csv.DictReader(f, delimiter="\t")}
        metrics = json.loads((proj / f"results/quant/{sample}/small_rna/metrics.json").read_text())
        total = sum(metrics["alignment_summary"].values())
        assert metrics["num_reads_processed"] == total, \
            f"{sample}: count conservation failed, {metrics['num_reads_processed']} reads in vs {total} kept+lost"
        pub = {}
        with gzip.open(data / f"{run['accession']}.published_counts.txt.gz", "rt") as f:
            for line in f:
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 2 and fields[-1].lstrip("-").isdigit():
                    pub[fields[0]] = int(fields[-1])
        published[sample] = pub

    pools = [run["pool"] for run in manifest["runs"]]
    assert len(pools) == 2, f"animal fixture expects exactly two pools, got {pools}"
    a, b = pools
    ratio_diffs = []
    for feature, exp_ratio in expected.items():
        ca, cb = observed[a].get(feature, 0), observed[b].get(feature, 0)
        if ca >= min_reads and cb >= min_reads:
            ratio_diffs.append(abs(math.log2(ca / cb) - math.log2(exp_ratio)))
    assert ratio_diffs, "no Table S2 sequence reached the read-depth floor in both pools"
    median_diff = statistics.median(ratio_diffs)
    assert median_diff <= thresholds["max_median_abs_log2_ratio_diff"], \
        f"median |log2 ratio diff| {median_diff:.3f} exceeds {thresholds['max_median_abs_log2_ratio_diff']} " \
        f"over {len(ratio_diffs)} qualifying sequences"

    for sample in pools:
        shared = sorted(set(observed[sample]) & set(published[sample]))
        assert len(shared) >= 10, f"{sample}: only {len(shared)} IDs shared with published per-run counts"
        rho = spearman([observed[sample][k] for k in shared], [published[sample][k] for k in shared])
        assert rho >= thresholds["min_spearman_vs_published_counts"], \
            f"{sample}: Spearman {rho:.3f} vs published counts below {thresholds['min_spearman_vs_published_counts']} ({len(shared)} shared IDs)"

    print(f"small-RNA animal public check: {len(ratio_diffs)} sequences at >= {min_reads} reads/pool, "
          f"median |log2 ratio diff| {median_diff:.3f}, "
          + ", ".join(f"{s} Spearman {spearman([observed[s][k] for k in sorted(set(observed[s]) & set(published[s]))], [published[s][k] for k in sorted(set(observed[s]) & set(published[s]))]):.3f}" for s in pools)
          + " passed")

else:
    manifest = json.loads((repo / "tests/fixtures/plant_small_rna_public.json").read_text())
    thresholds = manifest["acceptance_thresholds"]

    with open(proj / "results/small_rna_loci.tsv") as f:
        loci = list(csv.DictReader(f, delimiter="\t"))
    assert loci, "no ShortStack loci produced"
    dicer_calls = {r["DicerCall"] for r in loci}
    missing_calls = set(thresholds["required_dicer_calls"]) - dicer_calls
    assert not missing_calls, f"DicerCall values missing from the cohort: {sorted(missing_calls)}"

    for family in thresholds["known_mirna_families"]:
        hit = any(family in r.get("Known_miRNAs", "") or (r["MIRNA"] == "Y" and family in r["Name"]) for r in loci)
        assert hit, f"{family}: not annotated MIRNA and not matched in Known_miRNAs"

    total = kept = 0
    for run in manifest["runs"]:
        sample = run["accession"]
        with gzip.open(proj / f"results/quant/{sample}/{sample}.trim.fastq.gz", "rt") as f:
            for i, line in enumerate(f):
                if i % 4 == 1:
                    total += 1
                    kept += 20 <= len(line.rstrip("\n")) <= 24
    assert total, "no trimmed reads found"
    fraction = kept / total
    assert fraction >= thresholds["min_trimmed_20_24nt_fraction"], \
        f"20-24 nt fraction {fraction:.3f} below {thresholds['min_trimmed_20_24nt_fraction']} ({kept}/{total} reads)"

    print(f"small-RNA plant public check: {len(loci)} loci, DicerCall {sorted(dicer_calls)}, "
          f"20-24 nt fraction {fraction:.3f} ({kept}/{total}), repeat run identical, passed")
PY
echo "Small-RNA public project ($mode): $proj"

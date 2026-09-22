#!/usr/bin/env bash
# V09: WGCNA must actually run, and behave the way its metrics claim.
#
# Master plan 12: "WGCNA | Memory behavior, sample/gene filtering, documented
# parameters and stable output identity". Until now Stage 5 had only ever taken
# its skip gate, so none of that had been observed. The fixture here has 20
# samples (over thresholds.wgcna.min_samples) and 2000 genes carrying four
# planted co-expression modules.
#
# Three runs, because the claims need more than one:
#   a, b   the same small memory allocation, run independently -> output identity
#   big    enough memory for a single block                    -> memory behaviour
#
# mem_mb is what sizes the blocks. R16 was "WGCNA forced one huge block": the
# fix derives maxBlockSize from the allocation, and the only way to see that is
# to vary the allocation and watch the block count move.
#
# Slow by the standards of this suite, about 90 seconds, because it runs
# blockwiseModules four times. Exits 77 (skipped) without snakemake and the R
# stack.
#
# Run:  bash tests/check_network.sh
set -uo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
FIXTURE="$REPO/tests/fixtures/tiny_project.py"
. "$REPO/tests/_gate.sh"
gate_runtime r
if ! Rscript -e 'q(status = if (requireNamespace("WGCNA", quietly = TRUE)) 0 else 1)' \
        >/dev/null 2>&1; then
    echo "needs R package: WGCNA"
    exit 77
fi

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

# Fresh build per run rather than a copy: copying rewrites mtimes, and Snakemake
# then decides the reference index is out of date and tries to rebuild it with a
# salmon that is not here. The fixture is seeded, so the inputs are identical.
one_run() {  # $1 = dir, $2 = mem_mb
    python3 "$FIXTURE" build "$1" "$REPO" rnaseq wgcna \
        || { echo "FAIL: fixture build" >&2; return 1; }
    SMK_EXTRA="--set-resources wgcna:mem_mb=$2" \
        run_workflow "$1" "$REPO" "$1.log" results/wgcna_done.flag
    local rc=$?
    if [ $rc -ne 0 ]; then
        echo "FAIL: stage 5 did not complete at mem_mb=$2 (exit $rc)" >&2
        tail -30 "$1.log" >&2
        return 1
    fi
}

one_run "$TMP/a" 20    || exit 1
one_run "$TMP/b" 20    || exit 1
one_run "$TMP/big" 16000 || exit 1

python3 "$FIXTURE" network "$TMP/a" "$TMP/b" "$TMP/big" || exit 1
echo "check_network.sh: WGCNA ran, and its memory, parameters and identity hold"

#!/usr/bin/env bash
# V04: bulk gene aggregation must agree numerically with an independent route.
#
# Master plan 12: "Bulk gene/transcript aggregation | Numerical agreement with
# an independent trusted route within declared tolerances". check_stages.sh runs
# the tag branch, where countsFromAbundance is "no" and gene counts are summed
# reads. This runs the BULK branch, where it is "lengthScaledTPM" and the
# arithmetic is not obvious by inspection.
#
# The route under test is Snakemake -> 02_aggregate.R -> tx2gene from txdbmaker
# -> tximport. The independent route is tests/fixtures/tiny_project.py: the GTF
# parsed as text, and lengthScaledTPM evaluated from its definition, in another
# language. Both read the same quant.sf files, so a disagreement is a wiring
# defect in the stage, not a difference of opinion about the data.
#
# The fixture makes the correction observable: G001 holds its read depth and
# moves it from a 500 bp isoform to a 5000 bp one between the groups. Raw counts
# cannot see that; lengthScaledTPM must. For a single-transcript gene of
# constant effective length the correction cancels exactly, so without that
# switch the two branches would agree everywhere and neither would be tested.
#
# Exits 77 (skipped) without snakemake and the R stack.
#
# Run:  bash tests/check_counts.sh
set -uo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
FIXTURE="$REPO/tests/fixtures/tiny_project.py"
. "$REPO/tests/_gate.sh"
gate_runtime r

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
PROJ="$TMP/proj"

python3 "$FIXTURE" build "$PROJ" "$REPO" rnaseq \
    || { echo "FAIL: fixture build" >&2; exit 1; }

# Stage 2 alone: V04 is about aggregation and count semantics.
run_workflow "$PROJ" "$REPO" "$TMP/run.log" results/counts.tsv
rc=$?
if [ $rc -ne 0 ]; then
    echo "FAIL: stage 2 did not complete on the bulk branch (exit $rc)" >&2
    tail -40 "$TMP/run.log" >&2
    exit 1
fi

python3 "$FIXTURE" qualify "$PROJ" || exit 1
echo "check_counts.sh: bulk counts agree with the independent route"

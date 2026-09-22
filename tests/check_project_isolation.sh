#!/usr/bin/env bash
# Checks that two projects sharing one pipeline checkout stay out of each
# other's state (R17d), and that nothing writes into the checkout itself.
#
# What must not break:
#   - submit.sh runs in the project directory, not the repo, and refuses to run
#     inside the repo
#   - each project's gates/decisions.log records only its own run
#   - the repo's tracked files are untouched by a run
#   - a project without config/thresholds.yaml fails with a clear message rather
#     than silently using the repo's copy (the Snakefile would read the project's)
#   - setup.py refuses --project-dir pointing at the checkout
#
# snakemake is not needed: the launcher gets as far as exec'ing it, which is
# past every path decision this checks.
#
# Run:  bash tests/check_project_isolation.sh
set -uo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

fail() { echo "FAIL: $*" >&2; exit 1; }

mkproject() {  # $1 = dir, $2 = seq_type, $3 = n samples
    # A complete, preflight-clean project: submit.sh runs preflight before it
    # submits anything, so a half-filled config would be rejected there and this
    # check would never reach the isolation behaviour it is testing.
    mkdir -p "$1/config"
    cat > "$1/config/config.yaml" <<EOF
project: {name: $(basename "$1"), description: isolation fixture, output_dir: "results/"}
organism: {common_name: grape, scientific_name: Vitis vinifera, tax_id: 29760,
  kegg_code: vvi, orgdb_package: org.Vvinifera.eg.db}
reference: {accession: GCF_1.1, assembly_name: ASM1}
samples: {sheet: config/samples.tsv, seq_type: $2}
model: {fixed_effects: "~ treatment", random_effects: null, primary_factor: treatment}
downstream: {run_go: false, run_kegg: false, run_wgcna: false}
orgdb: {strategy: skip}
# No storage_budget_gb: the platform has no default cap (master plan 10.1).
hpc: {delete_fastq_after_quant: true, samples_in_flight: null}
EOF
    cp "$REPO/config/thresholds.yaml" "$1/config/thresholds.yaml"
    # Two levels with >=2 replicates each, so the design is estimable.
    { echo -e "sample_id\ttreatment"
      for i in $(seq 1 "$3"); do
          if [ $((i % 2)) -eq 0 ]; then lvl=treated; else lvl=control; fi
          echo -e "$(basename "$1")_S$i\t$lvl"
      done
    } > "$1/config/samples.tsv"
}

# Snapshot the checkout so we can prove a run did not write into it.
repo_state() { (cd "$REPO" && git status --porcelain 2>/dev/null | sort); }
BEFORE=$(repo_state)

# --- two projects, different seq types and sizes -----------------------
A="$TMP/study_a"; B="$TMP/study_b"
mkproject "$A" tagseq 6            # -> small tier
mkproject "$B" rnaseq_paired 40    # -> medium tier
BEFORE_A_CFG=$(cat "$A/config/config.yaml")

"$REPO/submit.sh" -d "$A" > "$TMP/a.log" 2>&1
"$REPO/submit.sh" -d "$B" > "$TMP/b.log" 2>&1

grep -q "project dir:   $A" "$TMP/a.log" || fail "A did not run in its own dir; see $TMP/a.log"
grep -q "project dir:   $B" "$TMP/b.log" || fail "B did not run in its own dir"
grep -q "tier:          small"  "$TMP/a.log" || fail "A: expected small tier"
grep -q "tier:          medium" "$TMP/b.log" || fail "B: expected medium tier"
# The workflow itself is read from the checkout, not copied per project.
grep -q "workflow:      $REPO" "$TMP/a.log" || fail "A did not read the workflow from the repo"

# --- each project logged only its own decision ------------------------
[ -f "$A/gates/decisions.log" ] || fail "A wrote no gates/decisions.log in its project dir"
[ -f "$B/gates/decisions.log" ] || fail "B wrote no gates/decisions.log in its project dir"
grep -q "STRATEGY: small"  "$A/gates/decisions.log" || fail "A's log lacks its own decision"
grep -q "STRATEGY: medium" "$B/gates/decisions.log" || fail "B's log lacks its own decision"
grep -q "medium" "$A/gates/decisions.log" && fail "B's decision leaked into A's log"
grep -q "small"  "$B/gates/decisions.log" && fail "A's decision leaked into B's log"

# --- A's config survived B's run untouched ---------------------------
[ "$(cat "$A/config/config.yaml")" = "$BEFORE_A_CFG" ] \
    || fail "running B modified A's config"

# --- the checkout is unchanged ---------------------------------------
AFTER=$(repo_state)
[ "$BEFORE" = "$AFTER" ] || {
    echo "--- repo changed during the runs:" >&2
    diff <(echo "$BEFORE") <(echo "$AFTER") >&2
    fail "a run wrote into the pipeline checkout"
}
[ ! -e "$REPO/gates/decisions.log.new" ] || fail "unexpected file in repo"

# --- refuses to run inside the repo ----------------------------------
out=$("$REPO/submit.sh" -d "$REPO" 2>&1); rc=$?
[ "$rc" -ne 0 ] || fail "running inside the repo was allowed"
echo "$out" | grep -q "Refusing to run inside the pipeline repo" \
    || fail "wrong error for repo-as-project: $out"

# --- missing thresholds is a clear failure, not a silent fallback -----
C="$TMP/study_c"; mkproject "$C" tagseq 4; rm "$C/config/thresholds.yaml"
out=$("$REPO/submit.sh" -d "$C" 2>&1); rc=$?
[ "$rc" -ne 0 ] || fail "missing thresholds.yaml did not stop the run"
echo "$out" | grep -q "Missing $C/config/thresholds.yaml" \
    || fail "wrong error for missing thresholds: $out"

# --- setup.py refuses the checkout as a project ----------------------
out=$(python3 "$REPO/scripts/setup.py" --project-dir "$REPO" 2>&1); rc=$?
[ "$rc" -ne 0 ] || fail "setup.py accepted the checkout as a project dir"
echo "$out" | grep -q "must not be the pipeline checkout" \
    || fail "wrong error from setup.py: $out"

# --- helper scripts are invoked by absolute path ----------------------
# The working directory is the project, so a repo-relative script path in a
# shell: block does not resolve. Guard against one creeping back in.
bad=$(grep -rn "bash workflow/\|python3\? workflow/\|Rscript scripts/\|source workflow/" \
        "$REPO"/workflow/rules/*.smk || true)
[ -z "$bad" ] || fail "repo-relative script path in a shell block:
$bad"
grep -q 'REPO_DIR = Path(config.get("repo_dir"' "$REPO/Snakefile" \
    || fail "Snakefile no longer derives REPO_DIR from config"
grep -q -- '--config repo_dir="\$REPO"' "$REPO/submit.sh" \
    || fail "submit.sh no longer passes repo_dir, so shell rules cannot find their scripts"

# --- a nonexistent project dir is rejected ---------------------------
out=$("$REPO/submit.sh" -d "$TMP/does_not_exist" 2>&1); rc=$?
[ "$rc" -ne 0 ] || fail "nonexistent project dir was accepted"
echo "$out" | grep -q "Project directory not found" || fail "wrong error: $out"

echo "check_project_isolation.sh: all assertions passed"

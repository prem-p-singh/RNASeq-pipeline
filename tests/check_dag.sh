#!/usr/bin/env bash
# Checks the two Snakemake behaviours R17d depends on, by building a real DAG.
#
# WORKING_PLAN 3.2 item 3 named these the unverified assumptions behind the
# per-project layout:
#
#   1. `include:` resolves against the repo's Snakefile when the working
#      directory is a project, not against the working directory.
#   2. `--config repo_dir=<checkout>` reaches the rule shell commands, so a
#      helper script is invoked by absolute path.
#
# Both were reasoned from documentation and neither had been executed. If the
# first were false the workflow could not load its own rule modules; if the
# second were false every shell rule would fail to find its script.
#
# Also covers R14 (`rule all` is the default target and the dry run reaches all
# requested terminal outputs) and confirms the content-keyed reference cache
# path from R03/R17d appears in the resolved commands.
#
# Needs snakemake, pandas and pyyaml. Skips cleanly when snakemake is absent,
# because that is the normal state until the runtime of item 3 exists.
#
# Run:  bash tests/check_dag.sh
set -uo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)

SMK="${SNAKEMAKE:-}"
if [ -z "$SMK" ]; then
    if command -v snakemake >/dev/null 2>&1; then SMK=$(command -v snakemake); fi
fi
if [ -z "$SMK" ] || ! "$SMK" --version >/dev/null 2>&1; then
    echo "check_dag.sh: SKIP (no snakemake; set SNAKEMAKE=/path/to/snakemake)"
    exit 0
fi

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
fail() { echo "FAIL: $*" >&2; exit 1; }

PROJ="$TMP/proj"
mkdir -p "$PROJ/config"
cat > "$PROJ/config/config.yaml" <<EOF
project: {name: dagcheck, description: dag probe, output_dir: "results/"}
organism: {common_name: grape, scientific_name: Vitis vinifera, tax_id: 29760,
  kegg_code: vvi, orgdb_package: org.Vvinifera.eg.db}
reference: {accession: GCF_1.1, assembly_name: ASM1, cache_dir: "$PROJ/refcache"}
samples: {sheet: config/samples.tsv, seq_type: tagseq}
model: {fixed_effects: "~ treatment", random_effects: null, primary_factor: treatment}
downstream: {run_go: false, run_kegg: false, run_wgcna: false}
orgdb: {strategy: skip}
hpc: {delete_fastq_after_quant: true, samples_in_flight: null}
EOF
cp "$REPO/config/thresholds.yaml" "$PROJ/config/thresholds.yaml"
{ echo -e "sample_id\tfastq_url\ttreatment"
  for i in 1 2 3; do echo -e "S$i\t/data/S$i.fq.gz\tcontrol"; done
  for i in 4 5 6; do echo -e "S$i\t/data/S$i.fq.gz\ttreated"; done
} > "$PROJ/config/samples.tsv"

cd "$PROJ"
OUT="$TMP/dag.txt"
"$SMK" --snakefile "$REPO/Snakefile" \
       --configfile config/config.yaml \
       --config repo_dir="$REPO" \
       --dry-run --printshellcmds > "$OUT" 2>&1 || {
    echo "--- snakemake output:" >&2; tail -30 "$OUT" >&2
    fail "dry run did not succeed"
}

# --- 1. include: resolved against the repo Snakefile ------------------
# If it had resolved against the working directory, no rule module would have
# loaded and none of these rules would exist.
for rule in fetch_transcriptome fetch_gtf salmon_index qc_quant_sample \
            aggregate_counts differential_expression enrichment wgcna qc_report; do
    grep -q "^rule ${rule}:\|^ *${rule} " "$OUT" \
        || fail "rule ${rule} absent; include: did not resolve from the repo Snakefile"
done

# --- 2. rule all is the default target and reaches the terminal outputs (R14)
grep -q "^rule all:" "$OUT" || fail "rule all was not selected as the default target"
for target in "results/counts.tsv" "results/de_done.flag" \
              "results/enrichment_done.flag" "results/wgcna_done.flag" \
              "results/qc_report/qc_report.done"; do
    grep -q "$target" "$OUT" || fail "terminal output $target not in the DAG"
done

# --- 3. repo_dir reached the shell commands ---------------------------
# Every helper must be invoked by an absolute path inside the checkout, and
# that path must exist. A relative path here is the R17d regression.
found=0
while read -r script; do
    [ -f "$script" ] || fail "shell command references a missing script: $script"
    case "$script" in
        "$REPO"/*) found=$((found + 1)) ;;
        *) fail "script path is outside the checkout: $script" ;;
    esac
done < <(grep -oE "${REPO}/workflow/scripts/[A-Za-z0-9_]+\.(sh|py)" "$OUT" | sort -u)
[ "$found" -ge 4 ] || fail "expected several repo-absolute helper paths, found $found"

# no helper may be invoked by a repo-relative path
grep -qE "(bash|python3?|Rscript|source) +workflow/scripts/" "$OUT" \
    && fail "a helper is invoked by a repo-relative path; it would not resolve"

# --- 4. outputs are project-relative, not repo-relative ---------------
grep -q "results/quant/S1/quant.sf" "$OUT" \
    || fail "per-sample outputs are not project-relative"
grep -q "${REPO}/results/" "$OUT" \
    && fail "an output path points into the checkout"

# --- 5. the content-keyed reference cache is in use (R03/R17d) --------
grep -qE "refcache/GCF_1\.1-[0-9a-f]{12}/" "$OUT" \
    || fail "reference cache path is not content-keyed"

# --- 6. R15 cleanup flag is passed through ----------------------------
grep -q -- "--delete-intermediates" "$OUT" \
    || fail "storage cleanup flag not passed to the per-sample rule"

echo "check_dag.sh: all assertions passed ($(grep -c '^rule ' "$OUT") rules resolved)"

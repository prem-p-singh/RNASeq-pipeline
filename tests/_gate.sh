# Shared skip gate for the checks that need a Snakemake runtime and the R stack.
#
# Sourced, not run. Sets SMK, or prints a reason and exits 77, which run_all.sh
# reports as SKIP. It lives here because duplicating the probe once already let
# a skipped check be counted as a pass, and this one is twelve R packages long.
#
#   . "$(dirname "$0")/_gate.sh"
#   gate_runtime          # needs snakemake only
#   gate_runtime r        # also needs the R stack for stages 2-5

gate_runtime() {
    SMK="${SNAKEMAKE:-}"
    [ -z "$SMK" ] && command -v snakemake >/dev/null 2>&1 && SMK=$(command -v snakemake)
    if [ -z "$SMK" ] || ! "$SMK" --version >/dev/null 2>&1; then
        echo "needs snakemake (set SNAKEMAKE=/path/to/snakemake)"
        exit 77
    fi
    [ "${1:-}" = "r" ] || return 0

    if ! command -v Rscript >/dev/null 2>&1; then
        echo "needs Rscript on PATH"
        exit 77
    fi
    # 03_de.R attaches variancePartition unconditionally, so it is needed even
    # on the limma-voom path these fixtures take.
    local missing
    missing=$(Rscript -e 'cat(paste(commandArgs(TRUE)[
        !sapply(commandArgs(TRUE), requireNamespace, quietly = TRUE)],
        collapse = " "))' \
        tximport GenomicFeatures edgeR limma variancePartition emmeans \
        clusterProfiler dplyr readr jsonlite tibble 2>/dev/null)
    if [ -n "$missing" ]; then
        echo "needs R packages: $missing"
        exit 77
    fi
}

# Runs the workflow in PROJ for the given targets. --scheduler greedy: the
# default ILP scheduler shells out to a CBC binary that pulp ships only for
# x86-64, so on arm64 it raises "Bad CPU type in executable". Greedy needs no
# solver and makes the run deterministic.
# SMK_EXTRA holds any further options, and is placed BEFORE --cores on purpose:
# --set-resources and --config take one-or-more values, so anything greedy must
# be followed by another flag rather than by the target list.
run_workflow() {  # $1 = project dir, $2 = repo, $3 = log, rest = targets
    local proj=$1 repo=$2 log=$3
    shift 3
    ( cd "$proj" && "$SMK" --snakefile "$repo/Snakefile" \
        --configfile config/config.yaml --config repo_dir="$repo" \
        ${SMK_EXTRA:-} --cores 2 --scheduler greedy "$@" ) > "$log" 2>&1
}

#!/usr/bin/env bash
# Runs every self-check and reports one line per check.
#
# Use --strict for release qualification: any skip fails the suite.
# In developer mode, checks whose interpreter or libraries are
# absent are reported as SKIP, with the reason, and do not count as passes:
# a skipped check is not a passed check.
#
# Use --n8n to run only the three orchestration checks.
# Run:  bash tests/run_all.sh
set -uo pipefail

strict=0
n8n_only=0
for option in "$@"; do
    case "$option" in
        --strict) strict=1 ;;
        --n8n) n8n_only=1 ;;
        *) echo "Usage: bash tests/run_all.sh [--strict] [--n8n]" >&2; exit 2 ;;
    esac
done

cd "$(dirname "$0")/.."

pass=0; fail=0; skip=0
failed_names=()

have_r_pkgs() {  # $@ = package names
    Rscript -e "q(status = if (all(sapply(commandArgs(TRUE), requireNamespace, quietly = TRUE))) 0 else 1)" \
        "$@" >/dev/null 2>&1
}
have_py_mods() {
    # importlib.util must be imported explicitly; `import importlib` alone does
    # not bind the submodule, which silently reported every module missing.
    python3 -c "import importlib.util, sys
sys.exit(0 if all(importlib.util.find_spec(m) for m in sys.argv[1:]) else 1)" \
        "$@" >/dev/null 2>&1
}

run() {  # $1 = label, $2 = command
    local label=$1 cmd=$2 out rc
    out=$(eval "$cmd" 2>&1); rc=$?
    if [ $rc -eq 77 ]; then
        # 77 is "skipped", the automake convention. A check that gates itself
        # on a missing tool reports it this way, because exiting 0 would be
        # indistinguishable from having run and passed.
        skip "$label" "$(echo "$out" | tail -1)"
    elif [ $rc -eq 0 ]; then
        printf '  PASS  %-28s\n' "$label"
        pass=$((pass + 1))
    else
        printf '  FAIL  %-28s (exit %d)\n' "$label" "$rc"
        echo "$out" | sed 's/^/          /' | tail -15
        fail=$((fail + 1))
        failed_names+=("$label")
    fi
}

skip() {
    printf '  SKIP  %-28s %s\n' "$1" "$2"
    skip=$((skip + 1))
}

echo "Self-checks"

if [ "$n8n_only" -eq 0 ]; then

if command -v Rscript >/dev/null 2>&1; then
    if have_r_pkgs edgeR limma emmeans; then
        run check_de.R "Rscript tests/check_de.R"
    else
        skip check_de.R "needs edgeR, limma, emmeans"
    fi
    if have_r_pkgs WGCNA edgeR dplyr readr jsonlite tibble; then
        run check_wgcna.R "Rscript tests/check_wgcna.R"
    else
        skip check_wgcna.R "needs WGCNA, edgeR, dplyr, readr, jsonlite, tibble"
    fi
else
    skip check_de.R    "no Rscript on PATH"
    skip check_wgcna.R "no Rscript on PATH"
fi

if have_py_mods yaml; then
    run check_setup_helpers.py "python3 tests/check_setup_helpers.py"
else
    skip check_setup_helpers.py "needs pyyaml"
fi

run check_qc_report.py       "python3 tests/check_qc_report.py"
run check_reference_lock.py  "python3 tests/check_reference_lock.py"
run check_config_resolve.py  "python3 tests/check_config_resolve.py"
run check_artifacts.py       "python3 tests/check_artifacts.py"
run check_decoys.py          "python3 tests/check_decoys.py"
run check_reference_cache.py "python3 tests/check_reference_cache.py"
run check_metadata.py        "python3 tests/check_metadata.py"
run check_read_units.py      "python3 tests/check_read_units.py"
if have_py_mods yaml; then
    run check_raw_qc.py "python3 tests/check_raw_qc.py"
    run check_raw_qc_real "python3 tests/check_raw_qc.py --real"
else
    skip check_raw_qc.py "needs pyyaml"
    skip check_raw_qc_real "needs pyyaml and fastp"
fi
run check_resources.py       "python3 tests/check_resources.py"
run check_recommend.py       "python3 tests/check_recommend.py"
if have_py_mods openpyxl pandas yaml; then
    run check_intake.py "python3 tests/check_intake.py"
else
    skip check_intake.py "needs openpyxl, pandas and pyyaml"
fi

if have_py_mods yaml && command -v Rscript >/dev/null 2>&1 \
   && have_r_pkgs jsonlite; then
    run check_preflight.py "python3 tests/check_preflight.py"
else
    skip check_preflight.py "needs pyyaml, Rscript and jsonlite"
fi

# Integration checks require the full release environment; absent dependencies
# are explicit skips and fail strict release qualification.
run check_dag.sh             "bash tests/check_dag.sh"
run check_stages.sh          "bash tests/check_stages.sh"
run check_counts.sh          "bash tests/check_counts.sh"
run check_network.sh         "bash tests/check_network.sh"
run check_end_to_end.sh      "bash tests/check_end_to_end.sh"
run check_objectives.sh      "bash tests/check_objectives.sh"
run check_dream.sh           "bash tests/check_dream.sh"
run check_deseq2.sh          "bash tests/check_deseq2.sh"
run check_edger.sh           "bash tests/check_edger.sh"
run check_designs.sh         "bash tests/check_designs.sh"
run check_screen.py          "python3 tests/check_screen.py"
run check_screen_real        "python3 tests/check_screen.py --real"
run check_screen_workflow    "bash tests/check_screen_workflow.sh"
run check_star.py            "python3 tests/check_star.py"
run check_star_workflow.sh   "bash tests/check_star_workflow.sh"
if command -v Rscript >/dev/null 2>&1; then
    run check_enrichment.R "Rscript tests/check_enrichment.R"
else
    skip check_enrichment.R "needs Rscript"
fi

run check_storage_policy.sh  "bash tests/check_storage_policy.sh"

# submit.sh runs preflight, and preflight shells out to Rscript for the
# estimability check, so this needs R even though the isolation behaviour
# it tests does not.
if have_py_mods yaml && command -v Rscript >/dev/null 2>&1 \
   && have_r_pkgs jsonlite; then
    run check_project_isolation.sh "bash tests/check_project_isolation.sh"
else
    skip check_project_isolation.sh "needs pyyaml, Rscript and jsonlite (submit.sh runs preflight)"
fi

fi

# Orchestration is a separate source line; never run its release smoke against
# this schema-4 parent. --gate reports missing deployment prerequisites as 77.
echo "n8n integration gate"
if command -v node >/dev/null 2>&1; then
    run check_n8n_nodes.js "node tests/check_n8n_nodes.js"
else
    skip check_n8n_nodes.js "needs Node.js"
fi
if command -v python3 >/dev/null 2>&1 \
   && python3 -c "import fcntl, sys; sys.exit(sys.version_info < (3, 9))" >/dev/null 2>&1; then
    run check_n8n.py "python3 tests/check_n8n.py"
    run check_n8n_release.py "python3 tests/check_n8n_release.py --gate"
else
    skip check_n8n.py "needs Python 3.9+ and POSIX fcntl"
    skip check_n8n_release.py "needs Python 3.9+, Linux x86_64, reviewed schema-7 source and locked modules"
fi

echo
echo "  $pass passed, $fail failed, $skip skipped"
if [ $fail -gt 0 ]; then
    echo "  failed: ${failed_names[*]}"
    exit 1
fi
if [ "$strict" -eq 1 ] && [ "$skip" -gt 0 ]; then
    echo "  RELEASE BLOCKED: required checks were skipped"
    exit 1
fi
exit 0

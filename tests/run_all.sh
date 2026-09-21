#!/usr/bin/env bash
# Runs every self-check and reports one line per check.
#
# Exits 0 only if all of them pass. Checks whose interpreter or libraries are
# absent are reported as SKIP, with the reason, and do not count as passes:
# a skipped check is not a passed check.
#
# Run:  bash tests/run_all.sh
set -uo pipefail

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
    if [ $rc -eq 0 ]; then
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

if have_py_mods yaml && command -v Rscript >/dev/null 2>&1 \
   && have_r_pkgs jsonlite; then
    run check_preflight.py "python3 tests/check_preflight.py"
else
    skip check_preflight.py "needs pyyaml, Rscript and jsonlite"
fi

run check_storage_policy.sh  "bash tests/check_storage_policy.sh"

if have_py_mods yaml; then
    run check_project_isolation.sh "bash tests/check_project_isolation.sh"
else
    skip check_project_isolation.sh "needs pyyaml (submit.sh reads config with it)"
fi

echo
echo "  $pass passed, $fail failed, $skip skipped"
if [ $fail -gt 0 ]; then
    echo "  failed: ${failed_names[*]}"
    exit 1
fi
exit 0

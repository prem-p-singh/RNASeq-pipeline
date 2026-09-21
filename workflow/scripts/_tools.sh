# Sourced by any rule that shells out to salmon/fastp.
#
# Prefers whatever the activated conda environment provides, because
# environment.yml pins those versions and an unpinned cluster module silently
# changes results. Cluster modules are a fallback, not the default, and the
# fallback announces itself.
#
# Usage:  source workflow/scripts/_tools.sh
#         ensure_tools fastp salmon

ensure_tools() {
    local missing=0 t
    for t in "$@"; do
        command -v "$t" >/dev/null 2>&1 || missing=1
    done
    [ "$missing" -eq 0 ] && return 0

    # The system module init is not safe under `set -e`/`set -u`: inside an
    # active conda env it runs `conda deactivate`, which is undefined in this
    # child shell, exits 127, and would abort the caller before any work.
    set +euo pipefail
    if [ -f /etc/profile.d/modules.sh ]; then
        source /etc/profile.d/modules.sh
        # Assumes module names match tool names, which holds on FARM for
        # salmon and fastp.
        module load "$@"
    fi
    set -euo pipefail

    for t in "$@"; do
        command -v "$t" >/dev/null 2>&1 || {
            echo "required tool not found on PATH or as a module: $t" >&2
            return 1
        }
    done
    echo "[$(date -Iseconds)] NOTE: fell back to cluster module(s) for: $*" \
         "(environment.yml versions are NOT active)" >&2
}

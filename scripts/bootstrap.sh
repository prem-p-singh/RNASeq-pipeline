#!/usr/bin/env bash
# Prepare one immutable release environment using Conda's explicit lock support.
# Conda itself may be supplied by the HPC site. No system packages are modified.
set -euo pipefail
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
MODULE=core
if [ "$#" -gt 0 ]; then
    [ "$#" -eq 2 ] && [ "$1" = --module ] && [ "$2" = star ] || {
        echo "Usage: bootstrap.sh [--module star]" >&2; exit 2;
    }
    MODULE=star
fi
LOCK="$REPO/environments/linux-64.explicit.txt"
[ "$MODULE" = core ] || LOCK="$REPO/environments/star-linux-64.explicit.txt"
VERIFY_PYTHON=$(command -v python3 || true)
[ "$(uname -sm)" = "Linux x86_64" ] || {
    echo "Use a Linux x86_64 host for this release environment." >&2; exit 1;
}
[ -s "$LOCK" ] || { echo "Release lock missing: $LOCK" >&2; exit 1; }
command -v conda >/dev/null || {
    echo "Conda is required. Activate your installation or load your site's Conda module." >&2; exit 1;
}
digest=$(sha256sum "$LOCK" | cut -d ' ' -f1)
prefix=${RNASEQ_ENV_PREFIX:-${XDG_CACHE_HOME:-$HOME/.cache}/rnaseq/environments/$digest}
if [ "$MODULE" = star ]; then
    prefix=${RNASEQ_STAR_ENV_PREFIX:-${XDG_CACHE_HOME:-$HOME/.cache}/rnaseq/environments/$digest}
    [ -n "$VERIFY_PYTHON" ] || { echo "Activate the core runtime before installing STAR tools" >&2; exit 1; }
fi
mkdir -p "$(dirname "$prefix")"
# The prefix must not be renamed: installed entry points embed its absolute path.
# The lock prevents concurrent writes; verification is the publication boundary.
exec 9>"$prefix.build.lock"
flock -w 3600 9 || { echo "Timed out waiting for environment build" >&2; exit 1; }
if [ ! -d "$prefix/conda-meta" ]; then
    conda create --yes --prefix "$prefix" --file "$LOCK" >&2
fi
set +u
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$prefix"
set -u
export PYTHONNOUSERSITE=1 R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null
export R_LIBS_USER="$prefix/lib/R/library" R_LIBS_SITE="$prefix/lib/R/library"
[ "$MODULE" != core ] || VERIFY_PYTHON="$prefix/bin/python3"
"$VERIFY_PYTHON" "$REPO/scripts/environment_check.py" --prefix "$prefix" --module "$MODULE" --lock "$LOCK" \
    --out "$prefix/environment_report.json" >&2
printf '%s\n' "$prefix"

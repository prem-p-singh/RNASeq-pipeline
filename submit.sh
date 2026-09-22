#!/usr/bin/env bash
# =============================================================================
# submit.sh — launcher with auto-tier selection
#
# Counts samples in config/samples.tsv, picks the SMALL / MEDIUM / LARGE
# SLURM profile, logs the decision, and kicks off Snakemake.
#
# Runs IN a project directory and reads the workflow FROM this repo, so several
# projects can share one checkout without overwriting each other. Everything the
# run writes (results/, reference/, metrics/, gates/, logs/, .snakemake/) is
# relative to the project directory; the Snakefile, rules and profiles are read
# from the repo.
#
#   cd ~/rnaseq_projects/my_study && ~/RNASeq_pipeline/submit.sh
#   ~/RNASeq_pipeline/submit.sh -d ~/rnaseq_projects/my_study
# =============================================================================
set -euo pipefail

# Where this script, and therefore the workflow, lives.
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

# --- Activate pipeline conda env (FARM) -------------------------------
# Auto-activate if available; no-op on systems without it.
# (temp disable -u; system profile scripts use unbound vars like MANPATH)
set +u
if [ -f /etc/profile.d/modules.sh ]; then
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load conda/base 2>/dev/null || true
fi
if command -v conda >/dev/null 2>&1; then
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate rnaseq-pipeline 2>/dev/null || \
        echo "NOTE: conda env 'rnaseq-pipeline' not found; running with current env"
fi
set -u

# --- Arguments ---------------------------------------------------------
# The first positional used to double as the config path *and* get forwarded
# to snakemake, where it was interpreted as a build target. Options are now
# parsed explicitly and only genuine pass-through arguments reach snakemake.
usage() {
    cat <<'USAGE'
Usage: submit.sh [options] [-- snakemake args...]

  -d, --directory PATH    project directory to run in (default: current dir).
                          All run output lands here; the workflow is read from
                          the repo this script lives in.
  -c, --configfile PATH   config file, relative to the project directory
                          (default: config/config.yaml)
  -p, --profile NAME      force a tier profile (small|medium|large),
                          overriding the sample-count choice
  -n, --dry-run           show the plan without running anything
  -h, --help              this message

Anything after `--` is passed straight to snakemake, e.g.
  submit.sh -- --forcerun differential_expression
USAGE
}

PROJDIR="."
CONFIG="config/config.yaml"
FORCE_PROFILE=""
DRY_RUN=""
PASSTHRU=()

while [ $# -gt 0 ]; do
    case "$1" in
        -d|--directory)  PROJDIR="$2";       shift 2 ;;
        -c|--configfile) CONFIG="$2";        shift 2 ;;
        -p|--profile)    FORCE_PROFILE="$2"; shift 2 ;;
        -n|--dry-run)    DRY_RUN="--dry-run"; shift ;;
        -h|--help)       usage; exit 0 ;;
        --)              shift; PASSTHRU=("$@"); break ;;
        -*)              echo "Unknown option: $1" >&2; usage >&2; exit 1 ;;
        *)               echo "Unexpected argument: $1" >&2
                         echo "(config files go after -c; snakemake targets after --)" >&2
                         exit 1 ;;
    esac
done

# Everything below reads and writes relative paths, so move into the project
# first. This is what keeps two projects out of each other's state.
if [ ! -d "$PROJDIR" ]; then
    echo "Project directory not found: $PROJDIR" >&2
    exit 1
fi
cd "$PROJDIR"
PROJDIR=$(pwd)

if [ "$PROJDIR" = "$REPO" ]; then
    echo "Refusing to run inside the pipeline repo itself ($REPO)." >&2
    echo "Give the project its own directory, e.g." >&2
    echo "  $REPO/submit.sh -d ~/rnaseq_projects/<name>" >&2
    exit 1
fi

if [ ! -f "$CONFIG" ]; then
    echo "Config not found: $PROJDIR/$CONFIG" >&2
    echo "Run scripts/new_project.sh, or setup.py, to create it." >&2
    exit 1
fi

# Thresholds must live in the project, not the repo: the Snakefile opens this
# same relative path, so a fallback here would let the launcher and the workflow
# disagree about tier sizes and QC gates. setup.py copies the defaults in.
THRESH="config/thresholds.yaml"
if [ ! -f "$THRESH" ]; then
    echo "Missing $PROJDIR/$THRESH" >&2
    echo "Copy the defaults in:  cp $REPO/config/thresholds.yaml $PROJDIR/config/" >&2
    exit 1
fi

# --- Preflight ---------------------------------------------------------
# Validates config keys, the sample sheet, the assay/design/backend combination
# and model estimability before a single FASTQ is fetched. An unestimable design
# used to surface only in Stage 3, i.e. after the whole cohort was quantified.
echo "--- preflight ---"
if ! python3 "$REPO/scripts/preflight.py" -d "$PROJDIR" -c "$CONFIG"; then
    echo
    echo "Preflight found blocking problems; nothing was submitted." >&2
    echo "Full issue table: $PROJDIR/gates/preflight_issues.tsv" >&2
    exit 1
fi
echo "-----------------"

SHEET=$(python3 -c "import yaml; print(yaml.safe_load(open('$CONFIG'))['samples']['sheet'])")

# --- Count samples (strip comments + header) --------------------------
N=$(grep -v '^#' "$SHEET" | awk 'NR>1 && NF>0' | wc -l | tr -d ' ')

# --- Tier profile: cluster settings only ------------------------------
# Master plan 10.7: "Replace small/medium/large sample tiers as the primary
# scheduler logic with per-rule resources and dataset-aware concurrency.
# Friendly size labels may remain in reports." The tier therefore still selects
# the SLURM profile (account, partition, qos), but no longer decides how much
# work runs at once: that comes from the storage plan below.
read TIER_SMALL TIER_MEDIUM SEQ_TYPE <<< $(python3 - <<PYCFG
import yaml
t = yaml.safe_load(open("$THRESH"))["tiers"]
c = yaml.safe_load(open("$CONFIG"))
print(t["small_max_samples"], t["medium_max_samples"], c["samples"]["seq_type"])
PYCFG
)

if   [ "$N" -le "$TIER_SMALL"  ]; then TIER=small
elif [ "$N" -le "$TIER_MEDIUM" ]; then TIER=medium
else                                   TIER=large
fi

if [ -n "$FORCE_PROFILE" ]; then
    TIER="$FORCE_PROFILE"
    echo "NOTE: tier overridden to '$TIER' by --profile"
fi
PROFILE="$REPO/profiles/$TIER"
if [ ! -d "$PROFILE" ]; then
    echo "No such profile: profiles/$TIER (expected small|medium|large)" >&2
    exit 1
fi

# --- Storage plan: measured, per filesystem, no platform cap ----------
# Replaces hpc.storage_budget_gb and the per-assay FASTQ guesses (R15, U07,
# master plan 10.7). plan_resources.py exits non-zero when the run cannot fit,
# after trying a lower concurrency first (RS08/RS09).
PLAN_JSON="$PROJDIR/gates/resource_plan.json"
if ! python3 "$REPO/scripts/plan_resources.py" \
        --project-dir "$PROJDIR" --configfile "$CONFIG" \
        --profile "$PROFILE" --out "$PLAN_JSON"; then
    echo
    echo "Storage planning blocked this run; nothing was submitted." >&2
    exit 1
fi
MAX_CONC=$(python3 -c "import json;print(json.load(open('$PLAN_JSON'))['planned_concurrency'])")
PROFILE_JOBS=$(python3 -c "import yaml;print(yaml.safe_load(open('$PROFILE/config.yaml')).get('jobs',1))")

JOBS_FLAG=""
if [ "$MAX_CONC" -lt "$PROFILE_JOBS" ]; then
    JOBS_FLAG="--jobs $MAX_CONC"
    echo "NOTE: --jobs capped at $MAX_CONC (profile allows $PROFILE_JOBS) by the storage plan"
fi

# --- Log decision -----------------------------------------------------
mkdir -p gates
{
    echo "[$(date -Iseconds)] STRATEGY: $TIER"
    echo "    N=$N  seq_type=$SEQ_TYPE  max_concurrent=$MAX_CONC"
    echo "    storage plan: $PLAN_JSON"
    echo "    profile=$PROFILE"
} >> gates/decisions.log

echo "==================================================="
echo " Submitting RNA-Seq pipeline"
echo "   project dir:   $PROJDIR"
echo "   workflow:      $REPO"
echo "   samples:       $N"
echo "   seq_type:      $SEQ_TYPE"
echo "   tier:          $TIER"
echo "   profile:       $PROFILE"
echo "   max concurrent: $MAX_CONC"
echo "==================================================="

# --- Launch -----------------------------------------------------------
# Already chdir'd into the project, so snakemake's working directory is the
# project and every relative path in the workflow lands there. -s points at the
# repo's Snakefile; its `include:` paths resolve relative to that Snakefile.
#
# ${PASSTHRU[@]+"${PASSTHRU[@]}"} expands to nothing when the array is empty,
# which plain "${PASSTHRU[@]}" does not do safely under `set -u` on bash 3.2
# (the version macOS ships).
exec snakemake \
    --snakefile "$REPO/Snakefile" \
    --profile "$PROFILE" \
    --configfile "$CONFIG" \
    --config repo_dir="$REPO" \
    $JOBS_FLAG \
    $DRY_RUN \
    ${PASSTHRU[@]+"${PASSTHRU[@]}"}

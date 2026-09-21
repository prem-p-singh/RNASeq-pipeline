#!/usr/bin/env bash
# =============================================================================
# submit.sh — launcher with auto-tier selection
#
# Counts samples in config/samples.tsv, picks the SMALL / MEDIUM / LARGE
# SLURM profile, logs the decision, and kicks off Snakemake.
# =============================================================================
set -euo pipefail

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
Usage: ./submit.sh [options] [-- snakemake args...]

  -c, --configfile PATH   config file (default: config/config.yaml)
  -p, --profile NAME      force a tier profile (small|medium|large),
                          overriding the sample-count choice
  -n, --dry-run           show the plan without running anything
  -h, --help              this message

Anything after `--` is passed straight to snakemake, e.g.
  ./submit.sh -- --forcerun differential_expression
USAGE
}

CONFIG="config/config.yaml"
FORCE_PROFILE=""
DRY_RUN=""
PASSTHRU=()

while [ $# -gt 0 ]; do
    case "$1" in
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

if [ ! -f "$CONFIG" ]; then
    echo "Config not found: $CONFIG" >&2
    exit 1
fi

THRESH="config/thresholds.yaml"
SHEET=$(python3 -c "import yaml; print(yaml.safe_load(open('$CONFIG'))['samples']['sheet'])")

# --- Count samples (strip comments + header) --------------------------
N=$(grep -v '^#' "$SHEET" | awk 'NR>1 && NF>0' | wc -l | tr -d ' ')

# --- Read tier thresholds + fastq size estimate -----------------------
read TIER_SMALL TIER_MEDIUM RESERVED TAGSEQ_GB RNA_SE_GB RNA_PE_GB SEQ_TYPE <<< $(python3 - <<PY
import yaml
t = yaml.safe_load(open("$THRESH"))["tiers"]
c = yaml.safe_load(open("$CONFIG"))
print(
    t["small_max_samples"], t["medium_max_samples"], t["reserved_gb"],
    t["fastq_size_estimate_gb"]["tagseq"],
    t["fastq_size_estimate_gb"]["rnaseq_single"],
    t["fastq_size_estimate_gb"]["rnaseq_paired"],
    c["samples"]["seq_type"],
)
PY
)

case "$SEQ_TYPE" in
    tagseq)         FASTQ_GB=$TAGSEQ_GB ;;
    rnaseq_single)  FASTQ_GB=$RNA_SE_GB ;;
    rnaseq_paired)  FASTQ_GB=$RNA_PE_GB ;;
    *) echo "Unknown seq_type: $SEQ_TYPE" >&2; exit 1 ;;
esac

# --- Pick tier --------------------------------------------------------
if   [ "$N" -le "$TIER_SMALL"  ]; then TIER=small
elif [ "$N" -le "$TIER_MEDIUM" ]; then TIER=medium
else                                   TIER=large
fi

if [ -n "$FORCE_PROFILE" ]; then
    TIER="$FORCE_PROFILE"
    echo "NOTE: tier overridden to '$TIER' by --profile"
fi
PROFILE="profiles/$TIER"
if [ ! -d "$PROFILE" ]; then
    echo "No such profile: $PROFILE (expected small|medium|large)" >&2
    exit 1
fi

# --- Storage budget: cap concurrency, refuse impossible runs ----------
# Assigned to a variable first, NOT piped straight into `read`: under `set -e`
# a failed command substitution only aborts the script when it is the whole
# assignment, so `read ... <<< $(python3 ...)` would print REFUSED and carry on.
BUDGET=$(python3 - <<'PY' "$CONFIG" "$PROFILE" "$RESERVED" "$FASTQ_GB" "$N" "$SEQ_TYPE"
import sys, yaml

cfg_path, profile, reserved, fastq_gb, n, seq_type = sys.argv[1:7]
reserved, fastq_gb, n = float(reserved), float(fastq_gb), int(n)

hpc = yaml.safe_load(open(cfg_path))["hpc"]
budget   = float(hpc["storage_budget_gb"])
delete   = bool(hpc.get("delete_fastq_after_quant", True))
inflight = hpc.get("samples_in_flight")
prof_jobs = int(yaml.safe_load(open(f"{profile}/config.yaml")).get("jobs", 1))

avail = budget - reserved

# Peak on-disk working set for ONE sample. Paired-end holds raw and trimmed
# reads at the same time (01_qc_quant.sh writes fastp output to disk so
# --gcBias can re-read it), so its peak is twice the raw-pair estimate.
# Single-end streams fastp into salmon, so only the raw file is on disk.
peak = fastq_gb * (2 if seq_type == "rnaseq_paired" else 1)

def refuse(msg):
    sys.exit(f"REFUSED: {msg}")

if avail <= 0:
    refuse(f"storage_budget_gb={budget:g} leaves nothing after {reserved:g} GB reserved "
           f"for the index and outputs")
if peak > avail:
    refuse(f"one {seq_type} sample peaks at ~{peak:g} GB but only {avail:g} GB is usable "
           f"(budget {budget:g} minus {reserved:g} reserved); raise storage_budget_gb")
if not delete:
    total = peak * n
    if total > avail:
        refuse(f"delete_fastq_after_quant is false, so all {n} samples stay on disk "
               f"(~{total:g} GB) but only {avail:g} GB is usable; set it to true")

conc = max(1, int(avail // peak))
if inflight:
    conc = min(conc, int(inflight))
print(f"{budget:g}", min(conc, prof_jobs), prof_jobs)
PY
)
read BUDGET_GB MAX_CONC PROFILE_JOBS <<< "$BUDGET"

JOBS_FLAG=""
if [ "$MAX_CONC" -lt "$PROFILE_JOBS" ]; then
    JOBS_FLAG="--jobs $MAX_CONC"
    echo "NOTE: --jobs capped at $MAX_CONC (profile allows $PROFILE_JOBS) to stay under ${BUDGET_GB}GB"
fi

# --- Log decision -----------------------------------------------------
mkdir -p gates
{
    echo "[$(date -Iseconds)] STRATEGY: $TIER"
    echo "    N=$N  seq_type=$SEQ_TYPE  avg_fastq=${FASTQ_GB}GB"
    echo "    budget=${BUDGET_GB}GB  reserved=${RESERVED}GB  max_concurrent=$MAX_CONC"
    echo "    profile=$PROFILE"
} >> gates/decisions.log

echo "==================================================="
echo " Submitting RNA-Seq pipeline"
echo "   samples:       $N"
echo "   seq_type:      $SEQ_TYPE"
echo "   tier:          $TIER"
echo "   profile:       $PROFILE"
echo "   max concurrent: $MAX_CONC"
echo "==================================================="

# --- Launch -----------------------------------------------------------
# ${PASSTHRU[@]+"${PASSTHRU[@]}"} expands to nothing when the array is empty,
# which plain "${PASSTHRU[@]}" does not do safely under `set -u` on bash 3.2
# (the version macOS ships).
exec snakemake \
    --profile "$PROFILE" \
    --configfile "$CONFIG" \
    $JOBS_FLAG \
    $DRY_RUN \
    ${PASSTHRU[@]+"${PASSTHRU[@]}"}

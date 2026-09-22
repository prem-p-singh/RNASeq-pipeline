#!/bin/bash -l
# Submit from the checkout: sbatch scripts/validate_farm.sbatch
# Site defaults below are FARM-specific. Override with sbatch options elsewhere.
#SBATCH --job-name=rnaseq-release
#SBATCH --account=publicgrp
#SBATCH --partition=low
#SBATCH --qos=publicgrp-low-qos
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=03:00:00
set -eo pipefail
module load conda/latest
TASK_REPO="$SLURM_SUBMIT_DIR"
TASK_WORK=$(mktemp -d /tmp/rnaseq-release.XXXXXX)
TASK_EVIDENCE="$TASK_REPO/docs/evidence/farm-$SLURM_JOB_ID"
mkdir -p "$TASK_EVIDENCE"
printf 'Worker: %s\nJob: %s\nUTC: %s\n' "$(hostname)" "$SLURM_JOB_ID" "$(date -u +%FT%TZ)" > "$TASK_EVIDENCE/worker.txt"
export TMPDIR="$TASK_WORK" CONDA_PKGS_DIRS="$TASK_WORK/pkgs"
export RNASEQ_ENV_PREFIX="$TASK_WORK/runtime"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
prefix=$(bash "$TASK_REPO/scripts/bootstrap.sh")
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$prefix"
export PYTHONNOUSERSITE=1 R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null
export R_LIBS_USER="$prefix/lib/R/library" R_LIBS_SITE="$prefix/lib/R/library"
cp "$prefix/environment_report.json" "$TASK_EVIDENCE/"
cd "$TASK_REPO"
set -o pipefail
bash tests/run_all.sh --strict 2>&1 | tee "$TASK_EVIDENCE/tests.log"
RNASEQ_PUBLIC_PROJECT="$TASK_WORK/public" bash tests/check_public.sh 2>&1 | tee "$TASK_EVIDENCE/public.log"
# Reusing the verified environment must also work without network access.
CONDA_OFFLINE=true bash scripts/bootstrap.sh > "$TASK_EVIDENCE/offline-reuse.txt"
printf 'Temporary runtime and public project retained at %s\n' "$TASK_WORK"

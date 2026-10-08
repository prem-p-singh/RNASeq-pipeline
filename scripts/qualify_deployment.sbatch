#!/bin/bash -l
# Deployment qualification beyond validate_slurm.sbatch, which runs every rule
# inside one allocation. Submit from the checkout with Conda available:
#   RNASEQ_QUALIFY_ROOT=/shared/dir sbatch [site options] scripts/qualify_deployment.sbatch [distributed|container]
#
# distributed (default): this job is only the controller. submit.sh hands each
#   rule to SLURM workers, one run is interrupted and resumed, and two projects
#   build the same reference cache entry at the same time. RNASEQ_QUALIFY_ROOT
#   must be visible to every worker. Set RNASEQ_ENV_PREFIX to reuse a verified
#   runtime. Workers take the account and partition from the profile, or from
#   SBATCH_ACCOUNT / SBATCH_PARTITION when the site has no default.
# container: fresh installation and the Salmon route in a stock Miniforge image
#   under Apptainer, without host modules, home directory or environment.
#SBATCH --job-name=rnaseq-deploy
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
set -eo pipefail
MODE=${1:-distributed}
REPO="$SLURM_SUBMIT_DIR"
WORK=$(mktemp -d "${RNASEQ_QUALIFY_ROOT:?set RNASEQ_QUALIFY_ROOT to a shared directory}/rnaseq-deploy.XXXXXX")
EVIDENCE="$REPO/docs/evidence/deploy-$MODE-$SLURM_JOB_ID"
mkdir -p "$EVIDENCE"
printf 'Mode: %s\nController: %s\nJob: %s\nUTC: %s\nWork: %s\n' "$MODE" "$(hostname)" "$SLURM_JOB_ID" "$(date -u +%FT%TZ)" "$WORK" | tee "$EVIDENCE/controller.txt"
fail() { echo "FAIL: $*" >&2; exit 1; }
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

if [ "$MODE" = container ]; then
    IMAGE=${RNASEQ_QUALIFY_IMAGE:-docker://condaforge/miniforge3:26.7.2-0}
    export APPTAINER_CACHEDIR="$WORK/apptainer-cache" APPTAINER_TMPDIR="$WORK/apptainer-tmp"
    mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR" "$WORK/home" "$WORK/tmp"
    apptainer exec --containall --no-home --bind "$WORK:/work" --bind "$REPO:/src" "$IMAGE" bash -eo pipefail -s <<'INSIDE' 2>&1 | tee "$EVIDENCE/container.log"
export HOME=/work/home TMPDIR=/work/tmp XDG_CACHE_HOME=/work/cache CONDA_PKGS_DIRS=/work/pkgs
export RNASEQ_ENV_PREFIX=/work/runtime OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
grep PRETTY_NAME /etc/os-release; conda --version
cd /src
prefix=$(bash scripts/bootstrap.sh)
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$prefix"
export PYTHONNOUSERSITE=1 R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null
export R_LIBS_USER="$prefix/lib/R/library" R_LIBS_SITE="$prefix/lib/R/library"
cp "$prefix/environment_report.json" /work/environment_report.json
bash tests/check_end_to_end.sh
RNASEQ_PUBLIC_PROJECT=/work/public bash tests/check_public.sh
CONDA_OFFLINE=true bash scripts/bootstrap.sh
echo "CONTAINER PASS"
INSIDE
    cp "$WORK/environment_report.json" "$EVIDENCE/"
    exit 0
fi
[ "$MODE" = distributed ] || fail "mode must be distributed or container"

prefix=$(bash "$REPO/scripts/bootstrap.sh")
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$prefix"
export RNASEQ_ENV_PREFIX="$prefix" PYTHONNOUSERSITE=1 R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null
export R_LIBS_USER="$prefix/lib/R/library" R_LIBS_SITE="$prefix/lib/R/library"
cp "$prefix/environment_report.json" "$EVIDENCE/"
cd "$REPO"

new_project() { python3 tests/fixtures/public_project.py "$WORK/$1" "$REPO"; }
# Background commands of a script ignore SIGINT; restore it so the launcher can be interrupted.
launch() { python3 -c 'import os,signal,sys; signal.signal(signal.SIGINT, signal.SIG_DFL); os.execvp(sys.argv[1], sys.argv[1:])' \
    bash "$REPO/submit.sh" -d "$WORK/$1" --executor slurm; }
verify() { RNASEQ_PUBLIC_PROJECT="$WORK/$1" bash tests/check_public.sh --verify; }
workers() {  # $1 = launcher log: every rule job must have completed on a SLURM worker
    local ids
    ids=$(grep -o 'SLURM jobid [0-9]*' "$1" | awk '{print $3}' | sort -u | paste -sd, -)
    [ -n "$ids" ] || fail "no SLURM worker jobs in $1"
    sacct -X -n -P -j "$ids" -o JobID,JobName,State,NodeList,Elapsed | tee "$1.workers"
    ! grep -q "^$SLURM_JOB_ID|" "$1.workers" || fail "a rule ran in the controller allocation"
}

echo "== 1. launcher to SLURM workers"
new_project distributed
launch distributed > "$EVIDENCE/distributed.log" 2>&1 || { tail -60 "$EVIDENCE/distributed.log"; fail "distributed run"; }
verify distributed
workers "$EVIDENCE/distributed.log"
! grep -v '|COMPLETED|' "$EVIDENCE/distributed.log.workers" || fail "a worker job did not complete"

echo "== 2. interruption and resume"
new_project restart
launch restart > "$EVIDENCE/restart-1.log" 2>&1 &
pid=$!
until ls "$WORK"/restart/results/quant/*/metrics.json >/dev/null 2>&1; do
    kill -0 $pid 2>/dev/null || { tail -60 "$EVIDENCE/restart-1.log"; fail "run ended before it could be interrupted"; }
    sleep 5
done
kill -INT $pid
if wait $pid; then fail "interrupted run reported success"; fi
kept=$(ls "$WORK"/restart/results/quant/*/metrics.json)
before=$(stat -c '%n %Y' $kept)
echo "$before" | tee "$EVIDENCE/restart-kept.txt"
[ "$(echo "$kept" | wc -l)" -lt 6 ] || fail "all six samples finished before the interruption"
launch restart > "$EVIDENCE/restart-2.log" 2>&1 || { tail -60 "$EVIDENCE/restart-2.log"; fail "resume"; }
verify restart
[ "$(stat -c '%n %Y' $kept)" = "$before" ] || fail "resume recomputed samples that had finished"
workers "$EVIDENCE/restart-2.log"

echo "== 3. two projects, one reference cache entry, started together"
new_project shared-a
new_project shared-b
python3 - "$WORK" <<'PY'
import sys, yaml
from pathlib import Path
work = Path(sys.argv[1])
source = yaml.safe_load((work / 'shared-a/config/config.yaml').read_text())['reference']['transcriptome_fasta_url']
for name in ('shared-a', 'shared-b'):
    path = work / name / 'config/config.yaml'
    cfg = yaml.safe_load(path.read_text())
    # The cache key covers the source URL, so both projects name one copy.
    cfg['reference']['transcriptome_fasta_url'] = source
    cfg['reference']['cache_dir'] = str(work / 'reference-cache')
    path.write_text(yaml.safe_dump(cfg, sort_keys=False))
PY
launch shared-a > "$EVIDENCE/shared-a.log" 2>&1 & pid_a=$!
launch shared-b > "$EVIDENCE/shared-b.log" 2>&1 & pid_b=$!
wait $pid_a || { tail -60 "$EVIDENCE/shared-a.log"; fail "shared-a"; }
wait $pid_b || { tail -60 "$EVIDENCE/shared-b.log"; fail "shared-b"; }
verify shared-a
verify shared-b
workers "$EVIDENCE/shared-a.log"
workers "$EVIDENCE/shared-b.log"
ls -la "$WORK/reference-cache" | tee "$EVIDENCE/reference-cache.txt"
[ "$(ls "$WORK/reference-cache" | wc -l)" -eq 1 ] || fail "expected one cache entry and no leftover build lock"
python3 - "$WORK" <<'PY' | tee "$EVIDENCE/shared-counts.txt"
import csv, sys
from pathlib import Path
def totals(project):
    rows = list(csv.reader((Path(sys.argv[1]) / project / 'results/counts.tsv').open(), delimiter='\t'))
    return rows[0], [r[0] for r in rows[1:]], [sum(float(r[i]) for r in rows[1:]) for i in range(1, len(rows[0]))]
(ha, ga, ta), (hb, gb, tb) = totals('shared-a'), totals('shared-b')
assert ha == hb and ga == gb, 'sample or gene sets differ'
# Tolerance set before the run: Salmon totals within 1% between the two projects.
worst = max(abs(x - y) / max(x, 1) for x, y in zip(ta, tb))
assert worst < 0.01, worst
print(f'{len(ga)} genes, {len(ta)} samples, largest per-sample total difference {worst:.2e}')
PY
echo "DEPLOYMENT PASS"

#!/usr/bin/env bash
# =============================================================================
# start_new.sh <project_name> <metadata_file> [fastq_dir]
#
# MAC-SIDE wrapper. Uploads metadata (and optional FASTQs) into FARM's
# new_project_inbox, then SSHes over to kick off scripts/new_project.sh
# interactively.
#
# Example:
#   scripts/start_new.sh grape_2021 ~/Downloads/grape_meta.xlsx
#   scripts/start_new.sh cancer_study ~/data/meta.csv ~/data/fastq/
# =============================================================================
set -euo pipefail

if [ $# -lt 2 ]; then
    cat <<EOF
Usage: scripts/start_new.sh <project_name> <metadata_file> [fastq_dir]

  project_name     alphanumeric + underscore (no spaces)
  metadata_file    local path to .xlsx/.csv/.tsv
  fastq_dir        optional local dir with FASTQs; if omitted, you'll be
                   prompted interactively for an S3/HTTPS URL or remote path.

After uploading, you'll drop into an interactive SSH session on FARM
running new_project.sh — pick organism, confirm column guesses, done.
EOF
    exit 1
fi

PROJECT="$1"
META="$2"
FASTQ_DIR="${3:-}"

if [[ ! "$PROJECT" =~ ^[A-Za-z0-9_]+$ ]]; then
    echo "Project name must be alphanumeric/underscore only"
    exit 1
fi
if [ ! -f "$META" ]; then
    echo "Metadata file not found: $META"
    exit 1
fi

INBOX="~/new_project_inbox/$PROJECT"

# --- create inbox, upload metadata -----------------------------------
echo "Creating inbox on FARM: $INBOX"
ssh farm "mkdir -p $INBOX"
echo "Uploading metadata ($(basename "$META"))..."
scp -q "$META" "farm:$INBOX/"

# --- upload FASTQs if given ------------------------------------------
if [ -n "$FASTQ_DIR" ]; then
    if [ ! -d "$FASTQ_DIR" ]; then
        echo "FASTQ dir not found: $FASTQ_DIR"
        exit 1
    fi
    n=$(find "$FASTQ_DIR" -name "*.fastq.gz" 2>/dev/null | wc -l | tr -d ' ')
    if [ "$n" -eq 0 ]; then
        echo "No *.fastq.gz under $FASTQ_DIR — skipping"
    else
        echo "Uploading $n FASTQ files (may take a while)..."
        rsync -avz --progress --include="*.fastq.gz" --include="*/" --exclude="*" \
              "$FASTQ_DIR/" "farm:$INBOX/"
    fi
fi

# --- kick off new_project.sh on FARM ---------------------------------
echo
echo "Starting interactive setup on FARM..."
ssh -t farm "cd ~/RNASeq_pipeline && scripts/new_project.sh $PROJECT"

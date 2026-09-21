#!/usr/bin/env bash
# =============================================================================
# new_project.sh <project_name>
#
# Run on FARM. Expects ~/new_project_inbox/<project_name>/ to contain:
#     - a metadata file (.xlsx / .csv / .tsv)  [required]
#     - FASTQ files (*.fastq.gz) OR fastq_urls.txt  [optional — can point elsewhere]
#
# Auto-detects sample-ID and factor columns, prompts for an organism preset and
# a few minimal choices, writes setup_inputs.yaml, then runs the full setup +
# submit pipeline.
# =============================================================================
set -euo pipefail

# --- args -------------------------------------------------------------
if [ $# -lt 1 ]; then
    echo "Usage: scripts/new_project.sh <project_name>"
    exit 1
fi
PROJECT="$1"

# Sanitize project name (letters, digits, underscore)
if [[ ! "$PROJECT" =~ ^[A-Za-z0-9_]+$ ]]; then
    echo "Project name must be alphanumeric/underscore only: got '$PROJECT'"
    exit 1
fi

INBOX=~/new_project_inbox/"$PROJECT"
REPO=$(cd "$(dirname "$0")/.." && pwd)
PRESETS="$REPO/scripts/presets"
# Projects live outside the checkout, one directory each, so nothing a run
# writes can reach another project or the pipeline source.
PROJECTS_ROOT=${RNASEQ_PROJECTS_ROOT:-~/rnaseq_projects}
PROJECTS_ROOT=${PROJECTS_ROOT/#\~/$HOME}

if [ ! -d "$INBOX" ]; then
    echo "Inbox missing: $INBOX"
    echo "Create it and drop your metadata + FASTQs there first, then re-run."
    exit 1
fi

# --- activate env (temp disable -u; these system scripts use unbound vars) -
set +u
if [ -f /etc/profile.d/modules.sh ]; then
    source /etc/profile.d/modules.sh
    module load conda/base 2>/dev/null || true
fi
if command -v conda >/dev/null 2>&1; then
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate rnaseq-pipeline 2>/dev/null || true
fi
set -u

echo
echo "╭──────────────────────────────────────────────────────────────╮"
echo "│  New RNA-Seq project: $PROJECT"
echo "│  Inbox:   $INBOX"
echo "│  Project: $PROJECTS_ROOT/$PROJECT"
echo "│  Workflow: $REPO"
echo "╰──────────────────────────────────────────────────────────────╯"

# ====================================================================== 1. metadata
META=$(ls "$INBOX"/*.xlsx "$INBOX"/*.csv "$INBOX"/*.tsv 2>/dev/null | head -1 || true)
if [ -z "$META" ]; then
    echo "ERROR: no metadata file (.xlsx/.csv/.tsv) found in $INBOX"
    exit 1
fi
echo
echo "[1/6] Metadata file:  $META"
INFO_JSON=$(python3 "$REPO/scripts/detect_metadata.py" "$META")
N_SAMPLES=$(echo "$INFO_JSON" | python3 -c 'import sys,json; print(json.load(sys.stdin)["n_samples"])')
echo "      $N_SAMPLES samples detected"

# ====================================================================== 2. FASTQs
echo
echo "[2/6] FASTQ source:"
if [ -f "$INBOX/fastq_urls.txt" ]; then
    FASTQ_SOURCE=$(head -1 "$INBOX/fastq_urls.txt" | sed 's|/[^/]*$|/|')   # infer base dir
    N_FASTQS=$(grep -c . "$INBOX/fastq_urls.txt")
    echo "      using URL list:  $INBOX/fastq_urls.txt  ($N_FASTQS URLs)"
    FASTQ_SOURCE_FINAL="$INBOX/fastq_urls.txt"
elif compgen -G "$INBOX/*.fastq.gz" > /dev/null; then
    N_FASTQS=$(ls -1 "$INBOX"/*.fastq.gz | wc -l)
    echo "      using local files in inbox:  $N_FASTQS *.fastq.gz"
    FASTQ_SOURCE_FINAL="$INBOX"
else
    echo "      (no FASTQs in inbox and no fastq_urls.txt)"
    read -r -p "      Path or S3/HTTPS URL where FASTQs live: " FASTQ_SOURCE_FINAL
fi

# ====================================================================== 3. organism preset
echo
echo "[3/6] Pick organism:"
PS3="      Choose: "
shopt -s nullglob
PRESET_FILES=("$PRESETS"/*.yaml)
shopt -u nullglob
PRESET_NAMES=()
for f in "${PRESET_FILES[@]}"; do
    PRESET_NAMES+=("$(basename "$f" .yaml)")
done
PRESET_NAMES+=("other (manual taxID)")

select choice in "${PRESET_NAMES[@]}"; do
    [ -n "$choice" ] && break
done

if [ "$choice" = "other (manual taxID)" ]; then
    read -r -p "      NCBI taxID: " TAX_ID
    read -r -p "      Genus species (e.g. 'Zea mays'): " SCI_NAME
    ORGDB="org.$(echo "$SCI_NAME" | awk '{print substr($1,1,1) $2}').eg.db"
    KEGG_CODE=""
    echo "      → orgdb guess: $ORGDB (will be built if not on Bioconductor)"
    echo "      → kegg code not set — edit setup_inputs.yaml if needed"
else
    PRESET_FILE="$PRESETS/$choice.yaml"
    # Parse YAML cleanly in Python (avoids brittle sed)
    eval "$(python3 - "$PRESET_FILE" <<'PY'
import sys, yaml
d = yaml.safe_load(open(sys.argv[1]))
for k, v in d.items():
    if k in ("tax_id","scientific_name","orgdb_package","kegg_code","common_name"):
        print(f"PRESET_{k.upper()}={v!r}")
PY
)"
    TAX_ID="$PRESET_TAX_ID"
    SCI_NAME="$PRESET_SCIENTIFIC_NAME"
    ORGDB="$PRESET_ORGDB_PACKAGE"
    KEGG_CODE="$PRESET_KEGG_CODE"
    echo "      → $SCI_NAME (taxID $TAX_ID, OrgDb $ORGDB, KEGG $KEGG_CODE)"
fi

# ====================================================================== 4. sample_id column
echo
echo "[4/6] Sample-ID column:"
SAMPLE_ID_GUESS=$(echo "$INFO_JSON" | python3 -c '
import sys, json
d = json.load(sys.stdin)
c = d["sample_id_candidates"]
print(c[0] if c else d["columns"][0])')
echo "      all columns: $(echo "$INFO_JSON" | python3 -c 'import sys,json; print(", ".join(json.load(sys.stdin)["columns"]))')"
read -r -p "      Which column holds sample IDs? [$SAMPLE_ID_GUESS]: " SAMPLE_ID_COL
SAMPLE_ID_COL=${SAMPLE_ID_COL:-$SAMPLE_ID_GUESS}

# ====================================================================== 5. primary factor
echo
echo "[5/6] Primary factor (for DE contrasts):"
echo "      Candidates (columns with 2–10 levels):"
# Write JSON to a temp file so python can read it from argv
_tmp=$(mktemp)
echo "$INFO_JSON" > "$_tmp"
python3 - "$_tmp" <<'PY'
import sys, json
d = json.load(open(sys.argv[1]))
for i, c in enumerate(d["factor_candidates"], 1):
    levels = ", ".join(c["levels"][:6])
    print(f"        {i}. {c['name']}  ({c['n_unique']} levels: {levels})")
PY
rm -f "$_tmp"
FACTOR_GUESS=$(echo "$INFO_JSON" | python3 -c '
import sys, json
d = json.load(sys.stdin)
f = d["factor_candidates"]
print(f[0]["name"] if f else "")')
read -r -p "      Which column is your primary factor? [$FACTOR_GUESS]: " PRIMARY
PRIMARY=${PRIMARY:-$FACTOR_GUESS}

# ====================================================================== 6. seq type
echo
echo "[6/6] Sequencing type:"
PS3="      Choose: "
select SEQ_TYPE in tagseq rnaseq_single rnaseq_paired; do
    [ -n "$SEQ_TYPE" ] && break
done

# ====================================================================== write setup_inputs.yaml
# Each project gets its own directory. The previous version wrote one shared
# $REPO/setup_inputs.yaml and ran everything inside the checkout, so a second
# project silently overwrote the first one's config, samples and results.
PROJDIR="$PROJECTS_ROOT/$PROJECT"
if [ -e "$PROJDIR" ] && [ -n "$(ls -A "$PROJDIR" 2>/dev/null)" ]; then
    echo
    echo "Project directory already exists and is not empty:"
    echo "  $PROJDIR"
    read -r -p "Reuse it (existing config and results are kept)? [y/N] " reuse
    if [[ ! "$reuse" =~ ^[Yy]$ ]]; then
        echo "Stopped. Pick another project name, or move that directory aside."
        exit 1
    fi
fi
mkdir -p "$PROJDIR/config"

SETUP_YAML="$PROJDIR/setup_inputs.yaml"
cat > "$SETUP_YAML" <<EOF
# Auto-generated by scripts/new_project.sh on $(date -Iseconds)
project_name: "$PROJECT"
tax_id: $TAX_ID
seq_type: "$SEQ_TYPE"
fastq_source: "$FASTQ_SOURCE_FINAL"
metadata_file: "$META"
sample_id_column: "$SAMPLE_ID_COL"
fastq_pattern: null
primary_factor: "$PRIMARY"
model_fixed_effects: "~ $PRIMARY"
model_random_effects: null
storage_budget_gb: 20
EOF

echo
echo "────────────────────────────────────────────────────────────────"
echo " Wrote $SETUP_YAML"
echo "────────────────────────────────────────────────────────────────"
cat "$SETUP_YAML"
echo "────────────────────────────────────────────────────────────────"

# ====================================================================== confirm + run
echo
read -r -p "Proceed with setup + submit? [y/N] " go
if [[ ! "$go" =~ ^[Yy]$ ]]; then
    echo "Stopped. Review the YAML, edit if needed, then run manually:"
    echo "  python3 $REPO/scripts/setup.py --project-dir $PROJDIR $SETUP_YAML"
    echo "  $REPO/submit.sh -d $PROJDIR"
    exit 0
fi

echo
echo "╭── Running scripts/setup.py ─────────────────────────────────╮"
python3 "$REPO/scripts/setup.py" --project-dir "$PROJDIR" "$SETUP_YAML"

echo
echo "╭── Launching pipeline ───────────────────────────────────────╮"
"$REPO/submit.sh" -d "$PROJDIR"

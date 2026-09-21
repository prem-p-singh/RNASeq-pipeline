#!/usr/bin/env bash
# =============================================================================
# 01_qc_quant.sh — Per-sample QC + trim + Salmon pseudoalignment.
#
# Merge of ritu-farm's real-data-tested implementation (paired-end, non-streaming
# for --gcBias safety, NextSeq/NovaSeq poly-G trim) + handbook-driven additions
# (Salmon --seqBias/--posBias, strandedness verification gate).
#
# For paired-end, fastp writes cleaned reads to disk and salmon reads those
# files. We do NOT stream through pipes, because salmon's --gcBias re-reads the
# input and a pipe can only be read once (streaming + --gcBias deadlocks).
#
# Storage policy: locally-staged FASTQs are read in place, never copied and
# never deleted. Files this script creates itself (downloaded reads, trimmed
# reads) are collected in $owned and removed after quant.sf is verified, when
# --delete-intermediates true. Cleanup iterates that list, so a user-supplied
# path cannot be reached by it even by accident.
# =============================================================================
set -euo pipefail

# --- Tools --------------------------------------------------------------
# Prefer the activated environment (environment.yml pins salmon and fastp);
# cluster modules are only a fallback and announce themselves when used.
source "$(dirname "$0")/_tools.sh"
ensure_tools fastp salmon

# Recorded per sample so results carry the software that produced them.
# tr -d '"' keeps the value safe to embed in metrics.json.
salmon_version=$(salmon --version 2>&1 | head -1 | tr -d '"')
fastp_version=$(fastp --version 2>&1 | head -1 | tr -d '"')

# --- Parse args ---------------------------------------------------------
sample=""
url=""
url2=""
seq_type=""
index=""
outdir=""
threads=4
min_map_rate=0.60
expected_libtype=""              # optional; if set, warn on mismatch with Salmon auto-detect
delete_intermediates="false"     # default off: an unpassed flag must never delete data

while [[ $# -gt 0 ]]; do
    case $1 in
        --sample)                sample=$2;               shift 2 ;;
        --url)                   url=$2;                  shift 2 ;;
        --url2)                  url2=$2;                 shift 2 ;;
        --seq-type)              seq_type=$2;             shift 2 ;;
        --index)                 index=$2;                shift 2 ;;
        --outdir)                outdir=$2;               shift 2 ;;
        --threads)               threads=$2;              shift 2 ;;
        --min-map-rate)          min_map_rate=$2;         shift 2 ;;
        --expected-libtype)      expected_libtype=$2;     shift 2 ;;
        --delete-intermediates)  delete_intermediates=$2; shift 2 ;;
        *) echo "Unknown arg: $1" >&2; exit 1 ;;
    esac
done

mkdir -p "$outdir"

# Files this script creates. Cleanup only ever iterates this list.
owned=()

# --- Fetch helper (remote sources only) ---------------------------------
fetch() {
    local src=$1 dst=$2
    case $src in
        s3://*)     aws s3 cp "$src" "$dst" ;;
        http*|ftp*) curl -sSL -o "$dst" "$src" ;;
        *) echo "Unrecognized remote URL scheme: $src" >&2; exit 1 ;;
    esac
}

# --- Locate / pull R1 ---------------------------------------------------
# Absolute local paths are read in place (no copy, no deletion).
# Remote sources are downloaded into the work dir and kept.
if [[ "$url" == /* ]]; then
    fastq_local="$url"
    echo "[$(date -Iseconds)] R1 in place: $fastq_local"
else
    fastq_local="$outdir/${sample}.fastq.gz"
    echo "[$(date -Iseconds)] Downloading R1: $url"
    fetch "$url" "$fastq_local"
    owned+=("$fastq_local")
fi

# --- Locate / pull R2 (paired-end only) --------------------------------
if [[ "$seq_type" == "rnaseq_paired" ]]; then
    if [[ "$url2" == /* ]]; then
        fastq_local_r2="$url2"
        echo "[$(date -Iseconds)] R2 in place: $fastq_local_r2"
    else
        fastq_local_r2="$outdir/${sample}_R2.fastq.gz"
        echo "[$(date -Iseconds)] Downloading R2: $url2"
        fetch "$url2" "$fastq_local_r2"
        owned+=("$fastq_local_r2")
    fi
fi

# --- Bias-correction flags (handbook section 7.2.2) ---------------------
# --gcBias, --seqBias, --posBias are designed for standard RNA-seq where
# fragments span the transcript. They REMOVE technical bias there.
#
# For TAGseq, the library IS intentionally 3'-end-biased; applying positional
# bias correction would fight the assay's biology. Use only --gcBias (mild,
# safe across protocols).
if [[ "$seq_type" == "tagseq" ]]; then
    salmon_bias_flags="--gcBias"
else
    salmon_bias_flags="--gcBias --seqBias --posBias"
fi

# --- fastp clean -> salmon quant ---------------------------------------
echo "[$(date -Iseconds)] Running fastp + salmon for $sample ($seq_type)"
fastp_report="$outdir/fastp.json"

if [[ "$seq_type" == "rnaseq_paired" ]]; then
    # Paired-end: fastp cleans R1+R2 to disk, then salmon reads those files.
    #   --detect_adapter_for_pe -> adapter detection for paired reads
    #   --trim_poly_g           -> removes poly-G tails NextSeq/NovaSeq produce
    trim_r1="$outdir/${sample}_R1.trim.fastq.gz"
    trim_r2="$outdir/${sample}_R2.trim.fastq.gz"
    owned+=("$trim_r1" "$trim_r2")
    fastp \
        -i "$fastq_local" \
        -I "$fastq_local_r2" \
        -o "$trim_r1" \
        -O "$trim_r2" \
        --detect_adapter_for_pe \
        --trim_poly_g \
        --json "$fastp_report" \
        --html "$outdir/fastp.html" \
        --thread "$threads" \
        2> "$outdir/fastp.log"

    salmon quant \
        -i "$index" \
        -l A \
        -1 "$trim_r1" \
        -2 "$trim_r2" \
        -p "$threads" \
        --validateMappings \
        $salmon_bias_flags \
        -o "$outdir" \
        2> "$outdir/salmon.log"
else
    # Single-end (e.g. TAGseq, rnaseq_single): fastp streams trimmed reads
    # into salmon via stdout. --gcBias is single-pass on SE so streaming is OK.
    fastp \
        -i "$fastq_local" \
        --stdout \
        --trim_poly_g \
        --json "$fastp_report" \
        --html "$outdir/fastp.html" \
        --thread "$threads" \
        2> "$outdir/fastp.log" \
    | salmon quant \
        -i "$index" \
        -l A \
        -r - \
        -p "$threads" \
        --validateMappings \
        $salmon_bias_flags \
        -o "$outdir" \
        2> "$outdir/salmon.log"
fi

# --- Verify salmon produced usable output -------------------------------
# quant.sf is the primary output and is checked here rather than left to
# Snakemake, because the cleanup below must not run on a failed quantification.
if [[ ! -s "$outdir/quant.sf" ]]; then
    echo "ERROR: $outdir/quant.sf missing or empty — salmon failed; see salmon.log" >&2
    exit 3
fi

# Salmon writes aux_info/meta_info.json with mapping stats.
meta="$outdir/aux_info/meta_info.json"
if [[ ! -f "$meta" ]]; then
    echo "ERROR: $meta missing — salmon failed" >&2
    exit 3
fi

map_rate=$(python3 -c "
import json
with open('$meta') as f: d = json.load(f)
# Salmon reports percent_mapped as a percentage (0-100)
print(d.get('percent_mapped', 0) / 100.0)
")

num_processed=$(python3 -c "
import json
with open('$meta') as f: d = json.load(f)
print(d.get('num_processed', 0))
")

# Salmon's auto-detected library type (e.g. ISR for paired stranded-reverse,
# SR for single-end stranded-reverse, U for unstranded). Handbook section 4.1.3
# flags this as a routine source of silent ~50% count errors.
detected_libtype=$(python3 -c "
import json
with open('$meta') as f: d = json.load(f)
lt = d.get('library_types', [])
print(lt[0] if lt else '')
")

# --- Gate: mapping rate -----------------------------------------------
flagged="false"
if awk "BEGIN {exit !($map_rate < $min_map_rate)}"; then
    flagged="true"
    echo "[$(date -Iseconds)] WARN: mapping rate $map_rate < $min_map_rate — flagged"
fi

# --- Gate: strandedness mismatch --------------------------------------
libtype_mismatch="false"
if [[ -n "$expected_libtype" && -n "$detected_libtype" ]]; then
    exp_upper=$(echo "$expected_libtype" | tr '[:lower:]' '[:upper:]')
    det_upper=$(echo "$detected_libtype" | tr '[:lower:]' '[:upper:]')
    if [[ "$exp_upper" != "$det_upper" ]]; then
        libtype_mismatch="true"
        echo "[$(date -Iseconds)] WARN: Salmon detected libtype=$detected_libtype but config expected $expected_libtype"
    fi
fi

# --- Write metrics.json -----------------------------------------------
cat > "$outdir/metrics.json" <<EOF
{
  "sample": "$sample",
  "seq_type": "$seq_type",
  "num_reads_processed": $num_processed,
  "mapping_rate": $map_rate,
  "min_map_rate_threshold": $min_map_rate,
  "flagged_low_mapping": $flagged,
  "detected_libtype": "$detected_libtype",
  "expected_libtype": "$expected_libtype",
  "libtype_mismatch": $libtype_mismatch,
  "salmon_bias_flags": "$salmon_bias_flags",
  "salmon_version": "$salmon_version",
  "fastp_version": "$fastp_version",
  "fastp_report": "$fastp_report"
}
EOF

# --- Log decisions if flagged -----------------------------------------
mkdir -p gates
if [[ "$flagged" == "true" ]]; then
    echo "[$(date -Iseconds)] SAMPLE_QC: $sample flagged (map_rate=$map_rate < $min_map_rate)" \
        >> gates/decisions.log
fi
if [[ "$libtype_mismatch" == "true" ]]; then
    echo "[$(date -Iseconds)] LIBTYPE_MISMATCH: $sample detected=$detected_libtype expected=$expected_libtype" \
        >> gates/decisions.log
fi

# --- Release owned intermediates ---------------------------------------
# Reached only once quant.sf verified above, so a failed run keeps its
# downloads (which may be the only local copy) for a rerun. A low mapping rate
# is NOT a reason to keep them: quant.sf is still valid input for Stage 2, which
# decides inclusion (see sample_disposition in 02_aggregate.R).
if [[ "$delete_intermediates" == "true" ]]; then
    for f in ${owned[@]+"${owned[@]}"}; do
        rm -f "$f"
        echo "[$(date -Iseconds)] released intermediate: $f"
    done
fi

echo "[$(date -Iseconds)] Done: $sample (map_rate=$map_rate, reads=$num_processed)"

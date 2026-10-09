#!/usr/bin/env bash
# =============================================================================
# 01_qc_quant.sh — Per-sample QC + trim + Salmon pseudoalignment.
#
# Merge of Ritu's real-data-tested implementation (paired-end, non-streaming
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

# --- Parse args ---------------------------------------------------------
stage="all"
processing_json="{}"
sample=""
url=""
url2=""
reads_json=""
seq_type=""
index=""
outdir=""
threads=4
min_map_rate=0.60
expected_libtype=""              # optional; if set, warn on mismatch with Salmon auto-detect
delete_intermediates="false"     # default off: an unpassed flag must never delete data

while [[ $# -gt 0 ]]; do
    case $1 in
        --stage)                 stage=$2; shift 2 ;;
        --processing-json)       processing_json=$2; shift 2 ;;
        --sample)                sample=$2;               shift 2 ;;
        --url)                   url=$2;                  shift 2 ;;
        --url2)                  url2=$2;                 shift 2 ;;
        --reads-json)            reads_json=$2;           shift 2 ;;
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

case "$stage" in all|preprocess|quant) ;; *) echo "Unknown stage: $stage" >&2; exit 2 ;; esac
source "$(dirname "$0")/_tools.sh"
if [[ "$stage" == quant ]]; then ensure_tools salmon; else ensure_tools fastp; fi
if [[ "$stage" == all ]]; then ensure_tools salmon; fi
mkdir -p "$outdir"

# Files this script creates. Cleanup only ever iterates this list.
owned=()

# --- Fetch helper (remote sources only) ---------------------------------
fetch() {
    local src=$1 dst=$2
    local partial
    partial=$(mktemp "${dst}.partial.XXXXXX")
    case $src in
        s3://*)     aws s3 cp "$src" "$partial" ;;
        http*|ftp*) curl --fail --retry 3 -sSL -o "$partial" "$src" ;;
        *) echo "Unrecognized remote URL scheme: $src" >&2; exit 1 ;;
    esac
    gzip -t "$partial"
    mv "$partial" "$dst"
}

if [[ "$stage" != quant ]]; then
# --- Locate / pull R1 ---------------------------------------------------
if [[ -n "$reads_json" ]]; then
    python3 "$(dirname "$0")/../../scripts/prepare_reads.py" \
        --manifest-json "$reads_json" --outdir "$outdir"
    url="$outdir/prepared_R1.fastq.gz"
    owned+=("$url")
    if [[ "$seq_type" == "rnaseq_paired" ]]; then
        url2="$outdir/prepared_R2.fastq.gz"
        owned+=("$url2")
    fi
fi
# Absolute local paths are read in place (no copy, no deletion).
# Remote sources are downloaded into the work dir and kept.
if [[ "$url" != *://* || "$url" == file://* ]]; then
    fastq_local="${url#file://}"
    echo "[$(date -Iseconds)] R1 in place: $fastq_local"
else
    fastq_local="$outdir/${sample}.fastq.gz"
    echo "[$(date -Iseconds)] Downloading R1: $url"
    fetch "$url" "$fastq_local"
    owned+=("$fastq_local")
fi

# --- Locate / pull R2 (paired-end only) --------------------------------
if [[ "$seq_type" == "rnaseq_paired" ]]; then
    if [[ "$url2" != *://* || "$url2" == file://* ]]; then
        fastq_local_r2="${url2#file://}"
        echo "[$(date -Iseconds)] R2 in place: $fastq_local_r2"
    else
        fastq_local_r2="$outdir/${sample}_R2.fastq.gz"
        echo "[$(date -Iseconds)] Downloading R2: $url2"
        fetch "$url2" "$fastq_local_r2"
        owned+=("$fastq_local_r2")
    fi
fi

python3 "$(dirname "$0")/../../scripts/preprocessing.py" \
    --r1 "$fastq_local" --r2 "${fastq_local_r2:-}" --outdir "$outdir" \
    --sample "$sample" --threads "$threads" --policy-json "$processing_json"
fi
fastp_version=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["fastp_version"])' "$outdir/preprocessing.json")
if [[ "$stage" == preprocess ]]; then
    # Downloads/merged copies have no downstream consumer. Trimmed reads do.
    if [[ "$delete_intermediates" == true ]]; then
        for f in ${owned[@]+"${owned[@]}"}; do rm -f "$f"; done
    fi
    exit 0
fi
salmon_version=$(salmon --version 2>&1 | head -1 | tr -d '"')

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

# Preprocessing is an independent producer; Salmon only consumes its trims.
if [[ "$seq_type" == rnaseq_paired ]]; then
    trim_r1="$outdir/${sample}_R1.trim.fastq.gz"
    trim_r2="$outdir/${sample}_R2.trim.fastq.gz"
    owned+=("$trim_r1" "$trim_r2")
    salmon quant -i "$index" -l A -1 "$trim_r1" -2 "$trim_r2" \
        -p "$threads" --validateMappings $salmon_bias_flags -o "$outdir" 2> "$outdir/salmon.log"
else
    trim_se="$outdir/${sample}.trim.fastq.gz"
    owned+=("$trim_se")
    salmon quant -i "$index" -l A -r "$trim_se" -p "$threads" \
        --validateMappings $salmon_bias_flags -o "$outdir" 2> "$outdir/salmon.log"
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

# Parse values as data, never as shell/Python source. Missing metrics fail closed.
qc_values=$(python3 - "$meta" "$outdir" "$sample" "$seq_type" "$min_map_rate" \
    "$expected_libtype" "$salmon_bias_flags" "$salmon_version" "$fastp_version" "${fastq_local:-}" "${fastq_local_r2:-}" <<'PYQC'
import csv, hashlib, json, math, sys
from pathlib import Path
meta, outdir, sample, seq_type, min_map, expected, bias, salmon_ver, fastp_ver = sys.argv[1:10]
out = Path(outdir)
with open(meta) as f:
    data = json.load(f)
rate = float(data["percent_mapped"]) / 100
processed = int(data["num_processed"])
if not math.isfinite(rate) or not 0 <= rate <= 1 or processed <= 0:
    raise SystemExit("Invalid Salmon mapping metrics")
with open(out / "quant.sf") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))
if not rows or not {"Name", "NumReads", "TPM", "EffectiveLength"}.issubset(rows[0]):
    raise SystemExit("Salmon quant.sf is not a valid quantification table")
for row in rows:
    if any(not math.isfinite(float(row[k])) or float(row[k]) < 0
           for k in ("NumReads", "TPM", "EffectiveLength")):
        raise SystemExit("Nonfinite or negative quantification value")
libtypes = data.get("library_types", [])
detected = libtypes[0] if libtypes else ""
flagged = rate < float(min_map)
mismatch = bool(expected and detected and expected.upper() != detected.upper())
metrics = dict(sample=sample, seq_type=seq_type, num_reads_processed=processed,
    mapping_rate=rate, min_map_rate_threshold=float(min_map), flagged_low_mapping=flagged,
    detected_libtype=detected, expected_libtype=expected, libtype_mismatch=mismatch,
    salmon_bias_flags=bias, salmon_version=salmon_ver, fastp_version=fastp_ver,
    fastp_report=str(out / "fastp.json"))
preprocessing = json.loads((out / "preprocessing.json").read_text())
metrics["input_files"] = preprocessing["input_files"]
metrics["preprocessing"] = {k: preprocessing[k] for k in
    ("policy", "fragments_before", "fragments_after", "fraction_lost", "warnings")}
(out / "metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
print(rate, processed, str(flagged).lower(), str(mismatch).lower(), detected)
PYQC
)
read -r map_rate num_processed flagged libtype_mismatch detected_libtype <<< "$qc_values"

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

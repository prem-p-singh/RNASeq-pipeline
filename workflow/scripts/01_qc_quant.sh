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
# Use the verified release environment; missing tools stop the stage.
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

# --- Locate / pull R1 ---------------------------------------------------
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
    # Keep a seekable input for every Salmon bias/optimization pass.
    trim_se="$outdir/${sample}.trim.fastq.gz"
    owned+=("$trim_se")
    fastp \
        -i "$fastq_local" \
        -o "$trim_se" \
        --trim_poly_g \
        --json "$fastp_report" \
        --html "$outdir/fastp.html" \
        --thread "$threads" \
        2> "$outdir/fastp.log"
    salmon quant \
        -i "$index" \
        -l A \
        -r "$trim_se" \
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

# Parse values as data, never as shell/Python source. Missing metrics fail closed.
qc_values=$(python3 - "$meta" "$outdir" "$sample" "$seq_type" "$min_map_rate" \
    "$expected_libtype" "$salmon_bias_flags" "$salmon_version" "$fastp_version" "$fastq_local" "${fastq_local_r2:-}" <<'PYQC'
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
metrics["input_files"] = []
for raw in sys.argv[10:]:
    if raw:
        with open(raw, "rb") as f:
            digest = hashlib.file_digest(f, "sha256").hexdigest()
        metrics["input_files"].append(dict(path=raw, bytes=Path(raw).stat().st_size, sha256=digest))
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

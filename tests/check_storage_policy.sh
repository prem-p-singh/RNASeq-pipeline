#!/usr/bin/env bash
# Checks the storage policy in workflow/scripts/01_qc_quant.sh:
#   - files the script creates itself (downloaded reads, trimmed reads) are
#     released once quant.sf exists
#   - a user-supplied FASTQ read in place is NEVER deleted
#   - nothing is released when quant.sf is missing, so a rerun still has input
#   - nothing is released when --delete-intermediates false
#
# fastp, salmon and curl are stubbed, so this runs anywhere in a second.
# Run:  bash tests/check_storage_policy.sh
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
SCRIPT="$ROOT/workflow/scripts/01_qc_quant.sh"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

# --- Stubs ------------------------------------------------------------
mkdir -p "$TMP/bin"

cat > "$TMP/bin/fastp" <<'STUB'
#!/usr/bin/env bash
out=""; out2=""; json=""; html=""; to_stdout=0
while [ $# -gt 0 ]; do
  case $1 in
    -o) out=$2; shift 2 ;;
    -O) out2=$2; shift 2 ;;
    --json) json=$2; shift 2 ;;
    --html) html=$2; shift 2 ;;
    --stdout) to_stdout=1; shift ;;
    --version) echo "fastp 0.0.0-stub"; exit 0 ;;
    *) shift ;;
  esac
done
[ -n "$json" ] && echo '{}' > "$json"
[ -n "$html" ] && echo '<html></html>' > "$html"
[ -n "$out" ]  && printf 'trimmed' | gzip > "$out"
[ -n "$out2" ] && printf 'trimmed' | gzip > "$out2"
[ "$to_stdout" -eq 1 ] && printf 'trimmed'
exit 0
STUB

# STUB_NO_QUANT=1 makes salmon "fail" by writing meta_info.json but no quant.sf,
# which is how the script learns the result is unusable.
cat > "$TMP/bin/salmon" <<'STUB'
#!/usr/bin/env bash
o=""
while [ $# -gt 0 ]; do
  case $1 in
    -o) o=$2; shift 2 ;;
    --version) echo "salmon 0.0.0-stub"; exit 0 ;;
    *) shift ;;
  esac
done
mkdir -p "$o/aux_info"
printf '{"percent_mapped": 91.0, "num_processed": 1000, "library_types": ["SR"]}\n' \
    > "$o/aux_info/meta_info.json"
if [ "${STUB_NO_QUANT:-0}" != "1" ]; then
    printf 'Name\tLength\tEffectiveLength\tTPM\tNumReads\nt1\t100\t80\t5.0\t10\n' > "$o/quant.sf"
fi
exit 0
STUB

cat > "$TMP/bin/curl" <<'STUB'
#!/usr/bin/env bash
dst=""
while [ $# -gt 0 ]; do
  case $1 in
    -o) dst=$2; shift 2 ;;
    *) shift ;;
  esac
done
printf 'downloaded' | gzip > "$dst"
exit 0
STUB

chmod +x "$TMP/bin"/*
export PATH="$TMP/bin:$PATH"

fail() { echo "FAIL: $*" >&2; exit 1; }
gone()   { [ ! -e "$1" ] || fail "$2: expected released, still present: $1"; }
exists() { [   -e "$1" ] || fail "$2: expected kept, is missing: $1"; }

# --- Case 1: downloaded single-end, delete=true ------------------------
# Single-end streams fastp into salmon, so the only owned file is the download.
cd "$TMP"
d="$TMP/c1"
bash "$SCRIPT" --sample S1 --url "http://example/S1.fastq.gz" --url2 "" \
    --seq-type tagseq --index "$TMP/idx" --outdir "$d" --threads 1 \
    --min-map-rate 0.6 --expected-libtype "" --delete-intermediates true \
    > "$TMP/c1.log" 2>&1 || fail "case 1: script exited nonzero (see $TMP/c1.log)"
exists "$d/quant.sf"            "case 1"
gone   "$d/S1.fastq.gz"         "case 1 download"
grep -q "released intermediate" "$TMP/c1.log" || fail "case 1: no release logged"

# --- Case 2: in-place paired-end, delete=true -------------------------
# The two source FASTQs belong to the user and must survive; the two trimmed
# files are ours and must go.
src1="$TMP/src/S2_R1.fastq.gz"; src2="$TMP/src/S2_R2.fastq.gz"
mkdir -p "$TMP/src"; printf 'raw' | gzip > "$src1"; printf 'raw' | gzip > "$src2"
d="$TMP/c2"
bash "$SCRIPT" --sample S2 --url "$src1" --url2 "$src2" \
    --seq-type rnaseq_paired --index "$TMP/idx" --outdir "$d" --threads 1 \
    --min-map-rate 0.6 --expected-libtype "" --delete-intermediates true \
    > "$TMP/c2.log" 2>&1 || fail "case 2: script exited nonzero (see $TMP/c2.log)"
exists "$src1" "case 2 USER SOURCE R1"
exists "$src2" "case 2 USER SOURCE R2"
gone   "$d/S2_R1.trim.fastq.gz" "case 2 trimmed R1"
gone   "$d/S2_R2.trim.fastq.gz" "case 2 trimmed R2"

# --- Case 3: delete=false keeps everything ----------------------------
d="$TMP/c3"
bash "$SCRIPT" --sample S3 --url "http://example/S3.fastq.gz" --url2 "" \
    --seq-type tagseq --index "$TMP/idx" --outdir "$d" --threads 1 \
    --min-map-rate 0.6 --expected-libtype "" --delete-intermediates false \
    > "$TMP/c3.log" 2>&1 || fail "case 3: script exited nonzero"
exists "$d/S3.fastq.gz" "case 3 download with delete=false"

# --- Case 4: no quant.sf fails the sample and keeps the download ------
d="$TMP/c4"
rc=0
STUB_NO_QUANT=1 bash "$SCRIPT" --sample S4 --url "http://example/S4.fastq.gz" --url2 "" \
    --seq-type tagseq --index "$TMP/idx" --outdir "$d" --threads 1 \
    --min-map-rate 0.6 --expected-libtype "" --delete-intermediates true \
    > "$TMP/c4.log" 2>&1 || rc=$?
[ "$rc" -eq 3 ] || fail "case 4: expected exit 3 on missing quant.sf, got $rc"
exists "$d/S4.fastq.gz" "case 4 download when quant.sf absent"
grep -q "quant.sf missing or empty" "$TMP/c4.log" || fail "case 4: no failure reason logged"

echo "check_storage_policy.sh: all assertions passed"

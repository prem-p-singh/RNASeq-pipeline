#!/usr/bin/env bash
# =============================================================================
# fetch_reference_file.sh — atomic, verified download of one reference file.
#
# R03 requires atomic transfer and format checks. The previous form,
#   curl -sSL "$url" | gunzip -c > "$out"
# wrote straight onto the final path, so an interrupted transfer left a
# truncated file sitting at the name everything else reads, and nothing checked
# that the payload was the format it claimed to be.
#
# Here the download lands beside the target, the archive is verified, the
# decompressed content is sniffed for its declared format, and only then is it
# renamed into place. Within one filesystem that rename is atomic, so a reader
# sees the final name only once the content is whole (master plan 5.3:
# "temporary build directory, validation, and atomic publication").
#
# A sha256 of the delivered bytes is written next to the file so the reference
# lock can record content identity without re-reading a large archive.
#
# Usage:
#   fetch_reference_file.sh --url URL --out PATH --kind fasta|gtf
# =============================================================================
set -euo pipefail

url=""; out=""; kind=""
while [ $# -gt 0 ]; do
    case "$1" in
        --url)  url=$2;  shift 2 ;;
        --out)  out=$2;  shift 2 ;;
        --kind) kind=$2; shift 2 ;;
        *) echo "Unknown arg: $1" >&2; exit 1 ;;
    esac
done
[ -n "$url" ] && [ -n "$out" ] && [ -n "$kind" ] || {
    echo "fetch_reference_file.sh: --url, --out and --kind are required" >&2
    exit 1
}

# A wildcard here cannot be expanded by curl and is a guaranteed failure.
# setup.py refuses to emit one; catch a hand-edited config too.
case "$url" in
    *'*'*) echo "Reference URL contains a wildcard, which curl cannot expand: $url" >&2
           exit 1 ;;
esac

outdir=$(dirname "$out")
mkdir -p "$outdir"

# Staging beside the destination, so the final move stays within one filesystem
# and is therefore atomic.
tmp_archive=$(mktemp "$outdir/.fetch.XXXXXX")
tmp_plain=$(mktemp "$outdir/.fetch.XXXXXX")
cleanup() { rm -f "$tmp_archive" "$tmp_plain"; }
trap cleanup EXIT

echo "[$(date -Iseconds)] fetching $url"
# --fail turns an HTTP error page into a non-zero exit instead of a file whose
# contents are HTML. --retry covers transient transport failures only.
curl -sSL --fail --retry 3 --retry-delay 2 -o "$tmp_archive" "$url"

if [ ! -s "$tmp_archive" ]; then
    echo "Downloaded file is empty: $url" >&2
    exit 1
fi

# --- Verify the archive before trusting its contents ---------------------
if gzip -t "$tmp_archive" 2>/dev/null; then
    gunzip -c "$tmp_archive" > "$tmp_plain"
else
    # Not gzip. Accept a plain file, but only if it is not an error page.
    case "$(head -c 512 "$tmp_archive" | tr -d '\0')" in
        *'<html'*|*'<HTML'*|*'<?xml'*)
            echo "Server returned markup, not a reference file: $url" >&2
            exit 1 ;;
    esac
    cp "$tmp_archive" "$tmp_plain"
fi

if [ ! -s "$tmp_plain" ]; then
    echo "Decompressed file is empty: $url" >&2
    exit 1
fi

# --- Sniff the declared format -------------------------------------------
# Cheap structural checks only: enough to catch a wrong URL or a truncated
# transfer, not a substitute for the identifier compatibility RF03 requires.
first=$(grep -m1 -v '^#' "$tmp_plain" || true)
case "$kind" in
    fasta)
        case "$first" in
            '>'*) : ;;
            *) echo "Expected FASTA (first record should start with '>'), got: ${first:0:80}" >&2
               exit 1 ;;
        esac
        ;;
    gtf)
        # A GTF data line has 9 tab-separated fields.
        n=$(printf '%s' "$first" | awk -F'\t' '{print NF}')
        if [ "${n:-0}" -lt 9 ]; then
            echo "Expected GTF (9 tab-separated fields), got $n in: ${first:0:80}" >&2
            exit 1
        fi
        ;;
    *) echo "Unknown --kind: $kind" >&2; exit 1 ;;
esac

# --- Record content identity, then publish atomically --------------------
if command -v sha256sum >/dev/null 2>&1; then
    sum=$(sha256sum "$tmp_plain" | awk '{print $1}')
else
    sum=$(shasum -a 256 "$tmp_plain" | awk '{print $1}')
fi
printf '%s  %s\n' "$sum" "$(basename "$out")" > "$out.sha256"

mv "$tmp_plain" "$out"
echo "[$(date -Iseconds)] $out  sha256=$sum  bytes=$(wc -c < "$out" | tr -d ' ')"

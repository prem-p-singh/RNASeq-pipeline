# Tools must come from the activated release runtime. No module fallback.
ensure_tools() {
    local missing=0 t
    for t in "$@"; do
        command -v "$t" >/dev/null 2>&1 || missing=1
    done
    [ "$missing" -eq 0 ] && return 0

    for t in "$@"; do
        command -v "$t" >/dev/null 2>&1 || {
            echo "Required tool absent from the activated release environment: $t. Run scripts/bootstrap.sh." >&2
            return 1
        }
    done
}

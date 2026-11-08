#!/usr/bin/env bash
# Build dist/rnaseq-pipeline_<version>_amd64.deb from the committed tree.
# Needs dpkg-deb (Debian/Ubuntu, or `brew install dpkg` on macOS).
set -euo pipefail
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
VERSION=$(sed -n 's/^version: "\(.*\)"/\1/p' "$REPO/CITATION.cff")
[ -n "$VERSION" ] || { echo "No version in CITATION.cff" >&2; exit 1; }

stage=$(mktemp -d)
trap 'rm -rf "$stage"' EXIT
lib=$stage/usr/lib/rnaseq-pipeline
mkdir -p "$lib" "$stage/usr/bin" "$stage/DEBIAN" "$stage/usr/share/doc/rnaseq-pipeline"

# Committed files only, minus README artwork, CI and packaging itself.
git -C "$REPO" archive HEAD -- . ':!assets' ':!.github' ':!packaging' | tar -x -C "$lib"
cp "$REPO/LICENSE" "$stage/usr/share/doc/rnaseq-pipeline/copyright"

# Wrappers rather than symlinks: submit.sh locates the workflow from its own path.
# The packaged Conda is used only when the caller has none of their own.
wrapper() {
    printf '%s\n' '#!/bin/sh' \
        'command -v conda >/dev/null 2>&1 || PATH=/opt/rnaseq-pipeline/conda/bin:$PATH' \
        'export PATH' \
        "exec bash /usr/lib/rnaseq-pipeline/$2 \"\$@\"" > "$stage/usr/bin/$1"
}
wrapper rnaseq-pipeline submit.sh
wrapper rnaseq-new-project scripts/new_project.sh

sed "s/@VERSION@/$VERSION/" "$REPO/packaging/control.in" > "$stage/DEBIAN/control"
install -m 755 "$REPO/packaging/postinst" "$REPO/packaging/postrm" "$stage/DEBIAN/"
chmod 755 "$stage"/usr/bin/*
chmod -R go-w "$stage"

mkdir -p "$REPO/dist"
out="$REPO/dist/rnaseq-pipeline_${VERSION}_amd64.deb"
dpkg-deb --root-owner-group --build "$stage" "$out" >&2
echo "$out"

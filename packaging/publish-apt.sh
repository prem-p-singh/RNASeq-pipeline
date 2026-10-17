#!/usr/bin/env bash
# publish-apt.sh <site-dir> <deb>... : add .debs to a flat apt repo in <site-dir>
# and sign it. Needs apt-utils, gnupg, and the signing key in the gpg keyring.
# Users: curl -fsSL <site>/rnaseq-pipeline.gpg | sudo tee /usr/share/keyrings/rnaseq-pipeline.gpg >/dev/null
#        echo "deb [signed-by=/usr/share/keyrings/rnaseq-pipeline.gpg] <site> ./" | sudo tee /etc/apt/sources.list.d/rnaseq-pipeline.list
set -euo pipefail
site=$1; shift
cp "$@" "$site/"
cd "$site"
apt-ftparchive packages . > Packages
gzip -kf Packages
apt-ftparchive release . > Release
gpg --batch --yes --clearsign -o InRelease Release
gpg --batch --yes -abs -o Release.gpg Release
gpg --export > rnaseq-pipeline.gpg

#!/usr/bin/env bash
# ON THE PUBLISHING MACHINE: publication of a package directory downloaded from the build server.
#   1. Zenodo: new version of the record (archive + checksum); the DOI goes to package/zenodo.json
#   2. GitHub Pages: the summary page with the new DOI, pushed to the gh-pages branch
#
#   ./publish_local.sh ~/Downloads/package              – production (zenodo.org)
#   SANDBOX=1 ./publish_local.sh ~/Downloads/package    – test (sandbox.zenodo.org)
#   NEW=1 ./publish_local.sh ~/Downloads/package        – first version of a new record (once)
#   NO_ZENODO=1 ./publish_local.sh ...                  – web page only (e.g. already published)
#   NO_SITE=1 ./publish_local.sh ...                    – Zenodo only, the web page is not pushed
#
# With SANDBOX=1 or NO_SITE=1 the web page is only built locally (site_preview/), never pushed, so
# a test DOI never reaches the public page.
#
# Run from the repository directory (zenodo.yaml, site.yaml and zenodo_state.json are here).
# Requires Python with pyyaml and pandas, git, and write access to the GitHub repository.
set -euo pipefail
PKG=$(cd "${1:?give the package directory}" && pwd)
cd "$(dirname "$0")"
[[ -f $PKG/package.json && -f $PKG/site/summary.json ]] || { echo "Not a package directory: $PKG" >&2; exit 1; }
ZEN=(); [[ -n ${SANDBOX:-} ]] && ZEN=(--sandbox)
[[ -n ${NEW:-} ]] && ZEN+=(--new)
log() { echo "[$(date '+%F %T')] $*"; }

if [[ -z ${NO_ZENODO:-} ]]; then
    log "Zenodo${SANDBOX:+ (sandbox)}"
    python3 -m argeement zenodo --package "$PKG" "${ZEN[@]}"
fi

if [[ -n ${SANDBOX:-} || -n ${NO_SITE:-} ]]; then
    rm -rf site_preview
    python3 -m argeement site --package "$PKG" --out site_preview --config site.yaml
    log "web page built locally, not pushed: $(pwd)/site_preview/index.html"
    exit 0
fi

SITE=$(mktemp -d)
python3 -m argeement site --package "$PKG" --out "$SITE" --config site.yaml
REMOTE=$(git remote get-url origin)
RUN=$(python3 -c "import json,sys;print(json.load(open(sys.argv[1]))['run'])" "$PKG/package.json")
(
    cd "$SITE"
    git init -q && git checkout -q -b gh-pages && git add -A
    git -c user.name="ARGeement" -c user.email="argeement@users.noreply.github.com" \
        commit -qm "Release $RUN"
    git push -qf "$REMOTE" gh-pages
)
rm -rf "$SITE"
DOI=$(python3 -c "import json,sys,os;p=sys.argv[1]+'/zenodo.json';print(json.load(open(p))['doi'] if os.path.exists(p) else '–')" "$PKG")
log "done: $RUN (DOI: $DOI); the page updates within a few minutes"

#!/usr/bin/env bash
# ON THE BUILD SERVER: test suite, update, and if any database has a new version, a new release with
# the package needed for publication (export archive, metadata, web page summary, figures).
# If the test suite fails, nothing is built.
#
#   THREADS=36 ./server_release.sh
#   FORCE=1 THREADS=36 ./server_release.sh      – new release even if no version has changed
#   PACKAGE_ONLY=1 ./server_release.sh          – no build, package the latest release
#   NO_DIFF=1 PACKAGE_ONLY=1 ./server_release.sh – the same, without the changes against the
#                                                  previous build (first published version)
#
# Cron (weekly, early Monday morning):
#   0 3 * * 1  cd /path/to/ARGeement && ARGEEMENT_ENV=/path/to/conda/envs/argeement THREADS=36 ./server_release.sh >> logs/release.log 2>&1
#
# At the end the location of the package directory is printed: download it to the machine used
# for publication and publish it there with publish_local.sh.
set -euo pipefail
cd "$(dirname "$0")"
ARGEEMENT_ENV=${ARGEEMENT_ENV:-$(dirname "$(dirname "$(command -v python)")")}
export PATH="$ARGEEMENT_ENV/bin:$PATH"
THREADS=${THREADS:-8}
PROFILE_FILE=${PROFILE_FILE:-config.yaml}
mkdir -p logs
log() { echo "[$(date '+%F %T')] $*"; }
PROFILE=$(python - "$PROFILE_FILE" << 'PY'
import sys, yaml
base = yaml.safe_load(open("config.yaml")); prof = yaml.safe_load(open(sys.argv[1])) or {}
print(prof.get("profile") or base.get("profile") or "current")
PY
)
BASE=releases/$PROFILE
PREV=$(readlink "$BASE/latest" 2>/dev/null || true)
if [[ -z ${PACKAGE_ONLY:-} ]]; then
    # the test suites (synthetic data; the whole Snakemake chain with two profiles) run before every
    # build; on failure no release is made
    TLOG=logs/tests_$(date +%Y%m%d_%H%M%S).log
    log "tests (log: $TLOG)"
    if ! { tests/run_test.sh && tests/run_pipeline_test.sh; } > "$TLOG" 2>&1; then
        log "TESTS FAILED, no release is built; see $TLOG"; tail -20 "$TLOG"; exit 1
    fi
    log "update ($PROFILE)"
    THREADS=$THREADS ./update.sh "$PROFILE_FILE"
fi
[[ -e $BASE/latest ]] || { log "no release"; exit 1; }
REL=$(readlink -f "$BASE/latest"); NEW=$(basename "$REL")
if [[ $NEW == "$PREV" && -z ${PACKAGE_ONLY:-} ]]; then
    log "no new release, nothing to do"; exit 0
fi
log "release: $REL"

if [[ ! -f $REL/validation_synthetic/summary.json ]]; then
    log "synthetic validation"
    python -m argeement validate --release "$REL" --n 500 --seed 1 --threads "$THREADS"
fi

FIG=()
if [[ -f manuscript/make_assets.py && -f manuscript/figures.R ]] && command -v Rscript > /dev/null; then   # optional
    log "figures for the web page"
    FT=$(mktemp -d)
    python manuscript/make_assets.py --release "$REL" --out "$FT" --lang en > /dev/null \
        && Rscript manuscript/figures.R "$FT" en > /dev/null 2>&1 \
        && { cp "$FT/en/szamok.json" "$FT/figures_en/macros.json" 2>/dev/null \
             || log "no szamok.json (old make_assets.py): figure captions without numbers"; } \
        && FIG=(--figures "$FT/figures_en") || log "figures not created (the page will have no figures)"
fi
ND=(); [[ -n ${NO_DIFF:-} ]] && ND=(--no-diff)
python -m argeement package --release "$REL" "${FIG[@]}" "${ND[@]}"
rm -rf "${FT:-}"
log "done. Download: $(pwd)/$BASE/$NEW/package/"

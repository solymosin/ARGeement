#!/usr/bin/env bash
# Update for one profile: download/copy, build, and if any database version has changed,
# a diff report against the previous release of the profile.
#
#   ./update.sh                       – default profile (config.yaml): the latest databases
#   ./update.sh profiles/aug2026.yaml – pinned profile
#   FORCE=1 ./update.sh ...           – new release even if no version has changed
#   CHECK_ONLY=1 ./update.sh ...      – only check whether a database has a new version (download
#                                       and parse, no build); exit code 0: no new version,
#                                       10: new version, run the update (e.g. server_release.sh)
#
# Releases: releases/<profile>/<run>; releases/<profile>/latest points to the newest one.
# Do not run two profiles at the same time (Snakemake locks the working directory).
set -euo pipefail
cd "$(dirname "$0")"
PROFILE_FILE=${1:-config.yaml}
THREADS=${THREADS:-4}
RUN_ID=$(date +%Y-%m-%d_%H%M%S)

PROFILE=$(python3 - "$PROFILE_FILE" << 'PY'
import sys, yaml
base = yaml.safe_load(open("config.yaml"))
prof = yaml.safe_load(open(sys.argv[1])) or {}
print(prof.get("profile") or base.get("profile") or "current")
PY
)
BASE="releases/$PROFILE"
mkdir -p "$BASE"
PREV=$(readlink "$BASE/latest" 2>/dev/null || true)

EXTRA=()
[[ "$PROFILE_FILE" != "config.yaml" ]] && EXTRA=(--configfile config.yaml "$PROFILE_FILE")
# 1. download and parse (cheap), then compare versions: if nothing has changed, the expensive
#    build is not started
if [[ -n "${CHECK_ONLY:-}" && -z "$PREV" ]]; then
    echo "[$PROFILE] No release yet: the first update builds one."
    exit 10
fi
if [[ -n "$PREV" && ( -z "${FORCE:-}" || -n "${CHECK_ONLY:-}" ) ]]; then
    snakemake --cores "$THREADS" "${EXTRA[@]}" --config run_id="$RUN_ID" threads="$THREADS" \
        --until parse inventory_pointfinder
    if python3 - "$BASE/$PREV/manifest.json" "work/$PROFILE/$RUN_ID/std" << 'PY'
import json, sys
from pathlib import Path
old = json.load(open(sys.argv[1]))["database_versions"]
std = Path(sys.argv[2])
new = {f.stem: f.read_text().strip() for f in std.glob("*.version")}
inv = std / "pointfinder.inventory.json"
if inv.exists():
    new["pointfinder"] = json.load(open(inv)).get("version", "")
# older manifests may lack a key (e.g. pointfinder): that is not treated as a change
same = all(old.get(k, v) == v for k, v in new.items()) and set(old) <= set(new)
print("versions:", ", ".join(f"{k} {v}" for k, v in sorted(new.items())))
for k in sorted(set(old) | set(new)):
    if old.get(k, new.get(k)) != new.get(k, old.get(k)):
        print(f"  changed: {k}: {old.get(k, '-')} -> {new.get(k, '-')}")
sys.exit(0 if same else 1)
PY
    then
        echo "[$PROFILE] No new database version, the build is not started."
        rm -rf "work/$PROFILE/$RUN_ID"
        exit 0
    fi
    if [[ -n "${CHECK_ONLY:-}" ]]; then
        echo "[$PROFILE] New database version; run the update to build a release."
        rm -rf "work/$PROFILE/$RUN_ID"
        exit 10
    fi
fi

# 2. full run (Snakemake does not repeat the steps already done)
snakemake --cores "$THREADS" "${EXTRA[@]}" --config run_id="$RUN_ID" threads="$THREADS"
NEW="$BASE/$RUN_ID"

if [[ -n "$PREV" && -z "${FORCE:-}" ]] && python3 - "$BASE/$PREV" "$NEW" << 'PY'
import json, sys
a, b = (json.load(open(f"{d}/manifest.json"))["database_versions"] for d in sys.argv[1:])
sys.exit(0 if a == b else 1)
PY
then
    echo "[$PROFILE] No new database version, the new run is discarded."
    rm -rf "$NEW" "work/$PROFILE/$RUN_ID"
    exit 0
fi

if [[ -n "$PREV" ]]; then
    python3 -m argeement diff "$BASE/$PREV" "$NEW" --out "$NEW/diff"
    echo "[$PROFILE] Diff report: $NEW/diff/diff_summary.md"
fi
ln -sfn "$RUN_ID" "$BASE/latest"
echo "[$PROFILE] New release: $NEW"

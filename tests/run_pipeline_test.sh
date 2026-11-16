#!/usr/bin/env bash
# Test of the whole Snakemake chain from local sources, with two profiles, without downloads:
# pinned profile (local directories + git ref), rerun without changes,
# a second profile, diff between profiles, and the version check.
set -euo pipefail
SRC="$(cd "$(dirname "$0")/.." && pwd)"
T=$(mktemp -d)
# copy the code only: on a build server the directory also holds work/, releases/ and local data
mkdir -p "$T/proj"
tar -C "$SRC" -cf - --exclude=./work --exclude=./releases --exclude=./.snakemake --exclude=./logs \
    --exclude=./manuscript --exclude=./docs --exclude=./saliva_dbs --exclude=./concordance/genomes \
    --exclude=./concordance/calls --exclude=./concordance/db --exclude='*/__pycache__' --exclude=./.git . \
    | tar -C "$T/proj" -xf -
cd "$T/proj"
python3 tests/make_test_data.py "$T/v1" > /dev/null
python3 tests/make_test_data.py "$T/v2" --v2 > /dev/null

# ResFinder as a git repository with two commits: the ref must select the older one
R="$T/resfinder_git"; mkdir -p "$R"
cp "$T"/v1/resfinder/*.fsa "$R"/ && echo "2.4.0" > "$R/VERSION"
git -C "$R" init -q && git -C "$R" add -A && git -C "$R" -c user.name=t -c user.email=t@t commit -qm v1
OLD=$(git -C "$R" rev-parse --short HEAD)
cp "$T"/v2/resfinder/*.fsa "$R"/ && echo "2.5.0" > "$R/VERSION"
git -C "$R" add -A && git -C "$R" -c user.name=t -c user.email=t@t commit -qm v2

cat > "$T/pinned.yaml" << YAML
profile: pinned
mmseqs: $(command -v mmseqs)
databases:
  card: {local: $T/v1/card, expect_version: "4.0.1"}
  amrfinder: {local: $T/v1/amrfinder, expect_version: "2026-08-15.1"}
  resfinder: {url: $R, ref: "$OLD", expect_version: "$OLD"}
  pointfinder: {local: $T/v1/pointfinder, expect_version: "4.1.1"}
  disinfinder: {local: $T/v1/disinfinder, expect_version: "09fc831"}
YAML
cat > "$T/newer.yaml" << YAML
profile: newer
mmseqs: $(command -v mmseqs)
databases:
  card: {local: $T/v2/card}
  amrfinder: {local: $T/v2/amrfinder}
  resfinder: {local: $T/v2/resfinder}
  pointfinder: {enabled: false}
  disinfinder: {enabled: false}
YAML

echo "== 1. pinned profile"
THREADS=2 ./update.sh "$T/pinned.yaml" > "$T/log1" 2>&1 || { cat "$T/log1"; exit 1; }
test -f releases/pinned/latest/manifest.json
grep -q '"profile": "pinned"' releases/pinned/latest/manifest.json
grep -q "2.4.0 (git $OLD)" releases/pinned/latest/manifest.json
test -f releases/pinned/latest/config.yaml
test -f releases/pinned/latest/inventory_pointfinder.json
grep -q "PointFinder" releases/pinned/latest/DATABASES.md
grep -q "1 AMR_DNA-\* files" releases/pinned/latest/DATABASES.md          # copied from the local source
grep -q '"disinfinder": "2.0.1 (git 09fc831)"' releases/pinned/latest/manifest.json
test ! -f releases/newer/latest/inventory_pointfinder.json 2>/dev/null || true
python3 - releases/pinned/latest/config.yaml << 'PY'
import sys, yaml
c = yaml.safe_load(open(sys.argv[1]))
# the profile overrides, but the rest of the base configuration is kept
assert c["profile"] == "pinned"
assert c["params"]["same_family_pident"] == 90.0, "params lost"
assert "AMR-SUSCEPTIBLE" in c["amrfinder"]["exclude_subtypes"], "amrfinder filters lost"
assert c["databases"]["card"]["url"].startswith("https://"), "card.url lost (non-recursive merge)"
assert c["databases"]["card"]["local"].endswith("/v1/card")
PY

echo "== 2. rerun without changes"
sleep 1
THREADS=2 ./update.sh "$T/pinned.yaml" 2>&1 | grep -q "No new database version"
test "$(ls releases/pinned | grep -vc latest)" -eq 1

echo "== 3. second profile with its own release series"
THREADS=2 ./update.sh "$T/newer.yaml" > "$T/log3" 2>&1 || { cat "$T/log3"; exit 1; }
test -f releases/newer/latest/manifest.json
if grep -q '"disinfinder"' releases/newer/latest/manifest.json; then echo "ERROR: disabled DisinFinder in the release"; exit 1; fi
test "$(readlink releases/pinned/latest)" != "$(readlink releases/newer/latest)"

echo "== 4. diff between profiles"
python3 -m argeement diff releases/pinned/latest releases/newer/latest --out "$T/xdiff" > /dev/null
grep -q "profile: pinned" "$T/xdiff/diff_summary.md"
grep -q "entry_added" "$T/xdiff/diff.tsv"

echo "== 4b. package with diff summary (second release of the profile)"
sleep 1
THREADS=2 FORCE=1 ./update.sh "$T/pinned.yaml" > "$T/log4b" 2>&1 || { cat "$T/log4b"; exit 1; }
python3 -m argeement package --release releases/pinned/latest > /dev/null 2>&1
grep -q "Changes since the previous release" releases/pinned/latest/package/RELEASE_NOTES.md
ls releases/pinned/latest/package/argeement-pinned-*.tar.gz > /dev/null

echo "== 5. lookup: prints the versions of the release"
python3 -m argeement lookup --release releases/pinned/latest --from resfinder --to card \
    --names blaTEM-1B 2>&1 >/dev/null | grep -q "\[pinned\]"

echo "== 6. different version: the run must stop"
sed 's/expect_version: "4.0.1"/expect_version: "9.9.9"/' "$T/pinned.yaml" > "$T/bad.yaml"
if THREADS=2 FORCE=1 ./update.sh "$T/bad.yaml" > "$T/log6" 2>&1; then
    echo "ERROR: the version check did not stop the run"; exit 1
fi
grep -q "but the profile expects" "$T/log6"

echo "== 7. ref with a plain (non-git) directory: the run must stop"
cat > "$T/badref.yaml" << YAML
profile: badref
mmseqs: $(command -v mmseqs)
databases:
  card: {local: $T/v1/card}
  amrfinder: {local: $T/v1/amrfinder}
  resfinder: {local: $T/v1/resfinder, ref: "$OLD"}
  pointfinder: {enabled: false}
  disinfinder: {enabled: false}
YAML
if THREADS=2 ./update.sh "$T/badref.yaml" > "$T/log7" 2>&1; then
    echo "ERROR: ref with a non-git directory did not stop the run"; exit 1
fi
grep -q "is not a git repository" "$T/log7"

rm -rf "$T"
echo "Pipeline test passed."

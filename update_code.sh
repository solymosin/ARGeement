#!/usr/bin/env bash
# ON THE BUILD SERVER: update the code from GitHub (or from a zip), keeping the local settings.
# The server needs no git repository and no GitHub credentials; publication runs elsewhere.
#
#   ./update_code.sh                  – latest code of the main branch on GitHub
#   ./update_code.sh v1.0.0           – a tag, branch or commit on GitHub
#   ./update_code.sh ARGeement.zip    – a local zip of the repository
#   DRY_RUN=1 ./update_code.sh ...    – only list what would change
#
# - config.yaml and existing profiles/*.yaml are never overwritten: the new version is written next
#   to them as <file>.new, and the difference of the settings (without comments) is shown.
# - Every other file of the repository is replaced; files not in the repository (work/, releases/,
#   manuscript/, concordance data, local profiles, local databases) are left alone.
# - Replaced files are backed up first, to ../argeement_code_backup_<date>.tar.gz.
set -euo pipefail

REPO_URL=${ARGEEMENT_REPO:-https://github.com/solymosin/ARGeement}
PROTECT_RE='^(config\.yaml|profiles/[^/]+\.yaml)$'

main() {
    local src_arg=${1:-main} dest dry=${DRY_RUN:-0}
    dest=$(cd "$(dirname "$0")" && pwd)
    tmp=$(mktemp -d)
    trap 'rm -rf "$tmp"' EXIT

    if [[ -f "$src_arg" ]]; then
        cp "$src_arg" "$tmp/code.zip"
    else
        local url=$src_arg
        [[ "$url" == http* ]] || url="$REPO_URL/archive/$src_arg.zip"
        echo "Download: $url"
        curl -fsSL -o "$tmp/code.zip" "$url" || { echo "Download failed: $url" >&2; exit 1; }
    fi
    unzip -q "$tmp/code.zip" -d "$tmp/x" 2>/dev/null || { echo "Not a zip archive: $src_arg" >&2; exit 1; }
    # the repository root is the directory that holds argeement/__init__.py
    local init; init=$(find "$tmp/x" -maxdepth 3 -path '*/argeement/__init__.py' | head -1)
    [[ -n "$init" ]] || { echo "Not an ARGeement repository archive." >&2; exit 1; }
    local src; src=$(dirname "$(dirname "$init")")

    local files=() f
    mapfile -t files < <(cd "$src" && find . -type f ! -path '*/__pycache__/*' | sed 's|^\./||' | sort)

    local changed=() kept=() new=() same=0
    for f in "${files[@]}"; do
        if [[ -e "$dest/$f" ]] && cmp -s "$src/$f" "$dest/$f"; then same=$((same + 1)); continue; fi
        if [[ -e "$dest/$f" && "$f" =~ $PROTECT_RE ]]; then kept+=("$f"); continue; fi
        if [[ -e "$dest/$f" ]]; then changed+=("$f"); else new+=("$f"); fi
    done
    echo "Target: $dest"
    echo "  new: ${#new[@]}, replaced: ${#changed[@]}, unchanged: $same, local settings kept: ${#kept[@]}"
    for f in "${changed[@]}"; do echo "  replace: $f"; done
    for f in "${new[@]}"; do echo "  new:     $f"; done
    for f in "${kept[@]}"; do echo "  keep:    $f (new version: $f.new)"; done
    if [[ "$dry" == 1 ]]; then echo "(DRY_RUN: nothing changed)"; return; fi

    if (( ${#changed[@]} + ${#kept[@]} )); then
        local backup; backup="$(dirname "$dest")/argeement_code_backup_$(date +%Y%m%d_%H%M%S).tar.gz"
        (cd "$dest" && tar czf "$backup" "${changed[@]}" "${kept[@]}")
        echo "Backup: $backup"
    fi
    for f in "${changed[@]}" "${new[@]}"; do
        mkdir -p "$dest/$(dirname "$f")"
        cp -p "$src/$f" "$dest/$f"
    done
    for f in "${kept[@]}"; do cp -p "$src/$f" "$dest/$f.new"; done
    find "$dest/argeement" "$dest/tests" -name __pycache__ -type d -prune -exec rm -rf {} + 2>/dev/null || true
    chmod +x "$dest"/*.sh "$dest"/tests/*.sh "$dest"/concordance/*.sh 2>/dev/null || true

    local nc='s/[[:space:]]*#.*$//; /^[[:space:]]*$/d'
    for f in "${kept[@]}"; do
        echo
        if diff -q <(sed -e "$nc" "$dest/$f") <(sed -e "$nc" "$dest/$f.new") >/dev/null; then
            echo "== $f: only the comments differ. To take the new comments: mv $f.new $f"
        else
            echo "== $f: the settings differ (< kept, > new; without comments):"
            diff <(sed -e "$nc" "$dest/$f") <(sed -e "$nc" "$dest/$f.new") || true
            echo "   Copy any new keys by hand, then: rm $f.new"
        fi
    done
    echo
    echo "Done. Check with: tests/run_test.sh && tests/run_pipeline_test.sh"
}

# the whole script is read before it runs, so it may replace itself safely
main "$@"
exit

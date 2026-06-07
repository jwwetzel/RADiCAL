#!/usr/bin/env bash
# organize_data.sh — consolidate the local RADiCAL data into the canonical
#   data/2023/{raw,reduced/<BUILD>} layout, using data/2023/MANIFEST.csv.
#
# Modes (run from the repo root):
#   ./Analysis/organize_data.sh            # default: symlink canonical -> legacy
#                                          #   (non-breaking; no duplication)
#   ./Analysis/organize_data.sh --copy     # copy legacy -> canonical (keeps legacy)
#   ./Analysis/organize_data.sh --move     # move legacy -> canonical (frees legacy;
#                                          #   only after the macros use DataPaths.h)
#   ./Analysis/organize_data.sh --check    # report which files are present, no changes
#
#   ./Analysis/organize_data.sh --legacy-links
#       For a NEWCOMER who downloaded data/2023/ from CERNBox: create the legacy
#       paths (Data/, reduced/, output/<E>GeV/) as symlinks INTO the canonical
#       files, so EVERY macro works (migrated or not) without moving anything.
#
# To prepare the CERNBox upload after linking:
#   rsync -aL data/2023/raw data/2023/reduced /path/to/upload/   # -L follows links
set -euo pipefail
cd "$(dirname "$0")/.."                      # repo root
MAN="data/2023/MANIFEST.csv"
MODE="${1:---link}"
[ -f "$MAN" ] || { echo "manifest not found: $MAN"; exit 1; }

n_ok=0; n_miss=0; n_done=0
while IFS=',' read -r build energy kind canonical legacy format; do
  case "$build" in '#'*|build|'') continue;; esac      # skip comments/header/blank
  if [ "$MODE" = "--legacy-links" ]; then
    # newcomer: canonical files exist (downloaded); create legacy -> canonical links
    if [ ! -e "$canonical" ]; then echo "  MISSING canonical: $canonical"; n_miss=$((n_miss+1)); continue; fi
    if [ -e "$legacy" ] || [ -L "$legacy" ]; then n_ok=$((n_ok+1)); continue; fi
    mkdir -p "$(dirname "$legacy")"
    ln -s "$(cd "$(dirname "$canonical")" && pwd)/$(basename "$canonical")" "$legacy"; n_done=$((n_done+1)); continue
  fi
  if [ -e "$canonical" ] || [ -L "$canonical" ]; then n_ok=$((n_ok+1)); continue; fi
  if [ ! -e "$legacy" ]; then
    echo "  MISSING legacy: $legacy   (-> $canonical)"; n_miss=$((n_miss+1)); continue
  fi
  mkdir -p "$(dirname "$canonical")"
  case "$MODE" in
    --check) echo "  would map: $legacy -> $canonical";;
    --copy)  cp    "$legacy" "$canonical"; n_done=$((n_done+1));;
    --move)  mv    "$legacy" "$canonical"; n_done=$((n_done+1));;
    --link)  ln -s "$(cd "$(dirname "$legacy")" && pwd)/$(basename "$legacy")" "$canonical"; n_done=$((n_done+1));;
    *) echo "unknown mode: $MODE"; exit 1;;
  esac
done < "$MAN"

echo "----"
echo "mode=$MODE   already-present=$n_ok   newly-mapped=$n_done   missing-legacy=$n_miss"
[ "$MODE" = "--check" ] && echo "(dry run — no changes made)"
echo "canonical tree now under data/2023/{raw,reduced}/"

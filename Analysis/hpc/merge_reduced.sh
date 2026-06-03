#!/bin/bash
# ============================================================================
# merge_reduced.sh — login-node hadd of the per-run reduced files into one
# compact per-ENERGY file per config, placed in the repo's reduced/ tree ready
# to copy home.  Fast (the reduced files are tiny); no SGE needed.
#
#   export RAD_CONFIG=LUAG
#   bash Analysis/hpc/merge_reduced.sh
#
# Output:  reduced/<RAD_CONFIG>/<label>.root   (e.g. reduced/LUAG/150GeV.root)
# ============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/env.sh"
setup_root
: "${RAD_CONFIG:?set RAD_CONFIG to the capillary config name (e.g. DSB1, LUAG)}"

src="$RAD_WORK/reduced/${RAD_CONFIG}/byrun"
dst="$RAD_REPO/reduced/${RAD_CONFIG}"
[ -d "$src" ] || { echo "no per-run output at $src — run sge_reduce first"; exit 1; }
mkdir -p "$dst"

for d in "$src"/*/; do
    [ -d "$d" ] || continue
    label=$(basename "$d")
    n=$(ls "$d"/RUN*.root 2>/dev/null | wc -l | tr -d ' ')
    [ "$n" -eq 0 ] && { echo "  $label: no runs, skip"; continue; }
    echo "  $label: merging $n runs -> $dst/${label}.root"
    hadd -f "$dst/${label}.root" "$d"/RUN*.root >/dev/null
done

echo "merged config '${RAD_CONFIG}' -> $dst"
ls -lh "$dst"
echo
echo "Copy home:  rsync -av <argon>:$dst/  ./reduced/${RAD_CONFIG}/"

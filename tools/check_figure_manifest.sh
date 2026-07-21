#!/usr/bin/env bash
# check_figure_manifest.sh — figure-drift checker (REORG_PLAN Phase B8).
#
# Verifies that every file listed in papers/timing/figs/MANIFEST.md still has the
# md5 recorded there. A mismatch means the committed manuscript figures have drifted
# from the recorded provenance — treat it as a FINDING, not a formality: either a
# gate legitimately regenerated the figure (then update MANIFEST.md in the same
# commit) or something touched a frozen file (dist.png / optimization.png must
# never change silently).
#
# Usage: tools/check_figure_manifest.sh          (from the repo root)
# Exit:  0 = all match; 1 = drift or missing file; 2 = environment problem.
set -u

MANIFEST="papers/timing/figs/MANIFEST.md"
FIGDIR="papers/timing/figs"

if [[ ! -f "$MANIFEST" ]]; then
  echo "ERROR: $MANIFEST not found — run from the repo root." >&2
  exit 2
fi

# md5 first-12-hex, portable across macOS (md5) and Linux (md5sum)
md5_12() {
  if command -v md5 >/dev/null 2>&1; then
    md5 -q "$1" | cut -c1-12
  else
    md5sum "$1" | cut -c1-12
  fi
}

fail=0
checked=0

# Parse table rows: | `file` | `md5(12)` | ...
while IFS= read -r line; do
  file=$(printf '%s' "$line" | sed -n 's/^| *`\([^`]*\)` *| *`\([0-9a-f]\{12\}\)`.*/\1/p')
  want=$(printf '%s' "$line" | sed -n 's/^| *`\([^`]*\)` *| *`\([0-9a-f]\{12\}\)`.*/\2/p')
  [[ -z "$file" || -z "$want" ]] && continue
  checked=$((checked + 1))
  path="$FIGDIR/$file"
  if [[ ! -f "$path" ]]; then
    echo "MISSING  $file  (manifest lists $want)"
    fail=1
    continue
  fi
  have=$(md5_12 "$path")
  if [[ "$have" == "$want" ]]; then
    echo "OK       $file  $have"
  else
    echo "DRIFT    $file  manifest=$want  tree=$have"
    fail=1
  fi
done < "$MANIFEST"

if [[ "$checked" -eq 0 ]]; then
  echo "ERROR: no manifest rows parsed — table format changed?" >&2
  exit 2
fi

echo "---"
if [[ "$fail" -eq 0 ]]; then
  echo "PASS: all $checked manifest entries match the tree."
else
  echo "FAIL: drift detected. If a gate legitimately regenerated a figure, update"
  echo "MANIFEST.md in the same commit; otherwise investigate before committing."
fi
exit $fail

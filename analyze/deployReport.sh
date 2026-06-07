#!/usr/bin/env bash
# ============================================================================
# deployReport.sh — publish the generated analysis report into the site/ folder
# (site/report/, served at <pages-url>/report/ via the Pages Actions workflow).
#
# output/ is a gitignored build dir, so Pages can't serve from it.
# This copies ONLY the assets report.html references into site/report/, and
# downscales the oversized multi-page renders so we don't bloat git history.
#
# Run from the repo root, after makeReport.py:   bash analyze/deployReport.sh
# ============================================================================
set -euo pipefail

SRC="output"
DST="site/report"

[ -f "$SRC/report.html" ] || { echo "No $SRC/report.html — run 'python3 analyze/makeReport.py' first."; exit 1; }

rm -rf "$DST"
mkdir -p "$DST"
cp "$SRC/report.html" "$DST/index.html"

# Copy only the Summary/ and report_images/ files the report actually references,
# as-is — the renderer's PNGs are already well-compressed (re-encoding bloats them).
grep -oE '(src|href)="(Summary|report_images)/[^"]+"' "$SRC/report.html" \
  | sed -E 's/.*"([^"]+)"/\1/' | sort -u | while read -r f; do
    [ -f "$SRC/$f" ] || { echo "  WARN missing asset: $f"; continue; }
    mkdir -p "$DST/$(dirname "$f")"
    cp "$SRC/$f" "$DST/$f"
  done

echo "Deployed $(find "$DST" -type f | wc -l | tr -d ' ') files ($(du -sh "$DST" | cut -f1)) to $DST/"
echo "Local preview: open $DST/index.html   |   Live: <pages-url>/RADiCAL/report/"

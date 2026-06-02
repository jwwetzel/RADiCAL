#!/usr/bin/env bash
# ============================================================================
# deployReport.sh — publish the generated analysis report to the committed
# Pages folder (report/ at the repo root, served at <pages>/RADiCAL/report/).
#
# Analysis/Output/ is a gitignored build dir, so Pages can't serve from it.
# This copies ONLY the assets report.html references into report/, and
# downscales the oversized multi-page renders so we don't bloat git history.
#
# Run from the repo root, after makeReport.py:   bash Analysis/deployReport.sh
# ============================================================================
set -euo pipefail

SRC="Analysis/Output"
DST="report"

[ -f "$SRC/report.html" ] || { echo "No $SRC/report.html — run 'python3 Analysis/makeReport.py' first."; exit 1; }

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

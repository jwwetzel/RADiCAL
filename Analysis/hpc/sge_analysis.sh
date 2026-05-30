#!/bin/bash
# ============================================================================
# sge_analysis.sh (SGE) — run the full downstream analysis on the merged
# per-energy ntuples.  SKIP_PROCESS=1 makes runAll.sh skip raw->ntuple Step 1
# (already done by the array job + merge) and run every analysis macro + report.
#
# Single sequential job (the macros loop the merged ntuples; no ACLiC race
# since Step 1's parallel section is skipped).  Give it generous time/memory.
#
# Needs python3 + a PDF rasterizer (pdftoppm/poppler) for the HTML report's PNG
# embedding.  If those aren't on the nodes, the macros still emit all PDFs and
# results.json; regenerate report.html elsewhere with makeReport.py.
# ============================================================================
#$ -N rad_analysis
#$ -cwd
#$ -j y
set -uo pipefail
source Analysis/hpc/env.sh
setup_root

export SKIP_PROCESS=1
bash Analysis/runAll.sh
rc=$?
echo "sge_analysis.sh: runAll exit=$rc"
echo "Outputs under: $(readlink -f Analysis/Output)"
exit $rc

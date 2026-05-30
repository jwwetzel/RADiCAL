#!/bin/bash
# ============================================================================
# compile.sh (SGE) — one-shot setup: point Analysis/Output at scratch and
# pre-compile processRun.C so the 252 array tasks all LOAD the same .so instead
# of racing to build it (the ACLiC race documented in runAll.sh).
# Submitted first; the array job -hold_jid's on it.
# ============================================================================
#$ -N rad_compile
#$ -cwd
#$ -j y
set -euo pipefail
source Analysis/hpc/env.sh
setup_root

mkdir -p "$RAD_WORK/Output" "$RAD_WORK/logs"
# Redirect all analysis output to scratch (transparent to every macro, which
# all write under Analysis/Output/).  Only replaces a symlink, never a real dir.
if [ -e Analysis/Output ] && [ ! -L Analysis/Output ]; then
    echo "ERROR: Analysis/Output exists as a real directory; move/remove it first." >&2
    exit 1
fi
ln -sfn "$RAD_WORK/Output" Analysis/Output

root -l -b -q -e '.L Analysis/processRun.C+' || true
ls -l Analysis/processRun_C*.so
echo "compile.sh done."

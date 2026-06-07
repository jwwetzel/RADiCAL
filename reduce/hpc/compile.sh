#!/bin/bash
# ============================================================================
# compile.sh (SGE) — prebuild the canonical reducer .so ONCE, so the reduce
# array tasks LOAD it instead of each racing to ACLiC-compile it (the ACLiC
# race). Submit this first; the reduce array -hold_jid's on rad_compile.
# ============================================================================
#$ -N rad_compile
#$ -cwd
#$ -j y
set -euo pipefail
source reduce/hpc/env.sh
setup_root

mkdir -p "$RAD_WORK/logs"

# Canonical config-driven reducer (reduce/reduceRun.C -> reduce/Reducer.C).
# lib/ headers must be on the include path (the four domain subdirs).
ROOT_INCLUDE_PATH="$RAD_REPO/lib/waveform:$RAD_REPO/lib/io:$RAD_REPO/lib/physics:$RAD_REPO/lib/viz" \
  root -l -b -q -e '.L reduce/reduceRun.C+'

ls -l reduce/reduceRun_C*.so
echo "compile.sh done."

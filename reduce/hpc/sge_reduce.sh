#!/bin/bash
# ============================================================================
# sge_reduce.sh (SGE ARRAY job) — config-AGNOSTIC raw->reduced, one task per run.
# Task $SGE_TASK_ID reads that line of $TASKLIST and runs reduceRaw.C on its raw
# file, writing a compact per-run reduced ntuple under
#   $RAD_WORK/reduced/$RAD_CONFIG/byrun/<label>/RUN<run>.root
#
# Unlike sge_process.sh (DSB1-specific processRun), this needs NO channel map:
# reduceRaw stores every slot's pulse features + the (config-invariant) WC track
# and both MCPs.  So the SAME job handles DSB1, LUAG, mixed, timing+energy —
# you just point it at that config's tasklist via $MANIFEST/$RAD_TASKLIST and
# tag the output with $RAD_CONFIG.
#
#   export RAD_CONFIG=LUAG
#   export MANIFEST=$RAD_REPO/reduce/hpc/manifest_LUAG.csv
#   export RAD_TASKLIST=$RAD_WORK/tasks_LUAG.txt
#   bash reduce/hpc/discover_tasklist.sh
#   qsub -hold_jid rad_compile -t 1-$(wc -l < $RAD_TASKLIST) reduce/hpc/sge_reduce.sh
# ============================================================================
#$ -N rad_reduce
#$ -cwd
#$ -j y
set -euo pipefail
source reduce/hpc/env.sh
setup_root

: "${RAD_CONFIG:?set RAD_CONFIG to the capillary config name (e.g. DSB1, LUAG)}"

line=$(sed -n "${SGE_TASK_ID}p" "$TASKLIST")
[ -z "$line" ] && { echo "no task on line $SGE_TASK_ID"; exit 1; }
run=$(   printf '%s' "$line" | cut -f1)
label=$( printf '%s' "$line" | cut -f2)
energy=$(printf '%s' "$line" | cut -f3)
raw=$(   printf '%s' "$line" | cut -f4)

outdir="$RAD_WORK/reduced/${RAD_CONFIG}/byrun/${label}"
mkdir -p "$outdir"
out="${outdir}/RUN${run}.root"
echo "[reduce ${SGE_TASK_ID}] cfg=${RAD_CONFIG} run=${run} E=${energy} GeV -> ${out}"
echo "  raw: ${raw}"

# CANONICAL config-driven reduction (reduce/Reducer): proper role-resolved
# hg_cfd05 etc. for EVERY build, via the build's channel map. .so prebuilt by
# compile.sh; lib/ headers on the include path.
cfg="$RAD_REPO/data/${RAD_YEAR}/configs/${RAD_CONFIG}.json"
[ -f "$cfg" ] || { echo "config not found: $cfg (expected data/${RAD_YEAR}/configs/${RAD_CONFIG}.json)"; exit 1; }
ROOT_INCLUDE_PATH="$RAD_REPO/lib/waveform:$RAD_REPO/lib/io:$RAD_REPO/lib/physics:$RAD_REPO/lib/viz" \
  root -l -b -q "reduce/reduceRun.C+(\"${cfg}\", \"${raw}\", ${energy}., \"${out}\")"
echo "[reduce ${SGE_TASK_ID}] done."

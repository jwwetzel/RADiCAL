#!/bin/bash
# ============================================================================
# submit_reduce.sh — reduce ONE capillary config on Argon, end to end:
#   resolve files -> reduce[array] -> merge per energy.
# (Run compile.sh ONCE first; it prebuilds reduceRaw.so for every config.)
#
# Usage (from repo root on the Argon LOGIN node):
#   qsub Analysis/hpc/compile.sh                 # once, wait for it to finish
#   python3 Analysis/hpc/build_manifest.py "$LOGBOOK" --capillary "LUAG" \
#           -o Analysis/hpc/manifest_LUAG.csv
#   RAD_CONFIG=LUAG RAD_MANIFEST=$PWD/Analysis/hpc/manifest_LUAG.csv \
#       bash Analysis/hpc/submit_reduce.sh
#
# Output: reduced/<RAD_CONFIG>/<E>GeV.root   (copy home when merge finishes)
# ============================================================================
set -euo pipefail
: "${RAD_CONFIG:?set RAD_CONFIG, e.g. RAD_CONFIG=LUAG}"
: "${RAD_MANIFEST:?set RAD_MANIFEST=/abs/path/manifest_<config>.csv}"
[ -f "$RAD_MANIFEST" ] || { echo "manifest not found: $RAD_MANIFEST"; exit 1; }

source Analysis/hpc/env.sh                     # MANIFEST<-RAD_MANIFEST (now respected)
# Give each config its own tasklist so concurrent configs never collide.
export RAD_TASKLIST="${RAD_TASKLIST:-$RAD_WORK/tasks_${RAD_CONFIG}.txt}"

echo "config=$RAD_CONFIG"
echo "manifest=$MANIFEST"
echo "tasklist=$RAD_TASKLIST"
echo "rec_dir=$REC_DIR"

# 1) resolve manifest runs -> raw file paths (login node, no ROOT needed)
bash Analysis/hpc/discover_tasklist.sh
N=$(wc -l < "$RAD_TASKLIST")
[ "$N" -gt 0 ] || { echo "ERROR: 0 tasks resolved — check REC_DIR and the manifest."; exit 1; }
mkdir -p "$RAD_WORK/logs"

# 2) submit reduce array + dependent merge, passing the env through with -V
COMMON=(-V -cwd -q "$RAD_QUEUE" -o "$RAD_WORK/logs" -e "$RAD_WORK/logs")
RES=(); [ -n "${RAD_MEM:-}" ] && RES+=(-l "mem_free=$RAD_MEM")
        [ -n "${RAD_HRT:-}" ] && RES+=(-l "h_rt=$RAD_HRT")
strip(){ echo "$1" | cut -d. -f1; }

# -hold_jid rad_compile: waits if compile.sh is still queued/running; if it has
# already finished (no such job) SGE proceeds immediately — so you can fire
# `qsub compile.sh` and this back-to-back without timing them by hand.
jid_r=$(qsub -terse "${COMMON[@]}" -N "rad_reduce_${RAD_CONFIG}" -hold_jid rad_compile \
             -t 1-"$N" "${RES[@]}" Analysis/hpc/sge_reduce.sh)
echo "reduce : $jid_r   ($N tasks)"
jid_m=$(qsub -terse "${COMMON[@]}" -N "rad_merge_${RAD_CONFIG}" \
             -hold_jid "$(strip "$jid_r")" Analysis/hpc/merge_reduced.sh)
echo "merge  : $jid_m"

echo
echo "Monitor : qstat -u \"$USER\""
echo "Result  : reduced/${RAD_CONFIG}/<E>GeV.root   (after merge finishes)"
echo "Copy home: rsync -av <argon>:$RAD_REPO/reduced/${RAD_CONFIG}/  ./reduced/${RAD_CONFIG}/"

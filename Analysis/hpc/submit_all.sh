#!/bin/bash
# ============================================================================
# submit_all.sh — submit the whole DSB1 pipeline to Argon (SGE) with the right
# job dependencies:  compile -> process[array] -> merge -> analysis.
# Run from the repo root on the Argon LOGIN node, AFTER discover_tasklist.sh:
#     bash Analysis/hpc/submit_all.sh
# ============================================================================
set -euo pipefail
source Analysis/hpc/env.sh

[ -s "$TASKLIST" ] || { echo "ERROR: $TASKLIST missing/empty. Run: bash Analysis/hpc/discover_tasklist.sh" >&2; exit 1; }
N=$(wc -l < "$TASKLIST")
mkdir -p "$RAD_WORK/logs"

# Common qsub flags.  -V passes the current (exported) environment so the jobs
# see RAD_REPO/REC_DIR/RAD_WORK/... -o/-e send logs to scratch.
COMMON=(-V -cwd -q "$RAD_QUEUE" -o "$RAD_WORK/logs" -e "$RAD_WORK/logs")
# Optional resource requests (blank them in env.sh if Argon rejects these names).
PROC_L=(); [ -n "${RAD_MEM:-}" ] && PROC_L+=(-l "mem_free=$RAD_MEM"); [ -n "${RAD_HRT:-}" ] && PROC_L+=(-l "h_rt=$RAD_HRT")
ANA_L=(-l "h_rt=12:00:00"); [ -n "${RAD_MEM:-}" ] && ANA_L+=(-l "mem_free=16G")

strip() { echo "$1" | cut -d. -f1; }   # array -terse id "123.1-252:1" -> "123"

jid_c=$(qsub -terse "${COMMON[@]}" Analysis/hpc/compile.sh)
echo "compile  : $jid_c"

jid_p=$(qsub -terse "${COMMON[@]}" -hold_jid "$(strip "$jid_c")" \
             -t 1-"$N" "${PROC_L[@]}" Analysis/hpc/sge_process.sh)
echo "process  : $jid_p   ($N array tasks)"

jid_m=$(qsub -terse "${COMMON[@]}" -hold_jid "$(strip "$jid_p")" Analysis/hpc/sge_merge.sh)
echo "merge    : $jid_m"

jid_a=$(qsub -terse "${COMMON[@]}" -hold_jid "$(strip "$jid_m")" "${ANA_L[@]}" Analysis/hpc/sge_analysis.sh)
echo "analysis : $jid_a"

echo
echo "Submitted. Monitor with:  qstat -u \"$USER\""
echo "Logs:      $RAD_WORK/logs/"
echo "Results:   $RAD_WORK/Output/  (report.html, Summary/, per-energy PDFs)"

#!/bin/bash
# ============================================================================
# sge_process.sh (SGE ARRAY job) — one task per run.  Task $SGE_TASK_ID reads
# that line of $TASKLIST and runs processRun on its raw file, writing a per-run
# ntuple to  Analysis/Output/byrun/<label>/RUN<run>/ntuple.root.
#
# Submit with the array range set to the number of tasks, e.g.
#   qsub -t 1-$(wc -l < $TASKLIST) Analysis/hpc/sge_process.sh
# (submit_all.sh does this for you).  Embarrassingly parallel.
# ============================================================================
#$ -N rad_proc
#$ -cwd
#$ -j y
set -euo pipefail
source Analysis/hpc/env.sh
setup_root

line=$(sed -n "${SGE_TASK_ID}p" "$TASKLIST")
[ -z "$line" ] && { echo "no task on line $SGE_TASK_ID"; exit 1; }
run=$(  printf '%s' "$line" | cut -f1)
label=$(printf '%s' "$line" | cut -f2)
energy=$(printf '%s' "$line" | cut -f3)
raw=$(  printf '%s' "$line" | cut -f4)

outlabel="byrun/${label}/RUN${run}"
echo "[task ${SGE_TASK_ID}] run=${run} E=${energy} GeV  -> Analysis/Output/${outlabel}"
echo "  raw: ${raw}"

# .so was built by compile.sh; ACLiC sees it current and just loads it.
root -l -b -q "Analysis/processRun.C+(\"${raw}\", ${energy}., \"${outlabel}\")"
echo "[task ${SGE_TASK_ID}] done."

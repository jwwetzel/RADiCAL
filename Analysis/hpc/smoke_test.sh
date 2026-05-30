#!/bin/bash
# ============================================================================
# smoke_test.sh — run the WHOLE pipeline on a tiny subset (default 2 runs per
# energy), serially, end-to-end, BEFORE committing the full 221-run submission.
#
# Validates the three things a subset is uniquely good at catching:
#   1. processRun actually runs on the real Argon `rec` files (raw "pulse" tree
#      schema matches) and produces a valid ntuple.
#   2. hadd merges the per-run ntuples cleanly.
#   3. every analysis macro compiles under Argon's ROOT (6.24) and makeReport runs.
#
# Writes to a SEPARATE output dir ($RAD_WORK/Output_smoke) so it never clobbers
# the real run, and restores the Output symlink at the end.
#
# RUN ON A COMPUTE NODE, not the login node:
#     qlogin                                  # grab an interactive node
#     cd <repo> && bash Analysis/hpc/smoke_test.sh   [N_runs_per_energy]
# ============================================================================
set -uo pipefail
source Analysis/hpc/env.sh
setup_root

N=${1:-2}                                  # runs per energy
[ -s "$TASKLIST" ] || { echo "ERROR: $TASKLIST missing — run discover_tasklist.sh first" >&2; exit 1; }

if [ -e Analysis/Output ] && [ ! -L Analysis/Output ]; then
    echo "ERROR: Analysis/Output is a real directory; move/remove it first." >&2; exit 1
fi

SMOKE_OUT="$RAD_WORK/Output_smoke"
SMOKE_TASKS="$RAD_WORK/tasks_smoke.txt"
mkdir -p "$SMOKE_OUT"
ln -sfn "$SMOKE_OUT" Analysis/Output

# First N runs of each energy (awk exact-match on the label column).
: > "$SMOKE_TASKS"
for lab in $(cut -f2 "$TASKLIST" | sort -u); do
    awk -F'\t' -v l="$lab" '$2==l' "$TASKLIST" | head -n "$N" >> "$SMOKE_TASKS"
done
echo "=== smoke subset: $(wc -l < "$SMOKE_TASKS") runs ($N per energy) ==="
cat "$SMOKE_TASKS"

echo "=== [1/4] compile processRun ==="
root -l -b -q -e '.L Analysis/processRun.C+' || true
ls -l Analysis/processRun_C*.so || { echo "compile FAILED"; exit 1; }

echo "=== [2/4] processRun on the subset (serial) ==="
while IFS=$'\t' read -r run label energy raw; do
    [ -z "$run" ] && continue
    echo "  run=$run E=$energy GeV"
    root -l -b -q "Analysis/processRun.C+(\"$raw\", ${energy}., \"byrun/${label}/RUN${run}\")" \
        || { echo "processRun FAILED on run $run — STOP, schema mismatch likely"; exit 1; }
done < "$SMOKE_TASKS"

echo "=== [3/4] merge per energy ==="
for lab in $(cut -f2 "$SMOKE_TASKS" | sort -u); do
    shopt -s nullglob; files=( Analysis/Output/byrun/${lab}/RUN*/ntuple.root ); shopt -u nullglob
    [ ${#files[@]} -eq 0 ] && continue
    mkdir -p "Analysis/Output/${lab}"
    hadd -f "Analysis/Output/${lab}/ntuple.root" "${files[@]}" || { echo "hadd FAILED for $lab"; exit 1; }
done

echo "=== [4/4] full analysis + report (SKIP_PROCESS) ==="
SKIP_PROCESS=1 bash Analysis/runAll.sh
rc=$?

# restore the symlink to the REAL output target for the production run
ln -sfn "$RAD_WORK/Output" Analysis/Output

echo
echo "=== SMOKE TEST exit=$rc ==="
echo "Inspect: $SMOKE_OUT/report.html , $SMOKE_OUT/Summary/"
echo "If clean -> production run:  bash Analysis/hpc/submit_all.sh"
exit $rc

#!/bin/bash
# ============================================================================
# discover_tasklist.sh — resolve the manifest run numbers to actual files on
# Argon and emit the per-run task list the array job indexes by $SGE_TASK_ID.
#
# Run once on the Argon LOGIN node (fast, no ROOT needed):
#     bash Analysis/hpc/discover_tasklist.sh
#
# Output: $TASKLIST  (one line per existing run file)
#     <run>\t<label>\t<energy>\t<absolute_raw_path>
# Runs in the manifest with no matching file are reported and skipped.
# ============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/env.sh"

mkdir -p "$RAD_WORK"
: > "$TASKLIST"
missing=0; found=0

# skip header; manifest columns: run,energy_GeV,label,v_up,v_down,beam,capillary,status
tail -n +2 "$MANIFEST" | while IFS=, read -r run energy label vup vdn beam cap status; do
    [ -z "$run" ] && continue
    f=""
    # Try the most likely names first: bare, 4-digit and 5-digit zero-padded.
    for cand in \
        "$REC_DIR/RUN${run}.root" \
        "$REC_DIR/RUN$(printf '%04d' "$run").root" \
        "$REC_DIR/RUN$(printf '%05d' "$run").root"; do
        if [ -f "$cand" ]; then f="$cand"; break; fi
    done
    # Fallback: a unique RUN<run> file with an energy/suffix (e.g. RUN1034_125_GeV.root).
    if [ -z "$f" ]; then
        match=$(ls "$REC_DIR"/RUN0*"${run}".root "$REC_DIR"/RUN"${run}"_*.root \
                   "$REC_DIR"/RUN"${run}".root 2>/dev/null | sort -u || true)
        if [ "$(printf '%s\n' "$match" | grep -c .)" = "1" ]; then f="$match"; fi
    fi
    if [ -n "$f" ]; then
        printf '%s\t%s\t%s\t%s\n' "$run" "$label" "$energy" "$f" >> "$TASKLIST"
        found=$((found+1))
    else
        echo "WARN: no file for run $run ($label)" >&2
        missing=$((missing+1))
    fi
done

n=$(wc -l < "$TASKLIST")
echo "Wrote $TASKLIST : $n tasks"
echo "Per energy:"; cut -f2 "$TASKLIST" | sort | uniq -c
[ "$missing" -gt 0 ] && echo "($missing manifest runs had no file — see WARN lines above)"
echo
echo "Next:  qsub -t 1-$n Analysis/hpc/sge_process.sh    (or: bash Analysis/hpc/submit_all.sh)"

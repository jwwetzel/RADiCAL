#!/bin/bash
# ============================================================================
# sge_merge.sh (SGE) — hadd the per-run ntuples of each energy into the single
# per-energy ntuple the analysis macros expect:
#     Analysis/Output/byrun/<label>/RUN*/ntuple.root  ->  Analysis/Output/<label>/ntuple.root
# The "rad" trees concatenate cleanly (identical schema; each carries its run #).
# ============================================================================
#$ -N rad_merge
#$ -cwd
#$ -j y
set -euo pipefail
source Analysis/hpc/env.sh
setup_root

for label in $(cut -f2 "$TASKLIST" | sort -u); do
    outdir="Analysis/Output/${label}"
    mkdir -p "$outdir"
    shopt -s nullglob
    files=( Analysis/Output/byrun/${label}/RUN*/ntuple.root )
    shopt -u nullglob
    if [ ${#files[@]} -eq 0 ]; then
        echo "WARN: no per-run ntuples for ${label} — skipping" >&2
        continue
    fi
    echo "[${label}] merging ${#files[@]} run-ntuples -> ${outdir}/ntuple.root"
    # -k: skip corrupt/zombie inputs (e.g. a processRun task killed mid-write by a
    # full disk) instead of aborting the whole merge.
    hadd -k -f -j 4 "${outdir}/ntuple.root" "${files[@]}"
done
echo "sge_merge.sh done."
echo "Per-energy event counts:"
for label in $(cut -f2 "$TASKLIST" | sort -u); do
    f="Analysis/Output/${label}/ntuple.root"
    [ -f "$f" ] && echo "  ${label}: $(root -l -b -q -e "TFile f(\"$f\"); printf(\"%lld\\n\",((TTree*)f.Get(\"rad\"))->GetEntries());" 2>/dev/null | tail -1)"
done

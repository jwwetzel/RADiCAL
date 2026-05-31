#!/bin/bash
# ============================================================================
# env.sh — shared configuration + ROOT environment for the Argon HPC pipeline.
# Every job script sources this.  EDIT the paths and the ROOT block for your
# Argon account (see README.md "Step 0: find ROOT").
# ============================================================================

# --- paths -----------------------------------------------------------------
# Repository root (where the Analysis/ directory lives), derived from THIS
# file's location (Analysis/hpc/env.sh -> ../.. = repo root).  Works wherever
# the repo is cloned and in -cwd batch jobs; override by exporting RAD_REPO.
_ENV_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
# Always re-derived (NOT ${RAD_REPO:-...}) so re-sourcing self-heals any stale
# value left exported in your shell by a previous source of this file.
export RAD_REPO="$(cd "$_ENV_DIR/../.." && pwd)"

# Raw reconstructed run files (READ-only input): <REC_DIR>/RUN<run>.root.
# These live on the group LSS allocation.
export REC_DIR="${REC_DIR:-/Shared/lss_yonel/jwwetzel/RADiCAL_CERN_May2023/rec/rec}"

# Work area for ALL outputs (per-run ntuples, merged ntuples, plots, report).
# Put this on roomy, fast, SHARED storage — Argon's /nfsscratch (303 TB) — NOT on
# the group LSS allocation, which is small and fills up (the raw data already
# lives there).  /nfsscratch is auto-purged periodically, so copy the final
# report.html + Summary/ + merged ntuples to permanent storage when done.
export RAD_WORK="${RAD_WORK:-/nfsscratch/$USER/RADiCAL_CERN_May2023/analysis_work}"

# Manifest (committed in repo) and the resolved per-run task list (made on Argon).
# Derived directly so re-sourcing self-heals stale exported values.  Override the
# task list for a subset run by exporting RAD_TASKLIST (env.sh never sets that,
# so it can't self-poison): e.g.  RAD_TASKLIST=$RAD_WORK/tasks_smoke.txt bash submit_all.sh
export MANIFEST="$RAD_REPO/Analysis/hpc/manifest_dsb1.csv"
export TASKLIST="${RAD_TASKLIST:-$RAD_WORK/tasks.txt}"

# SGE queue / resources (adjust to Argon + your group's allocation).
export RAD_QUEUE="${RAD_QUEUE:-UI}"          # e.g. UI (free) or an investor queue
export RAD_MEM="${RAD_MEM:-8G}"              # per-task memory request
export RAD_HRT="${RAD_HRT:-2:00:00}"         # per-task hard runtime

# --- ROOT environment ------------------------------------------------------
# Make `root` available in batch jobs.  Pick ONE of the options below and fill
# it in after running the probe in README.md "Step 0".  This MUST work in a
# non-interactive shell (batch jobs do not source your ~/.bashrc by default).
setup_root() {
    # Argon: `module load root` (also loads python/3.9.9 for makeReport.py).
    # In a non-interactive SGE batch shell the `module` function is often not
    # defined yet, so initialise the module system first if needed.
    if ! command -v module >/dev/null 2>&1; then
        for init in /etc/profile.d/lmod.sh /etc/profile.d/z00_lmod.sh \
                    /etc/profile.d/00-modulepath.sh /etc/profile.d/modules.sh \
                    /opt/apps/lmod/lmod/init/bash; do
            [ -f "$init" ] && source "$init" && break
        done
    fi
    module load root
}

# Two Layer-1 macros (averageWaveforms.C, drs4TimeBase.C) read RAW waveforms and
# locate files via kRuns[].inFiles, i.e. Data/RUN<n>_<E>_GeV.root.  Build a Data/
# symlink farm so those names resolve to the HPC rec files.  If a specific kRuns
# run is missing, fall back to the first available run of that energy.
setup_data_links() {
    mkdir -p "$RAD_REPO/Data"
    grep -oE 'Data/RUN[0-9]+_[0-9]+_GeV\.root' "$RAD_REPO/Analysis/ChannelConfig.h" | sort -u | \
    while read -r name; do
        local run lab src alt
        run=$(printf '%s' "$name" | sed -E 's#Data/RUN([0-9]+)_.*#\1#')
        lab=$(printf '%s' "$name" | sed -E 's#Data/RUN[0-9]+_([0-9]+)_GeV\.root#\1GeV#')
        src="$REC_DIR/RUN${run}.root"
        if [ ! -f "$src" ] && [ -s "$TASKLIST" ]; then
            alt=$(awk -F'\t' -v l="$lab" '$2==l{print $4; exit}' "$TASKLIST")
            [ -n "$alt" ] && src="$alt"
        fi
        if [ -f "$src" ]; then
            ln -sfn "$src" "$RAD_REPO/$name"
        else
            echo "WARN: no raw file for $name — Layer-1 macros will skip $lab" >&2
        fi
    done
    echo "Data/ raw-file symlinks:"; ls -l "$RAD_REPO/Data" 2>/dev/null | grep -c '\->' | xargs echo "  links:"
}

# Compile-cache dir so concurrent ACLiC builds never collide and the cache
# persists between stages.
export RAD_ACLIC="${RAD_ACLIC:-$RAD_WORK/aclic_cache}"

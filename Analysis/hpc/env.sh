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

# Raw reconstructed run files: <REC_DIR>/RUN<run>.root
export REC_DIR="${REC_DIR:-/Shared/lss_yonel/jwwetzel/RADiCAL_CERN_May2023/rec/rec}"

# Work/scratch area for all outputs (per-run ntuples, merged ntuples, plots,
# report).  Should be on fast, large storage — NOT your home dir.  Argon: use
# your group's /Shared/lss_yonel space or local node scratch staged back.
export RAD_WORK="${RAD_WORK:-/Shared/lss_yonel/jwwetzel/RADiCAL_CERN_May2023/analysis_work}"

# Manifest (committed in repo) and the resolved per-run task list (made on Argon).
# Derived directly so re-sourcing self-heals stale exported values.
export MANIFEST="$RAD_REPO/Analysis/hpc/manifest_dsb1.csv"
export TASKLIST="$RAD_WORK/tasks.txt"

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

# Compile-cache dir so concurrent ACLiC builds never collide and the cache
# persists between stages.
export RAD_ACLIC="${RAD_ACLIC:-$RAD_WORK/aclic_cache}"

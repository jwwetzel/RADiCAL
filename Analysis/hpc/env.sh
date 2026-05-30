#!/bin/bash
# ============================================================================
# env.sh — shared configuration + ROOT environment for the Argon HPC pipeline.
# Every job script sources this.  EDIT the paths and the ROOT block for your
# Argon account (see README.md "Step 0: find ROOT").
# ============================================================================

# --- paths -----------------------------------------------------------------
# Repository root (where the Analysis/ directory lives) on Argon.
export RAD_REPO="${RAD_REPO:-$HOME/RADiCAL}"

# Raw reconstructed run files: <REC_DIR>/RUN<run>.root
export REC_DIR="${REC_DIR:-/Shared/lss_yonel/jwwetzel/RADiCAL_CERN_May2023/rec/rec}"

# Work/scratch area for all outputs (per-run ntuples, merged ntuples, plots,
# report).  Should be on fast, large storage — NOT your home dir.  Argon: use
# your group's /Shared/lss_yonel space or local node scratch staged back.
export RAD_WORK="${RAD_WORK:-/Shared/lss_yonel/jwwetzel/RADiCAL_CERN_May2023/analysis_work}"

# Manifest (committed in repo) and the resolved per-run task list (made on Argon).
export MANIFEST="${MANIFEST:-$RAD_REPO/Analysis/hpc/manifest_dsb1.csv}"
export TASKLIST="${TASKLIST:-$RAD_WORK/tasks.txt}"

# SGE queue / resources (adjust to Argon + your group's allocation).
export RAD_QUEUE="${RAD_QUEUE:-UI}"          # e.g. UI (free) or an investor queue
export RAD_MEM="${RAD_MEM:-8G}"              # per-task memory request
export RAD_HRT="${RAD_HRT:-2:00:00}"         # per-task hard runtime

# --- ROOT environment ------------------------------------------------------
# Make `root` available in batch jobs.  Pick ONE of the options below and fill
# it in after running the probe in README.md "Step 0".  This MUST work in a
# non-interactive shell (batch jobs do not source your ~/.bashrc by default).
setup_root() {
    # (A) Environment module:
    #     module load root            # or e.g. `module load stack/root 6.x`
    #
    # (B) CVMFS / LCG view (if /cvmfs is mounted on the compute nodes):
    #     source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc13-opt/setup.sh
    #
    # (C) Container (Singularity/Apptainer): handled in the job script via
    #     `apptainer exec <image> ...` — leave this function empty in that case.
    #
    # (D) A local ROOT install:
    #     source /path/to/root/bin/thisroot.sh
    : # <-- REPLACE this line with your chosen option above
}

# Compile-cache dir so concurrent ACLiC builds never collide and the cache
# persists between stages.
export RAD_ACLIC="${RAD_ACLIC:-$RAD_WORK/aclic_cache}"

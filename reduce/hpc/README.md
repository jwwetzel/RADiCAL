# Running the full DSB1 analysis on Argon (UIowa)

Scales the local RADiCAL analysis to the **full DSB1 electron dataset** (all 6
energies, 252 runs at the nominal 42.25 V bias) on the Argon HPC (SGE/`qsub`).

```
 logbook.csv ──build_manifest.py──▶ manifest_dsb1.csv   (252 runs: run,energy,label,bias,…)
 manifest + REC_DIR ──discover_tasklist.sh──▶ tasks.txt  (one line per existing RUN*.root)

 SGE pipeline (submit_all.sh chains the dependencies):
   compile.sh        pre-build processRun.so, point output → scratch     [1 job]
        │ -hold_jid
   sge_process.sh    processRun per run  → Output/byrun/<E>/RUN<run>/ntuple.root  [252 array tasks]
        │ -hold_jid
   sge_merge.sh      hadd per energy     → Output/<E>/ntuple.root                 [1 job]
        │ -hold_jid
   sge_analysis.sh   SKIP_PROCESS=1 runAll.sh  → Output/Summary/, report.html     [1 job]
```

Only `processRun` is heavy and it is **embarrassingly parallel** (one run per
task). Merge + analysis are cheap. The validated physics code is unchanged.

---

## Step 0 — find ROOT on Argon (do this first)

Batch jobs run in a non-interactive shell that does **not** read your `~/.bashrc`,
so ROOT must be set up explicitly. On a login node, probe what's available:

```bash
module avail 2>&1 | grep -i root          # is there a root module?
ls /cvmfs/sft.cern.ch/lcg/views 2>/dev/null   # is CVMFS/LCG mounted?
which root thisroot.sh 2>/dev/null         # already on PATH?
which apptainer singularity 2>/dev/null    # container runtime?
```

Then edit **`env.sh`** → `setup_root()` and uncomment the matching option, e.g.
`module load root`, or `source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc13-opt/setup.sh`.
Verify in a fresh shell: `bash -lc 'source reduce/hpc/env.sh; setup_root; root-config --version'`.

> If ROOT only comes from a container, tell me and I'll switch the job scripts
> to `apptainer exec <image> root …`.

---

## Step 1 — put the repo + paths on Argon

```bash
# from your Argon home (or wherever you keep code):
git clone <your repo>  RADiCAL          # or rsync your local repo up
cd RADiCAL
```

Edit `reduce/hpc/env.sh`:
- `REC_DIR`  — raw files (default: `/Shared/lss_yonel/jwwetzel/RADiCAL_CERN_May2023/rec/rec`)
- `RAD_WORK` — outputs; use roomy SHARED scratch (Argon `/nfsscratch`, default), **not** the
  group LSS (it's small/full — raw data lives there). `/nfsscratch` is auto-purged, so copy the
  final `report.html` + `Summary/` + merged ntuples to permanent storage when done.
- `RAD_QUEUE`, `RAD_MEM`, `RAD_HRT` — your group's queue/limits
- `setup_root()` — from Step 0

## Step 2 — manifest (already built; rebuild only if needed)

`manifest_dsb1.csv` is committed. To regenerate (e.g. include the 125 GeV bias
scan, or a different config):

```bash
python3 reduce/hpc/build_manifest.py /path/to/logbook.csv -o reduce/hpc/manifest_dsb1.csv
#   --all-bias          keep every bias (default: 42.25 V only)
#   --capillary LUAG    a different config (NB: needs its own channel map first)
```

## Step 3 — resolve run numbers → files (login node, fast)

```bash
bash reduce/hpc/discover_tasklist.sh
```
Globs `REC_DIR`, handles zero-padding (`RUN1034.root`, `RUN0100.root`, …), writes
`$RAD_WORK/tasks.txt`, and reports any manifest run with no file.

## Step 3.5 — smoke test (recommended before the full run)

Run the whole pipeline on a few runs per energy, end-to-end, on a compute node.
This catches the things a subset uniquely reveals: that `processRun` works on the
real `rec` files (raw `pulse`-tree schema matches), that `hadd` merges, and that
every analysis macro compiles under Argon's ROOT + `makeReport` runs — before you
commit 221 array tasks.

```bash
qlogin                                       # interactive compute node (NOT the login node)
cd <repo>
bash reduce/hpc/smoke_test.sh 2            # 2 runs/energy; writes to Output_smoke/
```
Inspect `$RAD_WORK/Output_smoke/report.html` + `Summary/`. If it's clean, proceed.
(The script uses a separate output dir and restores the `Output` symlink at the end.)

## Step 4 — submit

```bash
bash reduce/hpc/submit_all.sh        # chains compile → process[array] → merge → analysis
qstat -u "$USER"                       # monitor
```

Outputs land in `$RAD_WORK/Output/`: `report.html`, `Summary/` (hero figures,
`results.json`), and per-energy PDFs — identical layout to the local run.

---

## Notes & gotchas

- **Data volume.** 252 runs × ~30k events × ~2 GB raw ≈ 0.5 TB read. Per-run
  `processRun` is ~1–5 min, so the array clears in minutes-to-an-hour on a busy
  queue. Merged per-energy ntuples are small (the compact `rad` tree).
- **Bias.** The manifest uses only 42.25 V so every energy shares one SiPM gain;
  the 125 GeV 41.25 V scan (runs 100, 1001–1033) is excluded on purpose.
- **`-l` resource names.** `mem_free`/`h_rt` are common SGE resources but Argon
  may name them differently or allocate whole nodes — blank `RAD_MEM`/`RAD_HRT`
  in `env.sh` (or edit `submit_all.sh`) if `qsub` rejects them.
- **Report step** needs `python3` + `pdftoppm` (poppler). Without them the job
  still produces every PDF + `results.json`; regenerate `report.html` elsewhere.
- **Re-running one stage.** Each `sge_*.sh` can be `qsub`'d alone once its inputs
  exist (e.g. re-merge: `qsub -V -cwd -q $RAD_QUEUE reduce/hpc/sge_merge.sh`).
- **The other 3 configs** (LUAG, 2×DSB1+2×LuAG, 3×DSB1+1×Energy) are *separate
  analyses* — each needs its own channel map in `ChannelConfig.h` before its
  data is meaningful. The manifest builder already supports `--capillary` to
  list their runs when you're ready.
```

# Re-reduce ALL builds to the canonical schema (Argon)

The HPC reduce pipeline now produces the ONE canonical `rad::RadEvent` schema for
every build via the config-driven `radcore` Reducer. This replaces the old
two-path split (`processRun.C` for DSB1 + `reduceRaw.C` for the configs) with a
single path, and crucially gives **proper role-resolved `hg_cfd05`** (MCP-
referenced, per the build's channel map) for **all** builds — fixing the
low-energy collapse that the legacy slot-derived `cfd05` had for LUAG/MIXED/TENERGY.

The reducer is **bit-identical to `processRun.C`** (`reduce/validateReduce.C`
gate, worst diff = 0), so re-reducing DSB1 reproduces the published numbers
exactly; the others gain a rigorous `hg_cfd05` for the first time.

## Prerequisites
- Repo cloned on Argon, on `main`: `git pull`
- `REC_DIR` in `reduce/hpc/env.sh` points at the raw rec files
- ROOT available in batch (`setup_root` → `module load root`)

## 1. One-time compile (prebuilds `reduce/reduceRun_C.so`)
```
qsub reduce/hpc/compile.sh        # wait for rad_compile to finish (qstat)
```

## 2. Reduce each build (resolve → reduce[array] → merge per energy)
Run from the repo root on the **login node**. Note the manifest filename casing
(`manifest_dsb1.csv` is lowercase; the others are upper):
```
RAD_CONFIG=DSB1    RAD_MANIFEST=$PWD/reduce/hpc/manifest_dsb1.csv    bash reduce/hpc/submit_reduce.sh
RAD_CONFIG=LUAG    RAD_MANIFEST=$PWD/reduce/hpc/manifest_LUAG.csv    bash reduce/hpc/submit_reduce.sh
RAD_CONFIG=MIXED   RAD_MANIFEST=$PWD/reduce/hpc/manifest_MIXED.csv   bash reduce/hpc/submit_reduce.sh
RAD_CONFIG=TENERGY RAD_MANIFEST=$PWD/reduce/hpc/manifest_TENERGY.csv bash reduce/hpc/submit_reduce.sh
```
`RAD_CONFIG` must match the config JSON name (`data/2023/configs/<RAD_CONFIG>.json`).
Each call submits a reduce array job + a dependent merge. Monitor: `qstat -u $USER`.

**Output (canonical):** `data/2023/reduced/<BUILD>/<E>GeV.root`

The canonical files are a bit larger than the old reduced ones (they carry the
full per-cap feature set + the generic `s_peak/s_cfd05/s_charge[36]` arrays). The
default `RAD_MEM=8G`, `RAD_HRT=2:00:00` per run-task are ample.

## 3. Validate on Argon (or after copying home)
```
ROOT_INCLUDE_PATH=lib/waveform:lib/io:lib/physics:lib/viz root -l -b -q 'analyze/sigmaT.C+("data/2023/configs/DSB1.json",150)'
ROOT_INCLUDE_PATH=lib/waveform:lib/io:lib/physics:lib/viz root -l -b -q analyze/timingLadder.C+
```
- DSB1 must still reproduce the published ladder (47.1 … 27.4 ps).
- LUAG/MIXED/TENERGY now give physical `sigma_t` at **every** energy (no more
  sub-floor collapse) — the cross-build ladder becomes paper-grade.

## 4. Copy home (if analysing on your Mac)
```
rsync -av <argon>:<repo>/data/2023/reduced/  ./data/2023/reduced/
```

---
Once validated, the old DSB1-only `sge_process.sh` (processRun) path is redundant
and can be retired (P6). For a new year, set `RAD_YEAR=2024`, drop the raw files
under `data/2024/raw/`, author `data/2024/configs/<BUILD>.json`, and run
the same flow.

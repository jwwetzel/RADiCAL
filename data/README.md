# RADiCAL datasets — multi-year, one repository

This directory holds **one subdirectory per test-beam campaign** (`2023/`, `2024/`,
`2025/`, `2026/`, …). The detector and DAQ are largely the same across years; what
differs is mostly the **channel map / build configs** and a few per-dataset tweaks.
Those differences live here as *configuration*; the analysis *code* stays shared in the
repo-root `lib/`, `analyze/`, and `reduce/`. The organizing principle and the full
new-dataset / new-paper playbooks are in **`docs/PUBLICATION_ENGINE.md`**.

## Why one repo (not one per year/run/report)

- **Shared framework, fixed once.** One `lib/` + `analyze/` + `reduce/` serves every
  year; a fix or method improvement propagates everywhere instead of being copy-pasted.
- **Reduction is config-agnostic.** `reduce/reduceRaw.C` stores all 36 DRS slots
  `(peak, CFD-5% time, charge)` plus the wire-chamber track and both MCPs — *no channel
  map needed at reduction time*. The map only matters when you interpret slots as
  capillaries (analysis), so a new year is a new **config**, not a code fork.
- **Cross-year synthesis is trivial here.** Comparing campaigns apples-to-apples wants a
  single tree; the analysis reads any campaign via `RAD_YEAR`.
- **Data size is irrelevant to topology.** Raw/reduced data is gitignored.

## Layout

```
data/
  <year>/
    configs/             # the ONLY files the code reads (tracked)
      <BUILD>.json       #   per-build: DRS4 [drs,grp,ch] map, materials/roles, geometry, run list
      <BUILD>.hglg       #   per-build HG/LG calibration sidecar (optional)
    metadata/            # human logbook, decorative — not read by code (tracked)
      dataset.yaml       #   campaign metadata + any per-dataset overrides
      channel_map.yaml   #   slot -> capillary wiring, mirror of the JSON
      runs.csv           #   run manifest / logbook (run, energy, build, bias, status)
    raw/                 # raw 'pulse' .root files, ~2 GB/run                (gitignored)
    reduced/<BUILD>/     # reduceRaw.C output 'rad' tree, ~12 MB/run         (gitignored)
  _TEMPLATE/             # copy this to start a new campaign
```

`configs/` and `metadata/` are tracked; `raw/` and `reduced/` are gitignored payloads.

## Switching campaigns — `RAD_YEAR`

Every macro resolves data through `lib/io/DataPaths.h` (the single resolver — no legacy
fallbacks). Two env vars parameterize it:

- **`RAD_YEAR`** picks the campaign: `export RAD_YEAR=2024` makes every `radReduced()` /
  `radConfig()` / `radRaw()` point at `data/2024/` — **no recompile**. Default is 2023.
- **`RAD_DATA`** overrides the base directory (defaults to the repo root) if your data
  payloads live elsewhere.

## Add a new campaign (summary; full steps in `docs/PUBLICATION_ENGINE.md`)

```
cp -r data/_TEMPLATE data/2027
# 1. fill data/2027/metadata/{dataset.yaml,channel_map.yaml,runs.csv}  (the logbook)
# 2. write data/2027/configs/<BUILD>.json per build  (copy a 2023 one;
#    VERIFY every [drs,grp,ch] vs this year's wiring + hg_sat_mV vs its DC offset)
# 3. drop a few sample runs into data/2027/raw/ and smoke-test:
#    RAD_YEAR=2027 root -l -b -q 'analyze/sigmaT.C+("DSB1",150)'
# 4. reduce on the cluster, then run the standard analysis with RAD_YEAR=2027
```

The per-build JSON (`BuildConfig`) is authoritative: it fully describes a build (DRS4
`[drs,grp,ch]` map, per-corner material/role, geometry, run list). The universal flat
index math lives in `lib/waveform/ChannelConfig.h`; the per-campaign wiring is what you
verify in each `configs/<BUILD>.json`.

## Cluster (HPC) workflow

Build the task list from `runs.csv`, `reduceRaw.C` each raw run into `reduced/`, merge,
then run the shared analysis. `reduce/hpc/` already parameterizes by manifest and respects
`RAD_YEAR` — point it at the campaign's runs. See `reduce/hpc/PRESCRIPTION.md`.

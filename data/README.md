# RADiCAL datasets — multi-year, one repository

This directory holds **one subdirectory per test-beam campaign** (`2023/`, `2024/`,
`2025/`, `2026/`, …). The experimental setup and DAQ are largely the same across
years; what differs between campaigns is mostly the **channel map** and a few
per-dataset tweaks. Those differences live here as *configuration*; the analysis
*code* stays shared in the repo-root `Analysis/`.

## Why one repo (not one per year/run/report)

- **Shared framework, fixed once.** One `Analysis/` serves every year; a fix or a
  method improvement propagates everywhere instead of being copy-pasted N times.
- **Reduction is already config-agnostic.** `Analysis/reduceRaw.C` stores all 36
  DRS slots `(peak, CFD-5% time, charge)` plus the wire-chamber track and both
  MCPs — *no channel map needed at reduction time*. The map only matters when you
  interpret slots as capillaries (analysis), so a new year is a new **config**,
  not a code fork.
- **Cross-year synthesis is trivial here, painful across many repos.** Comparing
  campaigns apples-to-apples (as the 2023 capillary study did within one year)
  wants a single tree. The synthesis report draws from every `data/<year>/`.
- **Data size is irrelevant to topology.** Raw/reduced data is gitignored and
  never enters git either way (see `.gitignore`).
- **Per-run / per-report isolation** = a subdirectory, a branch, or a git tag —
  not a separate repository.

## Layout

```
data/
  <year>/
    config/
      channel_map.yaml   # slot -> capillary wiring for THIS campaign (the thing that differs)
      dataset.yaml       # campaign metadata + any per-dataset analysis overrides
      runs.csv           # run manifest / logbook (run, energy, config, bias, status, ...)
    raw/                 # drop sample raw .root files here              (gitignored)
    reduced/             # reduceRaw.C output, ~12 MB/run                (gitignored)
    output/              # analysis ntuples / scratch                   (gitignored)
    figures/             # committed key figures for this campaign
  _TEMPLATE/             # copy this to start a new campaign
```

Only `config/` and `figures/` are tracked in git; `raw/`, `reduced/`, `output/`
are data payloads and are ignored (each keeps a `.gitkeep` so the empty dir
persists).

## Add a new campaign

```
cp -r data/_TEMPLATE data/2027
# 1. edit data/2027/config/dataset.yaml      (year, energies, beam, overrides)
# 2. edit data/2027/config/channel_map.yaml  (VERIFY every (drs,grp,ch) vs this year's wiring)
# 3. fill data/2027/config/runs.csv          (the run logbook)
# 4. drop a few sample runs into data/2027/raw/ to smoke-test before the cluster
```

The shared framework reads the campaign by name (planned: a `--dataset <year>`
argument / `RAD_DATASET` env that points the macros at `data/<year>/config/`;
today the 2023 map is compiled into `Analysis/ChannelConfig.h` — see *Status* below).

## Cluster (HPC) workflow

Per campaign: build the task list from `runs.csv`, `reduceRaw.C` each raw run on
the cluster into `reduced/`, merge, then run the shared analysis. The existing
`reduce/hpc/` scripts already parameterize by manifest — point them at the
campaign's `runs.csv`. (See `reduce/hpc/PRESCRIPTION.md`.)

## Data centralization (2023)

All 2023 data now resolves through **one** place — `Analysis/DataPaths.h`:

- **Canonical layout:** `data/2023/{raw, reduced/<BUILD>}` (mirrored on CERNBox).
- **`data/2023/MANIFEST.csv`** maps every file to its canonical and legacy path.
- **`Analysis/organize_data.sh`** consolidates local data into the canonical tree
  (link / copy / move), or — for a newcomer who downloaded `data/2023/` from
  CERNBox — bridges the old in-repo paths with `--legacy-links`.
- The resolver tries the canonical path first and **falls back** to the legacy
  locations (`Data/`, `reduced/`, `output/`), so existing checkouts and the
  live site keep working during the transition.

See `data/2023/README.md` to get started. The 2023 channel map is still compiled
into `Analysis/ChannelConfig.h` (mirrored in `config/channel_map.yaml`); teaching the
framework to load `channel_map.yaml` per campaign is the remaining multi-year step.

**Code-migration status:** the data-reading macros are being moved onto
`DataPaths.h` (`radRaw()` / `radReduced()`). Until all are migrated, the
`--legacy-links` bridge makes every macro work against a CERNBox download.

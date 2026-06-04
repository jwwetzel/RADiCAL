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
  wants a single tree. The synthesis report draws from every `datasets/<year>/`.
- **Data size is irrelevant to topology.** Raw/reduced data is gitignored and
  never enters git either way (see `.gitignore`).
- **Per-run / per-report isolation** = a subdirectory, a branch, or a git tag —
  not a separate repository.

## Layout

```
datasets/
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
cp -r datasets/_TEMPLATE datasets/2027
# 1. edit datasets/2027/config/dataset.yaml      (year, energies, beam, overrides)
# 2. edit datasets/2027/config/channel_map.yaml  (VERIFY every (drs,grp,ch) vs this year's wiring)
# 3. fill datasets/2027/config/runs.csv          (the run logbook)
# 4. drop a few sample runs into datasets/2027/raw/ to smoke-test before the cluster
```

The shared framework reads the campaign by name (planned: a `--dataset <year>`
argument / `RAD_DATASET` env that points the macros at `datasets/<year>/config/`;
today the 2023 map is compiled into `Analysis/ChannelConfig.h` — see *Status* below).

## Cluster (HPC) workflow

Per campaign: build the task list from `runs.csv`, `reduceRaw.C` each raw run on
the cluster into `reduced/`, merge, then run the shared analysis. The existing
`Analysis/hpc/` scripts already parameterize by manifest — point them at the
campaign's `runs.csv`. (See `Analysis/hpc/PRESCRIPTION.md`.)

## Status / migration note

The **2023 campaign is the live reference**: its analysis, the six-layer
`report/`, and `paper.html` currently sit at the repo root and on GitHub Pages,
and its channel map is compiled into `Analysis/ChannelConfig.h`. To avoid breaking
the live site, 2023 has **not** been physically moved under `datasets/2023/` yet —
that directory documents the 2023 config and points at the legacy locations. The
clean next step is to (a) teach the framework to load `channel_map.yaml` per
campaign, and (b) optionally migrate 2023's data/config under `datasets/2023/`.
New campaigns (2024+) start clean under this structure from day one.

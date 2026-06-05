# RADiCAL 2023 data — get started

Everything you need to analyse the CERN SPS H2 (May 2023) data lives under this one
folder. The files are large, so they are **not** in git — download them from CERNBox
(below) into this directory tree.

## Layout (mirror this on CERNBox)

```
datasets/2023/
  raw/                       # raw DRS4 waveforms ('pulse' tree), ~2 GB/run  [optional]
    RUN1211_25_GeV.root … RUN1261_150_GeV.root          (DSB1 runs)
  reduced/                   # small, analysis-ready ntuples ('rad' tree)    [START HERE]
    DSB1/    25GeV.root … 150GeV.root                   (rad-processRun format)
    LUAG/    50GeV.root … 150GeV.root                   (rad-reduceRaw format)
    MIXED/   50GeV.root … 150GeV.root
    TENERGY/ 50GeV.root … 150GeV.root
  config/    channel_map.yaml · dataset.yaml · runs.csv
  MANIFEST.csv               # every file, its canonical path, and its legacy path
```

The **reduced/** files (~12 MB each) are all you need for the full timing / energy /
position analysis. The **raw/** files are only needed for waveform-level studies or
re-reduction.

## Download from CERNBox

> **CERNBox:** <PASTE_CERNBOX_SHARE_LINK_HERE>  *(public share — folder structure
> matches the tree above)*

Download (at minimum) `reduced/` into `datasets/2023/reduced/` here. Then either:

**A. Point the code at it (recommended).** The analysis resolves data through
`Analysis/DataPaths.h`; by default the base is the repo root, so files placed at
`datasets/2023/reduced/<BUILD>/<E>GeV.root` are found automatically. (Set
`export RAD_DATA=/path/to/checkout` if your data lives elsewhere.)

**B. Bridge the legacy paths.** Some macros still use the old in-repo locations
(`reduced/`, `Analysis/Output/`). One command makes them all work against your
download, without moving anything:

```
./Analysis/organize_data.sh --legacy-links
```

## Two ntuple formats (both have a `rad` tree)

- **DSB1** (`reduced/DSB1/*`): the `processRun.C` output — per-capillary branches
  `hg_cfd05[8]`, `lg_peak[8]`, `sum_lg`, `wc_ok`, `x_trk`, `y_trk`, `mcp_peak`, …
- **LuAG / MIXED / TENERGY** (`reduced/<BUILD>/*`): the config-agnostic `reduceRaw.C`
  output — all 36 slots in `s_peak[36]`, `s_cfd05[36]`, `s_charge[36]`, plus
  `mcp1_time`, `mcp2_time`, `mcp1_peak`, the WC track, and `run`.

The channel map (`config/channel_map.yaml`) tells you which slot is which capillary.

## First analysis

Try the cross-build comparison once `reduced/` is in place:

```
ROOT_INCLUDE_PATH=Analysis root -l 'Analysis/configResolutionFull.C+'   # σ_t(E), σ_E(E) per build
ROOT_INCLUDE_PATH=Analysis root -l 'Analysis/fourDdemo.C+'              # time+energy+position
```

For a guided, hands-on first analysis on the raw data, see the
[Data Lab](../../data_lab.html).

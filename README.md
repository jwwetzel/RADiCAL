# RADiCAL — CERN May 2023 Test-Beam Analysis

Analysis of the RADiCAL radiation-hard shashlik calorimeter prototype, taken in
the CERN SPS H2 electron beam (25–150 GeV) in May 2023. The pipeline goes from
raw CAEN DT5742 (DRS4) waveforms to energy- and timing-resolution results and a
single browsable HTML report.

**Headline result:** σ_t ≈ 37 ps at 150 GeV (MCP-free (DW−UP)/2 estimator).

---

## Quick Start

```bash
# 1. Install prerequisites (Homebrew, ROOT, optional Python stack)
bash Analysis/setup.sh
#    then, in the current shell:
source "$(brew --prefix root)/bin/thisroot.sh"

# 2. Get the raw data into Data/  (~18 GB, NOT in this repo — see "Data" below)
#    Data/RUN1211_25_GeV.root, ... RUN1258-1261_150_GeV.root

# 3. Run the full chain (process all 6 energies → all plots → report)
bash Analysis/runAll.sh

# 4. View the report
open Analysis/Output/report.html
```

Run a single step instead of the whole chain, e.g. one energy:
```bash
ROOT_INCLUDE_PATH=Analysis root -l -b -q \
  'Analysis/processRun.C+("Data/RUN1211_25_GeV.root", 25., "25GeV")'
ROOT_INCLUDE_PATH=Analysis root -l -b -q 'Analysis/timingEnergyBins.C+'
```
> Tip: macros `#include` the shared headers in `Analysis/`, so set
> `ROOT_INCLUDE_PATH=Analysis` (or run from a shell where `runAll.sh` has set it).
> The report's PNG conversion needs `pdftoppm` (poppler) or `gs` (ghostscript).

---

## Repository map

```
Analysis/            ← the codebase (everything below runs from here)
  *.h                shared headers (config, cuts, utils, style)
  *.C                analysis macros (one concern each)
  runAll.sh          one-command pipeline (process → analyse → report)
  setup.sh           install ROOT + deps (macOS/Homebrew + cluster notes)
  makeReport.py      assemble Analysis/Output/report.html
  GAP_CLOSING_PLAN.md analysis roadmap / status
  Output/            generated results (gitignored; created by runAll)
Data/                raw DT5742 .root files (gitignored; fetch externally)
legacy/              superseded first-pass code (do not use — see legacy/README)
README.md            this file
```

### Shared headers (the foundation)
| Header | Role |
|--------|------|
| `ChannelConfig.h` | detector geometry, channel map, **run list `kRuns[]`**, capillaries `kCap[]` |
| `SelectionCuts.h` | **single source of truth** for every cut threshold |
| `WaveformUtils.h` | DRS4 pulse extraction (pedestal, peak, CFD, charge) |
| `PlotUtils.h` | Gaussian/Crystal-Ball fits, `ScanRunCenters`, plot helpers |
| `RADiCALStyle.h` | unified ROOT plotting style |
| `DRS4Calibration.h` | DRS4 stop-cell recovery + timing correction |

### Pipeline macros (run in this order by `runAll.sh`)
| Macro | Produces |
|-------|----------|
| `processRun.C` | raw waveforms → compact per-energy ntuple (`Output/<E>/ntuple.root`) |
| `analyzeResolution.C` | energy resolution, linearity, σ_t-vs-E summary |
| `timingResolution.C` | channel-combination timing (best-single / corners / A²-weighted) |
| `timingMethods.C` | CFD-fraction / LED / walk-correction method comparison |
| `qualityPlots.C` | per-energy quality report + cross-energy quality summary |
| `timingEnergyBins.C` | **headline** energy-binned (DW−UP)/2 σ_t |
| `compareEnergies.C` | cross-energy collated plots (beam, containment, channels, timing) |
| `investigatePbGlass.C` | EM / hadronic / halo / off-axis event classification |
| `mcpJitter.C` | MCP reference jitter + quadrature subtraction |
| `timingContainmentScan.C` | containment-cut optimisation (σ_t vs yield) |
| `channelCombinationScan.C` | brute-force best-N-channel combination scan |
| `positionCorrection.C` | charge-sharing position reconstruction + timing correction |
| `channelCalibration.C` | inter-channel relative gain calibration |
| `drs4TimeBase.C` | DRS4 time-base verification + stop-cell correction validation |
| `drs4Diagnostics.C` | DRS4 hardware diagnostics (noise, saturation, spikes) |
| `channelIntegrity.C` | per-channel integrity (active fraction, HG/LG, cross-talk) |
| `uniformityScan.C` | spatial σ_t(x,y) uniformity |
| `systematicUncertainties.C` | cut-variation + fit-model systematic budget |
| `makeReport.py` | converts all PDFs → PNGs and builds `report.html` |

### Adding things (scalability)
- **New beam energy** → add one line to `kRuns[]` in `ChannelConfig.h`.
- **Change a cut** → edit `SelectionCuts.h` only (no magic numbers elsewhere).
- **New macro** → `#include "ChannelConfig.h"` + `"PlotUtils.h"`, call
  `ApplyRADiCALStyle()` then `ScanRunCenters()`, follow the existing template.

---

## Data

The raw files are **not** in the repository (~18 GB). Copy them into `Data/`:

| Location | Path |
|----------|------|
| Iowa Argon cluster | `/Shared/lss_yonel/jwwetzel/RADiCAL_CERN_May2023/rec/rec` |
| CERN lxplus (EOS) | `/eos/experiment/iucmsdata/Data/2023/H2TB-202305-RADiCAL/RADiCAL_2023May/rec` |

Expected files: `RUN1211_25_GeV.root`, `RUN1148_50_GeV.root`, `RUN1112_75_GeV.root`,
`RUN1075_100_GeV.root`, `RUN1034_125_GeV.root`, and `RUN1258–1261_150_GeV.root`
(four files chained for 150 GeV).

---

## Detector & DAQ (background)

Two CAEN **DT5742** digitisers (16+1 channels each) record waveforms via the
**DRS4** switched-capacitor chip. Each DT5742 stores its channels in two groups
of 8 signal + 1 trigger, so the reconstructed data has channels 0–8 (8 = trigger)
and 9–17 (17 = trigger) per board — 18 "channels" because the trigger is
duplicated. Each channel saves 1024 time samples per event.

The detector is a compact (14 mm) W/LYSO shashlik module read out by 8 capillary
fibres (4 corners × Up/Down): **HG** channels (shower-max, → timing) and **LG**
channels (full-length, → energy), plus MCP timing references and PbGlass leakage
counters. Full channel map: `Analysis/ChannelConfig.h`.

## Raw data format (reference)

Per event, all channels are stored end-to-end in two flat arrays in the `pulse`
tree: `amplitude` (waveforms) and `timevalue` (time axes). Indexing:

```cpp
// drs ∈ {0,1}, group ∈ {0,1}, ch ∈ 0..8
channel_index = (1024*9*2)*drs + (1024*9)*group + 1024*ch
time_index    = (1024*2)*drs   + 1024*group
```
i.e. channel 0's waveform is samples `[0,1023]`, channel 1 is `[1024,2047]`, etc.
These index helpers are wrapped as `chanOff()`/`timeOff()` in `ChannelConfig.h`.

Reading the tree:
```cpp
TTreeReader reader("pulse", drs_file);
TTreeReaderValue<int>   run(reader, "run"), event0(reader, "event0"),
                        event1(reader, "event1"), trigger(reader, "trigger");
TTreeReaderArray<float> timevalue(reader, "timevalue");
TTreeReaderArray<float> amplitude(reader, "amplitude");
```
`processRun.C` turns these raw waveforms into the compact analysis ntuple that
all other macros read.

---

## Provenance

Method builds on A. Ledovskoy's prescription:
<https://github.com/ledovsk/RADiCAL_TB_May2023>. The original single-file
analysis lives in `legacy/` (superseded — do not use).

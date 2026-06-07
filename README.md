# RADiCAL — CERN May 2023 Test-Beam Analysis

Analysis of the RADiCAL radiation-hard W/LYSO:Ce shashlik calorimeter prototype,
taken in the CERN SPS H2 electron beam (25–150 GeV) in May 2023. The pipeline
goes from raw CAEN DT5742 (DRS4) waveforms to energy- and timing-resolution
results and a browsable HTML report.

**Headline result:** σ_t ≈ 27 ps at 150 GeV (DSB1, MCP-free `(DW−UP)/2`
shower-depth estimator) — a ~1 X₀ shower-depth timing floor.

---

## Quick start (local analysis)

```bash
source setup.sh                      # sets ROOT_INCLUDE_PATH + RAD_DATA, once per shell

# verify the reduced data is present + clean (all builds Y/Y/Y)
root -l -b -q 'reduce/verify.C+("data/2023/reduced")'

# headline timing resolution for a build (reads data/2023/reduced/<BUILD>/)
root -l -b -q 'analyze/timingHeadline.C+("DSB1")'
```

That's it — all **reduced** data lives in `data/2023/reduced/` (in the repo,
multi-GB but local for fast analysis). You only need Argon to **re-reduce** from
raw (see *Reducing from raw* below).

---

## Repository map

```
lib/                 ← shared library (all headers). The single include root.
  waveform/   WaveformUtils ChannelConfig DRS4Calibration   (raw-pulse domain)
  io/         MiniJson BuildConfig Schema RadView DataPaths  (config + canonical IO)
  physics/    SelectionCuts RadTiming                        (cuts + timing primitives)
  viz/        PlotUtils RADiCALStyle                         (plotting)
reduce/              ← raw → reduced pipeline
  Reducer.C reduceRun.C calibHGLG.C validateReduce.C verify.C
  hpc/        Argon SGE pipeline (env.sh, submit_reduce.sh, ...)
analyze/             ← reduced → results
  timingHeadline.C sigmaT.C timingLadder.C slopeVsE.C        (curated headline analyses)
  makeReport.py deployReport.sh                              (report build/publish)
  studies/    ~80 one-off investigation macros (kept, runnable, quarantined)
data/                ← ALL data, per year
  2023/{raw, reduced, configs, metadata}/                    (configs = <BUILD>.json + .hglg)
figures/             ← analysis output plots (PNG/PDF)
site/                ← published GitHub Pages website (self-contained)
  *.html  report/  assets/                                   (deployed via Actions; see site/README.md)
student/             ← self-contained student lab module (.C macros + expected figs)
docs/                ← design notes, plans, paper drafts
archive/             ← legacy code (superseded; kept for provenance)
setup.sh             ← source once per shell
```

**Three rules that keep it DRY + navigable:**
1. **One include root** — `source setup.sh` puts `lib/{waveform,io,physics,viz}`
   on `ROOT_INCLUDE_PATH`; every macro `#include`s headers bare.
2. **One data resolver** — `lib/io/DataPaths.h`: `radReduced(build,E)`,
   `radRaw(name)`, `radConfig(build)` all resolve under `data/<year>/`.
3. **Dependency direction** — `reduce/` and `analyze/` depend on `lib/`, never
   the reverse; `lib/` knows nothing about specific studies.

### The library (`lib/`)
| Header | Role |
|--------|------|
| `io/Schema.h` | the canonical reduced-ntuple branch layout (`rad` tree) |
| `io/BuildConfig.h` | per-build JSON config loader + channel-map resolution (+ `.hglg` sidecar) |
| `io/RadView.h` | typed reader over a reduced file (`hg_peak`, `timeOf(c,src)`, beam center) |
| `io/DataPaths.h` | **the** path resolver — `radReduced/radRaw/radConfig/radHglg` |
| `io/MiniJson.h` | minimal JSON parser for the build configs |
| `waveform/WaveformUtils.h` | DRS4 pulse extraction (pedestal, peak, CFD, charge) |
| `waveform/ChannelConfig.h` | detector geometry + `chanOff()`/`timeOff()` raw indexing |
| `waveform/DRS4Calibration.h` | DRS4 stop-cell recovery + timing correction |
| `physics/SelectionCuts.h` | single source of truth for every cut threshold |
| `physics/RadTiming.h` | `timingBestBin()` — the headline (DW−UP)/2 σ_t primitive |
| `viz/PlotUtils.h` | fits, `ScanRunCenters`, plot helpers |
| `viz/RADiCALStyle.h` | unified ROOT plotting style |

---

## Builds & energies

Four detector configurations ("builds"), identical channel map:

| Build | Material | Energies (GeV) |
|-------|----------|----------------|
| `DSB1` | LYSO:Ce, high light (headline) | 25, 50, 75, 100, 125, 150 |
| `LUAG` | LuAG:Ce, lower light | 50–150 |
| `MIXED` | NE+SW = DSB1, NW+SE = LuAG | 50–150 (valid era 2975–3044) |
| `TENERGY` | 3×DSB1 + 1×Energy capillary | 50–150 |

Each build has a config at `data/2023/configs/<BUILD>.json` and an HG/LG gain
calibration sidecar `<BUILD>.hglg` (written by `reduce/calibHGLG.C`).

---

## Reducing from raw (Argon HPC)

Reduction (raw waveforms → reduced `rad` ntuple) runs on the Iowa **Argon**
cluster. See `reduce/hpc/README.md` and `reduce/hpc/REREDUCE.md` for the full
recipe; the short version, per build:

```bash
# on the Argon login node, from the repo root
qsub reduce/hpc/compile.sh                                  # prebuild reducer .so (once)
RAD_CONFIG=DSB1 RAD_MANIFEST=$PWD/reduce/hpc/manifest_DSB1.csv \
  bash reduce/hpc/submit_reduce.sh                          # reduce[array] -> merge
# then copy home:  data/2023/reduced/<BUILD>/<E>GeV.root
```

The reduced files store HG times already MCP-referenced, plus `hg_lgcfd`
(LG-predicted-true-peak CFD) and the TR0 trigger-scintillator times. Run
`reduce/verify.C` after a reduction to confirm every file carries the new
branches before shipping.

---

## Data layout & where to get raw

```
data/<year>/raw/        RUN<n>_<E>_GeV.root   (raw 'pulse' tree; subset kept locally)
data/<year>/reduced/    <BUILD>/<E>GeV.root   (analysis 'rad' tree; ALL kept locally)
data/<year>/configs/    <BUILD>.json + .hglg  (build configs)
data/<year>/metadata/   channel_map.yaml, dataset.yaml, runs.csv
```

Raw files (~18 GB/build) are gitignored. Sources:

| Location | Path |
|----------|------|
| Iowa Argon cluster | `/Shared/lss_yonel/jwwetzel/RADiCAL_CERN_May2023/rec/rec` |
| CERN lxplus (EOS) | `/eos/experiment/iucmsdata/Data/2023/H2TB-202305-RADiCAL/RADiCAL_2023May/rec` |

---

## Detector & DAQ (background)

Two CAEN **DT5742** digitisers (16+1 channels each) record waveforms via the
**DRS4** switched-capacitor chip. Each DT5742 stores its channels in two groups
of 8 signal + 1 trigger, so reconstructed data has channels 0–8 (8 = trigger) and
9–17 (17 = trigger) per board. Each channel saves 1024 time samples per event.

The module (14 mm, W/LYSO shashlik) is read out by 8 capillary fibres (4 corners ×
Up/Down): **HG** channels (shower-max, → timing, negative-going, clip ~820 mV) and
**LG** channels (full-length, → energy, unsaturated), plus 2 MCP timing references
and the TR0 trigger scintillator (ch 8 of each DRS0 group). Channel map:
`lib/waveform/ChannelConfig.h` and the per-build `data/2023/configs/<BUILD>.json`.

### Raw data format (reference)

Per event, all channels are stored end-to-end in two flat arrays in the `pulse`
tree: `amplitude` (waveforms) and `timevalue` (time axes):

```cpp
// drs ∈ {0,1}, group ∈ {0,1}, ch ∈ 0..8
channel_index = (1024*9*2)*drs + (1024*9)*group + 1024*ch
time_index    = (1024*2)*drs   + 1024*group
```

These index helpers are wrapped as `chanOff()`/`timeOff()` in `ChannelConfig.h`.
`reduce/reduceRun.C` (→ `reduce/Reducer.C`) turns these raw waveforms into the
compact analysis ntuple that everything in `analyze/` reads.

---

## Provenance

Method builds on A. Ledovskoy's prescription:
<https://github.com/ledovsk/RADiCAL_TB_May2023>. The original single-file analysis
and superseded first-pass code live in `archive/` (kept for provenance — not used).

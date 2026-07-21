# RADiCAL — CERN May 2023 Test-Beam Analysis

Analysis of the RADiCAL radiation-hard W/LYSO:Ce shashlik shower-maximum calorimeter
prototype, CERN SPS H2 electron beam (25–150 GeV), May 2023. The pipeline goes from raw
CAEN DT5742 (DRS4) waveforms to timing-, energy-, and position-resolution results and a
journal manuscript, with every published number traced to a committed, audited generator.

**Headline result (timing paper, `papers/timing/`):**
σ_t = **25.7 ± 0.6 ps** at 150 GeV for the brightest-1000 showers of the high-light DSB1
build (saturation-recovered CFD, reference-free (DW−UP)/2), **≈50 ps over the full fiducial
sample**; the fitted floor **18.8 ± 0.8 ps confirms** the previously published 17.5 ps
(NIM A 1068 (2024) 169737). Across wavelength-shifter builds the stochastic term tracks
detected light — a = 203±6 → 440±18 ps·√GeV from the high-light to the low-light build —
while the same-shower MIXED control (σ_DSB1/σ_LuAG = 1.04 ± 0.05) disfavors intrinsic WLS
re-emission kinetics as the dominant cause, within this geometry and at these light levels.

**The trust path (start here):**
1. [`ANALYSIS_GUIDE.md`](ANALYSIS_GUIDE.md) — the top-down waterfall: raw → reduced →
   selection → estimator → gated results → manuscript, with the figure/number provenance table.
2. [`papers/timing/radical_timing.pdf`](papers/timing/radical_timing.pdf) — the manuscript.
3. [`papers/scripts/`](papers/scripts/) — one directory per published result ("gate"), each
   with a pre-registered `AUDIT.md`, the generator macro, and its committed result log.
4. [`papers/memory_claims_and_forbidden_language.md`](papers/memory_claims_and_forbidden_language.md)
   — the wording law every public number obeys; and
   [`papers/CAMPAIGN_SNAPSHOT_2026-06-09.md`](papers/CAMPAIGN_SNAPSHOT_2026-06-09.md) — the
   append-only campaign record.
5. [`chain_of_evidence.html`](chain_of_evidence.html) — browsable raw-to-result evidence for
   all four builds.

---

## Quick start (local analysis)

```bash
source setup.sh                      # ROOT_INCLUDE_PATH + RAD_DATA, once per shell

# one-command reproduction of the headline chain, verdict included:
tools/repro.sh --check-data          # full protocol + env pinning: REPRODUCE.md

# or by hand:
root -l -b -q 'reduce/verify.C+("data/2023/reduced")'          # data present + clean
root -l -b -q 'papers/scripts/timing_fit_summary/timingFitSummary.C'   # money plot
```

⚠ Gate macros regenerate their committed logs/tables/figures in place — run them to *verify*
(outputs should be identical), not casually; `git diff` afterward is the check (known
benign residual: one PDF timestamp — see [`REPRODUCE.md`](REPRODUCE.md)). Per-gate runbook:
[`papers/scripts/INDEX.md`](papers/scripts/INDEX.md); guarded runs: `tools/run_gate.sh`.
Data integrity: `shasum -a 256 -c data/2023/MANIFEST.sha256`.

**Data:** the reduced ntuples (~13 GB) are **not in git**. Collaborators: mirror from Argon
with `tools/pull_argon_data.sh`, or copy `data/2023/reduced/` from an existing checkout.
A public archival deposit (DOI) is a pre-submission task tracked in
[`WORKSPACE_REORG_PLAN_2026-06-10.md`](WORKSPACE_REORG_PLAN_2026-06-10.md) §5. Raw waveforms
live on Argon/EOS (paths below); a small local subset under `data/2023/raw/` feeds the
waveform-level figures.

---

## Repository map

```
ANALYSIS_GUIDE.md    ← THE spine: waterfall + provenance table (read first)
lib/                 ← shared header-only library. The single include root.
  waveform/   WaveformUtils ChannelConfig DRS4Calibration    (raw-pulse domain)
  io/         MiniJson BuildConfig Schema RadView DataPaths FigPaths (config + IO)
  physics/    SelectionCuts RadTiming                        (cuts + timing primitives)
  viz/        PlotUtils RADiCALStyle                         (plotting)
reduce/              ← raw → reduced pipeline (Reducer.C; hpc/ = Argon SGE recipe → REREDUCE.md)
analyze/             ← curated exploratory analyses + report builder (README inside)
  studies/    ~120 investigation macros (exploratory record, indexed in studies/INDEX.md;
              the PUBLISHED results live in papers/scripts/)
papers/              ← THE SCIENTIFIC RECORD
  timing/     Paper 1 manuscript (elsarticle) + figs/ + circulation/review documents
  scripts/    gated result generators: one dir per published result, AUDIT.md + logs
  tables/     generated single-source-of-truth result tables
  figures/    gate-script figure outputs
  memory_*.md + CAMPAIGN_SNAPSHOT_*.md ← claims law, gate history, dataset inventory
  energy_position/  Paper 2 skeleton (companion: energy + position + 4D)
data/                ← per-year data + configs + metadata (ntuples gitignored)
  2023/{raw, reduced, configs, metadata}/   metadata/DATASET_NOTES.md = dataset compendium
docs/                ← APPARATUS_2023.md, METHODS_2023.md, ARCHITECTURE.md, CONVENTIONS.md,
                       STATS_CONVENTIONS.md, design notes, outlines
figures/             ← analysis figure outputs (year-namespaced; narrative/ feeds the paper)
tools/               ← repro.sh (one-command reproduction), run_gate.sh (guarded gate runs),
                       check_figure_manifest.sh (figure drift), pull_argon_data.sh
REPRODUCE.md         ← fresh clone → published numbers, env pinning + expected residual
site/                ← GitHub Pages site (outreach/report; historical pages carry banners)
student/             ← self-contained student lab module
archive/             ← superseded code kept for provenance (do not use)
output/              ← transient per-run scratch (gitignored)
setup.sh             ← source once per shell
```

**Three rules that keep it navigable:**
1. **One include root** — `source setup.sh`; every macro `#include`s headers bare.
2. **One data resolver** — `lib/io/DataPaths.h`: `radReduced(build,E)`, `radRaw(name)`,
   `radConfig(build)`, `radHglg(build)` all resolve under `data/<year>/`.
3. **Dependency direction** — `lib ← reduce ← analyze ← papers/scripts`; never the reverse.

### The library (`lib/`)
| Header | Role |
|--------|------|
| `io/Schema.h` | canonical reduced-ntuple branch layout (`rad` tree, 39 branches) |
| `io/BuildConfig.h` | per-build JSON config + channel map (+ `.hglg` calibration sidecar) |
| `io/RadView.h` | typed reader over a reduced file (`hg_peak`, `timeOf(c,src)`, beam center) |
| `io/DataPaths.h` | **the** path resolver |
| `io/MiniJson.h` | minimal JSON parser for the build configs |
| `waveform/WaveformUtils.h` | DRS4 pulse extraction (pedestal, peak, CFD, charge) |
| `waveform/ChannelConfig.h` | detector geometry + `chanOff()`/`timeOff()` raw indexing |
| `waveform/DRS4Calibration.h` | DRS4 stop-cell recovery + timing correction |
| `physics/SelectionCuts.h` | single source of truth for every cut threshold (with WHY notes) |
| `physics/RadTiming.h` | production primitives: `timingBrightestK()` + robust `tebSigma()` (the retired `timingBestBin()` is kept for the historical record) |
| `viz/PlotUtils.h` | fits, `ScanRunCenters`, plot helpers |
| `viz/RADiCALStyle.h` | unified ROOT plotting style |

---

## Builds & energies

Four builds — LYSO:Ce + W common to all; **the WLS capillary is the variable**:

| Build | WLS capillaries | Light | Energies (GeV) |
|-------|-----------------|-------|----------------|
| `DSB1` | 4 × DSB1 organic (headline; the build+dataset of NIM A 1068 (2024) 169737) | high | 25–150 |
| `LUAG` | 4 × LuAG:Ce ceramic | ~3× less | 50–150 |
| `MIXED` | NE+SW = DSB1, NW+SE = LuAG (same-shower control; valid era runs 2975–3044) | mixed | 50–150 |
| `TENERGY` | 3 × DSB1 timing + 1 full-length energy capillary | high | 50–150 |

Configs: `data/2023/configs/<BUILD>.json` + `.hglg` gain sidecar (`reduce/calibHGLG.C`).
Apparatus: `docs/APPARATUS_2023.md`. Dataset compendium: `data/2023/metadata/DATASET_NOTES.md`.
Methods compendium: `docs/METHODS_2023.md`.

---

## Reducing from raw (Argon HPC)

Reduction (raw waveforms → reduced `rad` ntuple) runs on the Iowa **Argon** cluster.
Follow **`reduce/hpc/REREDUCE.md`** (the current recipe; `reduce/hpc/README.md` describes the
superseded legacy pipeline retained in `archive/hpc_legacy/`):

```bash
# on the Argon login node, from the repo root
qsub reduce/hpc/compile.sh                                  # prebuild reducer .so (once)
RAD_CONFIG=DSB1 RAD_MANIFEST=$PWD/reduce/hpc/manifest_DSB1.csv \
  bash reduce/hpc/submit_reduce.sh                          # reduce[array] -> merge
# then copy home:  data/2023/reduced/<BUILD>/<E>GeV.root
```

Reduced files carry MCP-referenced HG times, `hg_lgcfd` (the srCFD source branch), and TR0
trigger times. Run `reduce/verify.C` after any reduction; `reduce/validateReduce.C` enforces
bit-identity against the legacy reducer.

---

## Data layout & where to get raw

```
data/<year>/raw/        RUN<n>_<E>_GeV.root   (raw 'pulse' tree; subset kept locally)
data/<year>/reduced/    <BUILD>/<E>GeV.root   (analysis 'rad' tree; ~13 GB, gitignored)
data/<year>/configs/    <BUILD>.json + .hglg  (build configs + gain calibration)
data/<year>/metadata/   DATASET_NOTES.md, channel_map.yaml, dataset.yaml, runs.csv, logbook
```

Raw files (~18 GB/build) are gitignored. Sources:

| Location | Path |
|----------|------|
| Iowa Argon cluster | `/Shared/lss_yonel/jwwetzel/RADiCAL_CERN_May2023/rec/rec` |
| CERN lxplus (EOS) | `/eos/experiment/iucmsdata/Data/2023/H2TB-202305-RADiCAL/RADiCAL_2023May/rec` |

---

## Detector & DAQ (background)

Two CAEN **DT5742** digitisers (16+1 channels each) record waveforms via the **DRS4**
switched-capacitor chip: **DRS0 at 5 GS/s** (HG timing channels + MCP reference copies),
**DRS1 at 1 GS/s** (LG amplitude channels + wire chambers). Each DT5742 stores its channels
in two groups of 8 signal + 1 trigger; 1024 samples per channel per event.

The module (14×14 mm², W/LYSO shashlik) is read out by 8 capillary ends (4 corners ×
UP/DW): **HG** copies (→ timing, clip ~820 mV) and **LG** copies (→ amplitude, linear),
plus the MCP-PMT timing reference (split into both DRS0 groups) and the TR0 trigger
scintillator. Channel map: `lib/waveform/ChannelConfig.h`, the per-build
`data/2023/configs/<BUILD>.json`, and `docs/APPARATUS_2023.md`.

### Raw data format (reference)

Per event, all channels are stored end-to-end in two flat arrays in the `pulse` tree:
`amplitude` (waveforms) and `timevalue` (time axes):

```cpp
// drs ∈ {0,1}, group ∈ {0,1}, ch ∈ 0..8
channel_index = (1024*9*2)*drs + (1024*9)*group + 1024*ch
time_index    = (1024*2)*drs   + 1024*group
```

---

## Papers

- **Paper 1 — Timing** (`papers/timing/`): at internal coauthor circulation. σ_t(E) for four
  WLS builds, srCFD saturation recovery, same-shower MIXED control, systematics budget,
  depth-asymmetry drift. Index: `papers/README.md`.
- **Paper 2 — Energy + Position (+4D)** (`papers/energy_position/`): companion, in preparation.

Workspace-release status (license, citation, data DOI, tags):
`WORKSPACE_REORG_PLAN_2026-06-10.md` §5. Full panel audit: `WORKSPACE_AUDIT_2026-06-10.md`.

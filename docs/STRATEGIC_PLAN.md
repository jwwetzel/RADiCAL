# RADiCAL Analysis — Strategic Plan
> *Positioning RADiCAL as a standard candle for precision calorimetric timing.*

---

## Executive Context

We have confirmed ~37–56 ps timing resolution across 25–150 GeV with three complementary
methods. The path to understanding where the remaining resolution goes — and to demonstrating
this detector as a reproducible platform — requires answering five questions simultaneously:

1. **What are we actually selecting?** The PbGlass scatter shows structure that our current
   cuts may not fully exploit.
2. **Which levers matter most?** CFD fraction, energy binning, walk correction, fiducial
   tightening, shower containment — we need these compared on the same axis.
3. **Can we find it faster?** The per-energy folder structure makes trend-spotting hard.
   Seeing all six energies side-by-side immediately reveals whether a problem is a single
   run artifact or a systematic.
4. **How do we hand this off?** A future collaborator must be able to plug in a new crystal
   geometry and rerun everything without reading source code.
5. **How do we present it?** A static web report replaces the PDF hunt and lets anyone
   navigate straight to the plot they need.

---

## Phased Roadmap

```
Phase 1 — Reorganize & Reveal       (immediate, ~1-2 days coding)
  ├── 1A  Cross-energy collated plots
  ├── 1B  PbGlass dual-band investigation
  └── 1C  Static web report (HTML + JS, no server needed)

Phase 2 — Deepen the Timing Study   (next, ~3-5 days)
  ├── 2A  Containment cut scan (does tighter PbGlass cut improve σ_t?)
  ├── 2B  Position-resolved timing (σ_t vs x,y — lateral shower position effect)
  ├── 2C  MCP jitter decomposition (subtract MCP self-resolution from σ_t)
  └── 2D  Channel combination optimization (which subset of 8 gives best σ_t?)

Phase 3 — Platform Architecture     (medium-term, ~1 week)
  ├── 3A  Detector config file (JSON) — decouple geometry from code
  ├── 3B  Run config file (JSON) — decouple run list from code
  ├── 3C  N-channel generalization — remove hardcoded 8-channel arrays
  └── 3D  runAll.sh → unified pipeline script with config injection

Phase 4 — Standard Candle Deliverable  (publication-ready)
  ├── 4A  Final combined plots for paper figures
  ├── 4B  Systematic uncertainty table
  └── 4C  Reproducibility package (one-command rerun from raw data)
```

---

## Phase 1A — Cross-Energy Collated Plots

### Problem
Currently each energy lives in its own folder. To see whether σ_t at 100 GeV is anomalous
you must open six PDFs and flip between them. Trends are invisible.

### Solution: `Analysis/compareEnergies.C`

A new macro that reads the per-energy ntuples and summary ROOT files and produces one PDF
per *topic*, with all six energies on each page (2×3 or 3×2 grid). It writes to
`Analysis/Output/Summary/cross_energy_*.pdf`.

### Plot Groups and Rationale

Each group answers a specific physics question about timing resolution.

**Group 1 — Beam & Selection Purity**
*Question: are we selecting the same population at each energy?*
```
cross_energy_beam_quality.pdf
  Page 1: Hit maps (2×3 grid, all energies) — is the beam centred consistently?
  Page 2: MCP1 amplitude distributions (all energies overlaid) — saturation trend?
  Page 3: sum_lg energy distributions (all energies overlaid) — linearity visible?
  Page 4: Fiducial efficiency vs radius (all energies) — does 3 mm always work?
```

**Group 2 — Shower Containment**
*Question: what fraction of showers are truly contained, and does it change with energy?*
```
cross_energy_containment.pdf
  Page 1: sum_pb/sum_lg ratio all energies overlaid — does hadronic fraction grow?
  Page 2: PbGlass scatter (6-panel) — scaling of the two-population structure
  Page 3: Containment fraction vs r_beam (6-panel) — edge effects vs energy?
  Page 4: Contained fraction vs energy (single plot) — headline containment number
```

**Group 3 — Channel Performance**
*Question: are all 8 capillaries contributing equally, and does that change with energy?*
```
cross_energy_channels.pdf
  Page 1: σ_t per channel at each energy (8 lines on one plot) — which dominates?
  Page 2: HG amplitude per channel (6-panel, each panel = one energy, 8 channels)
  Page 3: Walk correction slope k_i per channel vs energy — is correction stable?
  Page 4: HG/LG correlation r^2 per channel vs energy — proxy for shower quality
```

**Group 4 — Timing Methods Comparison**
*Question: which method wins at each energy, and by how much?*
```
cross_energy_timing_methods.pdf
  Page 1: σ_t vs energy, one line per method (A, B, C + paper reference)
  Page 2: CFD fraction comparison (10/20/30/50%) vs energy — is 20% optimal?
  Page 3: Method A vs C improvement factor vs energy — is walk correction worth it?
  Page 4: σ_t vs 1/√E with stochastic fits (the headline physics plot)
```

**Group 5 — Optimization Scan Summary**
*Question: which cuts are we most sensitive to?*
```
cross_energy_cut_scans.pdf
  Page 1: σ_t vs MCP amplitude cut at all energies (multi-line)
  Page 2: σ_t vs fiducial radius at all energies
  Page 3: σ_t vs HG threshold at all energies
  Page 4: σ_t vs containment cut (kPb_maxRatio) at all energies [NEW — Phase 2]
```

### Design Rules for `compareEnergies.C`
- Pure reader: reads only existing ROOT files, never re-processes raw ntuples
- One function per plot group; main() calls them in sequence
- Plot style via the existing `StylePad()` from PlotUtils.h — zero new style code
- All energy-loop logic uses the same `kRuns[]` array from RunConfig.h

---

## Phase 1B — PbGlass Dual-Band Investigation

### What We're Seeing

The 25 GeV scatter plot (and likely all energies) shows two distinct populations
separated by the containment cut line. Looking at the image:

**Band 1 (horizontal, near sum_pb ≈ 0):** Events that deposited energy in RADiCAL and
very little in PbGlass. These are our signal — fully contained EM showers.

**Band 2 (near-vertical, sum_lg ≈ 0–300 mV):** Events that deposited little or nothing
in RADiCAL. These are beam halo particles that missed the crystal bundle entirely but
still created a WC track. They may or may not hit PbGlass.

**Above the cut line (scattered):** Events where sum_pb > 30% × sum_lg. Likely pion
contamination or edge showers. At 25 GeV the SPS H2/H4 beam has higher hadronic
contamination than at 150 GeV.

### Why It Matters for Timing

If hadronic events or halo events pass our fiducial cut (they can, if the pion travels
through the fiducial cylinder before showering), they will dilute the timing sample.
Hadrons have longer shower-development time and non-uniform longitudinal profiles —
both of which degrade CFD precision. Our 30% containment cut addresses this, but we
have never measured HOW MUCH it helps.

### Solution: `Analysis/investigatePbGlass.C`

This macro classifies every event into four populations and studies each separately:

```
Population A — "Good EM shower"
  sum_lg > kSumLG_centroid  AND  sum_pb < kPb_maxRatio × sum_lg
  → our signal; should have the best timing

Population B — "Punch-through / hadronic"
  sum_lg > kSumLG_centroid  AND  sum_pb >= kPb_maxRatio × sum_lg
  → rejected by current cut; quantify timing degradation

Population C — "Beam halo (missed RADiCAL)"
  sum_lg < kSumLG_centroid  AND  r_beam < kFiducial_r_timing
  → within timing fiducial but no RADiCAL signal; these must be cut
  NOTE: does the current kSumLG_centroid cut actually catch all of these?

Population D — "Off-axis / noise"
  r_beam >= kFiducial_r_timing
  → already excluded by fiducial; sanity check
```

**Output: `Analysis/Output/Summary/pbglass_investigation.pdf`**

```
Page 1: 2D scatter (all WC-ok) with populations colour-coded by label A/B/C/D
Page 2: Hit maps for A, B, C side-by-side — do hadronic events come from the edge?
Page 3: Timing distributions for A vs B at each energy — how much does B hurt?
Page 4: σ_t for population A vs B vs (A+B) vs energy — quantify the containment benefit
Page 5: sum_pb/sum_lg ratio for fiducial events vs energy — is pion fraction growing?
Page 6: sum_lg spectra for A and B overlaid — do they overlap?
```

**Key questions this answers:**
- Is our 30% containment cut in the right place?
- At what energy does hadronic contamination become negligible?
- Should we tighten kPb_maxRatio to 20%? What does that cost in efficiency?
- Does population C actually pass the current fiducial cut? (Could be a residual cut gap)

---

## Phase 1C — Static Web Report

### Why a Web Report, Not More PDFs

PDFs require sequential navigation. The web report lets you jump directly from
"I want to see the 100 GeV timing" to that exact plot in one click. It also allows
side-by-side comparison of plots that are physically in different PDFs.

### Implementation: `Analysis/makeReport.py`

A Python script (stdlib only — no external dependencies) that:
1. Converts all PDF pages to PNG using `ghostscript` (already installed with ROOT environments)
2. Generates a single static `Analysis/Output/report.html` with all plots embedded
3. Uses vanilla JavaScript for the sticky sidebar TOC and smooth scroll
4. No server, no npm, no build step — `open Analysis/Output/report.html` is all it takes

### Report Structure and Navigation

```
RADiCAL Analysis Report                        [sidebar TOC — sticky, scrollable]
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━            ┌─────────────────────────────────────┐
                                          │ 0  Executive Summary                │
  0  EXECUTIVE SUMMARY                   │ 1  Beam & Data Quality              │
     Best timing: 37 ps @ 150 GeV        │    1.1  Hit maps                    │
     Stochastic term: a = 181 ps√GeV     │    1.2  MCP amplitude               │
     Constant term:   b = 34.9 ps        │    1.3  Energy distributions        │
     Status vs paper: [table]            │    1.4  Channel integrity            │
                                         │ 2  Shower Containment               │
  1  BEAM & DATA QUALITY                 │    2.1  PbGlass scatter              │
     [cross-energy hit maps]             │    2.2  Dual-band analysis           │
     [MCP amplitude all energies]        │    2.3  Containment efficiency       │
     [energy distributions]              │ 3  Timing Methods                   │
                                         │    3.1  Method A/B/C vs energy      │
  2  SHOWER CONTAINMENT                  │    3.2  CFD fraction comparison     │
     [PbGlass scatter 6-panel]           │    3.3  Walk correction benefit     │
     [dual-band classification]          │    3.4  Channel-by-channel σ_t      │
     [containment vs radius]             │ 4  Best Timing Resolution           │
                                         │    4.1  Energy-binned estimator     │
  3  TIMING METHODS                      │    4.2  σ_t vs 1/√E (stochastic)   │
     [method comparison table]           │    4.3  vs published result         │
     [σ_t vs energy all methods]         │ 5  Optimization Studies             │
     [CFD fraction comparison]           │    5.1  MCP cut scan                │
                                         │    5.2  Fiducial scan               │
  4  BEST TIMING RESOLUTION              │    5.3  HG amplitude scan           │
     [headline plot: σ_t vs 1/√E]        │    5.4  Containment cut scan       │
     [comparison to arXiv:2401.01747]    │ 6  Per-Energy Detail               │
                                         │    6.1  25 GeV                      │
  5  OPTIMIZATION STUDIES                │    6.2  50 GeV                      │
     [cut scan overlays]                 │    ...                              │
                                         └─────────────────────────────────────┘
  6  PER-ENERGY DETAIL
     [click-through to individual energy plots]
```

### Technical Spec for `makeReport.py`

```python
# Guiding principles:
#   - One function per section; main() is a flat list of section() calls
#   - PNG conversion via subprocess call to gs (ghostscript)
#   - HTML template as a raw string with {placeholder} substitution
#   - TOC auto-generated from section list — no manual HTML editing
#   - Each plot has a required caption string (forces documentation)

def convert_pdf_pages_to_png(pdf_path, output_dir, dpi=150) -> list[Path]
def make_section(title, anchor, plots: list[PlotEntry]) -> str
def make_toc(sections: list[Section]) -> str
def make_report(sections, output_path)

class PlotEntry:
    png_path: Path
    caption: str          # required — if empty, script fails with clear error
    width_pct: int        # 33, 50, or 100 — controls layout
```

**Invocation:**
```bash
python3 Analysis/makeReport.py
# → Analysis/Output/report.html
open Analysis/Output/report.html
```

**runAll.sh integration:** `makeReport.py` is the last step in `runAll.sh`, so a full
rerun always regenerates the report.

---

## Phase 2 — Deepen the Timing Study

### 2A — Containment Cut Scan (`Analysis/timingContainmentScan.C`)

The most important unanswered question: **does tighter shower containment improve σ_t?**

We apply the PbGlass cut for purity but we've never measured its timing benefit.
This macro scans `kPb_maxRatio` from 0.10 to 0.50 in steps of 0.05, computes σ_t for
Method A at each cut value, and overlays curves for all six energies. Output:
`Analysis/Output/Summary/timing_containment_scan.pdf`

Expected result: σ_t should improve as the cut tightens (purer EM sample), then plateau
when hadronic contamination is negligible (probably at ~0.20 for 50+ GeV).

### 2B — Position-Resolved Timing (`Analysis/timingVsPosition.C`)

Does σ_t depend on where in the crystal face the shower starts?

A 2D profile histogram `TProfile2D` of σ_t vs (x_trk, y_trk) per energy. The key
physics question: is there a position-dependent timing bias from the fibre geometry?
If crystal boundaries (edges between capillaries) degrade CFD performance, that will
appear as high-σ bands at the crystal boundaries. This directly motivates tighter
fiducial cuts or per-position calibration.

### 2C — MCP Jitter Decomposition

Our measured σ_t includes the MCP timing reference jitter. The MCP self-resolution
can be estimated from the MCP1–MCP2 time difference distribution (which we already
store as `mcp2_time − mcp_time`). If we assume σ_MCP1 ≈ σ_MCP2:

  σ_MCP_single = σ(MCP1−MCP2) / √2

Then the true crystal timing is: σ_crystal = √(σ_measured² − σ_MCP²)

This should already be possible from existing ntuples. Add a page to `analyzeResolution.C`
or create `Analysis/mcpJitter.C`. This directly tells us how much headroom we gain from
a better time reference.

### 2D — Channel Combination Optimization

We always use all 8 channels. But channel 7 (SW-U) uses MCP2 as reference and has
historically had different performance. Questions:
- What is σ_t if we drop ch7 entirely?
- What is σ_t for (DW − UP) where DW = {ch0,ch1,ch2,ch3} and UP = {ch4,ch5,ch6,ch7}?
- What is σ_t for each individual capillary as a function of energy?
- Is there an optimal subset? (brute-force 255 subsets is feasible; add to timingMethods.C)

---

## Phase 3 — Platform Architecture

### The Standard Candle Vision

Future users will plug in a different crystal (LYSO instead of scintillating fibre,
different geometry, 4 or 16 channels instead of 8). The analysis chain must accept a
config file and produce the same report without touching a line of C++.

### 3A — Detector Config (`Analysis/config/detector.json`)

```json
{
  "name":        "RADiCAL-v1",
  "description": "8-channel scintillating-fibre calorimeter",
  "n_channels":  8,
  "channels": [
    { "id": 0, "name": "NE-D",  "role": "downstream", "use_mcp2": false },
    { "id": 1, "name": "NW-D",  "role": "downstream", "use_mcp2": false },
    { "id": 2, "name": "SE-D",  "role": "downstream", "use_mcp2": false },
    { "id": 3, "name": "SW-D",  "role": "downstream", "use_mcp2": false },
    { "id": 4, "name": "NE-U",  "role": "upstream",   "use_mcp2": false },
    { "id": 5, "name": "NW-U",  "role": "upstream",   "use_mcp2": false },
    { "id": 6, "name": "SE-U",  "role": "upstream",   "use_mcp2": false },
    { "id": 7, "name": "SW-U",  "role": "upstream",   "use_mcp2": true  }
  ],
  "reference_detectors": ["MCP1", "MCP2"],
  "geometry_mm": { "crystal_pitch": 2.5, "n_rows": 2, "n_cols": 4 }
}
```

### 3B — Run Config (`Analysis/config/runs.json`)

```json
{
  "dataset": "CERN_May_2023",
  "runs": [
    { "label": "25GeV",  "energy_GeV": 25,  "file": "run_25GeV.root" },
    { "label": "50GeV",  "energy_GeV": 50,  "file": "run_50GeV.root" },
    ...
  ]
}
```

This replaces `RunConfig.h`. The C++ reads it via a small `ConfigLoader.h` (JSON
parsed with ROOT's built-in `TFile`/TObject or a tiny nlohmann/json include).

### 3C — N-Channel Generalization

Every hardcoded `[8]` array becomes `[kNChannels]` where `kNChannels` is loaded from
`detector.json` at startup. The "downstream" and "upstream" roles replace the implicit
even/odd indexing. This is a mechanical but important refactor.

### 3D — Pipeline Script (`runAll.sh` → `run.sh`)

```bash
./run.sh --config Analysis/config/runs.json --detector Analysis/config/detector.json
```

All steps (processRun → qualityPlots → analyzeResolution → timingMethods →
timingResolution → timingEnergyBins → compareEnergies → investigatePbGlass →
makeReport) run in sequence. Errors abort with a clear message. The report appears
at `Analysis/Output/report.html`.

---

## Code Quality Standards

All code in this project — new and edited — must satisfy these rules.
They are listed in order of importance.

### 1. Single Responsibility
Every function does exactly one thing. A function that fills a histogram AND draws it
AND fits it is three functions. `FitGaussCore`, `DrawFitOverlay`, `ScanRunCenters` are
the existing model — keep to that pattern.

### 2. Named Constants, Never Magic Numbers
Every threshold lives in `SelectionCuts.h` or the JSON config. No `> 50.f` buried in
a loop. If it appears twice it should be a constant; if it appears once it should still
be a named constant with a comment explaining the physics.

### 3. Explicit Over Clever
Variable names describe physics, not implementation: `sigma_timing_ps` not `s`,
`n_fiducial_events` not `n2`. Array indices use named enums where possible.
A physicist who has never seen the codebase should be able to read any function in
under two minutes.

### 4. Comments Explain Why, Not What
```cpp
// BAD: fill the histogram
hSumLG->Fill(sum_lg);

// GOOD: use_mcp2 channels carry MCP2-referenced times; without this check,
// a saturated MCP2 would bias the (DW+UP)/2 method toward early crossing times.
if (kCap[i].use_mcp2 && mcp2_peak > kMCP2_maxPeak) continue;
```

### 5. Self-Contained Functions
Functions take explicit parameters; no reading global state. Every analysis macro
can be called standalone with a path argument. No macro assumes another has already
run (check for the ntuple; print a clear message if absent).

### 6. Fail Loudly
If a file is missing, a fit diverges, or a branch is absent, print a descriptive
message and return — never silently produce an empty or wrong plot.

---

## Timing Resolution Lever Summary

For reference when prioritizing investigations, here is the current understanding of
what limits σ_t and what we can pull:

| Lever | Current Status | Expected Gain | Phase |
|-------|---------------|---------------|-------|
| MCP saturation cut | Done (kMCP1_maxPeak=750mV) | ~1 ps | 1 |
| Walk correction (M7) | Implemented | ~5–10 ps at 25 GeV | done |
| Energy binning | Implemented | ~15 ps at 150 GeV | done |
| Shower containment cut | Implemented but not scanned | TBD | 2A |
| Fiducial tightening | Scanned per energy | ~2–5 ps | done |
| CFD fraction | Scanned (10/20/30/50%) | ~2–3 ps | done |
| MCP jitter subtraction | NOT YET DONE | ~5–10 ps | 2C |
| Channel combination | NOT YET DONE | ~3–8 ps | 2D |
| Position-resolved correction | NOT YET DONE | ~2–5 ps | 2B |
| Better time reference (SiPM MCP) | Hardware change | ~10–20 ps | future |

The two highest-value investigations that require zero hardware changes are:
**MCP jitter subtraction (2C)** and **channel combination optimization (2D)**.
Both should move to Phase 1 if time permits.

---

## Immediate Next Steps (ordered)

1. `Analysis/compareEnergies.C` — Groups 1–4 above (1 day)
2. `Analysis/makeReport.py` — PNG conversion + HTML assembly (1 day)
3. `Analysis/investigatePbGlass.C` — four-population classification (1 day)
4. `Analysis/timingContainmentScan.C` — containment vs timing (half day)
5. Add PNG export (`SaveAs`) to all existing macros so the report stays current
6. Wire everything into `runAll.sh`
7. JSON config files (detector + runs) — removes last hardcoded detector specifics

---

*Last updated: May 2026*
*Authors: RADiCAL analysis group*

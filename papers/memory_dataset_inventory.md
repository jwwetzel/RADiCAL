# Dataset Inventory — what exists, where, and what every branch means

Updated: 2026-06-09. Branch list verified by direct inspection of all four builds' 150GeV files
(identical 39-branch schema). Authoritative channel map: `data/2023/configs/<BUILD>.json` +
memory `dataset_2023.md` (outside repo, ~/.claude/.../memory/).

## Reduced data (analysis-ready, 'rad' tree)
| build | energies (GeV) | entries @150 | WLS configuration |
|---|---|---|---|
| DSB1 | 25,50,75,100,125,150 | 1.83M | 4× DSB1 organic T-type (the published build) |
| LUAG | 50,75,100,125,150 | 1.67M | 4× LuAG:Ce ceramic T-type |
| MIXED | 50,75,100,125,150 | 264k | NE,SW = DSB1; NW,SE = LuAG (map: GATE 1 pending logbook) |
| TENERGY | 50,75,100,125,150 | 807k | 3× DSB1 T-type + 1× E-type (energy capillary) |

Paths: `data/2023/reduced/<BUILD>/<E>GeV.root` · configs: `data/2023/configs/<BUILD>.json` (+ `.hglg`
calibration files) · manifests: `data/2023/MANIFEST.csv`, `reduce/hpc/manifest_*.csv`.

## Verified branch schema (39 branches, all builds)
- Event/beam: `run/I`, `event/I`, `beam_energy/F`
- Tracking: `wc_ok/O`, `x_trk/F`, `y_trk/F`, `in_fiducial/O`, `wc_peak[4]/F`
- References: `mcp1_peak,mcp1_time`, `mcp2_peak,mcp2_time`, `tr0a_*`, `tr0b_*` (all /F)
- DRS metadata: `stopcell[4]/I` (per-chip stop cell — enables stop-cell timing correction, GATE 4)
- HG (timing) per capillary end [8]: `hg_peak`, `hg_ped_rms`, `hg_saturated/O`, `hg_spike/O`,
  crossings `hg_cfd03/05` (`hg_cfd`=cfd20) `/10/30/50`, `hg_led`, `hg_lgcfd`, `hg_tot`, `hg_charge`
- LG (energy) [8]: `lg_peak`, `lg_charge`; sums: `sum_lg`, `sum_pb`; `pb_peak[4]`
- Legacy/raw-slot: `s_peak[36]`, `s_cfd05[36]`, `s_charge[36]`

Channel-index convention (kCap, universal across builds):
0=NW-D 1=NE-D 2=SE-D 3=SW-D (downstream) · 4=NW-U 5=NE-U 6=SE-U 7=SW-U (upstream).
Timing sources enum (RadView): kCFD03,05,10,20,30,50, kLED (20 mV fixed), kLGCFD (srCFD: CFD at 15% of
LG-predicted true peak — the clipped-pulse recovery). Times in ns; HG clips ~820 mV (DT5742 channel clip).

## Raw data (local; full set on Argon — ssh -p 40, setup_root)
Local `data/2023/raw/`: DSB1 all 6 energies (RUN1034/1075/1112/1148/1211/1258-61); MIXED 5 energies
(RUN2941/2913/2995/2985/2787/2975); LUAG 50,150 only (RUN2389/2186); TENERGY 50,150 only (RUN2722/2530).
Raw 'pulse' tree: `timevalue[]`,`amplitude[]` 1024-sample DRS4 waveforms.

## DAQ (from dataset_2023 memory; verify against configs before printing)
2× CAEN DT5742. DRS0 @5 GS/s: 8 HG timing + both MCP copies + 4 Pb-glass + trigger scints
(same-group wiring NW/NE/SE + MCP1 in G0 is DELIBERATE — differential cancellation).
DRS1 @1 GS/s: 8 LG energy + wire chambers. SiPM bias 42.25 V everywhere.

## Selection working points (lib/physics/SelectionCuts.h, post-fix)
MCP1 window 200–750 mV · kHG_minPeak 20 mV · TimingFiducialR: 2.5 mm (≤100 GeV) → 3.0 mm (≥125),
linear ramp between (continuous since 2026-06 fix) · in-event channel-consistency veto
kTimingChanConsistency_ns = 2.0 · containment: sum_pb < 0.30·sum_lg (report-era; verify).

## Analysis products supporting current claims
- Chain of evidence: `chain_of_evidence.html` (waveforms→reduction→fits→timing, all builds)
- Production estimator (post-fix): `lib/physics/RadTiming.h` (tebSigma + eventDWUP veto)
- Survey/compare/fit macros: `analyze/studies/{timingAllMethods,methodCompare,methodDist}.C`
- MIXED: `analyze/studies/mixedSeparate.C`, kill-shot `figures/2023/narrative/mixed_h2h.png`
- Fix evidence: `figures/2023/narrative/monotonicity_{fix,evidence}.png`, `analyze/studies/{sigmaProbe,outlierPeek,timingRegression}.C`
- Reduction QA: `figures/2023/narrative/reduction_<BUILD>.png`
- Depth dial (GATE 3): `papers/scripts/depth_dial/`, `papers/figures/depth_dial/`

## Known data caveats
- MIXED per-energy files merge 8–34 runs with drifting beam centroid; per-energy N is 16k–53k fiducial.
- 25 GeV exists for DSB1 ONLY — every cross-build figure states "50–150 GeV" for non-DSB1.
- UPSE SiPM had lower gain in the parent's data (parent Fig 14) — inter-calibration absorbed it; SE-U/SW-U
  run 10–15% low in-fiducial (WLS path length) per the report-era non-uniformity scan.
- `hg_saturated`/`hg_spike` flags exist per end — use them rather than re-deriving clip status ad hoc.

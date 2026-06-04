# RADiCAL follow-on paper — inventory, gap analysis & plan

Companion to the published paper **arXiv:2401.01747** (*"Study of time and energy
resolution of an ultra-compact sampling calorimeter (RADiCAL) module at EM shower
maximum, 25–150 GeV"*, RADiCAL Collaboration; Perez-Lara, Wetzel, Akgun et al.).
This file plans the **next** paper from the same data.

## 1. What the published paper established (don't re-derive)

- **Module:** 29 LYSO:Ce (1.5 mm) + 28 W (2.5 mm), 14 × 14 × 135 mm³ = **25 X₀**
  (X₀ = 5.4 mm, R_M = 13.7 mm). Contains longitudinally; ~1 R_M wide → leaks
  transversely (a 3×3 array is needed for full containment).
- **Readout:** **4 T-type capillaries** (DSB1 WLS filament at shower max), each read
  by a SiPM at **both ends** → 8 channels (DW/UP × NW,NE,SE,SW), each with low-gain
  (energy) + high-gain (timing) amplification. Hamamatsu HDR2 SiPMs. CAEN DT5742.
- **Reference:** Hamamatsu R3809U-50 MCP, 10–20 ps.
- **Timing method:** fixed threshold on the HG leading edge vs the MCP. Best
  estimator = **BestMinus = (Δt_DW − Δt_UP)/2**, evaluated in **measured-energy bins
  6–8**. Six methods compared (Fig 27). Documented a **15 %, ~0.2 ns satellite peak**
  from fixed-threshold + 1-DRS-sample quantization (cancels in the DW−UP difference).
- **Timing result:** σ_t = 256/√E ⊕ 17.5 ps; **27 ps @ 150 GeV.**
- **Energy:** sum of 8 LG amplitudes = **shower-max-localized** energy (≈3 LYSO
  tiles), explicitly NOT full-module. σ_E/E = 52.04 %/√E ⊕ 31.62 %/E ⊕ 9.31 %.
- **Section 7 roadmap (their stated next steps):**
  1. *Data analysis — EM shower localization* (capillary energy × wire-chamber position).
  2. *Data analysis — comparison of DSB1 and LuAG:Ce WLS* (T-type LuAG:Ce vs the DSB1
     data presented, 25 GeV steps, 25–150 GeV). **← our follow-on.**
  3. *Instrumentation — enhanced modular design* (new hardware; future).

## 2. Corrections made to our own internal docs (done)

Reading the PDF corrected two errors we had been propagating:
- **Depth:** we called this a "~4 X₀ / 14 mm short, leakage-dominated prototype." It
  is **25 X₀ deep**; the "14 mm" is transverse. The modest σ_E is a **localized
  shower-max measurement**, not longitudinal leakage. Fixed in `paper.html`,
  `report/index.html`, `geant4_lab.html`.
- **Attribution:** the (DW−UP)/2 corner estimator (= published "BestMinus") and the
  energy-binned best-bin are the **published method**, not ours. Our contribution is
  the rigour on top (CFD-5%, OOS, timebase study). Fixed the framing in `paper.html`.
- Also corrected the "eight capillaries" overcount → **four capillaries read at both
  ends.**

## 3. Inventory — our findings vs roadmap item 2

| Needed for the LuAG:Ce-vs-DSB1 paper | Status | Source |
|---|---|---|
| LuAG / MIXED / TENERGY + DSB1 reduced data | ✅ have | `reduced/`, local |
| Per-config σ_t & σ_E/E, all-fiducial, all E | ✅ have | `configResolution.C` |
| In-event confound-free DSB1 vs LuAG per-capillary | ✅ have (cleanest result) | `mixedHeadToHead.C` |
| OOS best-bin σ_t **@150 only** | ⚠️ partial | `compareConfigsPlot.C` |
| CFD-5% (removes their 0.2 ns satellite) | ✅ have (improvement over fixed threshold) | `WaveformUtils.h` |
| OOS cross-validation + DRS4 timebase study | ✅ have (rigour they lacked) | `timingEnergyBins.C`, `drs4TimeBase.C` |
| **σ_t vs light yield** (DSB1+LuAG on one curve) | ❌ **missing** — the key physics plot | — |
| Per-config σ_t(E) OOS at **all** energies + fits | ❌ missing | extend `configResolution`/OOS |
| Per-config σ_E(E) with a/√E⊕b/E⊕c fits | ❌ missing | extend |
| Per-config systematics | ❌ to run | `systematicUncertainties.C` |

Key cross-config numbers in hand (all-fiducial / OOS-best-bin @150 ps): DSB1
58.4 / 29.4, LuAG 49.8 / 36.6, MIXED 54.5 / 36.7, TENERGY 57.9 / 35.7.

## 4. Gap analysis — are we complete? **No (~70 %).**

Missing before drafting:
1. **σ_t(E) OOS best-bin at all six energies, per config**, + a/√E⊕b fits → the
   headline DSB1-vs-LuAG curve (today we only have @150).
2. **σ_E/E(E) per config** with a/√E⊕b/E⊕c fits.
3. **The σ_t-vs-light-yield plot** — DSB1 and LuAG capillaries on ONE curve
   (light yield from the LG amplitude as a relative proxy; no absolute p.e.
   calibration needed for the single-curve argument). This is the paper's thesis
   and the measured analog of the published GEANT4 Fig 8. **Highest priority.**
4. **Per-config systematics** (cut-variation, already scripted).
5. **Statistics:** DSB1/LuAG/TENERGY ample; **MIXED thin at high E** (~18 k fiducial
   @150 → best-bin may miss the ≥500/bin OOS cut). Either reduce more MIXED runs
   from Argon (analyzed 190/252) or present MIXED all-fiducial.

**Data sufficiency:** the data exist (DSB1+LuAG+MIXED+TENERGY reduced ntuples; raw on
Argon). No new beam time needed. Only MIXED may want a few more reduced runs.

## 5. Paper scoping decision

**Recommendation: ONE focused paper now — on TIMING — with energy as support;
keep EM-shower-localization as a separate later paper.**

Reasoning:
- The **strong, novel result is timing**: per capillary DSB1 ≈ LuAG:Ce, light yield
  (not species) drives timing, so a radiation-hard crystal that is "bright enough"
  times just as well. This is roadmap item 2 and addresses the program's
  radiation-hardness motivation. It carries a paper on its own.
- The **energy resolution here is a supporting null** — shower-max-localized
  (~14–17 %), build-independent. It is **not strong enough to anchor its own paper**
  (a full-module energy resolution is impossible from a single ~1 R_M-wide module —
  that needs the future 3×3 array). So **do not split timing from energy**; fold the
  per-material energy characterization in as a supporting section. Splitting it off
  would create a weak second paper and dilute both.
- **EM shower localization (roadmap item 1)** is a genuinely different observable
  (position, not time/energy). Bundling it would dilute the timing story. **Keep it
  as a separate, later paper.**
- Method refinements (CFD-5%, OOS, timebase, de-saturation): CFD-5% goes in the
  methods of the timing paper (it removes the published 0.2 ns satellite — a concrete
  improvement); OOS + timebase + de-saturation as a short methods/appendix.

So: **Paper A (now)** — *Comparison of DSB1 and LuAG:Ce wavelength shifters for
precision timing in RADiCAL modules* (timing headline + supporting per-material
energy + improved method). **Paper B (later)** — *EM shower localization at shower
max*. Not two papers split by timing-vs-energy; split by **story/observable**.

## 6. Ordered task plan → draft

1. ✅ Fix the depth & attribution errors (done).
2. Extend `configResolution`/OOS to produce **σ_t(E) OOS best-bin + σ_E(E), per
   config, all energies, with fits.**
3. Build the **σ_t-vs-light-yield** plot (DSB1 + LuAG on one curve). *Highest value.*
4. Run **per-config systematics**; check MIXED stats, reduce more from Argon if needed.
5. Draft **Paper A** outline (NIM A; same collaboration; cite 2401.01747 + LuAG:Ce
   Hu et al. NIM A 954 (2020) 161723).

## 7. Target

NIM A, RADiCAL Collaboration. Frames as the realization of arXiv:2401.01747
Section 7 item 2. Radiation hardness of LuAG:Ce is the motivation (ref Hu et al.).

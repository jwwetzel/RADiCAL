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

**FINAL STRUCTURE — two papers, split by observable/story (full outlines in
`PAPER_OUTLINES.md`):**

**Paper 1 — TIMING.** *DSB1 vs LuAG:Ce wavelength shifters for precision timing.*
Configs DSB1 / LuAG / MIXED (all T-type). Light yield, not crystal species, drives
timing → a radiation-hard crystal that is bright enough times just as well (roadmap
item 2; radiation-hardness motivation). Strong, nearly ready.

**Paper 2 — ENERGY + POSITION (+ the 4D combination).** *Shower-max energy and
position measurement with WLS capillaries.* Energy (shower-max-localized σ_E/E per
build) **and** position (the 4-corner light center-of-gravity, ~1.5 mm vs the wire
chamber) come from the *same* shower-max sampling — the published paper says so
explicitly. Anchored by the **E-type (full-length WLS) energy capillary** and the
**TENERGY module** as the simultaneous time+energy+position ("4D") demonstration.
Covers roadmap item 1 (localization) + the energy/E-type direction.

Why this split:
- Two thin halves (energy alone = published null; localization alone = a position
  measurement) combine into one solid paper because they are the **same measurement**
  (shower-max light → energy *and* position).
- **TENERGY is the combination's home (Paper 2)**, not the timing paper: its three
  timing caps are just DSB1 (already in Paper 1); its novel element is the E-type
  *energy* cap, and it demonstrates time+energy+position in one module.
- Method refinements (CFD-5% removing the published 0.2 ns satellite, OOS, timebase,
  de-saturation) → methods/appendix of Paper 1.

## 6. Ordered task plan → draft

1. ✅ Fix the depth & attribution errors (done).
2. ✅ **σ_t(E) OOS best-bin + σ_E(E), per config, all energies, with fits**
   (`configResolutionFull.C` → `capillary_figs/config_sigmat_vs_E.png`,
   `config_sigmaE_vs_E.png`).
3. ✅ **σ_t-vs-light-yield** plot, DSB1 vs LuAG (`lightYieldTiming.C` →
   `capillary_figs/sigmat_vs_lightyield.png`) — first pass; intrinsic
   (reference-subtracted) version is the paper-grade refinement still TODO.
4. ✅ **Per-config systematics** (`configSystematics.C`, on the stable all-fiducial
   σ_t; the jumpy best-bin systematic is selection-jitter-dominated — use OOS).
5. ⏳ Draft **Paper A** outline (NIM A; cite 2401.01747 + Hu et al. NIM A 954 (2020) 161723).

## 7. Results — first pass (June 2026)

**σ_t (OOS best-bin, (DW−UP)/2, CFD-5%, 5-fold run-folded, FitGaussCore) [ps]:**

| E (GeV) | DSB1 | LuAG | MIXED | TENERGY |
|---|---|---|---|---|
| 25  | 50.4 | — | — | — |
| 50  | 42.9 | 64.7 | 46.0 | 51.4 |
| 75  | 35.3 | 51.9 | 44.2 | 44.8 |
| 100 | 36.1 | 50.7 | 42.9 | 41.2 |
| 125 | 34.9 | 41.5 | 39.8 | 36.4 |
| 150 | **29.0** | 37.4 | 38.2 | 38.0 |

Fits σ_t = a/√E ⊕ b: DSB1 a=223, b=26.5; LuAG a=455, b=11.6; MIXED a=224,
b=34.6; TENERGY a=299, b=28.4 (ps). **DSB1 lowest at every energy.** DSB1 @150 =
29.0 reproduces the official 27.4 within ~1.6 ps (method-consistency check). σ_E/E
flat ~13.5–15 % across builds (build-independent); DSB1 σ_E fit a≈52 %/√GeV matches
the published 52.04.

**Per-config systematic (cut variation on stable all-fid σ_t) [ps]:** DSB1 1.3,
LuAG 1.1, MIXED 1.7, TENERGY 1.5 (≈ the published 2.7 total; fiducial dominates).

**Caveats:** MIXED is thin at high E (N=16 k @150) → its OOS/systematic are noisier;
consider reducing more MIXED runs from Argon (190/252 analyzed). The light-yield
plot's single-channel σ_t is reference-limited (~150–280 ps, matching the paper's
per-channel values); DSB1 & LuAG **overlap with no species separation** at equal
light yield (the claim), but the sharp 1/√LY trend needs the reference-subtracted
intrinsic σ_t (pairwise method) — paper-grade refinement.

## 7b. Capillary TYPE finding — TENERGY contains an E-type energy cap (June 2026)

T-type = WLS at shower max (timing); E-type = WLS full-length (energy).
**TENERGY ("3×DSB1, 1×Energy") has one E-type cap at NW** — logbook ("3xDSB1, 1x
Energy"; "E,52,54,57"; "full LuAG Energy Cap … NW") + data (NW ⟨LG⟩≈110 mV vs
≈480 for the DSB1 caps; `configCapDiag.C`). Our (DW−UP)/2 was averaging it as a
timing channel. Excluding NW (`tenergyClean.C`) → clean 3×DSB1 timing **≈ DSB1**
(150 GeV: 37.7 contaminated → **34.9** excl-NW vs DSB1 35.6).

**Reorganization (decided):** the clean timing-material comparison is
**DSB1 / LuAG / MIXED** (all T-type). **TENERGY is dropped from the timing
headline**; its E-type cap belongs in the energy / future-work discussion.
DSB1 (4 T-type) and MIXED (clean 2 DSB1 + 2 LuAG T-type) verified clean from the
per-cap data. **OPEN — needs author confirmation:** LuAG config's NW capillary
type (4 T-type, or does it also carry the full-length LuAG energy cap?). This
gates whether the LuAG timing curve and the σ_t-vs-LY plot are fully clean.

## 8. Target

NIM A, RADiCAL Collaboration. Frames as the realization of arXiv:2401.01747
Section 7 item 2. Radiation hardness of LuAG:Ce is the motivation (ref Hu et al.).

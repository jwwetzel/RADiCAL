# Depth-Dial Internal Review (referee-style) — GATE 3 hardening

Date: 2026-06-09. Reviewed objects: `depthDial.C` result (slope −33.6 ± 2.9 ps/e-fold) +
`depthDialDiag.C` stability suite. Figures: `papers/figures/depth_dial/depth_dial.png`,
`depth_dial_diag.png`, `diagnostics/depth_dial_stability.png`. Logs: `depth_dial_result.log`,
`diagnostics/depth_dial_stability.log`. Pre-registration: `AUDIT.md` (written before plotting).

---
## 1. Observable definition

**What exactly is t_DW, t_UP?** For each event, t_DW = the mean of the valid downstream capillary-end
crossing times (kCap channels 0–3: NW-D, NE-D, SE-D, SW-D) and t_UP = the mean over upstream ends
(channels 4–7), where each end's crossing time comes from ONE named estimator (primary: `hg_lgcfd` =
srCFD, the CFD at 15% of the LG-predicted true amplitude). In this analysis ALL timing ends are
required valid (peak ≥ 20 mV, valid crossing, |t − event median| < 2 ns), so the means are always over
the same four ends each — the channel composition cannot drift with energy.

**Same-group or cross-group?** All eight HG timing channels are digitized in the SAME DT5742 module
(DRS0 @ 5 GS/s) with deliberate same-group wiring; t_diff = (t_DW − t_UP)/2 is a **same-board,
same-event difference of raw crossing times** — no MCP reference enters, no stop-cell or inter-group
correction is applied or needed. MCP1 appears in the event selection only (200–750 mV window).

**Arbitrary offset?** Yes — absolutely. t_diff carries the per-channel constants (cable lengths, SiPM
transit, channel skew) which set the absolute value (≈ −0.2 to −0.3 ns depending on estimator). These
constants are ENERGY-INDEPENDENT. The analysis therefore reports Δ⟨t_diff⟩ relative to the lowest
energy; the offset cancels exactly, and only the energy DEPENDENCE — the slope — is physical. This is
why the publication figure's y-axis is Δ⟨t_diff⟩, not ⟨t_diff⟩.

## 2. Sign convention

Shower max migrates DEEPER (downstream, toward larger z) as ln(E): t_max ≈ ln(E/E_c) + C in units of
X₀. Deeper emission point ⇒ shorter light path to the DOWNSTREAM ends (t_DW decreases by Δz/v_eff) and
longer path to the UPSTREAM ends (t_UP increases by Δz/v_eff). Hence
Δ[(t_DW − t_UP)/2] = −Δz/v_eff < 0: **the dial must move NEGATIVE with increasing energy.**
Measured: ⟨t_diff⟩ goes −195.7 → −263.0 ps from 25 → 150 GeV — monotonically negative at every step.
**The sign matches the expectation in the current convention.** (Cross-check of the convention itself:
channels 0–3 are verified D-ends in the config/channel map; the same convention produced the published
(DW−UP)/2 distributions.)

## 3. Magnitude check (order-of-magnitude derivation)

- Migration: Δz ≈ X₀ · ln(E/E₀), with X₀ = 5.4 mm (effective, this stack).
- Dial: Δ⟨t_diff⟩ ≈ −Δz/v_eff. With the naive axial group velocity in fused quartz at 495 nm
  (n_g ≈ 1.47 ⇒ v_g ≈ 205 mm/ns): predicted slope = −X₀/v_g = −5.4/0.205 ≈ **−26.3 ps per e-fold**;
  total over ln(150/25) = 1.79 e-folds ≈ −47 ps.
- Measured: **−33.6 ± 2.9 ps per e-fold** (total −67.4 ps). Ratio measured/naive = 1.28.
- Interpretation of the excess: meridional/skew ray paths in the capillary travel 1/cosθ longer than
  axial rays; the acceptance-weighted v_eff is generically BELOW the axial v_g (v_eff ≈ 160 mm/ns
  reproduces −33.6 exactly). A WLS re-emission isotropy + capture-cone calculation or the Stage-2
  optical simulation will pin v_eff.
- **Why this is meaningful but NOT a calibrated depth resolution:** the agreement establishes that the
  observable responds to longitudinal shower position at the right ORDER (ps-per-X₀ scale). It does
  not provide (i) a per-event z estimate, (ii) a ps→mm constant better than ~30% (v_eff unknown), or
  (iii) a depth RESOLUTION — the event-by-event spread of t_diff mixes photostatistics with genuine
  depth fluctuations in unknown proportion until the G4 truth regression exists. All three remain
  forbidden claims.

## 4. Method dependence

| estimator | slope (ps/e-fold) | interpretation | acceptable as depth observable? |
|---|---|---|---|
| **srCFD (lgcfd)** | **−33.6 ± 2.9** | constant fraction of RECOVERED (unclipped-equivalent) amplitude → amplitude- and clip-independent; tracks the light centroid | **YES — primary** |
| cfd05 | −41.6 ± 13.7 full-range; ≈ −23 restricted to 50–150 | fraction of RECORDED peak; the 25→50 GeV step is the clip-fraction turn-on (4%→71%) — pure time-walk; past it, concurs with srCFD | conditionally (unclipped or clip-stable regimes only) |
| LED (20 mV fixed) | −5.5 ± 1.5 | fixed threshold fires on the FIRST-arriving photons; the crossing does not track the emission centroid → the dial is suppressed ×6 | **NO** |
| cfd20 | −85.8 ± 23.4 | fraction of recorded (clipped) peak: crossing slides down the true edge as clip deepens with E — walk-dominated, not depth | NO (on clipped data) |
| cfd30 | −105.5 ± 28.6 | same mechanism, worse (higher fraction = deeper into the clip distortion) | NO (on clipped data) |

**The depth dial is a CFD-family observable — specifically a saturation-recovered-CFD observable on
clipped data — and is NOT an LED-threshold observable.** This is itself a physics result (the
centroid-vs-first-photon distinction) and fixes the estimator choice for any future depth use.

## 5. Stability checks (all run; `diagnostics/depth_dial_stability.{png,log}`)

| variant | slope (ps/e-fold) | verdict |
|---|---|---|
| baseline: r<2.5, all 6 E | −33.6 ± 2.9 | reference |
| drop 150 GeV | −33.1 ± 4.1 | stable (the 125→150 re-steepening does not drive the slope) |
| drop 125 + 150 GeV | −36.4 ± 4.6 | stable |
| tight fiducial r<1.5 | −35.3 ± 5.6 | stable |
| per-capillary NW / NE / SE / SW | −41.6±15.5 / −33.1±18.4 / −45.6±28.1 / −46.9±14.4 | all four negative & mutually consistent — not a single-channel artifact |
| per-run (26–58 runs per energy, N>700) | run-to-run spread 4.8–8.8 ps at every energy | not a run-period artifact; the 125→150 local step (−10 ps) is comparable to the run spread (8.8 ps @150) |
| brightness quartiles Q1→Q4 (dim→bright) | −94.7 / −45.1 / −27.9 / −17.7 | **systematic gradient — see below** |

**The brightness-quartile gradient (the one genuinely conditional finding).** At fixed beam energy the
light acceptance of the fixed-position, ~15 mm WLS filament depends on where the shower deposits —
brightness selection therefore SCULPTS the depth distribution. The bright quartile is pinned to the
filament sweet spot (slow drift, −17.7); the dim quartile increasingly samples the migrating deep tail
(fast drift, −94.7). Three consequences:
1. This is CORROBORATING evidence for the depth mechanism (an electronics drift would not know about
   per-event brightness at fixed E), and is the population-level shadow of the per-event
   depth–amplitude correlation (the next unlock step).
2. **The quoted slope is acceptance-conditional**: −33.6 is the full-fiducial-population value;
   a brightness-selected sample would measure −18 to −45. Manuscript language must carry the
   "full-fiducial acceptance" qualifier.
3. Paper-2 relevance: brightest-K timing samples have a depth-biased MEAN — harmless for σ (a width),
   but t_diff means from selected samples must never be compared across selections.

**Chi-square honesty:** every fit has χ²/ndf ≫ 1 against 0.2–0.4 ps statistical errors — the drift is
not exactly logarithmic (local slope −47 → −19 → −56 ps/e-fold across the scan). The pre-registered
filament-edge compression explains the mid-range flattening qualitatively; the slope uncertainties
quoted everywhere are already inflated by √(χ²/ndf). The log-linear fit is a one-parameter summary of
a monotonic curve, not a shape claim.

## 6. Verdict

**STABLE in sign, monotonicity, and order of magnitude** across energy subsets, fiducial radius,
individual capillaries, and run periods. **CONDITIONAL in the numerical slope**, which is an
acceptance-weighted population quantity (quartile gradient ×5) with a ~30% scale uncertainty from
v_eff. **GATE 3 therefore remains PASSED at the structural level** — the depth dial exists, with the
predicted sign and scale — while the slope VALUE is downgraded to "conditional on the stated
(full-fiducial) acceptance" pending the per-event corroboration and the G4 v_eff/truth calibration.
No claim beyond the three-tier language in `papers/memory_claims_and_forbidden_language.md` is
licensed by this result.

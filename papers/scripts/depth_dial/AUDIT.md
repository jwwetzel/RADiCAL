# Depth-Dial Audit Note (pre-analysis) — GATE 3

Date: 2026-06-09. Author: analysis session (pre-registered before plotting).

## Observable and prediction
Depth dial = energy dependence of the MEAN of t_diff ≡ (t_DW − t_UP)/2.
Mechanism: shower max migrates downstream ~1 X₀ per e-fold of E (t_max ≈ ln(E/E_c)+C). The light-emission
centroid along the capillary moves with it; propagation to the DW end shortens while UP lengthens, so
⟨t_diff⟩ shifts by −Δz/v_g per e-fold.
Prediction: slope = −X₀/v_g ≈ −5.4 mm / 205 mm·ns⁻¹ ≈ **−26 ps per e-fold**, total ≈ −47 ps over
25→150 GeV (ln 6 = 1.79 e-folds). PASS band (pre-registered): −(26 ± 10) ps/e-fold, negative sign required.

## Branch identification (verified by direct inspection, 2026-06-09 — not assumed)
All four builds' reduced 'rad' trees carry the identical 39-branch schema (dump in session log; summary
in `papers/memory_dataset_inventory.md`):
- Downstream end times: `hg_lgcfd[0..3]`, `hg_led[0..3]`, `hg_cfd05[0..3]` (ns; kCap order NW,NE,SE,SW-D)
- Upstream end times: same arrays, indices 4..7 (NW,NE,SE,SW-U)
- MCP reference: `mcp1_time`, `mcp1_peak` (selection only — t_diff is reference-free by construction)
- Energy: `beam_energy` (per-file constant; files are per-energy) · light: `sum_lg`
- Quality: `wc_ok`, `x_trk`,`y_trk` (fiducial), `hg_peak[8]` (channel validity ≥ 20 mV),
  `hg_saturated[8]` (clip flags), timing-channel mask via BuildConfig/RadView `is_timing()`
  (TENERGY has only 6 timing ends — 3 T-type capillaries; the E-type capillary is not a timing channel).

## Pre-registered method choices (and why)
1. **Source = lgcfd (srCFD), primary.** A constant-fraction discriminator is amplitude-independent, so
   the E-dependent amplitude growth cannot fake a drift. LED (fixed 20 mV) walks earlier as amplitude
   grows — and the DW/UP amplitude asymmetry is itself depth-dependent — so LED drift conflates walk
   with depth. cfd05 on clipped pulses time-walks with clip fraction. Both are shown as CROSS-CHECKS in
   the diagnostic figure; method-consistency of the slope is part of the pass criterion.
2. **Fixed fiducial r = 2.5 mm at ALL energies.** The production 2.5→3.0 mm step at 112 GeV changes the
   accepted shower population exactly between 100 and 125 GeV — a potential fake step in ⟨t_diff⟩.
3. **Fixed channel composition: require ALL timing ends valid** (is_timing, hg_peak ≥ 20 mV, in-event
   |t−median| < 2 ns veto). Per-channel cable/offset constants differ; if the contributing-channel mix
   shifted with E, ⟨t_diff⟩ would shift without any physics.
4. **Full-fiducial events, no brightness selection.** Brightness correlates with shower depth at fixed E;
   a brightest-K cut could sculpt the depth distribution.
5. **Robust (2.5σ-truncated, iterated) mean**; uncertainty = truncated σ/√N. N ~ 10⁴–10⁵ per point →
   sub-ps statistical errors; systematics dominate (hence the method cross-checks).
6. **Report Δ⟨t_diff⟩ relative to the lowest energy** (25 GeV for DSB1; 50 GeV otherwise). The absolute
   offset encodes cable lengths and channel constants — physically meaningless.
7. Fit: Δ⟨t_diff⟩ = p₁ · ln(E/E_ref); p₁ in ps per e-fold, compared to −X₀/v_g.

## Known caveats (pre-registered)
- The WLS filament is ~15 mm long, positioned for FTBF-era 20–30 GeV shower max (parent Fig 7 caption:
  position "adequate (although not optimized)" for high-E showers at layers 11–13). Total migration over
  25→150 GeV ≈ 9.7 mm is comparable to the filament length → the response may COMPRESS at the high-E end
  (emission centroid clamped by the filament edge). A slope magnitude somewhat below 26 ps/e-fold,
  or curvature at high E, is therefore expected physics, not failure — sign and order of magnitude decide.
- v_g in the fused-quartz core at 495 nm: n_g ≈ 1.47 → v_g ≈ 204 mm/ns; we quote 205 ± 5 mm/ns.
- WLS re-emission delay (τ = 3.5 ns DSB1) is common to both ends at fixed emission point — cancels in
  t_diff to first order; only the propagation asymmetry survives.
- MIXED: per-material corner subsets only (DSB1 diag lgcfd, LuAG diag led) — module-wide mixing is
  ill-posed (established this session). LuAG:Ce τ ~ 70 ns slows edges; its slope is a weak cross-check.

## Inputs
`data/2023/reduced/DSB1/{25,50,75,100,125,150}GeV.root` (primary);
`{LUAG,MIXED,TENERGY}/{50..150}GeV.root` (cross-checks). Script: `depthDial.C` (this directory).
Outputs: `papers/figures/depth_dial/depth_dial.png` (publication), `depth_dial_diag.png` (diagnostics),
slope table on stdout (logged to `papers/scripts/depth_dial/depth_dial_result.log`).

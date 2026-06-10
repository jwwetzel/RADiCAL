# Analysis Gates — ranked, specified, auditable

Updated: 2026-06-09. Each gate: inputs → branches → script → output → pass/fail → manuscript consequence.
Ranking = publication leverage × risk reduction. Status legend: OPEN / IN PROGRESS / PASSED / FAILED.

Common inputs: reduced trees `data/2023/reduced/<BUILD>/<E>GeV.root` ('rad' tree, 39 branches — see
`memory_dataset_inventory.md`); configs `data/2023/configs/<BUILD>.json`; helpers `lib/` (RadView,
RadTiming post-fix, SelectionCuts, BuildConfig).

---
## GATE 3 — Depth dial: mean (t_DW − t_UP)/2 vs ln(E)   [RANK 1 — cheapest, highest leverage]
STATUS: **PASSED (structural), slope value CONDITIONAL on acceptance** — hardened 2026-06-09 by the
referee-style review `papers/scripts/depth_dial/DEPTH_DIAL_REVIEW.md` + stability suite
(`depthDialDiag.C` → `papers/figures/depth_dial/diagnostics/depth_dial_stability.{png,log}`).
HARDENING SUMMARY: slope stable under drop-150 (−33.1±4.1), drop-125+150 (−36.4±4.6), tight fiducial
r<1.5 (−35.3±5.6); all four capillaries independently negative (−33 to −47); run-to-run spread
4.8–8.8 ps at all six energies (26–58 runs each) so the 125→150 local step (−10 ps) is within run
scatter; BUT brightness quartiles show a ×5 slope gradient (Q1 dim −94.7 → Q4 bright −17.7) —
depth-dependent light acceptance sculpts the depth distribution, so the quoted slope is an
acceptance-weighted population quantity (this is also corroborating evidence FOR the depth mechanism
and the seed of the per-event unlock). Extended method table: cfd20/cfd30 walk-dominated on clipped
data (−86/−105) — the dial is specifically a SATURATION-RECOVERED-CFD observable.

### RESULT (2026-06-09)
- Primary (DSB1, srCFD/lgcfd, fixed r=2.5 mm, all-ends-valid, full fiducial, robust mean):
  ⟨t_diff⟩ = −195.7 → −263.0 ps over 25→150 GeV (N = 27k–90k/point, stat err 0.2–0.4 ps).
  **Slope −33.6 ± 2.9 ps per e-fold** (err inflated by √(χ²/ndf)); predicted −X₀/v_g = −26.3;
  pre-registered PASS band −(26±10) → **PASS** (negative sign ✓, magnitude ✓, monotonic ✓).
- Curvature beyond log-linear (χ²/ndf=972/4 against sub-ps errors): local slope −47 (25→50) →
  −19 to −20 (75→125) → −56 (125→150) ps/e-fold. Mid-range compression was PRE-REGISTERED
  (15 mm filament edge vs ~10 mm total migration); the 125→150 re-steepening is UNEXPLAINED — flag
  for run-period check before publication.
- Slope excess vs naive prediction (×1.28) consistent with modal dispersion lowering effective v_g
  (~160 mm/ns vs axial 205); Stage-2 optical sim to pin (then the ps→mm calibration constant).
- METHOD DEPENDENCE (new physics input for P3 S8): the dial is read by the CFD family only —
  cfd05 mid-range (50–150, past its clip turn-on) ≈ −23 ps/e-fold, consistent with srCFD; the LED
  family is suppressed (DSB1 led −5.5±1.5, TENERGY led −2.3±3.3): fixed-threshold crossings fire on
  the first-arriving photons and do not track the emission centroid. LuAG led −13.9±1.1 (τ~70 ns).
  ⇒ the depth estimator MUST be CFD-based (srCFD).
- RMS panel: full-fiducial t_diff RMS 50–62 ps ≈ photostat ⊕ depth (1 X₀ ↔ 26 ps) — direct tie to
  the Paper-2 floor narrative.
- Products: `papers/figures/depth_dial/depth_dial.png` (+`_diag.png`),
  `papers/scripts/depth_dial/{depthDial.C, AUDIT.md, depth_dial_result.log}`.
- Consequences taken: P2 S9 drift panel AVAILABLE; P3 S8 fifth-coordinate section UNLOCKED at
  population level; **"toward 5D" language permitted** with population-level scoping (see claims memory).
- Remaining within this gate (the per-event step, panel unlock #2): t_diff vs ln(A_DW/A_UP) and
  t_diff vs Pb-glass leakage at fixed E; plus the 125→150 run-period check.
- Inputs: DSB1 reduced 25–150 GeV (6 energies; the only full scan); LUAG/MIXED/TENERGY 50–150 as
  cross-checks. Branches: `hg_lgcfd[8]` (primary source; CFD = amplitude-independent → immune to the
  LED time-walk confound), `hg_led[8]` + `hg_cfd05[8]` (cross-checks), `hg_peak[8]`, `wc_ok`,
  `x_trk`,`y_trk`, `mcp1_peak`, `sum_lg`. DW = ch 0–3, UP = ch 4–7 (kCap universal map).
- Script: `papers/scripts/depth_dial/depthDial.C` (NEW).
- Output: `papers/figures/depth_dial/depth_dial.png` (publication) + `depth_dial_diag.png` (diagnostic)
  + printed slope table.
- Method guards: FIXED fiducial r=2.5 mm at ALL energies (kill the 2.5→3.0 seam); require ALL 8 ends
  valid (fixed channel composition vs E — composition drift fakes depth drift); full-fiducial events
  (no brightness selection — brightness correlates with depth); robust (truncated) mean; plot
  Δ⟨t_diff⟩ relative to 25 GeV (absolute offset = cable/channel constants, physically meaningless).
- PASS: ⟨t_diff⟩ drifts NEGATIVE (deeper→downstream-earlier) and linearly in ln(E) with slope
  ≈ −X₀/v_fiber ≈ −(26±10) ps per e-fold (X₀=5.4 mm, v_g≈205 mm/ns), total ≈ −48 ps over 25→150;
  lgcfd and cfd05 slopes agree (method-independence); per-build consistency.
- FAIL: no drift, wrong sign, or strongly method-dependent slope.
- If PASS → P2 S9 gains the drift consistency panel ("the floor IS depth physics"); P3 S8 becomes the
  measured fifth coordinate; "toward 5D" language unlocked (with G4 calibration to mm).
- If FAIL → floor decomposition stays simulation-only; "toward 5D" forbidden; P3 S8 reduced to an
  amplitude-asymmetry-only discussion.

## GATE 1 — MIXED corner-map / logbook confirmation   [RANK 2 — gates the kill-shot; USER ACTION]
STATUS: **DATA HALF CONFIRMED 2026-06-09 (8/8); logbook half still open (user action).**
DATA-HALF RESULT: brightness-independent pulse-shape classification
(`papers/scripts/mixed_corner_map/{AUDIT.md, cornerDiscriminant.C, corner_discriminant.log}`,
figure `papers/figures/mixed_corner_map/corner_discriminant.png`): all 8 MIXED ends classified
against MEASURED templates from the labeled pure builds (D1 = lg_charge/lg_peak: DSB1 51.2±1.9 vs
LuAG 58.5±3.7; D2 = cfd50−cfd05 on unclipped HG: 2.02±0.36 vs 1.65±0.31 ns; D3 weak, ±0.2 only).
**8/8 ends match NE+SW = DSB1, NW+SE = LuAG; all four capillary D/U pairs internally consistent;
D1 and D2 agree in sign for every end.** Margins 1.1–2.6 in template-width units (solid, not
overwhelming) — the logbook remains the belt-and-suspenders confirmation, but the kill-shot no longer
rests on a brightness-circular inference. Manuscript provenance sentence available: "the corner
assignment is confirmed by amplitude-independent pulse-shape discriminants validated on the
single-material modules."
- Inputs: 2023 beam-test logbook (assembly records, capillary serials); cross-check data: MIXED reduced,
  branches `hg_peak[8]`, `hg_tot[8]`, `hg_charge[8]`/`hg_peak` ratio (pulse-shape discriminant:
  DSB1 τ=3.5 ns vs LuAG:Ce τ~70 ns ⇒ charge/peak differs strongly).
- Script: `analyze/studies/cornerDiscriminant.C` (TO CREATE — independent pulse-shape ID per corner).
- Output: logbook scan/quote + discriminant figure (charge/peak per corner, bimodal by material).
- PASS: logbook AND discriminant both give NE+SW=DSB1, NW+SE=LuAG.
- If PASS → kill-shot publishable with provenance sentence. If logbook unavailable → the pulse-shape
  discriminant alone (NON-brightness-based, breaking the circularity) + explicit caveat.
- If FAIL (map wrong) → kill-shot re-derived with corrected map (analysis is map-parametric; rerun
  `mixedSeparate.C` + h2h with swapped indices).

## GATE 2 — Position arithmetic reconciliation   [RANK 3 — blocks every P3 mm claim]
STATUS: OPEN.
- Inputs: DSB1+TENERGY reduced all E. Branches: `lg_peak[8]` (4-corner light division), `x_trk`,`y_trk`,
  `wc_ok`, `wc_peak[4]`, `in_fiducial`.
- Script: `papers/scripts/position/positionReconcile.C` (TO CREATE): unbinned event residuals (estimator −
  WC), quadrature decomposition σ²_resid = σ²_est ⊕ σ²_WC, WC σ re-derivation (3.6 vs 3.3 vs ~1 mm
  effective?), TLinearFitter train/test run-fold split, tracker-free split-estimator internal σ_x
  (disjoint corner subsets, per-event difference /√2).
- Output: reconciliation note + unbinned residual figure + decomposition table.
- PASS: a self-consistent picture (e.g. residual ≥ WC term, or WC σ re-derived smaller than 3.6).
- If PASS → S6 written with "comparator-limited upper bound" + the internal split-estimator number.
- If FAIL (irreconcilable) → drop mm numbers; S6 becomes light-division demonstration only.

## GATE 4 — t/z separation for the 4D capstone (T_abs)   [RANK 4]
STATUS: OPEN.
- Inputs: TENERGY reduced (+DSB1 check). Branches: `hg_led[8]`, `mcp1_time`, `mcp1_peak`, `stopcell[4]`
  (exists in reduced — reducer-level fix may NOT be needed; verify which DRS group each channel maps to
  via `data/2023/configs/*.json` channel_map before trusting stopcell indexing).
- Script: `papers/scripts/fourD/tAbs.C` (TO CREATE): T_abs = mean(hg_led DW+UP) − mcp1_time, stop-cell-
  corrected; σ(T_abs) vs E; compare same-group vs cross-group reference.
- PASS: σ(T_abs) ~50 ps class, stop-cell correction reduces the 71–74 ps inter-group term toward ~15 ps.
- If PASS → 4D capstone time axis = T_abs; t_diff freed as the depth observable (5D narrative coherent).
- If FAIL → documented fallback: dual-label wording on the existing 4D figure (panel-approved fallback).

## GATE 6 — Kill-shot ratio with bootstrap CI   [RANK 5]
STATUS: OPEN.
- Inputs: MIXED reduced 50–150. Branches: per-end timing (adopted per material: `hg_lgcfd` DSB1 corners /
  `hg_led` LuAG corners), `hg_peak[8]`, selection branches.
- Script: extend `analyze/studies/mixedSeparate.C` or new `mixedH2hBootstrap.C`: paired event-bootstrap
  (resample events, recompute per-material σ and ratio), run-jackknife, audit the χ²/ndf=0.4/5.
- Output: ratio ± CI; updated `mixed_h2h.png` lower panel.
- PASS: CI excludes large deviations from 1 (e.g. 0.99 ± few %) and jackknife stable across runs.
- If PASS → S8 claim ladder complete. If FAIL (run-dependent) → ratio quoted per run-period with spread.

## GATE 5 — Universal collapse: σ_t vs detected light   [RANK 6 — the graphical abstract]
STATUS: OPEN (data exists; figure to build).
- Inputs: all four builds, all E. Branches: `sum_lg` (light proxy), adopted-source timing.
- Script: `analyze/studies/universalCollapse.C` (TO CREATE): σ_t (brightest-K slices AND per-E points)
  vs ⟨sum_lg⟩ of the slice, all builds on one panel; fit σ = p0/√L ⊕ p1.
- PASS: single curve within errors (the thesis in one panel). FAIL: build-dependent offsets → report as
  per-build photostatistics with the Npe budget instead.

## GATE 7 — Energy + position response figures (P3 backbone)   [RANK 7]
STATUS: OPEN.
- Inputs: all builds. Branches: `lg_peak[8]`, `sum_lg`, `sum_pb`, `pb_peak[4]`, `beam_energy`, WC.
- Script: `papers/scripts/energy/energyBackbone.C` (TO CREATE): linearity + residuals; σ_E_SM/E(E) per
  build; convention reconciliation (29.4 vs ~55 mV/GeV); Eq-1 decomposition (noise from `hg_ped_rms`/
  summed LG pedestals; leakage from `sum_pb`); E-type slopes from TENERGY.
- PASS: numbers reproduce report-era values with post-fix conventions. Consequence: S4–S5 written.

---
## Rank summary
1. GATE 3 depth dial (in progress) → 2. GATE 1 logbook (user) → 3. GATE 2 position arithmetic →
4. GATE 4 T_abs → 5. GATE 6 bootstrap CI → 6. GATE 5 collapse → 7. GATE 7 energy backbone.

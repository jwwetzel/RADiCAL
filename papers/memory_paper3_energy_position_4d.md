# Paper 3 — Energy + Position + 4D capstone (working memory)

Updated: 2026-06-09. Manuscript: `papers/energy_position/radical_energy_position.tex` (skeleton).
Blueprint: `papers/EXPERT_PANEL_BLUEPRINT.md` §Paper 3; JSON `paper3_blueprint`.

## Title (panel-proposed)
"Shower-maximum energy response, transverse shower localization, and a simultaneous four-dimensional
(x, y, E, t) measurement in an ultra-compact W/LYSO:Ce (RADiCAL) calorimeter module"

## Thesis
The SAME four corner-capillary shower-max signals deliver three coordinates at once — E_SM (linear to
150 GeV, σ/E ≈ 14% @150, WLS-species-independent, scoped as sampling), intra-cell transverse position
(comparator-limited residual vs a ~3.6 mm wire chamber — a regime no shashlik occupies), and with T-type
dual-end timing a simultaneous 4D (x, y, E_SM, t) demonstration in one TENERGY event sample. Plus the
series hinge: **(t_DW − t_UP)/2 is a longitudinal depth observable** whose mean must track the ln(E)
shower-max migration (the depth dial) and whose spread IS Paper 2's floor.

## Known results (verify before writing — most are report-era, pre-fix)
- E_SM linearity all builds; parent normalization 29.4 mV/GeV vs our LG-sum convention ~55 mV/GeV —
  RECONCILE IN PRINT (B-list). σ_E_SM/E ≈ 13.5–15% @150 across builds (recorded; re-derive).
- E-type (TENERGY): slopes ~1.6 vs 6.7 mV/GeV (~4–5× dimmer than T-type), linear 50–150 GeV (recorded).
- Position — **GATE 2 RECONCILED 2026-06-09** (was: "1.5 vs 3.6 impossible" — the 3.6 was the
  t₀-inflated end-time-SUM bound, never the difference-mode comparator; the memory had it wrong).
  Re-derived, current schema, train/test held-out: x/y residual RMS **1.54/1.45 mm** (full ±6 mm
  window, 150 GeV; Gaussian cores 1.26/1.16), **0.91/0.88 mm in the r<2.5 beam core**, roughly
  energy-independent (25–150 GeV). Closure slope 0.698 → linear estimator saturates beyond |x|≈3 mm;
  S6 must restrict to the core or use a nonlinear estimator. Quote ONLY as a joint upper bound
  (capillary AND tracker); no intrinsic resolution. Products:
  `papers/scripts/position_reconciliation/`, `papers/figures/position_reconciliation/`.
- Co-registration: module center (6.58, 4.66) mm stable ±0.15; MCP–RADiCAL alignment 0.79 mm; beam spot
  2.9 mm (recorded in report-era memory; re-verify when writing S2).
- Containment: SumPb < 0.30·SumLG cut, 95.3% eff @150; punch-through 6.2% @150 (recorded).

## Section outline (S1–S12 condensed)
S1 promise-kept intro (§7.1) · S2 co-registration frame · S3 selection/containment (4 populations) ·
S4 E_SM (linearity+residuals, convention reconciliation, multi-build σ_E_SM/E, Eq-1 decomposition) ·
S5 E-type in-beam characterization (the measured closure of parent Fig 4) · S6 transverse localization
(train/test folds, unbinned residuals, WC band, tracker-free split-estimator inset) · S7 pattern
recognition (module-ID demo + two-shower event mixing, π0/γ reach labeled extrapolation) ·
S8 **the fifth coordinate** (depth dial: mean t_diff in mm vs lnE over G4 migration + per-event
two-estimator correlation) · S9 4D capstone (T_abs time axis, t_diff relabeled depth, E_SM, x/y;
correlation matrix; all five energies) · S10 shashlik-heritage table · S11 systematics/OOS · S12 conclusions.

## Money figures
1. Energy backbone 2-panel (linearity + σ_E_SM/E all builds, leakage inset). 2. Localization (estimator
vs WC truth + unbinned residual + WC band + split-estimator inset). 3. **THE DEPTH DIAL** (series money
figure — GATE 3): mean t_diff [mm] vs ln(E) over the G4 shower-max migration + per-event t_diff vs
ln(A_DW/A_UP) companion. 4. 4D capstone rebuilt (T_abs, t_diff-as-depth, E_SM, x/y + correlation matrix,
multi-energy table). 5. Graphical abstract: the 4D shower passport (one real 150 GeV event). 6. Pattern
recognition (module-ID map + two-shower efficiency). 7. E-type 2-panel + G4 kernel fold.

## Unresolved gates blocking this paper
- GATE 2 — **PASSED/RECONCILED 2026-06-09** (see Known results above + memory_analysis_gates.md).
  S6 unblocked with the joint-upper-bound language; remaining S6 work: nonlinear/log-weighted
  estimator for the full window, σ_WC(diff) direct measurement for any future unfolding.
- GATE 3 depth dial — **PASSED 2026-06-09 at population level**: slope −33.6±2.9 ps/e-fold vs
  predicted −26.3 (DSB1 srCFD, 25–150 GeV); S8 unlocked; depth estimator must be CFD-based (LED
  suppressed, −5.5). Per-event corroboration (t_diff vs ln(A_DW/A_UP), vs Pb-glass leakage) still
  open; 125→150 re-steepening needs a run-period check. Products: `papers/figures/depth_dial/`,
  `papers/scripts/depth_dial/`. Details: memory_analysis_gates.md GATE 3 RESULT.
- GATE 4 t/z separation (T_abs via stop-cell-corrected mean(hg_led) vs same-group MCP) — blocks the S9
  4D time axis (fallback: dual-label wording, documented).
- GATE 7 energy/position figure re-derivation with one script + post-fix conventions.

## Required simulation (Stage 1 shared w/ P2)
Sampling-fraction calibration in GeV (what the ~3-tile/1-R_M slice samples) + leakage attribution
(transverse ~1 R_M dominant) + 4-corner light-sharing transfer function → predicted intrinsic σ_x,y(E)
+ E-type kernel fold + depth-truth regression (ps-to-mm calibration; instrumental-vs-fluctuation split).

## Must-cite
Parent §5.1.3/§7.1 · position-sensitive shashlik heritage: KOPIO/Atoian NIM A 531 (2004), PHENIX,
LHCb/HERA-B, COMPASS · Awes et al. NIM A 311 (1992) (log-weighted centroid) · wire-chamber/Spanggaard
DWC documentation · HGCAL TDR + SCEPCAL (when using any 4D/5D adjacency) · G4 references as in P2.

## Next actions
1. GATE 2 reconciliation (unbinned residual, quadrature decomposition, WC σ re-derivation, fold split).
2. GATE 3 depth dial (this session).
3. GATE 4 T_abs reducer fix (stop-cell correction available: `stopcell[4]` branch exists in reduced data).
4. Single energy/position re-derivation script → S4–S6 figures with post-fix conventions.
5. 4D capstone rebuild at all five TENERGY energies.

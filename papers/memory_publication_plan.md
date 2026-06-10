# Publication Plan — RADiCAL Papers 2 & 3 (master memory)

Updated: 2026-06-09. Source of truth for the two-paper campaign. Derived from the ten-expert panel
(synthesis: `papers/radical_expert_panel_synthesis_2026-06-09.json`, individual reports:
`papers/radical_expert_panel_reports_2026-06-09.json`, prose: `papers/EXPERT_PANEL_BLUEPRINT.md`).
Companion memories: `memory_paper2_timing.md`, `memory_paper3_energy_position_4d.md`,
`memory_analysis_gates.md`, `memory_claims_and_forbidden_language.md`, `memory_dataset_inventory.md`.

## Series framing (locked)
- Parent: Perez-Lara et al., NIM A 1068 (2024) 169737. Its **Section 7 promises our follow-ups**:
  §7.1 → Paper 3 (shower localization vs wire chamber); §7.2 → Paper 2 (DSB1 vs LuAG:Ce, 25 GeV steps).
- Both papers open as "the study announced in Sec. 7 of [parent]" — promise-kept framing (RADiCAL II/III).
- Shared assets across the pair: nomenclature Table 1 (retire BestMinus/BestPlus → t_diff/t_sum + srCFD/LED),
  one σ-extraction convention (post-fix tebSigma), one Stage-1 GEANT4 campaign serving both papers,
  the depth analysis run once with the reading split (ps-consistency in P2, mm-calibrated 5th coordinate in P3).

## Current status (2026-06-09)
- Expert panel synthesis complete and preserved in-repo (this directory).
- Estimator bug found and FIXED in production (`lib/physics/RadTiming.h` robust window + in-event veto
  + tail guard; `lib/physics/SelectionCuts.h` continuous fiducial). All adopted-method σ_t(E) monotonic;
  published cfd05 numbers transparent to the fix (Δ ≤ 0.4 ps). Evidence:
  `figures/2023/narrative/monotonicity_fix.png`, `monotonicity_evidence.png`, `analyze/studies/timingRegression.C`.
- Full chain-of-evidence per build (waveforms → reduction → fits → σ_t(E)): `chain_of_evidence.html`.
- 8-source × 4-build timing survey runs in production: `analyze/studies/timingAllMethods.C`.
- MIXED resolved: per-material diagonals (no single module-wide method); same-shower width parity
  10–20% (srCFD 1.04±0.05, GATE 6 — the old 0.99 is retired; figure:
  `papers/figures/mixed_killshot_bootstrap/mixed_h2h_corrected.png`, plus `mixed_separate.png`).
- Manuscript skeletons exist: `papers/timing/radical_timing.tex` (NEEDS reframe: currently DSB1-solo-method
  center of gravity; must become the materials comparison), `papers/energy_position/radical_energy_position.tex`.
- NOTHING from this campaign is committed to git yet (lib fix + macros + figures all uncommitted).

## Gated work plan — claims by category

### A. Ready-to-write (supported by existing, verified analysis)
1. Light-yield ordering of the stochastic term: DSB1 a≈201 vs LuAG a≈455 (×2.3), brightest-1000 (DW−UP)/2,
   adopted sources. Support: `analyze/studies/methodCompare.C`, `figures/2023/narrative/method_compare_*.png`,
   `papers/timing/tab_methods.tex`.
2. Adopted-source rule by light regime (clipped→srCFD/lgcfd, dim→LED) + the 8-source survey.
   Support: `analyze/studies/timingAllMethods.C`, `figures/2023/narrative/timing_allmethods_*.png`.
3. The σ-estimator fix narrative (robust window, in-event veto, kurtosis 888→0.6, monotonicity restored,
   headline preserved). Support: `monotonicity_fix.png`, `monotonicity_evidence.png`, `sigmaProbe.C`, `outlierPeek.C`.
4. The saturation trap (HG-rank selects most-clipped events → σ rises; LG-sum ranking correct).
   Support: `monotonicity_evidence.png` panel (a).
5. cfd05 clipping limitation on brightest events = the documented motivation for srCFD.
   Support: `method_compare_DSB1.png`, regression tables in `output/logs/`.
6. Per-point fit transparency (every σ has a clean Gaussian behind it). Support: `method_dist_*.png`.
7. MIXED requires per-material treatment; heterogeneous module-wide averaging is ill-posed (clip-fraction
   table measured: DSB1 corners 33–55% clipped, LuAG corners 0–24%). Support: `mixedSeparate.C`, session logs.
8. (Paper 3) E-type linearity & T/E slope comparison from TENERGY reduced data — analysis exists in report;
   re-verify numbers before writing (B-list item 7 covers the re-verification).

### B. Needs verification (plausible; check scripts/logbooks/branches/calibrations BEFORE writing)
1. **MIXED corner material map** (NE+SW=DSB1, NW+SE=LuAG) — logbook confirmation + independent
   pulse-shape/decay discriminant. GATE 1 (see memory_analysis_gates.md).
2. ~~Kill-shot ratio 0.99 with defensible CI~~ **GATE 6 CONDITIONAL: 0.99 RETIRED** (brightness-label
   artifact — SE-D misassigned; χ²=0.4 was over-conservative errors). Replacement claim (B→A once
   `mixed_h2h.png` is regenerated): same-shower width parity within ~10–20%, srCFD 1.04±0.05.
3. ~~Position residual 1.5 mm vs WC σ≈3.6 mm arithmetic~~ **GATE 2 RECONCILED 2026-06-09**: the 3.6
   was the t₀-inflated sum-side bound (inapplicable); the residual is real, survives train/test
   (Δ=3.5 µm), and is 0.9 mm in the beam core / 1.5 mm full-window (saturation beyond |x|≈3 mm,
   closure slope 0.70). Joint-upper-bound language only — see claims memory.
4. TENERGY floor 240/22±3 and MIXED 233/35±2 — recompute with the POST-FIX estimator before quoting.
5. srCFD validity on 94%-clipped data — closure test on unclipped events (25 GeV / dim builds).
6. Satellite-removal (0.2 ns artifact, parent Fig 19) — currently asserted; needs the same-event
   cfd05-vs-srCFD before/after demonstration.
7. Paper-3 energy numbers (σ_E_SM/E ~14% @150, linearity slopes, 29.4 vs ~55 mV/GeV convention) —
   reconcile conventions and re-derive from reduced data with one script.
8. MCP1−MCP2 71–74 ps attribution (inter-group domino phase, not MCP intrinsic) — re-verify with
   stop-cell correction before printing the ~15 ps corrected number.

### C. Forbidden until unlocked (must NOT appear unless the named analysis succeeds)
See `memory_claims_and_forbidden_language.md` for exact wording. Summary:
- "5D" as achieved → requires the full depth-unlock chain (GATES 3+4 + G4 regression).
- "toward 5D" → requires GATE 3 (depth dial) to pass.
- Position *resolution* in mm → requires GATE 2 reconciliation (until then: "comparator-limited residual").
- Absolute time-of-arrival from t_diff → requires GATE 4 (t/z separation, T_abs vs same-group MCP).
- "Material-independent ~20 ps floor" as measurement → stays an extrapolation statement (LuAG best
  measured point is 42 ps).
- Any "revises 17.5 ps" wording → permanently forbidden ("confirms/consistent with" only).

## Next actions (ordered)
1. ~~GATE 3 depth dial~~ **PASSED 2026-06-09** (−33.6±2.9 ps/e-fold vs −26.3 predicted; "toward 5D"
   unlocked at population level; depth readout must be srCFD — see memory_analysis_gates.md).
   Remaining inside GATE 3: per-event corroboration (t_diff vs ln(A_DW/A_UP), vs sum_pb) +
   125→150 run-period check.
2. GATE 1 — **data half CONFIRMED 2026-06-09 (8/8 ends match the map, brightness-independent)**;
   logbook confirmation remains the user-action half (belt and suspenders, no longer blocking).
3. ~~GATE 2 position arithmetic reconciliation~~ **PASSED/RECONCILED 2026-06-09**.
4. ~~GATE 6 kill-shot bootstrap CI~~ **CONDITIONAL 2026-06-09 — the 0.99 is RETIRED** (brightness-label
   artifact; SE-D misassigned). Kill-shot rewritten as within-~20% same-shower width parity
   (srCFD 1.04±0.05 scatter; methods 0.80–1.16) vs the >×2 cross-build a-ratio. Corrected figure
   SHIPPED 2026-06-09 (`mixed_h2h_corrected.{png,pdf}`); old figure DEPRECATED. Remaining
   follow-ups: method-gain bootstrap CIs, position-coupling dilution bound (G4 or amp-binned pairs).
5. ~~Recompute all four builds' (a,b)~~ **DONE 2026-06-09**: `papers/tables/timing_fit_summary_2026-06-09.md`
   — DSB1 203±6/18.8±0.8 (σ150=25.7±0.6), LuAG 440±18/24.6±3.3, TENERGY 198±21/26.0±1.8 (√(4/3)
   penalty), MIXED 253±29/34.3±2.4 (ill-posed ref). a-ratio 2.16; floor narrative SOFTENED (19–26 ps
   span, see table addendum). Corrected kill-shot figure shipped: `mixed_h2h_corrected.{png,pdf}`;
   old figure deprecated. ← Paper-2 numbers are now internally consistent.
6. ~~Reframe `papers/timing/radical_timing.tex`~~ **DONE 2026-06-09** (structural pass; audit:
   `papers/timing/MANUSCRIPT_REFRAME_AUDIT_2026-06-09.md`).
6b. ~~Method-gain + systematics post-fix recompute~~ **DONE 2026-06-09**
   (`papers/scripts/{method_gain_postfix,systematics_postfix}/`): §5.3 = core 1.3 ps / tail 4.0 ps
   [3.4,4.7] split; tab_systematics regenerated (nominals cross-validate the authoritative table;
   veto-window rows added; DSB1 floor 18.8±0.8, fit-range +0.2); tab_methods regenerated.
   **Every active manuscript number is now post-fix/gated.**
6c. ~~Satellite-removal demonstration~~ **DONE 2026-06-09** (`papers/scripts/satellite_removal/`):
   f_tail 4.60→3.30% identical events @150 (Gaussian 1.24%); clipped-concentration + the dim-pulse
   reversal (srCFD 7.15% unclipped) → per-regime source rule demonstrated; §5.3 + fig installed.
   7 TODO-P2 remain, ALL writing/citation/placement.
6d. ~~Full prose + bibliography pass~~ **DONE 2026-06-09**: title decided (same-shower framing),
   abstract/intro/conclusions rewritten, satellite→appendix, depth kept (decisions recorded),
   engineering extrapolation + benchmark sentence added, bibliography 24 entries (4 marked
   TODO-P2-CITE). Manuscript builds; compliance greps clean. Remaining: money-plot restyle,
   benchmark table (cite-blocked), 4 citation verifications, author/funding blocks.
6e. ~~Circulation prep~~ **DONE 2026-06-09**: all 4 citations web-verified and installed
   (spacal/crilin/cantone2023/gundacker/lucchini fixed); benchmark table built + \input;
   money plot restyled (bands, MIXED-as-reference) → figs/thesis_postfix.{png,pdf};
   coauthor note + referee-risk memo written. **Manuscript is READY FOR INTERNAL COAUTHOR
   CIRCULATION.** ← NEXT: send to coauthors; in parallel (optional pre-submission): GATE-1
   logbook confirmation, Stage-1 G4 campaign, author/funding blocks, then Paper-3 gates
   (GATE 4 T_abs, GATE 7 energy backbone).
6f. ~~Production-formatting pass~~ **DONE 2026-06-09** (layout only, zero content change): 3 table
   overflows fixed; no-internal-super-title figure convention adopted (5 scripts re-run, values
   identical); legacy title strips cropped; overfulls 4→1 (cosmetic). Audit:
   `papers/timing/FORMAT_AUDIT_2026-06-09.md`. Manuscript circulates from this state.
7. Stage-1 GEANT4 campaign (serves both papers; ~days of CPU; on Argon) — now also pins the
   depth-dial ps→mm calibration (modal-dispersion v_eff, the ×1.28 slope excess).
8. Commit the whole campaign to git (lib fix + macros + figures + papers/memory_*).

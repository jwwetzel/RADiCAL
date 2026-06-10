# Campaign Snapshot — 2026-06-09

Frozen state of the publication campaign at the moment GATE 3 (depth dial) provisionally passed.
Purpose: reproducibility anchor for internal review; everything below was verified by direct command,
not recalled from memory.

## Repository state
- Branch: `main` · HEAD: `e7ac826` "Add methodDist.C — per-point (DW-UP)/2 distribution + Gaussian-fit diagnostic"
- **All campaign work is UNCOMMITTED** as of this snapshot.
- Modified (tracked): `analyze/studies/{methodCompare,methodDist,waveformProfiles}.C`,
  `lib/physics/{RadTiming,SelectionCuts}.h` (the estimator fix + continuous fiducial + veto constant),
  `lib/viz/PlotUtils.h` (DrawSuperTitle/DrawEnergyLegend/GridWithTitle), 6 regenerated narrative PNGs.
- Untracked (new): 9 analysis macros under `analyze/studies/` (adaptiveTiming, build_chain_html.py,
  mixedSeparate, monotonicityEvidence, monotonicityFix, outlierPeek, reductionQA, sigmaProbe,
  timingAllMethods, timingRegression), `chain_of_evidence.html`, 16 new narrative figures,
  `papers/` additions: EXPERT_PANEL_BLUEPRINT.md, two panel JSONs, six memory_*.md,
  `papers/scripts/depth_dial/` (AUDIT.md, depthDial.C, depth_dial_result.log),
  `papers/figures/depth_dial/` (depth_dial.png, depth_dial_diag.png).
- Suggested commit message (NOT yet committed):
  `Add RADiCAL paper memory infrastructure and depth-dial gate analysis`

## Environment
- ROOT 6.40.00 (Homebrew), macOS 15.7.7, arm64 (Apple Silicon).
- Env via `source setup.sh` (sets ROOT_INCLUDE_PATH to lib/{waveform,io,physics,viz}, RAD_DATA=repo root).
- RAD_YEAR unset → defaults to 2023 (DataPaths.h `radYear()`).

## Depth-dial run (GATE 3)
- Exact command:
  ```
  source setup.sh
  root -l -b -q -e '.L papers/scripts/depth_dial/depthDial.C+' -e 'depthDial()'
  ```
  (output teed to `papers/scripts/depth_dial/depth_dial_result.log`)
- Pre-registered method: `papers/scripts/depth_dial/AUDIT.md` (written BEFORE plotting).
- Inputs (reduced 'rad' trees, reduced on Argon Jun 7 2026):
  | file | bytes | mtime |
  |---|---|---|
  | data/2023/reduced/DSB1/25GeV.root | 501,953,492 | Jun 7 00:54 |
  | data/2023/reduced/DSB1/50GeV.root | 808,201,197 | Jun 7 00:55 |
  | data/2023/reduced/DSB1/75GeV.root | 596,740,047 | Jun 7 00:55 |
  | data/2023/reduced/DSB1/100GeV.root | 771,506,031 | Jun 7 00:51 |
  | data/2023/reduced/DSB1/125GeV.root | 723,530,216 | Jun 7 00:52 |
  | data/2023/reduced/DSB1/150GeV.root | 1,380,189,183 | Jun 7 00:53 |
  | + LUAG/TENERGY 50–150 GeV (cross-checks) | — | same reduction batch |
- Outputs: `papers/figures/depth_dial/depth_dial.png` (publication),
  `papers/figures/depth_dial/depth_dial_diag.png` (diagnostics),
  log `papers/scripts/depth_dial/depth_dial_result.log`.

## Scientific interpretation (one paragraph)
The mean dual-end timing asymmetry ⟨(t_DW − t_UP)/2⟩ of the DSB1 module, measured with the
amplitude-independent srCFD estimator on the full fiducial sample at fixed r = 2.5 mm and fixed
channel composition, drifts monotonically from −195.7 ps to −263.0 ps as the beam energy rises from
25 to 150 GeV — negative-going exactly as expected when shower maximum migrates downstream (closer to
the DW ends) by ≈1 X₀ per e-fold of energy. The fitted slope, −33.6 ± 2.9 ps per e-fold, agrees in
sign and order of magnitude with the geometric expectation −X₀/v_g ≈ −26 ps per e-fold; the ~28%
excess is consistent with modal dispersion lowering the effective light-propagation speed in the
quartz capillary, to be pinned by simulation. The dial is read by constant-fraction estimators only
(cfd05 mid-range slope ≈ −23 ps/e-fold concurs) while fixed-threshold LED crossings suppress it
(−5.5 ps/e-fold), establishing that the depth observable must be CFD-based. This is the first
measured, population-level demonstration that the RADiCAL dual-end timing carries longitudinal
shower-depth information — the physics hinge between the timing paper (the floor is depth) and the
energy/position paper (the fifth coordinate).

## Claim language in force at snapshot time
(authoritative copy: `papers/memory_claims_and_forbidden_language.md`)
- ALLOWED: "the measured mean of (t_DW−t_UP)/2 tracks the ln(E) shower-max migration, slope
  −33.6 ± 2.9 ps per e-fold, consistent with −X₀/v_g"; "toward 5D shower-maximum calorimetry"
  (population-level scoping mandatory; CFD-based caveat sentence required).
- FORBIDDEN (unchanged): per-event depth resolution / calibrated z; "5D" as achieved; t_diff as
  clock-referenced time-of-arrival; mm numbers as calibrated position resolution; any revision wording
  against the published 17.5 ps / 27 ps; E_SM-untagged energy-resolution numbers.

## UPDATE (later 2026-06-09): committed + GATE 2 closed
- Campaign committed: `45f15dc` "Add RADiCAL paper memory infrastructure and depth-dial gate analysis"
  (60 files). GATE 2 committed: `c8280f9` "GATE 2 passed: position arithmetic reconciled" (8 files).
  Working tree clean after both. Branch `main` is ahead of `origin/main` (not pushed).
- **GATE 2 outcome: PASSED (RECONCILED).** The 1.5 mm is a legitimate unbinned event-level residual
  RMS that survives a train/test split (Δ = 3.5 µm); the "3.6 mm comparator" was the t₀-inflated
  end-time-SUM upper bound and never applied to the position difference (memory mischaracterization,
  now corrected everywhere). Held-out residuals: 1.54/1.45 mm (x/y, ±6 mm window, 150 GeV),
  0.91/0.88 mm in the r<2.5 beam core, ~energy-independent. New finding: closure slope 0.698 —
  linear light-division estimator saturates beyond |x|≈3 mm. Joint-upper-bound language only; no
  intrinsic resolution; unfolding vs 3.6 demonstrated imaginary and refused.
  Products: `papers/scripts/position_reconciliation/`, `papers/figures/position_reconciliation/`.

## UPDATE 2 (later 2026-06-09): GATE 6 CONDITIONAL — the 0.99 kill-shot number RETIRED
- Pre-state verified: working tree clean at `c8280f9`; commits 45f15dc + c8280f9 confirmed.
- **GATE 6 outcome: CONDITIONAL.** The recorded MIXED ratio 0.99 (`mixedHeadToHead.C`) used
  brightness-threshold corner labels (amp > 0.65·max) that misassign SE-D (LuAG by GATE-1
  pulse-shape) into the DSB1 group at mid/high E — the 0.99 and χ²/ndf=0.4 are a label-mixing
  artifact and are retired from all claims. Replacement (common-event, GATE-1 labels, Down-layer
  same-material diagonal pairs, srCFD): per-capillary width parity within ~10–20%; 5-E mean ratio
  1.038, bootstrap 95% [1.033,1.043], honest scatter-based 1.04 ± 0.05; methods srCFD/cfd05/led =
  1.04/1.16/0.80; jackknife spread 0.003 (76 variants, no single-run dependence); swapped map
  inverts exactly; scrambled map collapses the amplitude control. Kill-shot claim rewritten
  (see claims memory): within-~20% same-shower parity vs ×2.3 cross-build → kinetics disfavored as
  dominant, with the position-coupling dilution caveat. `mixed_h2h.png` must be regenerated before
  any draft circulates. Products: `papers/scripts/mixed_killshot_bootstrap/`,
  `papers/figures/mixed_killshot_bootstrap/`.

## UPDATE 3 (later 2026-06-09): stale-claim purge + corrected kill-shot figure + post-fix (a,b) table
- Working tree was clean at `4c49888` before this pass.
- Stale 0.99/χ² language removed from all ACTIVE Paper-2 files (`papers/timing/THREAD.md`,
  `papers/memory_paper2_timing.md`); the panel blueprint + JSONs marked HISTORICAL with a
  supersession banner; `analyze/studies/mixedHeadToHead.C` carries a DEPRECATED header; the old
  figure renamed `figures/2023/narrative/mixed_h2h_DEPRECATED.png` + warning README.
- Official corrected figure: `papers/figures/mixed_killshot_bootstrap/mixed_h2h_corrected.{png,pdf}`
  (3 panels: widths / ratio with 1.04±0.05 scatter band / method dependence; approved caption in
  `makeMixedKillshotFigure_CAPTION.txt`).
- **POST-FIX authoritative four-build table:** `papers/tables/timing_fit_summary_2026-06-09.md` —
  DSB1 srCFD a=203±6, b=18.8±0.8, σ(150)=25.7±0.6 (HOLDS vs old 201/19.5; headline becomes 25.7);
  LuAG LED a=440±18, b=24.6±3.3 (a HOLDS; floor up from 19.8); TENERGY LED a=198±21, b=26.0±1.8
  (√(4/3) penalty); MIXED module-wide 253±29/34.3±2.4 (ill-posed ref). a-ratio 2.16 ("more than a
  factor of two"). **Floor narrative SOFTENED:** floors span 18.8–26.0 ps; "shared ~20 ps" retired
  without qualifiers; DSB1's 18.8±0.8 still CONFIRMS the published 17.5.

## Open flags at snapshot time
- 125→150 GeV local re-steepening of the dial (−56 ps/e-fold locally) — run-period check pending
  (addressed in DEPTH_DIAL_REVIEW.md diagnostics).
- GATE 1 (MIXED corner map): logbook confirmation pending (user action); data-only pulse-shape
  discriminant queued (`papers/scripts/mixed_corner_map/`).
- GATE 2 (position arithmetic), GATE 4 (T_abs), GATE 6 (bootstrap CI) open — see
  `papers/memory_analysis_gates.md`.

## UPDATE 4 (later 2026-06-09): manuscript reframe (structural pass)
- Working tree was clean at `56dda6d` before this pass. `radical_timing.tex` reframed around the
  gated story (audit: `papers/timing/MANUSCRIPT_REFRAME_AUDIT_2026-06-09.md`): new title
  ("Detected light yield governs…"), GATE-compliant abstract/conclusions, post-fix numbers
  (203±6/440±18, floors 19–26 qualified, 25.7±0.6, "more than a factor of two"), corrected MIXED
  mechanism, NEW same-shower-control section (§sec:mixed, fig mixed_h2h_corrected) and NEW depth-dial
  section (§sec:depth, acceptance-conditional), srCFD named at definition. Pre-fix method-gain and
  systematics numbers softened/quarantined behind TODO-P2 markers pending post-fix recompute.
  Builds clean under tectonic (PDF regenerated). 11 TODO-P2 markers tracked in memory_paper2_timing.

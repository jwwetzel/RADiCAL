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

## Open flags at snapshot time
- 125→150 GeV local re-steepening of the dial (−56 ps/e-fold locally) — run-period check pending
  (addressed in DEPTH_DIAL_REVIEW.md diagnostics).
- GATE 1 (MIXED corner map): logbook confirmation pending (user action); data-only pulse-shape
  discriminant queued (`papers/scripts/mixed_corner_map/`).
- GATE 2 (position arithmetic), GATE 4 (T_abs), GATE 6 (bootstrap CI) open — see
  `papers/memory_analysis_gates.md`.

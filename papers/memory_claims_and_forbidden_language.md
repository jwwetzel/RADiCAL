# Claims & Forbidden Language — the wording law for Papers 2 & 3

Updated: 2026-06-09. Non-negotiable. Source: ten-expert panel consensus
(`papers/radical_expert_panel_synthesis_2026-06-09.json` → `fourD_consensus`, `referee_defense`).

## ALLOWED claim language (exact formulations)

### Timing (Paper 2)
- "σ_t = 25.3 ps at 150 GeV (brightest-1000 selection), 27–29 ps over the full fiducial sample" —
  the MEASURED number always leads; full-fiducial always quoted beside any selection.
- "The data are **consistent with / confirm** the published floor of 17.5 ps" — the ONLY floor verbs.
- "Light yield, not WLS species, governs the timing" — the thesis sentence.
- "The stochastic term tracks detected light: a = 201 → 455 ps·√GeV from DSB1 to LuAG:Ce."
- "srCFD (saturation-recovered CFD): a constant-fraction discriminator applied at f of the LG-predicted
  true amplitude, recovering the edge of clipped pulses."
- "Sub-30 ps EM-shower timing from a 14×14 mm², 16-SiPM module" — compactness + channel economy framing;
  NEVER "fastest"/"best" (LHCb SPACAL W-GAGG is faster).
- "t_diff = (t_DW − t_UP)/2 is reference-free (MCP jitter, DRS4 inter-group phase, and common-mode edge
  cancel); t_sum is MCP-referenced" — both numbers always labeled by taxonomy.
- LuAG positioning: "the radiation-tolerant configuration, timing-viable at equal conditions" + the ×2
  light-recovery projection labeled "engineering extrapolation".

### Energy / position / 4D (Paper 3)
- "E_SM — shower-maximum sampled energy (≈3 LYSO tiles ≈ 1 R_M)" — the symbol defined once and attached
  to EVERY energy number including abstract and captions. "σ_E_SM/E ≈ 14% at 150 GeV."
- "Transverse residual of X mm against a wire-chamber comparator of σ ≈ 3.6 mm — a comparator-limited
  upper bound" (exact value pending GATE 2) — meets the parent's "few mm" goal AS STATED.
- 4D (allowed NOW, scope-locked): "a single ultra-compact module DEMONSTRATES a simultaneous
  four-dimensional (x, y, E_SM, t) measurement of EM showers at shower maximum" — demonstration, from
  one co-registered event sample, scope stated in the same sentence.
- Depth — GATE 3 PASSED structurally (2026-06-09), hardened by the stability review
  (`papers/scripts/depth_dial/DEPTH_DIAL_REVIEW.md`). THREE TIERS of allowed language:

  **Tier 1 — Conservative (always safe):**
  "The mean upstream–downstream timing asymmetry exhibits a monotonic energy dependence consistent
  in sign and order of magnitude with the expected logarithmic migration of shower maximum."

  **Tier 2 — Stronger (diagnostics support it; the acceptance qualifier is MANDATORY):**
  "The measured slope, −33.6 ± 2.9 ps per e-fold in beam energy for the full-fiducial event
  population, agrees at the order-of-magnitude level with the expected X₀/v_eff scale for
  shower-maximum migration; the slope is stable against energy-range, fiducial, per-capillary, and
  run-period variations." (If the brightness-acceptance dependence is mentioned: "the slope is an
  acceptance-weighted population quantity" — quartile range −18 to −95 ps/e-fold.)
  Required companion caveat: the depth readout is CFD-based (srCFD); fixed-threshold (LED) crossings
  do not track the emission centroid (dial suppressed to −5.5 ps/e-fold), and raw-fraction CFDs on
  clipped pulses are walk-dominated (cfd20/30: −86/−105).
  "toward 5D shower-maximum calorimetry" — population-level scoping mandatory.

  **Tier 3 — NOT YET ALLOWED (no analysis to date licenses these):**
  ✗ per-event longitudinal position resolution; ✗ a calibrated z coordinate (v_eff unknown to ~30%);
  ✗ achieved 5D calorimetry; ✗ clock-referenced time-of-arrival from this observable;
  ✗ "the depth slope alone proves shower-depth reconstruction" (it is population-level consistency;
  per-event corroboration + G4 truth regression are the unlocks).
- After FULL unlock chain (gates 3+4+G4 regression): "single-module 5D estimator (x, y, z_depth, E_SM, t)
  with σ_z ≈ 1 X₀ per event".

## FORBIDDEN (must not appear; each with its unlock if one exists)
1. ✗ "4D calorimeter meeting program goals" — demo numbers ≠ the 10 ps / 10%/√E ⊕ 0.7% / 1 mm targets;
   never put 52%/√E and 10%/√E on the same axis. (No unlock with this dataset.)
2. ✗ "5D calorimeter" / any per-event depth MEASUREMENT or resolution today — no independent depth truth
   exists; spread-matches-1-X₀ is consistency, not validation. UNLOCK: gates 3→4→G4 regression chain;
   even then only "toward 5D" / "5D estimator", never "5D imaging calorimeter".
3. ✗ 0.9 or 1.5 mm as a calibrated position RESOLUTION — comparator is ~3.6 mm, beam spot 2.9 mm.
   Only "comparator-limited residual / upper bound". UNLOCK: GATE 2 + intrinsic decomposition.
4. ✗ 25–27 ps as absolute time-of-arrival or ToF — that is the self-referenced differential; the
   clock-referenced number is ~50 ps class. UNLOCK: GATE 4 (T_abs), and even then quote both with taxonomy.
   Map every physics case (H→γγ vertexing, BIB windows) to the RIGHT observable.
5. ✗ "Material-independent ~20 ps floor" as a measurement — LuAG's best measured point is 42 ps;
   b=19.8±5.6 is a form-dependent extrapolation. Allowed: "the fitted floors are consistent" + caveat.
6. ✗ ANY wording that reads as a revision of the published 17.5 ps / 27 ps — verbs locked to
   "confirms / consistent with". (Permanent.)
7. ✗ HGCAL-style "imaging" language — one longitudinal sample at shower max is not imaging. No
   array-level, hadron, pileup, or photon-conversion performance claims (electrons, normal incidence,
   single module, 25–150 GeV).
8. ✗ "Radiation hardness demonstrated here" — zero irradiation data in this dataset; the best-timing
   build is organic DSB1. Allowed: component-level tolerance BY CITATION (parent refs 14–16) + measured
   timing-viability of the LuAG configuration + honest SiPM-fluence paragraph.
9. ✗ Quoting best-selection numbers first — full-fiducial/all-event numbers always lead; brightest-1000
   or best-bin always second with the selection named.
10. ✗ "crystal differs between builds" — INVARIANT WORDING RULE: LYSO:Ce + W are common to all builds;
    the WLS capillary is the variable. (House rule, predates panel.)

## Referee-defense pairs (abbreviated; full list in synthesis JSON)
- "1.5 mm vs 3.6 mm impossible" → reconciliation shipped in-paper (GATE 2 outputs).
- "(DW−UP)/2 cancels arrival time" → two-end algebra section + both numbers with taxonomy.
- "Floor is an extrapolation" → measured 25.3 headlined; dual-form fits with χ²; LuAG 42 ps stated plainly.
- "Method shopping" → ex-ante algorithmic regime rule + full 8-source survey published.
- "MIXED map circular (brightness-inferred)" → logbook + pulse-shape discriminant (GATE 1).
- "No competitive context" → merged benchmark table with object-type column (per-MIP vs per-shower).
- "52%/√E isn't competitive" → E_SM tag on every number; parent §5.1.3 scoping repeated.
- "Saturated-pulse timing on a SCA digitizer" → digitizer-forensics figure (stop-cell, cell-width, satellite cure).
- "Garden of forking paths" → locked-pipeline OOS table; deterministic brightest-fraction selection.
- "4D/5D inflation" → locked wording above; HGCAL/SCEPCAL cited wherever adjacency is implied.

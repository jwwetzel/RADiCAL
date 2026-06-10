# Paper 2 — Timing (working memory)

Updated: 2026-06-09. Manuscript: `papers/timing/radical_timing.tex` (NEEDS REFRAME — see below).
Full blueprint: `papers/EXPERT_PANEL_BLUEPRINT.md` §Paper 2; machine-readable:
`papers/radical_expert_panel_synthesis_2026-06-09.json` → `paper2_blueprint`.

## Title (panel-proposed)
"Comparison of DSB1 and LuAG:Ce wavelength shifters for precision EM-shower timing at shower maximum
in an ultra-compact W/LYSO:Ce sampling calorimeter (RADiCAL), 25–150 GeV"

## Thesis (one sentence, locked)
**Light yield — not WLS species — governs shower-max timing**: the stochastic term tracks detected light
(DSB1 a≈201 → LuAG a≈455 ps·√GeV, ×2.3) while the ~20 ps floor is shared among LYSO builds; the MIXED
in-module comparison (same showers, paired per-capillary widths agreeing at the 10–20% level —
σ_DSB1/σ_LuAG = 1.04 ± 0.05 srCFD, GATE-6 corrected) DISFAVORS the LuAG-kinetics confound as dominant
with zero cross-build systematics; srCFD (né lgcfd) is the ENABLING method, and the floor CONFIRMS the
published 17.5 ps as shower-depth geometry.

## Known results — POST-FIX AUTHORITATIVE TABLE (2026-06-09)
Single source of truth: `papers/tables/timing_fit_summary_2026-06-09.md` (+ figure
`papers/figures/timing_fit_summary/timing_fit_summary.png`). Supersedes ALL pre-fix (a,b).
| result | value | support |
|---|---|---|
| DSB1 srCFD brightest-1000 | **a=203±6, b=18.8±0.8**, MONO; **σ(150)=25.7±0.6 ps**; L150=5185 mV | `timingFitSummary.C` |
| LuAG LED | **a=440±18, b=24.6±3.3**, MONO; best measured 44.4 ps @150; floor = extrapolation | same |
| MIXED srCFD module-wide | a=253±29, b=34.3±2.4 (ILL-POSED ref row; paper uses per-material) | same |
| TENERGY LED | **a=198±21, b=26.0±1.8**, one NM point; √(4/3) 3-cap penalty → ~22.5 4-cap-equiv | same |
| cross-build a-ratio | **2.16** (was "~2.3"; narrative updated to "more than a factor of 2") | same |
| floors | span 18.8–26.0 — "shared ~20" SOFTENED, see table addendum | table addendum |
| Kill-shot | **0.99 RETIRED (GATE 6)**. Same-shower width parity 10–20%; srCFD paired 1.04±0.05 (scatter); methods 0.80–1.16 | `mixed_h2h_corrected.{png,pdf}` |
| Published baseline | 256/√E ⊕ 17.5 ps; 27 ps @150 (BestMinus, cfd05) — confirm, never revise | parent Eq. 2 |
| Estimator fix | adopted methods monotonic; cfd05 Δ≤0.4 ps | `monotonicity_fix.png` |
| MIXED clip fractions | DSB1 corners 33–55%, LuAG corners 0–24% (E-dependent) | `mixedSeparate.C` |
| MCP1−MCP2 | 71–74 ps flat vs E (inter-group domino phase; ~15 ps stop-cell-corrected — VERIFY) | memory dataset_2023 |

## Section outline (panel S1–S12, condensed)
S1 promise-kept intro + competitive context · S2 module/builds + invariant wording + channel map ·
S3 nomenclature Table 1 (t_diff/t_sum; LED/CFD-f/srCFD; what cancels) + σ-extraction box ·
S4 digitizer timebase forensics (clip-fraction table, stop-cell, MCP jitter attribution) ·
S5 srCFD as enabling method (hglg calibration, mechanism, Δσ(E) panel, satellite cure, closure test) ·
S6 8-source survey + ex-ante adoption rule + veto characterization ·
S7 money plot (4-build σ_t(E)) + universal collapse (σ_t vs detected light) + photostatistics closure ·
S8 kill-shot (per-material diagonals + in-event pairwise, claim ladder) ·
S9 floor physics (dual-form fits, drift-vs-lnE consistency panel, G4 band; "confirms 17.5") ·
S10 benchmark+requirements table + deployability box · S11 systematics/OOS validation · S12 conclusions.

## Money figures
1. `light_yield_thesis.png` UPGRADED: 4-build σ_t(E), **error bars restored (FATAL fix)**, (a,b) inset
   table, collider requirement bands. 2. Universal collapse: σ_t vs sum_lg, all builds/energies, one law,
   floor asymptote, same-shower parity inset (graphical abstract). 3. `mixed_h2h_corrected.png`
   (GATE-6: widths + ratio panel with the 1.04±0.05 scatter band — never the retired 0.99).
4. method mechanism (pulse slopes 102 vs 381 mV/ns) + Δσ(E). 5. satellite before/after. 6. floor figure
(dual-form fits + plateau + G4 band + **mean-t_diff-vs-lnE drift panel ← GATE 3 feeds this**).
7. MCP-jitter forensics. Tables: nomenclature; 8×4 survey (a,b,χ²); benchmark; OOS; ps-budget.

## Unresolved gates blocking this paper
- GATE 1 — data half CONFIRMED (8/8 pulse-shape); logbook half = user action (belt & suspenders).
- GATE 6 — **CONDITIONAL (2026-06-09): the 0.99 is RETIRED** (brightness-label artifact, SE-D
  misassigned). S8 must be REWRITTEN around the replacement result: same-shower per-capillary width
  parity within ~10–20% (srCFD paired 1.04±0.05 scatter-based; methods 0.80–1.16) vs the ×2.3
  cross-build a-ratio → kinetics DISFAVORED as dominant, with the position-coupling dilution caveat.
  Corrected figure SHIPPED 2026-06-09: `mixed_h2h_corrected.{png,pdf}` (old figure DEPRECATED with
  warning README; generator `mixedHeadToHead.C` carries a DEPRECATED header).
  Remaining: bootstrap CI on the srCFD-vs-cfd05 method-gain claims (separate, smaller task).
- GATE 3 depth dial — **PASSED 2026-06-09**: the S9 drift-consistency panel is now MEASURED
  (−33.6±2.9 ps/e-fold vs −X₀/v_g = −26.3; `papers/figures/depth_dial/depth_dial.png`). S9 can write
  "the floor IS depth physics — shown, not asserted." Bonus S5/S6 material: the LED family does not
  read the dial (−5.5 ps/e-fold) — independent evidence that srCFD samples the light centroid while
  fixed-threshold crossings sample first photons.
- B-list: post-fix recompute of all (a,b); srCFD closure test; satellite demonstration; MCP attribution check.

## Required simulation (Stage 1, shared with Paper 3)
Energy-deposition-only G4 (exact stack, 25–150 GeV e−, 2.9 mm spot, 1% dp/p): per-event light-weighted
depth → predicted floor band + σ_t(Npe) twin of the money plot. Stage 2 (optional): single-capillary
optical run → v_g ≈ 200–205 mm/ns validation. Reproducibility paragraph mandatory (G4 version,
FTFP_BERT+EMZ, cuts, Birks, counts).

## Must-cite (top of list)
Parent (Eq 2, §7.2, Fig 19, Fig 7) · Anderson et al. NIM A 794 (2015) 7 (prior-art LYSO timing —
mandatory, author-network exposure) · Gundacker/Cates CTR canon · Hu NIM A 954 (2020) + Lucchini JINST 8
(LuAG) · Gedcke–McDonald (CFD origin), Genat–Varner–Tang–Frisch (digital timing), Cartiglia (time-walk) ·
Ritt DRS4 + Stricker-Shaver IEEE TNS 61 (2014) · CMS MTD/HGCAL TDRs, ATLAS HGTD, LHCb FTDR + SPACAL ·
FCC CDRs, DRD6, RDC9 · PDG stats, Efron/Tibshirani · Ruchti Annu.Rev. 46 (1996).

## Manuscript state + next actions
- **REFRAMED 2026-06-09** (audit: `papers/timing/MANUSCRIPT_REFRAME_AUDIT_2026-06-09.md`).
  `radical_timing.tex` now: new title ("Detected light yield governs…"), GATE-compliant abstract,
  post-fix numbers everywhere (Table tab:builds = the authoritative table), corrected MIXED
  mechanism (ill-posed module-wide estimator, NOT "compromised LG readout"), **NEW §sec:mixed**
  (same-shower control: corrected map, per-E ratios, 1.04±0.05, method dependence, dilution caveat,
  fig mixed_h2h_corrected) and **NEW §sec:depth** (depth dial, acceptance-conditional, consistency
  reading only), rewritten Discussion/Conclusions (disfavors form). srCFD named at definition.
  Builds clean under tectonic. Figures copied into figs/: thesis_postfix, mixed_h2h_corrected,
  depth_dial.
- **METHOD GAIN + SYSTEMATICS RECOMPUTED POST-FIX (2026-06-09):**
  §5.3 now quotes the identical-event result (`papers/scripts/method_gain_postfix/`, table
  `papers/tables/method_gain_postfix_2026-06-09.md`, figure figs/method_postfix.png):
  **core-width gain 1.3 ps @150 (26.9→25.6); tail-sensitive (truncated-RMS) gain 4.0 ps,
  68% CI [3.4,4.7] paired bootstrap** (the split IS the satellite-removal story); Δ(E):
  +5.2 @25 (noise foot), ≈0 @50–75, +2.7/+3.6 @100/125; gain concentrated in saturated events
  (fully clipped @150: 27.4→25.3); floors 21.8±2.7 (cfd05, poor non-monotonic fit) → 19.9±0.8
  (srCFD). LED−srCFD @150: +4.2 [3.6,4.9] (tail-RMS). The OLD 3.1 ps / 22.1→19.5 are retired
  and NOT reproduced — the post-fix split (1.3 core / 4.0 tail) replaces them.
  `tab_systematics.tex` regenerated post-fix (`papers/scripts/systematics_postfix/`): nominals
  match the authoritative table exactly (25.7/30.3/39.6/44.4 — machinery cross-validated);
  totals ±1.0/±1.1/±0.9/±1.9; NEW veto-window rows (largest: TENERGY −2.3 @1.5 ns); DSB1 floor
  block b=18.8±0.8, fit-range shift +0.2. `tab_methods.tex` regenerated post-fix.
- **Every active number in the manuscript is now post-fix/gated.**
- **SATELLITE-REMOVAL DEMONSTRATION DONE (2026-06-09)** (`papers/scripts/satellite_removal/`,
  figure `papers/figures/satellite_removal/satellite_removal.{png,pdf}` + approved caption):
  identical events @150: f_tail (common 66 ps window) 4.60%→3.30% (46→33 of 1000; Gaussian 1.24%),
  robust IQR nearly identical (28.7 vs 28.3) → the cfd05 excess lives in the tails; production
  widths reproduce §5.3 exactly (26.9/25.6). 75 GeV full-fiducial split: cfd05 tail excess
  concentrates in clipped events (3.39 vs 1.95%); REVERSAL on dim unclipped pulses (srCFD 7.15% —
  noisy LG anchor) → each estimator biased outside its design regime = the adopted-source rule's
  basis. §5.3 paragraph + Fig fig:satellite installed; §4 TODO replaced with a pointer.
- **Open TODO-P2 markers (7, ALL writing/citation/placement — no analysis left):**
  title confirmation; bibliography build-out; money-plot final style (bands+error bars);
  §mixed citations (+ optional dilution quantification via G4, can remain a stated caveat);
  depth-dial main-vs-appendix; satellite figure main-vs-appendix; engineering-projection +
  benchmark table. Author block + funding placeholders.
- Next: the full prose + bibliography pass (no analysis prerequisites remain).

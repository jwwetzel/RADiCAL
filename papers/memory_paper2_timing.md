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
in-module comparison (same showers, per-capillary ratio ≈0.99) resolves the LuAG-kinetics confound with
zero cross-build systematics; srCFD (né lgcfd) is the ENABLING method, and the floor CONFIRMS the
published 17.5 ps as shower-depth geometry.

## Known results (with provenance)
| result | value | support |
|---|---|---|
| DSB1 lgcfd brightest-1000 | a≈201, b=19.5±1.1 (pre-fix) / 18.8±0.8 (post-fix) | `methodCompare.C`, `timingRegression.C` |
| LuAG led | a≈455, b=19.8±5.6 (best measured point 42 ps) | same |
| MIXED lgcfd (post-fix, monotonic) | b≈34.3 (module-wide; ill-posed — use per-material) | `timingRegression.C` |
| TENERGY led | 240/22±3 (RECORDED — recompute post-fix before quoting) | report/memory |
| Kill-shot in-event ratio | ≈0.99 (needs bootstrap CI — GATE 6) | `mixed_h2h.png` |
| Published baseline | 256/√E ⊕ 17.5 ps; 27 ps @150 (BestMinus, cfd05) | parent Eq. 2 |
| Estimator fix | adopted methods monotonic; cfd05 Δ≤0.4 ps | `monotonicity_fix.png` |
| MIXED clip fractions | DSB1 corners 33–55%, LuAG corners 0–24% (E-dependent) | session log, `mixedSeparate.C` |
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
   floor asymptote, 0.99 inset (graphical abstract). 3. `mixed_h2h.png` + ratio panel with bootstrap CI.
4. method mechanism (pulse slopes 102 vs 381 mV/ns) + Δσ(E). 5. satellite before/after. 6. floor figure
(dual-form fits + plateau + G4 band + **mean-t_diff-vs-lnE drift panel ← GATE 3 feeds this**).
7. MCP-jitter forensics. Tables: nomenclature; 8×4 survey (a,b,χ²); benchmark; OOS; ps-budget.

## Unresolved gates blocking this paper
- GATE 1 MIXED corner map logbook (blocks S8 kill-shot).
- GATE 6 bootstrap CI on 0.99 + method-gain claims.
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
- `radical_timing.tex` currently centers DSB1-method-improvement; MUST be reframed to the materials
  comparison (panel worklist #14). Existing sections to keep: systematics table, optimization appendix.
- Next: close GATES 1/6/3 → post-fix (a,b) table → reframe tex → satellite + closure analyses → S10 tables.

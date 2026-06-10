# Manuscript Reframe Audit — radical_timing.tex (2026-06-09)

Baseline: commit 56dda6d, working tree clean. Manuscript: 461 lines, elsarticle 5p twocolumn,
8 figures (figs/), 2 table includes (tab_methods unused?, tab_systematics input at L363).

## Current manuscript state

Structure (pre-reframe): 1 Introduction · 2 Experimental setup (module/builds; readout/beam) ·
3 HG/LG saturated-edge recovery · 4 Estimator + selection · 5 Results (5.1 distributions,
5.2 "Light yield sets the resolution" + Table tab:builds, 5.3 method gain, 5.4 systematics) ·
6 Discussion · 7 Conclusions · App A selection optimization. Bibliography: only 4 refs.

### Stale claims found (all PRE-FIX numbers / forbidden phrasing)
| location | stale content | fix |
|---|---|---|
| Title | "Light yield, **not** WLS species, **sets**…" | soften verb → "governs" (claims-law thesis verb) |
| Abstract | "floor b≈20 ps consistent across LYSO builds" unqualified; "25.4 ps"; "therefore governed by light yield, not by the WLS species: at equal light all materials approach the same floor" | qualified floor span 19–26; 25.7±0.6; disfavors-form scoping |
| §5.2 | a 201→455 "factor 2.3"; floors 20±1/22±3/20±6; MIXED "low-gain readout partly compromised" (WRONG mechanism) | 203±6→440±18 "more than a factor of two"; 18.8±0.8/26.0±1.8/24.6±3.3; MIXED = ill-posed module-wide estimator (mixed clip regimes) |
| Table tab:builds | all pre-fix (201/240/233/455; σ150 25.4/28.7/40.0/41.2) | post-fix authoritative (203/198/253/440; 25.7/30.3/39.6/44.4) + caveat notes |
| Fig thesis.png | pre-fix money plot | replaced by post-fix `timing_fit_summary.png` (copied to figs/thesis_postfix.png; TODO restyle) |
| §5.3 | 22.1→19.5 ps floors, a 205→201, "best 25.4" | pre-fix method-gain numbers: softened to qualitative + TODO-P2 recompute (the method-gain bootstrap is a known GATE-6 follow-up) |
| §5.4 + tab_systematics.tex | nominal 25.3/19.5/34.6/19.8 | TODO-P2: regenerate paperSystematics post-fix; floor sentence softened |
| Discussion | "times worse … by **exactly** the factor its light implies, not more" | overclaim → GATE-6 within-10–20% form |
| Conclusions | unqualified shared floor; 25.4 | qualified span; 25.7±0.6 |

### Missing vs the authoritative story
- **The MIXED same-shower control section does not exist in the tex** (the old THREAD bead pointed
  to the retired figure; nothing was integrated). → NEW Section 5 (see outline).
- No depth-dial section. → NEW short section with acceptance-conditional scoping; TODO placement.
- `mixed_h2h_corrected` and the post-fix four-build figure not referenced anywhere.
- tab_methods.tex (pre-fix 25.3/19.5) is NOT \input anywhere — left in repo; TODO regenerate or drop.
- Bibliography is 4 refs — far below the panel's must-cite list. TODO-P2 markers added.

### Good news
- The retired "0.99 / χ²/ndf=0.4 / statistically indistinguishable" never entered the tex.
- §3 (srCFD mechanism), §4 (estimator), App A survive structurally intact.

## New thesis (adopted verbatim, per instruction)
The timing performance of the RADiCAL shower-maximum module is governed primarily by detected
light yield and module response rather than by WLS material identity alone. Across the four-build
energy scan, the stochastic term differs by more than a factor of two between the DSB1 and LuAG
builds, while the same-shower MIXED comparison shows DSB1/LuAG per-capillary timing widths agreeing
at the 10–20% level under the primary srCFD estimator. This combination disfavors intrinsic WLS
re-emission kinetics as the dominant explanation for the cross-build timing difference, within this
geometry and at these light levels.

## Revised outline → disposition of existing sections
1. **Introduction** — KEPT, light edits (already frames the species-vs-light question; add the
   disentanglement sentence + TODO citations).
2. **Apparatus and datasets** — KEPT; MIXED re-described as the same-shower control with the
   pulse-shape-confirmed corner map (one added paragraph).
3. **Timing reconstruction** — KEPT (current §3+§4 merged conceptually): srCFD primary; LED/cfd05
   diagnostics; estimator taxonomy naming added (TODO full nomenclature table).
4. **Four-build timing response** — REWRITTEN (current §5.2): post-fix numbers, qualified floor
   language, corrected MIXED mechanism, post-fix table + figure.
   (5.1 distributions KEPT; 5.3 method gain KEPT with softened numbers + TODO recompute;
   5.4 systematics KEPT + TODO regenerate table post-fix.)
5. **MIXED same-shower control** — NEW (the GATE-6 result): corrected map, per-energy ratios,
   1.04±0.05 scatter, method dependence, position-coupling dilution caveat, mixed_h2h_corrected.
6. **Longitudinal timing asymmetry (depth dial)** — NEW, short, acceptance-conditional scoping;
   TODO-P2: decide main text vs appendix.
7. **Discussion** — REWRITTEN claims (disfavors form; light-budget engineering framing).
8. **Conclusions** — REWRITTEN, claims-law compliant.
App A — KEPT. tab_systematics — KEPT with TODO. Title — verb softened.

## Claim constraints enforced (numbers in prose/table now)
DSB1 203±6 / 18.8±0.8 / 25.7±0.6 · LuAG 440±18 / 24.6±3.3 / 44.4 · TENERGY 198±21 / 26.0±1.8
(√(4/3) caveat) / 30.3 · MIXED module-wide ref 253±29 / 34.3±2.4 / 39.6 · MIXED same-shower ratio
1.04±0.05 (scatter), per-E 0.940/1.058/0.929/1.072/1.194, methods 1.038/1.161/0.795 · depth dial
−33.6±2.9 ps/e-fold (acceptance-conditional). Forbidden list as in claims-law memory (0.99, χ²=0.4,
indistinguishable, unity, proves, unqualified shared-20, 25.3/201/455/×2.3).

## Risks before a full prose pass
- §5.3 method-gain and §5.4 systematics tables still carry pre-fix derivations (TODO-marked, softened
  in prose; must be recomputed before submission).
- Bibliography skeleton (4 refs) vs the must-cite list.
- thesis_postfix.png is the summary-figure style, not the final money-plot style (error bars +
  requirement bands per panel blueprint).
- Author block placeholder; funding text placeholder.

# RADiCAL Paper 1 — Definitive Outline (expert-panel synthesis)

Panel: revered instrumentation author (story/beauty) · brutal NIM A referee (rigor) ·
FCC calorimeter expert (physics). Synthesized by the lead author. 2026-06.

## (a) RECOMMENDED TITLE
**"Light yield, not material, governs the timing of RADiCAL sampling-calorimeter
modules: edge recovery and a four-build comparison."**
- Soften "sets/not species" → "governs / dominant lever." The verb must be defensible:
  LuAG's coincidence-time-resolution ratio (≈2.6× vs LYSO at fixed 511 keV) is suspiciously
  close to the light ratio (2.3×), so a *between-build* comparison alone cannot separate
  "less light" from "intrinsically slower scintillator." The softened claim survives
  **because** it is defended by the in-module kill-shot (§V.C), not by the four-curve ordering.
- Fallback: "Light yield is the dominant design lever for the timing of RADiCAL modules."

## (b) THROUGH-LINE (the story)
A calorimeter that times a shower to tens of ps lives or dies by the steepness of one rising
edge — and here that edge is the part the readout throws away, clipped flat at 820 mV. We
**recover the discarded edge** event-by-event from the low-gain copy meant only for energy;
then, timing every build on that same recovered edge, we show four WLS materials separate
**almost entirely by how much light they collect, not which crystal they are** — proven
cleanest **inside one MIXED module where LYSO and LuAG capillaries see the same showers and
time identically.** The residual ~20–25 ps floor is shared and consistent with the shower's
own depth fluctuation (~1 X₀) — where a depth-corrected Paper 2 begins.
- Nested architecture: **edge recovery is the engine; "light, not material" is the destination.**
- Two bookend images: the **recovered edge (3.7× steeper)** and the **shared floor (~20 ps ≈ 1 X₀)**.

## (c) WHAT'S NEW vs NIM A 1068 (honest box, in the intro)
**New:** (1) HG/LG edge recovery *for timing* (published used LG for energy only) — clean
same-events head-to-head: 3.1 ps @150, floor 22→20; (2) four-build material comparison (the
published §7.2 mandate); (3) **in-module confound-free proof** (MIXED LYSO vs LuAG capillaries,
same showers); (4) the floor given a physical identity (shower-depth ≈ 1 X₀) + the route below it.
**NOT new (say so, to disarm the referee):** 25.4 ps is *not* a new best — Ref. [2]'s own
brightest amplitude bins already reach 24–27 ps @150. We **confirm**, not revise, 17.5 ps & a/√E.

---

## I. INTRODUCTION
- **Purpose:** pose species-vs-light with stakes; hook on the discarded edge; state two deliverables.
- **A.** Hook: the edge you must time is the part the readout discards (clip → CFD on the shallow foot).
- **B.** FCC stakes, QUANTITATIVE (currently the weakest part): FCC-hh pile-up ≈1000/crossing,
  collision-time spread ~150–180 ps RMS; two specs — 4D vertex assoc ≈30 ps, aggressive PU 5–10 ps;
  say which RADiCAL meets.
- **C.** The design payoff (keep verbatim): choose material for rad-hardness/mechanics "without a
  separate timing penalty."
- **D.** Baseline in one sentence (Ref. [2]: DSB1, cfd05, 27 ps, a=256/b=17.5; its §7.2 mandated this).
- **Reviewer-proofing:** say up front this *confirms & extends* [2]; raise the LuAG-kinetics confound
  yourself and promise to bound it.
- **Figure:** none (or the cartoon teased).

## II. EXPERIMENTAL SETUP
- **Purpose:** let the reader hold the device; make (DW−UP)/2 and the depth-floor visually inevitable.
- **A.** Module + builds + a REAL light-yield number per build (p.e. or ΣLG, not "factor ~2"). DSB1
  25–150; others 50–150 (no 25 GeV).
- **B.** Readout/beam/reference; DRS4 group map honestly (which corners share a group; SW-U separate,
  ~70 ps inter-chip, and how the differential handles it).
- **Figure:** **[BUILD] apparatus 3D cutaway** — highest-value missing figure; W plates, 4 corner
  capillaries through full depth, UP/DW ends, 8 SiPMs, MCP, wire chamber, depth ≈25 X₀.
- **Reviewer-proofing:** fiducial stated once here; disclose min-channels ≥2/4.

## III. RECOVERING THE SATURATED EDGE (HG/LG METHOD) — the engine, subordinated
- **Purpose:** deliver the method, explicitly in service of the comparison.
- **A.** Open: "To compare materials fairly we must time them all on the same edge; saturation denies
  it, so we recover it."
- **B.** Clip problem. **C.** Recovery: HG_true = a_c + b_c·LG_peak [Eq.1]; CFD@15% on the steep edge.
- **D.** Slope payoff 3.7× — but state HONESTLY: by the slew budget (~2 ps), 3.7× alone can't make the
  3.1 ps gain; the gain is from escaping the corrupted clipped-foot crossing + low-pulse noise.
- **E.** Per-build source (lgcfd vs led); the fragility (cfd05 fails on dim builds) is itself a result.
- **Figures:** **[BUILD] HG/LG cartoon** (show the trick) → **[HAVE,shrink title] clip.png** (villain)
  → **[REBUILD] hglg.png** (the ONE unpublishable figure: garbled axes/debug text → collapse to one
  clean HG∝LG hockey-stick channel + Eq.1) → **[HAVE,enhance] pulse.png** (the jewel; add a "3.7×
  steeper here" annotated gap).
- **Reviewer-proofing:** disclose the Eq.1 fit (per-energy median, slope 1.5–8, n>300) + floor
  sensitivity to ±5% slope; state cfd05=5% where first used.
- **FCC note:** front-end saturation is generic at high occupancy; edge recovery is transferable.

## IV. ESTIMATOR & SELECTION
- **A.** t=½(t̄_DW−t̄_UP), mean over ≥2/4 each side; cancels MCP; 4-ch averaging cuts uncorrelated noise.
- **B.** Brightest-1000 by ΣLG: "fixing the count fixes the precision of every point... the comparison
  is between materials, not sample sizes."
- **C.** Quote the TYPICAL-shower number (~30 ps quantile) alongside the best-case 25.4 ps.
- **D.** N=1000, σ_stat≈σ/√(2N)≈0.7 ps.
- **Reviewer-proofing:** surface the OOS/selection-stability number (fold-by-run 26.2/25.6/26.4;
  equal-width argmin overfits → not used); note (DW+UP)/2 is reference-limited (~70 ps), cross-check only;
  note ΣLG is MIXED's compromised quantity.

## V. RESULTS — reordered for ascending power
### V.A Timing distributions
- Clean Gaussians; differential narrower than absolute; width ↓ with brightness.
- **Figure [HAVE, FIX] dist.png:** repair collapsed title; enlarge per-panel σ; color rows to match.
- **Reviewer-proofing (MANDATORY): the satellite** — add inset showing the ~0.2 ns satellite, its
  fractional area (~15% per [2]), and Δσ(ps) from in/excluding it; cite [2] Figs.18/19.
### V.B Light yield sets the stochastic term (money plot)
- a: 201(DSB1)→240(TEN)→233(MIX)→455(LuAG), ×2.3. "Line the builds up and they sort by brightness."
- **Figure [REBUILD] thesis as two-panel:** (a) σ vs E + a/√E⊕b fits **WITH ERROR BARS** (vertical
  σ/√2N AND horizontal energy-spread — the current fig has NONE; single most important referee fix);
  (b) σ vs *light yield* (cleaned sigmat_vs_lightyield) captioned HONESTLY as cross-config-confounded,
  pointing to §V.C as the clean version.
- **Reviewer-proofing (FATAL): the floor, honest.** Don't headline "material-independent 20 ps." Per
  build give BOTH forms with χ²/ndf/p: photostat 1/√E vs slew 1/E (DSB1 19.5±1.1 vs 26.1±1.2; LuAG
  19.8±5.6 vs 39.3±2.7). State plainly LuAG's brightest 150 GeV point is still **42 ps — it never
  reaches 20 ps in data**; "LuAG reaches the same floor" is an extrapolation. The discriminating power
  is almost all in the single DSB1 25 GeV point. Headline the *measured* 25.4 ps; floor carries ±3–6 ps
  form systematic. Report ρ(a,b)≈−0.88, PDG √(χ²/ndf) scaling. The surviving claim: "the floor is
  consistent within the large extrapolation uncertainty AND is independently argued to be shower-depth
  physics (§VI) — making material-independence a *prediction*, not a coincidence."
### V.C In-module head-to-head (NEW — THE KILL SHOT, currently absent)
- **Purpose:** confound-free proof that ordering is light, not crystal — every systematic cancelled.
- Inside MIXED: 2 corners DSB1(LYSO), 2 LuAG, same showers, same MCP, same DRS group → per-capillary
  σ_t(E) curves lie on top of each other 50–150 GeV. Only difference is the crystal — and it doesn't matter.
- **Figure [HAVE→PROMOTE] mixed_h2h.png.** Caption: "With every systematic held fixed by construction
  — same shower, same reference, same digitiser — a DSB1 and a LuAG capillary time the same."
- **Reviewer-proofing:** uses TIMING ONLY (no LG) → MIXED's compromised LG does NOT contaminate it (this
  is what makes the kill-shot clean). **Confront the LuAG-kinetics confound here:** the CTR ratio
  (≈2.6×) ≈ the a-ratio (2.3×), so between-build can't separate light from rise-time — the in-module
  same-shower identity IS the separation. This is the one place the thesis is defended against its
  strongest attack.
### V.D Method-gain coda (validation)
- cfd05 vs hg_lgcfd, same brightest-1000 DSB1: 3.1 ps @150, floor 22.1→19.5, energy-dependent, a~unchanged.
- **Figure [HAVE] method.png** — tighten caption (same events, only time differs); add paired error.
- **Reviewer-proofing:** quote 3.1±X ps with PAIRED (correlated) error; reconcile the four DSB1 floor
  numbers (17.5/19.5/20/22.1) into ONE consistent treatment + a table of which selection/source/form each is.

## VI. DISCUSSION
- **A.** Light = design lever; LuAG worse "by exactly the factor its light implies" — cite an ACTUAL
  LuAG irradiation measurement, not only the generic ceramic ref.
- **B.** Floor = shower-depth, with NUMBERS: v_fiber≈200 mm/ns; 20 ps→σ_z≈4 mm; 28 ps→5.6 mm; 1 X₀=5.4 mm.
  Reconcile fitted-20 vs depth-28 (slew-form floor ~26→σ_z 5.2 mm ≈1 X₀, strengthening the interpretation).
  Vivid: "a shower a centimetre deeper sends light to DW sooner, UP later; (DW−UP)/2 measures exactly that."
- **C.** Instrumental-floor budget table (σ_V≈1.31 mV, dV/dt, slew~2 ps, cell-width<1 ps, inter-chip) —
  bounds detector vs shower; the cleanest argument the floor is NOT instrumental (so more light can't move it).
- **D.** FCC landscape (currently absent): CMS MTD BTL ~30 ps, HGCAL EM ~20–30 ps, RADiCAL 25 ps. The
  differentiator: RADiCAL times from a *sampling EM calorimeter at shower max* — timing+energy in the SAME
  ultra-compact (~25 X₀, ~13.5 cm) rad-hard detector, vs MTD as a dedicated separate layer.
- **E.** Route below the floor → Paper 2 (depth-corrected estimator → toward the ~2 ps slew limit, the
  5–10 ps FCC regime). Include AT MOST one 150-GeV depth-corrected demonstration point as a forward
  pointer if cheap; else state the falsifiable prediction. Don't over-promise.
- **Figure [REBUILD] floor-form discriminator:** floor_model DSB1 (1/√E χ²≈2 vs 1/E ≈12) paired with LuAG
  (20 vs 39, data never reaches 20). Can be an inset on V.B. **[OPT] FCC comparison table + budget table.**

## VII. CONCLUSIONS
- Land the two bookends. Verdict matched to the softened title ("dominant lever," not "sets/not species").
  Confirms (not revises) 17.5 ps & a/√E; edge recovery is the enabling technique; best measured 25.4 ps
  (typical ~30 ps stated alongside). Close: recovered edge → shared floor → "the floor is the shower
  itself, and that's where Paper 2 begins."

## ADMIN (blocking): full 39-author roster + 11 affiliations from [2]; funding; add refs (CMS HGCAL,
CMS MTD BTL, an FCC-hh requirements doc, a LuAG *irradiation* measurement).

---

## (d) PRIORITIZED FIGURE LIST
1. **MIXED in-module head-to-head** (mixed_h2h) — [HAVE→PROMOTE §V.C] — the kill shot; confound-free.
2. **Thesis σ(E) two-panel WITH ERROR BARS** (rebuild thesis + light-yield companion) — [REBUILD, FATAL].
3. **pulse.png** (recovered edge) — [HAVE, add "3.7× steeper" gap] — the jewel.
4. **Apparatus 3D schematic** — [BUILD] — highest-value missing.
5. **method.png** (cfd05 vs lgcfd) — [HAVE, paired errors].
6. **clip.png** (820 mV pile-up) — [HAVE, shrink title] — the villain.
7. **dist.png** (Gaussians + satellite inset) — [HAVE, FIX title + ADD satellite].
8. **hglg.png** — [REBUILD] — only unpublishable figure (garbled axes/debug text → one clean channel).
9. **HG/LG cartoon** — [BUILD] — prime the reader.
10. **Floor-form discriminator** (floor_model DSB1+LuAG) — [REBUILD] — forecloses photostat-vs-slew honestly.
11. **FCC comparison + instrumental-budget tables** — [BUILD, text/table].

## (e) MUST-DO-BEFORE-SUBMISSION (by rejection risk)
**FATAL:** (1) error bars on every result point + Table 1; (2) honest floor (both forms per build w/
χ²/p; b with 1/√E↔1/E systematic; LuAG never reaches 20 ps; headline measured 25.4); (3) systematics
table (δσ¹⁵⁰,δa,δb vs fit-window, brightest-N, source, HG/LG slope ±5%, beam momentum, MCP/inter-chip).
**MANDATORY:** (4) satellite demo + Δσ; (5) resolve MIXED (clean in timing kill-shot; define "LYSO-based"
by an objective LG-reliability criterion, not by which floors agree); (6) N/point + σ_stat + OOS number;
(7) method-gain paired significance + reconcile the four floor numbers; (8) depth-floor quantified;
(9) confront LuAG kinetics confound (answered by the kill-shot); (10) FCC quantitative context + MTD/HGCAL.
**EDITORIAL:** (11) disclose all cuts; (12) instrumental budget; (13) LuAG irradiation ref; (14) author
block; (15) typical-shower ~30 ps alongside 25.4.

## PANEL TENSIONS — RESOLVED
- Bold "sets/not species" vs kinetics confound vs caution → **soften to "governs/dominant lever"**, but it
  survives because the **kill-shot** defends it.
- Light-yield figure (confounded) vs the clean test → **include it, captioned honestly; kill-shot is the clean version.**
- Floor 20 ps clean vs fit-form artifact → **referee wins presentation (form-dependent extrapolation w/
  systematic); author/FCC win interpretation (shower-depth → material-independence is a prediction).**
- Depth-corrected point in Paper 1 vs 2-paper plan → **defer to Paper 2; at most one forward-pointer point.**
- Method ~40% of paper vs thesis is headline → **compress §III, expand §V (kill-shot + honest floor).**
- "25.4 best result" vs "not new" → **never sell 25.4 as a record; headline the two genuine novelties.**

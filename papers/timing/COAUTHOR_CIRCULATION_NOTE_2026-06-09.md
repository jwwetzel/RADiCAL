# Coauthor Circulation Note — RADiCAL timing paper (2026-06-09)

**Title:** *A same-shower comparison of wavelength shifters and light-yield limits in a compact
RADiCAL electromagnetic-shower timing module*
**PDF:** `papers/timing/radical_timing.pdf` · **Commit:** see git log (`7cd8e1f` prose pass,
`8a44b35` circulation prep, + the formatting-polish commit following it). Every number traces to a
committed, audited analysis product (scripts under `papers/scripts/`, audits pre-registered, gate
history in `papers/memory_analysis_gates.md`).

**Production-formatting pass (2026-06-09, after `8a44b35`):** layout/readability only — no claim,
number, selection, or interpretation changed. All three table overflows fixed (the benchmark table
was running off the page; Table 1's last column was clipped; Table 2 was overprinting body text);
internal ROOT super-titles removed from the five campaign figures (the LaTeX captions carry the
titles — regeneration reproduced identical values, e.g. the MIXED ratio re-emitted 1.04±0.05);
the method-gain figure's axis-label collisions fixed; the two legacy figures' garbled title strips
cropped; "Appendix Appendix" typo fixed. An independent page-by-page verification fan-out was run
on the rebuilt PDF; its two serious findings (a clipped panel title in Fig. 7; appendix figures
illegible at single-column width) and several annotation collisions were fixed in a second round,
and its numbers-integrity check confirmed every authoritative value unchanged. Full issue
inventory + dispositions: `papers/timing/FORMAT_AUDIT_2026-06-09.md`.

**Appendix figure consistency (resolved 2026-06-10):** the verification fan-out had caught that
Fig. A.10's baked annotations were pre-fix (stale floor/systematics numbers from the retired
paperSystematics era). The figure has been **regenerated from the post-fix chain**
(`systematicsPostfix.C` now emits it from the same arrays that produce Table 2, so figure and
table agree by construction): nominals 25.7/30.3/39.6/44.4 ps, totals ±1.0/±1.1/±0.9/±1.9 ps,
in-event veto variations included, and no baked floor numbers — the caption points to Tables 1–2
as authoritative. No claim, table, or text number changed.

## Thesis (one paragraph)
The timing performance of the RADiCAL shower-maximum module is governed primarily by detected light
yield, module response, and estimator choice rather than by WLS material identity alone. Across the
four-build energy scan the DSB1 and LuAG:Ce stochastic terms differ by more than a factor of two
(203±6 → 440±18 ps·√GeV), while the same-shower MIXED comparison gives per-capillary timing widths
agreeing at the 10–20% level (σ_DSB1/σ_LuAG = 1.04 ± 0.05 under the primary srCFD estimator). This
combination disfavors intrinsic WLS re-emission kinetics as the dominant explanation for the
cross-build difference, within this geometry and at these light levels.

## What changed scientifically vs. older internal drafts
1. **The "0.99, χ²/ndf=0.4" MIXED ratio is retired.** It came from brightness-threshold corner
   labels that misassigned SE-D (the corner map is now confirmed by amplitude-independent
   pulse-shape discriminants). The replacement is the paired same-shower width comparison:
   1.04 ± 0.05 (energy-period scatter), estimator spread 0.80–1.16.
2. **All numbers are "post-fix":** a σ-estimator bug (fit window inflated by rare broken-timing
   outliers) was found and fixed in production; every σ(E) is now monotonic and the published cfd05
   values are unchanged. Headline: 25.7 ± 0.6 ps (was 25.3 pre-fix); a = 203±6 / 440±18 (was
   201/455); floors 18.8±0.8 / 24.6±3.3 / 26.0±1.8 (the "shared ~20 ps floor" is now the qualified
   "19–26 ps, mutually consistent with stated caveats").
3. **Method gain re-derived on identical events:** core-width gain 1.3 ps @150 GeV; tail-sensitive
   gain 4.0 ps [3.4–4.7] — the gain is clipping-induced tail suppression, demonstrated in
   Appendix B, including the honest dim-pulse reversal (srCFD is NOT universally better).
4. **New:** the upstream–downstream mean-asymmetry drift (−33.6 ± 2.9 ps/e-fold), consistent with
   shower-maximum migration — supporting evidence that the floor is longitudinal shower physics.

## Intentionally conservative claims (please don't strengthen)
"Disfavors … as the dominant cause, within this geometry and at these light levels" (never
"identical"/"proves"); floors as a qualified range (never unqualified "shared 20 ps"); "confirms
the published 17.5 ps" (never "revises"); the depth drift is an acceptance-conditional ensemble
statement (never a per-event z measurement); position/4D claims are reserved for the companion
paper. The wording law lives in `papers/memory_claims_and_forbidden_language.md`.

## Please inspect first
1. Fig. money plot (four-build σ_t(E)) + Table tab:builds — the paper's numerical spine.
2. Fig. mixed_h2h_corrected + §6 — the same-shower control (the paper's unique asset).
3. §5.3 + Appendix B — the srCFD method gain and its tail decomposition.
4. Table tab:syst — the post-fix selection systematics (incl. the new veto-window rows).
5. §7 — the depth drift and its deliberately limited reading.

## Known limitations
Electrons only, normal incidence, one module per build; 25 GeV exists for DSB1 only; MIXED
module-wide row is a reference (ill-posed; per-material treatment is the result); the MIXED ratio's
position-coupling dilution is stated but not yet quantified; LuAG floor is an extrapolation (best
measured point 44.4 ps); no irradiated-module beam data in this dataset (hardness statements are
component-level citations); GEANT4 campaign not yet run (floor decomposition is consistency, not
simulation-validated).

## Specific questions for coauthors
1. Title: keep the same-shower framing, or adopt series naming ("RADiCAL II: …")?
2. §6 same-shower control: is the 10–20%/1.04±0.05 presentation with the estimator-dependence
   panel persuasive, or should the per-material σ(E) curves carry the section?
3. Depth section: agreed to keep in main text (short), with calibration deferred to the companion?
4. Benchmark table vs narrative comparison (Discussion): preference?
5. Engineering extrapolation (LuAG light recovery → ~31–35 ps): include at this strength?
6. Anyone able to confirm the MIXED corner map from the 2023 logbook (belt-and-suspenders; the
   pulse-shape confirmation is in §6)?

## Please don't relitigate (unless you find a concrete error)
The estimator-fix and its validation (committed audit trail); the retirement of the 0.99 ratio
(label artifact, demonstrated); the choice of srCFD as primary with LED/cfd05 as diagnostics (the
per-regime rule is demonstrated in Appendix B); the brightest-1000 headline selection (locked,
OOS-validated, systematics quoted). All have pre-registered audits under `papers/scripts/`.

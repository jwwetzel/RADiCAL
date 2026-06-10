# Format Audit — radical_timing.pdf (production-formatting pass, 2026-06-09)

Scope per instruction: **layout, plot uniformity, table sizing, figure-title visibility, caption
spacing, PDF readability only.** No scientific claims, numbers, selections, or interpretation may
change. Baseline: commit `8a44b35`, clean tree, tectonic build OK (10 pages).

Build-log overfull boxes (tectonic --keep-logs):

```
Overfull \hbox (471.93pt too wide)  tab_benchmarks.tex lines 12--31      <- CRITICAL
Overfull \hbox ( 63.19pt too wide)  radical_timing.tex lines 321--333    <- Table 1 (tab:builds)
Overfull \hbox ( 32.98pt too wide)  tab_systematics.tex lines 18--39     <- Table 2 (tab:syst)
Overfull \hbox (  2.22pt too wide)  \output (page head)                  <- negligible, ignore
```

## Issues (T = table, F = figure, X = text), with disposition

### Tables

- **T1 — CRITICAL. Table 3 (`tab_benchmarks.tex`) overflows the page by 472 pt** (p. 8): the
  *Scale* and *Comparability* columns are entirely off the right page edge; *Characteristic
  result* is truncated mid-sentence ("30–60 ps per track (begin–end of…"). The comparability
  caveats — which must NOT be silently lost — are currently invisible.
  **Fix:** fixed-width `p{}` columns + `\footnotesize` + reduced `\tabcolsep` so all six columns
  wrap inside `\textwidth`. No cell content deleted; comparability text retained verbatim.
- **T2 — HIGH. Table 1 (`tab:builds`, manuscript lines 321–333) overflows its column by 63 pt**
  (p. 4): the $\sigma_t^{150}$ column is visually clipped at the right margin (values cut).
  **Fix:** `\footnotesize` + `\tabcolsep` 3 pt + `@{}` outer trim. Headers and all values unchanged.
- **T3 — HIGH. Table 2 (`tab_systematics.tex`) overflows by 33 pt** (p. 6): the LUAG column
  overprints the adjacent text column (its values render on top of body text — unreadable both).
  **Fix:** same treatment as T2.

### Figures — campaign scripts (`papers/scripts/`), regenerate after edit

- **F1 — Internal ROOT super-titles are illegibly small in the PDF and duplicate the LaTeX
  captions** on all five campaign figures: `thesis_postfix` (Fig. 5), `method_postfix` (Fig. 6),
  `mixed_h2h_corrected` (Fig. 7), `depth_dial` (Fig. 8), `satellite_removal` (Fig. B.11; its
  super-title additionally renders garbled/overlapping at the canvas top edge).
  **Fix (per instruction preference):** remove the `DrawSuperTitle` call from the paper-bound
  canvas in each script; the LaTeX caption carries the title. Plot content untouched.
- **F2 — `method_postfix` (Fig. 6, p. 5): bottom-panel annotation overlaps the y-axis labels.**
  **Fix:** move the annotation right/into the pad interior in `methodGainPostfix.C`.
- **F3 — `satellite_removal` (Fig. B.11, p. 10): panel-C pad title truncated** ("tail excess
  concentrates in clipped eve…"). **Fix:** shorten the pad label (caption keeps the full wording).
- **F4 — `depth_dial` (Fig. 8, p. 7): log-x axis labels collide ("90100").**
  **Fix:** disable `MoreLogLabels`/thin the label set on the x axis. Data and fit unchanged.
- **F5 — `mixed_h2h_corrected` (Fig. 7, p. 7) is a 3-panel figure squeezed into one column**;
  panel internals (legends, axis labels) are at the edge of legibility.
  **Fix:** promote the float to `figure*` (two-column span). Pure layout; same PNG content.

### Figures — legacy/diagnostic era (no committed paper-script generator; regeneration could
### silently change plotted content, so fix must be non-generative)

- **F6 — `hglg.png` (Fig. 2, p. 3): top title strip garbled** (overlapping diagnostic run-label
  text baked into the PNG). **Fix:** crop the title strip off the PNG (content panels untouched);
  caption carries the title. Recorded: this figure remains diagnostic-grade; a styled regeneration
  is post-circulation work (needs its original study script + data revalidation).
- **F7 — `dist.png` (Fig. 4, p. 4): same garbled tiny header strip.** Same crop fix.
- **F8 — `clip.png` (Fig. 1) and `pulse.png` (Fig. 3): diagnostic-style internal titles**
  (run/energy stamps). Legible and informative; LEFT AS-IS this pass (recorded for the eventual
  pre-submission figure refresh).
- **F9 — `optimization.png` (Fig. A.9) and `systematics.png` (Fig. A.10): small internal
  super-titles.** Legible at print size; captions already carry full descriptions. LEFT AS-IS
  this pass (same rationale as F8).

### Text

- **X1 — "Appendix Appendix A" duplication** (p. 6, manuscript line 411): `Appendix~\ref{app:opt}`
  where elsarticle's `\ref` already yields "Appendix A". **Fix:** drop the literal word.

### Style standardization decision

The shared style already lives in `lib/viz/PlotUtils.h` (`ApplyRADiCALStyle`); the paper-figure
convention adopted in this pass is: **no internal super-titles on paper-bound canvases — the LaTeX
caption is the title.** A separate `RadicalPaperStyle.h` was considered and skipped: it would
duplicate `ApplyRADiCALStyle` for five scripts that already share it. Convention recorded here and
in `papers/memory_paper2_timing.md`.

### Non-issues verified

Page 1 (title/abstract/footnotes) clean; caption spacing consistent (elsarticle defaults);
bibliography typesets cleanly (pp. 8–9); Fig. A.10 panels legible; money plot (Fig. 5) bands +
legend correct after the circulation-prep restyle; the 2.22 pt page-head overfull is cosmetic
noise.

## Fix order executed

T1 → T2 → T3 → X1 → F1–F4 (script edits + parallel re-run) → F5 (float promotion) → F6/F7 (crops)
→ rebuild → page-by-page re-inspection → compliance grep → docs/memory/snapshot → commit.

## Outcome (post-fix rebuild)

- **Overfull boxes: 4 → 1.** Only the 2.22 pt page-head box remains (accepted as cosmetic noise).
  The 472 pt benchmark, 63 pt builds-table, and 33 pt systematics overflows are gone; the benchmark
  table now wraps in six `p{}` columns at `\footnotesize` with every comparability caveat retained.
- **All five campaign figures regenerated content-identical** (spot-verified: MIXED ratio
  re-emitted 1.038 ± 0.049 with estimator spread 1.038/1.161/0.795; depth dial re-emitted
  −33.6 ± 2.9 ps/e-fold PASS; method-gain identical-event floors 19.9/21.8/20.3 ps unchanged).
  Super-titles removed (captions carry titles); method-gain y-title/label collision and the
  90/100 log-label collisions fixed; satellite panel-C label shortened to fit.
- **Legacy figures:** hglg.png and dist.png garbled title strips cropped (panels untouched);
  clip/pulse/optimization/systematics left as recorded (F8/F9).
- **Compliance grep clean:** no retired 0.99/χ²=0.4 language, no forbidden floor/ToA/5D phrasing;
  the "27.4→25.3" at manuscript line ~364 is the *clipped-subset identical-event* pair from the
  post-fix method-gain script — a legitimate value, not the retired pre-fix headline.
- Build: tectonic exit 0, 10 pages, ~523 KB. Independent page-by-page verification fan-out run on
  the rebuilt PDF before commit (see circulation note).

## Round 2 — independent verification fan-out (5 page-pair reviewers + integrity checker)

**Integrity: ALL authoritative numbers verified unchanged in the typeset PDF** (Tables 1–3,
§6 ratio, §7 slope) — zero mismatches. Formatting findings and dispositions:

FIXED (scripts re-run, round 2):
- SERIOUS, Fig. 7(A): pad title clipped at right edge → shortened to "A same-shower widths
  (srCFD)"; bottom stamp shortened ("pulse-shape map"); panel C frame ymax 1.32→1.45 so the
  cfd05 marker clears the annotation block.
- SERIOUS, Fig. A.10 (and A.9): 2×2 grids illegible at single-column print size → promoted to
  `figure*` at 0.8\textwidth (≈2.1× label size). Their internal super-titles/baked stamps remain
  (legacy generators; recorded under F9).
- Fig. 5: HGTD band label moved clear of the MIXED curve/marker; legend moved up-right
  (0.47–0.96 × 0.66–0.93) + frame ymax 80→85 so LuAG 75 GeV no longer intrudes.
- Fig. 6 bottom strip: y-range −1.5→−3 (50 GeV point was clipped by the frame), y-tick density
  halved (505), stamp moved off the LED trace and reworded to "core width" only — the old stamp
  paired the +1.3 ps core value with the tail-convention CI [3.4,4.7], reading as
  self-contradictory; the tail-sensitive gain + CI stay in §5.3 where the conventions are defined.
- Fig. 8: legend moved to the empty upper-right (data descend lower-right).
- Fig. B.11: panel-A legend moved off the histogram peak; tail-window label shortened ("±66 ps",
  caption defines it) and moved above the red line's reach; panel-C x-title dropped (bin labels
  self-identify; was colliding), annotation block moved to the empty upper-left, frame headroom
  1.35→1.6.
- Caption of Fig. 7 relabeled (A)/(B)/(C) to match the internal panel labels (was Left/Centre/Right).

ACCEPTED / RECORDED (not fixed this pass):
- Legacy Figs. 1–4 stamps and small per-panel fonts (F6–F8 disposition): hglg y-title clip
  ("HG peak (mV)" overruns pad top), stamps over the 820 mV band, pulse.png annotation/trace
  collisions, dist.png middle-row "−25." clipped tick — all in baked PNGs; styled regeneration is
  pre-submission work with data revalidation.
- Affiliation "RADiCAL Collaboration," trailing comma — part of the flagged PLACEHOLDER author
  block (completed by the collaboration before submission).
- Prose dash-style inconsistency (§2.1 spaced en-dashes), heading hyphenation ("sta-/bility"),
  page-boundary hyphenation, B.11 lone-float page, Fig. A.9 caption (a)–(d) vs unlabeled panels —
  cosmetic; deferred to the submission pass.
- Raster resolutions (~230–315 ppi) meet circulation needs; Elsevier ≥500 ppi line-art export is a
  submission-time item (vector PDF export already exists for thesis_postfix).

## CONTENT-CONSISTENCY DISCOVERY (out of this pass's scope — flagged, NOT fixed here)

**Fig. A.10's baked annotations are PRE-FIX.** The enlarged (figure*) rendering made them readable
enough to compare against Table 2: the figure annotates LuAG floor b = 19.8 ± 5.6 (the retired
shared-20-ps-era fit; the post-fix value is 24.6 ± 3.3), MIXED syst ±1.5 vs the post-fix ±0.9,
DSB1 b = 19.5 ± 1.1 vs 18.8 ± 0.8, and it lacks the in-event veto variation points added to the
post-fix budget. systematics.png (and likely optimization.png's exact curves) come from the
pre-fix paperSystematics-era generator and were missed by the stale-number purge, which targeted
text, tables, and the kill-shot figure. The fix — extending `systematicsPostfix.C` to emit the
A.10 stability figure from the post-fix shift table — changes figure numbers and therefore
belongs to a dedicated follow-up, not a formatting pass. Flagged in the circulation note so
coauthors read Table 2 (not Fig. A.10's annotations) as authoritative.

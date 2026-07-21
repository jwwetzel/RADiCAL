# analyze/ — the exploration record

This directory is the preserved research trail of the 2023 timing analysis: hypotheses,
diagnostics, dead ends, and the studies that *became* the published gates. It is kept
deliberately — the paper's credibility rests on being able to see the road, not just the
destination.

**Nothing in here is load-bearing for the paper except two figure generators:**
`studies/narrativeFigs.C` (Figs. 2 and 4) and `studies/hgLgPlot.C` (Fig. 3). Every other
published number and figure comes from a gate in `papers/scripts/<gate>/` (macro +
pre-registered AUDIT + committed log — see `papers/scripts/INDEX.md`). If a study here
contradicts a gated product, **the gate wins**; retirements and why are recorded in
`papers/memory_analysis_gates.md`.

## Layout

- `studies/` — ~120 single-purpose ROOT macros, indexed one line each in
  `studies/INDEX.md`. Three carry DEPRECATED banners pointing to their gate
  replacements (`methodCompare.C`, `mixedHeadToHead.C`, `paperSystematics.C`).
- `sigmaT.C` — **RETIRED** historical driver (Method A "best-bin", the 27.4 ps era);
  its header says so. The production headline lives in
  `papers/scripts/timing_fit_summary/timingFitSummary.C` via
  `rad::timingBrightestK` + `kLGCFD` (srCFD).
- `timingLadder.C`, `timingHeadline.C`, `slopeVsE.C` — exploratory drivers from the
  same era (the first two also drive the retired best-bin path).
- `makeReport.py` + `deployReport.sh` — static HTML report assembly for `site/`
  (the deployed pages carry HISTORICAL-PAGE banners; see `site/`).

## Running studies

`source setup.sh` at the repo root, then `root -l -b -q 'analyze/studies/<macro>.C+'`.
Most need `data/2023/reduced/` (not in git — `data/2023/README.md`). Studies write to
`output/` and `figures/` (both regenerable); they do **not** touch gated products.
Estimator background: `lib/physics/ESTIMATOR.md`; conventions: `docs/CONVENTIONS.md`.

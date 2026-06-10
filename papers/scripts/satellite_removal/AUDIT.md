# Satellite-Removal Demonstration Audit (pre-analysis)

Date: 2026-06-09. Written BEFORE the script.

## 1. Manuscript claim being supported (and nothing more)
"The srCFD estimator improves timing primarily by suppressing non-Gaussian satellite/tail structure
induced by clipping/saturation and fixed-fraction estimator bias, rather than by changing the
underlying detector response." This substantiates §5.3's two-convention split (core gain 1.3 ps vs
tail-sensitive gain 4.0 ps at 150 GeV) — the difference between those numbers IS the tail story.
NOT supported here: any new detector-performance claim, "saturation fully corrected", or universal
statements about CFD methods.

## 2. Event sample (identical to the §5.3 method-gain headline)
DSB1; identical-event intersection (all 8 ends valid + per-source median veto for EVERY estimator);
brightest-1000 by sum_lg selected ONCE; production cuts (wc_ok, MCP 200–750 mV, TimingFiducialR(E));
same channels. Anchor energy **150 GeV** (94% clip — where the recovery is most visible and where
the §5.3 numbers were installed). For the clipped/less-clipped split, 150 GeV's unclipped subset is
too small (N=46, established in the method-gain run), so the split panel uses **75 GeV** (per-channel
clip fraction ≈50% ⇒ both subsets populated); this fallback is pre-registered here, not improvised.

## 3. Pre-registered observable and metrics
Observable: the per-event **(DW−UP)/2** under each estimator — exactly the quantity whose width is
the production resolution (no visually-convenient surrogate). Distributions median-centered per
estimator. Metrics per estimator (cfd05, srCFD; LED tabulated for completeness):
- core Gaussian σ (FitGaussCore, 2σ core — the paper's width convention);
- production width (post-fix tebSigma);
- robust σ (IQR/1.349);
- **tail fraction** f_tail = fraction of events with |x − median| > W, with the COMMON window
  W = 2.5 × σ_core(srCFD@that energy) applied to BOTH estimators (estimator-independent;
  Gaussian expectation ≈ 1.2%);
- satellite identification: any resolved secondary lobe noted on the log-y distribution
  (the parent documented a ~0.2 ns satellite in MCP-referenced sums; in the reference-free
  difference it is expected to appear as tails/shoulders rather than a clean lobe — stated honestly);
- N (events).

Expected outcome (pre-registered): f_tail(cfd05) > f_tail(srCFD) at 150 GeV; in the 75 GeV split,
the tail excess concentrates in the clipped subset and largely vanishes in the less-clipped subset.

## 4. Outputs
`satelliteRemoval.C` (this dir) → `papers/figures/satellite_removal/satellite_removal.{png,pdf}`
(A: identical-event (DW−UP)/2 distributions, unit-normalized, log-y, common tail window marked;
B: metric comparison (core σ / production width / tail fraction);
C: 75 GeV clipped vs less-clipped tail fractions), sidecar caption
`satellite_removal_CAPTION.txt`, log `satellite_removal_result.txt` (this dir).
§5.3 patch replaces the TODO-P2 satellite marker.

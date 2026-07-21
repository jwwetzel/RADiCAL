# The timing estimator, end to end

(CODE_AUDIT_2026-07-21 deep item D3 — the 30-line walk of how a published σ_t is built.
Code: `RadTiming.h`; thresholds: `SelectionCuts.h`; stats: `../../docs/STATS_CONVENTIONS.md`.)

**0. Naming bridge.** The paper's **srCFD** is the reduced-tree branch **`hg_lgcfd`**,
exposed as `RadView::kLGCFD`: a CFD run on the measured HG waveform at a threshold of
0.15 × the *LG-predicted true* (unclipped) peak, referenced to the MCP. Adopted sources
are per-regime (paper Sec. 5.3): srCFD for the bright/clipped builds (DSB1, MIXED), LED
for the dim builds (LuAG, TENERGY); `cfd05` is the clipped-peak diagnostic class.

**1. Per-event kernel — `eventDWUP` (RadTiming.h:97–110).** Each valid end (timing role,
`hg_peak ≥ kHG_minPeak`) contributes a time from the chosen source. Ends within
`kTimingChanConsistency_ns` = 2 ns of the **event median** survive (the in-event
broken-timing veto; a no-op for clean events); survivors split into Down (ends 0–3) and
Up (ends 4–7), and the event's observable is **t = (mean_DW − mean_UP)/2**. Needs ≥1
survivor on each side. The D4U4 end order this hardcodes is guarded at config load
(BuildConfig.h:147–158).

**2. Reference cancellation.** Every end's time is measured against an MCP reference on
the *same DRS4 group*, so the reference and group clock cancel **in the DW−UP
difference** — (DW−UP)/2 is reference-free by construction, immune to MCP jitter and
inter-group sync drift. The one cross-group end (DSB1: SW-U, the sole end on group 1)
uses its own group's copy, `mcp_ref: 2` in `data/2023/configs/DSB1.json`, preserving the
cancellation; it is guarded only by the reducer validity floor plus the in-event veto
(SelectionCuts.h note). Studies `mcpJitter.C`/`trSync.C` measured what the common
reference would otherwise contribute.

**3. Width estimator — `tebSigma` (RadTiming.h:36–76).** A robust seed first: iterated
±2.5σ truncated RMS, de-biased by 0.9546 (single-pass Gaussian truncation factor). Then
a 120-bin **Gaussian-core fit** seeded from it — the production width. Guards: heavy
tails (|skew| > 3 or kurtosis > 20) or a fit outside 0.5–2× the robust seed fall back to
the robust value. Caveat (lines 47–52): the iterated window's fixed point is ±2.34σ, so
the *fallback* de-bias sits ~2% low — fallback path only, headline unaffected,
post-submission fix.

**4. Production headline — `timingBrightestK` (RadTiming.h:170–190).** Event selection:
`wc_ok`, MCP1 ∈ (200, 750) mV, track within `TimingFiducialR(E)` (2.5 mm ≤ 100 GeV →
3.0 mm ≥ 125, linear ramp; SelectionCuts.h:161–165) of the sum_lg-weighted beam
centroid, and a valid kernel event. Of those, the **brightest K = 1000 by `sum_lg`**
(nth_element) go into `tebSigma`. A bare `s > 12.0` ps sanity floor (line 188) rejects
sub-physical widths, mirroring the best-bin guard rationale at lines 138–144. Result:
**σ_t(150 GeV, DSB1, srCFD) = 25.7 ± 0.6 ps**; the same chain with no brightness cut is
the ≈50 ps full-fiducial companion (`papers/scripts/full_fiducial_check/`).

`timingBestBin` (RadTiming.h:112–161) is **RETIRED** — the historical Method A
equal-population best-bin, kept for the exploration record; the paper's numbers come
from `timingBrightestK` via `papers/scripts/timing_fit_summary/`.

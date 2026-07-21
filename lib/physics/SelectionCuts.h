// ============================================================================
// SelectionCuts.h — RADiCAL analysis selection cuts (single source of truth)
// ============================================================================
//
// All numerical thresholds used for event and channel selection are defined
// here.  Every analysis macro must reference these named constants — no magic
// numbers anywhere else in the analysis chain.
//
// To study systematic uncertainties from a cut, change the value here and
// re-run the affected step.  The comment beside each cut explains WHY the
// value was chosen, so reviewers know what to change and why.
//
// ── Calorimeter geometry ────────────────────────────────────────────────────
//
//   kCalo_x0 / kCalo_y0  Nominal center of the calorimeter face in wire-
//                         chamber coordinates [mm].  Used only for annotation;
//                         the actual fiducial center is computed per-run from
//                         the signal-weighted beam centroid in ScanRunCenters().
//
// ── Fiducial selections ─────────────────────────────────────────────────────
//
//   Two radii are defined because the optimum trade-off between purity and
//   statistics differs between the two measurement goals:
//
//   kFiducial_r_energy  Tight: minimises shower leakage and lateral
//                        non-uniformity that bias the energy sum.
//
//   kFiducial_r_timing  Looser: timing is insensitive to shower leakage
//                        (it depends on the leading edge, not the integral),
//                        so we accept more events for a better σ_t estimate.
//
//   Both radii are centred on the DATA-DERIVED centroid, not kCalo_x0/y0.
//
// ── MCP reference quality ───────────────────────────────────────────────────
//
//   kMCP1_minPeak / kMCP2_minPeak
//       Applied for TIMING analysis.  A low MCP amplitude means the particle
//       was only partially captured by the MCP; the CFD crossing is noisy
//       and the timing reference becomes unreliable.  100 mV is ≈3× the
//       typical MCP noise floor.
//
//   kMCP1_maxPeak / kMCP2_maxPeak
//       Applied for TIMING analysis (upper bound).  When the MCP signal
//       exceeds the DRS4 full-scale range the waveform clips at the ADC
//       rail.  ExtractPulse sees peak ≈ pedestal ≈ 800 mV and computes a
//       20% CFD threshold that is too LOW relative to the true (unknown,
//       larger) amplitude.  The crossing is found earlier on the rising
//       edge than the true 20% point, biasing mcp_cfd20 early and all
//       hg_cfd[i] = cfd20[i] − mcp_cfd20 late.  These events are excluded
//       from timing analysis only; they are valid beam particles for energy.
//       The 750 mV ceiling sits below the observed saturation pile-up
//       (~800 mV) and well above the physical signal distribution tail.
//
//   kMCP_minPeak_E
//       Applied for ENERGY analysis.  Looser cut — we only need to confirm
//       that a beam particle passed; we don't use the MCP time.
//
// ── Signal-weighted centroid pre-scan ───────────────────────────────────────
//
//   kSumLG_centroid  Minimum total LG signal required for an event to
//                     contribute to the beam centroid computation.  Rejects
//                     sub-threshold noise events and beam halo.
//
// ── HG timing channel quality ───────────────────────────────────────────────
//
//   kHG_minPeak       Per-channel minimum HG amplitude [mV] required before
//                      using a timing measurement.  20 mV is ≈4× the typical
//                      HG noise floor (~5 mV RMS).  Events below this cut
//                      have poor LED/CFD precision due to noise walk.
//
//   kHG_minPeak_prescan  Stricter threshold used only in the pre-scan pass
//                         that computes the CFD histogram window center.
//                         The tighter cut gives a cleaner mean estimate.
//
//   kHG_LED_thresh    Absolute voltage threshold for LED and TOT extraction
//                      [mV above pedestal].  Chosen to match kHG_minPeak so
//                      that every event passing the quality cut also passes
//                      the LED threshold — no efficiency loss.
//
// ── LG energy channel quality ───────────────────────────────────────────────
//
//   kLG_minPeak  Per-channel minimum LG amplitude [mV].  5 mV is just above
//                 the LG noise floor; used for the HG/LG ratio correction and
//                 to build the energy sum with good channels only.
//
// ── Wire chamber ────────────────────────────────────────────────────────────
//
//   kWC_minPeak  Minimum WC pulse amplitude [mV] for a plane to be counted
//                 as valid.  Applied in reduce/Reducer.C (wc_ok branch;
//                 historically processRun.C).
//                 Must have all four planes (R, L, D, U) valid for wc_ok.
//
// ── Summary of event selection — WHICH CHAIN APPLIES WHAT ───────────────────
//
//   The PRODUCTION TIMING chain (the paper numbers; RadTiming.h timingBrightestK)
//   applies: 1, 2 (MCP1 window only), 3 (TimingFiducialR), 5, 6, plus the
//   in-event median-consistency veto (kTimingChanConsistency_ns) and the
//   brightest-K selection (K=1000 by sum_lg). Steps 4, 7, 8 and the MCP2/energy
//   variants are ENERGY- or STUDY-chain only — NOT applied by production timing.
//
//   1. wc_ok == true                      [timing + energy]
//   2. mcp1_peak in (kMCP1_min, kMCP1_max) [timing]
//      mcp_peak > kMCP_minPeak_E          [energy chain, lower only]
//   3. r_beam  < TimingFiducialR(E)       [timing: 2.5 mm ≤100 GeV → 3.0 mm ≥125]
//      r_beam  < kFiducial_r_energy       [energy, tighter]
//      where r_beam is computed from the data-derived centroid, not kCalo_x0/y0
//   4. sum_pb < kPb_maxRatio * sum_lg     [energy/uniformity STUDIES only]
//      (only applied when sum_lg > kSumLG_centroid to avoid divide-near-zero)
//   5. hg_peak[i] > kHG_minPeak           [timing + energy, per channel]
//   6. hg_cfd[i]  > kNoTime sentinel      [timing, per channel: crossing found]
//   7. hg_tot[i]  > 0                     [RETIRED legacy TOT study ("M5")]
//   8. lg_peak[i] > kLG_minPeak           [RETIRED legacy LG study ("M7")]
//
// ============================================================================

#ifndef SELECTIONCUTS_H
#define SELECTIONCUTS_H

// ---------------------------------------------------------------------------
// Calorimeter geometry  [mm]
// ---------------------------------------------------------------------------
static const double kCalo_x0 = 6.6;  // face center x in WC coords (data-derived
                                     // from shashlik edges, moduleCenter.C;
                                     // combined 6.60, energy-stable 6.4-6.7)
static const double kCalo_y0 = 4.7;  // face center y in WC coords (moduleCenter.C
                                     // combined 4.70, energy-stable 4.6-4.75)

// ---------------------------------------------------------------------------
// Fiducial radii  [mm]
// ---------------------------------------------------------------------------
static const double kFiducial_r_energy = 2.0;  // tight fiducial for energy analysis
static const double kFiducial_r_timing = 3.0;  // nominal/loose timing fiducial (display + default)

// Per-energy OPTIMISED timing fiducial radius [mm], set by the OUT-OF-SAMPLE-
// validated fiducial-radius scan (fiducialTimingScan.C + the run-folded OOS in
// timingEnergyBins.C; documented in the Layer-3 fiducial-optimisation chapter):
//
//   25-100 GeV  ->  2.5 mm   The outer-ring position/time-walk degradation
//                            dominates; tightening robustly improves the OOS
//                            best-bin sigma_t (50/75/100 GeV better by 1.7-3.6 ps;
//                            25 GeV ~tied).
//   125-150 GeV ->  3.0 mm   At 125 GeV the tighter cut OVERFITS (in-sample
//                            27.1 ps but OOS 33.1 ps); at 150 GeV the 4x larger
//                            sample lets the single-best-E_meas-bin selection
//                            benefit from the looser pool (OOS 27.4 vs 28.1 ps).
//                            Both are OOS-validated at 3 mm.
//
// Physics-motivated, OOS-validated endpoints -- NOT a per-energy argmin (the
// per-energy best-bin curve is too jumpy, and tighter cuts overfit in-sample).
//
// CONTINUITY: the former hard 2.5->3.0 STEP at 112 GeV created a fiducial seam --
// 125/150 GeV abruptly admit the noisy 2.5-3.0 mm outer annulus that lower
// energies exclude, contributing a non-monotonic bump at the 100->125 boundary
// (sigma-monotonicity audit, 2026-06). We keep the two OOS-validated endpoints
// (2.5 mm at <=100 GeV, 3.0 mm at >=125 GeV) but RAMP linearly between them so the
// radius is continuous in energy. At every MEASURED energy this is identical to the
// old step (25-100 ->2.5; 125,150 ->3.0), so the published headline is unchanged;
// it only removes the discontinuity for interpolated energies. (The seam's residual
// effect on the measured points is removed upstream by the in-event broken-timing
// veto + robust estimator window, not by moving the validated radii.)
static inline double TimingFiducialR(double energy_GeV) {
    if (energy_GeV <= 100.) return 2.5;
    if (energy_GeV >= 125.) return 3.0;
    return 2.5 + 0.5*(energy_GeV-100.)/25.;    // linear 2.5->3.0 across 100..125 GeV
}

// ---------------------------------------------------------------------------
// MCP reference quality  [mV]
// ---------------------------------------------------------------------------
static const float kMCP1_minPeak    = 200.0f;  // MCP1 amplitude cut for timing (lower)
static const float kMCP1_maxPeak    = 750.0f;  // MCP1 saturation cut for timing (upper)
static const float kMCP2_minPeak    = 200.0f;  // MCP2 amplitude cut (NOT applied in production —
static const float kMCP2_maxPeak    = 750.0f;  //   no macro cuts mcp2_peak; SW-U's MCP2 reference is
                                               //   guarded only by the reducer's 30 mV validity floor
                                               //   + the in-event veto. Kept for symmetry/studies.
static const float kMCP_minPeak_E   =  50.0f;  // MCP amplitude cut for energy (lower only)

// ---------------------------------------------------------------------------
// Pre-scan centroid thresholds
// ---------------------------------------------------------------------------
static const float kSumLG_centroid   = 300.0f;  // min sum_lg for centroid inclusion [mV]
static const float kHG_minPeak_prescan = 30.0f; // stricter HG cut used in pre-scan only [mV]

// ---------------------------------------------------------------------------
// HG timing channel quality  [mV]
// ---------------------------------------------------------------------------
static const float kHG_minPeak    = 20.0f;  // min HG amplitude for timing use
static const float kHG_LED_thresh = 20.0f;  // LED/TOT absolute threshold (= kHG_minPeak)
static const float kHG_maxPeak = 950.0f;  // DRS4 near-saturation threshold [mV] (diagnostic studies
                                          // ONLY — not applied in production; the 2023 HG chain
                                          // physically clips at ~820 mV, see BuildConfig hg_sat_mV)

// In-event broken-timing veto [ns]: a capillary-end crossing that disagrees with
// the event's MEDIAN end-time by more than this is a wrong-feature / near-noise
// crossing (e.g. a channel that triggered ~30 ns off on a pre-pulse, or a near-
// threshold 20-100 mV pulse picked up because brightness was ranked on LG). Such
// ends are dropped from that event's (DW-UP)/2 (the other ends still time it).
// Legitimate end-to-end spread is ~0.5 ns (shower depth), so 2 ns is wide open for
// real pulses and only removes genuine outliers -> no-op for clean events. This is
// the deterministic veto that de-fuels the estimator-window distortion behind the
// non-monotonic sigma_t(E) (sigma-monotonicity audit, 2026-06).
static const float kTimingChanConsistency_ns = 2.0f;

// ---------------------------------------------------------------------------
// LG energy channel quality  [mV]
// ---------------------------------------------------------------------------
static const float kLG_minPeak = 5.0f;   // min LG amplitude for energy / ratio use

// ---------------------------------------------------------------------------
// PbGlass reference quality  [mV]
// ---------------------------------------------------------------------------
static const float kPb_minPeak  = 5.0f;   // min PbGlass amplitude per PMT

// Shower containment: reject events where PbGlass signal is a large fraction of
// the RADiCAL signal — these are punch-through or edge-shower events where the
// shower is not fully contained in the RADiCAL module.
// Physical model: well-centred, fully-contained showers should deposit negligible
// energy in the Pb Glass behind the module.  Events with sum_pb > kPb_maxRatio *
// sum_lg are excluded from timing and energy analysis.
// Tune by inspecting Page 10 of quality_report.pdf (sum_pb vs sum_lg scatter).
static const float kPb_maxRatio = 0.30f;  // max sum_pb / sum_lg for containment cut

// ---------------------------------------------------------------------------
// Wire chamber  [mV / mm]
// ---------------------------------------------------------------------------
static const float  kWC_minPeak = 20.0f;  // min WC pulse amplitude per plane

// Position histogram bin width [mm] for all 2D/1D WC track plots.
// Set to match the effective WC position resolution (~1 mm, determined by
// the delay-line wire pitch of the CERN SPS test-beam chambers).
// Using finer bins reveals DRS4 cell-width non-uniformity as periodic
// artifacts; using this value ensures one bin ≥ one wire pitch.
// Change here to propagate consistently to all position histograms.
static const double kWC_resBin = 1.0;  // mm per bin for WC position plots

#endif // SELECTIONCUTS_H

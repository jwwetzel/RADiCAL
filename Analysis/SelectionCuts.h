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
//                 as valid.  Applied in processRun.C (wc_ok branch).
//                 Must have all four planes (R, L, D, U) valid for wc_ok.
//
// ── Summary of event selection (applied sequentially) ───────────────────────
//
//   1. wc_ok == true                      wire-chamber track reconstructed
//   2. mcp_peak > kMCP1_minPeak           clean MCP1 reference (timing, lower)
//      mcp_peak < kMCP1_maxPeak           MCP1 not saturated (timing, upper)
//      mcp_peak > kMCP_minPeak_E          MCP confirmed (energy, lower only)
//   3. r_beam  < kFiducial_r_timing       beam within timing fiducial
//      r_beam  < kFiducial_r_energy       beam within energy fiducial (tighter)
//      where r_beam is computed from data-derived centroid, not kCalo_x0/y0
//   4. sum_pb < kPb_maxRatio * sum_lg     shower contained in RADiCAL
//      (only applied when sum_lg > kSumLG_centroid to avoid divide-near-zero)
//   5. hg_peak[i] > kHG_minPeak           per-channel: signal above noise
//   6. hg_cfd[i]  > kNoTime sentinel      per-channel: CFD crossing found
//   7. hg_tot[i]  > 0                     per-channel: TOT valid (M5 only)
//   8. lg_peak[i] > kLG_minPeak           per-channel: LG valid (M7 only)
//
// ============================================================================

#ifndef SELECTIONCUTS_H
#define SELECTIONCUTS_H

// ---------------------------------------------------------------------------
// Calorimeter geometry  [mm]
// ---------------------------------------------------------------------------
static const double kCalo_x0 = 6.5;  // nominal face center x in WC coords
static const double kCalo_y0 = 4.5;  // nominal face center y in WC coords

// ---------------------------------------------------------------------------
// Fiducial radii  [mm]
// ---------------------------------------------------------------------------
static const double kFiducial_r_energy = 2.0;  // tight fiducial for energy analysis
static const double kFiducial_r_timing = 3.0;  // loose fiducial for timing analysis

// ---------------------------------------------------------------------------
// MCP reference quality  [mV]
// ---------------------------------------------------------------------------
static const float kMCP1_minPeak    = 200.0f;  // MCP1 amplitude cut for timing (lower)
static const float kMCP1_maxPeak    = 750.0f;  // MCP1 saturation cut for timing (upper)
static const float kMCP2_minPeak    = 200.0f;  // MCP2 amplitude cut for timing (lower)
static const float kMCP2_maxPeak    = 750.0f;  // MCP2 saturation cut for timing (upper)
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
static const float kHG_maxPeak = 950.0f;  // DRS4 near-saturation threshold [mV] -- above this, waveform may clip

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

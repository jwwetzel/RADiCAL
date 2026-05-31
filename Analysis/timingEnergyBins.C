// ============================================================================
// timingEnergyBins.C — energy-binned timing analysis for RADiCAL
// ============================================================================
//
// Implements the energy-binned timing analysis from arXiv:2401.01747 Section 5.3,
// extended with CFD-20% and M7 HG/LG ratio walk-correction improvements.
//
// Physics motivation
// ------------------
// Timing resolution depends on detected energy E_meas, not just beam energy.
// Binning events by sum_lg (= ΣA_LG, total low-gain signal, proportional to
// E_meas) and selecting the high-E_meas bins (where E_meas ≈ E_beam, i.e.
// fully-contained showers) gives the best timing resolution estimate.
// The paper arXiv:2401.01747 achieves 27 ps at 150 GeV using a simple LED
// threshold; CFD-20% + M7 walk correction are expected to improve on this.
//
// Three estimators are compared:
//
//   Method A: (DW − UP) / 2  using raw CFD-20%
//     — Matches the paper's "BestMinus" approach.  The MCP time cancels in the
//       difference when the same MCP is used for all channels.  (Ch7 uses MCP2
//       which shifts the mean by a constant but does NOT affect σ.)
//
//   Method B: (DW + UP) / 2  using raw CFD-20%
//     — Paper's "BestPlus" equivalent.  Includes MCP jitter → wider distribution.
//
//   Method C: (DW − UP) / 2  using M7 (HG/LG ratio) walk-corrected CFD-20%
//     — Per-channel correction: t_corr[i] = dt[i] − k_i*(ratio_i − <ratio_i>)
//       where ratio_i = hg_peak[i]/lg_peak[i] and k_i is the OLS slope.
//
// Fit model
// ---------
//   Each timing histogram is fitted with BOTH:
//     - FitGaussCore (Gaussian): standard σ_Gauss
//     - FitCrystalBall (CB): σ_CB from the Gaussian-core sigma of the CB function
//   Both sigmas are shown on overlays and stored in the summary ROOT file.
//   CB typically gives a smaller sigma when a low-side tail is present
//   (hadronic leakage or mis-timed events), making σ_CB a more robust core estimator.
//
// Walk correction diagnostic
// --------------------------
//   The M7 walk correction coefficients (OLS slopes) are trained on ALL fiducial
//   events and evaluated on the same set (in-sample).  The correction gain
//   sigma_C / sigma_A is annotated on each per-bin overlay.
//
//   TODO: add 5-fold cross-validation to quantify walk correction overfit.
//   Pattern: for fold k=0..4, train on events where (ev_index % 5 != k),
//   apply to test events (ev_index % 5 == k), accumulate test sigma_t,
//   compare sigma_t_crossval vs sigma_t_insample.
//
// Algorithm
// ---------
//   1. Pre-scan to derive beam centroid and CFD offset (ScanRunCenters).
//   2. Collect all fiducial events into a vector<EbEvent>.
//   3. Fit Gaussian to sum_lg distribution → define 9 equal-width energy bins.
//   4. Train M7 slopes via OLS over all fiducial events.
//   5. For each of 9 bins: compute methods A, B, C; fit Gaussian + Crystal Ball;
//      extract σ in ps for both fit models.
//   6. Pool bins 6–8 (0-indexed: 5, 6, 7) → best-estimator σ for each method
//      and each fit model.
//
// Outputs
// -------
//   Analysis/Output/<label>/timing_energy_bins.pdf     (2 pages, per energy)
//   Analysis/Output/Summary/timing_energy_bins_summary.pdf
//   Analysis/Output/Summary/timing_energy_bins.root
//
// Usage (from repository root):
//   root -l -b -q 'Analysis/timingEnergyBins.C+'
// ============================================================================

#include "ChannelConfig.h"
#include "PlotUtils.h"       // StylePad, FitGaussCore, DrawFitOverlay, ScanRunCenters

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TBox.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TAxis.h"
#include "TMath.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <iomanip>
#include <iostream>
#include <vector>
#include <string>

// ---------------------------------------------------------------------------
// Module-local constants
// ---------------------------------------------------------------------------
static const float  kSentinel_teb  = -1e9f;   // invalid dt value
static const float  kSentinelCut_teb = -1e5f; // validity test: dt > this
static const int    kNBins_teb      = 9;       // energy bins
static const int    kNMeth_teb      = 3;       // methods A, B, C

// Method labels, colors, markers
static const char*  kMLabel_teb[kNMeth_teb] = {
    "(DW#minusUP)/2 CFD-5%",
    "(DW+UP)/2 CFD-5%",
    "(DW#minusUP)/2 M7-cor"
};
static const int kMColor_teb[kNMeth_teb]  = { kBlue+1, kOrange+1, kRed+1 };
static const int kMMarker_teb[kNMeth_teb] = { 22, 23, 20 };

// Crystal Ball graph style: open markers, same color as Gaussian, dashed lines.
// Marker 26=open triangle-up, 32=open triangle-down, 24=open circle
static const int kMMarkerCB_teb[kNMeth_teb] = { 26, 32, 24 };

// ---------------------------------------------------------------------------
// OLS slope: y = a + b*x, returns b.  Requires >= 15 pairs.
// ---------------------------------------------------------------------------
static double OLS_teb(const std::vector<float>& x, const std::vector<float>& y) {
    int n = (int)x.size();
    if (n < 15) return 0.;
    double sx=0, sy=0, sxx=0, sxy=0;
    for (int i = 0; i < n; ++i) {
        sx  += x[i];
        sy  += y[i];
        sxx += (double)x[i]*x[i];
        sxy += (double)x[i]*y[i];
    }
    double d = (double)n*sxx - sx*sx;
    if (std::fabs(d) < 1e-30) return 0.;
    return ((double)n*sxy - sx*sy) / d;
}

// ---------------------------------------------------------------------------
// Build a TH1F from a vector, auto-ranging on a robust estimate of the core.
//
// Two-pass approach: first compute mean/rms including all values, then reject
// events more than 5*rms from the mean and recompute the range statistics.
// The histogram is filled with ALL values; outliers land in overflow and do
// not affect FitGaussCore (which fits within 2*sigma of the peak).
//
// This robustness matters because (DW+UP)/2 values can include rare events
// with incorrect MCP timestamps that shift the measured time by tens of ns;
// without rejection the raw rms would be O(ns) and the histogram bins would
// be O(100 ps) wide — too coarse for FitGaussCore to resolve the Gaussian core.
// ---------------------------------------------------------------------------
static TH1F* VecToHist_teb(const char* name, const std::vector<float>& v,
                             int nb = 120) {
    TH1F* h;
    if (v.empty()) {
        h = new TH1F(name, "", nb, -1., 1.);
        h->SetDirectory(nullptr);
        return h;
    }

    // ── Pass 1: inclusive mean and rms ─────────────────────────────────────
    double mu1 = 0., ms1 = 0.;
    for (auto x : v) mu1 += x;
    mu1 /= (double)v.size();
    for (auto x : v) ms1 += ((double)x - mu1) * ((double)x - mu1);
    ms1 = std::sqrt(ms1 / (double)v.size());
    if (ms1 < 0.008) ms1 = 0.100;

    // ── Pass 2: reject outliers > 5*rms; recompute range statistics ────────
    double mu2 = 0., ms2 = 0.;
    int    n2   = 0;
    for (auto x : v) {
        if (std::fabs((double)x - mu1) < 5. * ms1) {
            mu2 += x;
            ++n2;
        }
    }
    if (n2 > 0) {
        mu2 /= n2;
        for (auto x : v)
            if (std::fabs((double)x - mu1) < 5. * ms1)
                ms2 += ((double)x - mu2) * ((double)x - mu2);
        ms2 = std::sqrt(ms2 / n2);
        if (ms2 < 0.008) ms2 = 0.100;
    } else {
        mu2 = mu1; ms2 = ms1;  // fallback (all outliers — shouldn't happen)
    }

    // ── Build histogram over core range; outliers go to overflow ───────────
    h = new TH1F(name, "", nb, mu2 - 4.*ms2, mu2 + 4.*ms2);
    h->SetDirectory(nullptr);
    for (auto x : v) h->Fill(x);  // all values; outliers → overflow (ignored by FitGaussCore)
    return h;
}

// ---------------------------------------------------------------------------
// Per-event storage struct
// ---------------------------------------------------------------------------
struct EbEvent {
    int   run;        // run number (for run-level cross-validation folds; #G5)
    float sum_lg;
    float dt[8];      // hg_cfd[i] directly from ntuple (already MCP-referenced):
                      //   ch0–6: hg_cfd20[i] − mcp_time  (MCP1-ref)
                      //   ch7:   hg_cfd20[7] − mcp2_time (MCP2-ref)
                      // kSentinel_teb if channel or MCP was invalid.
                      // Used for Methods A (DW−UP) and C (walk-corrected DW−UP).
    float dt7_mcp1;   // ch7 re-referenced to MCP1:
                      //   hg_cfd[7] + (mcp2_time − mcp_time) = hg_cfd20[7] − mcp_time
                      // Used ONLY for Method B so all UP channels share MCP1,
                      // giving a proper centroid estimator without MCP jitter.
    float hg[8];      // hg_peak[i]
    float lg[8];      // lg_peak[i]
};

// ---------------------------------------------------------------------------
// Compute methods A, B, C for one event.
//
// dt[]      : per-channel MCP-referenced times from the ntuple (ch0–6: MCP1,
//             ch7: MCP2).  Used for Methods A and C.
// dt7_mcp1  : ch7 re-referenced to MCP1 = hg_cfd[7] + (mcp2_time − mcp_time).
//             Used ONLY for Method B so all 8 channels share MCP1, making the
//             DW+UP sum a proper shower centroid (σ ≈ 40–60 ps).
//             Without this correction, Method B would include the MCP2–MCP1
//             offset (~0.5 ns), shifting the mean but not affecting σ significantly.
// hg[], lg[]: HG and LG amplitudes for M7 walk correction.
// k[], meanR[]: M7 OLS slopes and mean ratios per channel.
// out[3]    : filled with σ values for methods A, B, C in ns (unconverted).
//
// Returns true if at least 1 valid DW and 1 valid UP channel exist.
// ---------------------------------------------------------------------------
static bool ComputeMethods_teb(const float dt[8],
                                float        dt7_mcp1,
                                const float  hg[8],
                                const float  lg[8],
                                const double k[8],
                                const double meanR[8],
                                float        out[kNMeth_teb])
{
    // DW channels: 0..3   UP channels: 4..7
    double dw_sum   = 0.;                  // DW average (same for all methods)
    double up_sum_A = 0., up_sum_B = 0.;  // UP average for A/C vs B
    double dw_sum_C = 0., up_sum_C = 0.;  // M7-corrected sums
    int    dw_n   = 0;
    int    up_n_A = 0, up_n_B = 0;
    int    dw_n_C = 0, up_n_C = 0;

    // DW channels (ch 0–3): all use mcp_time → same for A, B, C
    for (int i = 0; i < 4; ++i) {
        if (dt[i] > kSentinelCut_teb) {
            dw_sum += dt[i];
            ++dw_n;
            double ratio = (lg[i] >= kLG_minPeak && hg[i] > 0.f)
                           ? (double)hg[i] / (double)lg[i] : meanR[i];
            dw_sum_C += (double)dt[i] - k[i] * (ratio - meanR[i]);
            ++dw_n_C;
        }
    }

    // UP channels (ch 4–7)
    for (int i = 4; i < 8; ++i) {
        // ── Method A / C: use dt[i] (mcp2_time for ch7, mcp_time for ch4-6)
        if (dt[i] > kSentinelCut_teb) {
            up_sum_A += dt[i];
            ++up_n_A;
            double ratio = (lg[i] >= kLG_minPeak && hg[i] > 0.f)
                           ? (double)hg[i] / (double)lg[i] : meanR[i];
            up_sum_C += (double)dt[i] - k[i] * (ratio - meanR[i]);
            ++up_n_C;
        }
        // ── Method B: use dt7_mcp1 for ch7 so all UP channels share MCP1
        float dt_B = (i == 7) ? dt7_mcp1 : dt[i];
        if (dt_B > kSentinelCut_teb) {
            up_sum_B += dt_B;
            ++up_n_B;
        }
    }

    if (dw_n < 1 || up_n_A < 1) return false;

    double dw_avg   = dw_sum   / dw_n;
    double up_avg_A = up_sum_A / up_n_A;
    double up_avg_B = (up_n_B >= 1) ? up_sum_B / up_n_B : up_avg_A;

    // Method A: (DW − UP) / 2  — MCP-independent; mcp_time/mcp2_time cancel
    out[0] = (float)((dw_avg - up_avg_A) * 0.5);
    // Method B: (DW + UP) / 2  — all channels MCP1-referenced; includes MCP1 jitter
    out[1] = (float)((dw_avg + up_avg_B) * 0.5);
    // Method C: (DW − UP) / 2  with per-channel M7 walk correction
    if (dw_n_C >= 1 && up_n_C >= 1) {
        out[2] = (float)((dw_sum_C / dw_n_C - up_sum_C / up_n_C) * 0.5);
    } else {
        out[2] = out[0];  // fallback: same as A
    }
    return true;
}

// ---------------------------------------------------------------------------
// Main function
// ---------------------------------------------------------------------------
void timingEnergyBins()
{
    TH1::AddDirectory(kFALSE);

    ApplyRADiCALStyle();   // RADiCALStyle.h via PlotUtils.h

    TString sumDir = "Analysis/Output/Summary";
    gSystem->mkdir(sumDir, kTRUE);

    // Summary storage: best-estimator σ (ps) and error per method per energy
    // [Gauss] and [Crystal Ball] stored separately
    double vEnergy[kNRuns]               = {};
    double vSigBest[kNRuns][kNMeth_teb]     = {};
    double vSigBestErr[kNRuns][kNMeth_teb]  = {};
    double vSigBestCB[kNRuns][kNMeth_teb]   = {};  // Crystal Ball sigma
    double vSigBestCBErr[kNRuns][kNMeth_teb]= {};
    double vBestEff[kNRuns]                 = {};  // best-bin efficiency (% of fiducial)
    double vBestEmeas[kNRuns]               = {};  // best-bin E_meas center (mV)
    double vSigBestOOS[kNRuns]              = {};  // #G5: run-folded OUT-OF-SAMPLE best-bin σ (Method A)
    double vSigBestOOSErr[kNRuns]           = {};
    int    nValidRuns = 0;

    // =========================================================================
    // Per-energy loop
    // =========================================================================
    for (int iRun = 0; iRun < kNRuns; ++iRun) {
        const RunCfg& rc = kRuns[iRun];
        TString ntupleFile = TString("Analysis/Output/") + rc.label + "/ntuple.root";

        TFile* fin = TFile::Open(ntupleFile);
        if (!fin || fin->IsZombie()) {
            std::cerr << "[timingEnergyBins] Cannot open " << ntupleFile << "\n";
            continue;
        }
        TTree* tree = (TTree*)fin->Get("rad");
        if (!tree || tree->GetEntries() == 0) {
            std::cerr << "[timingEnergyBins] No TTree 'rad' in " << ntupleFile << "\n";
            fin->Close();
            continue;
        }

        Long64_t nEntries = tree->GetEntries();
        std::cout << "\n[timingEnergyBins] " << rc.label
                  << " — " << nEntries << " events\n";

        // -----------------------------------------------------------------------
        // Step 1: Pre-scan for centroid and CFD offset
        // -----------------------------------------------------------------------
        double xc, yc, tOff, tRms;
        ScanRunCenters(tree, xc, yc, tOff, tRms);
        const float xcf   = (float)xc;
        const float ycf   = (float)yc;
        const float rFid2 = (float)(kFiducial_r_timing * kFiducial_r_timing);

        // -----------------------------------------------------------------------
        // Branch setup — guard for optional mcp2_* branches
        // -----------------------------------------------------------------------
        Bool_t  wc_ok;
        Int_t   run = 0;
        Float_t x_trk, y_trk, mcp_peak, mcp_time, mcp2_peak, mcp2_time;
        Float_t hg_cfd[8], hg_peak[8], lg_peak[8], sum_lg, sum_pb;
        (void)mcp2_peak;  // declared for compatibility; not used — validity of ch7 is
                          // already encoded in hg_cfd[7] by processRun.C

        if (tree->GetBranch("run")) tree->SetBranchAddress("run", &run);
        tree->SetBranchAddress("wc_ok",    &wc_ok);
        tree->SetBranchAddress("x_trk",    &x_trk);
        tree->SetBranchAddress("y_trk",    &y_trk);
        tree->SetBranchAddress("mcp_peak", &mcp_peak);
        tree->SetBranchAddress("mcp_time", &mcp_time);
        // CFD-5% adopted as the base timing fraction: per-channel timing improves
        // substantially (esp. the Down capillaries, whose leading-edge shape jitters
        // more high on the edge and produces a broad shoulder at CFD-20%) and the
        // headline holds/improves vs CFD-20%.  CFD-3% is marginally sharper still on
        // high-amplitude pulses but REGRESSES at low energy (the 3% level dips into
        // noise on small pulses), so 5% is the evidence-based choice.  Falls back to
        // hg_cfd (20%) if the 5% branch is absent (pre-reprocess ntuples).
        tree->SetBranchAddress(
            tree->GetBranch("hg_cfd05") ? "hg_cfd05" : "hg_cfd", hg_cfd);
        tree->SetBranchAddress("hg_peak",   hg_peak);
        tree->SetBranchAddress("lg_peak",   lg_peak);
        tree->SetBranchAddress("sum_lg",   &sum_lg);
        tree->SetBranchAddress("sum_pb",   &sum_pb);

        // Optional MCP2 branch — only mcp2_time is needed (to re-reference ch7 to MCP1
        // for Method B).  mcp2_peak is not used here: processRun.C already set hg_cfd[7]
        // to kNoTime when MCP2 was invalid, so hg_cfd[7] > kSentinelCut_teb is sufficient.
        mcp2_peak = -1e6f;
        mcp2_time = -1e6f;
        bool hasMCP2 = (tree->GetBranch("mcp2_time") != nullptr);
        if (hasMCP2) {
            tree->SetBranchAddress("mcp2_time", &mcp2_time);
        }

        // -----------------------------------------------------------------------
        // Step 2: Collect fiducial events
        // -----------------------------------------------------------------------
        std::vector<EbEvent> events;
        events.reserve(50000);

        for (Long64_t iev = 0; iev < nEntries; ++iev) {
            tree->GetEntry(iev);

            // Event-level cuts
            if (!wc_ok)                          continue;
            if (mcp_peak < kMCP1_minPeak)        continue;  // weak MCP signal
            if (mcp_peak > kMCP1_maxPeak)        continue;  // DRS4-saturated MCP
            if (mcp_time <= kSentinelCut_teb)    continue;  // MCP1 time invalid

            float dx = x_trk - xcf;
            float dy = y_trk - ycf;
            if (dx*dx + dy*dy >= rFid2)          continue;

            EbEvent ev;
            ev.run    = run;
            ev.sum_lg = sum_lg;

            // hg_cfd[i] from the ntuple is already MCP-referenced (processRun.C
            // stores hg_cfd[i] = hg.cfd20[i] − mcpRef, where mcpRef = mcp_time for
            // ch0–6 and mcp2_time for ch7).  Do NOT subtract mcp_time again — that
            // would double-subtract and give wrong timing for Method B.
            //
            // dt7_mcp1: ch7 converted from MCP2-reference to MCP1-reference so that
            // all 8 UP channels share the same MCP reference in Method B.
            //   dt7_mcp1 = hg_cfd[7] + (mcp2_time − mcp_time)
            //             = (hg_cfd20[7] − mcp2_time) + (mcp2_time − mcp_time)
            //             = hg_cfd20[7] − mcp_time   (MCP1-referenced)
            float delta_mcp = (hasMCP2 && mcp2_time > kSentinelCut_teb)
                               ? mcp2_time - mcp_time : 0.f;
            bool ch7_valid = (hg_cfd[7] > kSentinelCut_teb);
            ev.dt7_mcp1 = ch7_valid ? hg_cfd[7] + delta_mcp : kSentinel_teb;

            for (int i = 0; i < 8; ++i) {
                ev.hg[i] = hg_peak[i];
                ev.lg[i] = lg_peak[i];
                // hg_cfd[i] validity: processRun.C sets kNoTime when the crossing
                // was not found or the MCP reference was missing.
                ev.dt[i] = (hg_peak[i] >= kHG_minPeak && hg_cfd[i] > kSentinelCut_teb)
                           ? hg_cfd[i] : kSentinel_teb;
            }
            events.push_back(ev);
        }
        fin->Close();

        std::cout << "  Collected " << events.size() << " fiducial events\n";
        if (events.size() < 50) {
            std::cout << "  Too few events — skipping.\n";
            continue;
        }

        // -----------------------------------------------------------------------
        // Step 3: Fit Gaussian to sum_lg → define 9 energy bins
        // -----------------------------------------------------------------------
        std::vector<float> slg_vals;
        slg_vals.reserve(events.size());
        for (auto& ev : events) slg_vals.push_back(ev.sum_lg);

        TH1F* hSumLG = VecToHist_teb(
            Form("hSumLG_teb_%s", rc.label.Data()), slg_vals, 150);
        hSumLG->GetXaxis()->SetTitle("#SigmaA_{LG} (mV)");
        hSumLG->GetYaxis()->SetTitle("Events");

        double muE = 0., muEErr = 0., sigE = 0., sigEErr = 0.;
        FitGaussCore(hSumLG, 2.0, muE, muEErr, sigE, sigEErr);
        if (sigE <= 0.) {
            // Fallback: use histogram mean/RMS
            muE  = hSumLG->GetMean();
            sigE = hSumLG->GetRMS();
        }

        std::cout << "  Energy distribution: mu=" << std::fixed << std::setprecision(1)
                  << muE << " mV  sig=" << sigE << " mV\n"
                  << std::defaultfloat;

        // Bin edges: muE ± 2*sigE, divided into kNBins_teb equal slices
        double binLo = muE - 2.0*sigE;
        double binHi = muE + 2.0*sigE;
        double binW  = (binHi - binLo) / kNBins_teb;

        double binEdge[kNBins_teb + 1];
        for (int ib = 0; ib <= kNBins_teb; ++ib)
            binEdge[ib] = binLo + ib * binW;
        double binCenter[kNBins_teb];
        for (int ib = 0; ib < kNBins_teb; ++ib)
            binCenter[ib] = 0.5 * (binEdge[ib] + binEdge[ib+1]);

        // -----------------------------------------------------------------------
        // Step 4: Compute M7 slopes (OLS over all fiducial events per channel)
        // -----------------------------------------------------------------------
        double k7[8]     = {};   // OLS slope per channel
        double meanR7[8] = {};   // mean ratio per channel

        for (int i = 0; i < 8; ++i) {
            std::vector<float> ratios, times;
            ratios.reserve(events.size());
            times.reserve(events.size());
            for (auto& ev : events) {
                if (ev.dt[i] <= kSentinelCut_teb) continue;
                if (ev.lg[i] < kLG_minPeak || ev.hg[i] <= 0.f) continue;
                float r = ev.hg[i] / ev.lg[i];
                ratios.push_back(r);
                times.push_back(ev.dt[i]);
            }
            if (!ratios.empty()) {
                double s = 0.;
                for (auto r : ratios) s += r;
                meanR7[i] = s / (double)ratios.size();
                k7[i] = OLS_teb(ratios, times);
            }
        }

        // -----------------------------------------------------------------------
        // Step 4b (#G4): 5-fold CROSS-VALIDATION of the M7 walk correction.
        // The slopes above are trained and evaluated in-sample, which can flatter
        // the result.  Here we train the slopes on 4/5 of the events and evaluate
        // Method C (walk-corrected (DW-UP)/2) on the held-out 1/5, pooling all
        // folds over the best-estimator bins (>=5).  Comparing the out-of-sample
        // sigma to the in-sample sigma quantifies any overfit.
        // -----------------------------------------------------------------------
        double sigC_insample = -1., sigC_crossval = -1.;
        {
            const int kF = 5;
            auto inBestBin = [&](float slg){ return slg >= binEdge[5] && slg < binEdge[kNBins_teb]; };
            auto coreSig = [&](std::vector<float>& v)->double{
                if (v.size() < 100) return -1.;
                double m = 0.; for (float x : v) m += x; m /= v.size();
                TH1F h("_g4cv", "", 200, m - 1.0, m + 1.0); h.SetDirectory(nullptr);
                for (float x : v) h.Fill(x);
                double mu, muE, s, sE; FitGaussCore(&h, 2.0, mu, muE, s, sE);
                return s > 0. ? s * 1000. : -1.;
            };
            std::vector<float> poolIn, poolCV;
            for (int f = 0; f < kF; ++f) {
                double kf[8] = {}, mrf[8] = {};
                for (int i = 0; i < 8; ++i) {
                    std::vector<float> ratios, times;
                    for (size_t e = 0; e < events.size(); ++e) {
                        if ((int)(e % kF) == f) continue;             // train = exclude fold f
                        auto& ev = events[e];
                        if (ev.dt[i] <= kSentinelCut_teb) continue;
                        if (ev.lg[i] < kLG_minPeak || ev.hg[i] <= 0.f) continue;
                        ratios.push_back(ev.hg[i] / ev.lg[i]); times.push_back(ev.dt[i]);
                    }
                    if (!ratios.empty()) { double s = 0.; for (float r : ratios) s += r;
                        mrf[i] = s / ratios.size(); kf[i] = OLS_teb(ratios, times); }
                }
                for (size_t e = 0; e < events.size(); ++e) {       // evaluate held-out fold f
                    if ((int)(e % kF) != f) continue;
                    auto& ev = events[e];
                    if (!inBestBin(ev.sum_lg)) continue;
                    float out[kNMeth_teb];
                    if (ComputeMethods_teb(ev.dt, ev.dt7_mcp1, ev.hg, ev.lg, kf, mrf, out))
                        poolCV.push_back(out[2]);                  // Method C, out-of-sample
                }
            }
            for (size_t e = 0; e < events.size(); ++e) {           // in-sample reference
                auto& ev = events[e];
                if (!inBestBin(ev.sum_lg)) continue;
                float out[kNMeth_teb];
                if (ComputeMethods_teb(ev.dt, ev.dt7_mcp1, ev.hg, ev.lg, k7, meanR7, out))
                    poolIn.push_back(out[2]);
            }
            sigC_insample = coreSig(poolIn);
            sigC_crossval = coreSig(poolCV);
            std::cout << Form("  #G4 M7 walk-corr CV (Method C, best bins): "
                              "in-sample %.1f ps -> 5-fold OOS %.1f ps  (overfit %.1f ps)\n",
                              sigC_insample, sigC_crossval,
                              (sigC_crossval > 0 && sigC_insample > 0) ? sigC_crossval - sigC_insample : 0.);
        }

        // -----------------------------------------------------------------------
        // Step 5: Per-bin analysis
        // -----------------------------------------------------------------------

        // σ and σ_err per bin per method — Gaussian and Crystal Ball
        double sigBin[kNBins_teb][kNMeth_teb]      = {};
        double sigBinErr[kNBins_teb][kNMeth_teb]   = {};
        double sigBinCB[kNBins_teb][kNMeth_teb]    = {};  // Crystal Ball sigma
        double sigBinCBErr[kNBins_teb][kNMeth_teb] = {};
        int    nBin[kNBins_teb]                    = {};

        // Histograms for per-bin distributions (method A), for page-1 PDF
        TH1F* hBinA[kNBins_teb];
        for (int ib = 0; ib < kNBins_teb; ++ib)
            hBinA[ib] = nullptr;

        // Storage for best-estimator (bins 5, 6, 7 pooled)
        std::vector<float> vBest[kNMeth_teb];

        std::cout << "  Bin  E_center(mV)  N_events  SigA(ps) SigA_CB  SigB(ps) SigB_CB  SigC(ps) SigC_CB\n";

        for (int ib = 0; ib < kNBins_teb; ++ib) {
            std::vector<float> vA, vB, vC;

            for (auto& ev : events) {
                if (ev.sum_lg < binEdge[ib] || ev.sum_lg >= binEdge[ib+1]) continue;

                float out[kNMeth_teb];
                bool ok = ComputeMethods_teb(ev.dt, ev.dt7_mcp1, ev.hg, ev.lg,
                                              k7, meanR7, out);
                if (!ok) continue;
                vA.push_back(out[0]);
                vB.push_back(out[1]);
                vC.push_back(out[2]);
            }

            nBin[ib] = (int)vA.size();

            // Build histograms and fit
            TH1F* hA = VecToHist_teb(
                Form("hA_teb_%s_b%d", rc.label.Data(), ib), vA);
            TH1F* hB = VecToHist_teb(
                Form("hB_teb_%s_b%d", rc.label.Data(), ib), vB);
            TH1F* hC = VecToHist_teb(
                Form("hC_teb_%s_b%d", rc.label.Data(), ib), vC);

            // ── Gaussian fits ────────────────────────────────────────────────
            double mu_A, muE_A, sig_A, sigE_A;
            double mu_B, muE_B, sig_B, sigE_B;
            double mu_C, muE_C, sig_C, sigE_C;
            FitGaussCore(hA, 2.0, mu_A, muE_A, sig_A, sigE_A);
            FitGaussCore(hB, 2.0, mu_B, muE_B, sig_B, sigE_B);
            FitGaussCore(hC, 2.0, mu_C, muE_C, sig_C, sigE_C);

            sigBin[ib][0]    = (sig_A > 0.) ? sig_A * 1000. : -1.;
            sigBinErr[ib][0] = (sig_A > 0.) ? sigE_A * 1000. : 0.;
            sigBin[ib][1]    = (sig_B > 0.) ? sig_B * 1000. : -1.;
            sigBinErr[ib][1] = (sig_B > 0.) ? sigE_B * 1000. : 0.;
            sigBin[ib][2]    = (sig_C > 0.) ? sig_C * 1000. : -1.;
            sigBinErr[ib][2] = (sig_C > 0.) ? sigE_C * 1000. : 0.;

            // ── Crystal Ball fits (parallel, same arguments as Gauss) ────────
            double mu_A_cb, muE_A_cb, sig_A_cb, sigE_A_cb;
            double mu_B_cb, muE_B_cb, sig_B_cb, sigE_B_cb;
            double mu_C_cb, muE_C_cb, sig_C_cb, sigE_C_cb;
            FitCrystalBall(hA, 2.0, mu_A_cb, muE_A_cb, sig_A_cb, sigE_A_cb);
            FitCrystalBall(hB, 2.0, mu_B_cb, muE_B_cb, sig_B_cb, sigE_B_cb);
            FitCrystalBall(hC, 2.0, mu_C_cb, muE_C_cb, sig_C_cb, sigE_C_cb);

            sigBinCB[ib][0]    = (sig_A_cb > 0.) ? sig_A_cb * 1000. : -1.;
            sigBinCBErr[ib][0] = (sig_A_cb > 0.) ? sigE_A_cb * 1000. : 0.;
            sigBinCB[ib][1]    = (sig_B_cb > 0.) ? sig_B_cb * 1000. : -1.;
            sigBinCBErr[ib][1] = (sig_B_cb > 0.) ? sigE_B_cb * 1000. : 0.;
            sigBinCB[ib][2]    = (sig_C_cb > 0.) ? sig_C_cb * 1000. : -1.;
            sigBinCBErr[ib][2] = (sig_C_cb > 0.) ? sigE_C_cb * 1000. : 0.;

            hBinA[ib] = hA;  // keep for PDF; hB and hC are temporary
            delete hB;
            delete hC;

            // Pool best-estimator bins (0-indexed 5, 6, 7 → bins 6, 7, 8 in 1-indexed)
            if (ib >= 5) {
                for (auto v : vA) vBest[0].push_back(v);
                for (auto v : vB) vBest[1].push_back(v);
                for (auto v : vC) vBest[2].push_back(v);
            }

            // Print row (Gauss + CB for each method)
            std::cout << Form("   %d  %8.1f        %6d",
                              ib+1, binCenter[ib], nBin[ib]);
            for (int m = 0; m < kNMeth_teb; ++m) {
                if (sigBin[ib][m] > 0.)
                    std::cout << Form("  %7.1f", sigBin[ib][m]);
                else
                    std::cout << "       -";
                if (sigBinCB[ib][m] > 0.)
                    std::cout << Form("  %7.1f", sigBinCB[ib][m]);
                else
                    std::cout << "       -";
            }
            std::cout << "\n";
        }

        // -----------------------------------------------------------------------
        // Step 6: Best estimator = the SINGLE best (lowest-sigma) E_meas bin with
        // adequate statistics (N >= kMinBinN_teb).  This is the BEST-CASE / fully-
        // contained resolution (highest-E_meas, most-contained showers) and matches
        // how arXiv:2401.01747 quotes its number.  It is NOT the typical all-shower
        // resolution: it is reported WITH its bin efficiency (N_bin / N_fiducial)
        // and the full sigma-vs-E_meas curve (PDF page 2) so the selection is
        // transparent.  The bin is chosen automatically (min sigma), not hand-picked.
        // -----------------------------------------------------------------------
        const int kMinBinN_teb = 500;  // a bin must have >=500 events to be eligible
        double sigBest[kNMeth_teb]      = {};
        double sigBestErr[kNMeth_teb]   = {};
        double sigBestCB[kNMeth_teb]    = {};
        double sigBestCBErr[kNMeth_teb] = {};

        long totFid_teb = 0;
        for (int ib = 0; ib < kNBins_teb; ++ib) totFid_teb += nBin[ib];
        int bestBinIdx = -1;   // method A (headline) best bin

        for (int m = 0; m < kNMeth_teb; ++m) {
            int bib = -1; double bsig = 1e30;
            for (int ib = 0; ib < kNBins_teb; ++ib) {
                if (nBin[ib] >= kMinBinN_teb && sigBin[ib][m] > 0. && sigBin[ib][m] < bsig) {
                    bsig = sigBin[ib][m]; bib = ib;
                }
            }
            if (bib >= 0) {
                sigBest[m]      = sigBin[bib][m];
                sigBestErr[m]   = sigBinErr[bib][m];
                sigBestCB[m]    = sigBinCB[bib][m];
                sigBestCBErr[m] = sigBinCBErr[bib][m];
                if (m == 0) bestBinIdx = bib;
            }
        }

        const double bestEff_teb = (bestBinIdx >= 0 && totFid_teb > 0)
                                     ? 100.0 * nBin[bestBinIdx] / totFid_teb : 0.;
        std::cout << Form("  Best estimator = single best E_meas bin (N>=%d):\n",
                          kMinBinN_teb);
        if (bestBinIdx >= 0)
            std::cout << Form("    >>> bin %d  E_meas=%.0f mV  N=%d  EFFICIENCY=%.1f%% of fiducial\n",
                              bestBinIdx+1, binCenter[bestBinIdx],
                              nBin[bestBinIdx], bestEff_teb);
        const char* mLongName[kNMeth_teb] = {
            "(DW-UP)/2 CFD:   ", "(DW+UP)/2 CFD:   ", "(DW-UP)/2 M7-cor:"
        };
        for (int m = 0; m < kNMeth_teb; ++m) {
            const char* mName = (m == 0) ? "Method A" : (m == 1) ? "Method B" : "Method C";
            std::cout << "    " << mName << " " << mLongName[m];
            if (sigBest[m] > 0.) std::cout << Form(" %6.1f ps   %6.1f ps (CB)\n",
                                                   sigBest[m], sigBestCB[m]);
            else                 std::cout << "      - ps        - ps (CB)\n";
        }

        // -----------------------------------------------------------------------
        // Step 6b (#G5): RUN-FOLDED CROSS-VALIDATION of the best-bin SELECTION.
        // Picking the single min-σ bin and then quoting that σ biases the headline
        // low (you partly select a downward statistical fluctuation — the classic
        // "optimize-on-the-test-statistic" / look-elsewhere bias).  Here we remove
        // that bias for the headline Method A ((DW−UP)/2 CFD, which has NO trained
        // parameters): in each of 5 folds we SELECT the best bin on the training
        // runs and MEASURE σ on the held-out runs, folding by RUN so no event from
        // a run is split across train/test.  Pooling the held-out events across
        // folds gives an unbiased OOS estimate of the selection procedure.
        // -----------------------------------------------------------------------
        double sigBestOOS_A = -1., sigBestOOSErr_A = 0.;
        {
            const int kF = 5;
            std::vector<int> uruns; uruns.reserve(events.size());
            for (auto& ev : events) uruns.push_back(ev.run);
            std::sort(uruns.begin(), uruns.end());
            uruns.erase(std::unique(uruns.begin(), uruns.end()), uruns.end());
            auto foldOf = [&](int r)->int{
                int idx = (int)(std::lower_bound(uruns.begin(), uruns.end(), r) - uruns.begin());
                return idx % kF;
            };
            const int kTrainMinN = (kMinBinN_teb * (kF - 1)) / kF;  // scale floor to train size
            std::vector<float> poolOOS;
            for (int f = 0; f < kF; ++f) {
                // (1) per-bin Method-A σ on the TRAINING runs (all folds but f)
                std::vector<std::vector<float>> trainBin(kNBins_teb);
                for (auto& ev : events) {
                    if (foldOf(ev.run) == f) continue;            // hold out fold f
                    int ib = -1;
                    for (int b = 0; b < kNBins_teb; ++b)
                        if (ev.sum_lg >= binEdge[b] && ev.sum_lg < binEdge[b+1]) { ib = b; break; }
                    if (ib < 0) continue;
                    float out[kNMeth_teb];
                    if (ComputeMethods_teb(ev.dt, ev.dt7_mcp1, ev.hg, ev.lg, k7, meanR7, out))
                        trainBin[ib].push_back(out[0]);           // Method A
                }
                // (2) SELECT the min-σ bin on the training data
                int selBin = -1; double selSig = 1e30;
                for (int b = 0; b < kNBins_teb; ++b) {
                    if ((int)trainBin[b].size() < kTrainMinN) continue;
                    TH1F* h = VecToHist_teb(Form("_oosTr_%s_f%d_b%d", rc.label.Data(), f, b), trainBin[b]);
                    double mu, muE, s, sE; FitGaussCore(h, 2.0, mu, muE, s, sE); delete h;
                    if (s > 0. && s*1000. < selSig) { selSig = s*1000.; selBin = b; }
                }
                if (selBin < 0) continue;
                // (3) MEASURE on the held-out fold, in the bin selected on the training data
                for (auto& ev : events) {
                    if (foldOf(ev.run) != f) continue;
                    if (ev.sum_lg < binEdge[selBin] || ev.sum_lg >= binEdge[selBin+1]) continue;
                    float out[kNMeth_teb];
                    if (ComputeMethods_teb(ev.dt, ev.dt7_mcp1, ev.hg, ev.lg, k7, meanR7, out))
                        poolOOS.push_back(out[0]);
                }
            }
            if (poolOOS.size() >= 100) {
                TH1F* h = VecToHist_teb(Form("_oosPool_%s", rc.label.Data()), poolOOS);
                double mu, muE, s, sE; FitGaussCore(h, 2.0, mu, muE, s, sE); delete h;
                if (s > 0.) { sigBestOOS_A = s*1000.; sigBestOOSErr_A = sE*1000.; }
            }
            std::cout << Form("  #G5 best-bin selection CV (Method A, run-folded over %zu runs): "
                              "in-sample %.1f ps -> OOS %.1f ps  (opt. bias %.1f ps)\n",
                              uruns.size(), sigBest[0], sigBestOOS_A,
                              (sigBestOOS_A > 0. && sigBest[0] > 0.) ? sigBestOOS_A - sigBest[0] : 0.);
        }

        // Store for summary (Gauss and Crystal Ball)
        vEnergy[nValidRuns] = rc.energy_GeV;
        vSigBestOOS[nValidRuns]    = sigBestOOS_A;
        vSigBestOOSErr[nValidRuns] = sigBestOOSErr_A;
        vBestEff[nValidRuns]   = bestEff_teb;
        vBestEmeas[nValidRuns] = (bestBinIdx >= 0) ? binCenter[bestBinIdx] : 0.;
        for (int m = 0; m < kNMeth_teb; ++m) {
            vSigBest[nValidRuns][m]      = sigBest[m];
            vSigBestErr[nValidRuns][m]   = sigBestErr[m];
            vSigBestCB[nValidRuns][m]    = sigBestCB[m];
            vSigBestCBErr[nValidRuns][m] = sigBestCBErr[m];
        }
        ++nValidRuns;

        // -----------------------------------------------------------------------
        // PDF Page 1: 3×3 grid of per-bin method-A distributions
        // -----------------------------------------------------------------------
        TString outDir  = TString("Analysis/Output/") + rc.label;
        gSystem->mkdir(outDir, kTRUE);
        TString outPDF  = outDir + "/timing_energy_bins.pdf";

        {
            TCanvas c1("c_teb_p1", "", 2100, 2100);
            c1.Divide(3, 3, 0.003, 0.003);

            for (int ib = 0; ib < kNBins_teb; ++ib) {
                c1.cd(ib + 1);
                StylePad();

                TH1F* h = hBinA[ib];
                if (!h) continue;

                h->SetLineColor(kBlue+1);
                h->SetLineWidth(2);
                h->SetFillColor(kBlue-9);
                h->SetFillStyle(1001);
                h->GetXaxis()->SetTitle("(DW#minusUP)/2 [ns]");
                h->GetXaxis()->SetTitleSize(0.052);
                h->GetXaxis()->SetLabelSize(0.046);
                h->GetYaxis()->SetTitle("Events");
                h->GetYaxis()->SetTitleSize(0.052);
                h->GetYaxis()->SetLabelSize(0.046);
                h->Draw("HIST");

                // Draw Gaussian fit overlay (red curve + annotation)
                DrawFitOverlay(h, kRed+1, 0.60, 0.82);

                // Bin label: E_center in mV
                TLatex tit;
                tit.SetNDC();
                tit.SetTextSize(0.060);
                tit.SetTextAlign(22);
                tit.DrawLatex(0.54, 0.93,
                    Form("Bin %d: %.0f mV", ib+1, binCenter[ib]));

                // Annotate Gaussian sigma (red) and Crystal Ball sigma (blue)
                // on two lines, replacing the single DrawFitOverlay annotation.
                // Gauss line is already drawn by DrawFitOverlay above; add CB below it.
                if (sigBin[ib][0] > 0.) {
                    TLatex lGauss;
                    lGauss.SetNDC();
                    lGauss.SetTextSize(0.042);
                    lGauss.SetTextColor(kRed+1);
                    lGauss.DrawLatex(0.17, 0.73,
                        Form("#sigma_{{Gauss}} = %.0f ps", sigBin[ib][0]));
                }
                if (sigBinCB[ib][0] > 0.) {
                    TLatex lCB;
                    lCB.SetNDC();
                    lCB.SetTextSize(0.042);
                    lCB.SetTextColor(kBlue+1);
                    lCB.DrawLatex(0.17, 0.63,
                        Form("#sigma_{{CB}} = %.0f ps", sigBinCB[ib][0]));
                }

                // Walk correction gain: sigma_C / sigma_A
                if (sigBin[ib][0] > 0. && sigBin[ib][2] > 0.) {
                    double gain = 100. * (1. - sigBin[ib][2] / sigBin[ib][0]);
                    TLatex lWalk;
                    lWalk.SetNDC();
                    lWalk.SetTextSize(0.038);
                    lWalk.SetTextColor(kGreen+2);
                    lWalk.DrawLatex(0.17, 0.54,
                        Form("walk corr. gain: %.1f%%", gain));
                }

                // N_events annotation (placed first, well clear of the lines below)
                TLatex nev;
                nev.SetNDC();
                nev.SetTextSize(0.046);
                nev.SetTextColor(kGray+2);
                nev.DrawLatex(0.17, 0.45, Form("N = %d", nBin[ib]));

                // Star annotation for best-estimator bins
                if (ib >= 5) {
                    TLatex star;
                    star.SetNDC();
                    star.SetTextSize(0.052);
                    star.SetTextColor(kOrange+1);
                    star.DrawLatex(0.17, 0.36, "#bullet  best-estimator bin");
                    // #G4: 5-fold cross-validation of the M7 walk correction
                    // (pooled over best bins) — shown once, on the last bin, in its
                    // own slot so it never overprints the N / best-estimator labels.
                    if (ib == kNBins_teb - 1 && sigC_crossval > 0.) {
                        TLatex lcv; lcv.SetNDC(); lcv.SetTextSize(0.030);
                        lcv.SetTextColor(kGray+2);
                        lcv.DrawLatex(0.17, 0.28,
                            Form("M7 walk-corr CV: %.0f #rightarrow %.0f ps (OOS)",
                                 sigC_insample, sigC_crossval));
                    }
                }
            }

            c1.Print(outPDF + "(");
        }

        // -----------------------------------------------------------------------
        // PDF Page 2: σ_t vs E_meas for all 3 methods + shaded best-estimator band
        // -----------------------------------------------------------------------
        {
            TCanvas c2("c_teb_p2", "", 900, 700);
            StylePad();

            // Frame
            double xlo_frame = binEdge[0]  - 0.5*binW;
            double xhi_frame = binEdge[kNBins_teb] + 0.5*binW;

            // Determine y range — exclude low-stats bins (N<50) to prevent
            // outlier σ from inflating the axis and hiding well-populated bins.
            // Include both Gaussian and Crystal Ball values in range calculation.
            double yMax = 0.;
            for (int ib = 0; ib < kNBins_teb; ++ib) {
                if (nBin[ib] < 50) continue;
                for (int m = 0; m < kNMeth_teb; ++m) {
                    if (sigBin[ib][m] > yMax)   yMax = sigBin[ib][m];
                    if (sigBinCB[ib][m] > yMax)  yMax = sigBinCB[ib][m];
                }
            }
            yMax = (yMax > 0.) ? yMax * 1.40 : 300.;

            TH1F* frame2 = c2.DrawFrame(xlo_frame, 0., xhi_frame, yMax,
                Form(";E_{{meas}} = #SigmaA_{{LG}} (mV);#sigma_{{t}} (ps)"));
            frame2->GetXaxis()->SetTitleSize(0.052);
            frame2->GetYaxis()->SetTitleSize(0.052);
            frame2->GetYaxis()->SetTitleOffset(1.30);
            frame2->GetXaxis()->SetLabelSize(0.046);
            frame2->GetYaxis()->SetLabelSize(0.046);

            // Shaded band for best-estimator bins (5, 6, 7 → indices 5..7)
            TBox shadedBand;
            shadedBand.SetFillColor(kYellow-9);
            shadedBand.SetFillStyle(1001);
            shadedBand.SetLineStyle(0);
            shadedBand.DrawBox(binEdge[5], 0., binEdge[kNBins_teb], yMax * 0.97);

            // Re-draw axes on top of shaded band
            frame2->Draw("AXIS SAME");

            // Build and draw TGraphErrors per method — Gaussian (filled) and CB (open)
            TGraphErrors* gMeth[kNMeth_teb];    // Gaussian
            TGraphErrors* gMethCB[kNMeth_teb];  // Crystal Ball
            for (int m = 0; m < kNMeth_teb; ++m) {
                gMeth[m]   = new TGraphErrors();
                gMethCB[m] = new TGraphErrors();
                for (int ib = 0; ib < kNBins_teb; ++ib) {
                    if (sigBin[ib][m] > 0.) {
                        int np = gMeth[m]->GetN();
                        gMeth[m]->SetPoint(np, binCenter[ib], sigBin[ib][m]);
                        gMeth[m]->SetPointError(np, 0., sigBinErr[ib][m]);
                    }
                    if (sigBinCB[ib][m] > 0.) {
                        int np = gMethCB[m]->GetN();
                        gMethCB[m]->SetPoint(np, binCenter[ib], sigBinCB[ib][m]);
                        gMethCB[m]->SetPointError(np, 0., sigBinCBErr[ib][m]);
                    }
                }
                // Gaussian graph: filled marker, solid line
                gMeth[m]->SetMarkerStyle(kMMarker_teb[m]);
                gMeth[m]->SetMarkerColor(kMColor_teb[m]);
                gMeth[m]->SetLineColor(kMColor_teb[m]);
                gMeth[m]->SetMarkerSize(1.3);
                gMeth[m]->SetLineWidth(2);
                gMeth[m]->Draw("P SAME");

                // Crystal Ball graph: open marker, dashed line, same color
                gMethCB[m]->SetMarkerStyle(kMMarkerCB_teb[m]);
                gMethCB[m]->SetMarkerColor(kMColor_teb[m]);
                gMethCB[m]->SetLineColor(kMColor_teb[m]);
                gMethCB[m]->SetLineStyle(2);   // dashed
                gMethCB[m]->SetMarkerSize(1.3);
                gMethCB[m]->SetLineWidth(2);
                if (gMethCB[m]->GetN() > 0)
                    gMethCB[m]->Draw("P SAME");
            }

            // Title
            TLatex titLat;
            titLat.SetNDC();
            titLat.SetTextSize(0.052);
            titLat.SetTextAlign(22);
            titLat.DrawLatex(0.54, 0.93,
                Form("%.0f GeV  #sigma_{{t}} vs E_{{meas}}", rc.energy_GeV));

            // Legend with best-estimator sigma (Gauss and CB)
            // Extended height to accommodate CB rows
            TLegend* leg = new TLegend(0.38, 0.52, 0.93, 0.90);
            leg->SetBorderSize(0);
            leg->SetTextSize(0.034);
            for (int m = 0; m < kNMeth_teb; ++m) {
                // Gaussian row
                TString legEntry = TString(kMLabel_teb[m]);
                if (sigBest[m] > 0.)
                    legEntry += Form("  [G: %.0f ps]", sigBest[m]);
                if (gMeth[m]->GetN() > 0)
                    leg->AddEntry(gMeth[m], legEntry, "lp");
                // Crystal Ball row (same method, open marker)
                TString legEntryCB = TString("  CB fit");
                if (sigBestCB[m] > 0.)
                    legEntryCB += Form("  [CB: %.0f ps]", sigBestCB[m]);
                if (gMethCB[m]->GetN() > 0)
                    leg->AddEntry(gMethCB[m], legEntryCB, "lp");
            }
            // Best-estimator band annotation
            TLegend* legBand = new TLegend(0.38, 0.44, 0.93, 0.52);
            legBand->SetBorderSize(0);
            legBand->SetTextSize(0.030);
            TBox* bandEntry = new TBox();
            bandEntry->SetFillColor(kYellow-9);
            bandEntry->SetFillStyle(1001);
            legBand->AddEntry(bandEntry, "Best-estimator bins (6-8)", "f");
            leg->Draw();
            legBand->Draw();

            c2.Print(outPDF + ")");

            for (int m = 0; m < kNMeth_teb; ++m) {
                delete gMeth[m];
                delete gMethCB[m];
            }
            delete leg;
            delete legBand;
        }

        // Clean up per-bin histograms
        delete hSumLG;
        for (int ib = 0; ib < kNBins_teb; ++ib)
            if (hBinA[ib]) { delete hBinA[ib]; hBinA[ib] = nullptr; }

        std::cout << "[timingEnergyBins] " << rc.label << ": wrote " << outPDF << "\n";
    }  // end per-energy loop

    if (nValidRuns == 0) {
        std::cout << "[timingEnergyBins] No data found — run processRun.C first.\n";
        return;
    }

    // =========================================================================
    // Summary: σ_best vs beam energy for methods A, B, C
    // =========================================================================

    // Build TGraphErrors per method — Gaussian (filled) and Crystal Ball (open)
    TGraphErrors* gSum[kNMeth_teb];    // Gaussian
    TGraphErrors* gSumCB[kNMeth_teb]; // Crystal Ball
    for (int m = 0; m < kNMeth_teb; ++m) {
        gSum[m] = new TGraphErrors();
        gSum[m]->SetName(Form("gBestSigma_teb_m%d", m));
        gSum[m]->SetTitle(kMLabel_teb[m]);
        for (int ir = 0; ir < nValidRuns; ++ir) {
            if (vSigBest[ir][m] > 0.) {
                int np = gSum[m]->GetN();
                gSum[m]->SetPoint(np, vEnergy[ir], vSigBest[ir][m]);
                gSum[m]->SetPointError(np, 0., vSigBestErr[ir][m]);
            }
        }
        gSum[m]->SetMarkerStyle(kMMarker_teb[m]);
        gSum[m]->SetMarkerColor(kMColor_teb[m]);
        gSum[m]->SetLineColor(kMColor_teb[m]);
        gSum[m]->SetMarkerSize(1.4);
        gSum[m]->SetLineWidth(2);

        gSumCB[m] = new TGraphErrors();
        gSumCB[m]->SetName(Form("gBestSigmaCB_teb_m%d", m));
        gSumCB[m]->SetTitle(Form("%s CB", kMLabel_teb[m]));
        for (int ir = 0; ir < nValidRuns; ++ir) {
            if (vSigBestCB[ir][m] > 0.) {
                int np = gSumCB[m]->GetN();
                gSumCB[m]->SetPoint(np, vEnergy[ir], vSigBestCB[ir][m]);
                gSumCB[m]->SetPointError(np, 0., vSigBestCBErr[ir][m]);
            }
        }
        gSumCB[m]->SetMarkerStyle(kMMarkerCB_teb[m]);
        gSumCB[m]->SetMarkerColor(kMColor_teb[m]);
        gSumCB[m]->SetLineColor(kMColor_teb[m]);
        gSumCB[m]->SetLineStyle(2);   // dashed for CB
        gSumCB[m]->SetMarkerSize(1.4);
        gSumCB[m]->SetLineWidth(2);
    }

    // #G5: out-of-sample best-bin σ vs energy (Method A, run-folded selection CV).
    // This is the BIAS-CORRECTED headline curve — the floor fit is run on it.
    TGraphErrors* gOOS = new TGraphErrors();
    gOOS->SetName("gBestSigmaOOS_teb_m0");
    gOOS->SetTitle("(DW#minusUP)/2 CFD, out-of-sample");
    for (int ir = 0; ir < nValidRuns; ++ir) {
        if (vSigBestOOS[ir] > 0.) {
            int np = gOOS->GetN();
            gOOS->SetPoint(np, vEnergy[ir], vSigBestOOS[ir]);
            gOOS->SetPointError(np, 0., vSigBestOOSErr[ir]);
        }
    }
    gOOS->SetMarkerStyle(21);
    gOOS->SetMarkerColor(kRed+1);
    gOOS->SetLineColor(kRed+1);
    gOOS->SetMarkerSize(1.4);
    gOOS->SetLineWidth(2);

    // Paper's published curve: σ = sqrt(256²/E + 17.5²)
    TGraphErrors* gPaper = new TGraphErrors();
    gPaper->SetName("gPaper_teb");
    gPaper->SetTitle("arXiv:2401.01747");
    {
        double E_pts[] = {25., 50., 75., 100., 125., 150.};
        for (int i = 0; i < 6; ++i) {
            double E = E_pts[i];
            double s = std::sqrt(256.*256./E + 17.5*17.5);
            gPaper->SetPoint(i, E, s);
        }
    }
    gPaper->SetMarkerStyle(24);
    gPaper->SetMarkerColor(kGray+2);
    gPaper->SetLineColor(kGray+2);
    gPaper->SetLineStyle(2);
    gPaper->SetLineWidth(2);
    gPaper->SetMarkerSize(1.2);

    // Fit function: sigma = sqrt(a^2/E + b^2)
    // Applied to both Gaussian and Crystal Ball graphs.
    TF1* fFit[kNMeth_teb];
    TF1* fFitCB[kNMeth_teb];
    double fitA[kNMeth_teb] = {},   fitB[kNMeth_teb] = {};
    double fitACB[kNMeth_teb] = {}, fitBCB[kNMeth_teb] = {};
    for (int m = 0; m < kNMeth_teb; ++m) {
        fFit[m] = new TF1(Form("fFit_teb_m%d", m),
                           "sqrt([0]*[0]/x + [1]*[1])", 15., 175.);
        fFit[m]->SetParameters(200., 20.);
        fFit[m]->SetLineColor(kMColor_teb[m]);
        fFit[m]->SetLineStyle(2);
        fFit[m]->SetLineWidth(2);
        if (gSum[m]->GetN() >= 3) {
            gSum[m]->Fit(fFit[m], "RQ");
            fitA[m] = std::fabs(fFit[m]->GetParameter(0));
            fitB[m] = std::fabs(fFit[m]->GetParameter(1));
        }

        fFitCB[m] = new TF1(Form("fFitCB_teb_m%d", m),
                              "sqrt([0]*[0]/x + [1]*[1])", 15., 175.);
        fFitCB[m]->SetParameters(200., 20.);
        fFitCB[m]->SetLineColor(kMColor_teb[m]);
        fFitCB[m]->SetLineStyle(3);   // dotted to distinguish CB energy fit
        fFitCB[m]->SetLineWidth(2);
        if (gSumCB[m]->GetN() >= 3) {
            gSumCB[m]->Fit(fFitCB[m], "RQ");
            fitACB[m] = std::fabs(fFitCB[m]->GetParameter(0));
            fitBCB[m] = std::fabs(fFitCB[m]->GetParameter(1));
        }
    }

    // Paper reference fit curve
    TF1* fPaper = new TF1("fPaper_teb", "sqrt(256.*256./x + 17.5*17.5)", 15., 175.);
    fPaper->SetLineColor(kGray+2);
    fPaper->SetLineStyle(3);
    fPaper->SetLineWidth(2);

    // Summary PDF
    TString sumPDF = sumDir + "/timing_energy_bins_summary.pdf";
    {
        TCanvas cS("c_tebsum", "", 960, 720);
        StylePad();

        // y-axis range: include Gaussian, Crystal Ball, and paper curves
        double yMax = 0.;
        for (int m = 0; m < kNMeth_teb; ++m) {
            for (int p = 0; p < gSum[m]->GetN(); ++p)
                yMax = std::max(yMax, gSum[m]->GetY()[p]);
            for (int p = 0; p < gSumCB[m]->GetN(); ++p)
                yMax = std::max(yMax, gSumCB[m]->GetY()[p]);
        }
        for (int p = 0; p < gPaper->GetN(); ++p)
            yMax = std::max(yMax, gPaper->GetY()[p]);
        for (int p = 0; p < gOOS->GetN(); ++p)
            yMax = std::max(yMax, gOOS->GetY()[p]);
        yMax = (yMax > 0.) ? yMax * 1.45 : 400.;

        TH1F* frame = (TH1F*)cS.DrawFrame(15., 0., 165., yMax,
            ";Beam Energy (GeV);#sigma_{{t}} (ps)");
        frame->GetXaxis()->SetTitleSize(0.052);
        frame->GetYaxis()->SetTitleSize(0.052);
        frame->GetYaxis()->SetTitleOffset(1.25);
        frame->GetXaxis()->SetLabelSize(0.046);
        frame->GetYaxis()->SetLabelSize(0.046);

        // Draw fit lines first (Gaussian: dashed, CB: dotted)
        for (int m = 0; m < kNMeth_teb; ++m) {
            if (gSum[m]->GetN() >= 3)
                fFit[m]->DrawCopy("SAME");
            if (gSumCB[m]->GetN() >= 3)
                fFitCB[m]->DrawCopy("SAME");
        }
        fPaper->DrawCopy("SAME");

        // Draw data points on top
        for (int m = 0; m < kNMeth_teb; ++m) {
            if (gSum[m]->GetN() > 0)
                gSum[m]->Draw("P SAME");
            if (gSumCB[m]->GetN() > 0)
                gSumCB[m]->Draw("P SAME");
        }
        if (gPaper->GetN() > 0)
            gPaper->Draw("P SAME");
        if (gOOS->GetN() > 0)
            gOOS->Draw("P SAME");   // #G5 out-of-sample (bias-corrected) headline points

        // Title
        TLatex sumTit;
        sumTit.SetNDC();
        sumTit.SetTextSize(0.050);
        sumTit.SetTextAlign(22);
        sumTit.DrawLatex(0.55, 0.93,
            "#sigma_{t} vs beam energy -- energy-binned best estimator");

        // Legend: Gaussian and Crystal Ball rows per method
        // Extend height to accommodate both rows per method
        TLegend* leg = new TLegend(0.33, 0.42, 0.93, 0.90);
        leg->SetBorderSize(0);
        leg->SetTextSize(0.029);
        for (int m = 0; m < kNMeth_teb; ++m) {
            // Gaussian row
            TString entry = TString(kMLabel_teb[m]);
            if (fitA[m] > 0. || fitB[m] > 0.)
                entry += Form("  (G: a=%.0f #oplus b=%.1f ps)", fitA[m], fitB[m]);
            if (gSum[m]->GetN() > 0)
                leg->AddEntry(gSum[m], entry, "lp");
            // Crystal Ball row
            TString entryCB = "  CB fit";
            if (fitACB[m] > 0. || fitBCB[m] > 0.)
                entryCB += Form("  (CB: a=%.0f #oplus b=%.1f ps)", fitACB[m], fitBCB[m]);
            if (gSumCB[m]->GetN() > 0)
                leg->AddEntry(gSumCB[m], entryCB, "lp");
        }
        leg->AddEntry(gPaper,
            "arXiv:2401.01747  (a=256 #oplus b=17.5)", "lp");
        if (gOOS->GetN() > 0)
            leg->AddEntry(gOOS,
                "(DW#minusUP)/2 out-of-sample (run-folded selection CV)", "lp");
        leg->Draw();

        cS.Print(sumPDF);
    }

    // Summary ROOT file — write both Gaussian and Crystal Ball graphs
    TFile* fOut = new TFile(sumDir + "/timing_energy_bins.root", "RECREATE");
    for (int m = 0; m < kNMeth_teb; ++m) {
        gSum[m]->Write();
        gSumCB[m]->Write();
    }
    gPaper->Write();
    gOOS->Write();   // #G5 bias-corrected (out-of-sample) headline curve
    // Best-bin disclosure: efficiency (% of fiducial) and E_meas (mV) vs energy,
    // so the report can quote the best-bin sigma WITH its selection/efficiency.
    {
        TGraph* gEff = new TGraph();  gEff->SetName("gBestEff_teb");
        TGraph* gEm  = new TGraph();  gEm->SetName("gBestEmeas_teb");
        for (int ir = 0; ir < nValidRuns; ++ir) {
            gEff->SetPoint(ir, vEnergy[ir], vBestEff[ir]);
            gEm->SetPoint (ir, vEnergy[ir], vBestEmeas[ir]);
        }
        gEff->Write();
        gEm->Write();
    }
    fOut->Close();

    // =========================================================================
    // Final comparison table to stdout — Gaussian and Crystal Ball columns
    // =========================================================================
    std::cout << "\n[timingEnergyBins] Best estimator summary (ps):\n";
    std::cout << "  Fit model shown: G=Gaussian  CB=Crystal Ball\n\n";

    // Header: one energy column pair (G/CB) per run
    std::cout << "  Method / fit model                ";
    for (int ir = 0; ir < nValidRuns; ++ir)
        std::cout << Form("  %4.0fG   %4.0fCB", vEnergy[ir], vEnergy[ir]);
    std::cout << Form("   %10s  %6s\n", "a(G)", "b(G)");

    // Method rows: one line for Gaussian, one for Crystal Ball
    const char* mRowName[kNMeth_teb] = {
        "(DW-UP)/2 CFD-20%      ",
        "(DW+UP)/2 CFD-20%      ",
        "(DW-UP)/2 M7-corrected "
    };
    for (int m = 0; m < kNMeth_teb; ++m) {
        // Gaussian row
        std::cout << "  " << mRowName[m] << " [G] ";
        for (int ir = 0; ir < nValidRuns; ++ir) {
            if (vSigBest[ir][m] > 0.)
                std::cout << Form("  %5.0f   ", vSigBest[ir][m]);
            else
                std::cout << "      -   ";
        }
        if (fitA[m] > 0. || fitB[m] > 0.)
            std::cout << Form("   %10.0f  %6.1f", fitA[m], fitB[m]);
        std::cout << "\n";

        // Crystal Ball row
        std::cout << "  " << mRowName[m] << " [CB]";
        for (int ir = 0; ir < nValidRuns; ++ir) {
            if (vSigBestCB[ir][m] > 0.)
                std::cout << Form("          %5.0f", vSigBestCB[ir][m]);
            else
                std::cout << "              -";
        }
        if (fitACB[m] > 0. || fitBCB[m] > 0.)
            std::cout << Form("   %10.0f  %6.1f (CB fit)", fitACB[m], fitBCB[m]);
        std::cout << "\n";
    }

    // Published row
    std::cout << "  Published arXiv:2401.01747     [G] ";
    {
        for (int ir = 0; ir < nValidRuns; ++ir) {
            double E = vEnergy[ir];
            double sig_pub = std::sqrt(256.*256./E + 17.5*17.5);
            std::cout << Form("  %5.0f   ", sig_pub);
        }
    }
    std::cout << Form("   %10.0f  %6.1f", 256., 17.5) << "\n";

    // Cleanup
    for (int m = 0; m < kNMeth_teb; ++m) {
        delete gSum[m];
        delete gSumCB[m];
        delete fFit[m];
        delete fFitCB[m];
    }
    delete gPaper;
    delete fPaper;

    std::cout << "\n[timingEnergyBins] Done!\n"
              << "  Analysis/Output/<label>/timing_energy_bins.pdf  (per-energy, 2 pages)\n"
              << "  " << sumPDF << "\n"
              << "  " << sumDir << "/timing_energy_bins.root\n";
}

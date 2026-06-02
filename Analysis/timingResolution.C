// ============================================================================
// timingResolution.C — combined timing resolution analysis for RADiCAL
// ============================================================================
//
// Reads the compact ntuples produced by processRun.C and computes the timing
// resolution for five combination strategies, exploiting the fact that each of
// the four calorimeter corners has two readout capillaries: one "Down" (shower-
// max WLS, DRS0 channels 0-6) and one "Up" (DRS0 channels 4-7 + DRS0 G1 ch0).
//
// Corner layout  (corner index k = 0..3):
//   k=0 NW:  Down = hg_cfd[0] (NW-D),  Up = hg_cfd[4] (NW-U)
//   k=1 NE:  Down = hg_cfd[1] (NE-D),  Up = hg_cfd[5] (NE-U)
//   k=2 SE:  Down = hg_cfd[2] (SE-D),  Up = hg_cfd[6] (SE-U)
//   k=3 SW:  Down = hg_cfd[3] (SW-D),  Up = hg_cfd[7] (SW-U, uses MCP2)
//
// Note: hg_cfd[7] is referenced to MCP2 while hg_cfd[0..6] are referenced to
// MCP1.  A constant MCP1-MCP2 cable offset (δ_MCP ≈ few hundred ps) may shift
// the SW-corner mean slightly relative to the other corners.  This shifts the
// mean of any combined quantity but does NOT broaden the per-event Gaussian σ,
// so the measured timing resolution is unaffected.
//
// ── Combination methods ─────────────────────────────────────────────────────
//
//   M0  best1       Best individual capillary (channel with smallest σ_t)
//   M1  corner4     Mean of valid per-corner averages: t_k = (t_Dk + t_Uk) / 2
//                   → cancels longitudinal (shower-depth) jitter
//   M2  avg8        Arithmetic mean of all N_valid (0–8) channels
//   M3  wgt8        Amplitude²-weighted mean of all N_valid channels
//   M4  wgtCorner   Amplitude²-weighted mean of valid corner pairs
//
// ── Expected improvement ────────────────────────────────────────────────────
//
//   Combining N independent channels of equal σ_single gives σ = σ_single/√N.
//   For 8 independent channels: σ_avg8 ≈ σ_single/√8 ≈ σ_single/2.83.
//   For 4 corner pairs (each pair halves the single σ): σ_corner4 ≈ σ_single/√8.
//   Amplitude weighting exploits SNR differences between channels and can
//   approach (but not exceed) the theoretical √N limit.
//
// ── Quality cuts ─────────────────────────────────────────────────────────────
//
//   wc_ok == true          (wire-chamber track reconstructed)
//   mcp_peak > 100 mV      (clean MCP reference signal)
//   in_fiducial cut        (dynamic, centred on signal-weighted beam centroid)
//   hg_peak[i] > 20 mV     (channel-level: ensures well-above-noise pulse)
//   hg_cfd[i] > -1e5f      (channel-level: CFD crossing was found)
//
// ── Outputs ──────────────────────────────────────────────────────────────────
//
//   Per energy (Analysis/Output/<label>/):
//     combined_timing.pdf — method distributions + corner-pair distributions
//
//   Summary (Analysis/Output/Summary/):
//     timing_resolution_methods.pdf — σ_t vs E for all methods + √N guides
//     timing_summary.root           — TGraphErrors for offline inspection
//
// Usage (from repository root):
//   root -l -b -q 'Analysis/timingResolution.C+'
// ============================================================================

#include "ChannelConfig.h"
#include "PlotUtils.h"        // StylePad, FitGaussCore, DrawFitOverlay, ScanRunCenters
#include "DRS4Calibration.h"  // StopCellCorrection — DRS4 cell-width fixed-pattern removal

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLatex.h"
#include "TLine.h"
#include "TPad.h"
#include "TMath.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

// kHG_minPeak — provided by SelectionCuts.h via ChannelConfig.h

// ---------------------------------------------------------------------------
// Number of corner pairs (each corner has a Down and an Up capillary)
// ---------------------------------------------------------------------------
static const int kNC = 4;   // NW, NE, SE, SW

// ---------------------------------------------------------------------------
// Method indices
// ---------------------------------------------------------------------------
static const int kNM = 5;
static const char* kMName[kNM] = {
    "Best single channel",
    "4-corner (U+D)/2 mean",
    "Mean all 8 channels",
    "A^{2}-wgt all 8 channels",
    "A^{2}-wgt 4 corner pairs"
};
static const int kMColor[kNM] = {
    kBlue+1, kRed+1, kGreen+2, kMagenta+1, kOrange+2
};
static const int kMMarker[kNM] = { 20, 21, 22, 23, 29 };

static const int kCornerColor[kNC] = { kBlue+1, kRed+1, kGreen+2, kMagenta+1 };
static const char* kCornerName[kNC] = { "NW", "NE", "SE", "SW" };

// ScanRunCenters, FitGaussCore, DrawFitOverlay, StylePad — provided by PlotUtils.h.
// TR_FitGauss is an alias kept for readability at call sites.
#define TR_FitGauss FitGaussCore
// TR_ScanRunCenters and StylePadT are also aliased:
#define TR_ScanRunCenters ScanRunCenters
#define StylePadT StylePad

// ---------------------------------------------------------------------------
// Combination-timing helpers — the combo math is defined exactly ONCE here and
// reused by both the stop-cell training pre-pass and the main fill loop, so
// the two passes can never drift apart.
//
// Cell-width physics (verified in drs4TimeBase.C):
//   * M2 (mean-all) and M3 (A^2-weighted-all) carry a COMMON-MODE DRS0-G0
//     cell-width error that does not average down → stop-cell-correctable.
//   * M0 (best single) is noise-dominated → not corrected.
//   * M1 / M4 (corner (U+D) estimators) CANCEL the cell-width error in the
//     Down−Up structure → correcting them adds noise → not corrected.
// ---------------------------------------------------------------------------
static void TR_ComputeValid(const Float_t hg_peak[8], const Float_t hg_cfd[8],
                            bool hasMCP2, Float_t mcp2_peak, bool valid[8])
{
    for (int i = 0; i < 8; ++i) {
        const bool mcp_ok = kCap[i].use_mcp2
            ? (hasMCP2 && mcp2_peak >= kMCP2_minPeak && mcp2_peak <= kMCP2_maxPeak)
            : true;   // MCP1 already gated at event level
        valid[i] = mcp_ok && (hg_peak[i] > kHG_minPeak && hg_cfd[i] > -1e5f);
    }
}

// Mean of all valid channels (M2).  Returns false if no valid channel.
static bool TR_MeanAll(const Float_t hg_cfd[8], const bool valid[8], double& out)
{
    double sum = 0.; int nc = 0;
    for (int i = 0; i < 8; ++i) if (valid[i]) { sum += hg_cfd[i]; ++nc; }
    if (nc == 0) return false;
    out = sum / nc; return true;
}

// A^2-weighted mean of all valid channels (M3).  Returns false if no weight.
static bool TR_A2WeightedAll(const Float_t hg_cfd[8], const Float_t hg_peak[8],
                             const bool valid[8], double& out)
{
    double wsum = 0., wt = 0.;
    for (int i = 0; i < 8; ++i) if (valid[i]) {
        const double w2 = static_cast<double>(hg_peak[i]) * hg_peak[i];
        wsum += w2 * hg_cfd[i]; wt += w2;
    }
    if (wt <= 0.) return false;
    out = wsum / wt; return true;
}

// ===========================================================================
// Main
// ===========================================================================
void timingResolution()
{
    ApplyRADiCALStyle();   // RADiCALStyle.h via PlotUtils.h

    gSystem->mkdir("Analysis/Output/Summary", kTRUE);

    // Multi-energy storage  [method][energy_index] in ps
    std::vector<double> vE, vEErr;
    std::vector<double> vTRes[kNM],      vTResErr[kNM];
    std::vector<double> vTResCorner[kNC], vTResCornerErr[kNC];

    // ===========================================================================
    // Per-energy loop
    // ===========================================================================
    for (int iRun = 0; iRun < kNRuns; ++iRun)
    {
        const RunCfg& rc = kRuns[iRun];
        TString ntFile   = TString("Analysis/Output/") + rc.label + "/ntuple.root";

        TFile* fIn = TFile::Open(ntFile);
        if (!fIn || fIn->IsZombie()) {
            std::cout << "[timingResolution] Skipping " << rc.label
                      << " — ntuple not found (" << ntFile << ")\n";
            continue;
        }
        TTree* t = (TTree*)fIn->Get("rad");
        if (!t || t->GetEntries() == 0) { fIn->Close(); continue; }

        Long64_t nEntries = t->GetEntries();
        std::cout << "\n[timingResolution] " << rc.label
                  << " — " << nEntries << " events\n";

        // -----------------------------------------------------------------------
        // Pre-scan
        // -----------------------------------------------------------------------
        double x_center, y_center, t_cfd_offset, t_cfd_rms;
        TR_ScanRunCenters(t, x_center, y_center, t_cfd_offset, t_cfd_rms);

        // Histogram half-window [ns].  IMPORTANT: do NOT drive this from
        // t_cfd_rms — ScanRunCenters pools all 8 channels, so t_cfd_rms (~6.5 ns)
        // is dominated by the inter-channel cable-delay SPREAD (constant offsets),
        // not the per-event timing jitter.  Using 5×t_cfd_rms gave a ±32 ns range
        // → 161 ps/bin, far too coarse to resolve a ~60-100 ps combination core
        // (the Gaussian fit then biased the combo σ upward by ~20 ps).
        //
        // The true structure (cable spread ~1.1 ns + per-channel σ ≤ ~0.28 ns,
        // so ±5σ ≈ ±1.4 ns) fits comfortably in ±3 ns.  Fixed window → 6 ns/400
        // = 15 ps/bin, which resolves both the ~250 ps per-channel and the
        // ~60-95 ps combination cores.
        double t_win = 3.0;
        double t_lo  = t_cfd_offset - t_win;
        double t_hi  = t_cfd_offset + t_win;
        int    nBins = 400;

        // -----------------------------------------------------------------------
        // Book histograms
        // -----------------------------------------------------------------------

        // Individual channels (0-7)
        TH1F* hCh[8];
        for (int i = 0; i < 8; ++i)
            hCh[i] = new TH1F(Form("hCh_%s_%d", rc.label.Data(), i),
                               ";t_{HG}#minust_{MCP} (ns);Events",
                               nBins, t_lo, t_hi);

        // Corner Up+Down averages (require both channels valid)
        TH1F* hCorner[kNC];
        for (int k = 0; k < kNC; ++k)
            hCorner[k] = new TH1F(Form("hCorner_%s_%d", rc.label.Data(), k),
                                   ";t_{corner} (ns);Events",
                                   nBins, t_lo, t_hi);

        // Combined methods M1–M4
        // (M0 = best single channel, determined after fitting individual channels)
        TH1F* hM[kNM];
        hM[0] = nullptr;   // set after event loop once bestCh is known
        for (int m = 1; m < kNM; ++m)
            hM[m] = new TH1F(Form("hM_%s_%d", rc.label.Data(), m),
                              ";t_{HG}#minust_{MCP} (ns);Events",
                              nBins, t_lo, t_hi);

        // -----------------------------------------------------------------------
        // Branch addresses
        // -----------------------------------------------------------------------
        Bool_t  wc_ok;
        Float_t x_trk, y_trk, mcp_peak, mcp2_peak;
        Float_t hg_cfd[8], hg_peak[8];
        t->SetBranchAddress("wc_ok",    &wc_ok);
        t->SetBranchAddress("x_trk",    &x_trk);
        t->SetBranchAddress("y_trk",    &y_trk);
        t->SetBranchAddress("mcp_peak", &mcp_peak);
        // CFD-5% adopted as the base timing fraction (#6a; 3% regresses low-E, noise); falls back to
        // hg_cfd (20%) if the 3% branch is absent (pre-reprocess ntuples).
        t->SetBranchAddress(t->GetBranch("hg_cfd05") ? "hg_cfd05" : "hg_cfd", hg_cfd);
        t->SetBranchAddress("hg_peak",   hg_peak);

        // MCP2 amplitude — needed to validate the SW-U (channel 7) timing reference.
        // Default to pass-through so older ntuples without this branch still run.
        mcp2_peak = 9999.f;
        bool hasMCP2_tr = (t->GetBranch("mcp2_peak") != nullptr);
        if (hasMCP2_tr) t->SetBranchAddress("mcp2_peak", &mcp2_peak);

        // DRS4 stop cell (per group).  Present only after processRun.C has been
        // re-run with the stopcell branch; absent in older ntuples (correction
        // then degrades gracefully to a no-op — output is byte-identical).
        Int_t stopcell[4] = {0, 0, 0, 0};
        bool  hasStop = (t->GetBranch("stopcell") != nullptr);
        if (hasStop) t->SetBranchAddress("stopcell", stopcell);

        // -----------------------------------------------------------------------
        // Training pre-pass: build the stop-cell timing correction for the two
        // common-mode estimators (M2 mean-all, M3 A^2-weighted-all).  Keyed on
        // the DRS0 G0 stop cell (where channels 0-6 — the dominant weight — live).
        // Only runs when the stopcell branch exists.
        // -----------------------------------------------------------------------
        // 16 bins (not 64): the cell-width-vs-stop-cell pattern is smooth, and
        // fewer bins give ~500+ events/bin at 150 GeV for stable per-bin means.
        // Train on CORE events only (|combo − offset| < 1.5 ns): multi-ns combo
        // outliers from a mis-reconstructed channel would otherwise corrupt the
        // per-bin MEAN and inject huge spurious "pattern" (seen as ~200 ps).
        drs4::StopCellCorrection corrM2(16), corrM3(16);
        bool applyM2 = false, applyM3 = false;
        if (hasStop) {
            const double coreWin = 1.5;   // ns, around the pooled CFD offset
            for (Long64_t ev = 0; ev < nEntries; ++ev) {
                t->GetEntry(ev);
                if (!wc_ok || mcp_peak < kMCP1_minPeak || mcp_peak > kMCP1_maxPeak) continue;
                float dxt = x_trk - (float)x_center;
                float dyt = y_trk - (float)y_center;
                if (std::sqrt(dxt*dxt + dyt*dyt) >= (float)TimingFiducialR(rc.energy_GeV)) continue;

                bool vtrain[8];
                TR_ComputeValid(hg_peak, hg_cfd, hasMCP2_tr, mcp2_peak, vtrain);
                double m2, m3;
                if (TR_MeanAll(hg_cfd, vtrain, m2) &&
                    std::fabs(m2 - t_cfd_offset) < coreWin)
                    corrM2.Accumulate(stopcell[0], m2);
                if (TR_A2WeightedAll(hg_cfd, hg_peak, vtrain, m3) &&
                    std::fabs(m3 - t_cfd_offset) < coreWin)
                    corrM3.Accumulate(stopcell[0], m3);
            }
            corrM2.Finalize(/*minCount=*/50);
            corrM3.Finalize(/*minCount=*/50);
            // Self-protection: apply only if the learned pattern is physically
            // plausible.  The validated cell-width pattern is ~70 ps (drs4TimeBase);
            // a >150 ps "pattern" means training was noise/outlier-dominated, so
            // the correction is skipped rather than allowed to harm the result.
            applyM2 = (corrM2.OffsetRMS() < 0.150);
            applyM3 = (corrM3.OffsetRMS() < 0.150);
            std::cout << "  stop-cell correction: "
                      << Form("M2 pattern %.1f ps (%s), M3 pattern %.1f ps (%s)\n",
                              corrM2.OffsetRMS()*1000., applyM2 ? "applied" : "SKIPPED",
                              corrM3.OffsetRMS()*1000., applyM3 ? "applied" : "SKIPPED");
        }

        // -----------------------------------------------------------------------
        // Event loop
        // -----------------------------------------------------------------------
        long nFid = 0;
        for (Long64_t ev = 0; ev < nEntries; ++ev) {
            t->GetEntry(ev);

            // Global event cuts: WC track, MCP1 quality window (lower + upper)
            if (!wc_ok || mcp_peak < kMCP1_minPeak || mcp_peak > kMCP1_maxPeak) continue;
            float dx = x_trk - (float)x_center;
            float dy = y_trk - (float)y_center;
            if (std::sqrt(dx*dx + dy*dy) >= (float)TimingFiducialR(rc.energy_GeV)) continue;
            ++nFid;

            // Per-channel validity.
            // Channel 7 (SW-U) uses MCP2 as its timing reference — apply the MCP2
            // amplitude cut here to prevent weak-MCP2 events from contaminating SW-U.
            // All other channels are gated at the event level by kMCP1_minPeak above.
            bool valid[8];
            for (int i = 0; i < 8; ++i) {
                bool mcp_ok = kCap[i].use_mcp2
                              ? (hasMCP2_tr && mcp2_peak >= kMCP2_minPeak && mcp2_peak <= kMCP2_maxPeak)
                              : true;  // MCP1 already checked at event level
                valid[i] = mcp_ok && (hg_peak[i] > kHG_minPeak && hg_cfd[i] > -1e5f);
            }

            // Fill individual channel histograms
            for (int i = 0; i < 8; ++i)
                if (valid[i]) hCh[i]->Fill(hg_cfd[i]);

            // Per-corner Up+Down average (require both channels)
            float t_corner[kNC];
            float w_corner[kNC];  // A²_Down + A²_Up
            bool  c_ok[kNC];
            for (int k = 0; k < kNC; ++k) {
                bool dv = valid[k];
                bool uv = valid[k+4];
                if (dv && uv) {
                    t_corner[k] = 0.5f * (hg_cfd[k] + hg_cfd[k+4]);
                    w_corner[k] = hg_peak[k]*hg_peak[k] + hg_peak[k+4]*hg_peak[k+4];
                    c_ok[k]     = true;
                    hCorner[k]->Fill(t_corner[k]);
                } else {
                    t_corner[k] = 0.f;
                    w_corner[k] = 0.f;
                    c_ok[k]     = false;
                }
            }

            // M1: arithmetic mean of valid corner-pair averages
            {
                double sum = 0.; int nc = 0;
                for (int k = 0; k < kNC; ++k)
                    if (c_ok[k]) { sum += t_corner[k]; ++nc; }
                if (nc > 0) hM[1]->Fill(sum / nc);
            }

            // M2: arithmetic mean of all valid channels (stop-cell corrected).
            // corrM2.Offset() returns 0 when the correction was not trained
            // (no stopcell branch), so this is behaviour-preserving on old ntuples.
            {
                double m2;
                if (TR_MeanAll(hg_cfd, valid, m2))
                    hM[2]->Fill(m2 - (applyM2 ? corrM2.Offset(stopcell[0]) : 0.0));
            }

            // M3: amplitude²-weighted mean of all valid channels (stop-cell corrected)
            {
                double m3;
                if (TR_A2WeightedAll(hg_cfd, hg_peak, valid, m3))
                    hM[3]->Fill(m3 - (applyM3 ? corrM3.Offset(stopcell[0]) : 0.0));
            }

            // M4: amplitude²-weighted mean of valid corner-pair averages
            {
                double wsum = 0., wt = 0.;
                for (int k = 0; k < kNC; ++k)
                    if (c_ok[k]) { wsum += w_corner[k]*t_corner[k]; wt += w_corner[k]; }
                if (wt > 0.) hM[4]->Fill(wsum / wt);
            }
        } // event loop

        std::cout << "  Fiducial events: " << nFid
                  << Form(" / %lld (%.1f%%)\n", nEntries, 100.*nFid/nEntries);

        // -----------------------------------------------------------------------
        // Fit individual channels — find best single (M0)
        // -----------------------------------------------------------------------
        double sigCh[8] = {}, sigChErr[8] = {}, muCh[8] = {};
        int bestCh = 0;
        double bestSig = 1e9;
        for (int i = 0; i < 8; ++i) {
            double mu, muE;
            TR_FitGauss(hCh[i], 2.0, mu, muE, sigCh[i], sigChErr[i]);
            muCh[i] = mu;
            if (sigCh[i] > 0 && sigCh[i] < bestSig)
                { bestSig = sigCh[i]; bestCh = i; }
        }
        hM[0] = hCh[bestCh];   // alias; not double-deleted (hCh[] freed last)

        // -----------------------------------------------------------------------
        // Fit corner pairs
        // -----------------------------------------------------------------------
        double muCorn[kNC] = {}, sCorn[kNC] = {};
        for (int k = 0; k < kNC; ++k) {
            double mu, muE, s, sE;
            TR_FitGauss(hCorner[k], 2.0, mu, muE, s, sE);
            muCorn[k] = mu;  sCorn[k] = s;
            vTResCorner[k].push_back(s > 0 ? s*1000. : 0.);
            vTResCornerErr[k].push_back(s > 0 ? sE*1000. : 0.);
        }

        // -----------------------------------------------------------------------
        // Fit combined methods and store results
        // -----------------------------------------------------------------------
        vE.push_back(rc.energy_GeV);
        vEErr.push_back(0.);

        double muMeth[kNM] = {}, sMeth[kNM] = {};
        std::cout << "  Timing resolution summary @ " << rc.energy_GeV << " GeV:\n";
        for (int m = 0; m < kNM; ++m) {
            double mu, muE, s, sE;
            TR_FitGauss(hM[m], 2.0, mu, muE, s, sE);
            muMeth[m] = mu;  sMeth[m] = s;
            vTRes[m].push_back(s > 0 ? s*1000. : 0.);
            vTResErr[m].push_back(s > 0 ? sE*1000. : 0.);
            std::cout << "    " << std::left << std::setw(32) << kMName[m]
                      << ": " << std::fixed << std::setprecision(0)
                      << (s > 0 ? s*1000. : 0.) << " ps";
            if (m == 0) std::cout << " (" << kCap[bestCh].name << ")";
            std::cout << "\n";
        }
        std::cout << std::defaultfloat;

        // Also print per-corner
        for (int k = 0; k < kNC; ++k) {
            std::cout << "    " << kCornerName[k] << " corner (U+D)/2"
                      << std::setw(16) << ""
                      << ": " << std::fixed << std::setprecision(0)
                      << vTResCorner[k].back() << " ps\n";
        }
        std::cout << std::defaultfloat;

        // -----------------------------------------------------------------------
        // Per-energy plots
        // -----------------------------------------------------------------------
        TString outDir = TString("Analysis/Output/") + rc.label;
        gSystem->mkdir(outDir, kTRUE);

        // Open multi-page PDF for this energy
        TString pdfPath = outDir + "/combined_timing.pdf";

        // ── Page 1: 5 method distributions ─────────────────────────────────────
        {
            TCanvas c("c_meth", "", 1400, 900);
            c.Divide(3, 2, 0.01, 0.01);

            // Panels 1-5: methods M0-M4
            for (int m = 0; m < kNM; ++m) {
                c.cd(m+1);
                StylePadT();
                TH1F* h = hM[m];
                h->SetLineColor(kMColor[m]);
                h->SetLineWidth(2);
                h->GetXaxis()->SetTitle("t_{HG} #minus t_{MCP} (ns)");
                h->GetXaxis()->SetTitleSize(0.054);
                h->GetXaxis()->SetLabelSize(0.048);
                // Zoom x-axis to ±5.5σ around the fitted peak
                {
                    double muZ = (m == 0) ? muCh[bestCh] : muMeth[m];
                    double sZ  = (m == 0) ? sigCh[bestCh] : sMeth[m];
                    if (sZ > 0)
                        h->GetXaxis()->SetRangeUser(muZ - 5.5*sZ, muZ + 5.5*sZ);
                }
                h->Draw("HIST");
                DrawFitOverlay(h, kRed+1);

                // Method label at top of pad
                TLatex tit;
                tit.SetNDC();
                tit.SetTextSize(0.048);
                tit.SetTextColor(kMColor[m]);
                tit.DrawLatex(0.17, 0.90, kMName[m]);

                // For M0 note the best channel name, just below the method label
                if (m == 0) {
                    TLatex cname;
                    cname.SetNDC();
                    cname.SetTextSize(0.044);
                    cname.SetTextColor(kGray+2);
                    cname.DrawLatex(0.17, 0.82, Form("(%s)", kCap[bestCh].name));
                }
            }

            // Panel 6: text summary
            c.cd(6);
            gPad->SetLeftMargin(0.05);
            TLatex info;
            info.SetNDC();
            info.SetTextSize(0.075);
            info.DrawLatex(0.08, 0.90, Form("%.0f GeV", rc.energy_GeV));
            info.SetTextSize(0.050);
            info.DrawLatex(0.08, 0.82, "Timing resolution:");
            double dy = 0.;
            for (int m = 0; m < kNM; ++m) {
                double ps  = vTRes[m].back();
                double dps = vTResErr[m].back();
                info.SetTextColor(kMColor[m]);
                info.SetTextSize(0.038);
                info.DrawLatex(0.08, 0.72 - dy,
                    Form("%s:  %.0f#pm%.0f ps", kMName[m],
                         ps > 0 ? ps : 0., dps > 0 ? dps : 0.));
                dy += 0.115;
            }

            c.Print(pdfPath + "(");
        }

        // ── Page 2: 4 corner (U+D)/2 distributions in a 2×2 grid ─────────────
        {
            TCanvas c2("c_corn", "", 1000, 900);
            c2.Divide(2, 2, 0.01, 0.01);
            for (int k = 0; k < kNC; ++k) {
                c2.cd(k+1);
                StylePadT();
                hCorner[k]->SetLineColor(kCornerColor[k]);
                hCorner[k]->SetLineWidth(2);
                hCorner[k]->GetXaxis()->SetTitle("t_{corner} (ns)");
                hCorner[k]->GetXaxis()->SetTitleSize(0.054);
                // Zoom x-axis to ±5.5σ around the fitted peak
                if (sCorn[k] > 0)
                    hCorner[k]->GetXaxis()->SetRangeUser(muCorn[k] - 5.5*sCorn[k],
                                                          muCorn[k] + 5.5*sCorn[k]);
                hCorner[k]->Draw("HIST");
                DrawFitOverlay(hCorner[k], kRed+1);

                TLatex tit;
                tit.SetNDC();
                tit.SetTextSize(0.055);
                tit.SetTextColor(kCornerColor[k]);
                tit.DrawLatex(0.17, 0.90,
                    Form("%s corner  (U+D)/2", kCornerName[k]));
            }
            c2.cd(0);
            DrawPageTitle(Form("Corner (U+D)/2 timing distributions  --  %.0f GeV",
                               rc.energy_GeV));
            c2.Print(pdfPath + ")");
        }

        fIn->Close();
        // Note: hM[0] aliases hCh[bestCh]; delete hCh[] after hM[] to avoid
        // accessing freed memory.  (hM[1..4] are heap-allocated independently.)
        // ROOT owns histograms in the current directory — no manual delete needed
        // between runs when we close fIn; histogram memory is recycled.

    } // end per-energy loop

    if (vE.empty()) {
        std::cout << "[timingResolution] No data found — run processRun.C first.\n";
        return;
    }

    // ===========================================================================
    // Multi-energy summary
    // ===========================================================================
    int N = (int)vE.size();

    // Build TGraphErrors per method
    TGraphErrors* gM[kNM];
    for (int m = 0; m < kNM; ++m) {
        gM[m] = new TGraphErrors();
        gM[m]->SetName(Form("gTiming_M%d", m));
        gM[m]->SetTitle(kMName[m]);
        gM[m]->SetMarkerStyle(kMMarker[m]);
        gM[m]->SetMarkerSize(1.4);
        gM[m]->SetMarkerColor(kMColor[m]);
        gM[m]->SetLineColor(kMColor[m]);
        gM[m]->SetLineWidth(2);
        for (int i = 0; i < N; ++i)
            if ((int)vTRes[m].size() > i && vTRes[m][i] > 0)
                gM[m]->SetPoint(gM[m]->GetN(), vE[i], vTRes[m][i]);
    }

    // Build TGraphErrors per corner
    TGraphErrors* gCorner[kNC];
    for (int k = 0; k < kNC; ++k) {
        gCorner[k] = new TGraphErrors();
        gCorner[k]->SetName(Form("gCorner_%d", k));
        gCorner[k]->SetTitle(Form("%s corner (U+D)/2", kCornerName[k]));
        gCorner[k]->SetMarkerStyle(24+k);
        gCorner[k]->SetMarkerSize(1.3);
        gCorner[k]->SetMarkerColor(kCornerColor[k]);
        gCorner[k]->SetLineColor(kCornerColor[k]);
        gCorner[k]->SetLineStyle(2);
        gCorner[k]->SetLineWidth(2);
        for (int i = 0; i < N; ++i)
            if ((int)vTResCorner[k].size() > i && vTResCorner[k][i] > 0)
                gCorner[k]->SetPoint(gCorner[k]->GetN(), vE[i], vTResCorner[k][i]);
    }

    // Determine a representative single-channel σ at ~125 GeV for √N guide lines
    double sigRef = 0., eRef = 0.;
    for (int i = 0; i < gM[0]->GetN(); ++i) {
        double e = gM[0]->GetX()[i];
        if (std::fabs(e - 125.) < 40. && gM[0]->GetY()[i] > sigRef)
        { sigRef = gM[0]->GetY()[i]; eRef = e; }
    }
    if (sigRef <= 0 && gM[0]->GetN() > 0) {
        sigRef = gM[0]->GetY()[0];
        eRef   = gM[0]->GetX()[0];
    }

    // ── Summary page 1: all five methods ─────────────────────────────────────
    TString sumPDF = "Analysis/Output/Summary/timing_resolution_methods.pdf";

    {
        TCanvas* cS = new TCanvas("cTimMeth", "Timing resolution — methods", 960, 720);
        cS->SetLeftMargin(0.13);
        cS->SetBottomMargin(0.12);
        cS->SetRightMargin(0.05);
        cS->SetTopMargin(0.06);
        cS->SetTickx(1);
        cS->SetTicky(1);

        // Determine y range
        double yMax = 100.;
        for (int m = 0; m < kNM; ++m)
            for (int i = 0; i < gM[m]->GetN(); ++i)
                yMax = std::max(yMax, gM[m]->GetY()[i]);
        yMax *= 1.35;

        TH1F* frame = (TH1F*)cS->DrawFrame(0, 0, 180, yMax);
        frame->GetXaxis()->SetTitle("Beam Energy (GeV)");
        frame->GetYaxis()->SetTitle("#sigma_{t} (ps)");
        frame->GetXaxis()->SetTitleSize(0.048);
        frame->GetYaxis()->SetTitleSize(0.048);
        frame->GetYaxis()->SetTitleOffset(1.25);

        // √N guide lines at 125 GeV reference (N = 2, 4, 8 channels)
        if (sigRef > 0) {
            const int divNs[3] = {2, 4, 8};
            for (int iN = 0; iN < 3; ++iN) {
                int divN = divNs[iN];
                double sig_n = sigRef / std::sqrt((double)divN);
                TLine* ln = new TLine(20., sig_n, 170., sig_n);
                ln->SetLineColor(kGray+1);
                ln->SetLineStyle(3);
                ln->SetLineWidth(1);
                ln->Draw("SAME");
                TLatex ann;
                ann.SetTextSize(0.026);
                ann.SetTextColor(kGray+2);
                ann.DrawLatex(152., sig_n + yMax*0.010,
                    Form("#sigma_{0}/#sqrt{%d}=%.0fps", divN, sig_n));
            }
        }

        // ── Full parametric fit: σ_t = a/√E ⊕ b  for M1 (4-corner mean) ──────
        // This is the best combined estimator and the most comparable to the
        // paper's BestMinus result (arXiv:2401.01747, a=256, b=17.5 ps).
        TF1* fFull = nullptr;
        double fitA_tr = 0., fitB_tr = 0.;
        if (gM[1]->GetN() >= 3) {
            fFull = new TF1("fFullM1", "sqrt([0]*[0]/x + [1]*[1])", 15., 175.);
            fFull->SetParameters(200., 20.);
            fFull->SetParNames("a (#sqrt{GeV}#cdotp ps)", "b (ps)");
            gM[1]->Fit(fFull, "RQ");
            fitA_tr = std::fabs(fFull->GetParameter(0));
            fitB_tr = std::fabs(fFull->GetParameter(1));
            fFull->SetLineColor(kRed+1);
            fFull->SetLineStyle(2);
            fFull->SetLineWidth(2);
            fFull->DrawCopy("SAME");
        }

        // ── Paper's published curve (arXiv:2401.01747 BestMinus) ────────────
        // Their result: σ = 256/√E ⊕ 17.5 ps.  Drawn as reference.
        {
            TF1* fPub = new TF1("fPub_tr", "sqrt(256.*256./x + 17.5*17.5)", 15., 175.);
            fPub->SetLineColor(kGray+1);
            fPub->SetLineStyle(3);
            fPub->SetLineWidth(2);
            fPub->DrawCopy("SAME");
        }

        // ── Simple stochastic guide a/√E for M0 ─────────────────────────────
        if (sigRef > 0) {
            TF1* fSt = new TF1("fSt0", "[0]/sqrt(x)", 20., 160.);
            fSt->SetParameter(0, sigRef * std::sqrt(eRef));
            gM[0]->Fit(fSt, "RQ0");   // Q0: quiet, don't draw automatically
            fSt->SetLineColor(kBlue-7);
            fSt->SetLineStyle(2);
            fSt->SetLineWidth(1);
            fSt->DrawCopy("SAME");
        }

        // Draw method graphs (reverse order so M0 lands on top)
        for (int m = kNM-1; m >= 0; --m)
            if (gM[m]->GetN() > 0)
                gM[m]->Draw("P SAME");

        TLegend* leg = new TLegend(0.36, 0.60, 0.94, 0.93);
        leg->SetBorderSize(0);
        leg->SetTextSize(0.028);
        for (int m = 0; m < kNM; ++m)
            if (gM[m]->GetN() > 0)
                leg->AddEntry(gM[m], kMName[m], "lp");
        // Add fit annotation and paper reference
        if (fFull && fitA_tr > 0.) {
            TGraph* dummy1 = new TGraph(); dummy1->SetLineColor(kRed+1);
            dummy1->SetLineStyle(2); dummy1->SetLineWidth(2);
            leg->AddEntry(dummy1,
                Form("4-corner fit: %.0f/#sqrt{E}#oplus%.1f ps", fitA_tr, fitB_tr), "l");
        }
        {
            TGraph* dummy2 = new TGraph(); dummy2->SetLineColor(kGray+1);
            dummy2->SetLineStyle(3); dummy2->SetLineWidth(2);
            leg->AddEntry(dummy2, "arXiv:2401.01747: 256/#sqrt{E}#oplus17.5 ps", "l");
        }
        leg->Draw();

        cS->Print(sumPDF + "(");
    }

    // ── Summary page 2: per-corner breakdown ─────────────────────────────────
    {
        TCanvas* cC = new TCanvas("cTimCorners", "Timing resolution — corners", 960, 720);
        cC->SetLeftMargin(0.13);
        cC->SetBottomMargin(0.12);
        cC->SetRightMargin(0.05);
        cC->SetTopMargin(0.06);
        cC->SetTickx(1);
        cC->SetTicky(1);

        double yCMax = 100.;
        for (int k = 0; k < kNC; ++k)
            for (int i = 0; i < gCorner[k]->GetN(); ++i)
                yCMax = std::max(yCMax, gCorner[k]->GetY()[i]);
        yCMax *= 1.35;

        TH1F* frame2 = (TH1F*)cC->DrawFrame(0, 0, 180, yCMax);
        frame2->GetXaxis()->SetTitle("Beam Energy (GeV)");
        frame2->GetYaxis()->SetTitle("#sigma_{t} (ps)");
        frame2->GetXaxis()->SetTitleSize(0.048);
        frame2->GetYaxis()->SetTitleSize(0.048);
        frame2->GetYaxis()->SetTitleOffset(1.25);

        for (int k = 0; k < kNC; ++k)
            if (gCorner[k]->GetN() > 0)
                gCorner[k]->Draw("P SAME");

        TLegend* legC = new TLegend(0.55, 0.66, 0.94, 0.93);
        legC->SetBorderSize(0);
        legC->SetTextSize(0.032);
        for (int k = 0; k < kNC; ++k)
            if (gCorner[k]->GetN() > 0)
                legC->AddEntry(gCorner[k],
                    Form("%s corner (U+D)/2", kCornerName[k]), "lp");
        legC->Draw();

        // Overlay the 4-corner mean for comparison
        if (gM[1]->GetN() > 0) {
            TGraphErrors* gM1copy = (TGraphErrors*)gM[1]->Clone("_gM1cmp");
            gM1copy->SetMarkerStyle(34);
            gM1copy->SetMarkerSize(1.6);
            gM1copy->SetMarkerColor(kBlack);
            gM1copy->SetLineColor(kBlack);
            gM1copy->Draw("P SAME");
            legC->AddEntry(gM1copy, "4-corner mean (M1)", "lp");
        }

        cC->Print(sumPDF + ")");
    }

    // ===========================================================================
    // Write summary ROOT file
    // ===========================================================================
    TFile* fOut = new TFile("Analysis/Output/Summary/timing_summary.root", "RECREATE");
    for (int m = 0; m < kNM; ++m)  gM[m]->Write();
    for (int k = 0; k < kNC; ++k)  gCorner[k]->Write();
    fOut->Close();

    std::cout << "\n[timingResolution] Done!\n"
              << "  Analysis/Output/<label>/combined_timing.pdf"
                 "   (per-energy, 2 pages)\n"
              << "  Analysis/Output/Summary/timing_resolution_methods.pdf"
                 "  (2-page summary)\n"
              << "  Analysis/Output/Summary/timing_summary.root\n";
}

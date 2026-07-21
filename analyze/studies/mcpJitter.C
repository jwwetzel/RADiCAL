// ============================================================================
// mcpJitter.C — inter-group (mezzanine) timing-reference jitter
// ============================================================================
//
// HARDWARE FACT (verified from data: corr(MCP1,MCP2)=1.000, equal amplitude):
// MCP1 and MCP2 are NOT two independent MCPs — they are ONE MCP signal passively
// split into the two DT5742 readout groups (group 0 = ch 0-7, group 1 = ch 9-16,
// on separate mezzanines).  Each group digitises its own copy with its own
// free-running DRS4 domino wave / stop cell.  Channels are referenced to the MCP
// copy IN THEIR OWN GROUP (ch 0-6 → MCP1; SW-U, the only group-1 capillary →
// MCP2) so that same-group timing is free of the group-to-group jitter.
//
// WHAT σ(MCP1 - MCP2) ACTUALLY MEASURES:
//   t_MCPk = T + j_MCP + d_k + n_k   (k = group 0,1)
//     T      = beam arrival,  j_MCP = the MCP's own jitter (SAME in both copies),
//     d_k    = group-k DRS4 digitisation timing (stop cell / clock / cell width),
//     n_k    = CFD/baseline noise on that copy.
//   ⇒ MCP1 - MCP2 = (d_0 - d_1) + (n_0 - n_1)   —  j_MCP and T CANCEL.
// So σ(MCP1-MCP2) is the INTER-GROUP (inter-mezzanine) relative timing jitter,
// NOT the MCP's intrinsic resolution (which is common-mode and cancels here).
// The cable-length difference between the two copies is a fixed offset (mean of
// the difference), not jitter.
//
//   σ(MCP1-MCP2) = √(σ_d0² + σ_d1² + ...)  ≈  √2 × σ_pergroup
//   σ_pergroup   = σ(MCP1-MCP2)/√2  ≈ 71 ps   (per-group reference jitter)
//
// CONSEQUENCES:
//   * The corner estimators kept WITHIN group 0 (NW, NE, SE: both D and U on
//     MCP1) are EXACTLY free of both the MCP jitter and the inter-group jitter —
//     which is why the headline (DW−UP)/2 sits below the 71 ps per-group floor.
//     Only SW-U sits in group 1 (MCP2), so the SW corner carries a sub-dominant
//     inter-group residual; the 3 group-0 corners do not.
//   * CAVEAT: the σ_crystal = √(σ_meas² − σ_pergroup²) subtraction below is only
//     valid for estimators that actually MIX the two groups (e.g. mean-all-8 with
//     SW-U on MCP2).  For a SAME-GROUP single channel (t_ch − t_MCP1) the group-0
//     digitisation d_0 cancels, so the inter-group term is NOT present and must
//     not be subtracted.  Treat the per-channel "crystal" graphs with this caveat.
//
// ── What this tells us ───────────────────────────────────────────────────────
//
//   If σ_MCP_single is small compared to σ_measured, the correction is minor
//   and the MCP is not limiting our resolution.
//
//   If σ_MCP_single is a significant fraction of σ_measured (> ~30%), then
//   upgrading to a better reference (e.g. a SiPM MCP) is a high-leverage
//   hardware change.
//
//   The per-energy σ_MCP_single also tells us whether MCP jitter scales with
//   energy (it shouldn't for a good reference) or is constant (electronic
//   noise dominated).
//
// ── Output ───────────────────────────────────────────────────────────────────
//
//   output/Summary/mcp_jitter.pdf        (5 pages)
//   output/Summary/mcp_jitter.root
//     → gMCP_jitter      : σ_MCP_single [ps] vs energy
//     → gTEB_corrected   : corrected energy-binned σ_t vs energy (method m0)
//     → gMCPvsAmp        : σ_MCP_single [ps] vs MCP1 amplitude (150 GeV only)
//     → mcp_jitter_p0    : TParameter<Double_t> — noise term [ps·mV]
//     → mcp_jitter_p1    : TParameter<Double_t> — electronic floor [ps]
//
//   Page 1:  MCP1–MCP2 time difference distributions (all 6 energies overlaid)
//   Page 2:  σ_MCP_single vs beam energy (with fit)
//   Page 2b: σ_MCP_single vs MCP1 amplitude at 150 GeV (amplitude-dependent fit)
//   Page 3:  Corrected vs uncorrected σ_t vs energy (energy-binned m0)
//   Page 3c: Combination methods — measured (MCP-limited) vs detector-only
//            intrinsic σ_t = √(σ_combo² − σ_MCP²); MCP is common-mode in all
//            combos.  → gIntrinsic_M0..M4 in mcp_jitter.root
//   Page 4:  Summary table — MCP jitter, σ_MCP(A_mean), measured σ_t, corrected σ_t, fraction
//
// Usage (from repository root):
//   root -l -b -q 'Analysis/mcpJitter.C+'
// ============================================================================

#include "ChannelConfig.h"
#include "PlotUtils.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TAxis.h"
#include "TF1.h"

#include "TParameter.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <utility>
#include <vector>

// Energy colors from RADiCALStyle.h (kREnergyCols, set by ApplyRADiCALStyle)
static const int kEMark[6] = { 20, 21, 22, 23, 24, 25 };

// ---------------------------------------------------------------------------
// OpenNtuple
// ---------------------------------------------------------------------------
static TFile* OpenNtuple(int r, TTree*& tree)
{
    TString path = Form("output/%s/ntuple.root",
                        kRuns[r].label.Data());
    TFile* f = TFile::Open(path);
    if (!f || f->IsZombie()) {
        std::cerr << "mcpJitter: cannot open " << path << "\n";
        tree = nullptr; return nullptr;
    }
    tree = static_cast<TTree*>(f->Get("rad"));
    if (!tree) {
        std::cerr << "mcpJitter: tree 'rad' missing in " << path << "\n";
        f->Close(); delete f; tree = nullptr; return nullptr;
    }
    return f;
}

// ---------------------------------------------------------------------------
// PageTitle
// ---------------------------------------------------------------------------
static void PageTitle(const char* t)
{
    TLatex lat; lat.SetNDC(); lat.SetTextSize(0.030); lat.SetTextAlign(21);
    lat.DrawLatex(0.50, 0.987, t);
}

// ===========================================================================
// Main
// ===========================================================================
void mcpJitter()
{
    TH1::AddDirectory(kFALSE);
    ApplyRADiCALStyle();   // RADiCALStyle.h via PlotUtils.h

    gSystem->mkdir("output/Summary", kTRUE);

    std::cout << "mcpJitter: measuring MCP1-MCP2 timing jitter\n";

    // Results to fill
    double sigMCP[6] = {};     // σ_MCP_single [ps]
    double sigMCPErr[6] = {};

    // Per-energy MCP1-MCP2 timing difference histograms
    TH1F* hDiff[6] = {};

    // Amplitude-binned analysis (150 GeV only)
    // Pairs of (A_MCP1 [mV], delta_t [ns])
    std::vector<std::pair<float,float>> ampDiffPairs;
    ampDiffPairs.reserve(20000);

    // =========================================================================
    // Event loop — fill MCP1-MCP2 time difference per energy
    // =========================================================================
    for (int r = 0; r < kNRuns; ++r) {
        TTree* t = nullptr;
        TFile* fin = OpenNtuple(r, t);
        if (!fin) continue;

        // Pre-scan to get beam centroid (needed for fiducial cut)
        double xc, yc, tcfd, trms;
        ScanRunCenters(t, xc, yc, tcfd, trms);

        Float_t x_trk, y_trk, mcp_peak, mcp_time, mcp2_peak, mcp2_time;
        Bool_t  wc_ok;
        t->SetBranchAddress("x_trk",     &x_trk);
        t->SetBranchAddress("y_trk",     &y_trk);
        t->SetBranchAddress("wc_ok",     &wc_ok);
        t->SetBranchAddress(t->GetBranch("mcp1_peak")?"mcp1_peak":"mcp_peak",  &mcp_peak);
        t->SetBranchAddress(t->GetBranch("mcp1_time")?"mcp1_time":"mcp_time",  &mcp_time);
        t->SetBranchAddress("mcp2_peak", &mcp2_peak);
        t->SetBranchAddress("mcp2_time", &mcp2_time);

        // First pass: collect Δt = mcp_time - mcp2_time for good events
        std::vector<float> diffs;
        diffs.reserve(10000);

        Long64_t nEv = t->GetEntries();
        for (Long64_t ev = 0; ev < nEv; ++ev) {
            t->GetEntry(ev);
            if (!wc_ok) continue;

            // Both MCPs must have valid, non-saturated signals
            if (mcp_peak  < kMCP1_minPeak || mcp_peak  > kMCP1_maxPeak) continue;
            if (mcp2_peak < kMCP2_minPeak || mcp2_peak > kMCP2_maxPeak) continue;

            // Timing fiducial
            float dx = x_trk - static_cast<float>(xc);
            float dy = y_trk - static_cast<float>(yc);
            if (std::sqrt(dx*dx + dy*dy) >= static_cast<float>(TimingFiducialR(kRuns[r].energy_GeV)))
                continue;

            // Sentinel check: both times must be valid
            if (mcp_time  < -1e5f || mcp2_time  < -1e5f) continue;

            float dt = mcp_time - mcp2_time;
            diffs.push_back(dt);

            // Collect amplitude-delta_t pairs for 150 GeV amplitude analysis
            if (std::fabs(kRuns[r].energy_GeV - 150.) < 1.)
                ampDiffPairs.emplace_back(mcp_peak, dt);
        }

        t->ResetBranchAddresses();
        fin->Close();
        delete fin;

        if (diffs.size() < 50) {
            std::cout << "  " << kRuns[r].label
                      << ": too few events (" << diffs.size() << ") — skipping\n";
            continue;
        }

        // Build histogram auto-ranged on mean ± 4 RMS
        double mu = 0., rms2 = 0.;
        for (float v : diffs) mu += v;
        mu /= diffs.size();
        for (float v : diffs) rms2 += (v-mu)*(v-mu);
        double rms = std::sqrt(rms2 / diffs.size());
        if (rms < 0.010) rms = 0.100;

        hDiff[r] = new TH1F(Form("hDiff_%s", kRuns[r].label.Data()), "",
                             200, mu - 4.*rms, mu + 4.*rms);
        hDiff[r]->SetDirectory(nullptr);
        for (float v : diffs) hDiff[r]->Fill(v);

        // Gaussian core fit
        double fitMu, fitMuErr, fitSig, fitSigErr;
        FitGaussCore(hDiff[r], 2.0, fitMu, fitMuErr, fitSig, fitSigErr);

        if (fitSig > 0.) {
            // σ(MCP1 - MCP2) = fitSig [ns]
            // σ_MCP_single    = fitSig / √2
            sigMCP[r]    = fitSig * 1000. / std::sqrt(2.0);   // ps
            sigMCPErr[r] = fitSigErr * 1000. / std::sqrt(2.0);
        }
        std::cout << "  " << kRuns[r].label
                  << Form(": σ(MCP1-MCP2) = %.1f ps"
                          "  →  σ_MCP_single = %.1f ± %.1f ps\n",
                          fitSig * 1000., sigMCP[r], sigMCPErr[r]);
    }

    // =========================================================================
    // Amplitude-binned σ_MCP analysis (150 GeV only)
    // Bin (A_MCP1, delta_t) pairs into 10 bins over 200–750 mV
    // =========================================================================

    TGraphErrors* gMCPvsAmp = new TGraphErrors();
    TF1* fMCPvA = nullptr;
    double mcp_p0 = 0., mcp_p1 = 0.;   // fit parameters to store

    if (!ampDiffPairs.empty()) {
        const int    kNAmpBins   = 10;
        const double kAmpLo      = static_cast<double>(kMCP1_minPeak);   // 200 mV
        const double kAmpHi      = static_cast<double>(kMCP1_maxPeak);   // 750 mV
        const double kBinWidth   = (kAmpHi - kAmpLo) / kNAmpBins;

        std::vector<std::vector<float>> binDt(kNAmpBins);
        std::vector<double>             binAmp(kNAmpBins, 0.);
        std::vector<int>                binCount(kNAmpBins, 0);

        for (auto& pr : ampDiffPairs) {
            int ib = static_cast<int>((pr.first - kAmpLo) / kBinWidth);
            if (ib < 0 || ib >= kNAmpBins) continue;
            binDt[ib].push_back(pr.second);
            binAmp[ib]  += static_cast<double>(pr.first);
            binCount[ib]++;
        }

        double meanSigForSeed = 0.; int nSeedPts = 0;

        for (int ib = 0; ib < kNAmpBins; ib++) {
            if (binDt[ib].size() < 20) continue;

            double aCenter = (binCount[ib] > 0)
                             ? binAmp[ib] / binCount[ib]
                             : kAmpLo + (ib + 0.5) * kBinWidth;

            // Compute RMS of delta_t in this bin (units: ns)
            double mu_b = 0.;
            for (float v : binDt[ib]) mu_b += v;
            mu_b /= binDt[ib].size();
            double rms2_b = 0.;
            for (float v : binDt[ib]) rms2_b += (v - mu_b) * (v - mu_b);
            double rms_b = std::sqrt(rms2_b / binDt[ib].size());

            // σ_MCP_single = rms / sqrt(2), convert ns -> ps
            double sigBin    = rms_b * 1000. / std::sqrt(2.0);
            double sigBinErr = sigBin / std::sqrt(2.0 * static_cast<double>(binDt[ib].size()) - 2.0);

            int np = gMCPvsAmp->GetN();
            gMCPvsAmp->SetPoint(np, aCenter, sigBin);
            gMCPvsAmp->SetPointError(np, kBinWidth * 0.5, sigBinErr);

            meanSigForSeed += sigBin; ++nSeedPts;
        }

        if (nSeedPts > 0) meanSigForSeed /= nSeedPts;

        // Fit: f(A) = sqrt(p0^2/A^2 + p1^2)
        //   p0 ~ mean_sigma * 400  (noise term, units ps*mV)
        //   p1 ~ 50 ps             (electronic floor)
        fMCPvA = new TF1("fMCPvA", "sqrt([0]*[0]/(x*x) + [1]*[1])",
                         kAmpLo, kAmpHi);
        fMCPvA->SetParameters(meanSigForSeed * 400., 50.);
        fMCPvA->SetParNames("p0_noise", "p1_floor");

        if (gMCPvsAmp->GetN() >= 3) {
            gMCPvsAmp->Fit(fMCPvA, "RQ0");   // quiet, no draw
            mcp_p0 = std::fabs(fMCPvA->GetParameter(0));
            mcp_p1 = std::fabs(fMCPvA->GetParameter(1));
            std::cout << Form("  Amplitude fit (150 GeV): p0 = %.1f ps*mV"
                              "  p1 = %.1f ps\n", mcp_p0, mcp_p1);
        }
    } else {
        std::cout << "  No 150 GeV amplitude pairs collected -- skipping amp analysis\n";
    }

    // =========================================================================
    // Apply MCP jitter correction to the per-channel average timing and to
    // the combination methods that still carry MCP jitter.
    //
    // IMPORTANT: the (DW−UP)/2 estimator (timing_energy_bins.root method m0)
    // is MCP-JITTER-FREE by construction:
    //
    //   t_DW = t_crystal_DW − t_MCP
    //   t_UP = t_crystal_UP − t_MCP
    //   (t_DW − t_UP)/2 → t_MCP drops out exactly
    //
    // Therefore σ_t from timing_energy_bins.root already represents the pure
    // crystal timing.  The MCP correction does NOT improve that estimate.
    //
    // The correction IS meaningful for:
    //   (a) Per-channel timing  (σ ≈ 230–280 ps from summary.root)
    //   (b) Absolute combination methods (mean-all-8, A²-wgt-all-8 from
    //       timing_summary.root) where t_MCP does not cancel
    //
    // We apply the correction to (a) via gTimingResolution from summary.root.
    // =========================================================================

    std::cout << "\n  NOTE: (DW-UP)/2 estimator is MCP-jitter-free — "
                 "correction applies to per-channel absolute timing.\n\n";

    TString tSumPath  = "output/Summary/summary.root";
    TFile*  fSum      = TFile::Open(tSumPath);
    if (!fSum || fSum->IsZombie()) {
        std::cerr << "mcpJitter: cannot open " << tSumPath
                  << " — run analyzeResolution.C first\n";
        return;
    }
    TGraphErrors* gAvgCh = static_cast<TGraphErrors*>(
        fSum->Get("gTimingResolution"));
    if (!gAvgCh) {
        std::cerr << "mcpJitter: gTimingResolution not found\n";
        fSum->Close(); delete fSum; return;
    }

    // gCorrected = per-channel average corrected for MCP jitter
    // gUncorrected = same uncorrected
    TGraphErrors* gCorrected   = new TGraphErrors();
    TGraphErrors* gUncorrected = new TGraphErrors();

    for (int p = 0; p < gAvgCh->GetN(); ++p) {
        double E       = gAvgCh->GetX()[p];
        double sigMeas = gAvgCh->GetY()[p];     // [ps]
        double sigMeasErr = gAvgCh->GetEY()[p];

        int n = gUncorrected->GetN();
        gUncorrected->SetPoint(n, E, sigMeas);
        gUncorrected->SetPointError(n, 0., sigMeasErr);

        int ir = -1;
        for (int r = 0; r < kNRuns; ++r)
            if (std::fabs(kRuns[r].energy_GeV - E) < 1.) { ir = r; break; }

        if (ir < 0 || sigMCP[ir] <= 0.) continue;

        double sig2_meas = sigMeas * sigMeas;
        double sig2_mcp  = sigMCP[ir] * sigMCP[ir];

        if (sig2_meas <= sig2_mcp) {
            std::cout << "  [WARN] σ_MCP >= σ_measured at E=" << E
                      << " GeV (unexpected for per-channel timing)\n";
            continue;
        }

        double sigCorr    = std::sqrt(sig2_meas - sig2_mcp);
        double sigCorrErr = sigMeasErr * sigMeas / sigCorr;

        int nc = gCorrected->GetN();
        gCorrected->SetPoint(nc, E, sigCorr);
        gCorrected->SetPointError(nc, 0., sigCorrErr);

        std::cout << Form("  E=%.0f GeV: σ_ch_meas=%.1f  σ_MCP=%.1f"
                          "  σ_crystal=%.1f ps  (gain: %.1f ps)\n",
                          E, sigMeas, sigMCP[ir], sigCorr, sigMeas - sigCorr);
    }

    // =========================================================================
    // Output PDF
    // =========================================================================
    TString outPDF  = "output/Summary/mcp_jitter.pdf";
    TString outROOT = "output/Summary/mcp_jitter.root";
    TCanvas c("c_mj", "", 960, 720);

    // ── Page 1: MCP1-MCP2 Δt distributions (overlaid, normalised) ────────────
    c.Clear(); c.cd(); StylePad();
    double yMaxD = 0.;
    for (int r = 0; r < kNRuns; ++r) {
        if (!hDiff[r]) continue;
        double integral = hDiff[r]->Integral();
        if (integral > 0.) hDiff[r]->Scale(1.0 / integral);
        yMaxD = std::max(yMaxD, hDiff[r]->GetMaximum());
    }

    bool first = true;
    TLegend leg1(0.64, 0.50, 0.93, 0.78);
    leg1.SetBorderSize(0); leg1.SetFillStyle(0); leg1.SetTextSize(0.040);
    for (int r = 0; r < kNRuns; ++r) {
        if (!hDiff[r]) continue;
        hDiff[r]->SetLineColor(kREnergyCols[r]);
        hDiff[r]->SetLineWidth(2);
        hDiff[r]->GetXaxis()->SetTitle("t_{MCP1} - t_{MCP2}  (ns)");
        hDiff[r]->GetYaxis()->SetTitle("Events (normalised)");
        hDiff[r]->GetYaxis()->SetRangeUser(0., 1.20 * yMaxD);
        if (first) { hDiff[r]->Draw("HIST"); first = false; }
        else          hDiff[r]->Draw("HIST SAME");
        leg1.AddEntry(hDiff[r],
                      Form("%.0f GeV  (#sigma_{1ch}=%.1f ps)",
                           kRuns[r].energy_GeV, sigMCP[r]), "l");
    }
    leg1.Draw();
    {
        TLatex ann; ann.SetNDC(); ann.SetTextSize(0.036); ann.SetTextColor(kGray+1);
        ann.DrawLatex(0.18, 0.20,
            "Both MCP quality cuts applied (200 < peak < 750 mV)");
        ann.DrawLatex(0.18, 0.14,
            "#sigma_{MCP,single} = #sigma(MCP1#minusMCP2) / #sqrt{2}");
    }
    PageTitle("MCP1 #minus MCP2 time difference  (normalised to unit area)");
    c.Print(outPDF + "(");

    // ── Page 2: σ_MCP_single vs beam energy ──────────────────────────────────
    c.Clear(); c.cd(); StylePad();

    TGraphErrors* gMCPsig = new TGraphErrors();
    for (int r = 0; r < kNRuns; ++r) {
        if (sigMCP[r] <= 0.) continue;
        int n = gMCPsig->GetN();
        gMCPsig->SetPoint(n, kRuns[r].energy_GeV, sigMCP[r]);
        gMCPsig->SetPointError(n, 0., sigMCPErr[r]);
    }
    gMCPsig->SetLineColor(kBlue+1); gMCPsig->SetMarkerColor(kBlue+1);
    gMCPsig->SetMarkerStyle(20); gMCPsig->SetMarkerSize(1.3);
    gMCPsig->SetLineWidth(2);
    gMCPsig->Draw("APL");
    gMCPsig->GetXaxis()->SetTitle("Beam energy (GeV)");
    gMCPsig->GetYaxis()->SetTitle("#sigma_{MCP,single}  (ps)");
    gMCPsig->GetXaxis()->SetLimits(15., 165.);

    // A flat MCP jitter vs energy is expected (intrinsic timing, not shower-dependent)
    {
        double mean = 0.; int nn = 0;
        for (int r = 0; r < kNRuns; ++r)
            if (sigMCP[r] > 0.) { mean += sigMCP[r]; ++nn; }
        if (nn > 0) {
            mean /= nn;
            TLine* lMean = new TLine(15., mean, 165., mean);
            lMean->SetLineColor(kGray+1); lMean->SetLineStyle(2); lMean->SetLineWidth(2);
            lMean->Draw("SAME");
            TLatex ann; ann.SetNDC(); ann.SetTextSize(0.042); ann.SetTextColor(kGray+1);
            ann.DrawLatex(0.18, 0.78, Form("mean = %.1f ps", mean));
        }
    }
    PageTitle("#sigma_{MCP,single} = #sigma(MCP1#minusMCP2)/#sqrt{2} vs beam energy");
    c.Print(outPDF);

    // ── Page 2b: σ_MCP_single vs MCP1 amplitude at 150 GeV ───────────────────
    c.Clear(); c.cd(); StylePad();

    if (gMCPvsAmp->GetN() > 0) {
        gMCPvsAmp->SetLineColor(kBlue+1); gMCPvsAmp->SetMarkerColor(kBlue+1);
        gMCPvsAmp->SetMarkerStyle(20); gMCPvsAmp->SetMarkerSize(1.2);
        gMCPvsAmp->SetLineWidth(2);
        gMCPvsAmp->Draw("APL");
        gMCPvsAmp->GetXaxis()->SetTitle("MCP1 amplitude (mV)");
        gMCPvsAmp->GetYaxis()->SetTitle("#sigma_{MCP,single}  (ps)");
        gMCPvsAmp->GetXaxis()->SetLimits(static_cast<double>(kMCP1_minPeak) - 20.,
                                         static_cast<double>(kMCP1_maxPeak) + 20.);

        if (fMCPvA && gMCPvsAmp->GetN() >= 3) {
            fMCPvA->SetLineColor(kRed+1);
            fMCPvA->SetLineStyle(2);
            fMCPvA->SetLineWidth(2);
            fMCPvA->Draw("SAME");

            TLatex ann2b; ann2b.SetNDC(); ann2b.SetTextSize(0.040);
            ann2b.SetTextColor(kRed+1);
            ann2b.DrawLatex(0.18, 0.80,
                Form("Fit: #sqrt{p_{0}^{2}/A^{2} + p_{1}^{2}}"));
            ann2b.DrawLatex(0.18, 0.73,
                Form("p_{0} = %.1f ps mV  (noise term)", mcp_p0));  // #cdot is not a valid ROOT TLatex token
            ann2b.DrawLatex(0.18, 0.66,
                Form("p_{1} = %.1f ps  (electronic floor)", mcp_p1));

            TLatex note2b; note2b.SetNDC();
            note2b.SetTextSize(0.032); note2b.SetTextColor(kGray+1);
            note2b.DrawLatex(0.18, 0.20,
                "Used to apply per-event MCP correction in amplitude-binned analyses");
        }
    } else {
        TLatex nodata; nodata.SetNDC(); nodata.SetTextSize(0.050);
        nodata.SetTextColor(kGray+1); nodata.SetTextAlign(22);
        nodata.DrawLatex(0.50, 0.50, "No 150 GeV data available for amplitude analysis");
    }
    PageTitle("#sigma_{MCP,single} vs MCP1 amplitude at 150 GeV");
    c.Print(outPDF);

    // ── Page 3: Corrected vs uncorrected σ_t (energy-binned m0) ──────────────
    c.Clear(); c.cd(); StylePad();

    double yMax3 = 0.;
    for (int p = 0; p < gUncorrected->GetN(); ++p)
        yMax3 = std::max(yMax3, gUncorrected->GetY()[p]);
    for (int p = 0; p < gCorrected->GetN(); ++p)
        yMax3 = std::max(yMax3, gCorrected->GetY()[p]);
    if (yMax3 < 10.) yMax3 = 100.;

    TH1F* frame3 = static_cast<TH1F*>(
        c.DrawFrame(15., 0., 165., 1.40 * yMax3,
                    ";Beam energy (GeV);#sigma_{t} (ps)"));
    frame3->GetXaxis()->SetTitleSize(0.050);
    frame3->GetYaxis()->SetTitleSize(0.050);
    frame3->GetYaxis()->SetTitleOffset(1.30);

    gUncorrected->SetLineColor(kRed+1); gUncorrected->SetMarkerColor(kRed+1);
    gUncorrected->SetMarkerStyle(20); gUncorrected->SetMarkerSize(1.3);
    gUncorrected->SetLineWidth(2);
    gUncorrected->Draw("PL SAME");

    gCorrected->SetLineColor(kGreen+2); gCorrected->SetMarkerColor(kGreen+2);
    gCorrected->SetMarkerStyle(21); gCorrected->SetMarkerSize(1.3);
    gCorrected->SetLineWidth(2); gCorrected->SetLineStyle(2);
    gCorrected->Draw("PL SAME");

    // Also plot MCP jitter for context
    gMCPsig->SetLineColor(kBlue+1); gMCPsig->SetMarkerColor(kBlue+1);
    gMCPsig->SetMarkerStyle(22); gMCPsig->SetMarkerSize(1.1);
    gMCPsig->SetLineWidth(1); gMCPsig->SetLineStyle(3);
    gMCPsig->Draw("PL SAME");

    TLegend leg3(0.65, 0.65, 0.93, 0.88);
    leg3.SetBorderSize(0); leg3.SetFillStyle(0); leg3.SetTextSize(0.042);
    leg3.AddEntry(gUncorrected, "#sigma_{measured}  (avg over 8 HG channels)", "lp");
    leg3.AddEntry(gCorrected,
                  "#sigma_{crystal} = #sqrt{#sigma_{meas}^{2} #minus #sigma_{MCP}^{2}}", "lp");
    leg3.AddEntry(gMCPsig, "#sigma_{MCP,single}", "lp");
    leg3.Draw();

    PageTitle("MCP jitter subtraction: measured vs crystal-only #sigma_{t}");
    c.Print(outPDF);

    // ── Page 3c: combination methods, MCP-subtracted intrinsic timing ────────
    //
    // The MCP-referenced combination methods (timing_summary.root gTiming_M0-M4)
    // carry the FULL MCP jitter as a common-mode term: every channel shares the
    // same t_MCP1, which factors out of the weighted mean, so
    //     σ_combo² = σ_combo,intrinsic² + σ_MCP1² .
    // Subtracting σ_MCP1 (amplitude-averaged σ_MCP_single, pinned by the
    // amplitude-dependent fit on Page 2b) reveals the detector-only combination
    // timing.  This is delicate when σ_MCP ≈ σ_combo: errors are propagated and
    // points with σ_MCP ≥ σ_combo are flagged rather than silently dropped.
    TGraphErrors* gIntrM[5] = {};
    {
        TString tsPath = "output/Summary/timing_summary.root";
        TFile*  fTS    = TFile::Open(tsPath);
        if (!fTS || fTS->IsZombie()) {
            std::cout << "  [Page 3c] " << tsPath
                      << " not found — run timingResolution.C first; skipping.\n";
        } else {
            static const char* mLab[5] = {
                "Best single", "4-corner (U+D)/2", "Mean all 8",
                "A^{2}-wgt all 8", "A^{2}-wgt corners"
            };
            static const int mCol[5] = { kRPurple, kRTeal, kROrange, kRRed, kRGreen };

            c.Clear(); c.cd(); StylePad(false, true);
            double yMaxC = 0.;
            // Build measured + intrinsic graphs per method
            TGraphErrors* gMeasM[5] = {};
            int nFlagged = 0;   // σ_MCP ≥ σ_combo: subtraction impossible
            int nDelicate = 0;  // σ_MCP > 0.85·σ_combo: combo is MCP-limited,
                                // intrinsic value highly sensitive to σ_MCP model
            for (int m = 0; m < 5; ++m) {
                TGraphErrors* gM = static_cast<TGraphErrors*>(
                    fTS->Get(Form("gTiming_M%d", m)));
                if (!gM) continue;
                gMeasM[m] = new TGraphErrors();
                gIntrM[m] = new TGraphErrors();
                for (int p = 0; p < gM->GetN(); ++p) {
                    const double E   = gM->GetX()[p];
                    const double sM  = gM->GetY()[p];
                    const double sME = gM->GetEY()[p];
                    // σ_MCP for this energy
                    int ir = -1;
                    for (int r = 0; r < kNRuns; ++r)
                        if (std::fabs(kRuns[r].energy_GeV - E) < 1.) { ir = r; break; }
                    if (ir < 0 || sigMCP[ir] <= 0.) continue;
                    const double sMcp  = sigMCP[ir];
                    const double sMcpE = sigMCPErr[ir];

                    const int nMeas = gMeasM[m]->GetN();
                    gMeasM[m]->SetPoint(nMeas, E, sM);
                    gMeasM[m]->SetPointError(nMeas, 0., sME);
                    yMaxC = std::max(yMaxC, sM);

                    if (sM <= sMcp) { ++nFlagged; continue; }   // unphysical: flag, skip
                    if (sMcp > 0.85 * sM) ++nDelicate;          // MCP-limited regime
                    const double sIntr = std::sqrt(sM*sM - sMcp*sMcp);
                    const double sIntrE = (sIntr > 0.)
                        ? std::sqrt(sM*sM*sME*sME + sMcp*sMcp*sMcpE*sMcpE) / sIntr
                        : 0.;
                    const int nI = gIntrM[m]->GetN();
                    gIntrM[m]->SetPoint(nI, E, sIntr);
                    gIntrM[m]->SetPointError(nI, 0., sIntrE);
                }
            }
            if (yMaxC < 10.) yMaxC = 200.;

            TH1F* frC = static_cast<TH1F*>(
                c.DrawFrame(15., 0., 165., 1.30 * yMaxC,
                            ";Beam energy (GeV);#sigma_{t} (ps)"));
            frC->GetXaxis()->SetTitleSize(0.048);
            frC->GetYaxis()->SetTitleSize(0.048);
            frC->GetYaxis()->SetTitleOffset(1.30);

            TLegend* Lc = MakeLegend(6);
            for (int m = 0; m < 5; ++m) {
                if (gMeasM[m] && gMeasM[m]->GetN() > 0) {
                    gMeasM[m]->SetLineColor(mCol[m]);   gMeasM[m]->SetMarkerColor(mCol[m]);
                    gMeasM[m]->SetMarkerStyle(24);      gMeasM[m]->SetMarkerSize(1.0);
                    gMeasM[m]->SetLineStyle(2);         gMeasM[m]->SetLineWidth(1);
                    gMeasM[m]->Draw("PL SAME");
                }
                if (gIntrM[m] && gIntrM[m]->GetN() > 0) {
                    gIntrM[m]->SetLineColor(mCol[m]);   gIntrM[m]->SetMarkerColor(mCol[m]);
                    gIntrM[m]->SetMarkerStyle(20);      gIntrM[m]->SetMarkerSize(1.2);
                    gIntrM[m]->SetLineWidth(2);
                    gIntrM[m]->Draw("PL SAME");
                    Lc->AddEntry(gIntrM[m], mLab[m], "lp");
                }
            }
            Lc->Draw();
            {
                // Compact, non-overlapping annotations.  This pad uses a right-side
                // sidebar legend (StylePad(false,true)), so all text is kept short
                // and left of x≈0.66, and the bottom caveats are merged into two
                // short lines (previously four long lines overlapped each other).
                TLatex a; a.SetNDC();
                // Headline (top-left, clear of the sidebar legend) — short form.
                if (gMeasM[3] && gMeasM[3]->GetN() > 0 &&
                    gIntrM[3] && gIntrM[3]->GetN() > 0) {
                    const int lm = gMeasM[3]->GetN() - 1;
                    const int li = gIntrM[3]->GetN() - 1;
                    a.SetTextColor(kRRed); a.SetTextSize(0.033);
                    a.DrawLatex(0.15, 0.93, Form("A^{2}-combo @ %.0f GeV:  "
                        "%.0f #rightarrow %.0f#pm%.0f ps (MCP-subtracted)",
                        gMeasM[3]->GetX()[lm], gMeasM[3]->GetY()[lm],
                        gIntrM[3]->GetY()[li], gIntrM[3]->GetEY()[li]));
                }
                a.SetTextColor(kGray+2); a.SetTextSize(0.028);
                a.DrawLatex(0.15, 0.27, "open = measured   filled = intrinsic (MCP subtracted)");
                a.DrawLatex(0.15, 0.22,
                    "#sigma_{intr} = #sqrt{#sigma_{combo}^{2} #minus #sigma_{MCP}^{2}}   (MCP common-mode)");
                if (nFlagged > 0 || nDelicate > 0) {
                    a.SetTextColor(kRRed);
                    a.DrawLatex(0.15, 0.15,
                        "Combos are MCP-limited (#sigma_{MCP} #approx #sigma_{combo}):");
                    a.DrawLatex(0.15, 0.10,
                        "intrinsic unreliable -- quote the MCP-free (DW#minusUP)/2 result.");
                }
            }
            PageTitle("Combination timing: measured (MCP-limited) vs detector-only intrinsic");
            c.Print(outPDF);
            fTS->Close(); delete fTS;
        }
    }

    // ── Page 4: Summary table ─────────────────────────────────────────────────
    c.Clear(); c.cd();
    gPad->SetLeftMargin(0.05); gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.08); gPad->SetBottomMargin(0.05);

    // Pre-compute mean MCP1 amplitude per energy for the table column
    // (use the same event selection as the amplitude analysis, all energies)
    double meanAmpPerRun[6] = {};
    {
        // Re-open ntuples briefly to collect mean MCP1 amplitude per run
        for (int r = 0; r < kNRuns; ++r) {
            TTree* ta = nullptr;
            TFile* fa = OpenNtuple(r, ta);
            if (!fa) continue;
            double xc2, yc2, tcfd2, trms2;
            ScanRunCenters(ta, xc2, yc2, tcfd2, trms2);
            Float_t mp, mt, mp2, mt2, xt, yt; Bool_t wok;
            ta->SetBranchAddress(ta->GetBranch("mcp1_peak")?"mcp1_peak":"mcp_peak",  &mp);
            ta->SetBranchAddress(ta->GetBranch("mcp1_time")?"mcp1_time":"mcp_time",  &mt);
            ta->SetBranchAddress("mcp2_peak", &mp2);
            ta->SetBranchAddress("mcp2_time", &mt2);
            ta->SetBranchAddress("x_trk",     &xt);
            ta->SetBranchAddress("y_trk",     &yt);
            ta->SetBranchAddress("wc_ok",     &wok);
            double ampSum = 0.; long nA = 0;
            Long64_t nEa = ta->GetEntries();
            for (Long64_t ev = 0; ev < nEa; ++ev) {
                ta->GetEntry(ev);
                if (!wok) continue;
                if (mp  < kMCP1_minPeak || mp  > kMCP1_maxPeak) continue;
                if (mp2 < kMCP2_minPeak || mp2 > kMCP2_maxPeak) continue;
                float dx = xt - static_cast<float>(xc2);
                float dy = yt - static_cast<float>(yc2);
                if (std::sqrt(dx*dx + dy*dy) >= static_cast<float>(TimingFiducialR(kRuns[r].energy_GeV))) continue;
                if (mt < -1e5f || mt2 < -1e5f) continue;
                ampSum += static_cast<double>(mp); ++nA;
            }
            ta->ResetBranchAddresses();
            fa->Close(); delete fa;
            if (nA > 0) meanAmpPerRun[r] = ampSum / nA;
        }
    }

    TLatex tab; tab.SetNDC();
    // Header (6 columns, tighter spacing)
    tab.SetTextSize(0.040); tab.SetTextColor(kBlack);
    tab.DrawLatex(0.04, 0.88, "E (GeV)");
    tab.DrawLatex(0.19, 0.88, "#sigma_{MCP} (ps)");
    tab.DrawLatex(0.36, 0.88, "#sigma_{MCP}(A) (ps)");
    tab.DrawLatex(0.53, 0.88, "#sigma_{meas} (ps)");
    tab.DrawLatex(0.68, 0.88, "#sigma_{crys} (ps)");
    tab.DrawLatex(0.84, 0.88, "Frac (%)");

    TLine* lHead = new TLine(0.04, 0.85, 0.96, 0.85);
    lHead->SetLineColor(kGray+1); lHead->SetLineWidth(1); lHead->Draw();

    // Data rows
    double yRow = 0.79;
    for (int r = 0; r < kNRuns; ++r) {
        tab.SetTextSize(0.038); tab.SetTextColor(kBlack);
        tab.DrawLatex(0.04, yRow, Form("%.0f", kRuns[r].energy_GeV));

        if (sigMCP[r] > 0.)
            tab.DrawLatex(0.19, yRow, Form("%.1f #pm %.1f", sigMCP[r], sigMCPErr[r]));
        else
            tab.DrawLatex(0.19, yRow, "#minus");

        // Amplitude-dependent σ_MCP at mean amplitude for this energy
        if (fMCPvA && meanAmpPerRun[r] > 0.) {
            double sigA = fMCPvA->Eval(meanAmpPerRun[r]);
            tab.SetTextColor(kBlue+1);
            tab.DrawLatex(0.36, yRow, Form("%.1f", sigA));
            tab.SetTextColor(kBlack);
        } else {
            tab.DrawLatex(0.36, yRow, "#minus");
        }

        // Find matching measured σ_t
        double sigMeas = -1.;
        for (int p = 0; p < gUncorrected->GetN(); ++p)
            if (std::fabs(gUncorrected->GetX()[p] - kRuns[r].energy_GeV) < 1.)
                { sigMeas = gUncorrected->GetY()[p]; break; }

        double sigCryst = -1.;
        for (int p = 0; p < gCorrected->GetN(); ++p)
            if (std::fabs(gCorrected->GetX()[p] - kRuns[r].energy_GeV) < 1.)
                { sigCryst = gCorrected->GetY()[p]; break; }

        if (sigMeas > 0.) {
            tab.DrawLatex(0.53, yRow, Form("%.1f", sigMeas));
            if (sigCryst > 0.) {
                tab.SetTextColor(kGreen+2);
                tab.DrawLatex(0.68, yRow, Form("%.1f", sigCryst));
                // MCP fraction = σ_MCP / σ_measured (in quadrature)
                double frac = (sigMCP[r] / sigMeas) * 100.;
                tab.SetTextColor(sigMCP[r]/sigMeas > 0.3 ? kRed+1 : kGray+1);
                tab.DrawLatex(0.84, yRow, Form("%.1f", frac));
                tab.SetTextColor(kBlack);
            } else {
                tab.DrawLatex(0.68, yRow, "n/a");
                tab.DrawLatex(0.84, yRow, "n/a");
            }
        } else {
            tab.DrawLatex(0.53, yRow, "#minus");
            tab.DrawLatex(0.68, yRow, "#minus");
            tab.DrawLatex(0.84, yRow, "#minus");
        }
        yRow -= 0.090;
    }

    TLine* lFoot = new TLine(0.04, yRow + 0.02, 0.96, yRow + 0.02);
    lFoot->SetLineColor(kGray+1); lFoot->SetLineWidth(1); lFoot->Draw();

    {
        TLatex note; note.SetNDC(); note.SetTextSize(0.030); note.SetTextColor(kGray+1);
        note.DrawLatex(0.04, yRow - 0.03,
            "Red fraction (>30%): MCP reference is limiting the measurement.");
        note.DrawLatex(0.04, yRow - 0.08,
            "#sigma_{MCP}(A): amplitude-dependent fit f(A)=#sqrt{p_{0}^{2}/A^{2}+p_{1}^{2}}"
            " evaluated at mean MCP1 amplitude.");
        note.DrawLatex(0.04, yRow - 0.13,
            "#sigma_{crystal} = #sqrt{#sigma_{meas}^{2} #minus #sigma_{MCP,single}^{2}} "
            "where #sigma_{MCP,single} = #sigma(MCP1#minusMCP2)/#sqrt{2}");
    }
    PageTitle("MCP jitter summary table");
    c.Print(outPDF + ")");

    // =========================================================================
    // Write summary ROOT file
    // =========================================================================
    TFile* fOut = new TFile(outROOT, "RECREATE");
    gMCPsig->Write("gMCP_jitter");
    gUncorrected->Write("gTEB_measured");
    gCorrected->Write("gTEB_corrected");
    gMCPvsAmp->Write("gMCPvsAmp");
    if (fMCPvA) {
        fMCPvA->Write("fMCPvA");
        (new TParameter<Double_t>("mcp_jitter_p0", mcp_p0))->Write();
        (new TParameter<Double_t>("mcp_jitter_p1", mcp_p1))->Write();
    }
    // MCP-subtracted intrinsic combination timing (detector-only, MCP common-mode removed)
    for (int m = 0; m < 5; ++m)
        if (gIntrM[m] && gIntrM[m]->GetN() > 0)
            gIntrM[m]->Write(Form("gIntrinsic_M%d", m));
    fOut->Close();
    delete fOut;

    // Cleanup
    for (int r = 0; r < kNRuns; ++r) delete hDiff[r];
    delete gMCPsig; delete gCorrected; delete gUncorrected;
    delete gMCPvsAmp;
    delete fMCPvA;
    fSum->Close(); delete fSum;

    std::cout << "mcpJitter: done -> " << outPDF << "\n";
    std::cout << "           root  -> " << outROOT << "\n";
}

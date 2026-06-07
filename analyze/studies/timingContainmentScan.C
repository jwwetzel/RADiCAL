// ============================================================================
// timingContainmentScan.C — containment cut optimisation
// ============================================================================
//
// Scans the PbGlass containment cut (sum_pb < ratio × sum_lg) over 7 threshold
// values and measures how the A²-weighted combination timing resolution σ_t
// changes at each beam energy.
//
// Motivation: investigatePbGlass.C found that at 150 GeV, 13.6% of in-fiducial
// events with a valid RADiCAL signal are hadronic punch-through.  This scan
// answers: does tightening the containment cut improve σ_t, and at what
// efficiency cost?
//
// Timing estimator: A²-weighted combination of 8 HG CFD-5% channels (M3).
//   t_combo = Σ(A_i² × hg_cfd[i]) / Σ(A_i²)
//
// Containment thresholds scanned:  0.10, 0.15, 0.20, 0.25, 0.30 (baseline),
//                                   0.40, 0.50
//
// ── Physics note ─────────────────────────────────────────────────────────────
//
//   The containment cut rejects events where sum_pb ≥ ratio × sum_lg, i.e.
//   the PbGlass sees ≥ (ratio×100)% of the RADiCAL LG signal.  This removes:
//
//     • Hadronic punch-through: long-range pions that leak past RADiCAL
//     • Edge showers: poorly-contained EM showers at the module boundary
//
//   Tighter (lower ratio) → purer EM sample, but smaller yield.
//   The optimal ratio minimises σ_t × 1/√N  (resolution × 1/√statistics).
//
//   A cut that is too tight removes genuine EM events (false positive rate ↑).
//   A cut that is too loose includes hadronic contamination (tail → σ↑).
//
// ── Output ───────────────────────────────────────────────────────────────────
//   output/Summary/timing_containment_scan.pdf   (4 pages)
//   output/Summary/timing_containment_scan.root
//     → gSigT_rXX   σ_t [ps] vs energy for each threshold XX
//     → gYield_rXX  surviving fraction vs energy for each threshold XX
//
//   Page 1: σ_t vs energy — all 7 thresholds overlaid
//   Page 2: σ_t vs threshold at each energy (6 curves, x=ratio)
//   Page 3: Surviving fraction vs threshold at each energy
//   Page 4: σ_t distributions at 150 GeV: tight (0.15), baseline (0.30), loose (0.50)
//
// Usage (from repository root):
//   root -l -b -q 'Analysis/timingContainmentScan.C+'
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

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

// ---------------------------------------------------------------------------
// Containment thresholds to scan
// ---------------------------------------------------------------------------
static const double kScanRatio[7] = {0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50};
static const int    kNScan        = 7;
static const int    kBaselineIdx  = 4;   // index of baseline 0.30

// Colors: hot=tight → cool=loose (red → violet)
static const int kScanCol[7]  = {
    kRed+1, kOrange+1, kOrange-3, kGreen+2, kCyan+2, kBlue+1, kViolet+1
};
static const int kScanMark[7] = { 20, 21, 22, 23, 24, 25, 26 };

// Energy colors from RADiCALStyle.h (kREnergyCols, set by ApplyRADiCALStyle)
static const int kEMark[6] = { 20, 21, 22, 23, 24, 25 };

// Sentinel for missing timing
static const float kNoTime = -1e9f;

// ---------------------------------------------------------------------------
// OpenNtuple — open per-energy ntuple, return the "rad" tree
// ---------------------------------------------------------------------------
static TFile* OpenNtuple(int r, TTree*& tree)
{
    TString path = Form("output/%s/ntuple.root", kRuns[r].label.Data());
    TFile* f = TFile::Open(path);
    if (!f || f->IsZombie()) {
        std::cerr << "timingContainmentScan: cannot open " << path << "\n";
        tree = nullptr; return nullptr;
    }
    tree = static_cast<TTree*>(f->Get("rad"));
    if (!tree) {
        std::cerr << "timingContainmentScan: tree 'rad' missing in " << path << "\n";
        f->Close(); delete f; tree = nullptr; return nullptr;
    }
    return f;
}

// ---------------------------------------------------------------------------
// A2WeightedCombo — A²-weighted mean of hg_cfd[8]
//
// Returns kNoTime if fewer than 2 channels pass the amplitude threshold.
// The 2-channel minimum guards against single-channel artifacts biasing
// the timing estimate in low-signal events.
// ---------------------------------------------------------------------------
static float A2WeightedCombo(const Float_t hg_peak[8], const Float_t hg_cfd[8])
{
    double sw = 0., swt = 0.;
    int    nValid = 0;
    for (int ic = 0; ic < kNCap; ++ic) {
        if (hg_peak[ic] < kHG_minPeak) continue;
        if (hg_cfd[ic]  < -1e5f)       continue;  // CFD-failed sentinel
        double w  = static_cast<double>(hg_peak[ic]) * hg_peak[ic];
        sw  += w;
        swt += w * hg_cfd[ic];
        ++nValid;
    }
    if (nValid < 2 || sw <= 0.) return kNoTime;
    return static_cast<float>(swt / sw);
}

// ---------------------------------------------------------------------------
// StatsFromVec — compute mean and population σ from a float vector
// ---------------------------------------------------------------------------
static void StatsFromVec(const std::vector<float>& v,
                          double& mean, double& sigma)
{
    mean = sigma = 0.;
    if (v.empty()) return;
    for (float x : v) mean += x;
    mean /= static_cast<double>(v.size());
    for (float x : v) sigma += (x - mean) * (x - mean);
    sigma = std::sqrt(sigma / static_cast<double>(v.size()));
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
void timingContainmentScan()
{
    TH1::AddDirectory(kFALSE);
    ApplyRADiCALStyle();   // RADiCALStyle.h via PlotUtils.h

    gSystem->mkdir("output/Summary", kTRUE);
    std::cout << "timingContainmentScan: scanning containment cut\n";

    // Results (ps) and yield fractions indexed [run][threshold]
    double sigT   [kNRuns][kNScan] = {};
    double sigTErr[kNRuns][kNScan] = {};
    double yield  [kNRuns][kNScan] = {};

    // Timing distributions at 150 GeV for the comparison page
    // Thresholds shown: 0.15 (k=1), 0.30 baseline (k=4), 0.50 (k=6)
    static const int kShow150[3] = {1, 4, 6};
    TH1F* hShow150[3] = {};

    // =========================================================================
    // Event loop — one pass per run, fill timing vectors for each threshold
    // =========================================================================
    for (int r = 0; r < kNRuns; ++r) {
        TTree* t = nullptr;
        TFile* fin = OpenNtuple(r, t);
        if (!fin) continue;

        double xc, yc, tcfd, trms;
        ScanRunCenters(t, xc, yc, tcfd, trms);

        Float_t x_trk, y_trk, mcp_peak, sum_lg, sum_pb;
        Float_t hg_peak[8], hg_cfd[8];
        Bool_t  wc_ok;
        t->SetBranchAddress("x_trk",    &x_trk);
        t->SetBranchAddress("y_trk",    &y_trk);
        t->SetBranchAddress("wc_ok",    &wc_ok);
        t->SetBranchAddress("mcp_peak", &mcp_peak);
        t->SetBranchAddress("sum_lg",   &sum_lg);
        t->SetBranchAddress("sum_pb",   &sum_pb);
        t->SetBranchAddress("hg_peak",   hg_peak);
        // CFD-5% (adopted headline fraction); guarded fallback to CFD-20%.
        t->SetBranchAddress(t->GetBranch("hg_cfd05") ? "hg_cfd05" : "hg_cfd", hg_cfd);

        // Per-threshold combo-timing accumulator and baseline count
        std::vector<float> tvec[kNScan];
        for (int k = 0; k < kNScan; ++k) tvec[k].reserve(20000);
        std::vector<float> tvecAll;   // all valid combos — used for histogram range
        tvecAll.reserve(20000);
        int nBase = 0;  // events with valid t_combo (denominator for yield)

        Long64_t nEv = t->GetEntries();
        for (Long64_t ev = 0; ev < nEv; ++ev) {
            t->GetEntry(ev);

            // Standard event quality cuts
            if (!wc_ok)                                               continue;
            if (mcp_peak < kMCP1_minPeak || mcp_peak > kMCP1_maxPeak) continue;
            float dx = x_trk - static_cast<float>(xc);
            float dy = y_trk - static_cast<float>(yc);
            if (std::sqrt(dx*dx + dy*dy) >= static_cast<float>(kFiducial_r_timing))
                continue;

            // Require valid A²-weighted combination timing
            float tcombo = A2WeightedCombo(hg_peak, hg_cfd);
            if (tcombo <= kNoTime + 1e5f) continue;

            ++nBase;
            tvecAll.push_back(tcombo);

            // Apply containment cut at each scan threshold.
            // The cut is only applied when there is enough RADiCAL signal
            // (sum_lg > kSumLG_centroid).  Events with too little signal are
            // excluded regardless — they are beam halo (see ClassifyEvent in
            // investigatePbGlass.C, population C).
            for (int k = 0; k < kNScan; ++k) {
                if (sum_lg <= kSumLG_centroid) continue;
                if (sum_pb >= static_cast<float>(kScanRatio[k]) * sum_lg) continue;
                tvec[k].push_back(tcombo);
            }
        }

        t->ResetBranchAddresses();
        fin->Close(); delete fin;

        if (nBase < 50 || tvecAll.empty()) {
            std::cout << "  " << kRuns[r].label << ": too few events\n";
            continue;
        }

        // Histogram range from the full (no containment cut) distribution
        double tMean, tSig;
        StatsFromVec(tvecAll, tMean, tSig);
        if (tSig < 0.020) tSig = 0.200;  // safety floor: 200 ps
        double histLo = tMean - 4. * tSig;
        double histHi = tMean + 4. * tSig;

        // Fit each threshold
        for (int k = 0; k < kNScan; ++k) {
            if (static_cast<int>(tvec[k].size()) < 50) {
                std::cout << "  " << kRuns[r].label
                          << " ratio=" << kScanRatio[k]
                          << ": too few events (" << tvec[k].size() << ")\n";
                continue;
            }

            TH1F h(Form("hCS_%d_%d", r, k), "",
                   200, histLo, histHi);
            h.SetDirectory(nullptr);
            for (float v : tvec[k]) h.Fill(v);

            double fitMu, fitMuErr, fitSig, fitSigErr;
            FitGaussCore(&h, 2.0, fitMu, fitMuErr, fitSig, fitSigErr);

            if (fitSig > 0.) {
                sigT[r][k]    = fitSig * 1000.;  // convert ns → ps
                sigTErr[r][k] = fitSigErr * 1000.;
            }
            yield[r][k] = static_cast<double>(tvec[k].size()) /
                          static_cast<double>(nBase);

            // Save the 150 GeV representative distributions
            if (r == kNRuns - 1) {
                for (int s = 0; s < 3; ++s) {
                    if (k == kShow150[s]) {
                        hShow150[s] = new TH1F(h);
                        hShow150[s]->SetName(Form("hShow150_%d", s));
                        hShow150[s]->SetDirectory(nullptr);
                    }
                }
            }
        }

        std::cout << "  " << kRuns[r].label
                  << Form(": baseline (30%%) σ_t = %.1f ps"
                          "  yield = %.3f  (%d events)\n",
                          sigT[r][kBaselineIdx],
                          yield[r][kBaselineIdx], nBase);
    }

    // =========================================================================
    // Build TGraphErrors for each threshold and each energy
    // =========================================================================

    // gSigT[k]  : σ_t [ps] vs energy for threshold k  (one per scan point)
    // gYield[k] : yield vs energy for threshold k
    TGraphErrors* gSigT [kNScan] = {};
    TGraphErrors* gYield[kNScan] = {};

    double yMaxSigT = 0., yMaxYield = 0.;
    for (int k = 0; k < kNScan; ++k) {
        gSigT[k]  = new TGraphErrors();
        gYield[k] = new TGraphErrors();
        for (int r = 0; r < kNRuns; ++r) {
            if (sigT[r][k] <= 0.) continue;
            int n = gSigT[k]->GetN();
            gSigT[k]->SetPoint(n, kRuns[r].energy_GeV, sigT[r][k]);
            gSigT[k]->SetPointError(n, 0., sigTErr[r][k]);
            yMaxSigT = std::max(yMaxSigT, sigT[r][k]);

            int m = gYield[k]->GetN();
            gYield[k]->SetPoint(m, kRuns[r].energy_GeV, yield[r][k]);
            gYield[k]->SetPointError(m, 0., 0.);
            yMaxYield = std::max(yMaxYield, yield[r][k]);
        }
        gSigT[k]->SetLineColor(kScanCol[k]);
        gSigT[k]->SetMarkerColor(kScanCol[k]);
        gSigT[k]->SetMarkerStyle(kScanMark[k]);
        gSigT[k]->SetMarkerSize(1.2);
        gSigT[k]->SetLineWidth(2);

        gYield[k]->SetLineColor(kScanCol[k]);
        gYield[k]->SetMarkerColor(kScanCol[k]);
        gYield[k]->SetMarkerStyle(kScanMark[k]);
        gYield[k]->SetMarkerSize(1.2);
        gYield[k]->SetLineWidth(2);
    }
    if (yMaxSigT < 10.) yMaxSigT = 400.;

    // gSigT_vs_ratio[r]  : σ_t vs threshold for energy r
    // gYield_vs_ratio[r] : yield vs threshold for energy r
    TGraphErrors* gSigTV[kNRuns]  = {};
    TGraphErrors* gYieldV[kNRuns] = {};
    for (int r = 0; r < kNRuns; ++r) {
        gSigTV[r]  = new TGraphErrors();
        gYieldV[r] = new TGraphErrors();
        for (int k = 0; k < kNScan; ++k) {
            if (sigT[r][k] <= 0.) continue;
            int n = gSigTV[r]->GetN();
            gSigTV[r]->SetPoint(n, kScanRatio[k], sigT[r][k]);
            gSigTV[r]->SetPointError(n, 0., sigTErr[r][k]);

            int m = gYieldV[r]->GetN();
            gYieldV[r]->SetPoint(m, kScanRatio[k], yield[r][k]);
            gYieldV[r]->SetPointError(m, 0., 0.);
        }
        gSigTV[r]->SetLineColor(kREnergyCols[r]);
        gSigTV[r]->SetMarkerColor(kREnergyCols[r]);
        gSigTV[r]->SetMarkerStyle(kEMark[r]);
        gSigTV[r]->SetMarkerSize(1.2);
        gSigTV[r]->SetLineWidth(2);

        gYieldV[r]->SetLineColor(kREnergyCols[r]);
        gYieldV[r]->SetMarkerColor(kREnergyCols[r]);
        gYieldV[r]->SetMarkerStyle(kEMark[r]);
        gYieldV[r]->SetMarkerSize(1.2);
        gYieldV[r]->SetLineWidth(2);
    }

    // =========================================================================
    // Output PDF
    // =========================================================================
    TString outPDF  = "output/Summary/timing_containment_scan.pdf";
    TString outROOT = "output/Summary/timing_containment_scan.root";
    TCanvas c("c_tcs", "", 960, 720);

    // ── Page 1: σ_t vs energy, 7 containment thresholds overlaid ─────────────
    // Use sidebar mode (right margin = 0.27) to keep the 7-entry legend
    // out of the data area.
    c.Clear(); c.cd(); StylePad(false, true);

    TH1F* frame1 = static_cast<TH1F*>(
        c.DrawFrame(15., 50., 165., 1.35 * yMaxSigT,
                    ";Beam energy (GeV);A^{2}-weighted combo #sigma_{t}  (ps)"));
    frame1->GetXaxis()->SetTitleSize(0.050);
    frame1->GetYaxis()->SetTitleSize(0.050);
    frame1->GetYaxis()->SetTitleOffset(1.30);

    TLegend* leg1 = MakeLegend(8);
    leg1->SetHeader("p_{pb}/p_{lg} threshold", "C");

    for (int k = 0; k < kNScan; ++k) {
        gSigT[k]->Draw("PL SAME");
        TString label = Form("%.2f", kScanRatio[k]);
        if (k == kBaselineIdx) label += "  (baseline)";
        leg1->AddEntry(gSigT[k], label, "lp");
    }
    leg1->Draw();

    // Guidance annotations
    {
        TLatex ann; ann.SetNDC(); ann.SetTextSize(0.034); ann.SetTextColor(kGray+1);
        ann.DrawLatex(0.18, 0.16,
            "Tighter cut (lower ratio) #rightarrow purer EM sample");
        ann.DrawLatex(0.18, 0.10,
            "Looser cut (higher ratio) #rightarrow more statistics");
    }
    PageTitle("A^{2}-weighted combo #sigma_{t} vs energy  |  containment threshold scan");
    c.Print(outPDF + "(");

    // ── Page 2: σ_t vs threshold at each energy ───────────────────────────────
    c.Clear(); c.cd(); StylePad();

    // Find y range across all energies and thresholds
    double yMin2 = 1e9, yMax2 = 0.;
    for (int r = 0; r < kNRuns; ++r) {
        for (int k = 0; k < kNScan; ++k) {
            if (sigT[r][k] > 0.) {
                yMin2 = std::min(yMin2, sigT[r][k]);
                yMax2 = std::max(yMax2, sigT[r][k]);
            }
        }
    }
    if (yMax2 < 10.) { yMin2 = 0.; yMax2 = 400.; }
    double yPad2 = 0.15 * (yMax2 - yMin2);

    TH1F* frame2 = static_cast<TH1F*>(
        c.DrawFrame(0.05, std::max(0., yMin2 - yPad2), 0.55,
                    yMax2 + 3. * yPad2,
                    ";Containment threshold (sum_{pb}/sum_{lg});#sigma_{t}  (ps)"));
    frame2->GetXaxis()->SetTitleSize(0.048);
    frame2->GetYaxis()->SetTitleSize(0.048);
    frame2->GetYaxis()->SetTitleOffset(1.30);

    // Dashed vertical line at the baseline (0.30)
    TLine* lBase = new TLine(kPb_maxRatio, std::max(0., yMin2 - yPad2),
                              kPb_maxRatio, yMax2 + 3. * yPad2);
    lBase->SetLineColor(kGray+1); lBase->SetLineStyle(2); lBase->SetLineWidth(2);
    lBase->Draw("SAME");
    {
        TLatex ann; ann.SetNDC(); ann.SetTextSize(0.034); ann.SetTextColor(kGray+1);
        ann.DrawLatex(0.58, 0.92, "dashed = current baseline (0.30)");
    }

    TLegend leg2(0.63, 0.45, 0.93, 0.78);
    leg2.SetBorderSize(0); leg2.SetFillStyle(0); leg2.SetTextSize(0.038);
    for (int r = 0; r < kNRuns; ++r) {
        gSigTV[r]->Draw("PL SAME");
        leg2.AddEntry(gSigTV[r], Form("%.0f GeV", kRuns[r].energy_GeV), "lp");
    }
    leg2.Draw();
    PageTitle("#sigma_{t} vs containment threshold at each beam energy");
    c.Print(outPDF);

    // ── Page 3: Surviving fraction vs threshold ───────────────────────────────
    c.Clear(); c.cd(); StylePad();

    TH1F* frame3 = static_cast<TH1F*>(
        c.DrawFrame(0.05, 0., 0.55, 1.10,
                    ";Containment threshold (sum_{pb}/sum_{lg});Surviving fraction"));
    frame3->GetXaxis()->SetTitleSize(0.048);
    frame3->GetYaxis()->SetTitleSize(0.048);
    frame3->GetYaxis()->SetTitleOffset(1.30);

    lBase = new TLine(kPb_maxRatio, 0., kPb_maxRatio, 1.10);
    lBase->SetLineColor(kGray+1); lBase->SetLineStyle(2); lBase->SetLineWidth(2);
    lBase->Draw("SAME");

    TLegend leg3(0.63, 0.15, 0.93, 0.48);
    leg3.SetBorderSize(0); leg3.SetFillStyle(0); leg3.SetTextSize(0.038);
    for (int r = 0; r < kNRuns; ++r) {
        gYieldV[r]->Draw("PL SAME");
        leg3.AddEntry(gYieldV[r], Form("%.0f GeV", kRuns[r].energy_GeV), "lp");
    }
    leg3.Draw();
    {
        TLatex ann; ann.SetNDC(); ann.SetTextSize(0.034); ann.SetTextColor(kGray+1);
        ann.DrawLatex(0.18, 0.88,
            "Denominator: events with wc_ok + MCP quality + fiducial + valid t_{combo}");
        ann.DrawLatex(0.18, 0.82,
            "(excludes beam halo with sum_{lg} #leq 300 mV)");
    }
    PageTitle("Surviving fraction vs containment threshold at each beam energy");
    c.Print(outPDF);

    // ── Page 4: σ_t distributions at 150 GeV — tight / baseline / loose ──────
    c.Clear(); c.cd(); StylePad();

    // Scale all three histograms to unit area for fair shape comparison
    const char* showLabel[3] = {"tight (0.15)", "baseline (0.30)", "loose (0.50)"};
    const int   showCol[3]   = { kRed+1, kCyan+2, kViolet+1 };

    double yMax4 = 0.;
    for (int s = 0; s < 3; ++s) {
        if (!hShow150[s]) continue;
        double intg = hShow150[s]->Integral();
        if (intg > 0.) hShow150[s]->Scale(1.0 / intg);
        yMax4 = std::max(yMax4, hShow150[s]->GetMaximum());
    }
    if (yMax4 <= 0.) yMax4 = 1.;

    TLegend leg4(0.57, 0.52, 0.93, 0.85);
    leg4.SetBorderSize(0); leg4.SetFillStyle(0); leg4.SetTextSize(0.042);

    bool first4 = true;
    for (int s = 0; s < 3; ++s) {
        if (!hShow150[s]) continue;
        hShow150[s]->SetLineColor(showCol[s]);
        hShow150[s]->SetLineWidth(2);
        hShow150[s]->GetXaxis()->SetTitle("t_{combo}  (ns)");
        hShow150[s]->GetYaxis()->SetTitle("Events (normalised)");
        hShow150[s]->GetYaxis()->SetRangeUser(0., 1.25 * yMax4);
        if (first4) { hShow150[s]->Draw("HIST"); first4 = false; }
        else           hShow150[s]->Draw("HIST SAME");

        // Annotate σ
        double mu4, muErr4, sig4, sigErr4;
        FitGaussCore(hShow150[s], 2.0, mu4, muErr4, sig4, sigErr4);
        leg4.AddEntry(hShow150[s],
                      Form("%s  #sigma=%.0f ps", showLabel[s], sig4 * 1000.), "l");
    }
    leg4.Draw();
    {
        TLatex ann; ann.SetNDC(); ann.SetTextSize(0.040);
        ann.DrawLatex(0.18, 0.88, "150 GeV A^{2}-weighted combo timing");
        ann.SetTextColor(kGray+1); ann.SetTextSize(0.034);
        ann.DrawLatex(0.18, 0.82, "(normalised to unit area)");
    }
    PageTitle("150 GeV: combo timing distribution  |  tight / baseline / loose containment cut");
    c.Print(outPDF + ")");

    // =========================================================================
    // Write summary ROOT file
    // =========================================================================
    TFile* fOut = new TFile(outROOT, "RECREATE");
    for (int k = 0; k < kNScan; ++k) {
        int ipct = static_cast<int>(kScanRatio[k] * 100. + 0.5);
        gSigT[k]->Write(Form("gSigT_r%02d", ipct));
        gYield[k]->Write(Form("gYield_r%02d", ipct));
    }
    for (int r = 0; r < kNRuns; ++r) {
        gSigTV[r]->Write(Form("gSigT_vs_ratio_%s",
                               kRuns[r].label.Data()));
        gYieldV[r]->Write(Form("gYield_vs_ratio_%s",
                                kRuns[r].label.Data()));
    }
    fOut->Close(); delete fOut;

    // Cleanup
    for (int k = 0; k < kNScan; ++k) { delete gSigT[k]; delete gYield[k]; }
    for (int r = 0; r < kNRuns; ++r) { delete gSigTV[r]; delete gYieldV[r]; }
    for (int s = 0; s < 3; ++s) delete hShow150[s];

    std::cout << "timingContainmentScan: done -> " << outPDF << "\n";
    std::cout << "                       root  -> " << outROOT << "\n";
}

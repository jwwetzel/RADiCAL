// ============================================================================
// channelCombinationScan.C — brute-force channel combination optimisation
// ============================================================================
//
// Tests all 255 non-empty subsets of the 8 HG capillary channels and measures
// the A²-weighted combination timing resolution σ_t for each subset at each
// beam energy.
//
// There are 2⁸ - 1 = 255 subsets.  The search is exhaustive; no pruning is
// applied so that unexpected combinations are not missed.
//
// ── Key question ─────────────────────────────────────────────────────────────
//
//   Channel 7 (SW-Up) uses MCP2 as its timing reference, while channels 0–6
//   all use MCP1.  When channels 0–6 and channel 7 are combined with A²-
//   weighting, channel 7 adds noise proportional to the MCP1–MCP2 jitter
//   (≈72 ps) in addition to its intrinsic crystal timing noise.
//
//   If ch7 hurts the combination, dropping it gives a 7-channel combination
//   that should outperform the naive 8-channel average.
//
//   This macro measures that effect directly at every beam energy.
//
// ── Timing estimator ─────────────────────────────────────────────────────────
//
//   For a given channel subset S:
//     t_combo(S) = Σ_{i∈S, valid} A_i² × hg_cfd[i] / Σ_{i∈S, valid} A_i²
//   where "valid" means hg_peak[i] > kHG_minPeak AND hg_cfd[i] not flagged.
//   A subset with no valid channels is skipped.
//
// ── Output ───────────────────────────────────────────────────────────────────
//   Analysis/Output/Summary/channel_combination_scan.pdf   (4 pages)
//   Analysis/Output/Summary/channel_combination_scan.root
//
//   Page 1: σ_t vs energy — best combination for each subset size (N=1..8)
//   Page 2: σ_t scatter at 150 GeV — all 255 subsets, colored by N channels
//           and whether SW-Up (ch7) is included
//   Page 3: SW-Up impact — (best-7 no ch7) vs (all-8) vs (best-7 any) vs energy
//   Page 4: Single-channel performance at each energy (σ_t bar chart)
//
// Usage (from repository root):
//   root -l -b -q 'Analysis/channelCombinationScan.C+'
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
// SubsetLabel — human-readable name for a channel bitmask
// ---------------------------------------------------------------------------
static TString SubsetLabel(int mask)
{
    TString s;
    for (int ic = 0; ic < kNCap; ++ic) {
        if (!(mask & (1 << ic))) continue;
        if (!s.IsNull()) s += "+";
        s += kCap[ic].name;
    }
    return s;
}

// ---------------------------------------------------------------------------
// NCh — number of bits set in a mask
// ---------------------------------------------------------------------------
static int NCh(int mask)
{
    int n = 0;
    for (int i = 0; i < 8; ++i) if (mask & (1<<i)) ++n;
    return n;
}

// ---------------------------------------------------------------------------
// OpenNtuple — open per-energy ntuple, return the "rad" tree
// ---------------------------------------------------------------------------
static TFile* OpenNtuple(int r, TTree*& tree)
{
    TString path = Form("Analysis/Output/%s/ntuple.root", kRuns[r].label.Data());
    TFile* f = TFile::Open(path);
    if (!f || f->IsZombie()) {
        std::cerr << "channelCombinationScan: cannot open " << path << "\n";
        tree = nullptr; return nullptr;
    }
    tree = static_cast<TTree*>(f->Get("rad"));
    if (!tree) {
        std::cerr << "channelCombinationScan: tree 'rad' missing in " << path << "\n";
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

// Color scheme: N=1..8 channels (warm → cool)
static const int kNChCol[9] = {
    0,               // [0] unused
    kRed+1,          // N=1
    kOrange+1,       // N=2
    kOrange-3,       // N=3
    kYellow+2,       // N=4
    kGreen+2,        // N=5
    kCyan+2,         // N=6
    kBlue+1,         // N=7
    kViolet+1        // N=8
};
static const int kNChMark[9] = { 0, 20, 21, 22, 23, 24, 25, 26, 27 };

// Energy colors from RADiCALStyle.h (kREnergyCols, set by ApplyRADiCALStyle)
static const int kEMark[6] = { 20, 21, 22, 23, 24, 25 };

// Bitmask of the 7-channel "no SW-Up" combination (all channels except bit 7)
static const int kMaskNoSWU = (1 << kNCap) - 1 - (1 << 7);  // = 0x7F = 127
static const int kMaskAll8  = (1 << kNCap) - 1;              // = 0xFF = 255

// ===========================================================================
// Main
// ===========================================================================
void channelCombinationScan()
{
    TH1::AddDirectory(kFALSE);
    ApplyRADiCALStyle();   // RADiCALStyle.h via PlotUtils.h

    gSystem->mkdir("Analysis/Output/Summary", kTRUE);
    std::cout << "channelCombinationScan: testing all 255 channel subsets\n";

    // sigT[mask][r]  — σ_t [ps] for subset `mask` at energy index `r`
    // mask=0 unused; valid range 1..255
    static double sigT   [256][kNRuns];
    static double sigTErr[256][kNRuns];
    for (int m = 0; m < 256; ++m)
        for (int r = 0; r < kNRuns; ++r)
            sigT[m][r] = sigTErr[m][r] = 0.;

    // =========================================================================
    // Main loop — one pass per energy run
    // =========================================================================
    for (int r = 0; r < kNRuns; ++r) {
        TTree* t = nullptr;
        TFile* fin = OpenNtuple(r, t);
        if (!fin) continue;

        double xc, yc, tcfd, trms;
        ScanRunCenters(t, xc, yc, tcfd, trms);

        // Histogram range: ±1.5 ns around the per-channel CFD mean.
        // This covers ±6σ for a single-channel distribution (σ ≈ 250 ps)
        // and ±15σ for the 8-channel combo (σ ≈ 100 ps).
        // 300 bins → 10 ps/bin.
        const double kHalfRange = 1.5;  // ns
        double histLo = tcfd - kHalfRange;
        double histHi = tcfd + kHalfRange;

        // Allocate 255 histograms (mask 1..255)
        TH1F* hCombo[256] = {};
        for (int mask = 1; mask < 256; ++mask) {
            hCombo[mask] = new TH1F(Form("hC_%d_%d", r, mask), "",
                                    300, histLo, histHi);
            hCombo[mask]->SetDirectory(nullptr);
        }

        // Set up branches
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
        // CFD-5% (adopted headline fraction); guarded fallback to CFD-20%.  This
        // brings the brute-force subset scan onto the same basis as the headline
        // and timingResolution.C (reconciling the former CFD-20% ~78 ps all-8 combo
        // with the CFD-5% ~62 ps one).
        t->SetBranchAddress(t->GetBranch("hg_cfd05") ? "hg_cfd05" : "hg_cfd", hg_cfd);

        Long64_t nEv = t->GetEntries();
        for (Long64_t ev = 0; ev < nEv; ++ev) {
            t->GetEntry(ev);

            // Standard event quality + containment cuts
            if (!wc_ok)                                               continue;
            if (mcp_peak < kMCP1_minPeak || mcp_peak > kMCP1_maxPeak) continue;
            float dx = x_trk - static_cast<float>(xc);
            float dy = y_trk - static_cast<float>(yc);
            if (std::sqrt(dx*dx + dy*dy) >= static_cast<float>(kFiducial_r_timing))
                continue;
            if (sum_lg > kSumLG_centroid &&
                sum_pb >= kPb_maxRatio * sum_lg) continue;

            // Pre-compute per-channel A²-weighted contributions for this event.
            // Channel with hg_peak below threshold or invalid CFD gets w=0
            // and is silently excluded from any subset's combination.
            double w [8] = {};
            double wt[8] = {};
            for (int ic = 0; ic < kNCap; ++ic) {
                if (hg_peak[ic] < kHG_minPeak) continue;
                if (hg_cfd[ic]  < -1e5f)       continue;  // CFD-failed
                w [ic] = static_cast<double>(hg_peak[ic]) * hg_peak[ic];
                wt[ic] = w[ic] * hg_cfd[ic];
            }

            // Fill all 255 subset histograms in a single inner loop
            for (int mask = 1; mask < 256; ++mask) {
                double sw = 0., swt = 0.;
                for (int ic = 0; ic < 8; ++ic) {
                    if (!(mask & (1 << ic))) continue;
                    sw  += w[ic];
                    swt += wt[ic];
                }
                if (sw <= 0.) continue;  // no valid channels in this subset
                hCombo[mask]->Fill(swt / sw);
            }
        }

        t->ResetBranchAddresses();
        fin->Close(); delete fin;

        // Fit all 255 histograms and store results
        int nGoodFits = 0;
        for (int mask = 1; mask < 256; ++mask) {
            double mu, muErr, sig, sigErr;
            FitGaussCore(hCombo[mask], 2.0, mu, muErr, sig, sigErr);
            if (sig > 0.) {
                sigT[mask][r]    = sig * 1000.;   // ns → ps
                sigTErr[mask][r] = sigErr * 1000.;
                ++nGoodFits;
            }
            delete hCombo[mask];
            hCombo[mask] = nullptr;
        }

        // Report baseline (all-8) and no-SW-Up (7-ch) results
        std::cout << "  " << kRuns[r].label
                  << Form(": all-8 σ_t = %.1f ps"
                          "  7ch(no SW-U) σ_t = %.1f ps"
                          "  (%d/%d fits ok)\n",
                          sigT[kMaskAll8][r],
                          sigT[kMaskNoSWU][r],
                          nGoodFits, 255);
    }

    // =========================================================================
    // Derived quantities: best-N combination per energy
    // =========================================================================
    // bestMask[N][r]  — bitmask of the best N-channel combination at energy r
    // bestSigT[N][r]  — its σ_t [ps]
    int    bestMask[9][kNRuns] = {};
    double bestSigT[9][kNRuns] = {};
    for (int N = 1; N <= 8; ++N)
        for (int r = 0; r < kNRuns; ++r)
            bestSigT[N][r] = 1e9;

    for (int mask = 1; mask < 256; ++mask) {
        int N = NCh(mask);
        for (int r = 0; r < kNRuns; ++r) {
            if (sigT[mask][r] > 0. && sigT[mask][r] < bestSigT[N][r]) {
                bestSigT[N][r] = sigT[mask][r];
                bestMask[N][r] = mask;
            }
        }
    }

    // Report best combinations at 150 GeV (r=5)
    std::cout << "\n  Best combinations at 150 GeV:\n";
    for (int N = 1; N <= 8; ++N) {
        if (bestMask[N][kNRuns-1] == 0) continue;
        std::cout << Form("  N=%d  σ_t = %.1f ps  [%s]\n",
                          N, bestSigT[N][kNRuns-1],
                          SubsetLabel(bestMask[N][kNRuns-1]).Data());
    }

    // =========================================================================
    // Build TGraphErrors for output
    // =========================================================================

    // Page 1: best σ_t vs energy for each N
    TGraphErrors* gBestN[9] = {};
    for (int N = 1; N <= 8; ++N) {
        gBestN[N] = new TGraphErrors();
        for (int r = 0; r < kNRuns; ++r) {
            if (bestSigT[N][r] >= 1e8) continue;
            int n = gBestN[N]->GetN();
            gBestN[N]->SetPoint(n, kRuns[r].energy_GeV, bestSigT[N][r]);
            gBestN[N]->SetPointError(n, 0., sigTErr[bestMask[N][r]][r]);
        }
        gBestN[N]->SetLineColor(kNChCol[N]);
        gBestN[N]->SetMarkerColor(kNChCol[N]);
        gBestN[N]->SetMarkerStyle(kNChMark[N]);
        gBestN[N]->SetMarkerSize(1.2);
        gBestN[N]->SetLineWidth(2);
    }

    // Page 3: SW-Up impact — all-8, 7-ch-no-SW-U, best-7-any vs energy
    TGraphErrors* gAll8   = new TGraphErrors();
    TGraphErrors* gNoSWU  = new TGraphErrors();
    TGraphErrors* gBest7  = new TGraphErrors();
    double yMaxP3 = 0.;
    for (int r = 0; r < kNRuns; ++r) {
        double E = kRuns[r].energy_GeV;
        auto addPt = [&](TGraphErrors* g, double sig, double sigE) {
            if (sig <= 0.) return;
            int n = g->GetN();
            g->SetPoint(n, E, sig);
            g->SetPointError(n, 0., sigE);
            yMaxP3 = std::max(yMaxP3, sig);
        };
        addPt(gAll8,  sigT[kMaskAll8][r],  sigTErr[kMaskAll8][r]);
        addPt(gNoSWU, sigT[kMaskNoSWU][r], sigTErr[kMaskNoSWU][r]);
        addPt(gBest7, bestSigT[7][r],
              sigTErr[bestMask[7][r]][r]);
    }
    auto styleG = [](TGraphErrors* g, int col, int mark, int lw = 2, int ls = 1) {
        g->SetLineColor(col); g->SetMarkerColor(col);
        g->SetMarkerStyle(mark); g->SetMarkerSize(1.3);
        g->SetLineWidth(lw); g->SetLineStyle(ls);
    };
    styleG(gAll8,  kBlue+1,   20);
    styleG(gNoSWU, kGreen+2,  21);
    styleG(gBest7, kRed+1,    22, 2, 2);
    if (yMaxP3 < 10.) yMaxP3 = 400.;

    // =========================================================================
    // Output PDF
    // =========================================================================
    TString outPDF  = "Analysis/Output/Summary/channel_combination_scan.pdf";
    TString outROOT = "Analysis/Output/Summary/channel_combination_scan.root";
    TCanvas c("c_ccs", "", 960, 720);

    // ── Page 1: best-N σ_t vs energy ─────────────────────────────────────────
    c.Clear(); c.cd(); StylePad();

    double yMaxP1 = 0.;
    for (int N = 1; N <= 8; ++N)
        for (int r = 0; r < kNRuns; ++r)
            if (bestSigT[N][r] < 1e8) yMaxP1 = std::max(yMaxP1, bestSigT[N][r]);
    if (yMaxP1 < 10.) yMaxP1 = 400.;

    TH1F* frame1 = static_cast<TH1F*>(
        c.DrawFrame(15., 50., 165., 1.35 * yMaxP1,
                    ";Beam energy (GeV);Best A^{2}-combo #sigma_{t} for N channels  (ps)"));
    frame1->GetXaxis()->SetTitleSize(0.048);
    frame1->GetYaxis()->SetTitleSize(0.048);
    frame1->GetYaxis()->SetTitleOffset(1.35);

    TLegend leg1(0.65, 0.48, 0.93, 0.88);
    leg1.SetBorderSize(0); leg1.SetFillStyle(0); leg1.SetTextSize(0.038);
    leg1.SetHeader("Best combination of N channels", "C");
    for (int N = 1; N <= 8; ++N) {
        if (!gBestN[N] || gBestN[N]->GetN() == 0) continue;
        gBestN[N]->Draw("PL SAME");
        TString entryLabel = Form("N = %d", N);
        if (N == 8) entryLabel += " (all)";
        leg1.AddEntry(gBestN[N], entryLabel, "lp");
    }
    leg1.Draw();
    PageTitle("Best A^{2}-combo #sigma_{t} vs energy by number of channels");
    c.Print(outPDF + "(");

    // ── Page 2: σ_t scatter at 150 GeV (all 255 subsets) ────────────────────
    c.Clear(); c.cd(); StylePad();

    // Determine y range from 150 GeV data
    double yMin2 = 1e9, yMax2 = 0.;
    for (int mask = 1; mask < 256; ++mask) {
        if (sigT[mask][kNRuns-1] > 0.) {
            yMin2 = std::min(yMin2, sigT[mask][kNRuns-1]);
            yMax2 = std::max(yMax2, sigT[mask][kNRuns-1]);
        }
    }
    if (yMax2 < 10.) { yMin2 = 50.; yMax2 = 400.; }

    // Use x range 0–9 so ROOT places integer tick labels at 1, 2, …, 8.
    // (Range 0.5–8.5 with 8 divisions would put labels at non-integer positions.)
    TH1F* frame2 = static_cast<TH1F*>(
        c.DrawFrame(0., std::max(50., yMin2 * 0.90),
                    9., yMax2 * 1.20,
                    ";Number of channels in combination;#sigma_{t}  (ps)  [150 GeV]"));
    frame2->GetXaxis()->SetNdivisions(9, kFALSE);
    frame2->GetXaxis()->SetTitleSize(0.050);
    frame2->GetYaxis()->SetTitleSize(0.050);
    frame2->GetYaxis()->SetTitleOffset(1.30);

    // Draw all 255 subsets as markers; color by N channels, open if ch7 excluded
    for (int mask = 1; mask < 256; ++mask) {
        if (sigT[mask][kNRuns-1] <= 0.) continue;
        int    N   = NCh(mask);
        bool   has7 = (mask & (1 << 7));
        int    col  = kNChCol[N];
        int    mark = has7 ? 20 : 24;   // filled = has ch7, open = no ch7
        double xjit = (N + (has7 ? +0.15 : -0.15));  // tiny x-jitter for readability

        TGraph gPt(1, &xjit, &sigT[mask][kNRuns-1]);
        gPt.SetMarkerColor(col);
        gPt.SetMarkerStyle(mark);
        gPt.SetMarkerSize(0.7);
        gPt.DrawClone("P SAME");
    }

    // Highlight key masks
    auto markKey = [&](int mask, int col, int mark, double xOff) {
        if (sigT[mask][kNRuns-1] <= 0.) return;
        double xp = NCh(mask) + xOff;
        double yp = sigT[mask][kNRuns-1];
        TGraph g1(1, &xp, &yp);
        g1.SetMarkerColor(col); g1.SetMarkerStyle(mark); g1.SetMarkerSize(2.0);
        g1.DrawClone("P SAME");
    };
    markKey(kMaskAll8,  kBlue+1,   29, +0.30);  // all-8
    markKey(kMaskNoSWU, kGreen+2,  29, -0.30);  // no SW-U

    TLegend leg2(0.53, 0.65, 0.93, 0.88);
    leg2.SetBorderSize(0); leg2.SetFillStyle(0); leg2.SetTextSize(0.034);
    leg2.AddEntry((TObject*)nullptr,
                  "Filled circles: subset includes SW-Up (ch7, MCP2 ref)", "");
    leg2.AddEntry((TObject*)nullptr,
                  "Open circles:   subset excludes SW-Up", "");
    leg2.AddEntry((TObject*)nullptr,
                  Form("Large #bullet (blue):  all-8  %.1f ps",
                       sigT[kMaskAll8][kNRuns-1]), "");
    leg2.AddEntry((TObject*)nullptr,
                  Form("Large #bullet (green): no SW-Up (7ch)  %.1f ps",
                       sigT[kMaskNoSWU][kNRuns-1]), "");
    leg2.Draw();

    PageTitle("150 GeV: #sigma_{t} for all 255 channel subsets (filled = incl. SW-Up)");
    c.Print(outPDF);

    // ── Page 3: SW-Up impact across all energies ──────────────────────────────
    c.Clear(); c.cd(); StylePad();

    TH1F* frame3 = static_cast<TH1F*>(
        c.DrawFrame(15., 50., 165., 1.35 * yMaxP3,
                    ";Beam energy (GeV);A^{2}-combo #sigma_{t}  (ps)"));
    frame3->GetXaxis()->SetTitleSize(0.050);
    frame3->GetYaxis()->SetTitleSize(0.050);
    frame3->GetYaxis()->SetTitleOffset(1.30);

    gAll8->Draw("PL SAME");
    gNoSWU->Draw("PL SAME");
    gBest7->Draw("PL SAME");

    TLegend leg3(0.65, 0.67, 0.93, 0.88);
    leg3.SetBorderSize(0); leg3.SetFillStyle(0); leg3.SetTextSize(0.040);
    leg3.AddEntry(gAll8,  "All 8 channels (baseline)", "lp");
    leg3.AddEntry(gNoSWU, "7 channels, no SW-Up  (ch7 removed)", "lp");
    leg3.AddEntry(gBest7, "Best 7-channel combination", "lp");
    leg3.Draw();

    // Annotate delta at 150 GeV — use colour matching the relevant graph
    {
        double d = sigT[kMaskAll8][kNRuns-1] - sigT[kMaskNoSWU][kNRuns-1];
        TLatex ann; ann.SetNDC(); ann.SetTextSize(0.042);
        ann.SetTextColor(d > 0. ? kGreen+2 : kRed+1);
        ann.DrawLatex(0.18, 0.20,
                      Form("150 GeV: removing SW-Up %s by %.1f ps",
                           d > 0. ? "improves" : "worsens",
                           std::fabs(d)));
    }
    PageTitle("SW-Up (ch7, MCP2 reference) impact on combination #sigma_{t}");
    c.Print(outPDF);

    // ── Page 4: Single-channel σ_t at all energies (bar chart overview) ──────
    c.Clear(); c.cd();
    gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.13);
    gPad->SetRightMargin(0.06); gPad->SetTopMargin(0.12);
    gPad->SetTickx(1); gPad->SetTicky(1);

    // Build per-channel graphs (σ_t vs energy, one per channel)
    TGraphErrors* gCh[8] = {};
    double yMaxP4 = 0.;
    for (int ic = 0; ic < kNCap; ++ic) {
        int mask = (1 << ic);
        gCh[ic] = new TGraphErrors();
        for (int r = 0; r < kNRuns; ++r) {
            if (sigT[mask][r] <= 0.) continue;
            int n = gCh[ic]->GetN();
            gCh[ic]->SetPoint(n, kRuns[r].energy_GeV, sigT[mask][r]);
            gCh[ic]->SetPointError(n, 0., sigTErr[mask][r]);
            yMaxP4 = std::max(yMaxP4, sigT[mask][r]);
        }
        gCh[ic]->SetLineColor(kREnergyCols[ic % kNRuns]);
        gCh[ic]->SetMarkerColor(kREnergyCols[ic % kNRuns]);
        gCh[ic]->SetMarkerStyle(20 + ic);
        gCh[ic]->SetMarkerSize(1.2);
        gCh[ic]->SetLineWidth(2);
    }
    if (yMaxP4 < 10.) yMaxP4 = 500.;

    TH1F* frame4 = static_cast<TH1F*>(
        c.DrawFrame(15., 50., 165., 1.35 * yMaxP4,
                    ";Beam energy (GeV);Single-channel #sigma_{t}  (ps)"));
    frame4->GetXaxis()->SetTitleSize(0.050);
    frame4->GetYaxis()->SetTitleSize(0.050);
    frame4->GetYaxis()->SetTitleOffset(1.20);

    TLegend leg4(0.65, 0.42, 0.93, 0.88);
    leg4.SetBorderSize(0); leg4.SetFillStyle(0); leg4.SetTextSize(0.038);
    for (int ic = 0; ic < kNCap; ++ic) {
        gCh[ic]->Draw("PL SAME");
        TString label = kCap[ic].name;
        if (ic == 7) label += " (MCP2 ref)";
        leg4.AddEntry(gCh[ic], label, "lp");
    }
    leg4.Draw();
    PageTitle("Single-channel A^{2}-weighted #sigma_{t} vs beam energy");
    c.Print(outPDF + ")");

    // =========================================================================
    // Write summary ROOT file
    // =========================================================================
    TFile* fOut = new TFile(outROOT, "RECREATE");
    for (int N = 1; N <= 8; ++N) {
        if (gBestN[N]) gBestN[N]->Write(Form("gBestN%d", N));
    }
    gAll8->Write("gAll8");
    gNoSWU->Write("gNoSWU");
    gBest7->Write("gBest7");
    for (int ic = 0; ic < kNCap; ++ic) {
        if (gCh[ic]) gCh[ic]->Write(Form("gCh%d_%s", ic, kCap[ic].name));
    }
    fOut->Close(); delete fOut;

    // Cleanup
    for (int N = 1; N <= 8; ++N) delete gBestN[N];
    delete gAll8; delete gNoSWU; delete gBest7;
    for (int ic = 0; ic < kNCap; ++ic) delete gCh[ic];

    std::cout << "channelCombinationScan: done -> " << outPDF << "\n";
    std::cout << "                        root  -> " << outROOT << "\n";
}

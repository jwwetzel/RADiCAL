// ============================================================================
// investigatePbGlass.C — PbGlass dual-band structure investigation
// ============================================================================
//
// The PbGlass vs ΣLG scatter plot shows two distinct populations:
//
//   Band 1 (horizontal, near sum_pb ≈ 0):
//     Fully contained EM showers — high ΣLG, negligible PbGlass signal.
//     These are our signal events.
//
//   Band 2 (near-vertical, sum_lg ≈ 0):
//     Beam halo particles that missed the RADiCAL crystal bundle entirely.
//     They may or may not hit the PbGlass.
//
//   Above the cut line (scattered):
//     Hadronic events (pion contamination) or edge showers where the shower
//     is not fully contained within RADiCAL.
//
// This macro classifies every WC-OK event into four exclusive populations
// and studies each separately to answer:
//   (1) Does the current 30% containment cut correctly isolate the signal?
//   (2) How much does hadronic contamination degrade σ_t?
//   (3) Do halo events pass our fiducial cut?
//   (4) Does the edge effect change with energy?
//
// ── Population definitions ───────────────────────────────────────────────────
//
//   A — "Good EM shower":      sum_lg > kSumLG_centroid  AND
//                               sum_pb < kPb_maxRatio × sum_lg
//       → Our signal.  Should have the best timing.
//
//   B — "Punch-through / hadronic":  sum_lg > kSumLG_centroid  AND
//                                     sum_pb >= kPb_maxRatio × sum_lg
//       → Rejected by the current containment cut.
//       → Quantify the timing degradation to justify the cut.
//
//   C — "Beam halo (missed RADiCAL)":  sum_lg <= kSumLG_centroid  AND
//                                       r_beam < kFiducial_r_timing
//       → Inside fiducial spatially but no RADiCAL signal.
//       → These events MUST be cut; question is whether they currently
//         slip through the containment cut.
//
//   D — "Off-axis":  r_beam >= kFiducial_r_timing
//       → Outside fiducial; included for spatial context only.
//
// ── Output ───────────────────────────────────────────────────────────────────
//
//   Analysis/Output/Summary/pbglass_investigation.pdf  (6 pages)
//
//   Page 1: Scatter (ΣPbGlass vs ΣLG) colour-coded by population — 2×3 grid
//   Page 2: Hit maps for populations A, B, C at 25 GeV and 150 GeV — 3×2 grid
//   Page 3: Timing distributions (A vs B) per energy — 2×3 grid
//   Page 4: σ_t for A / B / (A+B) vs energy — headline result
//   Page 5: Containment ratio distributions A vs B per energy — 2×3 grid
//   Page 6: ΣLG spectra A vs B overlaid per energy — 2×3 grid
//
// Usage (from repository root):
//   root -l -b -q 'Analysis/investigatePbGlass.C+'
// ============================================================================

#include "ChannelConfig.h"
#include "PlotUtils.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TEllipse.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TAxis.h"
#include "TF1.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

// ---------------------------------------------------------------------------
// Population identifiers and visual style
// ---------------------------------------------------------------------------
static const int kNPop = 4;

enum Pop { kA = 0, kB = 1, kC = 2, kD = 3 };

static const char* kPopName[kNPop] = {
    "A: Good EM",
    "B: Punch-through",
    "C: Beam halo",
    "D: Off-axis"
};
static const int kPopCol[kNPop] = {
    kGreen+2,    // A: good — green
    kRed+1,      // B: hadronic — red
    kOrange+1,   // C: halo — orange
    kGray+1      // D: off-axis — gray
};
static const int kPopMark[kNPop] = { 20, 24, 23, 25 };

// Max scatter-plot entries per population per energy (subsampling for file size)
static const int kMaxScatterPts = 3000;

// Sentinel for invalid timing
static const float kNoTimeSentinel = -1e5f;

// ---------------------------------------------------------------------------
// OpenNtuple — forward declaration (same as compareEnergies.C)
// ---------------------------------------------------------------------------
static TFile* OpenNtuple(int r, TTree*& tree)
{
    TString path = Form("Analysis/Output/%s/ntuple.root",
                        kRuns[r].label.Data());
    TFile* f = TFile::Open(path);
    if (!f || f->IsZombie()) {
        std::cerr << "investigatePbGlass: cannot open " << path << "\n";
        tree = nullptr;
        return nullptr;
    }
    tree = static_cast<TTree*>(f->Get("rad"));
    if (!tree) {
        std::cerr << "investigatePbGlass: tree 'rad' missing in " << path << "\n";
        f->Close(); delete f;
        tree = nullptr; return nullptr;
    }
    return f;
}

// ---------------------------------------------------------------------------
// ClassifyEvent — assign a population to one event.
//
// Returns kD if wc_ok is false (caller still needs to check wc_ok for the
// hit map; here it is used as a gate before calling this function).
// ---------------------------------------------------------------------------
static Pop ClassifyEvent(float sum_lg, float sum_pb, float r_beam)
{
    if (r_beam >= static_cast<float>(kFiducial_r_timing)) return kD;
    if (sum_lg <= kSumLG_centroid)                         return kC;
    if (sum_pb >= kPb_maxRatio * sum_lg)                   return kB;
    return kA;
}

// ---------------------------------------------------------------------------
// ComboTime — compute A²-weighted mean CFD-5% time for one event.
//
// Returns kNoTimeSentinel if fewer than two channels pass quality cuts.
// The returned value is in nanoseconds (= Δt relative to MCP).
// ---------------------------------------------------------------------------
static float ComboTime(const Float_t hg_peak[8], const Float_t hg_cfd[8])
{
    double sw = 0., swt = 0.;
    for (int ic = 0; ic < kNCap; ++ic) {
        if (hg_peak[ic] < kHG_minPeak)   continue;
        if (hg_cfd[ic]  < kNoTimeSentinel) continue;
        double w = static_cast<double>(hg_peak[ic]) * hg_peak[ic];
        sw  += w;
        swt += w * hg_cfd[ic];
    }
    if (sw <= 0.) return kNoTimeSentinel;
    return static_cast<float>(swt / sw);
}

// ---------------------------------------------------------------------------
// FitSigma — fit a Gaussian core to a vector of timing values.
// Returns σ [ps], or -1 on failure (< 50 entries).
// ---------------------------------------------------------------------------
static double FitSigma(const std::vector<float>& tv)
{
    if (static_cast<int>(tv.size()) < 50) return -1.;

    double mu = 0.;
    for (float v : tv) mu += v;
    mu /= tv.size();

    double rms = 0.;
    for (float v : tv) rms += (v - mu) * (v - mu);
    rms = std::sqrt(rms / tv.size());
    if (rms < 0.010) rms = 0.100;

    TH1F htmp("_ipb_tmp", "", 200, mu - 4.*rms, mu + 4.*rms);
    htmp.SetDirectory(nullptr);
    for (float v : tv) htmp.Fill(v);

    double m2, mE, sig, sigE;
    FitGaussCore(&htmp, 2.0, m2, mE, sig, sigE);
    return (sig > 0.) ? sig * 1000. : -1.;   // convert ns → ps
}

// ---------------------------------------------------------------------------
// PageTitle — centred title at top of canvas
// ---------------------------------------------------------------------------
static void PageTitle(const char* title)
{
    TLatex t; t.SetNDC(); t.SetTextSize(0.030); t.SetTextAlign(21);
    t.DrawLatex(0.50, 0.987, title);
}

// ===========================================================================
// Main
// ===========================================================================
void investigatePbGlass()
{
    TH1::AddDirectory(kFALSE);
    ApplyRADiCALStyle();   // RADiCALStyle.h via PlotUtils.h

    gSystem->mkdir("Analysis/Output/Summary", kTRUE);

    std::cout << "investigatePbGlass: classifying events into populations A/B/C/D\n";

    // ── Per-energy, per-population histograms ────────────────────────────────
    const double kHitXlo = -2., kHitXhi = 16.;
    const double kHitYlo = -4., kHitYhi = 14.;

    // Scatter (Page 1): TGraph per population per energy (subsampled)
    TGraph* gScat[6][kNPop] = {};

    // Hit maps (Page 2): populations A, B, C at 25 GeV (iRun=0) and 150 GeV (iRun=5)
    static const int kHitRuns[2] = { 0, 5 };   // 25 GeV, 150 GeV
    TH2F* hHit[2][kNPop] = {};                  // [iHitRun][pop]

    // Timing distributions (Page 3): A vs B per energy
    // Stored as vectors; histograms built after the loop
    std::vector<float> vTiming[6][2];  // [iRun][0=popA, 1=popB]

    // Containment ratio (Page 5): A vs B per energy
    TH1F* hRatio[6][2] = {};   // [iRun][0=popA, 1=popB]

    // ΣLG spectra (Page 6): A vs B per energy
    TH1F* hLG[6][2] = {};      // [iRun][0=popA, 1=popB]

    // σ_t for Page 4
    double sigA[6] = {}, sigB[6] = {}, sigAB[6] = {};
    long   nPopA[6] = {}, nPopB[6] = {}, nPopC[6] = {}, nPopD[6] = {};
    double xCen[6] = {}, yCen[6] = {};

    // =========================================================================
    // Per-energy event loop
    // =========================================================================
    for (int r = 0; r < kNRuns; ++r) {
        TTree* t = nullptr;
        TFile* fin = OpenNtuple(r, t);
        if (!fin) continue;

        double xc, yc, tcfd, trms;
        ScanRunCenters(t, xc, yc, tcfd, trms);
        xCen[r] = xc; yCen[r] = yc;

        double xmax_lg = 55. * kRuns[r].energy_GeV;
        double ymax_pb = 0.50 * xmax_lg;
        TString tag    = kRuns[r].label;

        // Book scatter graphs (we'll subsample below)
        for (int p = 0; p < kNPop; ++p)
            gScat[r][p] = new TGraph();

        // Book hit map histograms for the two representative energies
        for (int ih = 0; ih < 2; ++ih) {
            if (r == kHitRuns[ih]) {
                for (int p = 0; p < kNPop - 1; ++p) {  // skip D (off-axis)
                    // 1 mm/bin matches WC delay-line resolution (kWC_resBin)
                    int nHX = static_cast<int>(std::round((kHitXhi - kHitXlo) / kWC_resBin));
                    int nHY = static_cast<int>(std::round((kHitYhi - kHitYlo) / kWC_resBin));
                    hHit[ih][p] = new TH2F(
                        Form("hHit_%s_p%d", tag.Data(), p), "",
                        nHX, kHitXlo, kHitXhi, nHY, kHitYlo, kHitYhi);
                    hHit[ih][p]->SetDirectory(nullptr);
                }
            }
        }

        // Book ratio + LG histograms (A and B only)
        for (int ab = 0; ab < 2; ++ab) {
            hRatio[r][ab] = new TH1F(Form("hRatio_%s_p%d", tag.Data(), ab),
                                     "", 100, 0., 1.2);
            hRatio[r][ab]->SetDirectory(nullptr);
            hLG[r][ab]    = new TH1F(Form("hLG_%s_p%d",    tag.Data(), ab),
                                     "", 150, 0., xmax_lg);
            hLG[r][ab]->SetDirectory(nullptr);
        }

        // Subsample counters (scatter display only)
        long nDrawn[kNPop] = {};

        // Branch addresses
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

        Long64_t nEv = t->GetEntries();
        for (Long64_t ev = 0; ev < nEv; ++ev) {
            t->GetEntry(ev);
            if (!wc_ok) continue;

            float dx     = x_trk - static_cast<float>(xc);
            float dy     = y_trk - static_cast<float>(yc);
            float r_beam = std::sqrt(dx*dx + dy*dy);

            Pop pop = ClassifyEvent(sum_lg, sum_pb, r_beam);

            // Count populations
            switch (pop) {
                case kA: ++nPopA[r]; break;
                case kB: ++nPopB[r]; break;
                case kC: ++nPopC[r]; break;
                case kD: ++nPopD[r]; break;
            }

            // Scatter (subsample to kMaxScatterPts per population)
            if (nDrawn[pop] < kMaxScatterPts) {
                gScat[r][pop]->SetPoint(gScat[r][pop]->GetN(), sum_lg, sum_pb);
                ++nDrawn[pop];
            }

            // Hit maps for representative energies
            for (int ih = 0; ih < 2; ++ih) {
                if (r == kHitRuns[ih] && pop != kD)
                    hHit[ih][pop]->Fill(x_trk, y_trk);
            }

            // ---- Timing-only section: require MCP quality ----------------
            bool mcp_ok = (mcp_peak > kMCP1_minPeak && mcp_peak < kMCP1_maxPeak);
            if (!mcp_ok) continue;

            // Ratio and LG spectra for A and B (timing fiducial)
            if (pop == kA || pop == kB) {
                int ab = (pop == kA) ? 0 : 1;
                if (sum_lg > kSumLG_centroid) {
                    hRatio[r][ab]->Fill(sum_pb / sum_lg);
                    hLG[r][ab]->Fill(sum_lg);
                }
            }

            // Timing: compute combo time and store per population
            if (pop == kA || pop == kB) {
                float t_combo = ComboTime(hg_peak, hg_cfd);
                if (t_combo > kNoTimeSentinel) {
                    int ab = (pop == kA) ? 0 : 1;
                    vTiming[r][ab].push_back(t_combo);
                }
            }
        }

        t->ResetBranchAddresses();
        fin->Close();
        delete fin;

        long nTotal = nPopA[r] + nPopB[r] + nPopC[r] + nPopD[r];
        std::cout << "  " << kRuns[r].label
                  << Form(":  A=%ld (%.1f%%)  B=%ld (%.1f%%)"
                          "  C=%ld (%.1f%%)  D=%ld (%.1f%%)\n",
                          nPopA[r], 100.*nPopA[r]/nTotal,
                          nPopB[r], 100.*nPopB[r]/nTotal,
                          nPopC[r], 100.*nPopC[r]/nTotal,
                          nPopD[r], 100.*nPopD[r]/nTotal);
    }

    // Compute σ_t per population per energy
    for (int r = 0; r < kNRuns; ++r) {
        sigA[r]  = FitSigma(vTiming[r][0]);
        sigB[r]  = FitSigma(vTiming[r][1]);

        // Combined A+B
        std::vector<float> combined;
        combined.reserve(vTiming[r][0].size() + vTiming[r][1].size());
        for (float v : vTiming[r][0]) combined.push_back(v);
        for (float v : vTiming[r][1]) combined.push_back(v);
        sigAB[r] = FitSigma(combined);
    }

    // =========================================================================
    // Output PDF
    // =========================================================================
    const char* kSumDir = "Analysis/Output/Summary/";
    TString outPDF = Form("%spbglass_investigation.pdf", kSumDir);
    TCanvas c("c_ipb", "", 1200, 800);

    // ── Page 1: Scatter coloured by population (2×3 grid) ────────────────────
    c.Divide(3, 2, 0.003, 0.035);  // 3.5% top/bottom margin reserves space for PageTitle
    for (int r = 0; r < kNRuns; ++r) {
        c.cd(r + 1);
        gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
        gPad->SetRightMargin(0.04); gPad->SetTopMargin(0.10);
        gPad->SetTickx(1); gPad->SetTicky(1);

        double xmax_lg = 55. * kRuns[r].energy_GeV;
        double ymax_pb = 0.50 * xmax_lg;

        // Draw a blank frame so axes are set before the scatter points
        TH2F hFrame(Form("hFrame_%d", r), "",
                    1, 0., xmax_lg, 1, 0., ymax_pb);
        hFrame.GetXaxis()->SetTitle("#Sigma LG (mV)");
        hFrame.GetYaxis()->SetTitle("#Sigma PbGlass (mV)");
        hFrame.GetXaxis()->SetTitleSize(0.07); hFrame.GetXaxis()->SetLabelSize(0.06);
        hFrame.GetYaxis()->SetTitleSize(0.07); hFrame.GetYaxis()->SetLabelSize(0.06);
        hFrame.GetYaxis()->SetTitleOffset(0.90);
        hFrame.DrawCopy();

        // Containment cut line
        TLine* lcut = new TLine(0., 0., xmax_lg, kPb_maxRatio * xmax_lg);
        lcut->SetLineColor(kBlack); lcut->SetLineStyle(2); lcut->SetLineWidth(2);
        lcut->Draw("SAME");

        // Draw populations: D first (background), then C, B, A (foreground)
        static const Pop kDrawOrder[kNPop] = { kD, kC, kB, kA };
        for (int iord = 0; iord < kNPop; ++iord) {
            Pop p = kDrawOrder[iord];
            TGraph* g = gScat[r][p];
            if (!g || g->GetN() == 0) continue;
            g->SetMarkerColor(kPopCol[p]);
            g->SetMarkerStyle(7);   // small dot
            g->SetMarkerSize(0.4);
            g->Draw("P SAME");
        }

        // Legend (only in pad 1)
        if (r == 0) {
            TLegend* leg = new TLegend(0.63, 0.63, 0.97, 0.88);
            leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.065);
            for (int p = 0; p < kNPop; ++p)
                leg->AddEntry(gScat[r][p], kPopName[p], "p");
            leg->Draw();
        }

        TLatex lab; lab.SetNDC(); lab.SetTextSize(0.08); lab.SetTextFont(72);
        lab.DrawLatex(0.17, 0.85, Form("%.0f GeV", kRuns[r].energy_GeV));
        lab.SetTextFont(42); lab.SetTextSize(0.072);
        long nTot = nPopA[r]+nPopB[r]+nPopC[r]+nPopD[r];
        if (nTot > 0)
            lab.DrawLatex(0.17, 0.74,
                Form("B: %.1f%%", 100.*nPopB[r]/nTot));
    }
    c.cd(0);
    PageTitle("#Sigma_{PbGlass} vs #Sigma_{LG}  coloured by population  (WC-OK, dashed = 30% cut)");
    c.Print(outPDF + "(");

    // ── Page 2: Hit maps for populations A, B, C — 25 GeV and 150 GeV ───────
    // 3 columns (A, B, C), 2 rows (25 GeV, 150 GeV)
    c.Clear();
    c.Divide(3, 2, 0.003, 0.035);  // 3.5% top/bottom margin reserves space for PageTitle
    gStyle->SetPalette(kRust); TColor::InvertPalette();
    int padIdx = 1;
    for (int ih = 0; ih < 2; ++ih) {
        for (int p = 0; p < kNPop - 1; ++p) {  // A, B, C (skip D)
            c.cd(padIdx++);
            gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12);
            gPad->SetRightMargin(0.14); gPad->SetTopMargin(0.15);
            gPad->SetTickx(1); gPad->SetTicky(1);
            if (!hHit[ih][p]) continue;

            hHit[ih][p]->GetXaxis()->SetTitle("x_{trk} (mm)");
            hHit[ih][p]->GetYaxis()->SetTitle("y_{trk} (mm)");
            hHit[ih][p]->GetXaxis()->SetTitleSize(0.065);
            hHit[ih][p]->GetXaxis()->SetLabelSize(0.060);
            hHit[ih][p]->GetYaxis()->SetTitleSize(0.065);
            hHit[ih][p]->GetYaxis()->SetLabelSize(0.060);
            hHit[ih][p]->GetYaxis()->SetTitleOffset(0.90);
            hHit[ih][p]->Draw("COLZ");

            // Timing fiducial circle
            int er = kHitRuns[ih];
            TEllipse* ell = new TEllipse(xCen[er], yCen[er],
                                         kFiducial_r_timing, kFiducial_r_timing);
            ell->SetLineColor(kOrange+1); ell->SetLineWidth(2);
            ell->SetFillStyle(0); ell->Draw("SAME");

            TLatex lab; lab.SetNDC(); lab.SetTextSize(0.08); lab.SetTextFont(72);
            lab.SetTextColor(kPopCol[p]);
            lab.DrawLatex(0.17, 0.90, kPopName[p]);
            lab.SetTextFont(42); lab.SetTextColor(kBlack); lab.SetTextSize(0.065);
            lab.DrawLatex(0.17, 0.80,
                Form("%.0f GeV", kRuns[kHitRuns[ih]].energy_GeV));
        }
    }
    c.cd(0);
    PageTitle("Beam hit maps per population (top: 25 GeV  |  bottom: 150 GeV)");
    c.Print(outPDF);

    // ── Page 3: Timing distributions A vs B per energy (2×3 grid) ────────────
    // Each panel: A (green, solid) and B (red, solid) overlaid.
    // Normalised to unit area so shapes are directly comparable.
    TH1F* hTiming[6][2] = {};
    for (int r = 0; r < kNRuns; ++r) {
        for (int ab = 0; ab < 2; ++ab) {
            if (vTiming[r][ab].size() < 10) continue;
            // Auto-range on mean ± 4 RMS
            double mu = 0., rms2 = 0.;
            for (float v : vTiming[r][ab]) mu += v;
            mu /= vTiming[r][ab].size();
            for (float v : vTiming[r][ab]) rms2 += (v-mu)*(v-mu);
            double rms = std::sqrt(rms2 / vTiming[r][ab].size());
            if (rms < 0.010) rms = 0.100;

            hTiming[r][ab] = new TH1F(
                Form("hTiming_%s_p%d", kRuns[r].label.Data(), ab), "",
                120, mu - 4.*rms, mu + 4.*rms);
            hTiming[r][ab]->SetDirectory(nullptr);
            for (float v : vTiming[r][ab]) hTiming[r][ab]->Fill(v);
        }
    }

    c.Clear();
    c.Divide(3, 2, 0.003, 0.035);  // 3.5% top/bottom margin reserves space for PageTitle
    for (int r = 0; r < kNRuns; ++r) {
        c.cd(r + 1);
        gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);
        gPad->SetRightMargin(0.04); gPad->SetTopMargin(0.10);
        gPad->SetTickx(1); gPad->SetTicky(1);

        double yMaxT = 0.;
        for (int ab = 0; ab < 2; ++ab) {
            if (!hTiming[r][ab]) continue;
            double integral = hTiming[r][ab]->Integral();
            if (integral > 0.) hTiming[r][ab]->Scale(1.0 / integral);
            yMaxT = std::max(yMaxT, hTiming[r][ab]->GetMaximum());
        }
        if (yMaxT <= 0.) { c.cd(r+1); continue; }

        bool first = true;
        for (int ab = 0; ab < 2; ++ab) {
            if (!hTiming[r][ab]) continue;
            int p = (ab == 0) ? kA : kB;
            hTiming[r][ab]->SetLineColor(kPopCol[p]);
            hTiming[r][ab]->SetLineWidth(2);
            hTiming[r][ab]->GetXaxis()->SetTitle("t_{combo} (ns)");
            hTiming[r][ab]->GetXaxis()->SetTitleSize(0.070);
            hTiming[r][ab]->GetXaxis()->SetLabelSize(0.060);
            hTiming[r][ab]->GetYaxis()->SetTitleSize(0.070);
            hTiming[r][ab]->GetYaxis()->SetLabelSize(0.060);
            hTiming[r][ab]->GetYaxis()->SetTitle("Events (norm.)");
            hTiming[r][ab]->GetYaxis()->SetRangeUser(0., 1.20 * yMaxT);
            if (first) { hTiming[r][ab]->Draw("HIST"); first = false; }
            else          hTiming[r][ab]->Draw("HIST SAME");
        }

        TLatex lab; lab.SetNDC(); lab.SetTextSize(0.070); lab.SetTextFont(72);
        lab.DrawLatex(0.18, 0.86, Form("%.0f GeV", kRuns[r].energy_GeV));
        lab.SetTextFont(42); lab.SetTextSize(0.062);
        if (sigA[r] > 0.)
            lab.SetTextColor(kPopCol[kA]),
            lab.DrawLatex(0.18, 0.76, Form("#sigma_{A} = %.0f ps", sigA[r]));
        if (sigB[r] > 0.)
            lab.SetTextColor(kPopCol[kB]),
            lab.DrawLatex(0.18, 0.66, Form("#sigma_{B} = %.0f ps", sigB[r]));
    }
    c.cd(0);
    PageTitle("A²-wgt combo timing: Pop. A (green) vs Pop. B (red)  |  normalised to unit area");
    c.Print(outPDF);

    // ── Page 4: σ_t for A / B / (A+B) vs energy ─────────────────────────────
    c.Clear(); c.cd(); StylePad();

    TGraphErrors* gSigA  = new TGraphErrors();
    TGraphErrors* gSigB  = new TGraphErrors();
    TGraphErrors* gSigAB = new TGraphErrors();

    for (int r = 0; r < kNRuns; ++r) {
        double E = kRuns[r].energy_GeV;
        if (sigA[r]  > 0.) { int n = gSigA->GetN();  gSigA->SetPoint(n, E, sigA[r]);  }
        if (sigB[r]  > 0.) { int n = gSigB->GetN();  gSigB->SetPoint(n, E, sigB[r]);  }
        if (sigAB[r] > 0.) { int n = gSigAB->GetN(); gSigAB->SetPoint(n, E, sigAB[r]); }
    }

    auto styleG = [](TGraphErrors* g, int col, int mark, int width) {
        g->SetLineColor(col); g->SetMarkerColor(col);
        g->SetMarkerStyle(mark); g->SetMarkerSize(1.3);
        g->SetLineWidth(width);
    };
    styleG(gSigA,  kPopCol[kA], 20, 2);
    styleG(gSigB,  kPopCol[kB], 21, 2);
    styleG(gSigAB, kBlack,      22, 1);

    // Use only EM (A) and combined (A+B) to set y-range.
    // Population B Gaussian-core σ_t is unreliable for the non-Gaussian
    // hadronic timing distribution and would push the scale to thousands of ps,
    // making the EM population (our signal) invisible.
    double yMax4 = 0.;
    for (int r = 0; r < kNRuns; ++r) {
        if (sigA[r]  > 0.) yMax4 = std::max(yMax4, sigA[r]);
        if (sigAB[r] > 0.) yMax4 = std::max(yMax4, sigAB[r]);
    }
    if (yMax4 < 10.) yMax4 = 250.;

    TH1F* frame4 = static_cast<TH1F*>(
        c.DrawFrame(15., 50., 165., std::max(300., 1.50 * yMax4),
                    ";Beam energy (GeV);#sigma_{t} (ps)  [A^{2}-wgt combo]"));
    frame4->GetXaxis()->SetTitleSize(0.050);
    frame4->GetYaxis()->SetTitleSize(0.050);
    frame4->GetYaxis()->SetTitleOffset(1.30);

    if (gSigAB->GetN() > 0) gSigAB->Draw("PL SAME");
    if (gSigA->GetN()  > 0) gSigA->Draw("PL SAME");
    if (gSigB->GetN()  > 0) gSigB->Draw("PL SAME");  // clipped at axis top if off-scale

    TLegend leg4(0.65, 0.67, 0.93, 0.88);
    leg4.SetBorderSize(0); leg4.SetFillStyle(0); leg4.SetTextSize(0.042);
    if (gSigA->GetN()  > 0) leg4.AddEntry(gSigA,  "Pop. A #minus Good EM shower",        "lp");
    if (gSigB->GetN()  > 0) leg4.AddEntry(gSigB,  "Pop. B #minus Punch-through/hadronic","lp");
    if (gSigAB->GetN() > 0) leg4.AddEntry(gSigAB, "A + B combined",                      "lp");
    leg4.Draw();

    // Annotation: hadronic fraction at highest energy + note about Pop. B
    {
        int best = -1;
        for (int r = kNRuns-1; r >= 0; --r)
            if (sigA[r] > 0. && sigB[r] > 0.) { best = r; break; }
        if (best >= 0) {
            TLatex ann; ann.SetNDC(); ann.SetTextSize(0.038);
            ann.SetTextColor(kGray+1);
            ann.DrawLatex(0.18, 0.24,
                Form("containment cut removes %.0f%% of timing-fid. events",
                     100. * nPopB[best] /
                     std::max(1L, nPopA[best] + nPopB[best])));
        }
        TLatex note; note.SetNDC(); note.SetTextSize(0.033); note.SetTextColor(kRed+1);
        note.DrawLatex(0.16, 0.16,
            "Pop. B: hadronic timing non-Gaussian #minus Gaussian #sigma may exceed axis range");
    }
    PageTitle("#sigma_{t} per population vs beam energy  |  does the containment cut help?");
    c.Print(outPDF);

    // ── Page 5: Containment ratio A vs B per energy (2×3 grid) ───────────────
    c.Clear();
    c.Divide(3, 2, 0.003, 0.035);  // 3.5% top/bottom margin reserves space for PageTitle
    for (int r = 0; r < kNRuns; ++r) {
        c.cd(r + 1);
        gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);
        gPad->SetRightMargin(0.04); gPad->SetTopMargin(0.10);
        gPad->SetTickx(1); gPad->SetTicky(1);

        double yMaxR = 0.;
        for (int ab = 0; ab < 2; ++ab) {
            if (!hRatio[r][ab]) continue;
            double integral = hRatio[r][ab]->Integral();
            if (integral > 0.) hRatio[r][ab]->Scale(1.0 / integral);
            yMaxR = std::max(yMaxR, hRatio[r][ab]->GetMaximum());
        }
        if (yMaxR <= 0.) continue;

        bool first = true;
        for (int ab = 0; ab < 2; ++ab) {
            if (!hRatio[r][ab]) continue;
            int p = (ab == 0) ? kA : kB;
            hRatio[r][ab]->SetLineColor(kPopCol[p]);
            hRatio[r][ab]->SetLineWidth(2);
            hRatio[r][ab]->GetXaxis()->SetTitle("#Sigma_{PbGlass} / #Sigma_{LG}");
            hRatio[r][ab]->GetXaxis()->SetTitleSize(0.070);
            hRatio[r][ab]->GetXaxis()->SetLabelSize(0.060);
            hRatio[r][ab]->GetYaxis()->SetTitleSize(0.070);
            hRatio[r][ab]->GetYaxis()->SetLabelSize(0.060);
            hRatio[r][ab]->GetYaxis()->SetTitle("Events (norm.)");
            hRatio[r][ab]->GetYaxis()->SetRangeUser(0., 1.25 * yMaxR);
            if (first) { hRatio[r][ab]->Draw("HIST"); first = false; }
            else          hRatio[r][ab]->Draw("HIST SAME");
        }
        TLine* lcut = new TLine(kPb_maxRatio, 0., kPb_maxRatio, 1.15*yMaxR);
        lcut->SetLineColor(kBlack); lcut->SetLineStyle(2); lcut->SetLineWidth(2);
        lcut->Draw("SAME");   // pad takes ownership; do not delete

        TLatex lab; lab.SetNDC(); lab.SetTextFont(72); lab.SetTextSize(0.070);
        lab.DrawLatex(0.57, 0.86, Form("%.0f GeV", kRuns[r].energy_GeV));
    }
    c.cd(0);
    PageTitle("#Sigma_{PbGlass}/#Sigma_{LG} ratio for Pop. A (green) vs Pop. B (red)  |  dashed: 30% cut");
    c.Print(outPDF);

    // ── Page 6: ΣLG spectra A vs B per energy (2×3 grid) ─────────────────────
    c.Clear();
    c.Divide(3, 2, 0.003, 0.035);  // 3.5% top/bottom margin reserves space for PageTitle
    for (int r = 0; r < kNRuns; ++r) {
        c.cd(r + 1);
        gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);
        gPad->SetRightMargin(0.04); gPad->SetTopMargin(0.10);
        gPad->SetTickx(1); gPad->SetTicky(1);

        double yMaxL = 0.;
        for (int ab = 0; ab < 2; ++ab) {
            if (!hLG[r][ab]) continue;
            double integral = hLG[r][ab]->Integral();
            if (integral > 0.) hLG[r][ab]->Scale(1.0 / integral);
            yMaxL = std::max(yMaxL, hLG[r][ab]->GetMaximum());
        }
        if (yMaxL <= 0.) continue;

        bool first = true;
        for (int ab = 0; ab < 2; ++ab) {
            if (!hLG[r][ab]) continue;
            int p = (ab == 0) ? kA : kB;
            hLG[r][ab]->SetLineColor(kPopCol[p]);
            hLG[r][ab]->SetLineWidth(2);
            hLG[r][ab]->GetXaxis()->SetTitle("#Sigma LG (mV)");
            hLG[r][ab]->GetXaxis()->SetTitleSize(0.070);
            hLG[r][ab]->GetXaxis()->SetLabelSize(0.060);
            hLG[r][ab]->GetYaxis()->SetTitleSize(0.070);
            hLG[r][ab]->GetYaxis()->SetLabelSize(0.060);
            hLG[r][ab]->GetYaxis()->SetTitle("Events (norm.)");
            hLG[r][ab]->GetYaxis()->SetRangeUser(0., 1.25 * yMaxL);
            if (first) { hLG[r][ab]->Draw("HIST"); first = false; }
            else          hLG[r][ab]->Draw("HIST SAME");
        }
        TLatex lab; lab.SetNDC(); lab.SetTextFont(72); lab.SetTextSize(0.070);
        lab.DrawLatex(0.57, 0.86, Form("%.0f GeV", kRuns[r].energy_GeV));
    }
    c.cd(0);
    PageTitle("#Sigma_{LG} spectra for Pop. A (green) vs Pop. B (red)  |  normalised to unit area");
    c.Print(outPDF + ")");

    // =========================================================================
    // Persist headline numbers for the data-driven report (results.json harvest)
    //   punch-through fraction = Pop.B / (Pop.A + Pop.B), per energy
    // =========================================================================
    {
        TGraph gPunch;  // x = beam energy [GeV], y = punch-through fraction [%]
        for (int r = 0; r < kNRuns; ++r) {
            long denom = std::max(1L, nPopA[r] + nPopB[r]);
            gPunch.SetPoint(gPunch.GetN(), kRuns[r].energy_GeV,
                            100. * static_cast<double>(nPopB[r]) / denom);
        }
        TFile fout(Form("%spbglass_investigation.root", kSumDir), "RECREATE");
        gPunch.Write("gPunchThroughFrac");
        fout.Close();
        std::cout << "investigatePbGlass: wrote pbglass_investigation.root "
                     "(gPunchThroughFrac)\n";
    }

    // =========================================================================
    // Cleanup
    // =========================================================================
    for (int r = 0; r < kNRuns; ++r) {
        for (int p = 0; p < kNPop; ++p) delete gScat[r][p];
        for (int ab = 0; ab < 2; ++ab) {
            delete hRatio[r][ab];
            delete hLG[r][ab];
            delete hTiming[r][ab];
        }
    }
    for (int ih = 0; ih < 2; ++ih)
        for (int p = 0; p < kNPop - 1; ++p)
            delete hHit[ih][p];
    delete gSigA; delete gSigB; delete gSigAB;

    std::cout << "\ninvestigatePbGlass: done -> " << outPDF << "\n";
}

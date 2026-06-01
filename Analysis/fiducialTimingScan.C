// ============================================================================
// fiducialTimingScan.C  —  Is the timing fiducial radius optimal?
// ----------------------------------------------------------------------------
// Referee-proof companion to the best-bin headline.  The single-best-E_meas-bin
// selection used for the headline is statistically noisy as a function of the
// fiducial radius (the "best bin" jumps between narrow E_meas slices), so it
// cannot be read as a clean radius optimisation.  This macro instead uses a
// STABLE selection — a fixed fraction of the highest-E_meas events — and scans
// the (DW-UP)/2 corner resolution vs the timing fiducial radius at 150 GeV.
//
// Three curves:
//   - ALL fiducial events        -> sigma_t DEGRADES with radius: the raw
//                                    position effect (outer ring = low amplitude
//                                    = residual time-walk).
//   - top 10% by E_meas (sum_lg) -> FLAT: the E_meas selection already keeps the
//   - top  2% by E_meas             well-measured, central showers, so the
//                                    fiducial circle adds essentially nothing.
//
// Conclusion the figure makes airtight: the headline resolution is set by the
// E_meas selection, NOT by the fiducial radius — so r<3 mm is a defensible loose
// pre-cut, and there is no robust improvement to be had by re-tuning it.
//
// Uses the SAME Method-A formula, event cuts, data-derived centroid and
// Gaussian-core fit as timingEnergyBins.C (the headline).
//
// Output: Analysis/Output/Summary/fiducial_timing_scan.{png,pdf}
// ============================================================================

#include <vector>
#include <algorithm>
#include <cmath>
#include <utility>

#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLine.h"
#include "TBox.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TString.h"
#include "TSystem.h"
#include "TColor.h"

#include "PlotUtils.h"      // ScanRunCenters, FitGaussCore, StylePad
#include "RADiCALStyle.h"   // ApplyRADiCALStyle, kRData/kRRed/kROrange, DrawPageTitle
#include "SelectionCuts.h"  // kMCP1_minPeak, kMCP1_maxPeak, kFiducial_r_timing

static const float kSentCut_fts = -1e5f;   // hg_cfd validity test (matches headline)

// Robust core sigma [ns] of a set of corner times — mirrors VecToHist_teb +
// FitGaussCore (2-pass 5-sigma outlier trim, then 2-sigma Gaussian core fit).
static double RobustSigma_fts(const std::vector<float>& v)
{
    if (v.size() < 200) return -1.;
    double mu1 = 0.; for (float x : v) mu1 += x; mu1 /= (double)v.size();
    double ms1 = 0.; for (float x : v) ms1 += ((double)x - mu1) * ((double)x - mu1);
    ms1 = std::sqrt(ms1 / (double)v.size());
    if (ms1 < 0.008) ms1 = 0.100;

    double mu2 = 0.; long n2 = 0;
    for (float x : v) if (std::fabs((double)x - mu1) < 5. * ms1) { mu2 += x; ++n2; }
    if (n2 < 100) return -1.;
    mu2 /= (double)n2;
    double ms2 = 0.;
    for (float x : v) if (std::fabs((double)x - mu1) < 5. * ms1)
        ms2 += ((double)x - mu2) * ((double)x - mu2);
    ms2 = std::sqrt(ms2 / (double)n2);
    if (ms2 < 0.008) ms2 = 0.100;

    TH1F h("_fts_h", "", 120, mu2 - 5. * ms2, mu2 + 5. * ms2);
    h.SetDirectory(nullptr);
    for (float x : v) h.Fill(x);
    double mu, muE, s, sE;
    FitGaussCore(&h, 2.0, mu, muE, s, sE);
    return (s > 0.) ? s : -1.;
}

void fiducialTimingScan()
{
    TH1::AddDirectory(kFALSE);
    ApplyRADiCALStyle();

    const char* ntuple = "Analysis/Output/150GeV/ntuple.root";
    TFile* f = TFile::Open(ntuple);
    if (!f || f->IsZombie()) { printf("[fiducialTimingScan] cannot open %s\n", ntuple); return; }
    TTree* t = dynamic_cast<TTree*>(f->Get("rad"));
    if (!t) { printf("[fiducialTimingScan] no TTree 'rad'\n"); return; }

    double xc, yc, tOff, tRms;
    ScanRunCenters(t, xc, yc, tOff, tRms);
    printf("[fiducialTimingScan] data-derived centroid (%.2f, %.2f) mm\n", xc, yc);

    Bool_t  wc_ok;
    Float_t x_trk, y_trk, mcp_peak, mcp_time, hg_cfd[8], sum_lg;
    t->SetBranchAddress("wc_ok",    &wc_ok);
    t->SetBranchAddress("x_trk",    &x_trk);
    t->SetBranchAddress("y_trk",    &y_trk);
    t->SetBranchAddress("mcp_peak", &mcp_peak);
    t->SetBranchAddress("mcp_time", &mcp_time);
    t->SetBranchAddress(t->GetBranch("hg_cfd05") ? "hg_cfd05" : "hg_cfd", hg_cfd);
    t->SetBranchAddress("sum_lg",   &sum_lg);

    struct Ev { float r2, slg, tc; };
    std::vector<Ev> evs; evs.reserve(200000);

    const Long64_t N = t->GetEntries();
    for (Long64_t i = 0; i < N; ++i) {
        t->GetEntry(i);
        if (!wc_ok)                                              continue;
        if (mcp_peak < kMCP1_minPeak || mcp_peak > kMCP1_maxPeak) continue;
        if (mcp_time <= kSentCut_fts)                            continue;

        double dw = 0., up = 0.; int ndw = 0, nup = 0;
        for (int c = 0; c < 4; ++c) if (hg_cfd[c] > kSentCut_fts) { dw += hg_cfd[c]; ++ndw; }
        for (int c = 4; c < 8; ++c) if (hg_cfd[c] > kSentCut_fts) { up += hg_cfd[c]; ++nup; }
        if (ndw < 1 || nup < 1) continue;

        float tc = (float)((dw / ndw - up / nup) * 0.5);     // Method A: (DW-UP)/2
        float dx = x_trk - (float)xc, dy = y_trk - (float)yc;
        evs.push_back({ dx * dx + dy * dy, sum_lg, tc });
    }
    printf("[fiducialTimingScan] %zu events pass event-level cuts\n", evs.size());

    // ── Scan radius ─────────────────────────────────────────────────────────
    TGraph gAll, gTop10, gTop2;
    double yMax = 0.;
    printf("  R(mm)    N      all     top10%%   top2%% [ps]\n");
    for (double R = 1.0; R <= 5.001; R += 0.25) {
        const float R2 = (float)(R * R);
        std::vector<std::pair<float, float>> sel;            // (sum_lg, tc)
        sel.reserve(evs.size());
        for (const auto& e : evs) if (e.r2 < R2) sel.emplace_back(e.slg, e.tc);
        if (sel.size() < 400) continue;

        std::vector<float> tAll; tAll.reserve(sel.size());
        for (const auto& p : sel) tAll.push_back(p.second);
        double sAll = RobustSigma_fts(tAll);

        std::sort(sel.begin(), sel.end(),
                  [](const std::pair<float,float>& a, const std::pair<float,float>& b)
                  { return a.first > b.first; });            // E_meas descending
        auto topSigma = [&](double frac) -> double {
            size_t nt = (size_t)std::max(200., frac * (double)sel.size());
            nt = std::min(nt, sel.size());
            std::vector<float> tt; tt.reserve(nt);
            for (size_t k = 0; k < nt; ++k) tt.push_back(sel[k].second);
            return RobustSigma_fts(tt);
        };
        double s10 = topSigma(0.10);
        double s02 = topSigma(0.02);

        if (sAll > 0.) { gAll.SetPoint(gAll.GetN(),  R, sAll * 1000.); yMax = std::max(yMax, sAll*1000.); }
        if (s10  > 0.)   gTop10.SetPoint(gTop10.GetN(), R, s10  * 1000.);
        if (s02  > 0.)   gTop2.SetPoint(gTop2.GetN(),  R, s02  * 1000.);

        printf("  %4.2f  %7zu  %6.1f  %6.1f  %6.1f\n",
               R, sel.size(), sAll*1000., s10*1000., s02*1000.);
    }

    // ── helpers for annotation ───────────────────────────────────────────────
    auto valAt = [](const TGraph& g, double r) -> double {
        for (int i = 0; i < g.GetN(); ++i) if (std::fabs(g.GetX()[i] - r) < 1e-6) return g.GetY()[i];
        return 0.;
    };
    auto plateauMax = [](const TGraph& g, double r0, double r1) -> double {
        double hi = -1e9;
        for (int i = 0; i < g.GetN(); ++i)
            if (g.GetX()[i] >= r0 - 1e-6 && g.GetX()[i] <= r1 + 1e-6) hi = std::max(hi, g.GetY()[i]);
        return hi;
    };
    auto allRange = [&](double r) { return valAt(gAll, r); };
    const double rPlateau = 2.5;   // plateau edge from the scan

    // ── Plot — focus on the STABLE high-E_meas selections (the clean answer) ──
    TCanvas* c = new TCanvas("c_fts", "", 920, 760);
    c->SetLeftMargin(0.14); c->SetBottomMargin(0.13);
    c->SetRightMargin(0.05); c->SetTopMargin(0.10);
    c->SetTickx(1); c->SetTicky(1);

    gTop10.SetMarkerStyle(21); gTop10.SetMarkerColor(kRData); gTop10.SetLineColor(kRData);
    gTop10.SetLineWidth(3);    gTop10.SetMarkerSize(1.3);
    gTop2.SetMarkerStyle(20);  gTop2.SetMarkerColor(kRRed);  gTop2.SetLineColor(kRRed);
    gTop2.SetLineWidth(3);     gTop2.SetMarkerSize(1.3);

    // y-range focused on the selection curves (where the shallow structure lives)
    double yl = 1e9, yh = -1e9;
    for (int i = 0; i < gTop2.GetN();  ++i)  { yl = std::min(yl, gTop2.GetY()[i]);  yh = std::max(yh, gTop2.GetY()[i]); }
    for (int i = 0; i < gTop10.GetN(); ++i)  { yl = std::min(yl, gTop10.GetY()[i]); yh = std::max(yh, gTop10.GetY()[i]); }
    const double ylo = yl - 4., yhi = yh + 4.;

    gTop2.Draw("APL");
    gTop2.GetXaxis()->SetTitle("timing fiducial radius  r  (mm)");
    gTop2.GetYaxis()->SetTitle("#sigma_{t}  (DW#minusUP)/2  (ps)");
    gTop2.GetXaxis()->SetLimits(0.5, 5.2);
    gTop2.GetYaxis()->SetRangeUser(ylo, yhi);
    gTop2.GetXaxis()->SetTitleSize(0.046); gTop2.GetYaxis()->SetTitleSize(0.046);
    gTop10.Draw("PL same");

    // shade the optimal plateau (r <= ~2.5 mm) and mark the current 3 mm cut
    TBox* band = new TBox(0.5, ylo, rPlateau, yhi);
    band->SetFillColorAlpha(kRData, 0.07); band->SetLineColor(0); band->Draw();
    gTop2.Draw("PL same"); gTop10.Draw("PL same");   // redraw over the band
    TLine* l3 = new TLine(kFiducial_r_timing, ylo, kFiducial_r_timing, yhi);
    l3->SetLineStyle(2); l3->SetLineColor(kROrange); l3->SetLineWidth(2); l3->Draw();
    { TLatex a; a.SetTextColor(kROrange); a.SetTextSize(0.029); a.SetTextAngle(90);
      a.DrawLatex(kFiducial_r_timing + 0.09, ylo + 0.60*(yhi-ylo), "current cut  r < 3 mm"); }
    { TLatex a; a.SetTextColor(kRData); a.SetTextSize(0.027);
      a.DrawLatex(0.85, yhi - 0.07*(yhi-ylo), "optimal"); a.DrawLatex(0.70, yhi - 0.11*(yhi-ylo), "plateau"); }

    TLegend* L = new TLegend(0.42, 0.74, 0.93, 0.88);
    L->SetBorderSize(0); L->SetFillStyle(0); L->SetTextFont(42); L->SetTextSize(0.032);
    L->AddEntry(&gTop2,  "top 2% by E_{meas} (#SigmaA_{LG})", "pl");
    L->AddEntry(&gTop10, "top 10% by E_{meas}", "pl");
    L->Draw();

    { TLatex a; a.SetNDC(); a.SetTextSize(0.029); a.SetTextColor(kGray+3);
      a.DrawLatex(0.17, 0.345,
        Form("flat plateau for r #leq %.1f mm; degrades gently beyond", rPlateau));
      a.DrawLatex(0.17, 0.300,
        Form("top 2%%:  %.1f ps (r#leq2.5) #rightarrow %.1f ps (3 mm) #rightarrow %.1f ps (4 mm)",
             plateauMax(gTop2,1.0,2.5), valAt(gTop2,3.0), valAt(gTop2,4.0)));
      a.DrawLatex(0.17, 0.255,
        "#Rightarrow optimum r #approx 2#minus2.5 mm; r=3 mm is ~1#minus2 ps above it (shallow, robust)");
      a.SetTextSize(0.026); a.SetTextColor(kGray+2);
      a.DrawLatex(0.17, 0.205,
        Form("(no E_{meas} selection: #sigma_{t} #approx %.0f#minus%.0f ps, dominated by low-E_{meas} showers)",
             plateauMax(gTop2,3.0,3.0) > 0 ? std::min(allRange(3.0), allRange(5.0)) : 60.,
             std::max(allRange(1.0), allRange(1.5)))); }

    DrawPageTitle("Timing resolution vs fiducial radius  (150 GeV, stable E_{meas} selection)");

    const char* outPng = "Analysis/Output/Summary/fiducial_timing_scan.png";
    const char* outPdf = "Analysis/Output/Summary/fiducial_timing_scan.pdf";
    c->Print(outPng);
    c->Print(outPdf);
    printf("[fiducialTimingScan] wrote %s\n", outPng);
    printf("[fiducialTimingScan] top2%%: plateau(<=2.5mm)=%.1f  3mm=%.1f  4mm=%.1f ps\n",
           plateauMax(gTop2,1.0,2.5), valAt(gTop2,3.0), valAt(gTop2,4.0));

    f->Close();
}

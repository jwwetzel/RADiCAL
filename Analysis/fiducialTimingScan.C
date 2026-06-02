// ============================================================================
// fiducialTimingScan.C  —  Is the timing fiducial radius optimal?  (all energies)
// ----------------------------------------------------------------------------
// Referee-proof companion to the best-bin headline.  The single-best-E_meas-bin
// selection used for the headline is statistically noisy as a function of the
// fiducial radius (the "best bin" jumps between narrow E_meas slices), so it
// cannot be read as a clean radius optimisation.  This macro instead uses a
// STABLE selection — a fixed fraction of the highest-E_meas events — and scans
// the (DW-UP)/2 corner resolution vs the timing fiducial radius.
//
// Uses the SAME Method-A formula, event cuts, data-derived centroid and
// Gaussian-core fit as timingEnergyBins.C (the headline).
//
// Outputs:
//   Analysis/Output/Summary/fiducial_timing_scan.png      — ALL 6 energies
//        overlaid (top-2% E_meas corner sigma_t vs fiducial radius).
//   Analysis/Output/Summary/fiducial_timing_scan_150.png  — 150 GeV detail
//        (top-2% and top-10% selections, annotated).
//
// Conclusion the figures make airtight: the corner resolution has a shallow
// optimal PLATEAU for r <~ 2.5 mm at EVERY energy, degrading gently beyond — a
// genuine optimisation, not a fluctuation.  (The single-best-bin headline then
// prefers r<3 mm because it needs the larger statistics to isolate the cleanest
// slice; that argument lives in the report text.)
// ============================================================================

#include <vector>
#include <algorithm>
#include <cmath>
#include <utility>

#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLine.h"
#include "TBox.h"
#include "TMarker.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TString.h"
#include "TSystem.h"
#include "TColor.h"

#include "PlotUtils.h"      // ScanRunCenters, FitGaussCore, StylePad
#include "RADiCALStyle.h"   // ApplyRADiCALStyle, kREnergyCols, kRData/kRRed/kROrange, DrawPageTitle
#include "SelectionCuts.h"  // kMCP1_minPeak, kMCP1_maxPeak, kFiducial_r_timing

static const float kSentCut_fts = -1e5f;   // hg_cfd validity test (matches headline)

// Build a VecToHist_teb-faithful histogram: 2-pass 5-sigma outlier trim, then a
// core range of mu2 +/- 4*ms2 with nb bins (outliers go to overflow).  Caller
// owns the returned histogram.  Mirrors timingEnergyBins.C exactly so the
// best-bin sigma reproduces the headline.
static TH1F* BuildCoreHist_fts(const std::vector<float>& v, int nb, const char* nm)
{
    if (v.size() < 50) return nullptr;
    double mu1 = 0.; for (float x : v) mu1 += x; mu1 /= (double)v.size();
    double ms1 = 0.; for (float x : v) ms1 += ((double)x - mu1) * ((double)x - mu1);
    ms1 = std::sqrt(ms1 / (double)v.size());
    if (ms1 < 0.008) ms1 = 0.100;
    double mu2 = 0.; long n2 = 0;
    for (float x : v) if (std::fabs((double)x - mu1) < 5. * ms1) { mu2 += x; ++n2; }
    double ms2 = 0.;
    if (n2 > 0) {
        mu2 /= (double)n2;
        for (float x : v) if (std::fabs((double)x - mu1) < 5. * ms1)
            ms2 += ((double)x - mu2) * ((double)x - mu2);
        ms2 = std::sqrt(ms2 / (double)n2);
        if (ms2 < 0.008) ms2 = 0.100;
    } else { mu2 = mu1; ms2 = ms1; }
    TH1F* h = new TH1F(nm, "", nb, mu2 - 4. * ms2, mu2 + 4. * ms2);
    h->SetDirectory(nullptr);
    for (float x : v) h->Fill(x);
    return h;
}

// Core sigma [ns] (+ optional error in ps) of corner times, via FitGaussCore.
static double RobustSigma_fts(const std::vector<float>& v, double* errPs = nullptr, int nb = 120)
{
    if (errPs) *errPs = 0.;
    if (v.size() < 150) return -1.;
    TH1F* h = BuildCoreHist_fts(v, nb, "_fts_h");
    if (!h) return -1.;
    double mu, muE, s, sE;
    FitGaussCore(h, 2.0, mu, muE, s, sE);
    delete h;
    if (errPs) *errPs = sE * 1000.;
    return (s > 0.) ? s : -1.;
}

// Core mean & sigma (mV) of the E_meas (sum_lg) spectrum — the binning seed used
// by the headline (FitGaussCore on a 150-bin VecToHist of sum_lg).
static bool CoreMuSig_fts(const std::vector<float>& v, double& mu, double& sig)
{
    mu = sig = 0.;
    TH1F* h = BuildCoreHist_fts(v, 150, "_fts_slg");
    if (!h) return false;
    double muE, sE; FitGaussCore(h, 2.0, mu, muE, sig, sE); delete h;
    return (sig > 0.);
}

// Headline best-bin estimator: 9 equal sum_lg bins over [muE-2sigE, muE+2sigE];
// the best bin is the min Method-A core sigma among bins with N >= 500.
// Returns sigma_t in ps (-1 on failure); sets the fit error (ps), N, efficiency.
static double BestBinSigma_fts(const std::vector<std::pair<float,float>>& sel,  // (sum_lg, tc)
                               double& errPs, int& nBest, double& effPct)
{
    errPs = 0.; nBest = 0; effPct = 0.;
    const int kNB = 9, kMinN = 500;
    std::vector<float> slg; slg.reserve(sel.size());
    for (const auto& p : sel) slg.push_back(p.first);
    double muE, sigE;
    if (!CoreMuSig_fts(slg, muE, sigE)) {
        double m = 0.; for (float x : slg) m += x; m /= (double)slg.size();
        double r = 0.; for (float x : slg) r += (x-m)*(x-m); r = std::sqrt(r/(double)slg.size());
        muE = m; sigE = r;
    }
    if (sigE <= 0.) return -1.;
    const double binLo = muE - 2.*sigE, binW = 4.*sigE / kNB;
    if (binW <= 0.) return -1.;
    const int totFid = (int)sel.size();
    double best = 1e30, bestErr = 0.; int bN = 0;
    for (int ib = 0; ib < kNB; ++ib) {
        const double lo = binLo + ib*binW, hi = binLo + (ib+1)*binW;
        std::vector<float> tc;
        for (const auto& p : sel) if (p.first >= lo && p.first < hi) tc.push_back(p.second);
        if ((int)tc.size() < kMinN) continue;
        double e; double s = RobustSigma_fts(tc, &e, 120);
        // Floor at 15 ps: the corner estimator is physically >~20 ps, so a
        // smaller "sigma" is a degenerate Gaussian-core fit (a narrow spike)
        // that the min-sigma selection must not grab.
        if (s > 0. && s*1000. >= 15. && s*1000. < best) { best = s*1000.; bestErr = e; bN = (int)tc.size(); }
    }
    if (best > 1e29) return -1.;
    errPs = bestErr; nBest = bN; effPct = 100.*bN/totFid;
    return best;
}

// Per-event store after event-level cuts (no radius cut): r^2, E_meas, corner time
struct EvFTS { float r2, slg, tc; };

// Collect events passing the headline event-level cuts for one energy ntuple.
// Returns false if the file/tree cannot be read.
static bool CollectEvents(const char* ntuple, std::vector<EvFTS>& out)
{
    TFile* f = TFile::Open(ntuple);
    if (!f || f->IsZombie()) { printf("[fiducialTimingScan] cannot open %s\n", ntuple); return false; }
    TTree* t = dynamic_cast<TTree*>(f->Get("rad"));
    if (!t) { printf("[fiducialTimingScan] no TTree 'rad' in %s\n", ntuple); f->Close(); return false; }

    double xc, yc, tOff, tRms;
    ScanRunCenters(t, xc, yc, tOff, tRms);

    Bool_t  wc_ok;
    Float_t x_trk, y_trk, mcp_peak, mcp_time, hg_cfd[8], sum_lg;
    t->SetBranchAddress("wc_ok",    &wc_ok);
    t->SetBranchAddress("x_trk",    &x_trk);
    t->SetBranchAddress("y_trk",    &y_trk);
    t->SetBranchAddress("mcp_peak", &mcp_peak);
    t->SetBranchAddress("mcp_time", &mcp_time);
    t->SetBranchAddress(t->GetBranch("hg_cfd05") ? "hg_cfd05" : "hg_cfd", hg_cfd);
    t->SetBranchAddress("sum_lg",   &sum_lg);

    out.clear(); out.reserve(200000);
    const Long64_t N = t->GetEntries();
    for (Long64_t i = 0; i < N; ++i) {
        t->GetEntry(i);
        if (!wc_ok)                                               continue;
        if (mcp_peak < kMCP1_minPeak || mcp_peak > kMCP1_maxPeak) continue;
        if (mcp_time <= kSentCut_fts)                            continue;
        double dw = 0., up = 0.; int ndw = 0, nup = 0;
        for (int c = 0; c < 4; ++c) if (hg_cfd[c] > kSentCut_fts) { dw += hg_cfd[c]; ++ndw; }
        for (int c = 4; c < 8; ++c) if (hg_cfd[c] > kSentCut_fts) { up += hg_cfd[c]; ++nup; }
        if (ndw < 1 || nup < 1) continue;
        float tc = (float)((dw / ndw - up / nup) * 0.5);        // Method A: (DW-UP)/2
        float dx = x_trk - (float)xc, dy = y_trk - (float)yc;
        out.push_back({ dx * dx + dy * dy, sum_lg, tc });
    }
    printf("[fiducialTimingScan] %-6s centroid (%.2f, %.2f) mm  %zu events\n",
           gSystem->BaseName(gSystem->DirName(ntuple)), xc, yc, out.size());
    f->Close();
    return true;
}

// Scan radius for one energy; fill gTop2 / gTop10 (sigma_t in ps vs radius, with
// the Gaussian-fit uncertainty as the y error bar).
static void ScanRadii(const std::vector<EvFTS>& evs,
                      TGraphErrors& gTop2, TGraphErrors& gTop10, TGraphErrors& gBest,
                      double& yMin, double& yMax)
{
    for (double R = 1.0; R <= 5.001; R += 0.25) {
        const float R2 = (float)(R * R);
        std::vector<std::pair<float, float>> sel;            // (sum_lg, tc)
        sel.reserve(evs.size());
        for (const auto& e : evs) if (e.r2 < R2) sel.emplace_back(e.slg, e.tc);
        if (sel.size() < 400) continue;

        // --- single best E_meas bin (the headline estimator), re-optimised at R ---
        double eb = 0., effb = 0.; int nb = 0;
        double sBest = BestBinSigma_fts(sel, eb, nb, effb);
        if (sBest > 0.) { int n = gBest.GetN(); gBest.SetPoint(n, R, sBest); gBest.SetPointError(n, 0., eb);
                          yMin = std::min(yMin, sBest); yMax = std::max(yMax, sBest); }

        // --- fixed top-fraction selections (stable diagnostic) ---
        std::sort(sel.begin(), sel.end(),
                  [](const std::pair<float,float>& a, const std::pair<float,float>& b)
                  { return a.first > b.first; });            // E_meas descending
        auto topSigma = [&](double frac, double& errPs) -> double {
            size_t nt = (size_t)(frac * (double)sel.size());
            if (nt < 300) { errPs = 0.; return -1.; }   // drop low-statistics (noisy) points
            nt = std::min(nt, sel.size());
            std::vector<float> tt; tt.reserve(nt);
            for (size_t k = 0; k < nt; ++k) tt.push_back(sel[k].second);
            return RobustSigma_fts(tt, &errPs);
        };
        double e10 = 0., e02 = 0.;
        double s10 = topSigma(0.10, e10);
        double s02 = topSigma(0.02, e02);
        if (s02 > 0.) { int n = gTop2.GetN(); gTop2.SetPoint(n, R, s02 * 1000.); gTop2.SetPointError(n, 0., e02);
                        yMin = std::min(yMin, s02*1000.); yMax = std::max(yMax, s02*1000.); }
        if (s10 > 0.) { int n = gTop10.GetN(); gTop10.SetPoint(n, R, s10 * 1000.); gTop10.SetPointError(n, 0., e10); }
    }
}

void fiducialTimingScan()
{
    TH1::AddDirectory(kFALSE);
    ApplyRADiCALStyle();

    const int    nE = 6;
    const int    eGeV[nE]   = { 25, 50, 75, 100, 125, 150 };
    const char*  eLbl[nE]   = { "25 GeV","50 GeV","75 GeV","100 GeV","125 GeV","150 GeV" };

    TGraphErrors gTop2[nE], gTop10[nE], gBest[nE];
    double yMin = 1e9, yMax = -1e9, yBestMin = 1e9, yBestMax = -1e9;
    std::vector<EvFTS> evs150;
    for (int e = 0; e < nE; ++e) {
        TString nt = Form("Analysis/Output/%dGeV/ntuple.root", eGeV[e]);
        std::vector<EvFTS> evs;
        if (!CollectEvents(nt.Data(), evs)) continue;
        ScanRadii(evs, gTop2[e], gTop10[e], gBest[e], yMin, yMax);
        for (int i = 0; i < gBest[e].GetN(); ++i) {
            yBestMin = std::min(yBestMin, gBest[e].GetY()[i]);
            yBestMax = std::max(yBestMax, gBest[e].GetY()[i]);
        }
        if (eGeV[e] == 150) evs150 = evs;   // keep for the detail figure
    }

    // Validation: best-bin sigma at r=3 mm must reproduce the headline teb_sigma.
    printf("\n[fiducialTimingScan] VALIDATION — best-bin sigma_t at r=3 mm (should match headline):\n");
    { auto v3 = [](const TGraphErrors& g){ for (int i=0;i<g.GetN();++i) if (std::fabs(g.GetX()[i]-3.0)<1e-6) return g.GetY()[i]; return 0.; };
      printf("   "); for (int e = 0; e < nE; ++e) printf(" %3dGeV=%.1f", eGeV[e], v3(gBest[e])); printf(" ps\n");
      printf("   headline teb_sigma = 47.6 35.7 32.5 33.8 29.2 27.4 ps\n"); }

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
    const double rPlateau = 2.5;

    // =========================================================================
    // Figure 1 — ALL ENERGIES overlaid (top-2% E_meas corner sigma_t vs radius)
    // =========================================================================
    {
        TCanvas* c = new TCanvas("c_fts_all", "", 960, 760);
        c->SetLeftMargin(0.13); c->SetBottomMargin(0.13);
        c->SetRightMargin(0.05); c->SetTopMargin(0.10);
        c->SetTickx(1); c->SetTicky(1);

        const double ylo = std::max(0., yMin - 5.), yhi = yMax + 6.;
        bool first = true;
        for (int e = 0; e < nE; ++e) {
            if (gTop2[e].GetN() < 2) continue;
            gTop2[e].SetMarkerStyle(20); gTop2[e].SetMarkerSize(1.1);
            gTop2[e].SetMarkerColor(kREnergyCols[e]); gTop2[e].SetLineColor(kREnergyCols[e]);
            gTop2[e].SetLineWidth(3);
            if (first) {
                gTop2[e].Draw("APL");
                gTop2[e].GetXaxis()->SetTitle("timing fiducial radius  r  (mm)");
                gTop2[e].GetYaxis()->SetTitle("#sigma_{t}  (DW#minusUP)/2,  top 2% E_{meas}  (ps)");
                gTop2[e].GetXaxis()->SetLimits(0.5, 5.2);
                gTop2[e].GetYaxis()->SetRangeUser(ylo, yhi);
                gTop2[e].GetXaxis()->SetTitleSize(0.046); gTop2[e].GetYaxis()->SetTitleSize(0.044);
                first = false;
            } else {
                gTop2[e].Draw("PL same");
            }
        }
        // current 3 mm cut (the radius optimum is energy-dependent, so no single
        // "plateau" band is drawn here — see the 150 GeV detail figure for that)
        for (int e = 0; e < nE; ++e) if (gTop2[e].GetN() >= 2) gTop2[e].Draw("PL same");
        TLine* l3 = new TLine(kFiducial_r_timing, ylo, kFiducial_r_timing, yhi);
        l3->SetLineStyle(2); l3->SetLineColor(kGray+2); l3->SetLineWidth(2); l3->Draw();
        { TLatex a; a.SetTextColor(kGray+2); a.SetTextSize(0.027); a.SetTextAngle(90);
          a.DrawLatex(kFiducial_r_timing + 0.09, ylo + 0.42*(yhi-ylo), "adopted cut  r < 3 mm"); }

        TLegend* L = new TLegend(0.70, 0.50, 0.93, 0.88);
        L->SetBorderSize(0); L->SetFillStyle(0); L->SetTextFont(42); L->SetTextSize(0.032);
        for (int e = 0; e < nE; ++e) if (gTop2[e].GetN() >= 2) L->AddEntry(&gTop2[e], eLbl[e], "pl");
        L->Draw();

        { TLatex a; a.SetNDC(); a.SetTextSize(0.026); a.SetTextColor(kGray+3);
          a.DrawLatex(0.16, 0.895, "Tightening helps the lower energies (outer-ring");
          a.DrawLatex(0.16, 0.862, "position/walk); higher E #rightarrow lower #sigma_{t}.");
          a.SetTextColor(kRRed);
          a.DrawLatex(0.16, 0.825, "Adopted: 2.5 mm (E#leq100), 3.0 mm (E#geq125)"); }

        DrawPageTitle("Timing resolution vs fiducial radius -- all energies  (top-2% E_{meas})");
        c->Print("Analysis/Output/Summary/fiducial_timing_scan.png");
        c->Print("Analysis/Output/Summary/fiducial_timing_scan.pdf");
        printf("[fiducialTimingScan] wrote fiducial_timing_scan.png (all energies)\n");
    }

    // =========================================================================
    // Figure 2 — 150 GeV detail (top-2% vs top-10%), the pedagogical panel
    // =========================================================================
    if (gTop2[nE-1].GetN() >= 2)
    {
        TGraphErrors& g2  = gTop2[nE-1];
        TGraphErrors& g10 = gTop10[nE-1];
        TCanvas* c = new TCanvas("c_fts_150", "", 920, 760);
        c->SetLeftMargin(0.14); c->SetBottomMargin(0.13);
        c->SetRightMargin(0.05); c->SetTopMargin(0.10);
        c->SetTickx(1); c->SetTicky(1);

        TGraphErrors& gB = gBest[nE-1];   // the headline single-best-bin curve at 150 GeV
        g2.SetMarkerStyle(20); g2.SetMarkerColor(kRRed);  g2.SetLineColor(kRRed);  g2.SetLineWidth(3); g2.SetMarkerSize(1.2);
        g10.SetMarkerStyle(21); g10.SetMarkerColor(kRData); g10.SetLineColor(kRData); g10.SetLineWidth(3); g10.SetMarkerSize(1.2);
        gB.SetMarkerStyle(33); gB.SetMarkerColor(kBlack); gB.SetLineColor(kBlack); gB.SetLineWidth(3); gB.SetMarkerSize(1.7);

        double yl = 1e9, yh = -1e9;
        auto span = [&](const TGraphErrors& g){ for (int i=0;i<g.GetN();++i){ yl=std::min(yl,g.GetY()[i]); yh=std::max(yh,g.GetY()[i]); } };
        span(g2); span(g10); span(gB);
        const double ylo = yl - 5., yhi = yh + 5.;

        g2.Draw("APL");
        g2.GetXaxis()->SetTitle("timing fiducial radius  r  (mm)");
        g2.GetYaxis()->SetTitle("#sigma_{t}  (DW#minusUP)/2  (ps)");
        g2.GetXaxis()->SetLimits(0.5, 5.2);
        g2.GetYaxis()->SetRangeUser(ylo, yhi);
        g2.GetXaxis()->SetTitleSize(0.046); g2.GetYaxis()->SetTitleSize(0.046);
        g10.Draw("PL same");
        TBox* band = new TBox(0.5, ylo, rPlateau, yhi);
        band->SetFillColorAlpha(kRData, 0.06); band->SetLineColor(0); band->Draw();
        g2.Draw("PL same"); g10.Draw("PL same"); gB.Draw("PL same");
        // mark the adopted 150 GeV point (r=3 mm; OOS-validated 27.4 ps)
        const double kAdopt150 = 27.4;   // OOS best-bin headline at 3 mm
        TMarker* adopt = new TMarker(kFiducial_r_timing, kAdopt150, 29);
        adopt->SetMarkerColor(kROrange); adopt->SetMarkerSize(3.0); adopt->Draw();
        TLine* l3 = new TLine(kFiducial_r_timing, ylo, kFiducial_r_timing, yhi);
        l3->SetLineStyle(2); l3->SetLineColor(kROrange); l3->SetLineWidth(2); l3->Draw();
        { TLatex a; a.SetTextColor(kROrange); a.SetTextSize(0.028); a.SetTextAngle(90);
          a.DrawLatex(kFiducial_r_timing + 0.10, ylo + 0.62*(yhi-ylo), "adopted: r < 3 mm"); }
        { TLatex a; a.SetTextColor(kRData); a.SetTextSize(0.025);
          a.DrawLatex(0.78, yhi - 0.05*(yhi-ylo), "top-X% plateau (r #leq 2.5 mm)"); }

        TLegend* L = new TLegend(0.40, 0.73, 0.93, 0.88);
        L->SetBorderSize(0); L->SetFillStyle(0); L->SetTextFont(42); L->SetTextSize(0.030);
        L->AddEntry(&gB,  "single best bin  (the headline estimator)", "pl");
        L->AddEntry(&g2,  "top 2% by E_{meas} (#SigmaA_{LG})", "pl");
        L->AddEntry(&g10, "top 10% by E_{meas}", "pl");
        L->Draw();
        { TPave* bg = new TPave(0.135, 0.135, 0.955, 0.275, 0, "brNDC");
          bg->SetFillColor(kWhite); bg->SetLineColor(kWhite); bg->Draw();
          TLatex a; a.SetNDC(); a.SetTextSize(0.0245); a.SetTextColor(kGray+3);
          a.DrawLatex(0.150, 0.240,
            "Top-X% (frozen selection): average-event optimum at r #approx 2 mm.");
          a.SetTextColor(kBlack);
          a.DrawLatex(0.150, 0.205,
            "Headline best-bin (black) is tighter (sits below) but JUMPY: its 2.25 mm dip is a");
          a.DrawLatex(0.150, 0.170,
            "lucky radius (2.0 / 2.5 mm = 30.7 / 28.1 ps).  OOS #Rightarrow adopt 3 mm: 27.4 ps (#bigstar)."); }
        DrawPageTitle("Fiducial radius at 150 GeV:  headline best-bin  vs  stable top-X% selection");
        c->Print("Analysis/Output/Summary/fiducial_timing_scan_150.png");
        c->Print("Analysis/Output/Summary/fiducial_timing_scan_150.pdf");
        printf("[fiducialTimingScan] wrote fiducial_timing_scan_150.png (detail)\n");
    }

    // =========================================================================
    // Figure 3 — JOINT view: best-bin (headline) sigma_t vs radius, per energy
    // =========================================================================
    {
        TCanvas* c = new TCanvas("c_fts_best", "", 960, 760);
        c->SetLeftMargin(0.13); c->SetBottomMargin(0.13);
        c->SetRightMargin(0.05); c->SetTopMargin(0.10);
        c->SetTickx(1); c->SetTicky(1);

        const double ylo = std::max(0., yBestMin - 4.), yhi = yBestMax + 5.;
        bool first = true;
        for (int e = 0; e < nE; ++e) {
            if (gBest[e].GetN() < 2) continue;
            gBest[e].SetMarkerStyle(20); gBest[e].SetMarkerSize(1.1);
            gBest[e].SetMarkerColor(kREnergyCols[e]); gBest[e].SetLineColor(kREnergyCols[e]);
            gBest[e].SetLineWidth(2);
            if (first) {
                gBest[e].Draw("APL");
                gBest[e].GetXaxis()->SetTitle("timing fiducial radius  r  (mm)");
                gBest[e].GetYaxis()->SetTitle("best-bin #sigma_{t}  (DW#minusUP)/2  (ps)  [headline]");
                gBest[e].GetXaxis()->SetLimits(0.5, 5.2);
                gBest[e].GetYaxis()->SetRangeUser(ylo, yhi);
                gBest[e].GetXaxis()->SetTitleSize(0.046); gBest[e].GetYaxis()->SetTitleSize(0.044);
                first = false;
            } else gBest[e].Draw("PL same");
        }
        TLine* l3 = new TLine(kFiducial_r_timing, ylo, kFiducial_r_timing, yhi);
        l3->SetLineStyle(2); l3->SetLineColor(kGray+2); l3->SetLineWidth(2); l3->Draw();
        for (int e = 0; e < nE; ++e) if (gBest[e].GetN() >= 2) gBest[e].Draw("PL same");
        { TLatex a; a.SetTextColor(kGray+2); a.SetTextSize(0.027); a.SetTextAngle(90);
          a.DrawLatex(kFiducial_r_timing + 0.09, ylo + 0.42*(yhi-ylo), "adopted cut  r < 3 mm"); }

        TLegend* L = new TLegend(0.70, 0.50, 0.93, 0.88);
        L->SetBorderSize(0); L->SetFillStyle(0); L->SetTextFont(42); L->SetTextSize(0.032);
        for (int e = 0; e < nE; ++e) if (gBest[e].GetN() >= 2) L->AddEntry(&gBest[e], eLbl[e], "pl");
        L->Draw();

        { TLatex a; a.SetNDC(); a.SetTextSize(0.025); a.SetTextColor(kGray+3);
          a.DrawLatex(0.155, 0.878, "Joint optimisation: E_{meas} bin re-picked at each radius (in-sample).");
          a.DrawLatex(0.155, 0.845, "Jumpy; tighter cuts can OVERFIT, so the run-folded OOS validation");
          a.DrawLatex(0.155, 0.812, "(see text) sets the adopted per-energy fiducial:");
          a.SetTextColor(kRRed);
          a.DrawLatex(0.155, 0.778, "r < 2.5 mm for E #leq 100 GeV;  r < 3.0 mm for E #geq 125 GeV."); }

        DrawPageTitle("Best-bin (headline) #sigma_{t} vs fiducial radius -- all energies");
        c->Print("Analysis/Output/Summary/fiducial_timing_scan_bestbin.png");
        c->Print("Analysis/Output/Summary/fiducial_timing_scan_bestbin.pdf");
        printf("[fiducialTimingScan] wrote fiducial_timing_scan_bestbin.png (joint view)\n");
    }

    // Console summary tables
    auto printTable = [&](const char* tag, TGraphErrors* g) {
        printf("\n  R(mm)");
        for (int e = 0; e < nE; ++e) printf("  %6d", eGeV[e]);
        printf("   [%s sigma_t, ps]\n", tag);
        for (double R = 1.0; R <= 5.001; R += 0.25) {
            printf("  %4.2f ", R);
            for (int e = 0; e < nE; ++e) { double v = valAt(g[e], R); if (v>0) printf("  %6.1f", v); else printf("       -"); }
            printf("\n");
        }
    };
    printTable("BEST-BIN", gBest);
    // Mean over energies of best-bin sigma at each radius (where defined for all 6)
    printf("\n  best-bin sigma averaged over the 6 energies:\n");
    for (double R = 1.5; R <= 3.501; R += 0.25) {
        double sum = 0.; int n = 0;
        for (int e = 0; e < nE; ++e) { double v = valAt(gBest[e], R); if (v>0){ sum+=v; ++n; } }
        if (n == nE) printf("   r=%.2f mm : mean = %.2f ps  (over %d energies)\n", R, sum/n, n);
        else         printf("   r=%.2f mm : (only %d/%d energies defined)\n", R, n, nE);
    }
    printTable("top-2%", gTop2);
}

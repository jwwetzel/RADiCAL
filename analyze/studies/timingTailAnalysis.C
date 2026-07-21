// ============================================================================
// timingTailAnalysis.C — what are the tails of the (DW-UP)/2 timing distribution?
// ----------------------------------------------------------------------------
// The headline quotes the Gaussian-core sigma (~27 ps); the full RMS is ~31 ps,
// i.e. the distribution has non-Gaussian tails.  This macro characterises them
// for the high-E_meas regime the headline operates in (150 GeV, r<3 mm, top-2%
// by ΣA_LG), testing three hypotheses for the tail events:
//   (1) low amplitude  -> residual time-walk,
//   (2) few valid channels -> noisier per-event average,
//   (3) a single rogue channel -> CFD mis-reconstruction (large intra-group dev).
//
// Output: output/Summary/timing_tails.{png,pdf} + console stats.
// ============================================================================

#include <vector>
#include <algorithm>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TString.h"
#include "TSystem.h"

#include "PlotUtils.h"      // ScanRunCenters, FitGaussCore
#include "RADiCALStyle.h"   // ApplyRADiCALStyle, kRData/kRRed, DrawPageTitle
#include "SelectionCuts.h"  // kMCP1_minPeak/maxPeak, kFiducial_r_timing

static const float kSentT = -1e5f;

void timingTailAnalysis()
{
    TH1::AddDirectory(kFALSE);
    ApplyRADiCALStyle();

    TFile* f = TFile::Open("output/150GeV/ntuple.root");
    if (!f || f->IsZombie()) { printf("[tails] cannot open ntuple\n"); return; }
    TTree* t = dynamic_cast<TTree*>(f->Get("rad"));
    if (!t) { printf("[tails] no tree\n"); return; }
    double xc, yc, to, tr; ScanRunCenters(t, xc, yc, to, tr);

    Bool_t  wc_ok;
    Float_t x_trk, y_trk, mcp_peak, mcp_time, mcp2_time = -1e6f, hg_cfd[8], sum_lg;
    t->SetBranchAddress("wc_ok",    &wc_ok);
    t->SetBranchAddress("x_trk",    &x_trk);
    t->SetBranchAddress("y_trk",    &y_trk);
    t->SetBranchAddress(t->GetBranch("mcp1_peak")?"mcp1_peak":"mcp_peak", &mcp_peak);
    t->SetBranchAddress(t->GetBranch("mcp1_time")?"mcp1_time":"mcp_time", &mcp_time);
    t->SetBranchAddress(t->GetBranch("hg_cfd05") ? "hg_cfd05" : "hg_cfd", hg_cfd);
    t->SetBranchAddress("sum_lg",   &sum_lg);
    const bool hasM2 = (t->GetBranch("mcp2_time") != nullptr);
    if (hasM2) t->SetBranchAddress("mcp2_time", &mcp2_time);

    const float rFid2 = (float)(kFiducial_r_timing * kFiducial_r_timing);

    // ---- pass 1: ΣLG of fiducial events -> top-2% threshold ----
    std::vector<float> slgs;
    const Long64_t N = t->GetEntries();
    for (Long64_t i = 0; i < N; ++i) {
        t->GetEntry(i);
        if (!wc_ok || mcp_peak < kMCP1_minPeak || mcp_peak > kMCP1_maxPeak || mcp_time <= kSentT) continue;
        float dx = x_trk-(float)xc, dy = y_trk-(float)yc; if (dx*dx+dy*dy >= rFid2) continue;
        slgs.push_back(sum_lg);
    }
    std::sort(slgs.begin(), slgs.end());
    const float slgThr = slgs[(size_t)(0.98*slgs.size())];   // top 2% by E_meas (headline regime)

    // ---- pass 2: per-event timing + diagnostics for the high-E_meas sample ----
    struct Ev { float tc, slg, r, maxdev; int nch; };
    std::vector<Ev> evs; evs.reserve(40000);
    for (Long64_t i = 0; i < N; ++i) {
        t->GetEntry(i);
        if (!wc_ok || mcp_peak < kMCP1_minPeak || mcp_peak > kMCP1_maxPeak || mcp_time <= kSentT) continue;
        float dx = x_trk-(float)xc, dy = y_trk-(float)yc; float r2 = dx*dx+dy*dy;
        if (r2 >= rFid2) continue;
        if (sum_lg < slgThr) continue;

        double dw = 0., up = 0.; int ndw = 0, nup = 0;
        for (int c = 0; c < 4; ++c) if (hg_cfd[c] > kSentT) { dw += hg_cfd[c]; ++ndw; }
        for (int c = 4; c < 8; ++c) if (hg_cfd[c] > kSentT) { up += hg_cfd[c]; ++nup; }
        if (ndw < 1 || nup < 1) continue;
        const double dwa = dw/ndw, upa = up/nup;
        const float tc = (float)((dwa - upa) * 0.5);

        // worst single-channel deviation within the MCP1-referenced groups
        // (DW ch0-3, UP ch4-6; ch7 is MCP2-referenced so excluded from this test).
        double mdev = 0.;
        for (int c = 0; c < 4; ++c) if (hg_cfd[c] > kSentT) mdev = std::max(mdev, std::fabs(hg_cfd[c] - dwa));
        { double m = 0.; int n = 0; for (int c = 4; c < 7; ++c) if (hg_cfd[c] > kSentT) { m += hg_cfd[c]; ++n; }
          if (n > 0) { m /= n; for (int c = 4; c < 7; ++c) if (hg_cfd[c] > kSentT) mdev = std::max(mdev, std::fabs(hg_cfd[c]-m)); } }

        evs.push_back({ tc, sum_lg, std::sqrt(r2), (float)(mdev*1000.), ndw+nup });
    }
    f->Close();
    printf("[tails] %zu high-E_meas events (top 2%%, r<3mm, 150 GeV); ΣLG>%.0f mV\n", evs.size(), slgThr);

    // ---- core fit + tail classification ----
    std::vector<float> tcs; tcs.reserve(evs.size());
    for (auto& e : evs) tcs.push_back(e.tc);
    std::sort(tcs.begin(), tcs.end());
    const float med = tcs[tcs.size()/2];
    const float kWin = 0.25f;                       // +/-250 ps: core+shoulder window
    // bounded histogram: catastrophic (>250 ps) outliers go to overflow (excluded
    // from the fit and the 'shoulder RMS'), and are counted separately as mis-reco.
    TH1F* hAll = new TH1F("hAll", "", 200, med-kWin, med+kWin); hAll->SetDirectory(nullptr);
    long nCat = 0;
    for (float x : tcs) { hAll->Fill(x); if (std::fabs(x-med) > kWin) ++nCat; }
    double mu, muE, s, sE; FitGaussCore(hAll, 2.0, mu, muE, s, sE);
    const double sigPs = s*1000.;
    const double rmsPs = hAll->GetRMS()*1000.;       // shoulder RMS (within +/-250 ps)
    const double catFrac = 100.*nCat/tcs.size();      // catastrophic (>250 ps) fraction
    const double tailCut = 2.5*s;   // |tc-mu| > 2.5 sigma_core = (mild) tail
    printf("[tails] catastrophic (>250 ps) outliers: %ld  (%.2f%% -- mis-reconstruction)\n", nCat, catFrac);

    // ---- compare tail vs core ----
    TH1F* hSlgC = new TH1F("hSlgC","",60,slgThr,7000); TH1F* hSlgT = new TH1F("hSlgT","",60,slgThr,7000);
    TH1F* hDevC = new TH1F("hDevC","",60,0,400);       TH1F* hDevT = new TH1F("hDevT","",60,0,400);
    for (auto* h : {hSlgC,hSlgT,hDevC,hDevT}) h->SetDirectory(nullptr);
    long nC = 0, nT = 0; double full8C = 0, full8T = 0;
    for (auto& e : evs) {
        const bool tail = std::fabs(e.tc - mu) > tailCut;
        if (tail) { ++nT; hSlgT->Fill(e.slg); hDevT->Fill(e.maxdev); if (e.nch==8) ++full8T; }
        else      { ++nC; hSlgC->Fill(e.slg); hDevC->Fill(e.maxdev); if (e.nch==8) ++full8C; }
    }
    const double tailFrac = 100.*nT/(nC+nT);
    const double gausExp  = 100.*std::erfc(2.5/std::sqrt(2.0));   // Gaussian expectation beyond 2.5 sigma
    printf("[tails] sigma_core=%.1f ps  RMS=%.1f ps   tail (|.|>2.5sig)=%.2f%%  (pure-Gaussian would be %.2f%%)\n",
           sigPs, rmsPs, tailFrac, gausExp);
    printf("[tails] median ΣLG:  core=%.0f  tail=%.0f mV   |   full-8ch fraction: core=%.1f%% tail=%.1f%%\n",
           hSlgC->GetMean(), hSlgT->GetMean(), 100.*full8C/nC, 100.*full8T/nT);
    printf("[tails] mean worst-channel deviation:  core=%.1f ps   tail=%.1f ps\n",
           hDevC->GetMean(), hDevT->GetMean());

    // ====================== figure ======================
    TCanvas* c = new TCanvas("c_tails", "", 1280, 1080);
    c->Divide(2, 2, 0.005, 0.02);

    // (a) the distribution, log y, core fit overlaid
    c->cd(1); gPad->SetLogy(); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.12); gPad->SetTopMargin(0.10);
    hAll->SetLineColor(kRData); hAll->SetLineWidth(2); hAll->SetFillColorAlpha(kRData,0.20);
    hAll->GetXaxis()->SetTitle("(DW#minusUP)/2  (ns)"); hAll->GetYaxis()->SetTitle("events");
    hAll->GetXaxis()->SetTitleSize(0.05); hAll->GetYaxis()->SetTitleSize(0.05); hAll->SetMinimum(0.5);
    hAll->Draw("HIST");
    TF1* fg = new TF1("fg","[0]*exp(-0.5*((x-[1])/[2])^2)", mu-6*s, mu+6*s);
    fg->SetParameters(hAll->GetMaximum(), mu, s); fg->SetLineColor(kRRed); fg->SetLineWidth(3); fg->Draw("same");
    for (int sgn=-1; sgn<=1; sgn+=2){ TLine* l=new TLine(mu+sgn*tailCut, 0.5, mu+sgn*tailCut, hAll->GetMaximum());
        l->SetLineColor(kBlack); l->SetLineStyle(2); l->SetLineWidth(2); l->Draw(); }
    { TLatex a; a.SetNDC(); a.SetTextSize(0.045);
      a.SetTextColor(kRRed); a.DrawLatex(0.16,0.84,Form("#sigma_{core} = %.1f ps",sigPs));
      a.SetTextColor(kGray+2); a.DrawLatex(0.16,0.78,Form("RMS = %.1f ps",rmsPs));
      a.SetTextColor(kBlack); a.DrawLatex(0.16,0.72,Form("tail (>2.5#sigma) = %.1f%%",tailFrac));
      a.SetTextColor(kGray+3); a.SetTextSize(0.036); a.DrawLatex(0.16,0.66,Form("(Gaussian: %.1f%%)",gausExp)); }
    DrawPadTitle("(DW#minusUP)/2 with core fit (dashed = tail boundary)");

    // (b) amplitude: tail vs core (normalised)
    c->cd(2); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.12); gPad->SetTopMargin(0.10);
    hSlgC->Scale(1./hSlgC->Integral()); hSlgT->Scale(1./hSlgT->Integral());
    hSlgC->SetLineColor(kRData); hSlgC->SetLineWidth(3); hSlgT->SetLineColor(kRRed); hSlgT->SetLineWidth(3);
    hSlgC->GetXaxis()->SetTitle("#SigmaA_{LG}  (mV)"); hSlgC->GetYaxis()->SetTitle("fraction");
    hSlgC->GetXaxis()->SetTitleSize(0.05); hSlgC->GetYaxis()->SetTitleSize(0.05);
    hSlgC->SetMaximum(1.25*std::max(hSlgC->GetMaximum(),hSlgT->GetMaximum())); hSlgC->Draw("HIST"); hSlgT->Draw("HIST same");
    { TLegend* L=new TLegend(0.55,0.74,0.93,0.88); L->SetBorderSize(0); L->SetFillStyle(0); L->SetTextSize(0.042);
      L->AddEntry(hSlgC,"core events","l"); L->AddEntry(hSlgT,"tail events","l"); L->Draw(); }
    DrawPadTitle("Amplitude:  are tails low-E_{meas}?");

    // (c) worst-channel deviation: tail vs core
    c->cd(3); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.12); gPad->SetTopMargin(0.10);
    hDevC->Scale(1./hDevC->Integral()); hDevT->Scale(1./hDevT->Integral());
    hDevC->SetLineColor(kRData); hDevC->SetLineWidth(3); hDevT->SetLineColor(kRRed); hDevT->SetLineWidth(3);
    hDevC->GetXaxis()->SetTitle("worst single-channel deviation  (ps)"); hDevC->GetYaxis()->SetTitle("fraction");
    hDevC->GetXaxis()->SetTitleSize(0.05); hDevC->GetYaxis()->SetTitleSize(0.05);
    hDevC->SetMaximum(1.25*std::max(hDevC->GetMaximum(),hDevT->GetMaximum())); hDevC->Draw("HIST"); hDevT->Draw("HIST same");
    { TLegend* L=new TLegend(0.55,0.74,0.93,0.88); L->SetBorderSize(0); L->SetFillStyle(0); L->SetTextSize(0.042);
      L->AddEntry(hDevC,"core events","l"); L->AddEntry(hDevT,"tail events","l"); L->Draw(); }
    DrawPadTitle("Mis-reco:  do tails have a rogue channel?");

    // (d) text summary / verdict
    c->cd(4); gPad->SetLeftMargin(0.05);
    { TLatex a; a.SetTextSize(0.052); a.SetTextFont(42);
      a.DrawLatex(0.06,0.86,"Tail summary (150 GeV, top-2% E_{meas}):");
      a.SetTextSize(0.046); a.SetTextColor(kGray+3);
      a.DrawLatex(0.08,0.74,Form("#bullet tail fraction = %.1f%%  (Gaussian: %.1f%%)",tailFrac,gausExp));
      a.DrawLatex(0.08,0.65,Form("#bullet mean #SigmaA_{LG}: tail %.0f vs core %.0f mV",hSlgT->GetMean(),hSlgC->GetMean()));
      a.DrawLatex(0.08,0.56,Form("#bullet worst-ch dev: tail %.0f vs core %.0f ps",hDevT->GetMean(),hDevC->GetMean()));
      a.DrawLatex(0.08,0.47,Form("#bullet full-8-channel: tail %.0f%% vs core %.0f%%",100.*full8T/nT,100.*full8C/nC)); }
    DrawPadTitle("Verdict");

    c->cd(0);
    DrawPageTitle("Non-Gaussian tails of the (DW#minusUP)/2 timing  (150 GeV)");
    c->Print("output/Summary/timing_tails.png");
    c->Print("output/Summary/timing_tails.pdf");
    printf("[tails] wrote timing_tails.png\n");
}

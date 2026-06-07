// ===========================================================================
// layer6Summary.C  —  compact, Ledovskoy-clean "Layer 6: Systematic
//                     Uncertainties" heroes.
//
//   layer6_budget.png   H1  systematic budget breakdown at 150 GeV (horizontal
//                           bars): which cut variation dominates, + quadrature total
//   layer6_band.png     H2  A^2-weighted combo sigma_t vs energy with the
//                           stat (+) systematic uncertainty band
//
// NB: the systematic study (systematicUncertainties.C) is run on the A^2-weighted
// COMBINATION (nominal ~83 ps), which is the estimator most sensitive to the
// cuts.  The dominant term is the MCP-window variation -- which the MCP-free
// (DW-UP)/2 headline does NOT carry, so the headline is even more robust.
//
// Run:  ROOT_INCLUDE_PATH=Analysis root -l -b -q 'Analysis/layer6Summary.C+'
// ===========================================================================
#include <cmath>

#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TAxis.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TSystem.h"

#include "RADiCALStyle.h"

namespace {

const char* kSumDir = "output/Summary/";

double at150(TFile& f, const char* g)
{
    TGraph* G = dynamic_cast<TGraph*>(f.Get(g));
    if (!G) return 0.;
    double x, y; for (int i = 0; i < G->GetN(); ++i) { G->GetPoint(i, x, y); if (fabs(x - 150.) < 1.) return y; }
    return 0.;
}

// ── H1 — systematic budget breakdown at 150 GeV ────────────────────────────
void HeroBudget()
{
    TFile f(Form("%ssystematic_uncertainties.root", kSumDir));
    if (f.IsZombie()) { printf("[layer6Summary] no systematic_uncertainties.root — skip H1\n"); return; }

    // Each cut group = quadrature of its two-sided variations, so the group bars
    // reconcile EXACTLY with the quadrature total (gSystTotal).
    auto quad = [](double a, double b){ return std::sqrt(a*a + b*b); };
    const double mcp  = quad(at150(f, "gDelta_MCP_lo_p50mV"), at150(f, "gDelta_MCP_hi_m50mV"));
    const double fid  = quad(at150(f, "gDelta_Fid_p0d5mm"),   at150(f, "gDelta_Fid_m0d5mm"));
    const double cont = quad(at150(f, "gDelta_Contain_p0d05"),at150(f, "gDelta_Contain_m0d05"));
    const double hg   = fabs(at150(f, "gDelta_HG_thresh_p5mV"));
    const double tot  = at150(f, "gSystTotal");

    // 5 bins, ascending so TOTAL lands at the top of the horizontal-bar chart.
    const char*  lab[5] = {"HG threshold", "fiducial radius", "containment cut",
                           "MCP window", "TOTAL (quadrature)"};
    const double val[5] = {hg, fid, cont, mcp, tot};

    TH1F* hC = new TH1F("hSystC", "", 5, 0.5, 5.5);   // components
    TH1F* hT = new TH1F("hSystT", "", 5, 0.5, 5.5);   // total (different colour)
    hC->SetDirectory(nullptr); hT->SetDirectory(nullptr);
    for (int i = 0; i < 5; ++i) {
        hC->GetXaxis()->SetBinLabel(i + 1, lab[i]);
        if (i < 4) hC->SetBinContent(i + 1, val[i]);
        else       hT->SetBinContent(i + 1, val[i]);
    }

    TCanvas* c = NewSquareCanvas("c_l6_budget", 600, 224, 92, 54, 40);  // wide left for labels
    c->cd();
    hC->SetFillColor(kRBlue);  hC->SetLineColor(kRBlue);  hC->SetBarWidth(0.72); hC->SetBarOffset(0.14);
    hT->SetFillColor(kROrange);hT->SetLineColor(kROrange);hT->SetBarWidth(0.72); hT->SetBarOffset(0.14);
    hC->GetXaxis()->SetLabelSize(0.041);
    hC->GetYaxis()->SetTitle("systematic shift in #sigma_{t}  (ps)");
    hC->GetYaxis()->SetTitleSize(0.045);
    hC->GetYaxis()->SetRangeUser(0., tot * 1.18);
    hC->Draw("HBAR");
    hT->Draw("HBAR SAME");

    { TLatex t; t.SetNDC(); t.SetTextSize(0.028); t.SetTextColor(kGray + 3);
      t.DrawLatex(0.36, 0.34, "Evaluated on the A^{2}-weighted combo");
      t.DrawLatex(0.36, 0.295, "(the most cut-sensitive estimator);");
      t.DrawLatex(0.36, 0.250, "every term is only a few ps.");
      t.SetTextColor(kROrange);
      t.DrawLatex(0.36, 0.195, "The MCP-window term does NOT apply");
      t.DrawLatex(0.36, 0.155, "to the MCP-free (DW#minusUP)/2 headline."); }

    DrawPageTitle("Systematic budget at 150 GeV  (A^{2}-weighted combo)");
    c->Print(Form("%slayer6_budget.png", kSumDir));
    printf("[layer6Summary] wrote layer6_budget.png\n");
}

// ── H2 — A^2-combo sigma_t vs energy with stat (+) systematic band ─────────
void HeroBand()
{
    TFile f(Form("%ssystematic_uncertainties.root", kSumDir));
    if (f.IsZombie()) { printf("[layer6Summary] no systematic_uncertainties.root — skip H2\n"); return; }
    TGraphErrors* gNom = dynamic_cast<TGraphErrors*>(f.Get("gSigT_nominal"));
    TGraph*       gSys = dynamic_cast<TGraph*>(f.Get("gSystTotal"));
    if (!gNom || !gSys) { printf("[layer6Summary] graphs missing — skip H2\n"); return; }

    const int n = gNom->GetN();
    TGraph* band = new TGraph(2 * n);
    double maxhi = 0.;
    for (int i = 0; i < n; ++i) {
        double x, y; gNom->GetPoint(i, x, y);
        double sx, sy; gSys->GetPoint(i, sx, sy);
        band->SetPoint(i, x, y + sy);
        band->SetPoint(2 * n - 1 - i, x, y - sy);
        if (y + sy > maxhi) maxhi = y + sy;
    }

    TCanvas* c = NewSquareCanvas("c_l6_band");
    c->cd();
    band->SetFillColorAlpha(kROrange, 0.30); band->SetLineColor(0);
    gNom->SetMarkerStyle(20); gNom->SetMarkerSize(1.4); gNom->SetMarkerColor(kRData);
    gNom->SetLineColor(kRData); gNom->SetLineWidth(2);

    gNom->Draw("AP");
    band->Draw("F SAME");
    gNom->Draw("P SAME");   // points on top of the band
    gNom->GetXaxis()->SetTitle("beam energy (GeV)");
    gNom->GetYaxis()->SetTitle("A^{2}-combo #sigma_{t}  (ps)");
    gNom->GetXaxis()->SetLimits(0., 165.);
    gNom->GetYaxis()->SetRangeUser(0., maxhi * 1.12);
    gNom->GetXaxis()->SetTitleSize(0.046);
    gNom->GetYaxis()->SetTitleSize(0.046);

    TLegend* leg = new TLegend(0.42, 0.80, 0.95, 0.90);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextFont(42); leg->SetTextSize(0.033);
    leg->AddEntry(gNom, "nominal (stat. error bars)", "lp");
    leg->AddEntry(band, "systematic band", "f");
    leg->Draw();

    { TLatex t; t.SetNDC(); t.SetTextSize(0.028); t.SetTextColor(kGray + 3);
      t.DrawLatex(0.40, 0.40, "Systematic is a few ps at every energy");
      t.DrawLatex(0.40, 0.360, "(robust core #sigma, stable at low stats).");
      t.DrawLatex(0.40, 0.300, "The MCP-free (DW#minusUP)/2 headline");
      t.DrawLatex(0.40, 0.260, "is more robust still (Layer 5)."); }

    DrawPageTitle("A^{2}-combo #sigma_{t} with systematic band");
    c->Print(Form("%slayer6_band.png", kSumDir));
    printf("[layer6Summary] wrote layer6_band.png\n");
}

} // namespace

void layer6Summary()
{
    ApplyRADiCALStyle();
    gSystem->mkdir(kSumDir, kTRUE);

    HeroBudget();
    HeroBand();

    printf("\n[layer6Summary] Done — hero figures in %s\n", kSumDir);
}

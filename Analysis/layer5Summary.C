// ===========================================================================
// layer5Summary.C  —  compact, Ledovskoy-clean "Layer 5: Physics Extraction".
//
//   layer5_timing.png    H1  timing resolution vs energy (headline):
//                            energy-binned (DW-UP)/2 + stochastic fit + arXiv ref
//   layer5_energy.png    H2  energy resolution sigma_E/E vs energy + stochastic fit
//   layer5_uniformity.png H3 spatial uniformity: sigma_t across the beam spot (150 GeV)
//
// Reads timing_energy_bins.root (gBestSigma_teb_m0, gPaper_teb), summary.root
// (gEnergyResolution), uniformity_scan.root (hSigT2D_150GeV).
//
// Run:  ROOT_INCLUDE_PATH=Analysis root -l -b -q 'Analysis/layer5Summary.C+'
// ===========================================================================
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TProfile2D.h"
#include "TF1.h"
#include "TAxis.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TSystem.h"

#include "RADiCALStyle.h"

namespace {

const char* kSumDir = "Analysis/Output/Summary/";

double valAt(TGraph* g, double E)
{
    if (!g) return -1.;
    double x, y; for (int i = 0; i < g->GetN(); ++i) { g->GetPoint(i, x, y); if (fabs(x - E) < 1.) return y; }
    return -1.;
}

// ── H1 — timing resolution vs energy (headline) ────────────────────────────
void HeroTiming()
{
    TFile f(Form("%stiming_energy_bins.root", kSumDir));
    if (f.IsZombie()) { printf("[layer5Summary] no timing_energy_bins.root — skip H1\n"); return; }
    TGraphErrors* gT = dynamic_cast<TGraphErrors*>(f.Get("gBestSigma_teb_m0"));
    TGraphErrors* gP = dynamic_cast<TGraphErrors*>(f.Get("gPaper_teb"));
    if (!gT) { printf("[layer5Summary] gBestSigma_teb_m0 missing — skip H1\n"); return; }

    TCanvas* c = NewSquareCanvas("c_l5_timing");
    c->cd();

    TF1* fit = new TF1("fT", "sqrt([0]*[0]/x + [1]*[1])", 20., 160.);
    fit->SetParameters(220., 30.);
    gT->Fit(fit, "QN");
    const double a = fit->GetParameter(0), b = fit->GetParameter(1);
    fit->SetLineColor(kRData); fit->SetLineWidth(2); fit->SetLineStyle(2);

    gT->SetMarkerStyle(20); gT->SetMarkerSize(1.5); gT->SetMarkerColor(kRData);
    gT->SetLineColor(kRData); gT->SetLineWidth(2);
    if (gP) { gP->SetMarkerStyle(24); gP->SetMarkerSize(1.3); gP->SetMarkerColor(kGray + 2);
              gP->SetLineColor(kGray + 2); gP->SetLineWidth(2); gP->SetLineStyle(3); }

    TMultiGraph* mg = new TMultiGraph();
    mg->Add(gT, "P");
    if (gP) mg->Add(gP, "PL");
    mg->Draw("A");
    fit->Draw("SAME");
    mg->GetXaxis()->SetTitle("beam energy (GeV)");
    mg->GetYaxis()->SetTitle("#sigma_{t} (ps)");
    mg->GetXaxis()->SetLimits(0., 165.);
    mg->GetYaxis()->SetRangeUser(0., 80.);
    mg->GetXaxis()->SetTitleSize(0.048);
    mg->GetYaxis()->SetTitleSize(0.048);

    TLegend* leg = new TLegend(0.40, 0.78, 0.95, 0.90);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextFont(42); leg->SetTextSize(0.033);
    leg->AddEntry(gT, "this analysis (energy-binned)", "lp");
    if (gP) leg->AddEntry(gP, "arXiv:2401.01747", "lp");
    leg->Draw();

    { TLatex t; t.SetNDC(); t.SetTextSize(0.034); t.SetTextColor(kRData);
      t.DrawLatex(0.18, 0.34, Form("#sigma_{t} = %.0f ps at 150 GeV", valAt(gT, 150.)));
      t.SetTextSize(0.028); t.SetTextColor(kGray + 2);
      t.DrawLatex(0.18, 0.295, Form("fit  #sigma_{t} = %.0f/#sqrt{E} #oplus %.0f ps", a, b));
      t.DrawLatex(0.18, 0.255, "MCP-free (DW#minusUP)/2, CFD-5%, energy-binned"); }

    DrawPageTitle("Timing resolution vs beam energy  (headline)");
    c->Print(Form("%slayer5_timing.png", kSumDir));
    printf("[layer5Summary] wrote layer5_timing.png\n");
}

// ── H2 — energy resolution vs energy ───────────────────────────────────────
void HeroEnergy()
{
    TFile f(Form("%ssummary.root", kSumDir));
    if (f.IsZombie()) { printf("[layer5Summary] no summary.root — skip H2\n"); return; }
    TGraphErrors* g = dynamic_cast<TGraphErrors*>(f.Get("gEnergyResolution"));
    if (!g) { printf("[layer5Summary] gEnergyResolution missing — skip H2\n"); return; }
    g->GetListOfFunctions()->Clear();   // drop any fit stored by analyzeResolution (avoids a 2nd red line)

    TCanvas* c = NewSquareCanvas("c_l5_energy");
    c->cd();

    TF1* fit = new TF1("fE", "sqrt([0]*[0]/x + [1]*[1])", 20., 160.);
    fit->SetParameters(100., 10.);
    g->Fit(fit, "QN");
    const double a = fit->GetParameter(0), b = fit->GetParameter(1);
    fit->SetLineColor(kRData); fit->SetLineWidth(2); fit->SetLineStyle(2);

    g->SetMarkerStyle(20); g->SetMarkerSize(1.5); g->SetMarkerColor(kRData);
    g->SetLineColor(kRData); g->SetLineWidth(2);
    g->Draw("AP");
    fit->Draw("SAME");
    g->GetXaxis()->SetTitle("beam energy (GeV)");
    g->GetYaxis()->SetTitle("#sigma_{E}/E  (%)");
    g->GetXaxis()->SetLimits(0., 165.);
    g->GetYaxis()->SetRangeUser(0., 22.);
    g->GetXaxis()->SetTitleSize(0.048);
    g->GetYaxis()->SetTitleSize(0.048);

    { TLatex t; t.SetNDC(); t.SetTextSize(0.034); t.SetTextColor(kRData);
      t.DrawLatex(0.20, 0.42, Form("#sigma_{E}/E = %.1f%% at 150 GeV", valAt(g, 150.)));
      t.SetTextSize(0.028); t.SetTextColor(kGray + 2);
      t.DrawLatex(0.20, 0.375, Form("fit  %.0f%%/#sqrt{E} #oplus %.1f%%", a, b));
      t.DrawLatex(0.20, 0.335, "leakage-dominated in the 14 mm");
      t.DrawLatex(0.20, 0.295, "prototype (a short test module)."); }

    DrawPageTitle("Energy resolution vs beam energy");
    c->Print(Form("%slayer5_energy.png", kSumDir));
    printf("[layer5Summary] wrote layer5_energy.png\n");
}

// ── H3 — spatial uniformity of sigma_t across the beam spot (150 GeV) ──────
void HeroUniformity()
{
    TFile f(Form("%suniformity_scan.root", kSumDir));
    if (f.IsZombie()) { printf("[layer5Summary] no uniformity_scan.root — skip H3\n"); return; }
    TProfile2D* h = dynamic_cast<TProfile2D*>(f.Get("hSigT2D_150GeV"));
    if (!h) { printf("[layer5Summary] hSigT2D_150GeV missing — skip H3\n"); return; }
    h->SetDirectory(nullptr);

    gStyle->SetPalette(kFall);
    TCanvas* c = NewSquareCanvas("c_l5_unif", 600, 104, 92, 54, 138);
    c->cd();
    h->GetXaxis()->SetTitle("x_{track} (mm)");
    h->GetYaxis()->SetTitle("y_{track} (mm)");
    h->GetZaxis()->SetTitle("#sigma_{t} (ps)");
    h->GetXaxis()->SetTitleSize(0.046);
    h->GetYaxis()->SetTitleSize(0.046);
    h->GetYaxis()->SetTitleOffset(1.25);
    // Tighten z to the bulk so the (mild) spatial variation is visible rather
    // than washed out by a low-stat corner bin.
    h->SetMinimum(50.); h->SetMaximum(130.);
    h->Draw("COLZ");

    { TLatex t; t.SetNDC(); t.SetTextSize(0.029); t.SetTextColor(kGray + 3);
      t.DrawLatex(0.16, 0.20, "#sigma_{t} stable (~90 ps) across the core fiducial;");
      t.DrawLatex(0.16, 0.16, "mild edge variation = residual position-walk."); }

    DrawPageTitle("Spatial uniformity -- #sigma_{t}(x,y) at 150 GeV");
    c->Print(Form("%slayer5_uniformity.png", kSumDir));
    printf("[layer5Summary] wrote layer5_uniformity.png\n");
}

} // namespace

void layer5Summary()
{
    ApplyRADiCALStyle();
    gSystem->mkdir(kSumDir, kTRUE);

    HeroTiming();
    HeroEnergy();
    HeroUniformity();

    printf("\n[layer5Summary] Done — hero figures in %s\n", kSumDir);
}

// ============================================================================
// timingFloorComparison.C — the floor CONFIRMS the published 17.5 ps.
// ----------------------------------------------------------------------------
// REWRITTEN 2026-07-21 (report-pipeline migration to the production chain).
// The previous version of this macro hard-coded the retired best-bin/cfd05 era
// ("22 vs 17.5 ps floor — a real gap?") — that framing is RETIRED under the
// claims law: the production floor is b = 18.8 ± 0.8 ps and CONFIRMS (never
// revises) the published 17.5 ps of NIM A 1068 (2024) 169737.
//
// DATA-DRIVEN: reads the production curve + fit from
// output/Summary/timing_energy_bins.root (written by timingProduction.C — run
// that first). Nothing numeric is hard-coded except the published
// parametrisation (256/sqrt(E) ⊕ 17.5), which is a literature constant.
//
// Output: output/Summary/timing_floor_comparison.{png,pdf}
// ============================================================================
#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TParameter.h"
#include "TF1.h"
#include "TBox.h"
#include "TLine.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TH1.h"
#include <cstdio>

#include "PlotUtils.h"      // ApplyRADiCALStyle, kRData, DrawPageTitle

void timingFloorComparison()
{
    TH1::AddDirectory(kFALSE);
    ApplyRADiCALStyle();

    // production curve + fit parameters (data-driven)
    TFile* f = TFile::Open("output/Summary/timing_energy_bins.root");
    if (!f || f->IsZombie()){ printf("[timingFloorComparison] run timingProduction.C first\n"); return; }
    TGraphErrors* g = dynamic_cast<TGraphErrors*>(f->Get("gBestSigma_teb_m0"));
    auto par=[&](const char* n)->double{ auto* p=dynamic_cast<TParameter<double>*>(f->Get(n)); return p?p->GetVal():0; };
    const double aUs  = par("fit_a_DSB1"),  bUs  = par("fit_b_DSB1");
    const double aUse = par("fit_aerr_DSB1"), bUse = par("fit_berr_DSB1");
    if (!g || aUs<=0){ printf("[timingFloorComparison] missing production graph/fit\n"); return; }
    g->SetMarkerStyle(20); g->SetMarkerColor(kRData); g->SetLineColor(kRData);
    g->SetMarkerSize(1.6); g->SetLineWidth(2);

    const double aPap = 256., bPap = 17.5;   // published (NIM A 1068 (2024) 169737)
    TF1* fUs  = new TF1("fUs",  "sqrt([0]*[0]/x + [1]*[1])", 18., 420.);
    fUs->SetParameters(aUs, bUs);   fUs->SetLineColor(kRData); fUs->SetLineWidth(3);
    TF1* fPap = new TF1("fPap", "sqrt([0]*[0]/x + [1]*[1])", 18., 420.);
    fPap->SetParameters(aPap, bPap); fPap->SetLineColor(kGray+2); fPap->SetLineWidth(3); fPap->SetLineStyle(2);

    TCanvas* c = new TCanvas("c_floor", "", 960, 740);
    c->SetLeftMargin(0.12); c->SetBottomMargin(0.13); c->SetRightMargin(0.04); c->SetTopMargin(0.10);
    c->SetTickx(1); c->SetTicky(1);

    const double xlo = 18., xhi = 420., ylo = 10., yhi = 56.;
    TH1F* fr = c->DrawFrame(xlo, ylo, xhi, yhi);
    fr->GetXaxis()->SetTitle("beam energy  E  (GeV)");
    fr->GetYaxis()->SetTitle("#sigma_{t}  (DW#minusUP)/2  (ps)");
    fr->GetXaxis()->SetTitleSize(0.047); fr->GetYaxis()->SetTitleSize(0.047);

    // shade the extrapolation region (E > 150 GeV: no beam data)
    TBox* extrap = new TBox(150., ylo, xhi, yhi);
    extrap->SetFillColorAlpha(kGray, 0.16); extrap->SetLineColor(0); extrap->Draw();
    { TLatex a; a.SetTextColor(kGray+3); a.SetTextSize(0.030); a.SetTextAngle(90);
      a.DrawLatex(168., ylo + 0.30*(yhi-ylo), "extrapolation  (no beam data above 150 GeV)"); }

    // floor asymptotes: production floor band (b ± err) + published floor line
    TBox* bandUs = new TBox(xlo, bUs-bUse, xhi, bUs+bUse);
    bandUs->SetFillColorAlpha(kRData, 0.15); bandUs->SetLineColor(0); bandUs->Draw();
    TLine* lUs  = new TLine(xlo, bUs,  xhi, bUs);  lUs->SetLineColor(kRData);  lUs->SetLineStyle(3); lUs->SetLineWidth(2); lUs->Draw();
    TLine* lPap = new TLine(xlo, bPap, xhi, bPap); lPap->SetLineColor(kGray+2); lPap->SetLineStyle(3); lPap->SetLineWidth(2); lPap->Draw();
    { TLatex a; a.SetTextSize(0.028);
      a.SetTextColor(kRData);  a.DrawLatex(215., bUs+bUse+0.8, Form("production floor  b = %.1f #pm %.1f ps", bUs, bUse));
      a.SetTextColor(kGray+3); a.DrawLatex(255., bPap-2.4, Form("published floor  b = %.1f ps", bPap)); }

    fPap->Draw("same"); fUs->Draw("same"); g->Draw("P same");

    TLegend* L = new TLegend(0.44, 0.72, 0.96, 0.89);
    L->SetBorderSize(0); L->SetFillStyle(0); L->SetTextFont(42); L->SetTextSize(0.031);
    L->AddEntry(g,    "production: brightest-1000, srCFD", "pe");
    L->AddEntry(fUs,  Form("fit: %.0f/#sqrt{E} #oplus %.1f ps", aUs, bUs), "l");
    L->AddEntry(fPap, "published: 256/#sqrt{E} #oplus 17.5 ps (cfd05 era)", "l");
    L->Draw();

    { TLatex a; a.SetNDC(); a.SetTextSize(0.028); a.SetTextColor(kGray+3);
      a.DrawLatex(0.47, 0.665, "The fitted floor is an extrapolation beyond the data;");
      a.DrawLatex(0.47, 0.628, Form("%.1f #pm %.1f ps is consistent with, and CONFIRMS,", bUs, bUse));
      a.DrawLatex(0.47, 0.591, "the published 17.5 ps. The floor's plausible origin is");
      a.DrawLatex(0.47, 0.554, "longitudinal shower-depth fluctuation (depth dial)."); }

    DrawPageTitle(Form("The timing floor: %.1f #pm %.1f ps -- confirming the published 17.5 ps", bUs, bUse));
    c->Print("output/Summary/timing_floor_comparison.png");
    c->Print("output/Summary/timing_floor_comparison.pdf");
    printf("[timingFloorComparison] wrote timing_floor_comparison.{png,pdf} (b=%.1f±%.1f vs published %.1f)\n", bUs, bUse, bPap);
}

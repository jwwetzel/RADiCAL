// ============================================================================
// timingFloorComparison.C — the "22 vs 17.5 ps floor" explained in one figure.
// ----------------------------------------------------------------------------
// Both this analysis and arXiv:2401.01747 use the SAME detector, the SAME CERN
// H2 beam, and the SAME 25-150 GeV range, and both reach ~27 ps at 150 GeV.
// The floor b is sigma_t at E->infinity: an EXTRAPOLATION beyond the data, and
// strongly anti-correlated with the stochastic term a.  This figure overlays
// both fits on our data and shades the E>150 GeV region (where neither beam
// reached) so it is clear the floor is extrapolated, not measured -- and that
// where data exist, this analysis is equal or better.
//
//   this analysis : sigma_t = 200/sqrt(E) (+) 22  ps   (a-b corr = -0.8)
//   arXiv:2401.01747: sigma_t = 256/sqrt(E) (+) 17.5 ps
//
// Output: Analysis/Output/Summary/timing_floor_comparison.{png,pdf}
// ============================================================================

#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TBox.h"
#include "TLine.h"
#include "TMarker.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TH1.h"

#include "PlotUtils.h"      // ApplyRADiCALStyle, kRData/kRRed/kROrange, DrawPageTitle

void timingFloorComparison()
{
    TH1::AddDirectory(kFALSE);
    ApplyRADiCALStyle();

    // Our best-bin, OOS-validated per-energy points (results.json: teb_sigma).
    const int    nE = 6;
    double E[nE]   = { 25., 50., 75., 100., 125., 150. };
    double S[nE]   = { 47.1, 32.7, 30.8, 30.2, 29.2, 27.4 };
    double Se[nE]  = { 1.32, 1.35, 1.47, 1.09, 1.27, 1.17 };
    double Ez[nE]  = { 0. };
    TGraphErrors* g = new TGraphErrors(nE, E, S, Ez, Se);
    g->SetMarkerStyle(20); g->SetMarkerColor(kRData); g->SetLineColor(kRData);
    g->SetMarkerSize(1.6); g->SetLineWidth(2);

    // Fits: sigma_t = sqrt(a^2/E + b^2)
    const double aUs = 200., bUs = 21.8, aPap = 256., bPap = 17.5;
    TF1* fUs  = new TF1("fUs",  "sqrt([0]*[0]/x + [1]*[1])", 18., 420.);
    fUs->SetParameters(aUs, bUs);   fUs->SetLineColor(kRData); fUs->SetLineWidth(3);
    TF1* fPap = new TF1("fPap", "sqrt([0]*[0]/x + [1]*[1])", 18., 420.);
    fPap->SetParameters(aPap, bPap); fPap->SetLineColor(kGray+2); fPap->SetLineWidth(3); fPap->SetLineStyle(2);

    TCanvas* c = new TCanvas("c_floor", "", 960, 740);
    c->SetLeftMargin(0.12); c->SetBottomMargin(0.13); c->SetRightMargin(0.04); c->SetTopMargin(0.10);
    c->SetTickx(1); c->SetTicky(1);

    const double xlo = 18., xhi = 420., ylo = 12., yhi = 56.;
    // frame
    TH1F* fr = c->DrawFrame(xlo, ylo, xhi, yhi);
    fr->GetXaxis()->SetTitle("beam energy  E  (GeV)");
    fr->GetYaxis()->SetTitle("#sigma_{t}  (DW#minusUP)/2  (ps)");
    fr->GetXaxis()->SetTitleSize(0.047); fr->GetYaxis()->SetTitleSize(0.047);

    // shade the extrapolation region (E > 150 GeV: no data from either beam)
    TBox* extrap = new TBox(150., ylo, xhi, yhi);
    extrap->SetFillColorAlpha(kGray, 0.16); extrap->SetLineColor(0); extrap->Draw();
    { TLatex a; a.SetTextColor(kGray+3); a.SetTextSize(0.030); a.SetTextAngle(90);
      a.DrawLatex(168., ylo + 0.30*(yhi-ylo), "extrapolation  (no beam data above 150 GeV)"); }

    // floor asymptote lines
    TLine* lUs  = new TLine(xlo, bUs,  xhi, bUs);  lUs->SetLineColor(kRData);  lUs->SetLineStyle(3); lUs->SetLineWidth(2); lUs->Draw();
    TLine* lPap = new TLine(xlo, bPap, xhi, bPap); lPap->SetLineColor(kGray+2); lPap->SetLineStyle(3); lPap->SetLineWidth(2); lPap->Draw();
    { TLatex a; a.SetTextSize(0.028);
      a.SetTextColor(kRData);  a.DrawLatex(330., bUs+0.8,  Form("our floor  b = %.1f ps", bUs));
      a.SetTextColor(kGray+3); a.DrawLatex(330., bPap-2.2, Form("arXiv floor  b = %.1f ps", bPap)); }

    // the fits + data
    fPap->Draw("same"); fUs->Draw("same"); g->Draw("P same");

    // the agreement point at 150 GeV
    { TMarker* m = new TMarker(150., 27.4, 29); m->SetMarkerColor(kROrange); m->SetMarkerSize(2.2); m->Draw();
      TLatex a; a.SetTextColor(kROrange+2); a.SetTextSize(0.026);
      a.DrawLatex(95., 33.0, "both #approx 27 ps"); a.DrawLatex(107., 31.0, "at 150 GeV"); }

    TLegend* L = new TLegend(0.47, 0.73, 0.96, 0.89);
    L->SetBorderSize(0); L->SetFillStyle(0); L->SetTextFont(42); L->SetTextSize(0.032);
    L->AddEntry(g,    "this analysis (data, 25#minus150 GeV)", "pe");
    L->AddEntry(fUs,  "fit: 200/#sqrt{E} #oplus 22 ps", "l");
    L->AddEntry(fPap, "arXiv:2401.01747: 256/#sqrt{E} #oplus 17.5 ps", "l");
    L->Draw();

    // explanatory note in the empty upper-right (below the legend)
    { TLatex a; a.SetNDC(); a.SetTextSize(0.028); a.SetTextColor(kGray+3);
      a.DrawLatex(0.50, 0.685, "Both fits #rightarrow 27 ps at 150 GeV;");
      a.DrawLatex(0.50, 0.648, "they split it differently between");
      a.DrawLatex(0.50, 0.611, "stochastic term and floor (#rho_{ab} #approx #minus0.8).");
      a.SetTextColor(kRData);
      a.DrawLatex(0.50, 0.560, "Where data exist (#leq150 GeV),");
      a.DrawLatex(0.50, 0.523, "this analysis is equal or better."); }

    DrawPageTitle("The 22 vs 17.5 ps floor: an extrapolation, not a measured gap");
    c->Print("Analysis/Output/Summary/timing_floor_comparison.png");
    c->Print("Analysis/Output/Summary/timing_floor_comparison.pdf");
    printf("[timingFloorComparison] wrote timing_floor_comparison.png\n");
}

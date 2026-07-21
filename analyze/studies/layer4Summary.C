// ===========================================================================
// layer4Summary.C  —  compact, Ledovskoy-clean "Layer 4: Calibration" heroes.
//
// The story: the calibration + combination steps turn a ~180 ps single channel
// into the 37 ps headline.
//
//   layer4_ladder.png   H1  resolution ladder vs energy: best single channel
//                           -> A^2-weighted 8-channel combo -> energy-binned
//                           (DW-UP)/2 headline
//   layer4_walk.png     H2  single-channel timing: plain CFD-20% vs the
//                           CFD + HG/LG-ratio walk correction
//
// Reads persisted graphs: timing_summary.root (gTiming_M0..M4), timing_methods
// .root (gBest_m*), timing_energy_bins.root (gBestSigma_teb_m0).
//
// Run:  ROOT_INCLUDE_PATH=Analysis root -l -b -q 'Analysis/layer4Summary.C+'
// ===========================================================================
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TAxis.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TSystem.h"

#include "RADiCALStyle.h"

namespace {

const char* kSumDir = "output/Summary/";

TGraph* getG(const char* file, const char* name)
{
    TFile* f = TFile::Open(Form("%s%s", kSumDir, file));
    if (!f || f->IsZombie()) { if (f) delete f; return nullptr; }
    TGraph* g = dynamic_cast<TGraph*>(f->Get(name));
    return g ? dynamic_cast<TGraph*>(g->Clone()) : nullptr;   // detached clone (file closes)
}

// ── H1 — resolution ladder: single -> combo -> headline ────────────────────
void HeroLadder()
{
    TGraph* gSingle = getG("timing_summary.root", "gTiming_M0");     // best single channel
    TGraph* gCombo  = getG("timing_summary.root", "gTiming_M3");     // A^2-weighted all-8
    TGraph* gHead   = getG("timing_energy_bins.root", "gBestSigma_teb_m0");  // headline
    if (!gSingle || !gCombo || !gHead) { printf("[layer4Summary] missing ladder inputs — skip H1\n"); return; }

    TCanvas* c = NewSquareCanvas("c_l4_ladder");
    c->cd();

    gSingle->SetMarkerStyle(25); gSingle->SetMarkerSize(1.3); gSingle->SetMarkerColor(kRRed);
    gSingle->SetLineColor(kRRed);    gSingle->SetLineWidth(2); gSingle->SetLineStyle(2);
    gCombo->SetMarkerStyle(21);  gCombo->SetMarkerSize(1.3);  gCombo->SetMarkerColor(kROrange);
    gCombo->SetLineColor(kROrange);  gCombo->SetLineWidth(2);
    gHead->SetMarkerStyle(20);   gHead->SetMarkerSize(1.4);   gHead->SetMarkerColor(kRData);
    gHead->SetLineColor(kRData);     gHead->SetLineWidth(3);

    TMultiGraph* mg = new TMultiGraph();
    mg->Add(gSingle, "PL");
    mg->Add(gCombo, "PL");
    mg->Add(gHead, "PL");
    mg->Draw("A");
    mg->GetXaxis()->SetTitle("beam energy (GeV)");
    mg->GetYaxis()->SetTitle("#sigma_{t} (ps)");
    mg->GetXaxis()->SetLimits(0., 165.);
    mg->GetYaxis()->SetRangeUser(0., 210.);
    mg->GetXaxis()->SetTitleSize(0.048);
    mg->GetYaxis()->SetTitleSize(0.048);

    TLegend* leg = new TLegend(0.34, 0.74, 0.95, 0.90);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextFont(42); leg->SetTextSize(0.033);
    leg->AddEntry(gSingle, "best single channel", "lp");
    leg->AddEntry(gCombo,  "A^{2}-weighted 8-channel combo", "lp");
    leg->AddEntry(gHead,   "(DW#minusUP)/2 energy-binned (headline)", "lp");
    leg->Draw();

    { TLatex t; t.SetNDC(); t.SetTextSize(0.030); t.SetTextColor(kGray + 3);
      t.DrawLatex(0.17, 0.63, "Combining the 8 channels cuts #sigma_{t} ~3#times;");
      t.DrawLatex(0.17, 0.585, "energy-binning + the corner difference");
      t.DrawLatex(0.17, 0.540, "then reach the 37 ps headline."); }

    DrawPageTitle("Resolution ladder -- from single channel to the headline");
    c->Print(Form("%slayer4_ladder.png", kSumDir));
    printf("[layer4Summary] wrote layer4_ladder.png\n");
}

// ── H2 — walk correction: plain CFD vs CFD + HG/LG-ratio correction ────────
void HeroWalk()
{
    TGraph* gCFD  = getG("timing_methods.root", "gBest_m1");   // CFD 20%, no walk corr
    TGraph* gWalk = getG("timing_methods.root", "gBest_m7");   // CFD + HG/LG ratio walk corr
    if (!gCFD || !gWalk) { printf("[layer4Summary] missing walk inputs — skip H2\n"); return; }

    TCanvas* c = NewSquareCanvas("c_l4_walk");
    c->cd();

    gCFD->SetMarkerStyle(25);  gCFD->SetMarkerSize(1.3);  gCFD->SetMarkerColor(kRRed);
    gCFD->SetLineColor(kRRed);   gCFD->SetLineWidth(2); gCFD->SetLineStyle(2);
    gWalk->SetMarkerStyle(20); gWalk->SetMarkerSize(1.3); gWalk->SetMarkerColor(kRData);
    gWalk->SetLineColor(kRData); gWalk->SetLineWidth(2);

    TMultiGraph* mg = new TMultiGraph();
    mg->Add(gCFD, "PL");
    mg->Add(gWalk, "PL");
    mg->Draw("A");
    mg->GetXaxis()->SetTitle("beam energy (GeV)");
    mg->GetYaxis()->SetTitle("best single-channel #sigma_{t} (ps)");
    mg->GetXaxis()->SetLimits(0., 165.);
    mg->GetYaxis()->SetRangeUser(0., 240.);
    mg->GetXaxis()->SetTitleSize(0.046);
    mg->GetYaxis()->SetTitleSize(0.046);

    TLegend* leg = new TLegend(0.34, 0.795, 0.95, 0.90);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextFont(42); leg->SetTextSize(0.034);
    leg->AddEntry(gCFD,  "plain CFD (no walk correction)", "lp");
    leg->AddEntry(gWalk, "CFD + HG/LG-ratio walk correction", "lp");
    leg->Draw();

    { TLatex t; t.SetNDC(); t.SetTextSize(0.029); t.SetTextColor(kGray + 3);
      t.DrawLatex(0.17, 0.40, "The HG/LG-ratio correction removes the");
      t.DrawLatex(0.17, 0.355, "amplitude-dependent time-walk, improving");
      t.DrawLatex(0.17, 0.310, "single-channel timing by 20-50 ps.");
      t.SetTextSize(0.026); t.SetTextColor(kGray + 2);
      t.DrawLatex(0.17, 0.255, "(Headline uses srCFD #minus the 5% edge with LG");
      t.DrawLatex(0.17, 0.220, "clip recovery; lower edge = less walk.)"); }

    DrawPageTitle("Walk correction -- single-channel timing");
    c->Print(Form("%slayer4_walk.png", kSumDir));
    printf("[layer4Summary] wrote layer4_walk.png\n");
}

} // namespace

void layer4Summary()
{
    ApplyRADiCALStyle();
    gSystem->mkdir(kSumDir, kTRUE);

    HeroLadder();
    HeroWalk();

    printf("\n[layer4Summary] Done — hero figures in %s\n", kSumDir);
}

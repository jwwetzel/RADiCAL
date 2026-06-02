// ============================================================================
// idealUniform.C — "what if every capillary followed the trend?" projection.
// ----------------------------------------------------------------------------
// Two Up capillaries buck the single-channel sigma_t(E) trend: NW-U (the
// consistently weakest channel) drifts UP with energy, and SW-U (the lone
// MCP2-referenced channel, gain ~0.74) stays flat.  Both are averaged into the
// headline (DW-UP)/2, so they cap it slightly.
//
// IMPORTANT (measured, not modelled): dropping them makes the headline WORSE,
// not better — the 8-channel corner average benefits from every channel
// (full-ntuple reconstruction at 150 GeV: all-4-Up 67 ps -> drop-both 101 ps).
// So they are a modest handicap, not a liability to remove.
//
// This figure projects the *intrinsic detector potential*: replace NW-U & SW-U
// with the well-behaved Up trend (mean of NE-U, SE-U) and propagate to the
// headline via the effective combination model sigma_HL ~ (1/8) sqrt(sum_i
// sigma_i^2) — validated at high E against the reconstruction (69 vs 67 ps at
// 150).  The measured 27.4 ps headline STANDS; this is a clearly-labelled
// projection only.
//
// Output: Analysis/Output/Summary/ideal_uniform_projection.{png,pdf}
// ============================================================================
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TFile.h"
#include "PlotUtils.h"

void idealUniform()
{
    ApplyRADiCALStyle();
    const char* nm[8]={"NW-D","NE-D","SE-D","SW-D","NW-U","NE-U","SE-U","SW-U"};
    double E[6]={25,50,75,100,125,150};
    double real[6]={47.1,32.7,30.8,30.2,29.2,27.4};   // measured headline (best-bin)

    TFile f("Analysis/Output/Summary/summary.root");
    if (f.IsZombie()) { printf("[idealUniform] no summary.root — skip\n"); return; }
    double S[8][6];
    for(int i=0;i<8;++i){ auto g=(TGraphErrors*)f.Get(Form("gTRes_%s",nm[i]));
        for(int e=0;e<6;++e){ double v=0; if(g)for(int p=0;p<g->GetN();++p){double x,y;g->GetPoint(p,x,y);if(fabs(x-E[e])<2)v=y;} S[i][e]=v; } }

    double good[6], idealHL[6], dHL[6];
    for(int e=0;e<6;++e){ good[e]=0.5*(S[5][e]+S[6][e]);
        double sa=0,si=0; for(int i=0;i<8;++i){ double a=S[i][e],id=(i==4||i==7)?good[e]:a; sa+=a*a; si+=id*id; }
        idealHL[e]=real[e]*sqrt(si/sa); dHL[e]=real[e]-idealHL[e]; }

    TCanvas* c=new TCanvas("c_ideal","",1340,600); c->Divide(2,1,0.015,0.02);

    // LEFT: per-channel single-channel sigma_t with the ideal Up trend overlaid
    c->cd(1); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13); gPad->SetTopMargin(0.10);
    TH1F* fr=gPad->DrawFrame(0,140,165,400);
    fr->GetXaxis()->SetTitle("beam energy (GeV)"); fr->GetYaxis()->SetTitle("single-channel #sigma_{t} (ps)");
    fr->GetXaxis()->SetTitleSize(0.047); fr->GetYaxis()->SetTitleSize(0.047);
    for(int i=0;i<8;++i){ TGraph* g=new TGraph(6,E,S[i]); int col=kRChannelCols[i];
        g->SetLineColor(col); g->SetMarkerColor(col);
        bool bad=(i==4||i==7);
        g->SetMarkerStyle(bad?20:24); g->SetMarkerSize(bad?1.3:0.8);
        g->SetLineWidth(bad?2:1); g->SetLineStyle(bad?1:3); g->Draw("PL"); }
    { TGraph* gi=new TGraph(6,E,good); gi->SetLineColor(kRGreen+1); gi->SetLineStyle(2);
      gi->SetLineWidth(3); gi->Draw("L"); }
    { TLatex t; t.SetTextSize(0.029); t.SetTextColor(kGray+3);
      t.DrawLatex(64,378,"NW-U (#bullet) drifts up, SW-U (#bullet) stays flat;");
      t.DrawLatex(64,364,"the other six (#circ) fall with energy.");
      t.SetTextColor(kRGreen+2);
      t.DrawLatex(64,206,"#it{ideal}: put NW-U & SW-U on the");
      t.DrawLatex(64,192,"#it{good-Up trend (green dashed)}"); }
    DrawPadTitle("Single-channel #sigma_{t}: two capillaries buck the trend", 0.055);

    // RIGHT: measured headline vs ideal-uniform projection
    c->cd(2); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13); gPad->SetTopMargin(0.10);
    TH1F* fr2=gPad->DrawFrame(0,20,165,52);
    fr2->GetXaxis()->SetTitle("beam energy (GeV)"); fr2->GetYaxis()->SetTitle("headline #sigma_{t}  (DW#minusUP)/2  (ps)");
    fr2->GetXaxis()->SetTitleSize(0.047); fr2->GetYaxis()->SetTitleSize(0.047);
    TGraph* gR=new TGraph(6,E,real); gR->SetLineColor(kRData); gR->SetMarkerColor(kRData);
    gR->SetMarkerStyle(20); gR->SetMarkerSize(1.5); gR->SetLineWidth(3); gR->Draw("PL");
    TGraph* gI=new TGraph(6,E,idealHL); gI->SetLineColor(kRGreen+1); gI->SetMarkerColor(kRGreen+1);
    gI->SetMarkerStyle(24); gI->SetMarkerSize(1.4); gI->SetLineWidth(3); gI->SetLineStyle(2); gI->Draw("PL");
    TLegend* L=new TLegend(0.30,0.76,0.96,0.89); L->SetBorderSize(0); L->SetFillStyle(0); L->SetTextSize(0.033);
    L->AddEntry(gR,"measured headline (all 8 channels)","pl");
    L->AddEntry(gI,"ideal-uniform projection","pl"); L->Draw();
    { TLatex t; t.SetTextSize(0.029); t.SetTextColor(kGray+3);
      t.DrawLatex(40,25.2,Form("125 GeV: %.1f #rightarrow %.1f ps  (#minus%.1f)",real[4],idealHL[4],dHL[4]));
      t.DrawLatex(40,23.7,Form("150 GeV: %.1f #rightarrow %.1f ps  (#minus%.1f)",real[5],idealHL[5],dHL[5]));
      t.SetTextColor(kROrange+2); t.SetTextSize(0.025);
      t.DrawLatex(40,21.9,"projection only #minus the measured 27.4 ps stands as the result"); }
    DrawPadTitle("Projected headline if all channels were uniform", 0.055);

    c->cd(0);
    DrawPageTitle("Detector-potential projection: what an all-uniform calorimeter would time");
    c->Print("Analysis/Output/Summary/ideal_uniform_projection.png");
    c->Print("Analysis/Output/Summary/ideal_uniform_projection.pdf");
    printf("[idealUniform] wrote ideal_uniform_projection.png  (150 GeV: 27.4 -> %.1f ps)\n", idealHL[5]);
}

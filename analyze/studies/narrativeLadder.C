// ============================================================================
// narrativeLadder.C — the payoff figures, with the FIXED estimators (no equal-
// width tail-starvation bias). Produces TWO plots:
//   narrative_ladder.png  — cfd05 vs hg_lgcfd (DW-UP)/2, brightest-K estimator,
//                           a/sqrtE (+) b fits. Monotonic; preserves the
//                           best-achievable (~thin bright slice) magnitude.
//   narrative_methods.png — hg_lgcfd read TWO ways: quantile best-bin (broad,
//                           conservative) vs brightest-K (thin, best-case). The
//                           band between is the resolution's selection range.
//   source setup.sh; root -l -b -q 'analyze/studies/narrativeLadder.C+("DSB1")'
// ============================================================================
#include "RadView.h"
#include "RadTiming.h"
#include "DataPaths.h"
#include "PlotUtils.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include <vector>
#include <cstdio>
#include <cmath>
using namespace rad;

static double fitFloor(TGraph* g, int col, int sty, int wid){
    TF1* f=new TF1(Form("f%d",col),"sqrt([0]*[0]/x+[1]*[1])",20,160); f->SetParameters(200,20);
    g->Fit(f,"Q"); f->SetLineColor(col); f->SetLineStyle(sty); f->SetLineWidth(wid); f->Draw("SAME");
    return std::fabs(f->GetParameter(1));
}

void narrativeLadder(const char* build="DSB1", int K=1000){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    BuildConfig cfg=BuildConfig::Load(Form("data/2023/configs/%s.json",build));
    const double Es[]={25,50,75,100,125,150};
    TGraph *g5=new TGraph,*gLk=new TGraph,*gLq=new TGraph; // cfd05-K, lgcfd-K, lgcfd-quantile
    printf("\n%s  (DW-UP)/2 sigma_t [ps]:  E | cfd05(K) lgcfd(K) | lgcfd(quantile)\n",build);
    for(double E:Es){
        TFile* fp=TFile::Open(radReduced(build,E)); if(!fp||fp->IsZombie())continue;
        TTree* t=(TTree*)fp->Get("rad"); RadView v; v.attach(t,&cfg);
        double s5=timingBrightestK(v,E,RadView::kCFD05,K).sigma_ps;
        double sLk=timingBrightestK(v,E,RadView::kLGCFD,K).sigma_ps;
        double sLq=timingBestBin(v,E,RadView::kLGCFD).sigma_ps;
        printf("   %3.0f |  %5.1f  %5.1f  |  %5.1f\n",E,s5,sLk,sLq);
        if(s5>0)g5->SetPoint(g5->GetN(),E,s5); if(sLk>0)gLk->SetPoint(gLk->GetN(),E,sLk); if(sLq>0)gLq->SetPoint(gLq->GetN(),E,sLq);
        fp->Close();
    }
    // ---- Figure A: cfd05 vs lgcfd (brightest-K), monotonic ----
    { TCanvas* c=new TCanvas("cA","",880,640); c->SetLeftMargin(0.13); c->SetRightMargin(0.05); c->SetGridy();
      g5->SetTitle(""); g5->SetMarkerStyle(24); g5->SetMarkerColor(kRed+1); g5->SetMarkerSize(1.6);
      g5->GetXaxis()->SetTitle("beam energy (GeV)"); g5->GetYaxis()->SetTitle("(DW#minusUP)/2  #sigma_{t} (ps)");
      g5->GetXaxis()->SetLimits(0,160); g5->GetYaxis()->SetRangeUser(0,60); g5->Draw("AP");
      double b5=fitFloor(g5,kRed+1,2,2);
      gLk->SetMarkerStyle(20); gLk->SetMarkerColor(kAzure+2); gLk->SetMarkerSize(1.6); gLk->Draw("P SAME");
      double bL=fitFloor(gLk,kAzure+2,1,3);
      TLegend* lg=new TLegend(0.42,0.72,0.93,0.88); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.034);
      lg->AddEntry(g5,Form("cfd05 (clipped-peak foot): floor %.0f ps",b5),"pl");
      lg->AddEntry(gLk,Form("hg_lgcfd (true-peak edge): floor %.0f ps",bL),"pl"); lg->Draw();
      TLatex tt; tt.SetNDC(); tt.SetTextFont(62); tt.SetTextSize(0.038);
      tt.DrawLatex(0.13,0.94,Form("%s timing (de-biased, monotonic): lgcfd best %.1f ps @150",build,gLk->GetY()[gLk->GetN()-1]));
      gSystem->mkdir("figures/narrative",kTRUE); c->Print("figures/narrative/narrative_ladder.png"); }
    // ---- Figure B: the two estimators (selection band) ----
    { TCanvas* c=new TCanvas("cB","",880,640); c->SetLeftMargin(0.13); c->SetRightMargin(0.05); c->SetGridy();
      gLq->SetTitle(""); gLq->SetMarkerStyle(21); gLq->SetMarkerColor(kViolet+1); gLq->SetMarkerSize(1.5);
      gLq->GetXaxis()->SetTitle("beam energy (GeV)"); gLq->GetYaxis()->SetTitle("hg_lgcfd (DW#minusUP)/2  #sigma_{t} (ps)");
      gLq->GetXaxis()->SetLimits(0,160); gLq->GetYaxis()->SetRangeUser(0,60); gLq->Draw("AP");
      double bq=fitFloor(gLq,kViolet+1,2,2);
      TGraph* gLk2=new TGraph(*gLk); gLk2->SetMarkerStyle(20); gLk2->SetMarkerColor(kAzure+2); gLk2->SetMarkerSize(1.5); gLk2->Draw("P SAME");
      double bk=fitFloor(gLk2,kAzure+2,1,3);
      TLegend* lg=new TLegend(0.40,0.72,0.93,0.88); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.034);
      lg->AddEntry(gLq,Form("quantile best-bin (typical bright): floor %.0f ps",bq),"pl");
      lg->AddEntry(gLk2,Form("brightest-%d slice (best case): floor %.0f ps",K,bk),"pl"); lg->Draw();
      TLatex tt; tt.SetNDC(); tt.SetTextFont(62); tt.SetTextSize(0.036);
      tt.DrawLatex(0.13,0.94,"hg_lgcfd, two selections: both monotonic; the band is the resolution range");
      c->Print("figures/narrative/narrative_methods.png"); }
    printf("  wrote figures/narrative/narrative_ladder.png + narrative_methods.png\n");
}

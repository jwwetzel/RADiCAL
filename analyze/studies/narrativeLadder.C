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
    BuildConfig cfg=BuildConfig::Load(radConfig(build).Data());
    const double Es[]={25,50,75,100,125,150};
    TGraph *g5q=new TGraph,*gLq=new TGraph,*gLk=new TGraph; // cfd05-quantile, lgcfd-quantile, lgcfd-K
    printf("\n%s  (DW-UP)/2 sigma_t [ps]:  E | cfd05(quant) lgcfd(quant) | lgcfd(brightest-%d)\n",build,K);
    for(double E:Es){
        TFile* fp=TFile::Open(radReduced(build,E)); if(!fp||fp->IsZombie())continue;
        TTree* t=(TTree*)fp->Get("rad"); RadView v; v.attach(t,&cfg);
        double s5q=timingBestBin(v,E,RadView::kCFD05).sigma_ps;   // OOS-robust headline estimator
        double sLq=timingBestBin(v,E,RadView::kLGCFD).sigma_ps;
        double sLk=timingBrightestK(v,E,RadView::kLGCFD,K).sigma_ps;
        printf("   %3.0f |    %5.1f      %5.1f     |   %5.1f\n",E,s5q,sLq,sLk);
        if(s5q>0)g5q->SetPoint(g5q->GetN(),E,s5q); if(sLq>0)gLq->SetPoint(gLq->GetN(),E,sLq); if(sLk>0)gLk->SetPoint(gLk->GetN(),E,sLk);
        fp->Close();
    }
    // ---- Figure A: cfd05 vs lgcfd, BOTH quantile (OOS-robust headline) ----
    { TCanvas* c=new TCanvas("cA","",880,640); c->SetLeftMargin(0.13); c->SetRightMargin(0.05); c->SetGridy();
      g5q->SetTitle(""); g5q->SetMarkerStyle(24); g5q->SetMarkerColor(kRed+1); g5q->SetMarkerSize(1.6);
      g5q->GetXaxis()->SetTitle("beam energy (GeV)"); g5q->GetYaxis()->SetTitle("(DW#minusUP)/2  #sigma_{t} (ps)");
      g5q->GetXaxis()->SetLimits(0,160); g5q->GetYaxis()->SetRangeUser(0,60); g5q->Draw("AP");
      double b5=fitFloor(g5q,kRed+1,2,2);
      gLq->SetMarkerStyle(20); gLq->SetMarkerColor(kAzure+2); gLq->SetMarkerSize(1.6); gLq->Draw("P SAME");
      double bL=fitFloor(gLq,kAzure+2,1,3);
      TLegend* lg=new TLegend(0.40,0.72,0.93,0.88); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.034);
      lg->AddEntry(g5q,Form("cfd05 (clipped-peak foot): floor %.0f ps",b5),"pl");
      lg->AddEntry(gLq,Form("hg_lgcfd (true-peak edge): floor %.0f ps",bL),"pl"); lg->Draw();
      TLatex tt; tt.SetNDC(); tt.SetTextFont(62); tt.SetTextSize(0.036);
      tt.DrawLatex(0.13,0.94,Form("%s: equal-population best-bin (OOS-robust) -- lgcfd ~11%% below cfd05",build));
      gSystem->mkdir("figures/narrative",kTRUE); c->Print("figures/narrative/narrative_ladder.png"); }
    // ---- Figure B: lgcfd two selections (quantile robust vs brightest-K best-case) ----
    { TCanvas* c=new TCanvas("cB","",880,640); c->SetLeftMargin(0.13); c->SetRightMargin(0.05); c->SetGridy();
      gLq->SetTitle(""); gLq->SetMarkerStyle(21); gLq->SetMarkerColor(kViolet+1); gLq->SetMarkerSize(1.5);
      gLq->GetXaxis()->SetTitle("beam energy (GeV)"); gLq->GetYaxis()->SetTitle("hg_lgcfd (DW#minusUP)/2  #sigma_{t} (ps)");
      gLq->GetXaxis()->SetLimits(0,160); gLq->GetYaxis()->SetRangeUser(0,60); gLq->Draw("AP");
      double bq=fitFloor(gLq,kViolet+1,2,2);
      gLk->SetMarkerStyle(20); gLk->SetMarkerColor(kAzure+2); gLk->SetMarkerSize(1.5); gLk->Draw("P SAME");
      double bk=fitFloor(gLk,kAzure+2,1,3);
      TLegend* lg=new TLegend(0.34,0.72,0.94,0.88); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.033);
      lg->AddEntry(gLq,Form("quantile best-bin -- broad/typical bright: floor %.0f ps",bq),"pl");
      lg->AddEntry(gLk,Form("brightest slice -- best-contained showers: floor %.0f ps",bk),"pl"); lg->Draw();
      TLatex tt; tt.SetNDC(); tt.SetTextFont(62); tt.SetTextSize(0.034);
      tt.DrawLatex(0.13,0.94,"hg_lgcfd: both selections OOS-stable -- 26 ps (brightest) to 30 ps (typical)");
      c->Print("figures/narrative/narrative_methods.png"); }
    printf("  wrote narrative_ladder.png + narrative_methods.png\n");
}

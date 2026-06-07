// ============================================================================
// narrativeLadder.C — the payoff figure: (DW-UP)/2 best-bin sigma_t(E),
// cfd05 (clipped-peak foot) vs hg_lgcfd (LG-predicted true-peak edge), each fit
// to sigma = a/sqrt(E) (+) b. Uses the canonical rad::timingBestBin pipeline.
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
using namespace rad;

void narrativeLadder(const char* build="DSB1"){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    BuildConfig cfg = BuildConfig::Load(Form("data/2023/configs/%s.json",build));
    const double Es[]={25,50,75,100,125,150};
    TGraph* g5=new TGraph(); TGraph* gL=new TGraph(); double bestLG=1e9,bestE=0;
    printf("\n%s  (DW-UP)/2 best-bin sigma_t [ps]:  E   cfd05   lgcfd\n",build);
    for(double E:Es){
        TFile* fp=TFile::Open(radReduced(build,E)); if(!fp||fp->IsZombie()) continue;
        TTree* t=(TTree*)fp->Get("rad"); RadView v; v.attach(t,&cfg);
        TimingResult r5=timingBestBin(v,E,RadView::kCFD05);
        TimingResult rL=timingBestBin(v,E,RadView::kLGCFD);
        printf("   %3.0f   %5.1f   %5.1f\n",E,r5.sigma_ps,rL.sigma_ps);
        if(r5.sigma_ps>0) g5->SetPoint(g5->GetN(),E,r5.sigma_ps);
        if(rL.sigma_ps>0){ gL->SetPoint(gL->GetN(),E,rL.sigma_ps); if(rL.sigma_ps<bestLG){bestLG=rL.sigma_ps;bestE=E;} }
        fp->Close();
    }
    TF1* f5=new TF1("f5","sqrt([0]*[0]/x+[1]*[1])",20,160); f5->SetParameters(200,30);
    TF1* fL=new TF1("fL","sqrt([0]*[0]/x+[1]*[1])",20,160); fL->SetParameters(200,20);
    g5->Fit(f5,"Q"); gL->Fit(fL,"Q");
    double b5=std::fabs(f5->GetParameter(1)), bL=std::fabs(fL->GetParameter(1));

    TCanvas* c=new TCanvas("c_lad","",880,640); c->SetLeftMargin(0.13); c->SetRightMargin(0.05); c->SetGridy();
    g5->SetTitle(""); g5->SetMarkerStyle(24); g5->SetMarkerColor(kRed+1); g5->SetMarkerSize(1.6); g5->SetLineWidth(0);
    g5->GetXaxis()->SetTitle("beam energy (GeV)"); g5->GetYaxis()->SetTitle("(DW#minusUP)/2  best-bin  #sigma_{t} (ps)");
    g5->GetXaxis()->SetLimits(0,160); g5->GetYaxis()->SetRangeUser(0,60); g5->Draw("AP");
    f5->SetLineColor(kRed+1); f5->SetLineStyle(2); f5->SetLineWidth(2); f5->Draw("SAME");
    gL->SetMarkerStyle(20); gL->SetMarkerColor(kAzure+2); gL->SetMarkerSize(1.6); gL->Draw("P SAME");
    fL->SetLineColor(kAzure+2); fL->SetLineWidth(3); fL->Draw("SAME");
    // mark the best lgcfd bin
    TLatex bm; bm.SetTextFont(62); bm.SetTextSize(0.034); bm.SetTextColor(kAzure+3);
    bm.DrawLatex(bestE-22,bestLG-5.5,Form("best bin: %.1f ps @ %.0f GeV",bestLG,bestE));
    TLegend* lg=new TLegend(0.45,0.70,0.93,0.88); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.036);
    lg->AddEntry(f5,Form("cfd05  (clipped-peak foot): floor b = %.0f ps",b5),"pl");
    lg->AddEntry(fL,Form("hg_lgcfd  (true-peak edge): floor b = %.0f ps",bL),"pl");
    lg->Draw();
    TLatex tt; tt.SetNDC(); tt.SetTextFont(62); tt.SetTextSize(0.040);
    tt.DrawLatex(0.13,0.94,Form("%s timing: the HG/LG-ratio method drops the floor %.0f #rightarrow %.0f ps",build,b5,bL));
    gSystem->mkdir("figures/narrative",kTRUE); c->Print("figures/narrative/narrative_ladder.png");
    printf("  cfd05 floor b=%.1f ps ; lgcfd floor b=%.1f ps ; lgcfd best bin %.1f ps @%.0f\n",b5,bL,bestLG,bestE);
    printf("  wrote figures/narrative/narrative_ladder.png\n");
}

// ============================================================================
// floorFit.C — is the (DW-UP)/2 "floor" b physical or a fit artifact? Fit the
// quantile and brightest-slice lgcfd ladders (a) FREELY (a,b each) and (b) with
// a SHARED b (common depth floor, per-selection a). If the shared-floor fit is
// nearly as good, the 20-vs-23 split is the a<->b anti-correlation, not a real
// difference: the floors converge and the selection lives in the stochastic term.
//   source setup.sh; root -l -b -q 'analyze/studies/floorFit.C+("DSB1")'
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
#include "TLine.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include "FigPaths.h"
#include <vector>
#include <cmath>
#include <cstdio>
using namespace rad;

void floorFit(const char* build="DSB1", int K=1000){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    BuildConfig cfg=BuildConfig::Load(radConfig(build).Data());
    const double Es[]={25,50,75,100,125,150};
    TGraph* gq=new TGraph; TGraph* gb=new TGraph;
    for(double E:Es){ TFile* fp=TFile::Open(radReduced(build,E)); if(!fp||fp->IsZombie())continue;
        TTree* t=(TTree*)fp->Get("rad"); RadView v; v.attach(t,&cfg);
        double sq=timingBestBin(v,E,RadView::kLGCFD).sigma_ps;
        double sb=timingBrightestK(v,E,RadView::kLGCFD,K).sigma_ps;
        if(sq>0) gq->SetPoint(gq->GetN(),E,sq); if(sb>0) gb->SetPoint(gb->GetN(),E,sb); fp->Close(); }

    // (a) free fits
    TF1 fq("fq","sqrt([0]*[0]/x+[1]*[1])",20,160); fq.SetParameters(200,22);
    TF1 fb("fb","sqrt([0]*[0]/x+[1]*[1])",20,160); fb.SetParameters(200,20);
    gq->Fit(&fq,"Q0"); gb->Fit(&fb,"Q0");
    printf("\n%s lgcfd — FREE fits:\n",build);
    printf("  quantile : a=%.0f  b=%.1f  chi2/ndf=%.2f\n",fq.GetParameter(0),std::fabs(fq.GetParameter(1)),fq.GetChisquare()/fq.GetNDF());
    printf("  brightest: a=%.0f  b=%.1f  chi2/ndf=%.2f\n",fb.GetParameter(0),std::fabs(fb.GetParameter(1)),fb.GetChisquare()/fb.GetNDF());

    // (b) shared-b fits: scan a common b, refit a for each, total chi2
    double bestB=0, bestChi=1e18, aqAt=0, abAt=0;
    for(double bc=14; bc<=30; bc+=0.25){
        TF1 hq("hq","sqrt([0]*[0]/x+[1]*[1])",20,160); hq.FixParameter(1,bc); hq.SetParameter(0,200);
        TF1 hb("hb","sqrt([0]*[0]/x+[1]*[1])",20,160); hb.FixParameter(1,bc); hb.SetParameter(0,200);
        gq->Fit(&hq,"Q0"); gb->Fit(&hb,"Q0");
        double chi=hq.GetChisquare()+hb.GetChisquare();
        if(chi<bestChi){ bestChi=chi; bestB=bc; aqAt=hq.GetParameter(0); abAt=hb.GetParameter(0);} }
    int ndf=(gq->GetN()-1)+(gb->GetN()-1);   // each loses 1 dof (a); b is shared
    printf("  SHARED-b best: b=%.1f ps  (a_quant=%.0f, a_bright=%.0f)  chi2/ndf=%.2f\n",bestB,aqAt,abAt,bestChi/ndf);
    printf("  => if shared-b chi2/ndf ~ free, the floors CONVERGE (~%.0f ps); the 20-vs-23 split is the a<->b degeneracy.\n",bestB);

    // figure: both ladders + the shared-floor fits
    TCanvas* c=new TCanvas("cf","",880,640); c->SetLeftMargin(0.12); c->SetRightMargin(0.05); c->SetGridy();
    gq->SetTitle(""); gq->SetMarkerStyle(21); gq->SetMarkerColor(kViolet+1); gq->SetMarkerSize(1.5);
    gq->GetXaxis()->SetTitle("beam energy (GeV)"); gq->GetYaxis()->SetTitle("hg_lgcfd (DW#minusUP)/2  #sigma_{t} (ps)");
    gq->GetXaxis()->SetLimits(0,200); gq->GetYaxis()->SetRangeUser(0,55); gq->Draw("AP");
    gb->SetMarkerStyle(20); gb->SetMarkerColor(kAzure+2); gb->SetMarkerSize(1.5); gb->Draw("P SAME");
    TF1* sq=new TF1("sq","sqrt([0]*[0]/x+[1]*[1])",20,200); sq->SetParameters(aqAt,bestB); sq->SetLineColor(kViolet+1); sq->SetLineWidth(2); sq->Draw("SAME");
    TF1* sb=new TF1("sb","sqrt([0]*[0]/x+[1]*[1])",20,200); sb->SetParameters(abAt,bestB); sb->SetLineColor(kAzure+2); sb->SetLineWidth(3); sb->Draw("SAME");
    TLine* fl=new TLine(0,bestB,200,bestB); fl->SetLineColor(kGray+2); fl->SetLineStyle(2); fl->Draw();
    TLatex tx; tx.SetTextColor(kGray+3); tx.SetTextSize(0.034); tx.DrawLatex(150,bestB+1.5,Form("shared floor %.0f ps",bestB));
    TLegend* lg=new TLegend(0.40,0.74,0.93,0.89); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.033);
    lg->AddEntry(gq,Form("quantile (typical bright): a=%.0f",aqAt),"p");
    lg->AddEntry(gb,Form("brightest slice: a=%.0f",abAt),"p"); lg->Draw();
    TLatex tt; tt.SetNDC(); tt.SetTextFont(62); tt.SetTextSize(0.034);
    tt.DrawLatex(0.12,0.94,Form("Both selections share ONE depth floor (~%.0f ps); they differ only in light/slew (a)",bestB));
    gSystem->mkdir(Form("figures/%d/narrative",radYear()),kTRUE); c->Print(radFigP("figures/narrative/floor_fit.png"));
    printf("  wrote figures/narrative/floor_fit.png\n");
}

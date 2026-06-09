// ============================================================================
// methodCompare.C — clean HEAD-TO-HEAD of the timing extraction: cfd05 (CFD at 5%
// of the CLIPPED measured peak, the published method) vs hg_lgcfd (CFD at 15% of the
// LG-predicted TRUE peak, on the recovered steep edge). SAME build, SAME brightest-1000
// selection, SAME (DW-UP)/2 estimator and a/sqrt(E)(+)b fit -- only the per-channel
// time source differs. Isolates the gain from recovering the saturated edge.
//   source setup.sh; root -l -b -q 'analyze/studies/methodCompare.C+("DSB1")'
// ============================================================================
#include "RadView.h"
#include "RadTiming.h"
#include "DataPaths.h"
#include "PlotUtils.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include "FigPaths.h"
#include <vector>
#include <cmath>
#include <cstdio>
using namespace rad;

static void ladder(BuildConfig& cfg,int src,std::vector<double>&E,std::vector<double>&S,std::vector<double>&Se,std::vector<double>&ze){
    const double Es[]={25,50,75,100,125,150};
    for(double e:Es){ TFile* fp=TFile::Open(radReduced("DSB1",e)); if(!fp||fp->IsZombie())continue;
        TTree* t=(TTree*)fp->Get("rad"); RadView v; v.attach(t,&cfg);
        TimingResult r=timingBrightestK(v,e,src,1000);
        if(r.sigma_ps>0){ E.push_back(e); S.push_back(r.sigma_ps); Se.push_back(r.sigma_ps/std::sqrt(2000.0)); ze.push_back(0); }
        fp->Close(); }
}
static void fitAB(TGraphErrors* g,double&a,double&b,double&be,double&chi){
    TF1 f("f","sqrt([0]*[0]/x+[1]*[1])",20,160); f.SetParameters(200,18); g->Fit(&f,"QN");
    a=std::fabs(f.GetParameter(0)); b=std::fabs(f.GetParameter(1));
    chi=f.GetChisquare()/std::max(1,f.GetNDF()); be=f.GetParError(1); if(chi>1) be*=std::sqrt(chi);
}

void methodCompare(const char* build="DSB1"){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    BuildConfig cfg=BuildConfig::Load(radConfig(build).Data());
    std::vector<double> Ec,Sc,Sec,zc, El,Sl,Sel,zl;
    ladder(cfg,RadView::kCFD05,Ec,Sc,Sec,zc);
    ladder(cfg,RadView::kLGCFD,El,Sl,Sel,zl);
    TGraphErrors* gc=new TGraphErrors(Ec.size(),&Ec[0],&Sc[0],&zc[0],&Sec[0]);
    TGraphErrors* gl=new TGraphErrors(El.size(),&El[0],&Sl[0],&zl[0],&Sel[0]);
    double ac,bc,bec,chic, al,bl,bel,chil; fitAB(gc,ac,bc,bec,chic); fitAB(gl,al,bl,bel,chil);
    printf("\n=== %s: cfd05 (clipped foot) vs hg_lgcfd (recovered edge), brightest-1000, same selection ===\n",build);
    printf("   E      cfd05    lgcfd    gain\n");
    for(size_t i=0;i<Ec.size()&&i<El.size();++i) printf("  %3.0f   %6.1f   %6.1f   %+5.1f\n",Ec[i],Sc[i],Sl[i],Sc[i]-Sl[i]);
    printf("  cfd05 : a=%.0f  b=%.1f+-%.1f  chi2/ndf=%.1f\n",ac,bc,bec,chic);
    printf("  lgcfd : a=%.0f  b=%.1f+-%.1f  chi2/ndf=%.1f\n",al,bl,bel,chil);
    printf("  -> @150 GeV gain = %.1f ps ; stochastic a: %.0f -> %.0f (%.0f%%)\n",
           Sc.back()-Sl.back(),ac,al,100.0*(ac-al)/ac);

    TCanvas* c=new TCanvas("mc","",900,650);
    c->SetLeftMargin(0.12); c->SetRightMargin(0.04); c->SetTopMargin(0.06); c->SetBottomMargin(0.13); c->SetGridy();
    TH1F* fr=c->DrawFrame(0,0,165,std::max(50.0,Sc.empty()?50.0:Sc[0]*1.15));
    fr->SetTitle(";beam energy E (GeV);brightest-1000  (DW#minusUP)/2  #sigma_{t} (ps)");
    fr->GetYaxis()->SetTitleSize(0.044); fr->GetYaxis()->SetTitleOffset(1.25); fr->GetXaxis()->SetTitleSize(0.044);
    gc->SetMarkerStyle(24); gc->SetMarkerColor(kRed+1); gc->SetLineColor(kRed+1); gc->SetMarkerSize(1.5);
    gl->SetMarkerStyle(20); gl->SetMarkerColor(kAzure+2); gl->SetLineColor(kAzure+2); gl->SetMarkerSize(1.5);
    TF1* fc=new TF1("fc","sqrt([0]*[0]/x+[1]*[1])",20,160); fc->SetParameters(ac,bc); fc->SetLineColor(kRed+1); fc->SetLineWidth(2); fc->Draw("SAME");
    TF1* fl=new TF1("fl","sqrt([0]*[0]/x+[1]*[1])",20,160); fl->SetParameters(al,bl); fl->SetLineColor(kAzure+2); fl->SetLineWidth(2); fl->Draw("SAME");
    gc->Draw("P SAME"); gl->Draw("P SAME");
    TLegend* lg=new TLegend(0.40,0.68,0.95,0.90); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.034);
    lg->AddEntry(gc,Form("cfd05 (clipped foot):  a=%.0f, b=%.0f#pm%.0f ps",ac,bc,bec),"lp");
    lg->AddEntry(gl,Form("hg_lgcfd (recovered edge):  a=%.0f, b=%.0f#pm%.0f ps",al,bl,bel),"lp");
    lg->Draw();
    gSystem->mkdir(Form("figures/%d/narrative",radYear()),kTRUE); c->Print(radFigP(Form("figures/narrative/method_compare_%s.png",build)));
    printf("  wrote figures/narrative/method_compare_%s.png\n",build);
}

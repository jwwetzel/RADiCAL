// ============================================================================
// lightYieldThesis.C — Paper 1 thesis figure (published-consistent energy form).
// For all four builds: brightest-1000 (DW-UP)/2 sigma_t(E), each fit to the
// PUBLISHED form sigma_t = a/sqrt(E) (+) b. The STOCHASTIC term a tracks light
// yield (DSB1 high-light smallest, LuAG low-light largest); the FLOOR b is
// ~material-independent (~20 ps for the LYSO builds). => light yield, not species,
// sets the (stochastic) timing; the floor is shared detector physics. Adopts and
// extends NIM A 1068 (2024) 169737 (DSB1, b=17.5 ps); does NOT revise it.
//   source setup.sh; root -l -b -q 'analyze/studies/lightYieldThesis.C+'
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
#include <vector>
#include <cmath>
#include <cstdio>
#include <string>
using namespace rad;

static int robustSrc(const char* b){ std::string s=b; if(s=="LUAG"||s=="TENERGY") return RadView::kLED; return RadView::kLGCFD; }
static const char* srcName(int s){ return s==RadView::kLED?"led":(s==RadView::kLGCFD?"lgcfd":"cfd05"); }

void lightYieldThesis(){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    const char* builds[]={"DSB1","TENERGY","MIXED","LUAG"};
    const char* lab[]={"DSB1 (LYSO, high light)","TENERGY (3#timesLYSO+E)","MIXED (LYSO+LuAG)","LUAG (LuAG, low light)"};
    int col[4]={kAzure+2,kViolet+1,kOrange+8,kGreen+3};
    const double Es[]={25,50,75,100,125,150};

    TCanvas* c=new TCanvas("ly","",920,680);
    c->SetLeftMargin(0.12); c->SetRightMargin(0.04); c->SetTopMargin(0.09); c->SetBottomMargin(0.13); c->SetGridy();
    TH1F* fr=c->DrawFrame(0,0,165,70);
    fr->SetTitle(";beam energy E (GeV);brightest-1000  (DW#minusUP)/2  #sigma_{t} (ps)");
    fr->GetYaxis()->SetTitleSize(0.044); fr->GetYaxis()->SetTitleOffset(1.25); fr->GetXaxis()->SetTitleSize(0.044);
    TLegend* lg=new TLegend(0.45,0.62,0.95,0.89); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.032);
    printf("\n=== Paper 1 thesis: sigma_t = a/sqrt(E) (+) b per build ===\n");
    printf("build     src     a [ps.sqrt(GeV)]   floor b [ps]\n");
    for(int bi=0;bi<4;++bi){ const char* B=builds[bi]; int SRC=robustSrc(B);
        BuildConfig cfg=BuildConfig::Load(radConfig(B).Data()); if(!cfg.valid())continue;
        std::vector<double> E,S,Se,ze;
        for(double e:Es){ TFile* fp=TFile::Open(radReduced(B,e)); if(!fp||fp->IsZombie())continue;
            TTree* t=(TTree*)fp->Get("rad"); RadView v; v.attach(t,&cfg);
            TimingResult r=timingBrightestK(v,e,SRC,1000);
            if(r.sigma_ps>0){ E.push_back(e); S.push_back(r.sigma_ps); Se.push_back(r.sigma_ps/std::sqrt(2000.0)); ze.push_back(0); }
            fp->Close(); }
        if(E.size()<3) continue;
        TGraphErrors* g=new TGraphErrors(E.size(),&E[0],&S[0],&ze[0],&Se[0]);
        g->SetMarkerStyle(20+bi); g->SetMarkerColor(col[bi]); g->SetLineColor(col[bi]); g->SetMarkerSize(1.5);
        TF1* f=new TF1(Form("f%d",bi),"sqrt([0]*[0]/x+[1]*[1])",20,160); f->SetParameters(200,18);
        g->Fit(f,"QN"); double a=std::fabs(f->GetParameter(0)), b=std::fabs(f->GetParameter(1));
        double cn=f->GetChisquare()/std::max(1,f->GetNDF()); double be=f->GetParError(1); if(cn>1) be*=std::sqrt(cn);
        f->SetLineColor(col[bi]); f->SetLineWidth(2); f->SetLineStyle(1); f->Draw("SAME");
        g->Draw("P SAME");
        lg->AddEntry(g,Form("%s:  a=%.0f,  b=%.0f#pm%.0f ps",lab[bi],a,b,be),"lp");
        printf("%-8s  %-5s   %8.1f          %5.1f +- %.1f  (chi2/ndf=%.1f)\n",B,srcName(SRC),a,b,be,cn);
    }
    lg->Draw();
    gSystem->mkdir("figures/narrative",kTRUE); c->Print("figures/narrative/light_yield_thesis.png");
    printf("  wrote figures/narrative/light_yield_thesis.png\n");
}

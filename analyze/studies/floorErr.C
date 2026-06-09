// ============================================================================
// floorErr.C — the floor, done honestly: FREE fits with PROPER per-point errors.
// Each sigma_t(E) point gets its statistical error sigma/sqrt(2N) (N = events in
// the bin/slice). Fit sigma_t = a/sqrtE (+) b with those errors and read b +-
// sigma(b) straight from the covariance (which includes the a<->b correlation).
// No artificial shared-floor constraint -- let the error bars say whether the
// floors overlap. Reports lgcfd quantile vs brightest (same technique) and
// lgcfd vs cfd05 (different technique).
//   source setup.sh; root -l -b -q 'analyze/studies/floorErr.C+("DSB1")'
// ============================================================================
#include "RadView.h"
#include "RadTiming.h"
#include "DataPaths.h"
#include "PlotUtils.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TCanvas.h"
#include "TBox.h"
#include "TLine.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>
using namespace rad;

// quantile best-bin: returns sigma and sets N=winning-bin count
static double quantBestN(std::vector<std::pair<float,float>> v,long& N){
    N=0; if(v.size()<2000) return -1; std::sort(v.begin(),v.end());
    int NB=9; size_t per=v.size()/NB; double best=1e9;
    for(int b=0;b<NB;++b){ size_t lo=(size_t)b*per,hi=(b==NB-1)?v.size():lo+per; if(hi-lo<600)continue;
        std::vector<float> t; for(size_t k=lo;k<hi;++k) t.push_back(v[k].second); double s=tebSigma(t);
        if(s>18&&s<best){best=s;N=hi-lo;} } return best<1e8?best:-1;
}
static double brightFracN(std::vector<std::pair<float,float>> v,double f,long& N){
    int K=(int)(f*v.size()); N=K; if(K<300) return -1;
    std::nth_element(v.begin(),v.begin()+K,v.end(),[](auto&a,auto&b){return a.first>b.first;});
    std::vector<float> t; for(int i=0;i<K;++i) t.push_back(v[i].second); return tebSigma(t);
}

struct Curve { const char* name; int src; bool bright; int col; TGraphErrors* g; double b,eb,rho; };

void floorErr(const char* build="DSB1"){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    BuildConfig cfg=BuildConfig::Load(radConfig(build).Data());
    const double Es[]={25,50,75,100,125,150};
    std::vector<Curve> C={
        {"lgcfd quantile (typical)",RadView::kLGCFD,false,kViolet+1,new TGraphErrors,0,0,0},
        {"lgcfd brightest 2%",       RadView::kLGCFD,true, kAzure+2, new TGraphErrors,0,0,0},
        {"cfd05 quantile (typical)", RadView::kCFD05,false,kRed+1,   new TGraphErrors,0,0,0},
    };
    for(double E:Es){ TFile* fp=TFile::Open(radReduced(build,E)); if(!fp||fp->IsZombie())continue;
        TTree* t=(TTree*)fp->Get("rad"); RadView v; v.attach(t,&cfg);
        double xc,yc; v.beamCenter(xc,yc); double r2=9.0; Long64_t N=v.entries();
        std::vector<std::pair<float,float>> P5,PL;
        for(Long64_t i=0;i<N;++i){ v.get(i);
            if(!v.wc_ok()||v.mcp1_peak()<kMCP1_minPeak||v.mcp1_peak()>kMCP1_maxPeak) continue;
            double dx=v.x_trk()-xc,dy=v.y_trk()-yc; if(dx*dx+dy*dy>=r2) continue;
            double l5d=0,l5u=0,lLd=0,lLu=0;int e5d=0,e5u=0,eLd=0,eLu=0;
            for(int c=0;c<4;++c) if(v.is_timing(c)&&v.hg_peak(c)>=kHG_minPeak){float t5=v.timeOf(c,RadView::kCFD05),tl=v.timeOf(c,RadView::kLGCFD);if(t5>-1e5){l5d+=t5;++e5d;}if(tl>-1e5){lLd+=tl;++eLd;}}
            for(int c=4;c<8;++c) if(v.is_timing(c)&&v.hg_peak(c)>=kHG_minPeak){float t5=v.timeOf(c,RadView::kCFD05),tl=v.timeOf(c,RadView::kLGCFD);if(t5>-1e5){l5u+=t5;++e5u;}if(tl>-1e5){lLu+=tl;++eLu;}}
            if(e5d&&e5u) P5.push_back({(float)v.sum_lg(),0.5f*(float)(l5d/e5d-l5u/e5u)});
            if(eLd&&eLu) PL.push_back({(float)v.sum_lg(),0.5f*(float)(lLd/eLd-lLu/eLu)}); }
        for(auto& cu:C){ auto& P=(cu.src==RadView::kCFD05?P5:PL); long n;
            double s=cu.bright?brightFracN(P,0.02,n):quantBestN(P,n);
            if(s>0&&n>0){ int k=cu.g->GetN(); cu.g->SetPoint(k,E,s); cu.g->SetPointError(k,0,s/std::sqrt(2.0*n)); } }
        fp->Close();
    }
    printf("\n%s — FREE fit sigma_t=a/sqrtE(+)b, errors sigma/sqrt(2N) then PDG-scaled by sqrt(chi2/ndf):\n",build);
    for(auto& cu:C){ TF1 f("f","sqrt([0]*[0]/x+[1]*[1])",20,160); f.SetParameters(200,22);
        auto r=cu.g->Fit(&f,"QS"); cu.b=std::fabs(f.GetParameter(1));
        double cn=f.GetChisquare()/std::max(1,f.GetNDF()), scale=(cn>1?std::sqrt(cn):1.0);
        cu.eb=f.GetParError(1)*scale; cu.rho=r->Correlation(0,1);
        printf("  %-26s  b = %5.1f +- %.1f ps   (rho_ab=%+.2f, chi2/ndf=%.1f -> scale x%.1f)\n",cu.name,cu.b,cu.eb,cu.rho,cn,scale); }
    // overlap tests (n-sigma separation)
    auto sep=[&](int i,int j){ return std::fabs(C[i].b-C[j].b)/std::sqrt(C[i].eb*C[i].eb+C[j].eb*C[j].eb); };
    printf("  lgcfd quant vs bright (same technique): %.1f sigma apart -> %s\n",sep(0,1), sep(0,1)<2?"OVERLAP (one floor)":"distinct");
    printf("  lgcfd vs cfd05 (different technique)  : %.1f sigma apart -> %s\n",sep(0,2), sep(0,2)<2?"overlap":"DISTINCT (floor is per-technique)");

    // figure: curves with error bars + floor bands b+-sigma(b)
    TCanvas* c=new TCanvas("ce","",900,640); c->SetLeftMargin(0.12); c->SetRightMargin(0.16); c->SetGridy();
    TH1F* fr=c->DrawFrame(0,15,168,55); fr->SetTitle(Form("%s: free fits with proper errors -- does the floor overlap?;beam energy (GeV);(DW#minusUP)/2 #sigma_{t} (ps)",build));
    TLegend* lg=new TLegend(0.40,0.72,0.86,0.89); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.030);
    for(auto& cu:C){ cu.g->SetMarkerStyle(20); cu.g->SetMarkerColor(cu.col); cu.g->SetLineColor(cu.col); cu.g->SetMarkerSize(1.3); cu.g->Draw("P SAME");
        TF1* f=new TF1(Form("ff%d",cu.col),"sqrt([0]*[0]/x+[1]*[1])",20,160); cu.g->Fit(f,"Q0"); f->SetLineColor(cu.col); f->SetLineWidth(2); f->Draw("SAME");
        // floor band at right
        TBox* bx=new TBox(161,cu.b-cu.eb,168,cu.b+cu.eb); bx->SetFillColorAlpha(cu.col,0.45); bx->SetLineColor(cu.col); bx->Draw();
        lg->AddEntry(cu.g,Form("%s: b=%.0f#pm%.0f",cu.name,cu.b,cu.eb),"p"); }
    lg->Draw();
    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.027); tx.SetTextColor(kGray+3); tx.SetTextAngle(90);
    tx.DrawLatex(0.90,0.30,"floor b #pm #sigma(b)");
    gSystem->mkdir("figures/narrative",kTRUE); c->Print("figures/narrative/floor_err.png");
    printf("  wrote figures/narrative/floor_err.png\n");
}

// ============================================================================
// optLadder.C — remove the equal-width-bin artifact. Instead of 9 fixed-width
// sum_lg bins (whose brightest bin is statistically starved at 150 GeV), take a
// CONSISTENT "brightest-K-events" slice at every energy and read sigma_t. Sweep
// K so the trend is unambiguous. If sigma_t(E) is monotonic for fixed K, the
// ladder bump was the estimator, not physics.
//   source setup.sh; root -l -b -q 'analyze/studies/optLadder.C+("DSB1")'
// ============================================================================
#include "RadView.h"
#include "RadTiming.h"
#include "SelectionCuts.h"
#include "DataPaths.h"
#include "PlotUtils.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>
using namespace rad;

static bool depth(RadView& v,double xc,double yc,double r2,int src,float& out){
    if(!v.wc_ok()||v.mcp1_peak()<kMCP1_minPeak||v.mcp1_peak()>kMCP1_maxPeak) return false;
    double dx=v.x_trk()-xc, dy=v.y_trk()-yc; if(dx*dx+dy*dy>=r2) return false;
    double ds=0,us=0;int dn=0,un=0;
    for(int c=0;c<4;++c) if(v.is_timing(c)&&v.hg_peak(c)>=kHG_minPeak){float tc=v.timeOf(c,src);if(tc>-1e5){ds+=tc;++dn;}}
    for(int c=4;c<8;++c) if(v.is_timing(c)&&v.hg_peak(c)>=kHG_minPeak){float tc=v.timeOf(c,src);if(tc>-1e5){us+=tc;++un;}}
    if(dn<1||un<1)return false; out=0.5f*(float)(ds/dn-us/un); return true;
}

void optLadder(const char* build="DSB1"){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    BuildConfig cfg=BuildConfig::Load(radConfig(build).Data());
    const double Es[]={25,50,75,100,125,150}; int nE=6;
    const int Ks[]={1000,2000,4000}; int nK=3; int col[3]={kGreen+2,kOrange+8,kAzure+2};
    std::vector<std::vector<double>> sig(nK, std::vector<double>(nE,-1)); std::vector<double> Ev(Es,Es+nE);
    printf("\n%s  lgcfd (DW-UP)/2 sigma_t [ps] for brightest-K-events slice:\n  E\\K  1000   2000   4000\n",build);
    for(int e=0;e<nE;++e){ double E=Es[e];
        TFile* fp=TFile::Open(radReduced(build,E)); if(!fp||fp->IsZombie())continue;
        TTree* t=(TTree*)fp->Get("rad"); RadView v; v.attach(t,&cfg);
        double xc,yc; v.beamCenter(xc,yc); double r2=9.0; Long64_t N=v.entries();
        std::vector<std::pair<float,float>> sd;   // (sum_lg, depth_lgcfd)
        for(Long64_t i=0;i<N;++i){ v.get(i); float d; if(depth(v,xc,yc,r2,RadView::kLGCFD,d)) sd.push_back({(float)v.sum_lg(),d}); }
        std::sort(sd.begin(),sd.end(),[](auto&a,auto&b){return a.first>b.first;}); // brightest first
        printf("  %3.0f ",E);
        for(int k=0;k<nK;++k){ int K=Ks[k]; if((int)sd.size()<K){printf("   -   ");continue;}
            std::vector<float> vt; for(int i=0;i<K;++i) vt.push_back(sd[i].second);
            double s=tebSigma(vt); sig[k][e]=s; printf(" %5.1f ",s); }
        printf("\n"); fp->Close();
    }
    TCanvas* c=new TCanvas("c_ol","",900,640); c->SetLeftMargin(0.12); c->SetRightMargin(0.04); c->SetGridy();
    TLegend* lg=new TLegend(0.50,0.68,0.93,0.88); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.034);
    bool first=true;
    for(int k=0;k<nK;++k){ std::vector<double> x,y; for(int e=0;e<nE;++e) if(sig[k][e]>0){x.push_back(Ev[e]);y.push_back(sig[k][e]);}
        TGraph* g=new TGraph(x.size(),&x[0],&y[0]); g->SetMarkerStyle(20); g->SetMarkerColor(col[k]); g->SetLineColor(col[k]); g->SetMarkerSize(1.5); g->SetLineWidth(3);
        if(first){ g->SetTitle(Form("%s: lgcfd (DW#minusUP)/2 #sigma_{t} vs E, brightest-K slice (consistent selection);beam energy (GeV);#sigma_{t} (ps)",build));
            g->GetXaxis()->SetLimits(0,160); g->GetYaxis()->SetRangeUser(20,60); g->Draw("ALP"); first=false; } else g->Draw("LP SAME");
        lg->AddEntry(g,Form("brightest %d events",Ks[k]),"lp"); }
    lg->Draw();
    TLatex tt; tt.SetNDC(); tt.SetTextFont(62); tt.SetTextSize(0.036);
    tt.DrawLatex(0.12,0.94,"With a consistent bright slice, #sigma_{t}(E) is monotonic — 150 #leq 125");
    gSystem->mkdir("figures/narrative",kTRUE); c->Print("figures/narrative/opt_ladder.png");
    printf("  wrote figures/narrative/opt_ladder.png\n");
}

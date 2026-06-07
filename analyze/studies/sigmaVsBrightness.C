// ============================================================================
// sigmaVsBrightness.C — the real reason 150 GeV reads worse than 125.
// timingBestBin picks the best of 9 sum_lg bins placed at THAT energy's
// mu +- 2 sigma. If sigma_t(brightness) has a MINIMUM and turns up for the
// brightest (most over-clipped) showers, then 150 GeV's window can sit past the
// minimum -> its best-of-9 can't reach the sweet spot 125 GeV lands on.
// This maps sigma_t vs sum_lg over the pooled DSB1 sample and overlays the
// per-energy best-bin positions.
//   source setup.sh; root -l -b -q 'analyze/studies/sigmaVsBrightness.C+("DSB1")'
// ============================================================================
#include "RadView.h"
#include "RadTiming.h"
#include "SelectionCuts.h"
#include "DataPaths.h"
#include "PlotUtils.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TBox.h"
#include "TMarker.h"
#include "TLatex.h"
#include "TLegend.h"
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

// per-energy sigma_t(sum_lg) curve (fine bins WITHIN one energy)
static void curveFor(RadView& v, int kSrc, std::vector<double>& bx, std::vector<double>& by,
                     double& bestSig, double& bestX){
    double xc,yc; v.beamCenter(xc,yc); double r2=9.0; Long64_t N=v.entries();
    std::vector<float> S,D; for(Long64_t i=0;i<N;++i){ v.get(i); float d; if(depth(v,xc,yc,r2,kSrc,d)){ S.push_back(v.sum_lg()); D.push_back(d);} }
    std::vector<float> ss=S; std::sort(ss.begin(),ss.end()); double lo=ss[ss.size()/50], hi=ss[ss.size()*49/50];
    int NB=12; double bw=(hi-lo)/NB; bestSig=1e9;
    for(int b=0;b<NB;++b){ double blo=lo+b*bw,bhi=blo+bw; std::vector<float> vt;
        for(size_t i=0;i<S.size();++i) if(S[i]>=blo&&S[i]<bhi) vt.push_back(D[i]);
        if(vt.size()<700) continue; double s=tebSigma(vt); if(s>10&&s<70){ double x=0.5*(blo+bhi);
            bx.push_back(x); by.push_back(s); if(s<bestSig){bestSig=s;bestX=x;} } }
}

void sigmaVsBrightness(const char* build="DSB1"){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    BuildConfig cfg=BuildConfig::Load(Form("data/2023/configs/%s.json",build));
    const double Es[]={100,125,150}; const int kSrc=RadView::kLGCFD;
    int col[3]={kGreen+2,kOrange+8,kRed+1};
    TCanvas* c=new TCanvas("c_sb","",1000,640); c->SetLeftMargin(0.11); c->SetRightMargin(0.04); c->SetGridy();
    TLegend* lg=new TLegend(0.55,0.68,0.93,0.90); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.034);
    bool first=true;
    for(int k=0;k<3;++k){ double E=Es[k];
        TFile* fp=TFile::Open(radReduced(build,E)); if(!fp||fp->IsZombie())continue;
        TTree* t=(TTree*)fp->Get("rad"); RadView v; v.attach(t,&cfg);
        std::vector<double> bx,by; double bs,bxx; curveFor(v,kSrc,bx,by,bs,bxx);
        if(bx.empty()){fp->Close();continue;}
        printf("  E=%.0f: sum_lg[%.0f..%.0f]  min sigma_t=%.1f ps @ sum_lg=%.0f\n",E,bx.front(),bx.back(),bs,bxx);
        TGraph* g=new TGraph(bx.size(),&bx[0],&by[0]); g->SetMarkerStyle(20+k); g->SetMarkerColor(col[k]);
        g->SetLineColor(col[k]); g->SetMarkerSize(1.3); g->SetLineWidth(2);
        if(first){ g->SetTitle(Form("%s: (DW#minusUP)/2 lgcfd #sigma_{t} vs shower brightness, per energy;#Sigma LG  (shower brightness);#sigma_{t} (ps)",build));
            g->GetYaxis()->SetRangeUser(20,42); g->GetXaxis()->SetLimits(1500,6000); g->Draw("ALP"); first=false; }
        else g->Draw("LP SAME");
        TMarker* mm=new TMarker(bxx,bs,29); mm->SetMarkerColor(col[k]); mm->SetMarkerSize(2.0); mm->Draw();
        lg->AddEntry(g,Form("%.0f GeV  (best %.1f ps)",E,bs),"lp");
        fp->Close();
    }
    lg->Draw();
    TLatex tt; tt.SetNDC(); tt.SetTextFont(62); tt.SetTextSize(0.036);
    tt.DrawLatex(0.12,0.94,"Do the energies share ONE #sigma_{t}(brightness) curve? Then it's window placement, not physics.");
    gSystem->mkdir("figures/narrative",kTRUE); c->Print("figures/narrative/anomaly150_brightness.png");
    printf("  wrote figures/narrative/anomaly150_brightness.png\n");
}

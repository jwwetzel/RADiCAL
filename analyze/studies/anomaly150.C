// ============================================================================
// anomaly150.C — why is DSB1 150 GeV (75 chained runs) worse than 125 (35 runs)?
// Splits the MERGED reduced file by per-event run id and asks: is the pooled
// (DW-UP)/2 sigma_t inflated by run-to-run beam drift? Compares
//   (a) pooled, ONE global LG-weighted centroid  (what timingBestBin does)
//   (b) pooled, each run RE-CENTERED on its own centroid
//   (c) per-run sigma_t (each run alone, own centroid)
// Plots sigma_t(run) and the mean (DW-UP)/2 offset + beam center vs run.
//   source setup.sh; root -l -b -q 'analyze/studies/anomaly150.C+("DSB1",150)'
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
#include "TLine.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include <map>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>
using namespace rad;

// (DW-UP)/2 lgcfd for one event given a beam center; returns false if not usable
static bool depthVal(RadView& v, double xc, double yc, double r2, float& out){
    if(!v.wc_ok()||v.mcp1_peak()<kMCP1_minPeak||v.mcp1_peak()>kMCP1_maxPeak) return false;
    double dx=v.x_trk()-xc, dy=v.y_trk()-yc; if(dx*dx+dy*dy>=r2) return false;
    double ds=0,us=0; int dn=0,un=0;
    for(int c=0;c<4;++c) if(v.is_timing(c)&&v.hg_peak(c)>=kHG_minPeak){ float tc=v.timeOf(c,RadView::kLGCFD); if(tc>-1e5){ds+=tc;++dn;} }
    for(int c=4;c<8;++c) if(v.is_timing(c)&&v.hg_peak(c)>=kHG_minPeak){ float tc=v.timeOf(c,RadView::kLGCFD); if(tc>-1e5){us+=tc;++un;} }
    if(dn<1||un<1) return false; out=0.5f*(float)(ds/dn-us/un); return true;
}

void anomaly150(const char* build="DSB1", double E=150){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    BuildConfig cfg = BuildConfig::Load(radConfig(build).Data());
    TFile* fp=TFile::Open(radReduced(build,E)); if(!fp||fp->IsZombie()){printf("no file\n");return;}
    TTree* t=(TTree*)fp->Get("rad"); RadView v; v.attach(t,&cfg);
    double xcg,ycg; v.beamCenter(xcg,ycg); double r2=3.0*3.0;
    Long64_t N=v.entries();

    // pass 1: per-run LG-weighted centroid + event count
    std::map<int,double> sx,sy,sw; std::map<int,long> nrun;
    for(Long64_t i=0;i<N;++i){ v.get(i); int r=v.ev.run; double w=v.sum_lg();
        if(w<=0||!v.wc_ok()) continue; sx[r]+=w*v.x_trk(); sy[r]+=w*v.y_trk(); sw[r]+=w; nrun[r]++; }
    std::map<int,double> cx,cy; for(auto&kv:sw) if(kv.second>0){ cx[kv.first]=sx[kv.first]/kv.second; cy[kv.first]=sy[kv.first]/kv.second; }
    printf("[anomaly150] %s %.0f GeV: %lld events, %zu runs, global center (%.2f,%.2f)\n",
           build,E,(long long)N,nrun.size(),xcg,ycg);

    // pass 2: collect (DW-UP)/2 three ways
    std::vector<float> poolGlobal, poolRecentered;
    std::map<int,std::vector<float>> perRun;
    for(Long64_t i=0;i<N;++i){ v.get(i); int r=v.ev.run; float d;
        if(depthVal(v,xcg,ycg,r2,d)) poolGlobal.push_back(d);                 // (a) one global center
        if(cx.count(r) && depthVal(v,cx[r],cy[r],r2,d)){ poolRecentered.push_back(d); perRun[r].push_back(d); } // (b)/(c) own center
    }
    double sigA=tebSigma(poolGlobal), sigB=tebSigma(poolRecentered);
    printf("  (a) pooled, ONE global center : sigma_t = %.1f ps   (N=%zu)\n", sigA, poolGlobal.size());
    printf("  (b) pooled, per-run recentered: sigma_t = %.1f ps   (N=%zu)\n", sigB, poolRecentered.size());

    // per-run sigma + offset (runs with >=1500 fiducial events)
    std::vector<double> rN, rS, rOff, rCx, rCy; double sumS=0; int nS=0;
    for(auto&kv:perRun){ if(kv.second.size()<1500) continue; double s=tebSigma(kv.second);
        double m=0; for(float x:kv.second) m+=x; m/=kv.second.size();
        if(s>5&&s<120){ rN.push_back(kv.first); rS.push_back(s); rOff.push_back(m); rCx.push_back(cx[kv.first]); rCy.push_back(cy[kv.first]); sumS+=s; ++nS; } }
    double meanPerRun=nS?sumS/nS:0;
    printf("  (c) per-run sigma_t (N>=1500): mean over %d runs = %.1f ps   [min %.1f]\n",
           nS, meanPerRun, rS.empty()?0:*std::min_element(rS.begin(),rS.end()));
    printf("  >> if (a) >> (b)~(c), the 150 GeV inflation is run-to-run beam drift, not physics.\n");

    // ---- plot: sigma_t(run) [top] + offset & center drift [bottom] ----------
    TCanvas* c=new TCanvas("c_an","",1000,780); c->Divide(1,2,0.004,0.03);
    c->cd(1); gPad->SetGridy(); gPad->SetLeftMargin(0.10); gPad->SetRightMargin(0.04);
    TGraph* gS=new TGraph(rN.size(),&rN[0],&rS[0]); gS->SetMarkerStyle(20); gS->SetMarkerColor(kAzure+2); gS->SetMarkerSize(1.1);
    gS->SetTitle(Form("%s %.0f GeV: per-run (DW#minusUP)/2 lgcfd #sigma_{t};run number;per-run #sigma_{t} (ps)",build,E));
    gS->GetYaxis()->SetRangeUser(0,60); gS->Draw("AP");
    TLine* lA=new TLine(rN.front(),sigA,rN.back(),sigA); lA->SetLineColor(kRed+1); lA->SetLineWidth(2); lA->Draw();
    TLine* lB=new TLine(rN.front(),sigB,rN.back(),sigB); lB->SetLineColor(kGreen+2); lB->SetLineStyle(2); lB->SetLineWidth(2); lB->Draw();
    TLatex tx; tx.SetTextSize(0.045); tx.SetTextColor(kRed+1); tx.DrawLatex(rN.front()+1,sigA+1.5,Form("pooled, 1 center = %.1f ps",sigA));
    tx.SetTextColor(kGreen+3); tx.DrawLatex(rN.front()+1,sigB-3.5,Form("pooled, recentered = %.1f ps",sigB));
    c->cd(2); gPad->SetGridy(); gPad->SetLeftMargin(0.10); gPad->SetRightMargin(0.04);
    TGraph* gO=new TGraph(rN.size(),&rN[0],&rOff[0]); gO->SetMarkerStyle(21); gO->SetMarkerColor(kOrange+7); gO->SetMarkerSize(1.0);
    gO->SetTitle("run-to-run drift;run number;mean (DW#minusUP)/2 (ps)  [#bullet]  &  beam x (mm) [#circ]");
    gO->Draw("AP");
    TGraph* gC=new TGraph(rN.size(),&rN[0],&rCx[0]); gC->SetMarkerStyle(24); gC->SetMarkerColor(kViolet+1); gC->SetMarkerSize(1.0); gC->Draw("P SAME");
    gSystem->mkdir("figures/narrative",kTRUE); c->Print("figures/narrative/anomaly150_perrun.png");
    printf("  wrote figures/narrative/anomaly150_perrun.png\n");
    fp->Close();
}

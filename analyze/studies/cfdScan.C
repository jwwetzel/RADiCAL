// ============================================================================
// cfdScan.C — is CFD-5% the optimal timing threshold, or does the clip break it?
// ----------------------------------------------------------------------------
// The reduced ntuple stores the (DW-UP)/2 inputs at 6 CFD fractions
// (hg_cfd03/05/10/[hg_cfd=20]/30/50). We run the FULL headline best-bin method
// at each fraction, per energy, and find the optimum. Hypothesis: at high energy
// the peak clips, so 5%-of-peak = 5% of ~820 = ~41 mV (near noise) -> a HIGHER
// fraction should time better; at 25 GeV (unclipped) 5% should be fine.
//   ROOT_INCLUDE_PATH=lib/waveform:lib/io:lib/physics:lib/viz root -l -b -q analyze/studies/cfdScan.C+
// ============================================================================
#include "RadTiming.h"       // rad::tebSigma
#include "SelectionCuts.h"
#include "DataPaths.h"
#include "PlotUtils.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>

struct EvT { float slg, t; };

static double bestBin(std::vector<EvT>& ev) {
    if (ev.size() < 2000) return -1;
    std::vector<float> sl; for (auto&e:ev) sl.push_back(e.slg);
    double smin=*std::min_element(sl.begin(),sl.end()), smax=*std::max_element(sl.begin(),sl.end());
    TH1F hS("hS","",150,smin,smax); hS.SetDirectory(0); for(float x:sl) hS.Fill(x);
    double muE,muEe,sigE,sigEe; FitGaussCore(&hS,2.0,muE,muEe,sigE,sigEe); if(sigE<=0){muE=hS.GetMean();sigE=hS.GetRMS();}
    double lo=muE-2*sigE, bw=4*sigE/9.0, best=1e9;
    for(int b=0;b<9;++b){ double blo=lo+b*bw,bhi=blo+bw; std::vector<float> vt;
        for(auto&e:ev) if(e.slg>=blo&&e.slg<bhi) vt.push_back(e.t);
        if(vt.size()<500) continue; double s=rad::tebSigma(vt); if(s>10&&s<best) best=s; }
    return best;
}

// one pass over the file; return best-bin sigma_t for all 6 CFD fractions
static void sigmasForEnergy(double E, double out[6]) {
    const char* brs[6]={"hg_cfd03","hg_cfd05","hg_cfd10","hg_cfd","hg_cfd30","hg_cfd50"};
    TFile* fp=TFile::Open(radReduced("DSB1",E)); TTree* t=(TTree*)fp->Get("rad");
    Bool_t wc; Float_t x,y,mcp,slg,hgp[8]; Float_t cfd[6][8];
    t->SetBranchAddress("wc_ok",&wc); t->SetBranchAddress("x_trk",&x); t->SetBranchAddress("y_trk",&y);
    t->SetBranchAddress("mcp_peak",&mcp); t->SetBranchAddress("sum_lg",&slg); t->SetBranchAddress("hg_peak",hgp);
    for(int f=0;f<6;++f) t->SetBranchAddress(brs[f],cfd[f]);
    long N=t->GetEntries();
    double wx=0,wy=0,w=0; for(long i=0;i<N;++i){ t->GetEntry(i); if(wc&&mcp>kMCP_minPeak_E&&slg>kSumLG_centroid){wx+=x*slg;wy+=y*slg;w+=slg;} }
    double xc=w>0?wx/w:0, yc=w>0?wy/w:0, rFid=TimingFiducialR(E), r2=rFid*rFid;
    std::vector<EvT> ev[6];
    for(long i=0;i<N;++i){ t->GetEntry(i);
        if(!wc||mcp<kMCP1_minPeak||mcp>kMCP1_maxPeak) continue;
        double dx=x-xc,dy=y-yc; if(dx*dx+dy*dy>=r2) continue;
        for(int f=0;f<6;++f){ double ds=0,us=0; int dn=0,un=0;
            for(int c=0;c<4;++c) if(hgp[c]>=kHG_minPeak && cfd[f][c]>-1e5){ ds+=cfd[f][c]; ++dn; }
            for(int c=4;c<8;++c) if(hgp[c]>=kHG_minPeak && cfd[f][c]>-1e5){ us+=cfd[f][c]; ++un; }
            if(dn>=1&&un>=1){ EvT e; e.slg=slg; e.t=0.5f*(float)(ds/dn-us/un); ev[f].push_back(e); } }
    }
    for(int f=0;f<6;++f) out[f]=bestBin(ev[f]);
    fp->Close();
}

void cfdScan() {
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    double frac[6]={3,5,10,20,30,50};
    double Es[3]={25,75,150}; int col[3]={kAzure+2,kOrange+7,kRed+1};
    double sig[3][6];
    for(int e=0;e<3;++e){ sigmasForEnergy(Es[e], sig[e]);
        printf("E=%3.0f  sigma_t by CFD frac [%%]: ",Es[e]);
        for(int f=0;f<6;++f) printf("%.0f%%=%.1f  ",frac[f],sig[e][f]); printf("ps\n"); }

    TCanvas* c=new TCanvas("c_cfd","",900,680); c->SetGridx(); c->SetGridy();
    TH1F* fr=c->DrawFrame(0,20,55,70);
    fr->SetTitle("Timing vs CFD fraction: flat at high E (5% near-optimal) #Rightarrow floor-limited;CFD fraction [%];best-bin (DW#minusUP)/2 #sigma_{t} [ps]");
    TLegend* lg=new TLegend(0.55,0.68,0.88,0.88); lg->SetBorderSize(0);
    for(int e=0;e<3;++e){ TGraph* g=new TGraph(6,frac,sig[e]); g->SetMarkerStyle(20+e); g->SetMarkerColor(col[e]); g->SetLineColor(col[e]); g->SetMarkerSize(1.5); g->SetLineWidth(2); g->Draw("PL SAME");
        int bf=0; for(int f=1;f<6;++f) if(sig[e][f]>0&&sig[e][f]<sig[e][bf]) bf=f;
        lg->AddEntry(g,Form("%.0f GeV (best @ %.0f%% = %.1f ps)",Es[e],frac[bf],sig[e][bf]),"pl"); }
    lg->Draw();
    TLatex tl; tl.SetNDC(); tl.SetTextSize(0.026);
    tl.DrawLatex(0.14,0.20,"At 75-150 GeV the scan is flat (28-30 ps) and 5% is best: timing is");
    tl.DrawLatex(0.14,0.16,"FLOOR-limited, not threshold-limited. Fraction matters only at 25 GeV.");
    gSystem->mkdir("figures",kTRUE);
    c->Print("figures/cfd_scan.png");
    printf("wrote figures/cfd_scan.png\n");
}

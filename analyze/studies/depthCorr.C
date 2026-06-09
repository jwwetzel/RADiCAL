// ============================================================================
// depthCorr.C — break the shower-depth timing floor with a depth-INDEPENDENT
//               (and depth-CORRECTED) shower-time estimator.   (canonical data)
// ----------------------------------------------------------------------------
// Per group-0 corner k (NW,NE,SE: both ends in DRS0-G0, referenced to MCP1):
//     t_down = T + (L - z)/v + d0 ,  t_up = T + z/v + d0      (d0 = G0 domino phase)
//   A = (t_down - t_up)/2 = (L - 2z)/(2v)        -> DEPTH  (shower-depth floor, ~28 ps)
//   S = (t_down + t_up)/2 = T + L/(2v) + d0      -> z CANCELS (depth-independent)
//   B = S - mcp_time = (T - beam) + const        -> d0 CANCELS (MCP1 same group)
//                                                   limited by MCP-copy noise +
//                                                   SiPM photostatistics (LIGHT!)
//   C = B - k*A                                  -> remove any residual depth walk
//
// A is pinned to shower-depth fluctuation (does NOT improve with light).  B is
// depth-free and should IMPROVE with energy (more light) -- the physics the user
// expects.  Compare best-bin sigma_t(E) for A, B, C.
//
// ---------------------------------------------------------------------------
// RESULT (DSB1, canonical reduced data): A is clean (58 ps @25 -> 31 ps @150).
// B is NOT recoverable from the reduced data:  sigma(B) ~ 2.5 ns and is
//   * NOT stop-cell-indexed  (mean(hg-mcp) is FLAT vs stopcell0, -135 ns), and
//   * corr(hg_cfd05, mcp_time) = -0.017  -- the SiPM and MCP crossing times are
//     UNCORRELATED.
// The ~26 ns common SiPM domino-phase fluctuation cancels in down-up (so the
// headline works) but the MCP copy does NOT share it -> the absolute SiPM<->MCP
// timebase is incoherent in this reduction.  This is the SAME missing ingredient
// as the 71 ps inter-group jitter: a per-event GLOBAL DRS4 timebase (recorded
// reference clock / reconstructed domino phase).  The depth-independent estimator
// is sound but requires a REDUCER-LEVEL absolute-timebase fix + re-reduction on
// Argon -- it cannot be built as a reduced-data macro.  This macro stands as the
// diagnostic of record.
//   ROOT_INCLUDE_PATH=lib/waveform:lib/io:lib/physics:lib/viz root -l -b -q analyze/studies/depthCorr.C+
// ============================================================================
#include "RadTiming.h"       // rad::tebSigma
#include "DataPaths.h"
#include "SelectionCuts.h"
#include "PlotUtils.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "FigPaths.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>

struct Dev { float slg, A, B, Bc; int sc0; char par; };  // Bc = stop-cell-unwrapped B

static double medOf(std::vector<float> v){ if(v.empty())return 0; size_t k=v.size()/2;
    std::nth_element(v.begin(),v.begin()+k,v.end()); return v[k]; }

// The MCP-referenced B = (DW+UP)/2 - MCP carries a per-event DRS4 stop-cell offset
// (crossing times are measured from the stop cell, which cancels in down-up but NOT
// vs the MCP at a different cell). Train an OOS stop-cell offset table d(stopcell0)
// on par==0 and subtract it -> Bc unwraps the absolute SiPM-MCP timing.
static void unwrapStopCell(std::vector<Dev>& ev){
    const int NB=128; std::vector<std::vector<float>> bin(NB); std::vector<float> all;
    for(auto&e:ev) if(e.par==0){ int b=e.sc0*NB/1024; if(b<0)b=0;if(b>=NB)b=NB-1; bin[b].push_back(e.B); all.push_back(e.B); }
    std::vector<float> off(NB,0); double gm=medOf(all);
    for(int b=0;b<NB;++b) off[b]= bin[b].size()>20 ? medOf(bin[b])-gm : 0.0;
    for(auto&e:ev){ int b=e.sc0*NB/1024; if(b<0)b=0;if(b>=NB)b=NB-1; e.Bc=e.B-off[b]; }
}

// 9-bin best-bin core sigma. which: 0=A(all ev), 1=Bc, 2=Bc-k*A (B/C use par==1, OOS).
static double bestBin(std::vector<Dev>& ev, int which, double k){
    std::vector<float> slg, t;
    for(auto&e:ev){ if(which>0 && e.par!=1) continue;
        float v = which==0? e.A : (which==1? e.Bc : e.Bc - (float)k*e.A);
        slg.push_back(e.slg); t.push_back(v); }
    if(slg.size()<2000) return -1;
    double smin=*std::min_element(slg.begin(),slg.end()), smax=*std::max_element(slg.begin(),slg.end());
    TH1F hS("hS","",150,smin,smax); hS.SetDirectory(0); for(float x:slg) hS.Fill(x);
    double muE,muEe,sigE,sigEe; FitGaussCore(&hS,2.0,muE,muEe,sigE,sigEe); if(sigE<=0){muE=hS.GetMean();sigE=hS.GetRMS();}
    double lo=muE-2*sigE, bw=4*sigE/9.0, best=1e9;
    for(int b=0;b<9;++b){ double blo=lo+b*bw,bhi=blo+bw; std::vector<float> vt;
        for(size_t i=0;i<slg.size();++i) if(slg[i]>=blo&&slg[i]<bhi) vt.push_back(t[i]);
        if(vt.size()<300) continue; double s=rad::tebSigma(vt); if(s>3&&s<best) best=s; }
    return best;
}

// robust slope k = cov(Bc,A)/var(A) on the par==1 core
static double depthSlope(std::vector<Dev>& ev){
    std::vector<float> bb; for(auto&e:ev) if(e.par==1) bb.push_back(e.Bc); double mB=medOf(bb);
    double ma=0,mb=0; long n=0; for(auto&e:ev) if(e.par==1&&std::fabs(e.Bc-mB)<0.40){ma+=e.A;mb+=e.Bc;++n;}
    if(n<50) return 0; ma/=n; mb/=n;
    double sxy=0,sxx=0;
    for(auto&e:ev) if(e.par==1&&std::fabs(e.Bc-mB)<0.40){ sxy+=(e.A-ma)*(e.Bc-mb); sxx+=(e.A-ma)*(e.A-ma); }
    return sxx>0? sxy/sxx : 0.0;
}

void depthCorr(){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    const double Es[6]={25,50,75,100,125,150};
    double sA[6], sB[6], sC[6];
    const int G0[6]={0,1,2,4,5,6};   // NW-D,NE-D,SE-D, NW-U,NE-U,SE-U  (all DRS0-G0, MCP1 ref)

    for(int e=0;e<6;++e){
        TFile* fp=TFile::Open(radReduced("DSB1",Es[e])); TTree* t=(TTree*)fp->Get("rad");
        Bool_t wc; Float_t x,y,mcp,mt,slg,cfd[8],hgp[8]; Int_t sc[4];
        t->SetBranchAddress("wc_ok",&wc); t->SetBranchAddress("x_trk",&x); t->SetBranchAddress("y_trk",&y);
        t->SetBranchAddress("mcp_peak",&mcp); t->SetBranchAddress("mcp_time",&mt); t->SetBranchAddress("sum_lg",&slg);
        t->SetBranchAddress("hg_cfd05",cfd); t->SetBranchAddress("hg_peak",hgp); t->SetBranchAddress("stopcell",sc);
        long N=t->GetEntries(); double wx=0,wy=0,w=0;
        for(long i=0;i<N;++i){ t->GetEntry(i); if(wc&&mcp>kMCP_minPeak_E&&slg>kSumLG_centroid){wx+=x*slg;wy+=y*slg;w+=slg;} }
        double xc=w>0?wx/w:0,yc=w>0?wy/w:0,rFid=TimingFiducialR(Es[e]),r2=rFid*rFid;

        std::vector<Dev> ev; ev.reserve(N); long ie=0;
        for(long i=0;i<N;++i){ t->GetEntry(i);
            if(!wc||mcp<kMCP1_minPeak||mcp>kMCP1_maxPeak||mt<-1e4) continue;
            double dx=x-xc,dy=y-yc; if(dx*dx+dy*dy>=r2) continue;
            double dsum=0,usum=0; int dn=0,un=0;
            for(int j=0;j<3;++j){ int c=G0[j];   if(hgp[c]>=kHG_minPeak&&cfd[c]>-1e4){dsum+=cfd[c];++dn;} }
            for(int j=3;j<6;++j){ int c=G0[j];   if(hgp[c]>=kHG_minPeak&&cfd[c]>-1e4){usum+=cfd[c];++un;} }
            if(dn<1||un<1) continue;
            double td=dsum/dn, tu=usum/un;
            Dev d; d.slg=slg; d.A=0.5f*(float)(td-tu); d.B=0.5f*(float)(td+tu)-mt; d.Bc=d.B;
            d.sc0=sc[0]; d.par=(char)(ie++&1); ev.push_back(d);
        }
        unwrapStopCell(ev);                       // OOS stop-cell unwrap of the absolute SiPM-MCP timing
        double k=depthSlope(ev);
        sA[e]=bestBin(ev,0,0); sB[e]=bestBin(ev,1,0); sC[e]=bestBin(ev,2,k);
        printf("E=%3.0f  (DW-UP)/2 depth A=%.1f ps | depth-indep B=(DW+UP)/2-MCP1 (stopcell-unwrapped)=%.1f ps | depth-corr C=%.1f ps  (k=%.3f, N=%zu)\n",
               Es[e], sA[e], sB[e], sC[e], k, ev.size());
        fp->Close();
    }

    TCanvas* c=new TCanvas("c_dc","",900,680);
    TH1F* fr=c->DrawFrame(0,0,160,110);
    fr->SetTitle("Depth-independent shower-time vs the depth-pinned headline;beam energy (GeV);best-bin #sigma_{t} (ps)");
    TGraph* gA=new TGraph(6,Es,sA); gA->SetMarkerStyle(20); gA->SetMarkerColor(kAzure+1); gA->SetLineColor(kAzure+1); gA->SetMarkerSize(1.6); gA->SetLineWidth(2); gA->Draw("PL SAME");
    TGraph* gB=new TGraph(6,Es,sB); gB->SetMarkerStyle(21); gB->SetMarkerColor(kRed+1);  gB->SetLineColor(kRed+1);  gB->SetMarkerSize(1.6); gB->SetLineWidth(2); gB->Draw("PL SAME");
    TGraph* gC=new TGraph(6,Es,sC); gC->SetMarkerStyle(33); gC->SetMarkerColor(kGreen+2);gC->SetLineColor(kGreen+2);gC->SetMarkerSize(2.0); gC->SetLineWidth(2); gC->Draw("PL SAME");
    TLegend* lg=new TLegend(0.40,0.70,0.88,0.88); lg->SetBorderSize(0);
    lg->AddEntry(gA,"(DW#minusUP)/2  depth-pinned (headline)","pl");
    lg->AddEntry(gB,"(DW+UP)/2 #minus MCP1  depth-independent","pl");
    lg->AddEntry(gC,"depth-corrected  B #minus k#upointA","pl"); lg->Draw();
    gSystem->mkdir("figures",kTRUE); c->Print(radFigP("figures/depth_corr.png"));
    printf("wrote figures/depth_corr.png\n");
}

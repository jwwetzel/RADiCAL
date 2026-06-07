// ============================================================================
// absShowerTime.C — absolute shower time from the REDUCED data, done correctly.
// ----------------------------------------------------------------------------
// KEY: the Reducer stores every hg time ALREADY MCP-referenced (Reducer.C:122,
//   ev.hg_led[i] = hg.ledTime - ref,  ref = mcp1 (G0) or mcp2 (SW-U)).
// So hg_led IS already (SiPM - MCP).  The absolute shower time is therefore
//   T_abs = mean_c( hg_led[c] )            <-- NO second mcp_time subtraction
// (my earlier absTime/depthCorr double-subtracted -> 1.7 ns artifact.)
// Compare best-bin core sigma_t(E):
//   A = (DW-UP)/2  (depth, MCP cancels)         -- the 27 ps headline
//   T = mean(hg_led)  absolute shower time        -- does it improve with light?
//   Tw = T with per-channel LED walk correction (hg_led vs 1/lg_peak)
//   ROOT_INCLUDE_PATH=lib/waveform:lib/io:lib/physics:lib/viz root -l -b -q analyze/studies/absShowerTime.C+
// ============================================================================
#include "RadTiming.h"
#include "DataPaths.h"
#include "SelectionCuts.h"
#include "PlotUtils.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>

struct Ev { float slg, A, T, Tw; };
static double medOf(std::vector<float> v){ if(v.empty())return 0; size_t k=v.size()/2;
    std::nth_element(v.begin(),v.begin()+k,v.end()); return v[k]; }

static double bb(std::vector<Ev>&ev, int which){ // 0=A,1=T,2=Tw
    std::vector<float> slg,t; for(auto&e:ev){ slg.push_back(e.slg);
        t.push_back(which==0?e.A:(which==1?e.T:e.Tw)); }
    if(slg.size()<2000) return -1;
    double smin=*std::min_element(slg.begin(),slg.end()), smax=*std::max_element(slg.begin(),slg.end());
    TH1F hS("hS","",150,smin,smax); hS.SetDirectory(0); for(float x:slg) hS.Fill(x);
    double muE,muEe,sigE,sigEe; FitGaussCore(&hS,2.0,muE,muEe,sigE,sigEe); if(sigE<=0){muE=hS.GetMean();sigE=hS.GetRMS();}
    double lo=muE-2*sigE, bw=4*sigE/9.0, best=1e9;
    for(int b=0;b<9;++b){ double blo=lo+b*bw,bhi=blo+bw; std::vector<float> vt;
        for(size_t i=0;i<slg.size();++i) if(slg[i]>=blo&&slg[i]<bhi) vt.push_back(t[i]);
        if(vt.size()<300) continue; double m=medOf(vt); std::vector<float> c;
        for(float x:vt) if(std::fabs(x-m)<1.0) c.push_back(x);
        if(c.size()<200) continue; double s=rad::tebSigma(c); if(s>3&&s<best) best=s; }
    return best;
}

void absShowerTime(){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    const double Es[6]={25,50,75,100,125,150};
    const int G0[7]={0,1,2,3,4,5,6};   // all 7 pure-G0 timing channels (mcp1-referenced)
    double sA[6],sT[6],sTw[6];
    for(int e=0;e<6;++e){
        TFile* fp=TFile::Open(radReduced("DSB1",Es[e])); TTree* t=(TTree*)fp->Get("rad");
        Bool_t wc; Float_t x,y,mcp,slg,led[8],hgp[8],lgp[8];
        t->SetBranchAddress("wc_ok",&wc); t->SetBranchAddress("x_trk",&x); t->SetBranchAddress("y_trk",&y);
        t->SetBranchAddress("mcp_peak",&mcp); t->SetBranchAddress("sum_lg",&slg);
        t->SetBranchAddress("hg_led",led); t->SetBranchAddress("hg_peak",hgp); t->SetBranchAddress("lg_peak",lgp);
        long N=t->GetEntries(); double wx=0,wy=0,w=0;
        for(long i=0;i<N;++i){ t->GetEntry(i); if(wc&&mcp>kMCP_minPeak_E&&slg>kSumLG_centroid){wx+=x*slg;wy+=y*slg;w+=slg;} }
        double xc=w>0?wx/w:0,yc=w>0?wy/w:0,rFid=TimingFiducialR(Es[e]),r2=rFid*rFid;

        struct Raw{float slg; float led[8],lg[8],hgp[8];};
        std::vector<Raw> raw;
        for(long i=0;i<N;++i){ t->GetEntry(i);
            if(!wc||mcp<kMCP1_minPeak||mcp>kMCP1_maxPeak) continue;
            double dx=x-xc,dy=y-yc; if(dx*dx+dy*dy>=r2) continue;
            Raw r; r.slg=slg; for(int c=0;c<8;++c){r.led[c]=led[c];r.lg[c]=lgp[c];r.hgp[c]=hgp[c];} raw.push_back(r);
        }
        // per-channel LED walk vs 1/lg_peak
        double aa[8]={0},bb_[8]={0};
        for(int c=0;c<8;++c){ double sx=0,sy=0,sxx=0,sxy=0; long n=0;
            for(auto&r:raw){ if(r.hgp[c]<kHG_minPeak||r.led[c]<-1e4||r.lg[c]<10) continue; double X=1.0/r.lg[c]; sx+=X;sy+=r.led[c];sxx+=X*X;sxy+=X*r.led[c];++n; }
            if(n<200){bb_[c]=0;aa[c]=0;continue;} double d=n*sxx-sx*sx; if(std::fabs(d)<1e-12){bb_[c]=0;aa[c]=0;continue;}
            bb_[c]=(n*sxy-sx*sy)/d; aa[c]=(sy-bb_[c]*sx)/n; }

        std::vector<Ev> ev;
        for(auto&r:raw){
            double ds=0,us=0; int dn=0,un=0; double Ts=0,Tws=0; int tn=0;
            for(int j=0;j<7;++j){ int c=G0[j]; if(r.hgp[c]<kHG_minPeak||r.led[c]<-1e4) continue;
                Ts+=r.led[c]; double wcorr=(r.lg[c]>10)?(aa[c]+bb_[c]/r.lg[c]):0; Tws+=r.led[c]-wcorr; ++tn;
                if(c<4){ds+=r.led[c];++dn;} else {us+=r.led[c];++un;} }
            if(tn<4||dn<1||un<1) continue;
            Ev e; e.slg=r.slg; e.A=0.5f*(float)(ds/dn-us/un); e.T=(float)(Ts/tn); e.Tw=(float)(Tws/tn); ev.push_back(e);
        }
        sA[e]=bb(ev,0); sT[e]=bb(ev,1); sTw[e]=bb(ev,2);
        printf("E=%3.0f  depth (DW-UP)/2=%.1f ps | ABSOLUTE shower time mean(hg_led)=%.1f ps | walk-corrected=%.1f ps  (N=%zu)\n",
               Es[e], sA[e], sT[e], sTw[e], ev.size());
        fp->Close();
    }
    TCanvas* c=new TCanvas("c_as","",900,680);
    TH1F* fr=c->DrawFrame(0,0,160,140);
    fr->SetTitle("Absolute shower time (hg_led is already MCP-referenced) vs depth;beam energy (GeV);best-bin #sigma_{t} (ps)");
    TGraph* gA=new TGraph(6,Es,sA); gA->SetMarkerStyle(20); gA->SetMarkerColor(kAzure+1); gA->SetLineColor(kAzure+1); gA->SetMarkerSize(1.6); gA->SetLineWidth(2); gA->Draw("PL SAME");
    TGraph* gT=new TGraph(6,Es,sT); gT->SetMarkerStyle(24); gT->SetMarkerColor(kGray+2); gT->SetLineColor(kGray+2); gT->SetMarkerSize(1.5); gT->SetLineWidth(2); gT->Draw("PL SAME");
    TGraph* gW=new TGraph(6,Es,sTw);gW->SetMarkerStyle(21); gW->SetMarkerColor(kRed+1); gW->SetLineColor(kRed+1); gW->SetMarkerSize(1.6); gW->SetLineWidth(2); gW->Draw("PL SAME");
    TLegend* lg=new TLegend(0.38,0.70,0.88,0.88); lg->SetBorderSize(0);
    lg->AddEntry(gA,"(DW#minusUP)/2 depth (headline)","pl");
    lg->AddEntry(gT,"absolute shower time mean(hg_led)","pl");
    lg->AddEntry(gW,"absolute, LED-walk corrected","pl"); lg->Draw();
    gSystem->mkdir("figures",kTRUE); c->Print("figures/abs_shower_time.png");
    printf("wrote figures/abs_shower_time.png\n");
}

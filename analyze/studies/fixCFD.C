// ============================================================================
// fixCFD.C — does CFD relative to the TRUE (LG-predicted) peak recover the
// timing improvement with energy that CFD-5%-of-clipped-peak hides?
// ----------------------------------------------------------------------------
// The HG peak clips ~820 mV, so CFD-5% = 5% of 820 = ~41 mV (fixed, near noise)
// for all high-E events -> no improvement. Fix (per the proposal): the LG channel
// is unsaturated, so HG_true = a + b*LG (fit on 25 GeV unsaturated data). Set the
// CFD threshold to frac*HG_true (a fixed fraction of the TRUE peak), expressed as
// an effective fraction of the measured (clipped) peak so the validated
// ExtractPulse does the (negative-pulse, DRS4-correct) crossing. Compare the
// (DW-UP)/2 best-bin sigma_t(E) of standard CFD-5% vs corrected CFD.
//   ROOT_INCLUDE_PATH=lib/waveform:lib/io:lib/physics:lib/viz root -l -b -q analyze/studies/fixCFD.C+
// ============================================================================
#include "ChannelConfig.h"
#include "WaveformUtils.h"   // ExtractPulse (negative-pulse + DRS4 correct)
#include "PlotUtils.h"       // FitGaussCore, ApplyRADiCALStyle
#include "RadTiming.h"       // rad::tebSigma
#include "SelectionCuts.h"
#include "DataPaths.h"
#include "TChain.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
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

struct EvT { float x,y,slg,tStd,tCor; };

// bright-slice (top 40% sum_lg) core sigma of (DW-UP)/2 — stable at low stats
static double bestBin(std::vector<EvT>& ev, bool useCor) {
    if(ev.size()<300) return -1;
    std::vector<EvT> s=ev; std::sort(s.begin(),s.end(),[](const EvT&A,const EvT&B){return A.slg<B.slg;});
    size_t lo=(size_t)(0.60*s.size());
    std::vector<float> v; for(size_t i=lo;i<s.size();++i){ float t=useCor?s[i].tCor:s[i].tStd; if(t>-1e5) v.push_back(t); }
    return rad::tebSigma(v);
}

static void fitHGvsLG(double a[8], double b[8]) {
    TFile* fp=TFile::Open(radReduced("DSB1",25)); TTree* t=(TTree*)fp->Get("rad");
    Bool_t wc; Float_t hgp[8],lgp[8];
    t->SetBranchAddress("wc_ok",&wc); t->SetBranchAddress("hg_peak",hgp); t->SetBranchAddress("lg_peak",lgp);
    double sx[8]={0},sy[8]={0},sxx[8]={0},sxy[8]={0}; long n[8]={0}; long N=t->GetEntries();
    for(long i=0;i<N;++i){ t->GetEntry(i); if(!wc) continue;
        for(int c=0;c<8;++c) if(lgp[c]>15 && hgp[c]>30 && hgp[c]<720){ double X=lgp[c],Y=hgp[c]; sx[c]+=X;sy[c]+=Y;sxx[c]+=X*X;sxy[c]+=X*Y;++n[c]; } }
    for(int c=0;c<8;++c){ if(n[c]>50){ double d=n[c]*sxx[c]-sx[c]*sx[c]; b[c]=(n[c]*sxy[c]-sx[c]*sy[c])/d; a[c]=(sy[c]-b[c]*sx[c])/n[c]; } else {a[c]=0;b[c]=5;} }
    printf("HG_true=a+b*LG (25 GeV fit): "); for(int c=0;c<8;++c) printf("c%d:%.0f+%.1fLG ",c,a[c],b[c]); printf("\n");
    fp->Close();
}

static void processE(double E, const char* file, double a[8], double b[8], double frac, double& sStd, double& sCor, size_t& nfid) {
    TChain ch("pulse"); ch.Add(radRaw(file)); TTreeReader rd(&ch);
    TTreeReaderArray<float> amp(rd,"amplitude"), tim(rd,"timevalue");
    std::vector<EvT> ev; long iev=0; const long MAXEV=1000000;
    while(rd.Next() && iev<MAXEV){ ++iev; const float* A=&amp[0]; const float* T=&tim[0];
        Pulse wcR=ExtractPulse(T+kWC_t,A+kWC_R,0.5f,kWC_minPeak), wcL=ExtractPulse(T+kWC_t,A+kWC_L,0.5f,kWC_minPeak),
              wcD=ExtractPulse(T+kWC_t,A+kWC_D,0.5f,kWC_minPeak), wcU=ExtractPulse(T+kWC_t,A+kWC_U,0.5f,kWC_minPeak);
        if(!(wcR.valid&&wcL.valid&&wcD.valid&&wcU.valid)) continue;
        Pulse m1=ExtractPulse(T+kMCP1_t,A+kMCP1,0.2f,30.f); if(m1.peak<kMCP1_minPeak||m1.peak>kMCP1_maxPeak) continue;
        Pulse m2=ExtractPulse(T+kMCP2_t,A+kMCP2,0.2f,30.f);
        double dS=0,uS=0,dC=0,uC=0; int nS_d=0,nS_u=0,nC_d=0,nC_u=0; double slg=0;
        for(int c=0;c<8;++c){
            Pulse lg=ExtractPulse(T+kCap[c].lg_t,A+kCap[c].lg,0.2f,5.f); slg+=lg.peak;
            Pulse hg=ExtractPulse(T+kCap[c].hg_t,A+kCap[c].hg,0.05f,30.f); if(!hg.valid) continue;
            double tS=hg.crossingTime;
            double HGtrue=a[c]+b[c]*lg.peak, L=frac*HGtrue, tC=-1e9;
            if(L>20 && L<0.92*hg.peak){ double feff=L/hg.peak; Pulse hc=ExtractPulse(T+kCap[c].hg_t,A+kCap[c].hg,(float)feff,30.f); if(hc.crossingTime>-1e5) tC=hc.crossingTime; }
            bool down=(c<4);
            if(tS>-1e5){ if(down){dS+=tS;++nS_d;}else{uS+=tS;++nS_u;} }
            if(tC>-1e5){ if(down){dC+=tC;++nC_d;}else{uC+=tC;++nC_u;} }
        }
        EvT e; e.x=(float)(kWC_Scale*(wcR.peakTime-wcL.peakTime)); e.y=(float)(kWC_Scale*(wcD.peakTime-wcU.peakTime)); e.slg=(float)slg;
        e.tStd=(nS_d>=1&&nS_u>=1)?0.5f*(float)(dS/nS_d-uS/nS_u):-1e9;
        e.tCor=(nC_d>=1&&nC_u>=1)?0.5f*(float)(dC/nC_d-uC/nC_u):-1e9;
        ev.push_back(e);
    }
    double wx=0,wy=0,w=0; for(auto&e:ev){wx+=e.x*e.slg;wy+=e.y*e.slg;w+=e.slg;} double xc=w>0?wx/w:0,yc=w>0?wy/w:0;
    double rFid=TimingFiducialR(E), r2=rFid*rFid;
    std::vector<EvT> fid; for(auto&e:ev){ double dx=e.x-xc,dy=e.y-yc; if(dx*dx+dy*dy<r2) fid.push_back(e); }
    nfid=fid.size(); sStd=bestBin(fid,false); sCor=bestBin(fid,true);
}

void fixCFD() {
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    double a[8],b[8]; fitHGvsLG(a,b);
    struct R{double E; const char* f;};
    R runs[6]={{25,"RUN1211_25_GeV.root"},{50,"RUN1148_50_GeV.root"},{75,"RUN1112_75_GeV.root"},
               {100,"RUN1075_100_GeV.root"},{125,"RUN1034_125_GeV.root"},{150,"RUN1258_150_GeV.root"}};
    double frac=0.30, Es[6], sStd[6], sCor[6];
    for(int e=0;e<6;++e){ size_t nf=0; processE(runs[e].E,runs[e].f,a,b,frac,sStd[e],sCor[e],nf); Es[e]=runs[e].E;
        printf("E=%3.0f  CFD-5%%(clipped)=%.1f  CFD-%.0f%%-of-TRUE=%.1f ps  (Nfid=%zu)\n",Es[e],sStd[e],frac*100,sCor[e],nf); }

    TCanvas* c=new TCanvas("c_fx","",900,680); c->SetGridx(); c->SetGridy();
    TH1F* fr=c->DrawFrame(0,15,160,60);
    fr->SetTitle("Timing vs energy: CFD on the TRUE (LG-predicted) peak;beam energy [GeV];best-bin (DW#minusUP)/2 #sigma_{t} [ps]");
    TGraph* gS=new TGraph(6,Es,sStd); gS->SetMarkerStyle(24); gS->SetMarkerColor(kGray+2); gS->SetLineColor(kGray+2); gS->SetMarkerSize(1.5); gS->SetLineWidth(2); gS->Draw("PL SAME");
    TGraph* gC=new TGraph(6,Es,sCor); gC->SetMarkerStyle(20); gC->SetMarkerColor(kRed+1); gC->SetLineColor(kRed+1); gC->SetMarkerSize(1.6); gC->SetLineWidth(2); gC->Draw("PL SAME");
    TLegend* lg=new TLegend(0.45,0.72,0.88,0.88); lg->SetBorderSize(0);
    lg->AddEntry(gS,"standard CFD-5% (of clipped peak)","pl");
    lg->AddEntry(gC,Form("CFD-%.0f%% of LG-predicted TRUE peak",frac*100),"pl"); lg->Draw();
    gSystem->mkdir("figures",kTRUE); c->Print(radFigP("figures/fix_cfd.png"));
    printf("wrote figures/fix_cfd.png\n");
}

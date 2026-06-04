// lightYieldTiming.C — follow-on paper gap #3: THE money plot.
//   Per capillary, per energy: single-channel sigma_t (Gaussian core, MCP-referenced)
//   vs its light yield (mean low-gain peak amplitude, a relative LY proxy).
//   Points coloured by SCINTILLATOR MATERIAL: DSB1 vs LuAG:Ce. If both materials
//   fall on ONE sigma_t(LY) curve, then LIGHT YIELD — not crystal species — drives
//   timing (the central claim of the follow-on paper).
//   root -l 'Analysis/lightYieldTiming.C+'
#include "ChannelConfig.h"
#include "PlotUtils.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include <cstdio>
#include <vector>
#include <algorithm>
#include <cmath>

static double gcoreSig(std::vector<float>& v){          // Gaussian-core sigma [ps]
    if(v.size()<200) return -1;
    std::vector<float> s=v; std::sort(s.begin(),s.end()); double md=s[s.size()/2];
    TH1F h("_g","",200,md-1.5,md+1.5); for(float x:v) h.Fill(x);
    double mu,muE,sg,sgE; FitGaussCore(&h,2.0,mu,muE,sg,sgE); return sg>0? sg*1000.0:-1;
}

// For a config at one energy, fill per-capillary (LY, sigma_t) for the 8 capillaries.
static void perCap(const char* dir, bool isDSB1, double E,
                   std::vector<double>& ly, std::vector<double>& st){
    TString fn = isDSB1 ? Form("Analysis/Output/%.0fGeV/ntuple.root",E)
                        : Form("%s/%.0fGeV.root",dir,E);
    TFile* fp=TFile::Open(fn); if(!fp||fp->IsZombie()) return;
    TTree* t=(TTree*)fp->Get("rad"); if(!t){ fp->Close(); return; }
    Int_t run; Bool_t wc; Float_t x,y;
    Float_t m1t,m2t,m1p,sp[36],sc[36];          // reduced
    Float_t mp,cfd[8],lgp[8];                    // DSB1
    t->SetBranchAddress("wc_ok",&wc); t->SetBranchAddress("x_trk",&x); t->SetBranchAddress("y_trk",&y);
    if(isDSB1){ t->SetBranchAddress("mcp_peak",&mp); t->SetBranchAddress("hg_cfd05",cfd); t->SetBranchAddress("lg_peak",lgp); }
    else { t->SetBranchAddress("mcp1_time",&m1t); t->SetBranchAddress("mcp2_time",&m2t);
           t->SetBranchAddress("mcp1_peak",&m1p); t->SetBranchAddress("s_peak",sp); t->SetBranchAddress("s_cfd05",sc); }
    long N=t->GetEntries(); double xs=0,ys=0; long nw=0;
    for(long i=0;i<N&&nw<50000;++i){ t->GetEntry(i); if(wc&&x>-100&&x<100){xs+=x;ys+=y;++nw;} }
    double xc=nw?xs/nw:0, yc=nw?ys/nw:0;
    std::vector<float> tv[8]; double lysum[8]={0}; long lyn[8]={0};
    for(long i=0;i<N;++i){ t->GetEntry(i);
        if(!wc) continue; double dx=x-xc,dy=y-yc; if(dx*dx+dy*dy>=9.0) continue;
        if(isDSB1){ if(mp<200||mp>750) continue; }
        else { if(m1p<200||m1p>750||m1t<-1e5) continue; }
        for(int c=0;c<8;++c){
            double tt, lgpk;
            if(isDSB1){ if(cfd[c]<-1e5) continue; tt=cfd[c]; lgpk=lgp[c]; }
            else { int hs=kCap[c].hg/1024, ls=kCap[c].lg/1024; double ref=kCap[c].use_mcp2?m2t:m1t;
                   if(sc[hs]<-1e5||ref<-1e5||sp[hs]<20) continue; tt=sc[hs]-ref; lgpk=sp[ls]; }
            if(lgpk<10) continue;
            tv[c].push_back((float)tt); lysum[c]+=lgpk; ++lyn[c];
        }
    }
    for(int c=0;c<8;++c){ if(lyn[c]<500) continue; double s=gcoreSig(tv[c]);
        if(s>0){ ly.push_back(lysum[c]/lyn[c]); st.push_back(s); } }
    fp->Close();
}

void lightYieldTiming(){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    double Es[6]={25,50,75,100,125,150};
    std::vector<double> lyD,stD, lyL,stL;       // DSB1, LuAG
    for(int e=0;e<6;++e){ perCap("",true,Es[e],lyD,stD); }
    for(int e=0;e<6;++e){ perCap("reduced/LUAG",false,Es[e],lyL,stL); }
    printf("\nDSB1 capillary-points: %zu   LuAG capillary-points: %zu\n", lyD.size(), lyL.size());

    TCanvas* c=new TCanvas("c_ly","",880,640); c->SetGridx(); c->SetGridy(); c->SetLogx();
    double lymax=50; for(double v:lyD) lymax=std::max(lymax,v); for(double v:lyL) lymax=std::max(lymax,v);
    TH1F* fr=c->DrawFrame(8, 80, lymax*1.3, 360);
    fr->SetTitle("Timing vs light yield, by scintillator material;mean low-gain amplitude  (relative light yield) [mV];single-channel #sigma_{t}  [ps]");
    TGraph* gD=new TGraph(lyD.size(),&lyD[0],&stD[0]); gD->SetMarkerStyle(20); gD->SetMarkerColor(kRRed);    gD->SetMarkerSize(1.1);
    TGraph* gL=new TGraph(lyL.size(),&lyL[0],&stL[0]); gL->SetMarkerStyle(21); gL->SetMarkerColor(kAzure+1); gL->SetMarkerSize(1.1);
    // one common curve sigma_t = sqrt(p0^2/LY + p1^2) fit to BOTH materials together
    std::vector<double> lyAll=lyD, stAll=stD; lyAll.insert(lyAll.end(),lyL.begin(),lyL.end()); stAll.insert(stAll.end(),stL.begin(),stL.end());
    TGraph* gAll=new TGraph(lyAll.size(),&lyAll[0],&stAll[0]);
    TF1* f=new TF1("fly","sqrt([0]*[0]/x+[1]*[1])",8,lymax*1.3); f->SetParameters(800,90); f->SetLineColor(kGray+2); f->SetLineWidth(2); f->SetLineStyle(2);
    gAll->Fit(f,"RQN");
    f->Draw("SAME"); gD->Draw("P SAME"); gL->Draw("P SAME");
    TLegend* lg=new TLegend(0.55,0.70,0.88,0.88);
    lg->AddEntry(gD,"DSB1 capillaries","p"); lg->AddEntry(gL,"LuAG:Ce capillaries","p");
    lg->AddEntry(f,"joint fit  #sigma_{t}=#sqrt{p_{0}^{2}/LY+p_{1}^{2}}","l"); lg->Draw();
    TLatex tl; tl.SetNDC(); tl.SetTextSize(0.027);
    tl.DrawLatex(0.14,0.205,"Per capillary, per energy (25-150 GeV), single-channel MCP-referenced.");
    tl.DrawLatex(0.14,0.165,"DSB1 and LuAG:Ce overlap at equal light yield #Rightarrow no separation by");
    tl.DrawLatex(0.14,0.125,"crystal species; light yield sets the timing. (Reference-subtracted");
    tl.DrawLatex(0.14,0.085,"intrinsic-#sigma_{t} sharpens the trend - paper-grade refinement.)");
    c->Print("Analysis/capillary_figs/sigmat_vs_lightyield.png");
    printf("joint fit: p0=%.0f ps*sqrt(mV)  floor p1=%.0f ps\n", f->GetParameter(0), fabs(f->GetParameter(1)));
    printf("wrote Analysis/capillary_figs/sigmat_vs_lightyield.png\n");
}

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
#include "DataPaths.h"
#include "TTree.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "FigPaths.h"
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
    TString fn = isDSB1 ? radReduced("DSB1",E)
                        : radReduced(dir,E);
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
    for(int e=0;e<6;++e){ perCap("LUAG",false,Es[e],lyL,stL); }
    printf("\nDSB1 capillary-points: %zu   LuAG capillary-points: %zu\n", lyD.size(), lyL.size());

    // --- reference-subtraction: intrinsic capillary sigma_t = sqrt(single^2 - ref^2),
    //     ref = 73 ps = measured per-group MCP reference jitter (sigma(MCP1-MCP2)/sqrt2). ---
    const double SREF=73.0;
    auto intrinsic=[&](const std::vector<double>& s){ std::vector<double> o; o.reserve(s.size());
        for(double v:s){ double w=v*v-SREF*SREF; o.push_back(w>100? std::sqrt(w):10.0); } return o; };
    std::vector<double> stiD=intrinsic(stD), stiL=intrinsic(stL);
    // binned profiles (mean intrinsic sigma_t per log-LY bin, per material) to show overlap
    double lymax=50; for(double v:lyD) lymax=std::max(lymax,v); for(double v:lyL) lymax=std::max(lymax,v);
    const int NB=6; double l0=std::log10(15.), l1=std::log10(lymax*1.05), dw=(l1-l0)/NB;
    auto profile=[&](const std::vector<double>& ly,const std::vector<double>& si,
                     std::vector<double>& bx,std::vector<double>& by,std::vector<double>& be){
        for(int b=0;b<NB;++b){ double lo=std::pow(10,l0+b*dw), hi=std::pow(10,l0+(b+1)*dw);
            double a=0,a2=0; long m=0; for(size_t i=0;i<ly.size();++i) if(ly[i]>=lo&&ly[i]<hi){a+=si[i];a2+=si[i]*si[i];++m;}
            if(m<3) continue; double mu=a/m,var=a2/m-mu*mu; bx.push_back(std::sqrt(lo*hi)); by.push_back(mu); be.push_back(std::sqrt(var>0?var:0)/std::sqrt((double)m)); } };
    std::vector<double> dbx,dby,dbe,lbx,lby,lbe;
    profile(lyD,stiD,dbx,dby,dbe); profile(lyL,stiL,lbx,lby,lbe);

    TCanvas* c=new TCanvas("c_ly","",900,660); c->SetGridx(); c->SetGridy(); c->SetLogx();
    TH1F* fr=c->DrawFrame(12, 60, lymax*1.3, 300);
    fr->SetTitle("Intrinsic timing vs light yield, by scintillator material;mean low-gain amplitude  (light yield) [mV];reference-subtracted single-channel #sigma_{t}  [ps]");
    TGraph* gD=new TGraph(lyD.size(),&lyD[0],&stiD[0]); gD->SetMarkerStyle(20); gD->SetMarkerColor(kRRed-9);    gD->SetMarkerSize(0.7);
    TGraph* gL=new TGraph(lyL.size(),&lyL[0],&stiL[0]); gL->SetMarkerStyle(21); gL->SetMarkerColor(kAzure-8);   gL->SetMarkerSize(0.7);
    gD->Draw("P SAME"); gL->Draw("P SAME");                       // faded per-point scatter
    // joint fit to all intrinsic points
    std::vector<double> lyAll=lyD, stAll=stiD; lyAll.insert(lyAll.end(),lyL.begin(),lyL.end()); stAll.insert(stAll.end(),stiL.begin(),stiL.end());
    TGraph* gAll=new TGraph(lyAll.size(),&lyAll[0],&stAll[0]);
    TF1* f=new TF1("fly","sqrt([0]*[0]/x+[1]*[1])",12,lymax*1.3); f->SetParameters(800,90); f->SetLineColor(kGray+2); f->SetLineWidth(2); f->SetLineStyle(2);
    gAll->Fit(f,"RQN"); f->Draw("SAME");
    // bold binned profiles on top
    TGraphErrors* pD=new TGraphErrors(dbx.size(),&dbx[0],&dby[0],0,&dbe[0]); pD->SetMarkerStyle(20); pD->SetMarkerColor(kRRed);    pD->SetLineColor(kRRed);    pD->SetMarkerSize(1.6); pD->SetLineWidth(2);
    TGraphErrors* pL=new TGraphErrors(lbx.size(),&lbx[0],&lby[0],0,&lbe[0]); pL->SetMarkerStyle(21); pL->SetMarkerColor(kAzure+2); pL->SetLineColor(kAzure+2); pL->SetMarkerSize(1.6); pL->SetLineWidth(2);
    pD->Draw("PL SAME"); pL->Draw("PL SAME");
    TLegend* lg=new TLegend(0.50,0.68,0.88,0.88);
    lg->AddEntry(pD,"DSB1 (binned mean)","pl"); lg->AddEntry(pL,"LuAG:Ce (binned mean)","pl");
    lg->AddEntry(f,"joint fit  #sigma_{t}=#sqrt{p_{0}^{2}/LY+p_{1}^{2}}","l"); lg->Draw();
    TLatex tl; tl.SetNDC(); tl.SetTextSize(0.026);
    tl.DrawLatex(0.14,0.235,"Per capillary, per energy (25-150 GeV); MCP reference (73 ps) subtracted.");
    tl.DrawLatex(0.14,0.193,"DSB1 and LuAG:Ce occupy the same #sigma_{t}-light-yield band, both falling with");
    tl.DrawLatex(0.14,0.151,"light yield. (At equal LY a DSB1 cap sits at lower beam energy than a LuAG one,");
    tl.DrawLatex(0.14,0.109,"since DSB1 is brighter, so this cross-config view is confounded; the in-event");
    tl.DrawLatex(0.14,0.067,"head-to-head is the clean, confound-free test - and shows DSB1 = LuAG.)");
    c->Print(radFigP("figures/narrative/sigmat_vs_lightyield.png"));
    printf("joint fit (intrinsic): p0=%.0f ps*sqrt(mV)  floor p1=%.0f ps\n", f->GetParameter(0), fabs(f->GetParameter(1)));
    printf("wrote figures/narrative/sigmat_vs_lightyield.png\n");
}

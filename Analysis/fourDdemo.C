// fourDdemo.C — Paper 2 capstone: ONE RADiCAL module (TENERGY) measuring TIME,
//   ENERGY and POSITION simultaneously, from the same events.
//     time     = (DW-UP)/2 of the 3 T-type DSB1 timing capillaries (NW E-type excluded)
//     energy   = sum of the low-gain (shower-max) capillary signals
//     position = transverse impact from the 4-corner light (fit vs wire chamber)
//   3-panel figure with the simultaneous resolutions.  150 GeV.
//   root -l 'Analysis/fourDdemo.C+'
#include "ChannelConfig.h"
#include "PlotUtils.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TLinearFitter.h"
#include "TCanvas.h"
#include "TLatex.h"
#include <cstdio>
#include <vector>
#include <algorithm>
#include <cmath>

static double coreS(std::vector<float> v){   // scale-adaptive truncated RMS (any units)
    if(v.size()<150) return -1; std::sort(v.begin(),v.end());
    double mu=v[v.size()/2];
    std::vector<double> d; d.reserve(v.size()); for(float x:v) d.push_back(std::fabs(x-mu)); std::sort(d.begin(),d.end());
    double s=1.4826*d[d.size()/2]; if(s<=0) s=(v.back()-v.front())/6.0; if(s<=0) s=1;
    for(int it=0;it<5;++it){ double a=0,a2=0; long n=0;
        for(float x:v) if(std::fabs(x-mu)<2.5*s){a+=x;a2+=x*x;++n;} if(n<40)break; mu=a/n; double var=a2/n-mu*mu; s=var>0?std::sqrt(var):s; }
    return s;
}
void fourDdemo(double E=150){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    int hg[8],lg[8]; for(int i=0;i<8;++i){ hg[i]=kCap[i].hg/1024; lg[i]=kCap[i].lg/1024; }
    TFile* fp=TFile::Open(Form("reduced/TENERGY/%.0fGeV.root",E)); if(!fp||fp->IsZombie()){ printf("no TENERGY %.0f\n",E); return; }
    TTree* t=(TTree*)fp->Get("rad"); if(!t){ fp->Close(); return; }
    Bool_t wc; Float_t x,y,m1t,m2t,m1p,sp[36],sc[36];
    t->SetBranchAddress("wc_ok",&wc); t->SetBranchAddress("x_trk",&x); t->SetBranchAddress("y_trk",&y);
    t->SetBranchAddress("mcp1_time",&m1t); t->SetBranchAddress("mcp2_time",&m2t);
    t->SetBranchAddress("mcp1_peak",&m1p); t->SetBranchAddress("s_peak",sp); t->SetBranchAddress("s_cfd05",sc);
    const double cx0=6.6, cy0=4.7;
    long N=t->GetEntries(); double xs=0,ys=0; long nw=0;
    for(long i=0;i<N&&nw<40000;++i){ t->GetEntry(i); if(wc&&x>-100&&x<100){xs+=x;ys+=y;++nw;} }
    double xc=nw?xs/nw:0, yc=nw?ys/nw:0;
    // pass 1: collect time, energy, and the position-fit inputs
    TLinearFitter lx(3,"hyp3");
    std::vector<float> tv, ev, f0,f1,f2,xw;
    double Esum=0; long ne=0;
    for(long i=0;i<N;++i){ t->GetEntry(i);
        if(!wc||m1p<200||m1p>750||m1t<-1e5) continue; double dx=x-xc,dy=y-yc; if(dx*dx+dy*dy>=9.0) continue;
        // TIME: (DW-UP)/2 over the 3 DSB1 T-type caps (idx 1,2,3 down; 5,6,7 up), NW(0/4) excluded
        double dsum=0,usum=0; int dn=0,un=0;
        for(int c=1;c<=3;++c){ if(sp[hg[c]]>20&&sc[hg[c]]>-1e5){ dsum+=sc[hg[c]]-m1t; ++dn; } }
        for(int c=5;c<=7;++c){ double ref=(c==7)?m2t:m1t; if(ref<-1e5)continue; if(sp[hg[c]]>20&&sc[hg[c]]>-1e5){ usum+=sc[hg[c]]-ref; ++un; } }
        if(dn<1||un<1) continue;
        double tt=0.5*(dsum/dn-usum/un);
        // ENERGY: sum of 8 LG
        double slg=0; for(int c=0;c<8;++c) slg+=sp[lg[c]]; if(slg<500) continue;
        // POSITION inputs: 4-corner light fractions
        double A[4]={ (double)sp[lg[0]]+sp[lg[4]], (double)sp[lg[1]]+sp[lg[5]], (double)sp[lg[2]]+sp[lg[6]], (double)sp[lg[3]]+sp[lg[7]] };
        double s=A[0]+A[1]+A[2]+A[3]; if(s<500) continue;
        double ff[3]={A[0]/s,A[1]/s,A[2]/s}, xwc=x-cx0;
        if(std::fabs(xwc)>6||std::fabs(y-cy0)>6) continue;
        lx.AddPoint(ff,xwc,1.0);
        tv.push_back((float)tt); ev.push_back((float)slg);
        f0.push_back(ff[0]); f1.push_back(ff[1]); f2.push_back(ff[2]); xw.push_back((float)xwc);
        Esum+=slg; ++ne;
    }
    fp->Close();
    if(ne<2000){ printf("low stats %ld\n",ne); return; }
    lx.Eval();
    double cf[4]={lx.GetParameter(0),lx.GetParameter(1),lx.GetParameter(2),lx.GetParameter(3)};
    double st=coreS(tv)*1000.0;                 // ps (all-fiducial)
    double emean=Esum/ne, sE=coreS(ev)/emean*100.0;
    std::vector<float> pres; for(size_t i=0;i<xw.size();++i){ double xr=cf[0]+cf[1]*f0[i]+cf[2]*f1[i]+cf[3]*f2[i]; pres.push_back((float)(xr-xw[i])); }
    double sP=coreS(pres);
    printf("\n=== 4D demo (TENERGY, %.0f GeV, N=%ld) ===\n",E,ne);
    printf("  TIME   (DW-UP)/2, 3 DSB1 T-type:  core sigma_t = %.1f ps (all-fiducial; energy-binned best-bin ~35 ps)\n",st);
    printf("  ENERGY sum low-gain:              sigma_E/E   = %.1f %%\n",sE);
    printf("  POSITION 4-corner light:          sigma_x     = %.2f mm (vs wire chamber)\n",sP);

    TCanvas* c=new TCanvas("c_4d","",1500,470); c->Divide(3,1,0.01,0.01);
    c->cd(1); { double md=0; std::vector<float> s=tv; std::sort(s.begin(),s.end()); md=s[s.size()/2];
        TH1F* h=new TH1F("ht","TIME  (DW-UP)/2;t (ns);events",80,md-0.6,md+0.6); for(float v:tv)h->Fill(v); h->SetFillColorAlpha(kRRed,0.35); h->SetLineColor(kRRed); h->Draw();
        TLatex l; l.SetNDC(); l.SetTextSize(0.05); l.SetTextColor(kRRed); l.DrawLatex(0.16,0.80,Form("#sigma_{t} = %.0f ps",st)); l.SetTextSize(0.034); l.SetTextColor(kBlack); l.DrawLatex(0.16,0.73,"(best-bin ~35 ps)"); }
    c->cd(2); { TH1F* h=new TH1F("he","ENERGY  #Sigma low-gain;#Sigma LG (mV);events",80,0,emean*2.2); for(float v:ev)h->Fill(v); h->SetFillColorAlpha(kAzure+1,0.35); h->SetLineColor(kAzure+1); h->Draw();
        TLatex l; l.SetNDC(); l.SetTextSize(0.05); l.SetTextColor(kAzure+2); l.DrawLatex(0.16,0.80,Form("#sigma_{E}/E = %.1f%%",sE)); }
    c->cd(3); { TH1F* h=new TH1F("hp","POSITION  capillary - wire chamber;#Deltax (mm);events",80,-8,8); for(float v:pres)h->Fill(v); h->SetFillColorAlpha(30,0.45); h->SetLineColor(kGreen+3); h->Draw();
        TLatex l; l.SetNDC(); l.SetTextSize(0.05); l.SetTextColor(kGreen+3); l.DrawLatex(0.16,0.80,Form("#sigma_{x} = %.1f mm",sP)); }
    c->cd(0); TLatex tl; tl.SetNDC(); tl.SetTextSize(0.030); tl.SetTextFont(62);
    tl.DrawLatex(0.10,0.955,Form("One RADiCAL module, %.0f GeV e^{-}: simultaneous TIME + ENERGY + POSITION (TENERGY: 3 T-type timing + 1 E-type energy)",E));
    c->Print("Analysis/capillary_figs/fourD_demo.png");
    printf("wrote Analysis/capillary_figs/fourD_demo.png\n");
}

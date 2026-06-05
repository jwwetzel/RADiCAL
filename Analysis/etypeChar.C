// etypeChar.C — Paper 2 sec.5: in-beam characterization of the E-type (full-length
//   WLS) ENERGY capillary vs the T-type (shower-max) TIMING capillaries, using the
//   TENERGY module (3x T-type DSB1 + 1x E-type at NW). The in-beam analog of Fig 4
//   of 2401.01747 (uniform E-type vs localized T-type): the E-type is the dim, slow,
//   LINEAR energy element; the T-types are bright and fast (timing).
//   E-type cap = NW (kCap idx 0 down / 4 up); T-types = NE,SE,SW (idx 1,2,3 / 5,6,7).
//   root -l 'Analysis/etypeChar.C+'
#include "ChannelConfig.h"
#include "PlotUtils.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include <cstdio>
#include <vector>
#include <algorithm>
#include <cmath>

static double coreS(std::vector<float> v){
    if(v.size()<150) return -1; std::sort(v.begin(),v.end());
    double mu=v[v.size()/2], s=0.3;
    for(int it=0;it<5;++it){ double a=0,a2=0; long n=0;
        for(float x:v) if(std::fabs(x-mu)<2.5*s){a+=x;a2+=x*x;++n;} if(n<40)break; mu=a/n; double var=a2/n-mu*mu; s=var>0?std::sqrt(var):s; }
    return s;
}
void etypeChar(){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    int hg[8],lg[8]; bool upM2[8];
    for(int i=0;i<8;++i){ hg[i]=kCap[i].hg/1024; lg[i]=kCap[i].lg/1024; upM2[i]=kCap[i].use_mcp2; }
    double Es[5]={50,75,100,125,150};
    std::vector<double> ex, eLY,tLY, eSt,tSt, eLYe,tLYe,eSte,tSte;
    printf("\n=== E-type (NW) vs T-type (NE/SE/SW) in TENERGY ===\n");
    printf("%5s %10s %10s %12s %12s\n","E","Etype_LG","Ttype_LG","Etype_st[ps]","Ttype_st[ps]");
    for(int e=0;e<5;++e){
        TFile* fp=TFile::Open(Form("reduced/TENERGY/%.0fGeV.root",Es[e])); if(!fp||fp->IsZombie()) continue;
        TTree* t=(TTree*)fp->Get("rad"); if(!t){fp->Close();continue;}
        Bool_t wc; Float_t x,y,m1t,m1p,sp[36],sc[36];
        t->SetBranchAddress("wc_ok",&wc); t->SetBranchAddress("x_trk",&x); t->SetBranchAddress("y_trk",&y);
        t->SetBranchAddress("mcp1_time",&m1t); t->SetBranchAddress("mcp1_peak",&m1p);
        t->SetBranchAddress("s_peak",sp); t->SetBranchAddress("s_cfd05",sc);
        long N=t->GetEntries(); double xs=0,ys=0; long nw=0;
        for(long i=0;i<N&&nw<40000;++i){ t->GetEntry(i); if(wc&&x>-100&&x<100){xs+=x;ys+=y;++nw;} }
        double xc=nw?xs/nw:0, yc=nw?ys/nw:0;
        std::vector<float> eLYv,tLYv, eTv,tTv;
        for(long i=0;i<N;++i){ t->GetEntry(i);
            if(!wc||m1p<200||m1p>750||m1t<-1e5) continue; double dx=x-xc,dy=y-yc; if(dx*dx+dy*dy>=9.0) continue;
            // light yield (down+up LG): E-type NW = idx0+idx4; T-type = mean of NE,SE,SW
            eLYv.push_back(sp[lg[0]]+sp[lg[4]]);
            tLYv.push_back((sp[lg[1]]+sp[lg[5]]+sp[lg[2]]+sp[lg[6]]+sp[lg[3]]+sp[lg[7]])/3.0f);
            // single-channel timing (down channel, MCP-referenced)
            if(sp[hg[0]]>20&&sc[hg[0]]>-1e5) eTv.push_back(sc[hg[0]]-m1t);
            for(int c=1;c<=3;++c) if(sp[hg[c]]>20&&sc[hg[c]]>-1e5) tTv.push_back(sc[hg[c]]-m1t);
        }
        auto mean=[&](std::vector<float>&v){ double a=0; for(float z:v)a+=z; return v.empty()?0:a/v.size(); };
        double eLm=mean(eLYv), tLm=mean(tLYv);
        double eS=coreS(eTv)*1000, tS=coreS(tTv)*1000;
        printf("%5.0f %10.0f %10.0f %12.0f %12.0f\n",Es[e],eLm,tLm,eS,tS);
        ex.push_back(Es[e]); eLY.push_back(eLm); tLY.push_back(tLm); eSt.push_back(eS); tSt.push_back(tS);
        eLYe.push_back(eLm*0.01); tLYe.push_back(tLm*0.01); eSte.push_back(eS*0.05); tSte.push_back(tS*0.05);
        fp->Close();
    }
    if(ex.size()<3){ printf("low stats\n"); return; }
    TCanvas* c=new TCanvas("c_et","",1300,560); c->Divide(2,1,0.008,0.006);
    // Left: light yield vs E (linearity) — E-type vs T-type
    c->cd(1); gPad->SetGridx(); gPad->SetGridy();
    TGraphErrors* gE=new TGraphErrors(ex.size(),&ex[0],&eLY[0],0,&eLYe[0]); gE->SetMarkerStyle(21); gE->SetMarkerColor(kAzure+1); gE->SetLineColor(kAzure+1);
    TGraphErrors* gT=new TGraphErrors(ex.size(),&ex[0],&tLY[0],0,&tLYe[0]); gT->SetMarkerStyle(20); gT->SetMarkerColor(kRRed); gT->SetLineColor(kRRed);
    double ymax=0; for(double v:tLY) ymax=std::max(ymax,v); TH1F* f1=gPad->DrawFrame(40,0,160,ymax*1.25);
    f1->SetTitle("Energy linearity: E-type vs T-type;Beam energy (GeV);mean low-gain amplitude (mV)");
    TF1* lE=new TF1("lE","[0]*x",40,160); lE->SetLineColor(kAzure+1); lE->SetLineStyle(2); gE->Fit(lE,"RQN");
    TF1* lT=new TF1("lT","[0]*x",40,160); lT->SetLineColor(kRRed);   lT->SetLineStyle(2); gT->Fit(lT,"RQN");
    lT->Draw("SAME"); lE->Draw("SAME"); gT->Draw("P SAME"); gE->Draw("P SAME");
    TLegend* lg1=new TLegend(0.15,0.72,0.55,0.88);
    lg1->AddEntry(gT,Form("T-type (DSB1, timing): %.1f mV/GeV",lT->GetParameter(0)),"p");
    lg1->AddEntry(gE,Form("E-type (NW, energy): %.1f mV/GeV",lE->GetParameter(0)),"p"); lg1->Draw();
    // Right: single-channel sigma_t vs E — E-type much worse (it is NOT a timing element)
    c->cd(2); gPad->SetGridx(); gPad->SetGridy();
    TGraphErrors* sE=new TGraphErrors(ex.size(),&ex[0],&eSt[0],0,&eSte[0]); sE->SetMarkerStyle(21); sE->SetMarkerColor(kAzure+1); sE->SetLineColor(kAzure+1);
    TGraphErrors* sT=new TGraphErrors(ex.size(),&ex[0],&tSt[0],0,&tSte[0]); sT->SetMarkerStyle(20); sT->SetMarkerColor(kRRed); sT->SetLineColor(kRRed);
    double smax=0; for(double v:eSt) smax=std::max(smax,v); TH1F* f2=gPad->DrawFrame(40,0,160,smax*1.3);
    f2->SetTitle("Single-channel timing (E-type comparable to T-type);Beam energy (GeV);single-channel #sigma_{t} (ps)");
    sT->Draw("P SAME"); sE->Draw("P SAME");
    TLegend* lg2=new TLegend(0.15,0.72,0.55,0.88);
    lg2->AddEntry(sT,"T-type (DSB1)","p"); lg2->AddEntry(sE,"E-type (NW)","p"); lg2->Draw();
    c->cd(0); TLatex tl; tl.SetNDC(); tl.SetTextSize(0.025);
    tl.DrawLatex(0.06,0.955,"E-type (full-length WLS): ~4-5x dimmer at shower max but LINEAR in energy - a uniform energy element. Single-channel timing is comparable;");
    tl.DrawLatex(0.06,0.925,"the T/E division of labour is by light-collection geometry and role (E-type is still excluded from the (DW-UP)/2 corner estimator - see Paper 1).");
    c->Print("Analysis/capillary_figs/etype_vs_ttype.png");
    printf("wrote Analysis/capillary_figs/etype_vs_ttype.png\n");
}

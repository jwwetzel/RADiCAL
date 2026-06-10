// ============================================================================
// hgLgPlot.C — 8-channel HG_peak vs LG_peak panels for one build, with the
//              spike-cut + robust main-line fit drawn. Diagnostic for why a
//              build's HG/LG calibration succeeds or fails.
//   ROOT_INCLUDE_PATH=lib/waveform:lib/io:lib/physics:lib/viz root -l -b -q \
//     'analyze/studies/hgLgPlot.C+("MIXED","data/2023/raw/RUN2941.root",50)'
// ============================================================================
#include "BuildConfig.h"
#include "WaveformUtils.h"
#include "PlotUtils.h"
#include "TChain.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include "DataPaths.h"
#include "FigPaths.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>

struct P { float lg, hg, shp; };

void hgLgPlot(const char* build, const char* rawPath, double E=50, long maxev=120000){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0); gStyle->SetPalette(kBird);
    rad::BuildConfig cfg = rad::BuildConfig::Load(radConfig(build).Data());
    if(!cfg.valid()){ printf("config load failed\n"); return; }

    TH2F* h[8]; std::vector<P> pts[8];
    for(int i=0;i<8;++i){ h[i]=new TH2F(Form("h%d",i),"",110,0,800,110,0,900); h[i]->SetDirectory(0); }
    TChain ch("pulse"); ch.Add(rawPath); TTreeReader rd(&ch);
    TTreeReaderArray<float> amp(rd,"amplitude"), tim(rd,"timevalue");
    long cnt=0;
    while(rd.Next() && cnt<maxev){ ++cnt; const float* T=&tim[0]; const float* A=&amp[0];
        for(int i=0;i<cfg.nend;++i){ const rad::EndMap& c=cfg.end[i];
            Pulse hg=ExtractPulse(T+c.hg_t,A+c.hg,0.20f,5.f);
            Pulse lg=ExtractPulse(T+c.lg_t,A+c.lg,0.20f,5.f);
            h[i]->Fill(lg.peak,hg.peak);
            if(hg.peak>30.f&&hg.peak<700.f&&lg.peak>10.f) pts[i].push_back({lg.peak,hg.peak,hg.charge/hg.peak}); }
    }

    // robust fit per channel (spike cut 0.6*median shape, then 3-sigma trim)
    double A8[8]={0},B8[8]={0}; long N8[8]={0},Spk[8]={0};
    auto lfit=[&](std::vector<std::pair<float,float>>&w,double&a,double&b)->long{
        double sx=0,sy=0,sxx=0,sxy=0; long n=w.size(); if(n<50)return 0;
        for(auto&p:w){sx+=p.first;sy+=p.second;sxx+=p.first*p.first;sxy+=p.first*p.second;}
        double d=n*sxx-sx*sx; if(std::fabs(d)<1e-9)return 0; b=(n*sxy-sx*sy)/d; a=(sy-b*sx)/n; return n; };
    for(int i=0;i<cfg.nend;++i){ std::vector<float> sh; for(auto&p:pts[i]) sh.push_back(p.shp);
        double med=0; if(!sh.empty()){std::sort(sh.begin(),sh.end()); med=sh[sh.size()/2];}
        std::vector<std::pair<float,float>> v; for(auto&p:pts[i]){ if(p.shp>0.6*med) v.push_back({p.lg,p.hg}); else ++Spk[i]; }
        double a=0,b=0; long n=lfit(v,a,b);
        if(n>0){double s2=0;for(auto&p:v){double r=p.second-(a+b*p.first);s2+=r*r;}double rms=std::sqrt(s2/n);
            std::vector<std::pair<float,float>> k;for(auto&p:v)if(std::fabs(p.second-(a+b*p.first))<3*rms)k.push_back(p); lfit(k,a,b); n=k.size();}
        A8[i]=a;B8[i]=b;N8[i]=n; }

    TCanvas* c=new TCanvas("c_hl","",1400,720); c->Divide(4,2,0.001,0.001);
    for(int i=0;i<cfg.nend;++i){ c->cd(i+1); gPad->SetLogz();
        gPad->SetRightMargin(0.06); gPad->SetLeftMargin(0.15); gPad->SetTopMargin(0.12);
        h[i]->GetXaxis()->SetTitle("LG peak (mV)"); h[i]->GetYaxis()->SetTitle("HG peak (mV)");
        h[i]->GetXaxis()->SetTitleSize(0.05); h[i]->GetYaxis()->SetTitleSize(0.05);
        h[i]->GetYaxis()->SetTitleOffset(1.35);   // keep the rotated title inside the pad
        h[i]->Draw("COL");
        TF1* f=new TF1(Form("f%d",i),"[0]+[1]*x",0,800); f->SetParameters(A8[i],B8[i]);
        f->SetLineColor(kRed+1); f->SetLineWidth(2); f->Draw("SAME");
        bool bad = (B8[i]<1.0||B8[i]>8.0);
        TLatex tx; tx.SetNDC(); tx.SetTextFont(42); tx.SetTextSize(0.058);
        tx.SetTextColor(bad?kRed+1:kBlack);
        tx.DrawLatex(0.16,0.905,Form("%s  HG=%.0f+%.2f#upointLG", cfg.end[i].name.c_str(), A8[i], B8[i]));   // in the top-margin strip, clear of the frame
        tx.SetTextSize(0.045); tx.SetTextColor(kGray+3);
        tx.DrawLatex(0.18,0.72,Form("n=%ld  spikes=%ld%s", N8[i], Spk[i], bad?"  DEGENERATE":""));           // below the 820 mV saturation band
    }
    // paper convention (format pass 2026-06-09): no internal super-title; the LaTeX caption carries it
    gSystem->mkdir("figures",kTRUE);
    c->Print(radFigP(Form("figures/hglg_panel_%s.png",build)));
    printf("[%s] ", build); for(int i=0;i<cfg.nend;++i) printf("%s:%.2f ",cfg.end[i].name.c_str(),B8[i]); printf("\n");
    printf("wrote figures/hglg_panel_%s.png\n",build);
}

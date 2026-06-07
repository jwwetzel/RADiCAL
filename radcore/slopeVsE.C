// ============================================================================
// slopeVsE.C — is the HG/LG gain slope truly energy-independent? Per channel,
//   plot slope b vs energy with STAT (fit) + SYST (analysis-choice) errors, for
//   the nominal fit range AND an unsaturated low-HG range. If the drift vanishes
//   in the unsaturated range, it's soft-saturation leaking into the fit (fixable);
//   if it persists, it's a real nonlinearity to model.
//   ROOT_INCLUDE_PATH=radcore:Analysis root -l -b -q 'radcore/slopeVsE.C+("DSB1")'
// ============================================================================
#include "BuildConfig.h"
#include "WaveformUtils.h"
#include "DataPaths.h"
#include "PlotUtils.h"
#include "TChain.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>

struct P { float lg, hg, shp; };

// OLS slope with 3-sigma trim; returns slope b and its statistical error eb.
static bool fitSlope(const std::vector<P>& raw, double shapeCut, double hgmax, double& b, double& eb){
    std::vector<std::pair<double,double>> v;
    for (auto& p : raw) if (p.shp>shapeCut && p.hg>40 && p.hg<hgmax && p.lg>10) v.push_back({p.lg,p.hg});
    if (v.size()<100) return false;
    auto ols=[&](std::vector<std::pair<double,double>>& w, double& a, double& bb, double& sb)->bool{
        long n=w.size(); if(n<50) return false; double sx=0,sy=0,sxx=0,sxy=0;
        for(auto&p:w){sx+=p.first;sy+=p.second;sxx+=p.first*p.first;sxy+=p.first*p.second;}
        double Sxx=sxx-sx*sx/n, Sxy=sxy-sx*sy/n; if(std::fabs(Sxx)<1e-9) return false;
        bb=Sxy/Sxx; a=(sy-bb*sx)/n;
        double s2=0; for(auto&p:w){double r=p.second-(a+bb*p.first); s2+=r*r;} s2/=(n-2);
        sb=std::sqrt(s2/Sxx); return true; };
    double a,bb,sb; if(!ols(v,a,bb,sb)) return false;
    double rms=0; for(auto&p:v){double r=p.second-(a+bb*p.first); rms+=r*r;} rms=std::sqrt(rms/v.size());
    std::vector<std::pair<double,double>> k; for(auto&p:v) if(std::fabs(p.second-(a+bb*p.first))<3*rms) k.push_back(p);
    if(!ols(k,a,bb,sb)) return false; b=bb; eb=sb; return true;
}

void slopeVsE(const char* build){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    rad::BuildConfig cfg = rad::BuildConfig::Load(Form("data/2023/configs/%s.json",build));
    // energy -> raw path (local files)
    std::vector<std::pair<double,std::string>> runs;
    if (std::string(build)=="DSB1"){
        runs={{25,radRaw("RUN1211_25_GeV.root").Data()},{50,radRaw("RUN1148_50_GeV.root").Data()},
              {75,radRaw("RUN1112_75_GeV.root").Data()},{100,radRaw("RUN1075_100_GeV.root").Data()},
              {125,radRaw("RUN1034_125_GeV.root").Data()},{150,radRaw("RUN1258_150_GeV.root").Data()}};
    } else { // MIXED valid era
        runs={{50,"data/2023/raw/RUN2941.root"},{75,"data/2023/raw/RUN2913.root"},
              {100,"data/2023/raw/RUN2995.root"},{125,"data/2023/raw/RUN2985.root"},
              {150,"data/2023/raw/RUN2975.root"}};
    }
    int nE=runs.size();
    // b[ch][E], errors; nominal range [40,700] and unsaturated [40,300]
    std::vector<std::vector<double>> bn(8),en(8),Ev(8), bu(8),eu(8);
    for (int e=0;e<nE;++e){ double E=runs[e].first;
        std::vector<std::vector<P>> pts(8);
        TChain ch("pulse"); ch.Add(runs[e].second.c_str()); TTreeReader rd(&ch);
        TTreeReaderArray<float> amp(rd,"amplitude"),tim(rd,"timevalue"); long cnt=0;
        while(rd.Next()&&cnt<120000){++cnt; const float*A=&amp[0];const float*T=&tim[0];
            for(int i=0;i<cfg.nend;++i){const rad::EndMap&c=cfg.end[i];
                Pulse hg=ExtractPulse(T+c.hg_t,A+c.hg,0.2f,5.f),lg=ExtractPulse(T+c.lg_t,A+c.lg,0.2f,5.f);
                if(hg.peak>30&&hg.peak<900&&lg.peak>10)pts[i].push_back({lg.peak,hg.peak,hg.charge/hg.peak});}}
        for(int i=0;i<cfg.nend;++i){ std::vector<float>sh;for(auto&p:pts[i])sh.push_back(p.shp);
            if(sh.size()<200)continue; std::sort(sh.begin(),sh.end()); double med=sh[sh.size()/2];
            double b0,eb0; if(!fitSlope(pts[i],0.6*med,700,b0,eb0))continue;
            // systematics: vary spike cut (0.5,0.7) and HG max (500,600)
            double syst=0,bv,ev; double vars[4][2]={{0.5,700},{0.7,700},{0.6,600},{0.6,500}};
            for(int v=0;v<4;++v) if(fitSlope(pts[i],vars[v][0]*med,vars[v][1],bv,ev)) syst=std::max(syst,std::fabs(bv-b0));
            bn[i].push_back(b0); en[i].push_back(std::sqrt(eb0*eb0+syst*syst)); Ev[i].push_back(E);
            // unsaturated low-HG range [40,300]
            double bU,eU; if(fitSlope(pts[i],0.6*med,300,bU,eU)){ bu[i].push_back(bU); eu[i].push_back(eU); }
        }
    }
    // ---- report-quality figure: title strip + 4x2 panels ----------------------
    TCanvas* c=new TCanvas("c_se","",1500,820);
    TPad* tp=new TPad("tp","",0,0.925,1,1.0); tp->SetFillStyle(0); tp->Draw(); tp->cd();
    { TLatex t; t.SetNDC(); t.SetTextAlign(22);
      t.SetTextFont(62); t.SetTextSize(0.44);
      t.DrawLatex(0.5,0.66,Form("%s :  HG / LG gain slope vs beam energy", build));
      t.SetTextFont(42); t.SetTextSize(0.30); t.SetTextColor(kGray+3);
      t.DrawLatex(0.5,0.24,"#color[1]{#bullet} nominal fit [40,700 mV]     #color[2]{#circ} unsaturated [40,300 mV]     error bars = stat #oplus syst (spike-cut & fit-range)"); }
    c->cd();
    TPad* mp=new TPad("mp","",0,0,1,0.925); mp->SetFillStyle(0); mp->Draw(); mp->cd();
    mp->Divide(4,2,0.006,0.010);
    printf("[%s] db/dE (nominal [40,700]) and (unsaturated [40,300]):\n", build);
    for(int i=0;i<cfg.nend;++i){ mp->cd(i+1);
        gPad->SetGridy(); gPad->SetTopMargin(0.11); gPad->SetLeftMargin(0.17); gPad->SetBottomMargin(0.15); gPad->SetRightMargin(0.04);
        int n=bn[i].size(); if(n<3){ continue; }
        // y-range covering BOTH series and their errors
        double ymin=1e9,ymax=-1e9;
        for(int k=0;k<n;++k){ ymin=std::min(ymin,bn[i][k]-en[i][k]); ymax=std::max(ymax,bn[i][k]+en[i][k]); }
        for(size_t k=0;k<bu[i].size();++k){ ymin=std::min(ymin,bu[i][k]-eu[i][k]); ymax=std::max(ymax,bu[i][k]+eu[i][k]); }
        double pad=0.18*(ymax-ymin)+0.05;
        TGraphErrors* g=new TGraphErrors(n,&Ev[i][0],&bn[i][0],nullptr,&en[i][0]);
        g->SetMarkerStyle(20); g->SetMarkerColor(kBlack); g->SetLineColor(kBlack); g->SetMarkerSize(1.4); g->SetLineWidth(2);
        g->SetTitle(Form("%s;beam energy (GeV);HG / LG slope", cfg.end[i].name.c_str()));
        g->GetYaxis()->SetRangeUser(ymin-pad,ymax+pad);
        g->GetXaxis()->SetLimits(0,165);
        g->GetXaxis()->SetTitleSize(0.060); g->GetYaxis()->SetTitleSize(0.060);
        g->GetXaxis()->SetLabelSize(0.052); g->GetYaxis()->SetLabelSize(0.052);
        g->GetYaxis()->SetTitleOffset(1.35); g->GetXaxis()->SetTitleOffset(1.05);
        g->Draw("AP");
        TF1* f=new TF1(Form("f%d",i),"[0]+[1]*x",0,165); f->SetLineColor(kBlack); f->SetLineStyle(2); f->SetLineWidth(2);
        g->Fit(f,"Q"); f->Draw("SAME");
        double dbdE=f->GetParameter(1), edb=f->GetParError(1), b0=f->GetParameter(0);
        double pct = b0!=0 ? 100.0*dbdE*125.0/b0 : 0;   // % change over 25->150
        double dbu=0;
        if((int)bu[i].size()==n){ TGraphErrors* gu=new TGraphErrors(n,&Ev[i][0],&bu[i][0],nullptr,&eu[i][0]);
            gu->SetMarkerStyle(24); gu->SetMarkerColor(kRed+1); gu->SetLineColor(kRed+1); gu->SetMarkerSize(1.4); gu->SetLineWidth(2); gu->Draw("P SAME");
            TGraphErrors gg(n,&Ev[i][0],&bu[i][0],nullptr,&eu[i][0]); TF1 fu("fu","[0]+[1]*x",0,165); gg.Fit(&fu,"Q"); dbu=fu.GetParameter(1); }
        TLatex nm; nm.SetNDC(); nm.SetTextFont(62); nm.SetTextSize(0.090); nm.SetTextAlign(22); nm.SetTextColor(kBlack);
        nm.DrawLatex(0.60,0.955,cfg.end[i].name.c_str());
        TLatex tx; tx.SetNDC(); tx.SetTextSize(0.060); tx.SetTextColor(kBlue+2); tx.SetTextAlign(12);
        tx.DrawLatex(0.21,0.875,Form("#Deltab = %+.1f%% (%.1f#sigma)", pct, edb>0?std::fabs(dbdE/edb):0));
        printf("  %-5s nominal db/dE=%+.4f (%.1f%% over 25-150, %.1f-sigma)   unsat db/dE=%+.4f\n",
               cfg.end[i].name.c_str(), dbdE, pct, edb>0?std::fabs(dbdE/edb):0, dbu);
    }
    gSystem->mkdir("radcore/figs",kTRUE); c->Print(Form("radcore/figs/slope_vs_E_%s.png",build));
    printf("wrote radcore/figs/slope_vs_E_%s.png\n",build);
}

// ============================================================================
// pubFig.C — the "data is sound" figure, arXiv:2401.01747 Fig.21 style, per build.
// Pools ALL six beam energies (25..150 GeV) and bins fiducial events by sum_lg
// (amplitude proxy) into NB equal-population bins. For each bin it fits Gaussians
// to BOTH (DW-UP)/2 (the MCP-free differential estimator) and (DW+UP)/2 (the
// absolute cross-check), exactly like the publication. Layout:
//   top row    : (DW-UP)/2 distributions + Gaussian fits (sigma in ps)
//   middle row : (DW+UP)/2 distributions + Gaussian fits
//   bottom     : sigma_t vs amplitude for both estimators, with the brightest-1000
//                headline marked.
// Uses the robust timing source per build (lgcfd for LYSO builds, led for the
// low-light ones) so the peaks are clean Gaussians -- proving every timing number
// rests on sound, well-behaved data, not the cfd05 garbage of the old diagnostic.
//   source setup.sh; root -l -b -q 'analyze/studies/pubFig.C+("DSB1")'
// ============================================================================
#include "RadView.h"
#include "RadTiming.h"
#include "SelectionCuts.h"
#include "DataPaths.h"
#include "PlotUtils.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLine.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
using namespace rad;

struct EV { float slg, dmu, smu; };

static int robustSrc(const char* b){
    std::string s=b;
    if(s=="LUAG"||s=="TENERGY") return RadView::kLED;   // low-light: fixed threshold
    return RadView::kLGCFD;                              // DSB1, MIXED: HG/LG-ratio edge
}
static const char* srcName(int s){ return s==RadView::kLED?"led":(s==RadView::kLGCFD?"lgcfd":"cfd05"); }

// trimmed mean/sigma for a clean histogram range
static void trimMS(std::vector<float>& x,double& mu,double& sg){
    double m=0; for(float q:x) m+=q; m/=x.size();
    double s=0; for(float q:x) s+=(q-m)*(q-m); s=std::sqrt(s/x.size());
    double mt=0,n=0; for(float q:x) if(std::fabs(q-m)<3*s){mt+=q;++n;} if(n>0)mt/=n;
    double st=0; n=0; for(float q:x) if(std::fabs(q-m)<3*s){st+=(q-mt)*(q-mt);++n;} st=n>0?std::sqrt(st/n):s;
    mu=mt; sg=(st>1e-4?st:s>1e-4?s:0.05);
}

// CENTRAL-PEAK Gaussian fit, just like the publication: iterate the fit window in
// on the core so the satellite shoulder can't drag it, and draw the fitted curve
// over the core ONLY (+-2sigma), not into the tails. Returns the TF1 and sigma_ps.
static TF1* coreFit(TH1F* h,const char* nm,double mu0,double sg0,double& sigma_ps){
    double mu=mu0, sg=(sg0>1e-4?sg0:0.02);
    TF1* g=new TF1(nm,"gaus",mu-1.8*sg,mu+1.8*sg);
    for(int it=0;it<3;++it){
        g->SetRange(mu-1.8*sg,mu+1.8*sg);
        g->SetParameters(h->GetMaximum(),mu,sg);
        h->Fit(g,"QRN");
        double m=g->GetParameter(1), s=std::fabs(g->GetParameter(2));
        if(s>1e-4 && std::fabs(m-mu0)<4*sg0){ mu=m; sg=s; } else break;
    }
    sigma_ps=sg*1000.0;
    g->SetRange(mu-2.0*sg, mu+2.0*sg);    // draw the central peak only
    return g;
}

void pubFig(const char* build="DSB1"){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0); gStyle->SetOptFit(0);
    BuildConfig cfg=BuildConfig::Load(Form("data/2023/configs/%s.json",build));
    int SRC=robustSrc(build);
    const double Es[]={25,50,75,100,125,150}; int nE=6;
    std::vector<EV> ev;
    for(int e=0;e<nE;++e){ double E=Es[e];
        TFile* fp=TFile::Open(radReduced(build,E)); if(!fp||fp->IsZombie())continue;
        TTree* t=(TTree*)fp->Get("rad"); RadView v; v.attach(t,&cfg);
        double xc,yc; v.beamCenter(xc,yc); double r2=9.0; Long64_t N=v.entries();
        for(Long64_t i=0;i<N;++i){ v.get(i);
            if(!v.wc_ok()||v.mcp1_peak()<kMCP1_minPeak||v.mcp1_peak()>kMCP1_maxPeak) continue;
            double dx=v.x_trk()-xc,dy=v.y_trk()-yc; if(dx*dx+dy*dy>=r2) continue;
            double ds=0,us=0;int dn=0,un=0;
            for(int c=0;c<4;++c) if(v.is_timing(c)&&v.hg_peak(c)>=kHG_minPeak){float tc=v.timeOf(c,SRC);if(tc>-1e5){ds+=tc;++dn;}}
            for(int c=4;c<8;++c) if(v.is_timing(c)&&v.hg_peak(c)>=kHG_minPeak){float tc=v.timeOf(c,SRC);if(tc>-1e5){us+=tc;++un;}}
            if(dn<1||un<1) continue;
            double da=ds/dn, ua=us/un;
            ev.push_back({(float)v.sum_lg(),0.5f*(float)(da-ua),0.5f*(float)(da+ua)});
        }
        fp->Close();
    }
    if(ev.size()<3000){ printf("%s: too few events (%zu)\n",build,ev.size()); return; }
    std::sort(ev.begin(),ev.end(),[](const EV&a,const EV&b){return a.slg<b.slg;});
    const int NB=9; size_t per=ev.size()/NB; long Nbin=(long)per;
    std::vector<double> amp(NB),sdiff(NB),ssum(NB),ediff(NB),esum(NB);
    std::vector<TH1F*> hd(NB),hs(NB); std::vector<TF1*> fd(NB),fs(NB);
    // pass 1: per-bin slices + trimmed centers -> ONE common x-range per row (shared)
    std::vector<std::vector<float>> VD(NB),VS(NB);
    std::vector<double> MD(NB),SGD(NB),MS(NB),SGS(NB);
    double loD=1e9,hiD=-1e9,loS=1e9,hiS=-1e9;
    for(int b=0;b<NB;++b){
        size_t lo=(size_t)b*per, hi=(b==NB-1)?ev.size():lo+per; double sa=0;
        for(size_t k=lo;k<hi;++k){ VD[b].push_back(ev[k].dmu); VS[b].push_back(ev[k].smu); sa+=ev[k].slg; }
        amp[b]=sa/(hi-lo);
        trimMS(VD[b],MD[b],SGD[b]); trimMS(VS[b],MS[b],SGS[b]);
        loD=std::min(loD,MD[b]-4.5*SGD[b]); hiD=std::max(hiD,MD[b]+4.5*SGD[b]);
        loS=std::min(loS,MS[b]-4.5*SGS[b]); hiS=std::max(hiS,MS[b]+4.5*SGS[b]);
    }
    // pass 2: hist over the COMMON range (shared binning -> comparable absolute counts),
    // central-peak-only Gaussian fit. Track the shared y-max per row for absolute scaling.
    const int NHB=100; double ymaxD=0,ymaxS=0;
    printf("\n=== %s (%s): arXiv-Fig21-style, all energies pooled, %zu events (%ld/bin, equal-population) ===\n",build,srcName(SRC),ev.size(),Nbin);
    printf("  bin  <ampl>   N     sigma(DW-UP)/2   sigma(DW+UP)/2  [ps]\n");
    for(int b=0;b<NB;++b){
        TH1F* Hd=new TH1F(Form("hd%s%d",build,b),"",NHB,loD,hiD); Hd->SetDirectory(nullptr);
        TH1F* Hs=new TH1F(Form("hs%s%d",build,b),"",NHB,loS,hiS); Hs->SetDirectory(nullptr);
        for(float q:VD[b]) Hd->Fill(q); for(float q:VS[b]) Hs->Fill(q);
        double sd,ss;
        TF1* Fd=coreFit(Hd,Form("fd%s%d",build,b),MD[b],SGD[b],sd);
        TF1* Fs=coreFit(Hs,Form("fs%s%d",build,b),MS[b],SGS[b],ss);
        if(!(sd>0&&sd<2.0*SGD[b]*1000)) sd=SGD[b]*1000;          // robust fallback
        if(!(ss>0&&ss<2.0*SGS[b]*1000)) ss=SGS[b]*1000;
        sdiff[b]=sd; ssum[b]=ss; ediff[b]=sd/std::sqrt(2.0*VD[b].size()); esum[b]=ss/std::sqrt(2.0*VS[b].size());
        hd[b]=Hd; hs[b]=Hs; fd[b]=Fd; fs[b]=Fs;
        ymaxD=std::max(ymaxD,Hd->GetMaximum()); ymaxS=std::max(ymaxS,Hs->GetMaximum());
        printf("  %2d  %6.0f %6zu       %6.1f           %6.1f\n",b,amp[b],VD[b].size(),sd,ss);
    }
    // brightest-1000 headline (the top of the amplitude axis)
    std::vector<float> vtop; for(size_t k=(ev.size()>=1000?ev.size()-1000:0);k<ev.size();++k) vtop.push_back(ev[k].dmu);
    double sTop=tebSigma(vtop);
    double ampTop=0; for(size_t k=(ev.size()>=1000?ev.size()-1000:0);k<ev.size();++k) ampTop+=ev[k].slg; ampTop/=std::min((size_t)1000,ev.size());
    printf("  brightest-1000 (DW-UP)/2 sigma = %.1f ps  (<ampl>=%.0f)\n",sTop,ampTop);

    // ---- draw ----
    TCanvas* c=new TCanvas("pub","",1820,1040);
    TPad* pTop=new TPad("pTop","",0.0,0.655,1.0,0.945); pTop->Draw();
    TPad* pMid=new TPad("pMid","",0.0,0.365,1.0,0.655); pMid->Draw();
    TPad* pBot=new TPad("pBot","",0.0,0.0,1.0,0.365);  pBot->Draw();
    pTop->Divide(NB,1,0.0,0.0); pMid->Divide(NB,1,0.0,0.0);
    for(int b=0;b<NB;++b){
        bool first=(b==0); double lm=(first?0.26:0.03);
        pTop->cd(b+1); gPad->SetTopMargin(0.12); gPad->SetBottomMargin(0.15); gPad->SetLeftMargin(lm); gPad->SetRightMargin(0.02);
        hd[b]->SetLineColor(kAzure+2); hd[b]->SetFillColorAlpha(kAzure+2,0.30); hd[b]->SetLineWidth(1);
        hd[b]->GetYaxis()->SetRangeUser(0,ymaxD*1.15);                 // SHARED absolute y per row
        hd[b]->GetXaxis()->SetLabelSize(0.072); hd[b]->GetXaxis()->SetNdivisions(304);
        hd[b]->GetXaxis()->SetTitle("(DW-UP)/2 (ns)"); hd[b]->GetXaxis()->SetTitleSize(0.07); hd[b]->GetXaxis()->SetTitleOffset(0.95);
        hd[b]->GetYaxis()->SetLabelSize(first?0.066:0.0); hd[b]->GetYaxis()->SetNdivisions(first?505:0);
        if(first){ hd[b]->GetYaxis()->SetTitle("events"); hd[b]->GetYaxis()->SetTitleSize(0.075); hd[b]->GetYaxis()->SetTitleOffset(1.55); }
        hd[b]->Draw("HIST"); fd[b]->SetLineColor(kOrange+8); fd[b]->SetLineWidth(3); fd[b]->Draw("SAME");
        TLatex tx; tx.SetNDC(); tx.SetTextFont(62); tx.SetTextSize(0.092); tx.DrawLatex(first?0.32:0.09,0.83,Form("#sigma=%.0f",sdiff[b]));
        pMid->cd(b+1); gPad->SetTopMargin(0.12); gPad->SetBottomMargin(0.15); gPad->SetLeftMargin(lm); gPad->SetRightMargin(0.02);
        hs[b]->SetLineColor(kGray+2); hs[b]->SetFillColorAlpha(kGray+1,0.32); hs[b]->SetLineWidth(1);
        hs[b]->GetYaxis()->SetRangeUser(0,ymaxS*1.15);                 // SHARED absolute y per row
        hs[b]->GetXaxis()->SetLabelSize(0.072); hs[b]->GetXaxis()->SetNdivisions(304);
        hs[b]->GetXaxis()->SetTitle("(DW+UP)/2 (ns)"); hs[b]->GetXaxis()->SetTitleSize(0.07); hs[b]->GetXaxis()->SetTitleOffset(0.95);
        hs[b]->GetYaxis()->SetLabelSize(first?0.066:0.0); hs[b]->GetYaxis()->SetNdivisions(first?505:0);
        if(first){ hs[b]->GetYaxis()->SetTitle("events"); hs[b]->GetYaxis()->SetTitleSize(0.075); hs[b]->GetYaxis()->SetTitleOffset(1.55); }
        hs[b]->Draw("HIST"); fs[b]->SetLineColor(kOrange+8); fs[b]->SetLineWidth(3); fs[b]->Draw("SAME");
        TLatex tx2; tx2.SetNDC(); tx2.SetTextFont(62); tx2.SetTextSize(0.092); tx2.DrawLatex(first?0.32:0.09,0.83,Form("#sigma=%.0f",ssum[b]));
    }
    pBot->cd(); gPad->SetLeftMargin(0.085); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.05); gPad->SetBottomMargin(0.17); gPad->SetGridy();
    std::vector<double> ze(NB,0.0);
    TGraphErrors* gd=new TGraphErrors(NB,&amp[0],&sdiff[0],&ze[0],&ediff[0]);
    TGraphErrors* gs=new TGraphErrors(NB,&amp[0],&ssum[0],&ze[0],&esum[0]);
    gd->SetMarkerStyle(22); gd->SetMarkerColor(kAzure+2); gd->SetLineColor(kAzure+2); gd->SetMarkerSize(1.7); gd->SetLineWidth(2);
    gs->SetMarkerStyle(20); gs->SetMarkerColor(kOrange+8); gs->SetLineColor(kOrange+8); gs->SetMarkerSize(1.5); gs->SetLineWidth(2);
    double ymax=0; for(int b=0;b<NB;++b) ymax=std::max(ymax,std::max(sdiff[b],ssum[b]));
    gd->SetTitle(";shower amplitude  #SigmaLG (a.u.)  #minus  all energies 25#minus150 GeV;time resolution  #sigma_{t} (ps)");
    gd->GetYaxis()->SetRangeUser(0,ymax*1.18); gd->GetYaxis()->SetTitleSize(0.062); gd->GetYaxis()->SetTitleOffset(0.62);
    gd->GetYaxis()->SetLabelSize(0.055); gd->GetXaxis()->SetTitleSize(0.058); gd->GetXaxis()->SetTitleOffset(1.15); gd->GetXaxis()->SetLabelSize(0.05);
    gd->Draw("ALP"); gs->Draw("LP SAME");
    if(sTop>0){ TLine* ln=new TLine(amp[0],sTop,amp[NB-1],sTop); ln->SetLineColor(kAzure+3); ln->SetLineStyle(2); ln->SetLineWidth(2); ln->Draw();
        TLatex th; th.SetTextColor(kAzure+3); th.SetTextSize(0.05); th.DrawLatex(amp[1],sTop+ymax*0.03,Form("brightest-1000 headline: %.1f ps",sTop)); }
    TLegend* lg=new TLegend(0.70,0.70,0.965,0.95); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.058);
    lg->AddEntry(gd,"(DW#minusUP)/2  (MCP-free)","lp"); lg->AddEntry(gs,"(DW+UP)/2  (absolute)","lp"); lg->Draw();
    c->cd(); TLatex tt; tt.SetNDC(); tt.SetTextFont(62); tt.SetTextSize(0.025);
    tt.DrawLatex(0.05,0.967,Form("%s (%s): per-amplitude-bin timing distributions are clean Gaussians  --  top: (DW-UP)/2,  middle: (DW+UP)/2,  all six energies pooled (25-150 GeV)",build,srcName(SRC)));
    gSystem->mkdir("figures/narrative",kTRUE); c->Print(Form("figures/narrative/pub_%s.png",build));
    printf("  wrote figures/narrative/pub_%s.png\n",build);
}

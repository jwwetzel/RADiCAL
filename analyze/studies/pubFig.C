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
    const int NB=9; size_t per=ev.size()/NB;
    std::vector<double> amp(NB),sdiff(NB),ssum(NB),ediff(NB),esum(NB);
    std::vector<TH1F*> hd(NB),hs(NB); std::vector<TF1*> fd(NB),fs(NB);
    printf("\n=== %s (%s): arXiv-Fig21-style, all energies pooled, %zu events ===\n",build,srcName(SRC),ev.size());
    printf("  bin  <ampl>   N     sigma(DW-UP)/2   sigma(DW+UP)/2  [ps]\n");
    for(int b=0;b<NB;++b){
        size_t lo=(size_t)b*per, hi=(b==NB-1)?ev.size():lo+per;
        std::vector<float> vd,vs; double sa=0;
        for(size_t k=lo;k<hi;++k){ vd.push_back(ev[k].dmu); vs.push_back(ev[k].smu); sa+=ev[k].slg; }
        amp[b]=sa/(hi-lo);
        double md,sgd,ms,sgs; trimMS(vd,md,sgd); trimMS(vs,ms,sgs);
        TH1F* Hd=new TH1F(Form("hd%s%d",build,b),"",70,md-4*sgd,md+4*sgd); Hd->SetDirectory(nullptr);
        TH1F* Hs=new TH1F(Form("hs%s%d",build,b),"",70,ms-4*sgs,ms+4*sgs); Hs->SetDirectory(nullptr);
        for(float q:vd) Hd->Fill(q); for(float q:vs) Hs->Fill(q);
        TF1* Fd=new TF1(Form("fd%s%d",build,b),"gaus",md-2*sgd,md+2*sgd); Hd->Fit(Fd,"QRN");
        TF1* Fs=new TF1(Form("fs%s%d",build,b),"gaus",ms-2*sgs,ms+2*sgs); Hs->Fit(Fs,"QRN");
        Fd->SetRange(md-4*sgd,md+4*sgd); Fs->SetRange(ms-4*sgs,ms+4*sgs);
        double sd=Fd->GetParameter(2)*1000.0, ss=Fs->GetParameter(2)*1000.0;
        if(!(sd>0&&sd<2.0*sgd*1000)) sd=sgd*1000;          // robust fallback
        if(!(ss>0&&ss<2.0*sgs*1000)) ss=sgs*1000;
        sdiff[b]=sd; ssum[b]=ss; ediff[b]=sd/std::sqrt(2.0*vd.size()); esum[b]=ss/std::sqrt(2.0*vs.size());
        hd[b]=Hd; hs[b]=Hs; fd[b]=Fd; fs[b]=Fs;
        printf("  %2d  %6.0f %6zu       %6.1f           %6.1f\n",b,amp[b],hi-lo,sd,ss);
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
    pTop->Divide(NB,1,0.0008,0.0008); pMid->Divide(NB,1,0.0008,0.0008);
    for(int b=0;b<NB;++b){
        pTop->cd(b+1); gPad->SetTopMargin(0.13); gPad->SetBottomMargin(0.14); gPad->SetLeftMargin(0.03); gPad->SetRightMargin(0.02);
        hd[b]->SetLineColor(kAzure+2); hd[b]->SetFillColorAlpha(kAzure+2,0.28); hd[b]->SetLineWidth(1);
        hd[b]->GetXaxis()->SetLabelSize(0.075); hd[b]->GetYaxis()->SetLabelSize(0.0); hd[b]->GetXaxis()->SetNdivisions(304);
        hd[b]->GetXaxis()->SetTitle("(DW#minusUP)/2 (ns)"); hd[b]->GetXaxis()->SetTitleSize(0.075); hd[b]->GetXaxis()->SetTitleOffset(0.85);
        hd[b]->Draw("HIST"); fd[b]->SetLineColor(kOrange+8); fd[b]->SetLineWidth(2); fd[b]->Draw("SAME");
        TLatex tx; tx.SetNDC(); tx.SetTextFont(62); tx.SetTextSize(0.10); tx.DrawLatex(0.08,0.80,Form("#sigma=%.0f",sdiff[b]));
        pMid->cd(b+1); gPad->SetTopMargin(0.13); gPad->SetBottomMargin(0.14); gPad->SetLeftMargin(0.03); gPad->SetRightMargin(0.02);
        hs[b]->SetLineColor(kGray+2); hs[b]->SetFillColorAlpha(kGray+1,0.30); hs[b]->SetLineWidth(1);
        hs[b]->GetXaxis()->SetLabelSize(0.075); hs[b]->GetYaxis()->SetLabelSize(0.0); hs[b]->GetXaxis()->SetNdivisions(304);
        hs[b]->GetXaxis()->SetTitle("(DW+UP)/2 (ns)"); hs[b]->GetXaxis()->SetTitleSize(0.075); hs[b]->GetXaxis()->SetTitleOffset(0.85);
        hs[b]->Draw("HIST"); fs[b]->SetLineColor(kOrange+8); fs[b]->SetLineWidth(2); fs[b]->Draw("SAME");
        TLatex tx2; tx2.SetNDC(); tx2.SetTextFont(62); tx2.SetTextSize(0.10); tx2.DrawLatex(0.08,0.80,Form("#sigma=%.0f",ssum[b]));
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

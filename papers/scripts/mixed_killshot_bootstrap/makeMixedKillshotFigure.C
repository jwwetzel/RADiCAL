// ============================================================================
// makeMixedKillshotFigure.C — the OFFICIAL Paper-2 MIXED head-to-head figure,
// GATE-6 compliant (replaces the retired mixed_h2h.png and its "0.99, chi2=0.4").
// Recomputes everything from source (same machinery as mixedKillshotBootstrap.C):
//   Panel A: per-capillary timing widths vs E, DSB1 (NE,SW) vs LuAG (NW,SE),
//            GATE-1 map, srCFD, same-module/same-shower label, stat errors.
//   Panel B: ratio sigma_DSB1/sigma_LuAG vs E + the scatter-based 1.04+-0.05 band
//            (NOT the bootstrap-only CI). No 0.99, no chi2, no "indistinguishable".
//   Panel C: method dependence (srCFD headline; cfd05/LED diagnostic).
// Caption text: makeMixedKillshotFigure_CAPTION.txt (this directory).
//   source setup.sh
//   root -l -b -q 'papers/scripts/mixed_killshot_bootstrap/makeMixedKillshotFigure.C+'
// Output: papers/figures/mixed_killshot_bootstrap/mixed_h2h_corrected.{png,pdf}
// ============================================================================
#include "RadView.h"
#include "DataPaths.h"
#include "PlotUtils.h"
#include "SelectionCuts.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TGraphErrors.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>
using namespace rad;

static const int NS=3; static const char* SN[NS]={"srCFD","cfd05","LED"};
static const char* HB[NS]={"hg_lgcfd","hg_cfd05","hg_led"};

struct Pair { std::vector<float> dD,dL; };
struct W { double lo,hi; };
static W cw(std::vector<float>& v){ W w{0,0}; if(v.size()<300)return w;
    double mu=0; for(float x:v)mu+=x; mu/=v.size();
    double sd=0; for(float x:v)sd+=(x-mu)*(x-mu); sd=std::sqrt(sd/v.size()); if(sd<1e-4)sd=0.1;
    for(int it=0;it<6;++it){ double s=0,ss=0;long n=0; for(float x:v)if(std::fabs(x-mu)<2.5*sd){s+=x;ss+=x*x;++n;}
        if(n<100)break; mu=s/n; double q=ss/n-mu*mu; if(q>0)sd=std::sqrt(q); }
    w.lo=mu-2.5*sd; w.hi=mu+2.5*sd; return w; }
static double wrms(std::vector<float>& v,W w,long* nOut=nullptr){ double s=0,ss=0;long n=0;
    for(float x:v){ if(x<w.lo||x>w.hi)continue; s+=x;ss+=x*x;++n; }
    if(nOut)*nOut=n; if(n<50)return -1; double mu=s/n,q=ss/n-mu*mu; return q>0?std::sqrt(q):-1; }

static Pair gather(double E,int src){
    Pair P;
    TString p=radReduced("MIXED",E); if(gSystem->AccessPathName(p.Data())) return P;
    TFile* f=TFile::Open(p.Data()); if(!f||f->IsZombie()){if(f)f->Close();return P;}
    TTree* t=(TTree*)f->Get("rad"); if(!t){f->Close();return P;}
    Bool_t wc; Float_t x,y,mp,hp[8],ht[8];
    t->SetBranchStatus("*",0);
    for(const char* b:{"wc_ok","x_trk","y_trk","mcp1_peak","hg_peak"}) t->SetBranchStatus(b,1);
    t->SetBranchStatus(HB[src],1);
    t->SetBranchAddress("wc_ok",&wc); t->SetBranchAddress("x_trk",&x); t->SetBranchAddress("y_trk",&y);
    t->SetBranchAddress("mcp1_peak",&mp); t->SetBranchAddress("hg_peak",hp); t->SetBranchAddress(HB[src],ht);
    Long64_t N=t->GetEntries(); double sx=0,sy=0; long nc=0;
    for(Long64_t i=0;i<N&&nc<40000;++i){ t->GetEntry(i); if(wc&&x>-100&&x<100){sx+=x;sy+=y;++nc;} }
    if(nc<1000){f->Close();return P;}
    double xc=sx/nc,yc=sy/nc;
    for(Long64_t i=0;i<N;++i){ t->GetEntry(i);
        if(!wc||mp<kMCP1_minPeak||mp>kMCP1_maxPeak) continue;
        double dx=x-xc,dy=y-yc; if(dx*dx+dy*dy>=9.0) continue;
        bool ok=true; for(int s=0;s<7;++s) if(hp[s]<20||ht[s]<=-1e5f){ok=false;break;}
        if(!ok) continue;
        P.dD.push_back(ht[1]-ht[3]);   // DSB1 pair: NE-D - SW-D  (GATE-1 map)
        P.dL.push_back(ht[0]-ht[2]);   // LuAG pair: NW-D - SE-D
    }
    f->Close(); return P;
}

void makeMixedKillshotFigure(){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0); gStyle->SetGridStyle(3); gStyle->SetGridColor(kGray);
    gSystem->mkdir("papers/figures/mixed_killshot_bootstrap",kTRUE);
    const double Es[5]={50,75,100,125,150};
    double sD[5],sL[5],eD[5],eL[5],R[5],eR[5],ze[5]={0,0,0,0,0};
    double Rsrc[NS]={0,0,0};
    printf("\n===== corrected MIXED head-to-head figure (GATE-6 language) =====\n");
    for(int s=0;s<NS;++s){ double m=0; int n=0;
        for(int ie=0;ie<5;++ie){
            Pair P=gather(Es[ie],s); if(P.dD.size()<5000) continue;
            W wD=cw(P.dD), wL=cw(P.dL); long nD=0,nL=0;
            double xD=wrms(P.dD,wD,&nD), xL=wrms(P.dL,wL,&nL);
            if(xD<0||xL<0) continue;
            double r=xD/xL; m+=r; ++n;
            if(s==0){ sD[ie]=xD/std::sqrt(2.0)*1000.0; sL[ie]=xL/std::sqrt(2.0)*1000.0;
                eD[ie]=sD[ie]/std::sqrt(2.0*nD); eL[ie]=sL[ie]/std::sqrt(2.0*nL);
                R[ie]=r; eR[ie]=r*std::sqrt(1.0/(2*nD)+1.0/(2*nL));
                printf("  %3.0f GeV: sigma_DSB1=%5.0f+-%.0f ps  sigma_LuAG=%5.0f+-%.0f ps  R=%.3f+-%.3f (N=%ld)\n",
                       Es[ie],sD[ie],eD[ie],sL[ie],eL[ie],R[ie],eR[ie],nD); } }
        Rsrc[s]=n?m/n:-1; }
    double Rm=0; for(int i=0;i<5;++i)Rm+=R[i]; Rm/=5;
    double Rsc=0; for(int i=0;i<5;++i)Rsc+=(R[i]-Rm)*(R[i]-Rm); Rsc=std::sqrt(Rsc/4.0)/std::sqrt(5.0);
    printf("  summary: mean R=%.3f, scatter-based err=%.3f; methods srCFD/cfd05/LED = %.3f/%.3f/%.3f\n",
           Rm,Rsc,Rsrc[0],Rsrc[1],Rsrc[2]);

    TCanvas* c=new TCanvas("h2hc","",1640,620); c->Divide(3,1,0.005,0.005);
    // --- Panel A: widths vs E ---
    c->cd(1); gPad->SetLeftMargin(0.15);gPad->SetRightMargin(0.03);gPad->SetTopMargin(0.085);gPad->SetBottomMargin(0.13);gPad->SetGridy();
    { double ymx=0,ymn=1e9; for(int i=0;i<5;++i){ymx=std::max({ymx,sD[i],sL[i]});ymn=std::min({ymn,sD[i],sL[i]});}
      TH1F* fr=gPad->DrawFrame(40,ymn-40,160,ymx+90);
      fr->SetTitle(";beam energy E (GeV);per-capillary #sigma_{t} (ps)");
      fr->GetYaxis()->SetTitleSize(0.05);fr->GetYaxis()->SetTitleOffset(1.45);fr->GetXaxis()->SetTitleSize(0.05);
      TGraphErrors* gD=new TGraphErrors(5,Es,sD,ze,eD);
      gD->SetMarkerStyle(20);gD->SetMarkerColor(kRData);gD->SetLineColor(kRData);gD->SetMarkerSize(1.6);gD->SetLineWidth(2);gD->Draw("PL SAME");
      TGraphErrors* gL=new TGraphErrors(5,Es,sL,ze,eL);
      gL->SetMarkerStyle(21);gL->SetMarkerColor(kRGreen+1);gL->SetLineColor(kRGreen+1);gL->SetMarkerSize(1.6);gL->SetLineWidth(2);gL->Draw("PL SAME");
      TLegend* lg=new TLegend(0.38,0.74,0.96,0.90); lg->SetBorderSize(0);lg->SetFillStyle(0);lg->SetTextSize(0.038);
      lg->AddEntry(gD,"DSB1 capillaries (NE, SW)","lp"); lg->AddEntry(gL,"LuAG:Ce capillaries (NW, SE)","lp"); lg->Draw();
      TLatex tx; tx.SetNDC(); tx.SetTextSize(0.034); tx.SetTextColor(kGray+3);
      tx.DrawLatex(0.18,0.20,"same module, same showers, same DRS group;");
      tx.DrawLatex(0.18,0.155,"only the WLS capillary differs (map: pulse-shape confirmed)");
      DrawPadTitle("A  same-shower timing widths (srCFD)"); }
    // --- Panel B: ratio vs E + scatter band ---
    c->cd(2); gPad->SetLeftMargin(0.15);gPad->SetRightMargin(0.03);gPad->SetTopMargin(0.085);gPad->SetBottomMargin(0.13);gPad->SetGridy();
    { TH1F* fr=gPad->DrawFrame(40,0.80,160,1.32);
      fr->SetTitle(";beam energy E (GeV);#sigma_{DSB1} / #sigma_{LuAG}");
      fr->GetYaxis()->SetTitleSize(0.05);fr->GetYaxis()->SetTitleOffset(1.45);fr->GetXaxis()->SetTitleSize(0.05);
      TBox* band=new TBox(40,Rm-Rsc*std::sqrt(5.0),160,Rm+Rsc*std::sqrt(5.0));   // +-1 RMS scatter band
      band->SetFillColorAlpha(kAzure+1,0.18); band->Draw();
      TLine* l1=new TLine(40,1.0,160,1.0); l1->SetLineColor(kGray+2); l1->SetLineStyle(2); l1->SetLineWidth(2); l1->Draw();
      TGraphErrors* gR=new TGraphErrors(5,Es,R,ze,eR);
      gR->SetMarkerStyle(20);gR->SetMarkerColor(kAzure+2);gR->SetLineColor(kAzure+2);gR->SetMarkerSize(1.7);gR->Draw("P SAME");
      TLatex tx; tx.SetNDC(); tx.SetTextSize(0.038);
      tx.DrawLatex(0.18,0.86,Form("mean = %.2f #pm %.2f (energy-period scatter)",Rm,Rsc));
      tx.SetTextSize(0.032); tx.SetTextColor(kGray+3);
      tx.DrawLatex(0.18,0.80,"band: #pm1 RMS of the per-energy ratios;");
      tx.DrawLatex(0.18,0.755,"the spread across energies dominates the");
      tx.DrawLatex(0.18,0.71,"statistical (bootstrap) uncertainty");
      DrawPadTitle("B  width ratio per energy"); }
    // --- Panel C: method dependence ---
    c->cd(3); gPad->SetLeftMargin(0.15);gPad->SetRightMargin(0.03);gPad->SetTopMargin(0.085);gPad->SetBottomMargin(0.13);gPad->SetGridy();
    { TH1F* fr=gPad->DrawFrame(-0.5,0.70,2.5,1.32);
      fr->SetTitle(";timing estimator;5-energy mean  #sigma_{DSB1}/#sigma_{LuAG}");
      fr->GetYaxis()->SetTitleSize(0.05);fr->GetYaxis()->SetTitleOffset(1.45);fr->GetXaxis()->SetTitleSize(0.05);
      for(int s=0;s<NS;++s) fr->GetXaxis()->SetBinLabel(fr->GetXaxis()->FindBin(s),SN[s]);
      fr->GetXaxis()->SetLabelSize(0.065); fr->GetXaxis()->SetNdivisions(3);
      TLine* l1=new TLine(-0.5,1.0,2.5,1.0); l1->SetLineColor(kGray+2);l1->SetLineStyle(2);l1->SetLineWidth(2);l1->Draw();
      double xs[NS]={0,1,2};
      TGraphErrors* gh=new TGraphErrors(1,xs,Rsrc,ze,ze);
      gh->SetMarkerStyle(20);gh->SetMarkerColor(kAzure+2);gh->SetMarkerSize(2.0);gh->Draw("P SAME");
      TGraphErrors* gd=new TGraphErrors(2,xs+1,Rsrc+1,ze,ze);
      gd->SetMarkerStyle(24);gd->SetMarkerColor(kGray+2);gd->SetMarkerSize(1.8);gd->Draw("P SAME");
      TLatex tx; tx.SetNDC(); tx.SetTextSize(0.033); tx.SetTextColor(kGray+3);
      tx.DrawLatex(0.18,0.86,"srCFD (filled) = headline: least sensitive to");
      tx.DrawLatex(0.18,0.815,"clipping and fixed-threshold walk.");
      tx.DrawLatex(0.18,0.755,"cfd05: clip-walk penalizes DSB1 corners;");
      tx.DrawLatex(0.18,0.71,"LED: noise penalizes dim LuAG corners.");
      tx.DrawLatex(0.18,0.655,"#Rightarrow diagnostic, not headline, estimators.");
      DrawPadTitle("C  estimator dependence"); }
    c->cd(0); DrawSuperTitle("MIXED module head-to-head: DSB1 vs LuAG:Ce wavelength shifters reading the same showers (GATE-1 map, GATE-6 statistics)",0.019f);
    c->Print("papers/figures/mixed_killshot_bootstrap/mixed_h2h_corrected.png");
    c->Print("papers/figures/mixed_killshot_bootstrap/mixed_h2h_corrected.pdf");
    printf("  wrote papers/figures/mixed_killshot_bootstrap/mixed_h2h_corrected.{png,pdf}\n");
}

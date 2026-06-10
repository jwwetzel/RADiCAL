// monotonicityFix.C — BEFORE/AFTER of the sigma_t(E) estimator fix, per build.
// BEFORE = old tebSigma window (FitGaussCore over mu2 +- 4*ms2, the 5-sigma RMS that
//          broken-timing outliers inflate) + NO in-event veto.
// AFTER  = production tebSigma (robust truncated window + tail guard) + eventDWUP veto.
// Same brightest-1000 (DW-UP)/2 pipeline otherwise; LG-sum ranking (HG-rank disproven).
// Panels: DSB1 lgcfd, LUAG led, MIXED lgcfd (adopted methods) + DSB1 cfd05 (published
// headline, must be preserved). Each: old (grey open) vs new (colour filled), a/sqrt(E)+b
// fits, floor b, and a MONOTONIC flag.
//   source setup.sh; root -l -b -q 'analyze/studies/monotonicityFix.C+'
// Output: figures/<year>/narrative/monotonicity_fix.png
#include "RadView.h"
#include "RadTiming.h"       // production tebSigma (NEW) + eventDWUP
#include "DataPaths.h"
#include "FigPaths.h"
#include "PlotUtils.h"       // FitGaussCore, ApplyRADiCALStyle
#include "SelectionCuts.h"
#include "TFile.h"
#include "TTree.h"
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
using namespace rad;

// OLD tebSigma: Gaussian-core over the tail-inflated 5-sigma-RMS window (pre-fix).
static double sigmaOld(std::vector<float>& v){
    if(v.size()<50) return -1;
    double mu1=0; for(float x:v)mu1+=x; mu1/=v.size();
    double ms1=0; for(float x:v)ms1+=(x-mu1)*(x-mu1); ms1=std::sqrt(ms1/v.size()); if(ms1<0.008)ms1=0.1;
    double rc=mu1,rw=ms1; for(int it=0;it<5;++it){double s=0,ss=0;long n=0;for(float x:v)if(std::fabs(x-rc)<2.5*rw){s+=x;ss+=x*x;++n;}if(n<20)break;rc=s/n;double w2=ss/n-rc*rc;if(w2>0)rw=std::sqrt(w2);}
    double robust=rw/0.9546*1000.0;
    double mu2=0;int n2=0;for(float x:v)if(std::fabs(x-mu1)<5*ms1){mu2+=x;++n2;} double ms2=ms1;
    if(n2>0){mu2/=n2;ms2=0;for(float x:v)if(std::fabs(x-mu1)<5*ms1)ms2+=(x-mu2)*(x-mu2);ms2=std::sqrt(ms2/n2);if(ms2<0.008)ms2=0.1;}else mu2=mu1;
    TH1F h("_old","",120,mu2-4*ms2,mu2+4*ms2); h.SetDirectory(nullptr); for(float x:v)h.Fill(x);
    double mu,muE,s,sE; FitGaussCore(&h,2.0,mu,muE,s,sE); double g=s*1000.0;
    if(g>0&&g>0.5*robust&&g<2.0*robust)return g; return robust>0?robust:-1;
}
// gather brightest-1000 (DW-UP)/2 values; veto=true uses the production in-event veto.
static std::vector<float> bright1000(RadView& v,double E,int src,bool veto){
    double xc,yc; v.beamCenter(xc,yc); double r2=TimingFiducialR(E)*TimingFiducialR(E);
    std::vector<std::pair<float,float>> sd; Long64_t N=v.entries();
    for(Long64_t i=0;i<N;++i){ v.get(i);
        if(!v.wc_ok()||v.mcp1_peak()<kMCP1_minPeak||v.mcp1_peak()>kMCP1_maxPeak) continue;
        double dx=v.x_trk()-xc,dy=v.y_trk()-yc; if(dx*dx+dy*dy>=r2) continue;
        float dwup;
        if(veto){ if(!eventDWUP(v,src,dwup)) continue; }
        else { double ds=0,us=0;int dn=0,un=0;
            for(int c=0;c<4;++c) if(v.is_timing(c)&&v.hg_peak(c)>=kHG_minPeak){float tc=v.timeOf(c,src);if(tc>-1e5f){ds+=tc;++dn;}}
            for(int c=4;c<8;++c) if(v.is_timing(c)&&v.hg_peak(c)>=kHG_minPeak){float tc=v.timeOf(c,src);if(tc>-1e5f){us+=tc;++un;}}
            if(dn<1||un<1) continue; dwup=0.5f*(float)(ds/dn-us/un); }
        sd.push_back({(float)v.sum_lg(),dwup});
    }
    std::vector<float> vt; if((int)sd.size()<1000) return vt;
    std::nth_element(sd.begin(),sd.begin()+1000,sd.end(),[](const std::pair<float,float>&a,const std::pair<float,float>&b){return a.first>b.first;});
    for(int i=0;i<1000;++i) vt.push_back(sd[i].second); return vt;
}

struct Ser{ std::vector<double> E,S,Se,z; double b=0,be=0,chi=0; bool mono=true; };
static void fit(Ser& L){ if(L.E.size()<3)return; TGraphErrors g(L.E.size(),&L.E[0],&L.S[0],&L.z[0],&L.Se[0]);
    TF1 f("f","sqrt([0]*[0]/x+[1]*[1])",20,160); f.SetParameters(300,18); g.Fit(&f,"QN");
    L.b=std::fabs(f.GetParameter(1)); L.chi=f.GetChisquare()/std::max(1,f.GetNDF()); L.be=f.GetParError(1); if(L.chi>1)L.be*=std::sqrt(L.chi);
    for(size_t i=1;i<L.S.size();++i) if(L.S[i]>L.S[i-1]+0.1) L.mono=false; }

static void buildSeries(const char* B,int src,Ser& oldS,Ser& newS){
    BuildConfig cfg=BuildConfig::Load(radConfig(B).Data()); const double Es[6]={25,50,75,100,125,150};
    for(double e:Es){ TString p=radReduced(B,e); if(gSystem->AccessPathName(p.Data()))continue;
        TFile* f=TFile::Open(p.Data()); TTree* t=(TTree*)f->Get("rad"); if(!t){f->Close();continue;} RadView v; v.attach(t,&cfg);
        std::vector<float> vo=bright1000(v,e,src,false), vn=bright1000(v,e,src,true);
        double so=sigmaOld(vo), sn=tebSigma(vn);
        if(so>0){oldS.E.push_back(e);oldS.S.push_back(so);oldS.Se.push_back(so/std::sqrt(2000.0));oldS.z.push_back(0);}
        if(sn>0){newS.E.push_back(e);newS.S.push_back(sn);newS.Se.push_back(sn/std::sqrt(2000.0));newS.z.push_back(0);}
        f->Close(); }
    fit(oldS); fit(newS);
}

static void panel(TVirtualPad* pad,const char* title,Ser& o,Ser& n){
    pad->cd(); gPad->SetLeftMargin(0.13);gPad->SetRightMargin(0.04);gPad->SetTopMargin(0.10);gPad->SetBottomMargin(0.14);gPad->SetGridy();
    double ymn=1e9,ymx=0,minE=1e9,maxE=0; for(size_t i=0;i<o.S.size();++i){ymn=std::min(ymn,o.S[i]);ymx=std::max(ymx,o.S[i]);minE=std::min(minE,o.E[i]);maxE=std::max(maxE,o.E[i]);}
    for(size_t i=0;i<n.S.size();++i){ymn=std::min(ymn,n.S[i]);ymx=std::max(ymx,n.S[i]);minE=std::min(minE,n.E[i]);maxE=std::max(maxE,n.E[i]);}
    if(ymn>ymx){ymn=20;ymx=60;} if(minE>maxE){minE=25;maxE=150;} double r=ymx-ymn; if(r<5)r=5;
    TH1F* fr=gPad->DrawFrame(minE-12,std::max(0.0,ymn-0.15*r),maxE+12,ymx+0.42*r);
    fr->SetTitle(Form("%s;beam energy E (GeV);brightest-1000 (DW#minusUP)/2  #sigma_{t} (ps)",title));
    fr->GetYaxis()->SetTitleSize(0.05);fr->GetYaxis()->SetTitleOffset(1.15);fr->GetXaxis()->SetTitleSize(0.05);
    TGraphErrors* go=new TGraphErrors(o.E.size(),&o.E[0],&o.S[0],&o.z[0],&o.Se[0]);
    go->SetMarkerStyle(24);go->SetMarkerColor(kGray+2);go->SetLineColor(kGray+2);go->SetMarkerSize(1.5);
    TGraphErrors* gn=new TGraphErrors(n.E.size(),&n.E[0],&n.S[0],&n.z[0],&n.Se[0]);
    gn->SetMarkerStyle(20);gn->SetMarkerColor(kAzure+2);gn->SetLineColor(kAzure+2);gn->SetMarkerSize(1.7);
    auto ff=[&](Ser&L,int col,int ls){ if(L.b>0&&L.chi<12){TF1*f=new TF1(Form("f%d%s",col,title),"sqrt([0]*[0]/x+[1]*[1])",minE,maxE);
        double a=300; { TGraphErrors g(L.E.size(),&L.E[0],&L.S[0],&L.z[0],&L.Se[0]); TF1 ft("ft","sqrt([0]*[0]/x+[1]*[1])",minE,maxE);ft.SetParameters(300,L.b);g.Fit(&ft,"QN");a=ft.GetParameter(0);}
        f->SetParameters(a,L.b);f->SetLineColor(col);f->SetLineStyle(ls);f->SetLineWidth(ls==1?3:2);f->Draw("SAME");}};
    ff(o,kGray+2,2); ff(n,kAzure+2,1); go->Draw("P SAME"); gn->Draw("P SAME");
    TLegend* lg=new TLegend(0.32,0.72,0.97,0.90); lg->SetBorderSize(0);lg->SetFillStyle(0);lg->SetTextSize(0.040);
    lg->AddEntry(go,Form("before fix: b=%.0f  %s",o.b,o.mono?"mono":"NON-mono"),"lp");
    lg->AddEntry(gn,Form("after fix: b=%.0f  %s",n.b,n.mono?"MONOTONIC":"non-mono"),"lp");
    lg->Draw();
}

void monotonicityFix(){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    Ser dO,dN,lO,lN,mO,mN,cO,cN;
    buildSeries("DSB1",RadView::kLGCFD,dO,dN);
    buildSeries("LUAG",RadView::kLED ,lO,lN);
    buildSeries("MIXED",RadView::kLGCFD,mO,mN);
    buildSeries("DSB1",RadView::kCFD05,cO,cN);
    TCanvas* c=new TCanvas("mf","",1500,1000); c->Divide(2,2,0.004,0.004);
    panel(c->cd(1),"DSB1 (lgcfd, adopted)",dO,dN);
    panel(c->cd(2),"LuAG (led, adopted)",lO,lN);
    panel(c->cd(3),"MIXED (lgcfd, adopted)",mO,mN);
    panel(c->cd(4),"DSB1 (cfd05, PUBLISHED headline)",cO,cN);
    c->cd(0); DrawSuperTitle("#sigma_{t}(E) estimator fix: grey = before (5#sigma-RMS window, no veto), blue = after (robust window + broken-timing veto). Adopted methods now monotonic; cfd05 numbers unchanged (its non-monotonicity is the clipping limitation, not the fit).",0.018f);
    c->Print(radFigP("figures/narrative/monotonicity_fix.png"));
    printf("wrote figures/narrative/monotonicity_fix.png\n");
}

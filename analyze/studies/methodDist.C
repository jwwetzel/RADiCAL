// ============================================================================
// methodDist.C — show the (DW-UP)/2 distribution + Gaussian-core fit behind EVERY
// timing-method data point, so each sigma_t can be checked against the actual data
// (clean Gaussian? tails? bimodal bright/dim-corner split? a selection problem?).
// Grid per build: columns = energies (25..150), rows = the 3 methods (cfd05/led/lgcfd).
// Each panel: the brightest-1000 (DW-UP)/2 histogram, the FitGaussCore Gaussian
// overlaid, and the sigma_t actually used (Gaussian core, or robust truncated-RMS
// fallback) plus BOTH values (G=.. rob=..) so divergence (= non-Gaussian) is visible.
// This reproduces timingBrightestK + tebSigma exactly, just drawn.
//   source setup.sh
//   root -l -b -q 'analyze/studies/methodDist.C+("MIXED")'   (or DSB1 / LUAG)
// Output: figures/<year>/narrative/method_dist_<BUILD>.png
// ============================================================================
#include "RadView.h"
#include "DataPaths.h"
#include "FigPaths.h"
#include "PlotUtils.h"       // FitGaussCore, ApplyRADiCALStyle
#include "SelectionCuts.h"   // TimingFiducialR, kMCP1_*, kHG_minPeak
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>
using namespace rad;

static const int   SRCS[3]={ RadView::kCFD05, RadView::kLED, RadView::kLGCFD };
static const char* SLAB[3]={ "cfd05", "led", "hg_lgcfd" };

struct Row { float slg; float d[3]; bool ok[3]; };

void methodDist(const char* build="MIXED"){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    BuildConfig cfg=BuildConfig::Load(radConfig(build).Data());
    const double Es[6]={25,50,75,100,125,150};
    TCanvas* c=new TCanvas("md","",1860,940); c->Divide(6,3,0.002,0.006);
    printf("\n=== %s: brightest-1000 (DW-UP)/2 fit check (G=Gaussian-core, rob=truncated-RMS) ===\n",build);
    printf("  E   src        N    sigma   G     rob   used\n");
    for(int ei=0; ei<6; ++ei){ double E=Es[ei];
        TString p=radReduced(build,E);
        std::vector<Row> rows;
        double xc=0,yc=0;
        if(!gSystem->AccessPathName(p.Data())){
            TFile* fp=TFile::Open(p.Data());
            if(fp&&!fp->IsZombie()){ TTree* t=(TTree*)fp->Get("rad");
                if(t){ RadView v; v.attach(t,&cfg); v.beamCenter(xc,yc);
                    double r2=TimingFiducialR(E)*TimingFiducialR(E); Long64_t N=v.entries(); rows.reserve(N);
                    for(Long64_t i=0;i<N;++i){ v.get(i);
                        if(!v.wc_ok()||v.mcp1_peak()<kMCP1_minPeak||v.mcp1_peak()>kMCP1_maxPeak) continue;
                        double dx=v.x_trk()-xc,dy=v.y_trk()-yc; if(dx*dx+dy*dy>=r2) continue;
                        Row r; r.slg=(float)v.sum_lg();
                        for(int s=0;s<3;++s){ double ds=0,us=0;int dn=0,un=0;
                            for(int ch=0;ch<4;++ch) if(v.is_timing(ch)&&v.hg_peak(ch)>=kHG_minPeak){float tc=v.timeOf(ch,SRCS[s]);if(tc>-1e5f){ds+=tc;++dn;}}
                            for(int ch=4;ch<8;++ch) if(v.is_timing(ch)&&v.hg_peak(ch)>=kHG_minPeak){float tc=v.timeOf(ch,SRCS[s]);if(tc>-1e5f){us+=tc;++un;}}
                            r.ok[s]=(dn>=1&&un>=1); r.d[s]= r.ok[s]? 0.5f*(float)(ds/dn-us/un) : 0.f; }
                        rows.push_back(r);
                    } }
                fp->Close(); }
        }
        for(int si=0; si<3; ++si){ int pad=si*6+ei+1; c->cd(pad);
            gPad->SetLeftMargin(0.04); gPad->SetRightMargin(0.02); gPad->SetTopMargin(0.02); gPad->SetBottomMargin(0.10);
            // brightest-1000 of the events valid for THIS source (matches timingBrightestK)
            std::vector<std::pair<float,float>> sd; for(const Row& r:rows) if(r.ok[si]) sd.push_back({r.slg,r.d[si]});
            std::vector<float> vt; const int K=1000;
            if((int)sd.size()>=K){ std::nth_element(sd.begin(),sd.begin()+K,sd.end(),
                [](const std::pair<float,float>&a,const std::pair<float,float>&b){return a.first>b.first;});
                for(int i=0;i<K;++i) vt.push_back(sd[i].second); }
            else for(auto&q:sd) vt.push_back(q.second);
            if(vt.size()<50){ TLatex tx; tx.SetNDC(); tx.SetTextSize(0.12); tx.DrawLatex(0.1,0.5,Form("%s %.0f: n/a",SLAB[si],E)); continue; }
            // tebSigma internals, drawn
            double mu1=0; for(float x:vt)mu1+=x; mu1/=vt.size();
            double ms1=0; for(float x:vt)ms1+=(x-mu1)*(x-mu1); ms1=std::sqrt(ms1/vt.size()); if(ms1<0.008)ms1=0.1;
            double rc=mu1,rw=ms1; for(int it=0;it<5;++it){double s=0,ss=0;long n=0;for(float x:vt)if(std::fabs(x-rc)<2.5*rw){s+=x;ss+=x*x;++n;}if(n<20)break;rc=s/n;double w2=ss/n-rc*rc;if(w2>0)rw=std::sqrt(w2);}
            double robust=rw/0.9546*1000.0;
            double mu2=0;int n2=0;for(float x:vt)if(std::fabs(x-mu1)<5*ms1){mu2+=x;++n2;} double ms2=ms1;
            if(n2>0){mu2/=n2;ms2=0;for(float x:vt)if(std::fabs(x-mu1)<5*ms1)ms2+=(x-mu2)*(x-mu2);ms2=std::sqrt(ms2/n2);if(ms2<0.008)ms2=0.1;}else mu2=mu1;
            TH1F* h=new TH1F(Form("h%d%d",ei,si),"",120,mu2-4*ms2,mu2+4*ms2); h->SetDirectory(nullptr);
            for(float x:vt)h->Fill(x);
            double mu,muE,s,sE; FitGaussCore(h,2.0,mu,muE,s,sE); double gfit=s*1000.0;
            bool usedG=(gfit>0&&gfit>0.5*robust&&gfit<2.0*robust); double used=usedG?gfit:robust;
            h->SetLineColor(kAzure+2); h->SetFillColorAlpha(kAzure+2,0.30); h->SetLineWidth(1);
            h->GetXaxis()->SetLabelSize(0.06); h->GetYaxis()->SetLabelSize(0.0);
            h->GetXaxis()->SetTitle("(DW#minusUP)/2 [ns]"); h->GetXaxis()->SetTitleSize(0.07); h->GetXaxis()->SetTitleOffset(0.95);
            h->Draw("HIST");
            if(s>0){ TF1* g=new TF1(Form("g%d%d",ei,si),"gaus",mu-4*s,mu+4*s);
                g->SetParameters(h->GetMaximum()*0.95,mu,s); g->SetLineColor(kRed+1); g->SetLineWidth(2); g->Draw("SAME"); }
            TLatex tx; tx.SetNDC();
            tx.SetTextSize(0.085); tx.SetTextColor(kBlack); tx.DrawLatex(0.06,0.90,Form("%s  %.0f GeV",SLAB[si],E));
            tx.SetTextSize(0.085); tx.SetTextColor(usedG?(kAzure+3):(kRed+2));
            tx.DrawLatex(0.06,0.80,Form("#sigma=%.1f ps%s",used,usedG?"":"  (robust!)"));
            tx.SetTextSize(0.065); tx.SetTextColor(kGray+2);
            tx.DrawLatex(0.06,0.71,Form("G=%.0f  rob=%.0f  N=%zu",gfit,robust,vt.size()));
            printf("  %3.0f %-9s %5zu  %5.1f  %4.0f  %4.0f   %s\n",E,SLAB[si],vt.size(),used,gfit,robust,usedG?"Gauss":"ROBUST");
        }
    }
    c->cd(0); TLatex tt; tt.SetNDC(); tt.SetTextFont(62); tt.SetTextSize(0.026);
    tt.DrawLatex(0.02,0.975,Form("%s: brightest-1000 (DW#minusUP)/2 distribution + Gaussian-core fit behind every method/energy point (red = fitted Gaussian)",build));
    c->Print(radFigP(Form("figures/narrative/method_dist_%s.png",build)));
    printf("  wrote figures/narrative/method_dist_%s.png\n",build);
}

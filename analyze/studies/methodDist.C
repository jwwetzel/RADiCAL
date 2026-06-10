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

// brightest-1000 (DW-UP)/2 values for one (energy,source); fills the panel vector.
static std::vector<float> gatherVT(const char* build,double E,int src,BuildConfig& cfg){
    std::vector<float> vt; TString p=radReduced(build,E); if(gSystem->AccessPathName(p.Data())) return vt;
    TFile* fp=TFile::Open(p.Data()); if(!fp||fp->IsZombie()){if(fp)fp->Close();return vt;}
    TTree* t=(TTree*)fp->Get("rad"); if(!t){fp->Close();return vt;}
    RadView v; v.attach(t,&cfg); double xc,yc; v.beamCenter(xc,yc); double r2=TimingFiducialR(E)*TimingFiducialR(E);
    std::vector<std::pair<float,float>> sd; Long64_t N=v.entries();
    for(Long64_t i=0;i<N;++i){ v.get(i);
        if(!v.wc_ok()||v.mcp1_peak()<kMCP1_minPeak||v.mcp1_peak()>kMCP1_maxPeak) continue;
        double dx=v.x_trk()-xc,dy=v.y_trk()-yc; if(dx*dx+dy*dy>=r2) continue;
        double ds=0,us=0;int dn=0,un=0;
        for(int ch=0;ch<4;++ch) if(v.is_timing(ch)&&v.hg_peak(ch)>=kHG_minPeak){float tc=v.timeOf(ch,src);if(tc>-1e5f){ds+=tc;++dn;}}
        for(int ch=4;ch<8;++ch) if(v.is_timing(ch)&&v.hg_peak(ch)>=kHG_minPeak){float tc=v.timeOf(ch,src);if(tc>-1e5f){us+=tc;++un;}}
        if(dn>=1&&un>=1) sd.push_back({(float)v.sum_lg(),0.5f*(float)(ds/dn-us/un)});
    }
    const int K=1000;
    if((int)sd.size()>=K){ std::nth_element(sd.begin(),sd.begin()+K,sd.end(),
        [](const std::pair<float,float>&a,const std::pair<float,float>&b){return a.first>b.first;});
        for(int i=0;i<K;++i) vt.push_back(sd[i].second); }
    else for(auto&q:sd) vt.push_back(q.second);
    fp->Close(); return vt;
}
// robust core (rc,rw) of a value vector (matches the production tebSigma window)
static void robustCore(const std::vector<float>& vt,double& rc,double& rw){
    double mu1=0; for(float x:vt)mu1+=x; mu1/=vt.size();
    double ms1=0; for(float x:vt)ms1+=(x-mu1)*(x-mu1); ms1=std::sqrt(ms1/vt.size()); if(ms1<0.008)ms1=0.1;
    rc=mu1; rw=ms1; for(int it=0;it<5;++it){double s=0,ss=0;long n=0;for(float x:vt)if(std::fabs(x-rc)<2.5*rw){s+=x;ss+=x*x;++n;}if(n<20)break;rc=s/n;double w2=ss/n-rc*rc;if(w2>0)rw=std::sqrt(w2);}
}

void methodDist(const char* build="MIXED"){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    gStyle->SetGridStyle(3); gStyle->SetGridColor(kGray); gStyle->SetGridWidth(1);   // faint dotted gridlines
    BuildConfig cfg=BuildConfig::Load(radConfig(build).Data());
    const double EsAll[6]={25,50,75,100,125,150};
    std::vector<double> Es; for(double e:EsAll){ TString p=radReduced(build,e); if(!gSystem->AccessPathName(p.Data())) Es.push_back(e); }
    int nE=(int)Es.size(); if(nE==0){ printf("%s: no reduced data\n",build); return; }

    // PASS 1: gather every panel's brightest-1000 + the COMMON x-range (union of cores)
    std::vector<std::vector<float>> VT(nE*3);
    double gxlo=1e9,gxhi=-1e9;
    for(int ei=0; ei<nE; ++ei) for(int si=0; si<3; ++si){
        std::vector<float> vt=gatherVT(build,Es[ei],SRCS[si],cfg); VT[ei*3+si]=vt;
        if(vt.size()>=50){ double rc,rw; robustCore(vt,rc,rw); double half=4.0*rw; if(half<0.08)half=0.08; if(half>0.35)half=0.35;
            gxlo=std::min(gxlo,rc-half); gxhi=std::max(gxhi,rc+half); }
    }
    if(gxlo>gxhi){ gxlo=-0.40; gxhi=0.05; }

    // PASS 2: draw all panels on the COMMON [gxlo,gxhi] range, with faint gridlines
    TCanvas* c=new TCanvas("md","",310*nE,960);
    TPad* g=GridWithTitle(c,nE,3,Form("%s: brightest-1000 (DW#minusUP)/2 histogram + Gaussian-core fit (red) per method (rows) #times energy (cols)",build),0.003,0.010,0.05,0.020);
    printf("\n=== %s: brightest-1000 (DW-UP)/2 fit check (common x=[%.2f,%.2f] ns) ===\n",build,gxlo,gxhi);
    printf("  E   src        N    sigma   G     rob   used\n");
    for(int ei=0; ei<nE; ++ei){ double E=Es[ei];
        for(int si=0; si<3; ++si){ int pad=si*nE+ei+1; g->cd(pad);
            bool leftcol=(ei==0), botrow=(si==2);
            gPad->SetLeftMargin(leftcol?0.135:0.045); gPad->SetRightMargin(0.025); gPad->SetTopMargin(0.025); gPad->SetBottomMargin(botrow?0.18:0.10);
            gPad->SetGridx(); gPad->SetGridy();
            std::vector<float>& vt=VT[ei*3+si];
            if(vt.size()<50){ gPad->DrawFrame(gxlo,0,gxhi,1); TLatex tx; tx.SetNDC(); tx.SetTextAlign(22); tx.SetTextSize(0.10); tx.SetTextColor(kGray+1);
                tx.DrawLatex(0.5,0.55,Form("%s %.0f GeV",SLAB[si],E)); tx.SetTextSize(0.08); tx.DrawLatex(0.5,0.40,"(low N)"); continue; }
            double rc,rw; robustCore(vt,rc,rw); double robust=rw/0.9546*1000.0;
            // sigma is fit on the TIGHT per-panel window (so it matches the production
            // timingBrightestK), but the histogram is DRAWN on the common range.
            double half=4.0*rw; if(half<0.08)half=0.08; if(half>0.35)half=0.35;
            TH1F hf("hf","",50,rc-half,rc+half); hf.SetDirectory(nullptr); for(float x:vt)hf.Fill(x);
            double mu,muE,s,sE; FitGaussCore(&hf,2.0,mu,muE,s,sE); double gfit=s*1000.0;
            bool usedG=(gfit>0&&gfit>0.5*robust&&gfit<2.0*robust); double used=usedG?gfit:robust;
            TH1F* h=new TH1F(Form("h%d%d",ei,si),"",50,gxlo,gxhi); h->SetDirectory(nullptr);  // display histo, COMMON range
            for(float x:vt)h->Fill(x);
            h->SetMaximum(1.34*h->GetMaximum());                       // top headroom for annotations
            h->SetLineColor(kAzure+2); h->SetFillColorAlpha(kAzure+2,0.30); h->SetLineWidth(1);
            h->GetXaxis()->SetLabelSize(botrow?0.075:0.0); h->GetYaxis()->SetLabelSize(leftcol?0.060:0.0);
            h->GetXaxis()->SetNdivisions(505); h->GetYaxis()->SetNdivisions(505); h->GetYaxis()->SetMaxDigits(3);
            h->GetXaxis()->SetTitle(botrow?"(DW#minusUP)/2 (ns)":""); h->GetXaxis()->SetTitleSize(0.085); h->GetXaxis()->SetTitleOffset(1.05);
            h->GetYaxis()->SetTitle(leftcol?"entries":""); h->GetYaxis()->SetTitleSize(0.075); h->GetYaxis()->SetTitleOffset(0.85);
            h->Draw("HIST");
            if(s>0){ TF1* gg=new TF1(Form("gf%d%d",ei,si),"gaus",mu-4*s,mu+4*s);
                gg->SetParameters(h->GetMaximum()/1.34*0.97,mu,s); gg->SetLineColor(kRed+1); gg->SetLineWidth(2); gg->Draw("SAME"); }
            double lx = leftcol ? 0.175 : 0.075;
            TLatex tx; tx.SetNDC(); tx.SetTextSize(0.080); tx.SetTextColor(kBlack); tx.SetTextAlign(12); tx.DrawLatex(lx,0.93,Form("%s  %.0f GeV",SLAB[si],E));
            tx.SetTextColor(usedG?(kAzure+3):(kRed+2)); tx.DrawLatex(lx,0.84,Form("#sigma=%.1f ps%s",used,usedG?"":" *rob"));
            tx.SetTextSize(0.056); tx.SetTextColor(kGray+2); tx.SetTextAlign(32); tx.DrawLatex(0.945,0.90,Form("G=%.0f",gfit)); tx.DrawLatex(0.945,0.83,Form("rob=%.0f",robust)); tx.DrawLatex(0.945,0.76,Form("N=%zu",vt.size()));
            printf("  %3.0f %-9s %5zu  %5.1f  %4.0f  %4.0f   %s\n",E,SLAB[si],vt.size(),used,gfit,robust,usedG?"Gauss":"ROBUST");
        }
    }
    c->Print(radFigP(Form("figures/narrative/method_dist_%s.png",build)));
    printf("  wrote figures/narrative/method_dist_%s.png\n",build);
}

// ============================================================================
// adaptiveTiming.C — per-CHANNEL adaptive timing source: clipped -> hg_lgcfd
// (recover the saturated edge), un-clipped -> hg_led (robust leading edge).
// All three crossings (cfd05/led/lgcfd) and hg_peak are already in the reduced
// 'rad' tree, so this is a pure analysis-level recombination -- no re-reduction.
//
// Motivation: MIXED has bright DSB1 corners that slam the ~820 mV clip (want edge
// recovery) AND dim LuAG corners that never saturate (lgcfd extrapolates a peak
// that isn't there -> noisy). A single module-wide source mishandles half the
// corners, and WHICH corners clip drifts with energy -> a non-monotonic sigma_t(E)
// with tight error bars (a systematic, not statistics). Choosing the source per
// channel by its own clip state fixes it at the root and is a UNIFIED rule:
//   DSB1  (all clip)  -> all lgcfd  (= current headline)
//   LuAG  (none clip) -> all led    (= current headline)
//   MIXED (mixed)     -> lgcfd on DSB1 corners, led on LuAG corners (adaptive)
//
//   source setup.sh
//   root -l -b -q 'analyze/studies/adaptiveTiming.C+'            // clipThr=790
//   root -l -b -q 'analyze/studies/adaptiveTiming.C+(800)'       // scan the threshold
// Output: figures/<year>/narrative/adaptive_timing.png
// ============================================================================
#include "RadView.h"
#include "RadTiming.h"       // rad::tebSigma
#include "DataPaths.h"
#include "FigPaths.h"
#include "PlotUtils.h"
#include "SelectionCuts.h"   // TimingFiducialR, kMCP1_*, kHG_minPeak
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
#include <cstdio>
using namespace rad;

static int adoptedSrc(const char* b){ std::string s=b; return (s=="LUAG"||s=="TENERGY")?RadView::kLED:RadView::kLGCFD; }

// per-channel time: clipped (hg_peak >= clipThr) -> lgcfd, else -> led
static inline float adaptT(RadView& v, int c, float clipThr){
    return (v.hg_peak(c) >= clipThr) ? v.timeOf(c,RadView::kLGCFD) : v.timeOf(c,RadView::kLED);
}

// brightest-1000 (DW-UP)/2 sigma. mode 0 = single fixed source `src`; mode 1 = adaptive.
static double sigmaT(RadView& v, double E, int mode, int src, float clipThr, int K=1000){
    double xc,yc; v.beamCenter(xc,yc); double r2=TimingFiducialR(E)*TimingFiducialR(E);
    std::vector<std::pair<float,float>> sd; Long64_t N=v.entries();
    for(Long64_t i=0;i<N;++i){ v.get(i);
        if(!v.wc_ok()||v.mcp1_peak()<kMCP1_minPeak||v.mcp1_peak()>kMCP1_maxPeak) continue;
        double dx=v.x_trk()-xc,dy=v.y_trk()-yc; if(dx*dx+dy*dy>=r2) continue;
        double ds=0,us=0;int dn=0,un=0;
        for(int c=0;c<4;++c) if(v.is_timing(c)&&v.hg_peak(c)>=kHG_minPeak){ float tc=mode?adaptT(v,c,clipThr):v.timeOf(c,src); if(tc>-1e5f){ds+=tc;++dn;} }
        for(int c=4;c<8;++c) if(v.is_timing(c)&&v.hg_peak(c)>=kHG_minPeak){ float tc=mode?adaptT(v,c,clipThr):v.timeOf(c,src); if(tc>-1e5f){us+=tc;++un;} }
        if(dn<1||un<1) continue;
        sd.push_back({(float)v.sum_lg(), 0.5f*(float)(ds/dn-us/un)});
    }
    if((int)sd.size()<K) return -1;
    std::nth_element(sd.begin(),sd.begin()+K,sd.end(),[](const std::pair<float,float>&a,const std::pair<float,float>&b){return a.first>b.first;});
    std::vector<float> vt; for(int i=0;i<K;++i) vt.push_back(sd[i].second);
    return tebSigma(vt);
}

struct Ld { std::vector<double> E,S,Se,ze; double a=0,b=0,be=0,chi=0; };
static void fitAB(Ld& L){ if(L.E.size()<3) return;
    TGraphErrors g(L.E.size(),&L.E[0],&L.S[0],&L.ze[0],&L.Se[0]);
    TF1 f("f","sqrt([0]*[0]/x+[1]*[1])",20,160); f.SetParameters(300,18); g.Fit(&f,"QN");
    L.a=std::fabs(f.GetParameter(0)); L.b=std::fabs(f.GetParameter(1));
    L.chi=f.GetChisquare()/std::max(1,f.GetNDF()); L.be=f.GetParError(1); if(L.chi>1) L.be*=std::sqrt(L.chi); }

void adaptiveTiming(double clipThr=790){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    const char* builds[3]={"DSB1","LUAG","MIXED"};
    const double Es[6]={25,50,75,100,125,150};
    TCanvas* c=new TCanvas("at","",1500,540); c->Divide(3,1,0.005,0.005);
    printf("\n========== ADAPTIVE (clipped->lgcfd, else->led; clipThr=%.0f) vs single-source ==========\n",clipThr);
    for(int bi=0; bi<3; ++bi){ const char* B=builds[bi]; int adopt=adoptedSrc(B);
        BuildConfig cfg=BuildConfig::Load(radConfig(B).Data());
        Ld La, Ad;   // adopted single-source, Adaptive
        for(double e:Es){ TString p=radReduced(B,e); if(gSystem->AccessPathName(p.Data())) continue;
            TFile* fp=TFile::Open(p.Data()); if(!fp||fp->IsZombie()){if(fp)fp->Close();continue;}
            TTree* t=(TTree*)fp->Get("rad"); if(!t){fp->Close();continue;} RadView v; v.attach(t,&cfg);
            double sa=sigmaT(v,e,0,adopt,clipThr), sd=sigmaT(v,e,1,0,clipThr);
            if(sa>0){ La.E.push_back(e); La.S.push_back(sa); La.Se.push_back(sa/std::sqrt(2000.0)); La.ze.push_back(0); }
            if(sd>0){ Ad.E.push_back(e); Ad.S.push_back(sd); Ad.Se.push_back(sd/std::sqrt(2000.0)); Ad.ze.push_back(0); }
            fp->Close(); }
        fitAB(La); fitAB(Ad);
        printf("\n--- %s ---\n  E      adopted   adaptive\n",B);
        for(size_t k=0;k<Es+6-Es && k<6;++k){ double e=Es[k]; double a=-1,d=-1;
            for(size_t j=0;j<La.E.size();++j) if(std::fabs(La.E[j]-e)<1) a=La.S[j];
            for(size_t j=0;j<Ad.E.size();++j) if(std::fabs(Ad.E[j]-e)<1) d=Ad.S[j];
            if(a>0||d>0) printf("  %3.0f    %6.1f    %6.1f\n",e,a,d); }
        printf("  adopted(%s):  a=%.0f b=%.1f+-%.1f chi2/ndf=%.1f\n", adopt==RadView::kLED?"led":"lgcfd",La.a,La.b,La.be,La.chi);
        printf("  ADAPTIVE:       a=%.0f b=%.1f+-%.1f chi2/ndf=%.1f\n",Ad.a,Ad.b,Ad.be,Ad.chi);
        // panel
        c->cd(bi+1); gPad->SetLeftMargin(0.13); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.09); gPad->SetBottomMargin(0.14); gPad->SetGridy();
        double ymin=1e9,ymax=0; for(double s:La.S){ymin=std::min(ymin,s);ymax=std::max(ymax,s);} for(double s:Ad.S){ymin=std::min(ymin,s);ymax=std::max(ymax,s);}
        if(ymin>ymax){ymin=20;ymax=60;} double r=ymax-ymin; if(r<5)r=5;
        TH1F* fr=gPad->DrawFrame(0,std::max(0.0,ymin-0.12*r),165,ymax+0.42*r);
        fr->SetTitle(Form("%s;beam energy E (GeV);brightest-1000  (DW#minusUP)/2  #sigma_{t} (ps)",B));
        fr->GetYaxis()->SetTitleSize(0.045); fr->GetYaxis()->SetTitleOffset(1.25); fr->GetXaxis()->SetTitleSize(0.045);
        TGraphErrors* gA=new TGraphErrors(La.E.size(),&La.E[0],&La.S[0],&La.ze[0],&La.Se[0]);
        gA->SetMarkerStyle(24); gA->SetMarkerColor(kGray+2); gA->SetLineColor(kGray+2); gA->SetMarkerSize(1.4);
        TGraphErrors* gD=new TGraphErrors(Ad.E.size(),&Ad.E[0],&Ad.S[0],&Ad.ze[0],&Ad.Se[0]);
        gD->SetMarkerStyle(20); gD->SetMarkerColor(kAzure+2); gD->SetLineColor(kAzure+2); gD->SetMarkerSize(1.6);
        if(La.a>0&&La.chi<8){TF1* f=new TF1(Form("fa%d",bi),"sqrt([0]*[0]/x+[1]*[1])",20,160);f->SetParameters(La.a,La.b);f->SetLineColor(kGray+2);f->SetLineStyle(2);f->SetLineWidth(2);f->Draw("SAME");}
        if(Ad.a>0&&Ad.chi<8){TF1* f=new TF1(Form("fd%d",bi),"sqrt([0]*[0]/x+[1]*[1])",20,160);f->SetParameters(Ad.a,Ad.b);f->SetLineColor(kAzure+2);f->SetLineWidth(3);f->Draw("SAME");}
        gA->Draw("P SAME"); gD->Draw("P SAME");
        TLegend* lg=new TLegend(0.40,0.74,0.96,0.91); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.036);
        lg->AddEntry(gA,Form("single %s: b=%.0f#pm%.0f",adopt==RadView::kLED?"led":"lgcfd",La.b,La.be),"lp");
        lg->AddEntry(gD,Form("adaptive: b=%.0f#pm%.0f",Ad.b,Ad.be),"lp"); lg->Draw();
    }
    c->cd(0); TLatex t; t.SetNDC(); t.SetTextFont(62); t.SetTextSize(0.030);
    t.DrawLatex(0.05,0.96,Form("Per-channel adaptive timing (clipped#rightarrowlgcfd, else#rightarrowled, clip@%.0f mV) vs single module-wide source",clipThr));
    c->Print(radFigP("figures/narrative/adaptive_timing.png"));
    printf("\n  wrote figures/narrative/adaptive_timing.png\n");
}

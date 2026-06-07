// ============================================================================
// hgLgClean.C — what are the off-line populations in HG_peak vs LG_peak, and how
//               do we select only the clean (well-correlated) events?
// ----------------------------------------------------------------------------
// The HG-vs-LG plane shows the main gain line PLUS contamination. We robustly fit
// the main UNCLIPPED line, classify each unclipped event by its residual
// res = HG - (a + b*LG) into CORE / STEEP(high HG/LG) / CLOUD(low HG/LG, high LG),
// and profile each class on discriminators that reveal the cause:
//   fiducial distance, mcp_peak, hg_tot (pulse width), hg_charge/hg_peak (shape),
//   and light concentration frac = hg_peak[c]/sum_c hg_peak.
//   ROOT_INCLUDE_PATH=radcore:Analysis root -l -b -q 'radcore/hgLgClean.C+("DSB1",50,0)'
// ============================================================================
#include "RadView.h"
#include "SelectionCuts.h"
#include "PlotUtils.h"
#include "DataPaths.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TSystem.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>

using namespace rad;

struct Acc{ double s=0,ss=0; long n=0; void add(double x){s+=x;ss+=x*x;++n;}
    double mean(){return n?s/n:0;} double rms(){return n?std::sqrt(std::max(0.,ss/n-mean()*mean())):0;} };

void hgLgClean(const char* build, double E, int ch){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    BuildConfig cfg = BuildConfig::Load(Form("data/2023/configs/%s.json",build));
    TFile* fp=TFile::Open(radReduced(build,E)); TTree* t=(TTree*)fp->Get("rad");
    RadView v; v.attach(t,&cfg);
    double xc,yc; v.beamCenter(xc,yc);
    Long64_t N=v.entries();

    // robust main-line fit on UNCLIPPED events (HG in [40,700]), 2 passes 3-sigma
    double a=0,b=5;
    for(int pass=0; pass<3; ++pass){ double sx=0,sy=0,sxx=0,sxy=0; long n=0; double rms=1e9;
        if(pass>0){ Acc r; for(Long64_t i=0;i<N;++i){ v.get(i); double X=v.lg_peak(ch),Y=v.hg_peak(ch);
            if(Y>40&&Y<700&&X>10) r.add(Y-(a+b*X)); } rms=r.rms(); }
        for(Long64_t i=0;i<N;++i){ v.get(i); double X=v.lg_peak(ch),Y=v.hg_peak(ch);
            if(!(Y>40&&Y<700&&X>10)) continue; if(pass>0 && std::fabs(Y-(a+b*X))>2.5*rms) continue;
            sx+=X;sy+=Y;sxx+=X*X;sxy+=X*Y;++n; }
        if(n>50){ double d=n*sxx-sx*sx; b=(n*sxy-sx*sy)/d; a=(sy-b*sx)/n; } }
    printf("[%s ch%d %.0f GeV] main line: HG = %.1f + %.3f*LG\n", build, ch, E, a, b);

    // residual RMS of the core (for the band)
    Acc cr; for(Long64_t i=0;i<N;++i){ v.get(i); double X=v.lg_peak(ch),Y=v.hg_peak(ch);
        if(Y>40&&Y<700&&X>10 && std::fabs(Y-(a+b*X))<200) cr.add(Y-(a+b*X)); }
    double band = std::max(80.0, 2.5*cr.rms());
    printf("  core residual RMS=%.1f mV -> band=±%.0f mV\n", cr.rms(), band);

    // classify UNCLIPPED (HG<780) events + profile discriminators
    const char* cls[3]={"CORE ","STEEP","CLOUD"};
    Acc fid[3], mcp[3], tot[3], shp[3], frac[3]; long cnt[3]={0};
    TH2F* h2=new TH2F("h2",Form("%s ch%d %.0f GeV;LG peak (mV);HG peak (mV)",build,ch,E),150,0,800,150,0,900);
    for(Long64_t i=0;i<N;++i){ v.get(i);
        double X=v.lg_peak(ch),Y=v.hg_peak(ch); if(X<5) continue;
        h2->Fill(X,Y);
        if(Y>=780||Y<40) continue;                       // skip clipped + tiny
        double res=Y-(a+b*X); int k = (std::fabs(res)<band)?0 : (res>0?1:2);
        // discriminators
        double dx=v.x_trk()-xc, dy=v.y_trk()-yc; double fr2=std::sqrt(dx*dx+dy*dy);
        double sumhg=0; for(int c=0;c<8;++c) sumhg+=std::max(0.f,v.hg_peak(c));
        double conc = sumhg>0 ? Y/sumhg : 0;
        fid[k].add(fr2); mcp[k].add(v.mcp1_peak());
        if(v.ev.hg_tot[ch]>0) tot[k].add(v.ev.hg_tot[ch]);
        if(Y>0&&v.ev.hg_charge[ch]>0) shp[k].add(v.ev.hg_charge[ch]/Y);
        frac[k].add(conc); ++cnt[k];
    }
    long tot3=cnt[0]+cnt[1]+cnt[2];
    printf("\n  class   frac%%   fiducialR   mcp_peak   hg_tot   charge/peak   lightConc\n");
    for(int k=0;k<3;++k) printf("  %s  %5.1f   %8.2f   %8.1f   %6.2f   %9.2f   %7.3f\n",
        cls[k], 100.0*cnt[k]/std::max(1L,tot3), fid[k].mean(), mcp[k].mean(), tot[k].mean(), shp[k].mean(), frac[k].mean());
    printf("  (STEEP = above the line = high HG/LG; CLOUD = below = high LG/low HG)\n");

    TCanvas* c=new TCanvas("c_hl","",900,680); c->SetLogz(); c->SetRightMargin(0.13);
    h2->Draw("COLZ");
    gSystem->mkdir("radcore/figs",kTRUE); c->Print(Form("radcore/figs/hglg_%s_ch%d_%.0f.png",build,ch,E));
    printf("wrote radcore/figs/hglg_%s_ch%d_%.0f.png\n",build,ch,E);
    fp->Close();
}

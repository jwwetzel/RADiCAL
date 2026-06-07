// ============================================================================
// brightCoreScan.C — how low does (DW-UP)/2 sigma_t go for hg_lgcfd vs cfd05 as
// we tighten the transverse fiducial onto the bright core? Mirrors
// RadTiming::timingBestBin exactly (LG-weighted centroid, brightest sum_lg bin,
// tebSigma) but sweeps the fiducial radius. Establishes the genuine best-bin
// floor of the HG/LG-ratio method on full-stat canonical DSB1.
//   source setup.sh; root -l -b -q 'analyze/studies/brightCoreScan.C+("DSB1")'
// ============================================================================
#include "RadView.h"
#include "RadTiming.h"
#include "SelectionCuts.h"
#include "PlotUtils.h"
#include "DataPaths.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include <vector>
#include <algorithm>
#include <cstdio>
using namespace rad;

// best-bin (DW-UP)/2 sigma for a given source + fiducial radius (mirrors timingBestBin)
static double bestBinSigma(RadView& v, int src, double rFid, double& bestE, size_t& nfid){
    double xc,yc; v.beamCenter(xc,yc); double r2=rFid*rFid;
    std::vector<float> slg, tval; Long64_t N=v.entries();
    for(Long64_t i=0;i<N;++i){ v.get(i);
        if(!v.wc_ok()||v.mcp1_peak()<kMCP1_minPeak||v.mcp1_peak()>kMCP1_maxPeak) continue;
        double dx=v.x_trk()-xc, dy=v.y_trk()-yc; if(dx*dx+dy*dy>=r2) continue;
        double ds=0,us=0; int dn=0,un=0;
        for(int c=0;c<4;++c) if(v.is_timing(c)&&v.hg_peak(c)>=kHG_minPeak){ float tc=v.timeOf(c,src); if(tc>-1e5){ds+=tc;++dn;} }
        for(int c=4;c<8;++c) if(v.is_timing(c)&&v.hg_peak(c)>=kHG_minPeak){ float tc=v.timeOf(c,src); if(tc>-1e5){us+=tc;++un;} }
        if(dn<1||un<1) continue;
        slg.push_back(v.sum_lg()); tval.push_back(0.5f*(float)(ds/dn-us/un));
    }
    nfid=slg.size(); if(slg.size()<800) return -1;
    double smin=*std::min_element(slg.begin(),slg.end()), smax=*std::max_element(slg.begin(),slg.end());
    TH1F hS("_bcS","",150,smin,smax); hS.SetDirectory(nullptr); for(float x:slg) hS.Fill(x);
    double mu,sig,mue,sige; FitGaussCore(&hS,2.0,mu,mue,sig,sige); if(sig<=0){mu=hS.GetMean();sig=hS.GetRMS();}
    double binLo=mu-2*sig, binW=4*sig/9.0, best=1e9;
    for(int b=0;b<9;++b){ double blo=binLo+b*binW,bhi=blo+binW; std::vector<float> vt;
        for(size_t i=0;i<slg.size();++i) if(slg[i]>=blo&&slg[i]<bhi) vt.push_back(tval[i]);
        if(vt.size()<400) continue; double s=tebSigma(vt);
        if(s>12.0&&s<best){ best=s; bestE=0.5*(blo+bhi); } }
    return best;
}

void brightCoreScan(const char* build="DSB1"){
    BuildConfig cfg = BuildConfig::Load(Form("data/2023/configs/%s.json",build));
    const double Es[]={100,125,150}; const double Rs[]={3.0,2.5,2.0,1.5,1.0};
    printf("\n%s  (DW-UP)/2 best-bin sigma_t [ps] vs fiducial radius — cfd05 -> lgcfd\n", build);
    printf("  E\\R     3.0        2.5        2.0        1.5        1.0\n");
    for(double E:Es){
        TFile* fp=TFile::Open(radReduced(build,E)); if(!fp||fp->IsZombie()){printf("  %.0f: no file\n",E);continue;}
        TTree* t=(TTree*)fp->Get("rad"); RadView v; v.attach(t,&cfg);
        printf("  %3.0f  ", E);
        for(double R:Rs){ double be; size_t nf; double s=bestBinSigma(v,RadView::kLGCFD,R,be,nf);
            double s5=bestBinSigma(v,RadView::kCFD05,R,be,nf);
            printf("%4.1f->%4.1f  ", s5, s); }
        printf("\n"); fp->Close();
    }
    printf("(each cell: cfd05 -> lgcfd best-bin sigma_t at that radius)\n");
}

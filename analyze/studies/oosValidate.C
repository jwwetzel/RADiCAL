// ============================================================================
// oosValidate.C — held-out (fold-by-run) validation of the lgcfd timing, the
// same discipline the report used for cfd05 (in-sample best-bin OVERFITS; the
// honest number is measured on data the selection never saw). Splits runs into
// two folds; for each energy reports:
//   quantile best-bin: in-sample (pick best bin on fold A, measure on A)
//                  vs   OOS       (that SAME sum_lg cut measured on held-out B)
//   brightest-K:    fold A vs fold B (a FIXED selection -> no argmin bias; the
//                   two should agree, i.e. it is already OOS-honest)
//   source setup.sh; root -l -b -q 'analyze/studies/oosValidate.C+("DSB1")'
// ============================================================================
#include "RadView.h"
#include "RadTiming.h"
#include "SelectionCuts.h"
#include "DataPaths.h"
#include "PlotUtils.h"
#include "TFile.h"
#include "TTree.h"
#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>
using namespace rad;

struct Ev { int run; float slg, d; };

static bool depthOf(RadView& v,double xc,double yc,double r2,float& out){
    if(!v.wc_ok()||v.mcp1_peak()<kMCP1_minPeak||v.mcp1_peak()>kMCP1_maxPeak) return false;
    double dx=v.x_trk()-xc, dy=v.y_trk()-yc; if(dx*dx+dy*dy>=r2) return false;
    double ds=0,us=0;int dn=0,un=0;
    for(int c=0;c<4;++c) if(v.is_timing(c)&&v.hg_peak(c)>=kHG_minPeak){float tc=v.timeOf(c,RadView::kLGCFD);if(tc>-1e5){ds+=tc;++dn;}}
    for(int c=4;c<8;++c) if(v.is_timing(c)&&v.hg_peak(c)>=kHG_minPeak){float tc=v.timeOf(c,RadView::kLGCFD);if(tc>-1e5){us+=tc;++un;}}
    if(dn<1||un<1)return false; out=0.5f*(float)(ds/dn-us/un); return true;
}
static double sigOf(std::vector<float>& v){ return v.size()>=300 ? tebSigma(v) : -1; }

void oosValidate(const char* build="DSB1", int K=1000){
    BuildConfig cfg=BuildConfig::Load(Form("data/2023/configs/%s.json",build));
    const double Es[]={25,50,75,100,125,150};
    printf("\n%s lgcfd (DW-UP)/2 OOS validation (fold-by-run):\n",build);
    printf("  E  | quantile best-bin: in-samp -> OOS | brightest-%d: foldA  foldB\n",K);
    for(double E:Es){
        TFile* fp=TFile::Open(radReduced(build,E)); if(!fp||fp->IsZombie())continue;
        TTree* t=(TTree*)fp->Get("rad"); RadView v; v.attach(t,&cfg);
        double xc,yc; v.beamCenter(xc,yc); double r2=9.0; Long64_t N=v.entries();
        std::vector<Ev> ev; std::set<int> runs;
        for(Long64_t i=0;i<N;++i){ v.get(i); float d; if(depthOf(v,xc,yc,r2,d)){ ev.push_back({v.ev.run,(float)v.sum_lg(),d}); runs.insert(v.ev.run);} }
        // fold assignment: alternate sorted runs into A/B
        std::vector<int> rl(runs.begin(),runs.end()); std::map<int,int> fold;
        for(size_t i=0;i<rl.size();++i) fold[rl[i]]=i&1;
        std::vector<Ev> A,B; for(auto&e:ev){ (fold[e.run]?B:A).push_back(e); }
        auto bestBinCut=[&](std::vector<Ev>& S,double& slo,double& shi)->double{ // quantile best bin -> sum_lg cut + in-samp sigma
            std::sort(S.begin(),S.end(),[](const Ev&a,const Ev&b){return a.slg<b.slg;});
            int NB=9; size_t per=S.size()/NB; double best=1e9; slo=shi=0;
            for(int b=0;b<NB;++b){ size_t lo=(size_t)b*per,hi=(b==NB-1)?S.size():lo+per; if(hi-lo<300)continue;
                std::vector<float> vt; for(size_t k=lo;k<hi;++k) vt.push_back(S[k].d); double s=sigOf(vt);
                if(s>12&&s<best){best=s; slo=S[lo].slg; shi=S[hi-1].slg;} } return best; };
        auto sigInCut=[&](std::vector<Ev>& S,double slo,double shi)->double{ std::vector<float> vt;
            for(auto&e:S) if(e.slg>=slo&&e.slg<=shi) vt.push_back(e.d); return sigOf(vt); };
        auto brightK=[&](std::vector<Ev> S)->double{ if((int)S.size()<K)return -1;
            std::nth_element(S.begin(),S.begin()+K,S.end(),[](const Ev&a,const Ev&b){return a.slg>b.slg;});
            std::vector<float> vt; for(int i=0;i<K;++i) vt.push_back(S[i].d); return sigOf(vt); };
        // quantile: train on A -> cut -> measure on B (OOS), and vice versa; average
        double sloA,shiA,sloB,shiB;
        double inA=bestBinCut(A,sloA,shiA), inB=bestBinCut(B,sloB,shiB);
        double oosAB=sigInCut(B,sloA,shiA), oosBA=sigInCut(A,sloB,shiB);
        double inSamp=0.5*(inA+inB), oos=0.5*(oosAB+oosBA);
        double bkA=brightK(A), bkB=brightK(B);
        printf("  %3.0f |        %5.1f    ->  %5.1f      |       %5.1f  %5.1f\n",E,inSamp,oos,bkA,bkB);
        fp->Close();
    }
    printf("  (quantile OOS = honest; brightest-K foldA~foldB => already unbiased)\n");
}

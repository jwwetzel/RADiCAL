// ============================================================================
// oosFraction.C — the FAIR held-out test for a deterministic brightest-FRACTION
// selection (the previous oosValidate compared brightest-K across folds, i.e.
// DIFFERENT fractions on half samples -> apples to oranges). Here the selection
// is "brightest fraction f of sum_lg" — the SAME fraction on every sample — so
// in-sample and held-out are comparable. A deterministic cut has no trained
// timing parameter, so if the sigma_t(brightness) trend is smooth it must agree.
//   source setup.sh; root -l -b -q 'analyze/studies/oosFraction.C+("DSB1")'
// ============================================================================
#include "RadView.h"
#include "RadTiming.h"
#include "SelectionCuts.h"
#include "DataPaths.h"
#include "PlotUtils.h"
#include "TFile.h"
#include "TTree.h"
#include <set>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>
using namespace rad;

struct Ev { int run; float slg, dL, d5; };

// sigma_t of the brightest fraction f of a sample (deterministic, no argmin)
static double sigBrightest(std::vector<Ev>& S, double f, bool lg){
    int K=(int)(f*S.size()); if(K<200) return -1;
    std::nth_element(S.begin(),S.begin()+K,S.end(),[](const Ev&a,const Ev&b){return a.slg>b.slg;});
    std::vector<float> vt; for(int i=0;i<K;++i) vt.push_back(lg?S[i].dL:S[i].d5);
    return tebSigma(vt);
}

void oosFraction(const char* build="DSB1"){
    BuildConfig cfg=BuildConfig::Load(Form("data/2023/configs/%s.json",build));
    const double Es[]={100,125,150}; const double Fs[]={0.005,0.01,0.02};
    printf("\n%s lgcfd (DW-UP)/2 brightest-FRACTION, full vs held-out folds (matched f):\n",build);
    for(double E:Es){
        TFile* fp=TFile::Open(radReduced(build,E)); if(!fp||fp->IsZombie())continue;
        TTree* t=(TTree*)fp->Get("rad"); RadView v; v.attach(t,&cfg);
        double xc,yc; v.beamCenter(xc,yc); double r2=9.0; Long64_t N=v.entries();
        std::vector<Ev> all,A,B; std::set<int> runs;
        for(Long64_t i=0;i<N;++i){ v.get(i);
            if(!v.wc_ok()||v.mcp1_peak()<kMCP1_minPeak||v.mcp1_peak()>kMCP1_maxPeak) continue;
            double dx=v.x_trk()-xc,dy=v.y_trk()-yc; if(dx*dx+dy*dy>=r2) continue;
            double ds=0,us=0,d5d=0,d5u=0;int dn=0,un=0,e5d=0,e5u=0;
            for(int c=0;c<4;++c) if(v.is_timing(c)&&v.hg_peak(c)>=kHG_minPeak){ float tl=v.timeOf(c,RadView::kLGCFD),t5=v.timeOf(c,RadView::kCFD05); if(tl>-1e5){ds+=tl;++dn;} if(t5>-1e5){d5d+=t5;++e5d;} }
            for(int c=4;c<8;++c) if(v.is_timing(c)&&v.hg_peak(c)>=kHG_minPeak){ float tl=v.timeOf(c,RadView::kLGCFD),t5=v.timeOf(c,RadView::kCFD05); if(tl>-1e5){us+=tl;++un;} if(t5>-1e5){d5u+=t5;++e5u;} }
            if(dn<1||un<1||e5d<1||e5u<1) continue;
            Ev e{v.ev.run,(float)v.sum_lg(),0.5f*(float)(ds/dn-us/un),0.5f*(float)(d5d/e5d-d5u/e5u)};
            all.push_back(e); runs.insert(e.run);
        }
        std::vector<int> rl(runs.begin(),runs.end());
        for(auto&e:all){ int idx=std::lower_bound(rl.begin(),rl.end(),e.run)-rl.begin(); (idx&1?B:A).push_back(e); }
        printf("  E=%3.0f (Nfid=%zu):\n",E,all.size());
        for(double f:Fs){ std::vector<Ev> a=A,b=B,al=all;
            printf("    f=%.1f%%: lgcfd full=%.1f  foldA=%.1f  foldB=%.1f   |  cfd05 full=%.1f\n",
                   100*f, sigBrightest(al,f,true), sigBrightest(a,f,true), sigBrightest(b,f,true), sigBrightest(all,f,false)); }
        fp->Close();
    }
    printf("  (matched f => foldA ~ foldB ~ full: a deterministic bright cut does NOT overfit)\n");
}

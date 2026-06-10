// outlierPeek.C — confirm whether the kurt~900 tail in the brightest-1000 (DW-UP)/2
// is broken-timing events (a channel crossing on the wrong pulse feature) vs physics.
// For (build,E,cornerSet,src): take brightest-1000 by sum_lg, find the events whose
// (DW-UP)/2 is farthest from the median, and print each one's per-timing-channel time
// and hg_peak so the broken channel is visible. Also count events beyond sanity bands.
//   root -l -b -q -e '.L analyze/studies/outlierPeek.C+' -e 'outlierPeek("DSB1",125,0,7)'
#include "RadView.h"
#include "DataPaths.h"
#include "SelectionCuts.h"
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>
using namespace rad;
void outlierPeek(const char* build,double E,int cornerSet=0,int src=7){
    std::vector<int> dn,up; if(cornerSet==1){dn={1,3};up={5,7};}else if(cornerSet==2){dn={0,2};up={4,6};}else{dn={0,1,2,3};up={4,5,6,7};}
    BuildConfig cfg=BuildConfig::Load(radConfig(build).Data());
    TString p=radReduced(build,E); TFile* f=TFile::Open(p.Data()); TTree* t=(TTree*)f->Get("rad");
    RadView v; v.attach(t,&cfg); double xc,yc; v.beamCenter(xc,yc); double r2=TimingFiducialR(E)*TimingFiducialR(E);
    struct Ev{float slg,dwup;Long64_t idx;};
    std::vector<Ev> sd; Long64_t N=v.entries();
    for(Long64_t i=0;i<N;++i){ v.get(i);
        if(!v.wc_ok()||v.mcp1_peak()<kMCP1_minPeak||v.mcp1_peak()>kMCP1_maxPeak) continue;
        double dx=v.x_trk()-xc,dy=v.y_trk()-yc; if(dx*dx+dy*dy>=r2) continue;
        double ds=0,us=0;int nd=0,nu=0;
        for(int c:dn) if(v.is_timing(c)&&v.hg_peak(c)>=kHG_minPeak){float tc=v.timeOf(c,src);if(tc>-1e5f){ds+=tc;++nd;}}
        for(int c:up) if(v.is_timing(c)&&v.hg_peak(c)>=kHG_minPeak){float tc=v.timeOf(c,src);if(tc>-1e5f){us+=tc;++nu;}}
        if(nd<1||nu<1) continue;
        sd.push_back({(float)v.sum_lg(),0.5f*(float)(ds/nd-us/nu),i});
    }
    std::nth_element(sd.begin(),sd.begin()+std::min((size_t)1000,sd.size()),sd.end(),[](const Ev&a,const Ev&b){return a.slg>b.slg;});
    int K=std::min((size_t)1000,sd.size()); std::vector<Ev> br(sd.begin(),sd.begin()+K);
    std::vector<float> dv; for(auto&e:br)dv.push_back(e.dwup); std::sort(dv.begin(),dv.end());
    double med=dv[dv.size()/2];
    long n100=0,n300=0,n500=0; for(auto&e:br){double d=std::fabs(e.dwup-med);if(d>0.1)n100++;if(d>0.3)n300++;if(d>0.5)n500++;}
    printf("\n=== %s E=%.0f corners=%d src=%d : brightest-1000 (DW-UP)/2 outlier scan ===\n",build,E,cornerSet,src);
    printf("median=%.3f ns | beyond |d|>0.1ns:%ld  >0.3ns:%ld  >0.5ns:%ld  (of %d)\n",med,n100,n300,n500,K);
    std::sort(br.begin(),br.end(),[&](const Ev&a,const Ev&b){return std::fabs(a.dwup-med)>std::fabs(b.dwup-med);});
    const char* nm[8]={"NW-D","NE-D","SE-D","SW-D","NW-U","NE-U","SE-U","SW-U"};
    printf("--- 12 most extreme events (dwup far from median): per-channel time(ps)/peak(mV) ---\n");
    for(int k=0;k<12&&k<(int)br.size();++k){ v.get(br[k].idx);
        printf("dwup=%+8.3f ns  sum_lg=%.0f | ",br[k].dwup,br[k].slg);
        for(int c=0;c<8;++c){ float tc=v.timeOf(c,src); printf("%s:%.0f/%.0f ",nm[c],tc>-1e5f?tc*1000:-99999,v.hg_peak(c)); }
        printf("\n"); }
    f->Close();
}

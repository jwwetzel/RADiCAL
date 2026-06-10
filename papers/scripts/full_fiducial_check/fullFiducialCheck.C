// ============================================================================
// fullFiducialCheck.C — verify the FULL-FIDUCIAL companion to the brightest-1000
// headline (claims-law requirement; prior-record audit edit M1, 2026-06-09/10).
// Production post-fix chain, identical to systematicsPostfix.C nominal:
// r<3.0 mm, MCP in [min,750] mV, HG>=20 mV, in-event veto W=2.0 ns, tebSigma.
// K=1000 reproduces the headline as a cross-check; K=all is the companion number.
// NO selection or estimator change — this is a read-out of the locked chain.
//   source setup.sh && root -l -b -q 'papers/scripts/full_fiducial_check/fullFiducialCheck.C+'
// ============================================================================
#include "RadView.h"
#include "RadTiming.h"
#include "DataPaths.h"
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>
using namespace rad;

struct Ev { float t[8], pk[8], r2, slg, mcp; bool istim[8]; };

static std::vector<Ev> gatherE(const char* build,double E,int src,BuildConfig& cfg){
    std::vector<Ev> out;
    TString p=radReduced(build,E); if(gSystem->AccessPathName(p.Data())) return out;
    TFile* f=TFile::Open(p.Data()); if(!f||f->IsZombie()){if(f)f->Close();return out;}
    TTree* t=(TTree*)f->Get("rad"); if(!t){f->Close();return out;}
    RadView v; v.attach(t,&cfg); double xc,yc; v.beamCenter(xc,yc);
    Long64_t N=v.entries(); out.reserve(N/4);
    for(Long64_t i=0;i<N;++i){ v.get(i);
        if(!v.wc_ok()) continue;
        Ev e; e.mcp=v.mcp1_peak();
        if(e.mcp<kMCP1_minPeak||e.mcp>kMCP1_maxPeak) continue;
        double dx=v.x_trk()-xc,dy=v.y_trk()-yc; e.r2=(float)(dx*dx+dy*dy);
        if(e.r2>=3.0*3.0) continue;
        e.slg=(float)v.sum_lg(); bool any=false;
        for(int c=0;c<8;++c){ e.istim[c]=v.is_timing(c); e.pk[c]=v.hg_peak(c);
            float tc=e.istim[c]? v.timeOf(c,src) : -1e9f; e.t[c]=tc; if(e.istim[c]&&tc>-1e5f)any=true; }
        if(any) out.push_back(e);
    }
    f->Close(); return out;
}

// production (DW-UP)/2 per event with the in-event consistency veto (W=2.0 ns)
static bool dwup(Ev& e,float hgThr,float W,float& out){
    float tt[8]; bool okc[8]; int m=0; float a[8];
    for(int c=0;c<8;++c){ okc[c]=e.istim[c]&&e.pk[c]>=hgThr&&e.t[c]>-1e5f; if(okc[c]){tt[c]=e.t[c];a[m++]=e.t[c];} }
    if(m<2) return false;
    std::nth_element(a,a+m/2,a+m); float med=a[m/2];
    double ds=0,us=0;int dn=0,un=0;
    for(int c=0;c<8;++c){ if(!okc[c]||std::fabs(tt[c]-med)>=W) continue;
        if(c<4){ds+=tt[c];++dn;} else {us+=tt[c];++un;} }
    if(dn<1||un<1) return false;
    out=0.5f*(float)(ds/dn-us/un); return true;
}

void fullFiducialCheck(){
    BuildConfig cfg=BuildConfig::Load(radConfig("DSB1").Data());
    const double Es[6]={25,50,75,100,125,150};
    printf("\n===== FULL-FIDUCIAL companion check (DSB1, srCFD, production chain) =====\n");
    printf("  %-6s %8s %12s %14s\n","E","N_fid","sigma_full","sigma_K1000");
    for(double E:Es){
        std::vector<Ev> EV=gatherE("DSB1",E,RadView::kLGCFD,cfg);
        std::vector<std::pair<float,float>> sd; sd.reserve(EV.size());
        for(Ev& e:EV){ float v; if(dwup(e,20.f,2.0f,v)) sd.push_back({e.slg,v}); }
        if(sd.size()<100){ printf("  %3.0f: too few events (%zu)\n",E,sd.size()); continue; }
        std::vector<float> all; all.reserve(sd.size()); for(auto&p:sd) all.push_back(p.second);
        double sFull=tebSigma(all);
        double sK=-1;
        if((int)sd.size()>=1000){
            std::nth_element(sd.begin(),sd.begin()+1000,sd.end(),
                [](const std::pair<float,float>&a,const std::pair<float,float>&b){return a.first>b.first;});
            std::vector<float> vt; vt.reserve(1000); for(int i=0;i<1000;++i)vt.push_back(sd[i].second);
            sK=tebSigma(vt);
        }
        printf("  %3.0f %8zu %9.1f ps %11.1f ps\n",E,sd.size(),sFull,sK);
    }
    printf("  (sigma_K1000 at 150 must reproduce the 25.7 headline; sigma_full at 150 is the\n");
    printf("   claims-law companion number for manuscript edit M1.)\n");
}

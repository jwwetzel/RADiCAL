// ============================================================================
// sigmaT.C — build-agnostic timing resolution on the canonical schema
// ----------------------------------------------------------------------------
// THE template for the refactored analysis layer: one analysis, any build. It
// reads the canonical rad::RadEvent and uses the BuildConfig roles to select
// the TIMING capillaries automatically (energy E-type caps are excluded by
// role, not by hardcoded index), then computes the published (DW-UP)/2 corner
// estimator with the run-folded out-of-sample energy-binned best-bin method.
//
// For DSB1 @150 GeV this reproduces the published headline (~27 ps). Point it at
// any config and it just works:
//   ROOT_INCLUDE_PATH=radcore:Analysis root -l -b -q \
//     'radcore/sigmaT.C+("datasets/2023/configs/DSB1.json", 150)'
// ============================================================================
#include "BuildConfig.h"     // radcore
#include "Schema.h"          // radcore
#include "PlotUtils.h"       // Analysis: FitGaussCore
#include "DataPaths.h"       // Analysis: radReduced
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>

namespace { // file-local

struct Ev { int run; float slg; float t; };

double gcoreSig(std::vector<float>& v) {
    if (v.size() < 150) return -1;
    std::vector<float> s = v; std::sort(s.begin(), s.end()); double md = s[s.size()/2];
    TH1F h("_g","",200,md-1.0,md+1.0); for (float x : v) h.Fill(x);
    double mu,muE,sg,sgE; FitGaussCore(&h,2.0,mu,muE,sg,sgE); return sg>0 ? sg*1000.0 : -1;
}

// run-folded out-of-sample best energy-bin sigma_t [ps]
double oosBest(std::vector<Ev>& ev) {
    if (ev.size() < 2500) return -1;
    std::vector<float> sl; for (auto& e : ev) sl.push_back(e.slg); std::sort(sl.begin(), sl.end());
    double md = sl[sl.size()/2], a=0, a2=0; long n=0;
    for (float v : sl) if (std::fabs(v-md) < 0.5*md) { a+=v; a2+=v*v; ++n; }
    double mE=a/n, sE=std::sqrt(a2/n-mE*mE), lo=mE-2*sE, bw=4*sE/9.0;
    std::vector<int> ur; for (auto& e : ev) ur.push_back(e.run); std::sort(ur.begin(), ur.end());
    ur.erase(std::unique(ur.begin(), ur.end()), ur.end());
    const int kF=5; auto fold=[&](int r){ int idx=(int)(std::lower_bound(ur.begin(),ur.end(),r)-ur.begin()); return idx%kF; };
    const int kTrainMin=(500*(kF-1))/kF; std::vector<float> pool;
    for (int f=0; f<kF; ++f) {
        int bb=-1; double btr=1e9;
        for (int b=0; b<9; ++b) {
            double blo=lo+b*bw, bhi=blo+bw; std::vector<float> vt;
            for (auto& e : ev) if (fold(e.run)!=f && e.slg>=blo && e.slg<bhi) vt.push_back(e.t);
            if ((long)vt.size()<kTrainMin) continue; double s=gcoreSig(vt); if (s>0 && s<btr){btr=s;bb=b;}
        }
        if (bb<0) continue; double blo=lo+bb*bw, bhi=blo+bw;
        for (auto& e : ev) if (fold(e.run)==f && e.slg>=blo && e.slg<bhi) pool.push_back(e.t);
    }
    return pool.size()<150 ? -1 : gcoreSig(pool);
}

} // namespace

void sigmaT(const char* configPath, double energy = 150) {
    rad::BuildConfig cfg = rad::BuildConfig::Load(configPath);
    if (!cfg.valid()) { printf("config load failed: %s\n", cfg.error()); return; }

    // role of each canonical end (timing vs energy), from the per-corner config
    bool timingEnd[8] = {false};
    int  nTimeDown=0, nTimeUp=0;
    for (int i=0; i<cfg.nend; ++i) {
        std::string role = "timing";
        for (auto& c : cfg.caps) if (c.corner == cfg.end[i].corner) role = c.role;
        timingEnd[i] = (role == "timing");
        if (timingEnd[i]) { if (i<4) ++nTimeDown; else ++nTimeUp; }
    }
    printf("build %s @ %.0f GeV: %d timing ends (%d down / %d up); energy caps excluded by role\n",
           cfg.build.c_str(), energy, nTimeDown+nTimeUp, nTimeDown, nTimeUp);

    TFile* fp = TFile::Open(radReduced(cfg.build.c_str(), energy));
    if (!fp || fp->IsZombie()) { printf("no reduced file for %s @ %.0f\n", cfg.build.c_str(), energy); return; }
    TTree* t = (TTree*)fp->Get("rad");
    rad::RadEvent ev; ev.ConnectBranches(t);
    long N = t->GetEntries();

    double xs=0, ys=0; long nw=0;
    for (long i=0; i<N && nw<50000; ++i) { t->GetEntry(i); if (ev.wc_ok && ev.x_trk>-100 && ev.x_trk<100){xs+=ev.x_trk;ys+=ev.y_trk;++nw;} }
    double xc=xs/nw, yc=ys/nw;

    std::vector<Ev> events;
    for (long i=0; i<N; ++i) {
        t->GetEntry(i);
        if (!ev.wc_ok || ev.mcp1_peak<200 || ev.mcp1_peak>750) continue;
        double dx=ev.x_trk-xc, dy=ev.y_trk-yc; if (dx*dx+dy*dy >= 9.0) continue;
        double ds=0, us=0; int dn=0, un=0;
        for (int c=0; c<4; ++c) if (timingEnd[c] && ev.hg_cfd05[c]>-1e5) { ds+=ev.hg_cfd05[c]; ++dn; }
        for (int c=4; c<8; ++c) if (timingEnd[c] && ev.hg_cfd05[c]>-1e5) { us+=ev.hg_cfd05[c]; ++un; }
        if (dn<1 || un<1) continue;
        Ev e; e.run=ev.run; e.slg=ev.sum_lg; e.t=0.5f*(float)(ds/dn-us/un);
        events.push_back(e);
    }
    double s = oosBest(events);
    printf("  events=%zu   out-of-sample best-bin (DW-UP)/2 sigma_t = %.1f ps\n", events.size(), s);
}

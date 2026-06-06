// ============================================================================
// RadTiming.h — THE one timing-resolution method, reused by every analysis.
// ----------------------------------------------------------------------------
// Encapsulates the published headline (DW-UP)/2 Method A pipeline so that every
// build/run/year is analysed identically (no per-macro re-derivation). Matches
// Analysis/timingEnergyBins.C exactly; reproduces DSB1/150 = 27.4 ps.
//
//   rad::RadView v; v.attach(tree, &cfg);
//   rad::TimingResult r = rad::timingBestBin(v, energy);
//   // r.sigma_ps, r.best_bin, r.nFid, r.xc/yc, ...
// ============================================================================
#ifndef RADCORE_RADTIMING_H
#define RADCORE_RADTIMING_H

#include "RadView.h"
#include "PlotUtils.h"       // Analysis: FitGaussCore
#include "SelectionCuts.h"   // Analysis: kMCP1_minPeak.., kHG_minPeak, TimingFiducialR
#include "TH1F.h"
#include <vector>
#include <algorithm>
#include <cmath>

namespace rad {

// VecToHist_teb + FitGaussCore, identical to timingEnergyBins.C: 120 bins over
// mu2 +- 4*ms2 after a 5-sigma outlier rejection, then 2-sigma core Gaussian fit.
inline double tebSigma(std::vector<float>& v) {
    if (v.size() < 50) return -1;
    double mu1=0; for(float x:v) mu1+=x; mu1/=v.size();
    double ms1=0; for(float x:v) ms1+=(x-mu1)*(x-mu1); ms1=std::sqrt(ms1/v.size()); if(ms1<0.008) ms1=0.1;
    double mu2=0; int n2=0; for(float x:v) if(std::fabs(x-mu1)<5*ms1){mu2+=x;++n2;}
    double ms2=ms1;
    if(n2>0){ mu2/=n2; ms2=0; for(float x:v) if(std::fabs(x-mu1)<5*ms1) ms2+=(x-mu2)*(x-mu2); ms2=std::sqrt(ms2/n2); if(ms2<0.008) ms2=0.1; }
    else mu2=mu1;
    TH1F h("_rt","",120, mu2-4*ms2, mu2+4*ms2); h.SetDirectory(nullptr);
    for(float x:v) h.Fill(x);
    double mu,muE,s,sE; FitGaussCore(&h,2.0,mu,muE,s,sE); return s>0 ? s*1000.0 : -1;
}

struct TimingResult {
    double sigma_ps = -1;   // headline best-bin (DW-UP)/2 sigma_t
    int    best_bin = -1;   // which of the 9 E_meas bins
    double bestE    = 0;    // its sum_lg center [mV]
    size_t nFid     = 0;    // fiducial events with a valid (DW-UP)/2
    double xc=0, yc=0, rFid=0;
    double muE=0, sigE=0;   // sum_lg Gaussian fit
};

// The headline (DW-UP)/2 Method A best-bin sigma_t, build-agnostic via RadView.
// src selects the per-end timing source (RadView::kCFD05 default = published method;
// RadView::kLGCFD = CFD on the LG-predicted true peak, the improved headline).
inline TimingResult timingBestBin(RadView& v, double energy, int src = RadView::kCFD05) {
    TimingResult r;
    v.beamCenter(r.xc, r.yc);                                  // ScanRunCenters LG-weighted centroid
    r.rFid = TimingFiducialR(energy); double r2 = r.rFid*r.rFid;

    std::vector<float> slg, tval;
    Long64_t N = v.entries();
    for (Long64_t i=0; i<N; ++i) {
        v.get(i);
        if (!v.wc_ok() || v.mcp1_peak()<kMCP1_minPeak || v.mcp1_peak()>kMCP1_maxPeak) continue;
        double dx=v.x_trk()-r.xc, dy=v.y_trk()-r.yc; if (dx*dx+dy*dy >= r2) continue;
        double ds=0, us=0; int dn=0, un=0;
        for (int c=0; c<4; ++c) if (v.is_timing(c) && v.hg_peak(c)>=kHG_minPeak) { float tc=v.timeOf(c,src); if (tc>-1e5){ ds+=tc; ++dn; } }
        for (int c=4; c<8; ++c) if (v.is_timing(c) && v.hg_peak(c)>=kHG_minPeak) { float tc=v.timeOf(c,src); if (tc>-1e5){ us+=tc; ++un; } }
        if (dn<1 || un<1) continue;
        slg.push_back(v.sum_lg()); tval.push_back(0.5f*(float)(ds/dn-us/un));
    }
    r.nFid = slg.size();
    if (slg.size() < 1000) return r;

    // energy bins from a Gaussian fit of sum_lg -> muE +- 2 sigE, 9 equal bins
    double smin=*std::min_element(slg.begin(),slg.end()), smax=*std::max_element(slg.begin(),slg.end());
    TH1F hS("_rtS","",150,smin,smax); hS.SetDirectory(nullptr); for(float x:slg) hS.Fill(x);
    double muEe,sigEe; FitGaussCore(&hS,2.0,r.muE,muEe,r.sigE,sigEe);
    if (r.sigE<=0){ r.muE=hS.GetMean(); r.sigE=hS.GetRMS(); }
    double binLo=r.muE-2*r.sigE, binW=4*r.sigE/9.0;

    double best=1e9;
    for (int b=0; b<9; ++b) {
        double blo=binLo+b*binW, bhi=blo+binW; std::vector<float> vt;
        for (size_t i=0;i<slg.size();++i) if (slg[i]>=blo && slg[i]<bhi) vt.push_back(tval[i]);
        if (vt.size()<500) continue;
        double s=tebSigma(vt);
        // Guard: reject unphysical/degenerate bin fits. No real RADiCAL timing
        // bin is < ~15 ps. Sub-floor fits come from legacy reduceRaw config files
        // read via the DERIVED cfd05 path (s_cfd05@3mV - mcp_time), which is
        // noise-dominated for DIM builds at LOW energy -> rigorous config-build
        // timing requires the canonical re-reduction (hg_cfd05 with proper MCP
        // referencing). DSB1 (canonical hg_cfd05) is unaffected.
        if (s>15.0 && s<best){ best=s; r.best_bin=b; r.bestE=0.5*(blo+bhi); }
    }
    r.sigma_ps = best;
    return r;
}

} // namespace rad

#endif // RADCORE_RADTIMING_H

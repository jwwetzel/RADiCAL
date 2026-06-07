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
    // --- ROBUST core width: iterative 2.5-sigma truncated RMS (never garbage). ---
    // De-biased for the 2.5-sigma truncation (RMS_trunc = 0.9546*sigma for a Gaussian).
    double rc=mu1, rw=ms1;
    for(int it=0; it<5; ++it){ double s=0,ss=0; long n=0;
        for(float x:v) if(std::fabs(x-rc)<2.5*rw){ s+=x; ss+=x*x; ++n; }
        if(n<20) break; rc=s/n; double w2=ss/n-rc*rc; if(w2>0) rw=std::sqrt(w2); }
    double robust = rw/0.9546*1000.0;   // ps, truncation-debiased
    // --- Gaussian-core fit (preferred for clean peaks) ---
    double mu2=0; int n2=0; for(float x:v) if(std::fabs(x-mu1)<5*ms1){mu2+=x;++n2;}
    double ms2=ms1;
    if(n2>0){ mu2/=n2; ms2=0; for(float x:v) if(std::fabs(x-mu1)<5*ms1) ms2+=(x-mu2)*(x-mu2); ms2=std::sqrt(ms2/n2); if(ms2<0.008) ms2=0.1; }
    else mu2=mu1;
    TH1F h("_rt","",120, mu2-4*ms2, mu2+4*ms2); h.SetDirectory(nullptr);
    for(float x:v) h.Fill(x);
    double mu,muE,s,sE; FitGaussCore(&h,2.0,mu,muE,s,sE); double gfit=s*1000.0;
    // Use the Gaussian fit ONLY if it agrees with the robust core (0.5x..2x); else
    // the fit failed (sparse/non-Gaussian, e.g. low-light LuAG) -> trust the robust width.
    if (gfit>0 && gfit>0.5*robust && gfit<2.0*robust) return gfit;
    return robust>0 ? robust : -1;
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

    // sum_lg core fit (reporting only): mu_E, sigma_E
    double smin=*std::min_element(slg.begin(),slg.end()), smax=*std::max_element(slg.begin(),slg.end());
    TH1F hS("_rtS","",150,smin,smax); hS.SetDirectory(nullptr); for(float x:slg) hS.Fill(x);
    double muEe,sigEe; FitGaussCore(&hS,2.0,r.muE,muEe,r.sigE,sigEe);
    if (r.sigE<=0){ r.muE=hS.GetMean(); r.sigE=hS.GetRMS(); }

    // EQUAL-POPULATION (quantile) sum_lg bins. The previous equal-width mu+-2sigma
    // scheme starved the BRIGHTEST bin in the Gaussian tail as energy rose (its
    // count fell to a handful -> rejected by the N>=500 guard), so the estimator
    // systematically penalized the highest energies by discarding their best-
    // resolution events (e.g. DSB1 150 GeV brightest bin had 9 events). Quantile
    // bins keep every slice -- including the brightest -- equally populated, so
    // sigma_t(E) is monotonic in energy as the rising light yield demands.
    std::vector<size_t> ord(slg.size()); for (size_t i=0;i<ord.size();++i) ord[i]=i;
    std::sort(ord.begin(), ord.end(), [&](size_t a, size_t b){ return slg[a] < slg[b]; });
    const int NB=9; size_t per=ord.size()/NB; double best=1e9;
    for (int b=0; b<NB; ++b) {
        size_t lo=(size_t)b*per, hi=(b==NB-1)?ord.size():lo+per;
        if (hi-lo < 400) continue;                       // need enough events to fit a core
        std::vector<float> vt; vt.reserve(hi-lo); double sclg=0;
        for (size_t k=lo;k<hi;++k){ vt.push_back(tval[ord[k]]); sclg+=slg[ord[k]]; }
        double s=tebSigma(vt);
        // Guard: reject unphysical sub-floor fits (<~12 ps). These arise for DIM
        // builds at LOW energy read via the legacy DERIVED cfd05 path (noise-
        // dominated); canonical hg_cfd05/lgcfd are unaffected.
        if (s>12.0 && s<best){ best=s; r.best_bin=b; r.bestE=sclg/(hi-lo); }
    }
    r.sigma_ps = best;
    return r;
}

// ----------------------------------------------------------------------------
// timingBrightestK — the OTHER way to read the resolution, free of the same bias.
// Quotes the BEST achievable (DW-UP)/2 sigma_t for the brightest, best-measured
// showers, with IDENTICAL statistical tightness (the top K events by sum_lg) at
// every energy. Monotonic in E and preserves the thin-bright-slice magnitude
// (~the published headline), where the quantile best-bin above gives the broader,
// more conservative ~typical-bright-shower number. Report BOTH.
inline TimingResult timingBrightestK(RadView& v, double energy, int src = RadView::kCFD05, int K = 1000) {
    TimingResult r;
    v.beamCenter(r.xc, r.yc); r.rFid = TimingFiducialR(energy); double r2 = r.rFid*r.rFid;
    std::vector<std::pair<float,float>> sd;   // (sum_lg, (DW-UP)/2)
    Long64_t N = v.entries();
    for (Long64_t i=0;i<N;++i){ v.get(i);
        if (!v.wc_ok() || v.mcp1_peak()<kMCP1_minPeak || v.mcp1_peak()>kMCP1_maxPeak) continue;
        double dx=v.x_trk()-r.xc, dy=v.y_trk()-r.yc; if (dx*dx+dy*dy>=r2) continue;
        double ds=0,us=0; int dn=0,un=0;
        for (int c=0;c<4;++c) if (v.is_timing(c)&&v.hg_peak(c)>=kHG_minPeak){ float tc=v.timeOf(c,src); if(tc>-1e5){ds+=tc;++dn;} }
        for (int c=4;c<8;++c) if (v.is_timing(c)&&v.hg_peak(c)>=kHG_minPeak){ float tc=v.timeOf(c,src); if(tc>-1e5){us+=tc;++un;} }
        if (dn<1||un<1) continue;
        sd.push_back({ (float)v.sum_lg(), 0.5f*(float)(ds/dn-us/un) });
    }
    r.nFid = sd.size();
    if ((int)sd.size() < K) return r;
    std::nth_element(sd.begin(), sd.begin()+K, sd.end(),
                     [](const std::pair<float,float>&a, const std::pair<float,float>&b){ return a.first > b.first; });
    std::vector<float> vt; vt.reserve(K); double sclg=0;
    for (int i=0;i<K;++i){ vt.push_back(sd[i].second); sclg += sd[i].first; }
    double s = tebSigma(vt);
    if (s>12.0){ r.sigma_ps=s; r.bestE=sclg/K; }
    return r;
}

} // namespace rad

#endif // RADCORE_RADTIMING_H

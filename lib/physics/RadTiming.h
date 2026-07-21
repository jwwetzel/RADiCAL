// ============================================================================
// RadTiming.h — THE one timing-resolution method, reused by every analysis.
// ----------------------------------------------------------------------------
// PRODUCTION HEADLINE (the paper's numbers): timingBrightestK(v, E, kLGCFD)
//   = brightest-1000 (DW-UP)/2 with the in-event consistency veto (eventDWUP)
//   and the robust core width (tebSigma). kLGCFD is the paper's "srCFD"
//   (saturation-recovered CFD; reduced-tree branch hg_lgcfd). Reproduces the
//   published sigma_t(150 GeV, DSB1) = 25.7 +- 0.6 ps.
//
//   rad::RadView v; v.attach(tree, &cfg);
//   rad::TimingResult r = rad::timingBrightestK(v, energy, RadView::kLGCFD);
//   // r.sigma_ps, r.nFid, r.xc/yc, ...
//
// HISTORICAL: timingBestBin (below) is the RETIRED parent-paper-era Method A
// best-bin pipeline (cfd05 default; reproduced the published 27.4 ps era).
// Kept for the record; NOT the paper headline. See README.md (lib table) and
// ANALYSIS_GUIDE.md. Historical reference macro: analyze/studies/timingEnergyBins.C.
// ============================================================================
#ifndef RADCORE_RADTIMING_H
#define RADCORE_RADTIMING_H

#include "RadView.h"
#include "PlotUtils.h"       // lib/viz: FitGaussCore
#include "SelectionCuts.h"   // lib/physics: kMCP1_minPeak.., kHG_minPeak, TimingFiducialR
#include "TH1F.h"
#include <vector>
#include <algorithm>
#include <cmath>

namespace rad {

// The production width estimator: iterative robust core window, then a 2-sigma
// Gaussian-core fit inside it, with a tail guard. (The historical pre-fix version
// windowed on a 5-sigma RMS, which rare broken-timing outliers inflated — the
// sigma-monotonicity bug documented below and fixed 2026-06.)
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
    // NOTE (stats audit 2026-07): 0.9546 is the SINGLE-pass 2.5-sigma truncation
    // factor; the iterated fixed point converges slightly narrower (~0.938), so the
    // robust FALLBACK is ~-2% biased for a pure Gaussian. The published resolutions
    // use the Gaussian-core fit below (fallback only trips for pathological samples),
    // so no gated number is affected; correcting the constant is a post-submission
    // item requiring a full gate rerun (see CODE_AUDIT_2026-07-21.md).
    // --- Gaussian-core fit (preferred for clean peaks) ---
    // Build the 120-bin core window from the ROBUST truncated centre/width (rc,rw),
    // NOT a 5-sigma RMS. The 5-sigma RMS is inflated by the very broken-timing
    // outliers it is meant to fit around -- a single ~30-sigma |DW-UP|/2 event drives
    // the sample kurtosis to ~900 -- which shifted/scaled the window and drove the
    // core fit anomalously low/high, manufacturing a non-monotonic sigma_t(E).
    // FitGaussCore seeds itself from this window's mean/RMS, so a robust window fixes
    // the seed too. rc,rw (ns) already computed above. (sigma-monotonicity audit, 2026-06.)
    double half = 4.0*rw; if(half < 0.04) half = 0.04;   // ns; floor avoids a degenerate window
    TH1F h("_rt","",120, rc-half, rc+half); h.SetDirectory(nullptr);
    for(float x:v) h.Fill(x);
    double mu,muE,s,sE; FitGaussCore(&h,2.0,mu,muE,s,sE); double gfit=s*1000.0;
    // Tail guard: if the sample is STILL pathologically non-Gaussian (a residual
    // outlier the upstream in-event veto missed), the core fit is untrustworthy ->
    // fall back to the robust width. Clean cores have |skew|<~1, kurt<~15, so this
    // never trips for legitimate data and the published Gaussian-core headline stands.
    double sk=0,ku=0; { double m2=0,m3=0,m4=0; for(float x:v){double d=x-rc;m2+=d*d;m3+=d*d*d;m4+=d*d*d*d;}
        m2/=v.size(); m3/=v.size(); m4/=v.size(); if(m2>0){ sk=m3/std::pow(m2,1.5); ku=m4/(m2*m2)-3.0; } }
    if (std::fabs(sk)>3.0 || ku>20.0) return robust>0 ? robust : -1;
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

// Per-event (DW-UP)/2 with an IN-EVENT broken-timing veto. A capillary end whose
// crossing time disagrees with the event's median end-time by more than
// kTimingChanConsistency_ns is a wrong-feature / near-noise crossing (a channel that
// fired ~30 ns off on a pre-pulse, or a 20-100 mV near-threshold pulse admitted
// because brightness is ranked on LG, not on the HG timing pulse) and is dropped --
// the surviving ends still time the event. This removes the ~30-sigma |DW-UP|/2 tail
// that otherwise inflates the estimator window (-> non-monotonic sigma_t(E)). It is a
// NO-OP for clean events (all ends agree to ~0.5 ns). The event is assumed already
// positioned (v.get(i)) and to have passed the wc/MCP/fiducial cuts.
// Returns false (skip event) if fewer than one down AND one up end survive.
inline bool eventDWUP(RadView& v, int src, float& dwup) {
    float t[8]; bool ok[8]; int nok=0;
    for (int c=0;c<8;++c){ ok[c]=false;
        if (v.is_timing(c) && v.hg_peak(c)>=kHG_minPeak){ float tc=v.timeOf(c,src); if(tc>-1e5){ t[c]=tc; ok[c]=true; ++nok; } } }
    if (nok<2) return false;
    float tmp[8]; int m=0; for(int c=0;c<8;++c) if(ok[c]) tmp[m++]=t[c];
    std::nth_element(tmp, tmp+m/2, tmp+m); float med=tmp[m/2];   // event median end-time
    double ds=0,us=0; int dn=0,un=0;
    for (int c=0;c<4;++c) if(ok[c] && std::fabs(t[c]-med)<kTimingChanConsistency_ns){ ds+=t[c]; ++dn; }
    for (int c=4;c<8;++c) if(ok[c] && std::fabs(t[c]-med)<kTimingChanConsistency_ns){ us+=t[c]; ++un; }
    if (dn<1 || un<1) return false;
    dwup = 0.5f*(float)(ds/dn - us/un);
    return true;
}

// RETIRED (see file header) — historical parent-paper-era Method A best-bin
// pipeline; NOT the paper headline. The production headline is timingBrightestK
// below. src: RadView::kCFD05 (historical default) | RadView::kLGCFD (paper "srCFD").
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
        float dwup; if (!eventDWUP(v, src, dwup)) continue;   // in-event broken-timing veto
        slg.push_back(v.sum_lg()); tval.push_back(dwup);
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
        float dwup; if (!eventDWUP(v, src, dwup)) continue;   // in-event broken-timing veto
        sd.push_back({ (float)v.sum_lg(), dwup });
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

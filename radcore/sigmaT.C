// ============================================================================
// sigmaT.C — build-agnostic timing resolution, EXACT headline method
// ----------------------------------------------------------------------------
// The world-class timing analysis applied uniformly to any build/run via the
// canonical RadView layer. It reproduces Analysis/timingEnergyBins.C's headline
// (DW-UP)/2 Method A result (DSB1/150 -> 27.4 ps) by matching it exactly:
//   * beam center  : ScanRunCenters LG-weighted centroid  (RadView::beamCenter)
//   * fiducial     : r < TimingFiducialR(E)                (SelectionCuts.h)
//   * MCP quality  : kMCP1_minPeak..kMCP1_maxPeak
//   * per channel  : hg_peak >= kHG_minPeak, valid CFD-5%
//   * estimator    : (DW-UP)/2 over TIMING capillaries (energy caps excluded by role)
//   * energy bins  : Gaussian fit of sum_lg -> muE +- 2 sigE, 9 equal bins
//   * result       : single best E_meas bin (min sigma, N >= 500)
//   * sigma fit    : VecToHist_teb (120 bins, mu2 +- 4 ms2, 5-sigma reject) + FitGaussCore
//
//   ROOT_INCLUDE_PATH=radcore:Analysis root -l -b -q \
//     'radcore/sigmaT.C+("datasets/2023/configs/DSB1.json", 150)'
// ============================================================================
#include "RadView.h"         // radcore: format-agnostic view (Schema + BuildConfig)
#include "PlotUtils.h"       // Analysis: FitGaussCore
#include "SelectionCuts.h"   // Analysis: kMCP1_minPeak.., kHG_minPeak, TimingFiducialR
#include "DataPaths.h"       // Analysis: radReduced
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>

namespace { // file-local

// VecToHist_teb + FitGaussCore — identical to timingEnergyBins.C
double tebSigma(std::vector<float>& v) {
    if (v.size() < 50) return -1;
    double mu1=0; for(float x:v) mu1+=x; mu1/=v.size();
    double ms1=0; for(float x:v) ms1+=(x-mu1)*(x-mu1); ms1=std::sqrt(ms1/v.size()); if(ms1<0.008) ms1=0.1;
    double mu2=0; int n2=0; for(float x:v) if(std::fabs(x-mu1)<5*ms1){mu2+=x;++n2;}
    double ms2=ms1;
    if(n2>0){ mu2/=n2; ms2=0; for(float x:v) if(std::fabs(x-mu1)<5*ms1) ms2+=(x-mu2)*(x-mu2); ms2=std::sqrt(ms2/n2); if(ms2<0.008) ms2=0.1; }
    else mu2=mu1;
    TH1F h("_st","",120, mu2-4*ms2, mu2+4*ms2); for(float x:v) h.Fill(x);
    double mu,muE,s,sE; FitGaussCore(&h,2.0,mu,muE,s,sE); return s>0 ? s*1000.0 : -1;
}

} // namespace

void sigmaT(const char* configPath, double energy = 150) {
    rad::BuildConfig cfg = rad::BuildConfig::Load(configPath);
    if (!cfg.valid()) { printf("config load failed: %s\n", cfg.error()); return; }

    TFile* fp = TFile::Open(radReduced(cfg.build.c_str(), energy));
    if (!fp || fp->IsZombie()) { printf("no reduced file for %s @ %.0f\n", cfg.build.c_str(), energy); return; }
    TTree* t = (TTree*)fp->Get("rad");
    rad::RadView v; v.attach(t, &cfg);                 // format-agnostic (canonical or legacy slots)

    int nT=0; for (int i=0;i<cfg.nend;++i) if (v.is_timing(i)) ++nT;
    double xc, yc; v.beamCenter(xc, yc);               // ScanRunCenters LG-weighted centroid
    double rFid = TimingFiducialR(energy), rFid2 = rFid*rFid;

    // collect (sum_lg, (DW-UP)/2) over fiducial events
    std::vector<float> slg, tval;
    Long64_t N = v.entries();
    for (Long64_t i=0; i<N; ++i) {
        v.get(i);
        if (!v.wc_ok() || v.mcp1_peak()<kMCP1_minPeak || v.mcp1_peak()>kMCP1_maxPeak) continue;
        double dx=v.x_trk()-xc, dy=v.y_trk()-yc; if (dx*dx+dy*dy >= rFid2) continue;
        double ds=0, us=0; int dn=0, un=0;
        for (int c=0; c<4; ++c) if (v.is_timing(c) && v.hg_peak(c)>=kHG_minPeak) { float tc=v.cfd05(c); if (tc>-1e5){ ds+=tc; ++dn; } }
        for (int c=4; c<8; ++c) if (v.is_timing(c) && v.hg_peak(c)>=kHG_minPeak) { float tc=v.cfd05(c); if (tc>-1e5){ us+=tc; ++un; } }
        if (dn<1 || un<1) continue;
        slg.push_back(v.sum_lg()); tval.push_back(0.5f*(float)(ds/dn-us/un));
    }
    if (slg.size() < 1000) { printf("low stats (%zu)\n", slg.size()); fp->Close(); return; }

    // energy bins: Gaussian fit of sum_lg -> muE +- 2 sigE, 9 equal bins
    double smin=*std::min_element(slg.begin(),slg.end()), smax=*std::max_element(slg.begin(),slg.end());
    TH1F hS("hS","",150,smin,smax); for(float x:slg) hS.Fill(x);
    double muE,muEe,sigE,sigEe; FitGaussCore(&hS,2.0,muE,muEe,sigE,sigEe);
    if (sigE<=0){ muE=hS.GetMean(); sigE=hS.GetRMS(); }
    double binLo=muE-2*sigE, binW=4*sigE/9.0;

    // single best E_meas bin (min sigma, N >= 500)
    double best=1e9; int bb=-1; double bestE=0;
    for (int b=0; b<9; ++b) {
        double blo=binLo+b*binW, bhi=blo+binW; std::vector<float> vt;
        for (size_t i=0;i<slg.size();++i) if (slg[i]>=blo && slg[i]<bhi) vt.push_back(tval[i]);
        if (vt.size()<500) continue;
        double s=tebSigma(vt); if (s>0 && s<best){ best=s; bb=b; bestE=0.5*(blo+bhi); }
    }
    printf("build %-8s @ %.0f GeV [%-9s]  %d timing ends; centroid (%.2f,%.2f) rFid %.1f mm\n",
           cfg.build.c_str(), energy, v.named?"canonical":"slots", nT, xc, yc, rFid);
    printf("  fiducial=%zu  muE=%.0f sigE=%.0f mV  best bin %d (sum_lg~%.0f)\n", slg.size(), muE, sigE, bb, bestE);
    printf("  (DW-UP)/2 best-bin sigma_t = %.1f ps\n", best);
    fp->Close();
}

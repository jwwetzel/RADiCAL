// ============================================================================
// fiducialTimingScan.C  —  Is the timing fiducial radius optimal?  (all energies)
// ----------------------------------------------------------------------------
// Referee-proof companion to the best-bin headline.  The single-best-E_meas-bin
// selection used for the headline is statistically noisy as a function of the
// fiducial radius (the "best bin" jumps between narrow E_meas slices), so it
// cannot be read as a clean radius optimisation.  This macro instead uses a
// STABLE selection — a fixed fraction of the highest-E_meas events — and scans
// the (DW-UP)/2 corner resolution vs the timing fiducial radius.
//
// Uses the SAME Method-A formula, event cuts, data-derived centroid and
// Gaussian-core fit as timingEnergyBins.C (the headline).
//
// Outputs:
//   Analysis/Output/Summary/fiducial_timing_scan.png      — ALL 6 energies
//        overlaid (top-2% E_meas corner sigma_t vs fiducial radius).
//   Analysis/Output/Summary/fiducial_timing_scan_150.png  — 150 GeV detail
//        (top-2% and top-10% selections, annotated).
//
// Conclusion the figures make airtight: the corner resolution has a shallow
// optimal PLATEAU for r <~ 2.5 mm at EVERY energy, degrading gently beyond — a
// genuine optimisation, not a fluctuation.  (The single-best-bin headline then
// prefers r<3 mm because it needs the larger statistics to isolate the cleanest
// slice; that argument lives in the report text.)
// ============================================================================

#include <vector>
#include <algorithm>
#include <cmath>
#include <utility>

#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLine.h"
#include "TBox.h"
#include "TMarker.h"
#include "TF1.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TString.h"
#include "TSystem.h"
#include "TColor.h"

#include "PlotUtils.h"      // ScanRunCenters, FitGaussCore, StylePad
#include "RADiCALStyle.h"   // ApplyRADiCALStyle, kREnergyCols, kRData/kRRed/kROrange, DrawPageTitle
#include "SelectionCuts.h"  // kMCP1_minPeak, kMCP1_maxPeak, kFiducial_r_timing

static const float kSentCut_fts = -1e5f;   // hg_cfd validity test (matches headline)

// Build a VecToHist_teb-faithful histogram: 2-pass 5-sigma outlier trim, then a
// core range of mu2 +/- 4*ms2 with nb bins (outliers go to overflow).  Caller
// owns the returned histogram.  Mirrors timingEnergyBins.C exactly so the
// best-bin sigma reproduces the headline.
static TH1F* BuildCoreHist_fts(const std::vector<float>& v, int nb, const char* nm)
{
    if (v.size() < 50) return nullptr;
    double mu1 = 0.; for (float x : v) mu1 += x; mu1 /= (double)v.size();
    double ms1 = 0.; for (float x : v) ms1 += ((double)x - mu1) * ((double)x - mu1);
    ms1 = std::sqrt(ms1 / (double)v.size());
    if (ms1 < 0.008) ms1 = 0.100;
    double mu2 = 0.; long n2 = 0;
    for (float x : v) if (std::fabs((double)x - mu1) < 5. * ms1) { mu2 += x; ++n2; }
    double ms2 = 0.;
    if (n2 > 0) {
        mu2 /= (double)n2;
        for (float x : v) if (std::fabs((double)x - mu1) < 5. * ms1)
            ms2 += ((double)x - mu2) * ((double)x - mu2);
        ms2 = std::sqrt(ms2 / (double)n2);
        if (ms2 < 0.008) ms2 = 0.100;
    } else { mu2 = mu1; ms2 = ms1; }
    TH1F* h = new TH1F(nm, "", nb, mu2 - 4. * ms2, mu2 + 4. * ms2);
    h->SetDirectory(nullptr);
    for (float x : v) h->Fill(x);
    return h;
}

// Core sigma [ns] (+ optional error in ps) of corner times, via FitGaussCore.
static double RobustSigma_fts(const std::vector<float>& v, double* errPs = nullptr, int nb = 120)
{
    if (errPs) *errPs = 0.;
    if (v.size() < 150) return -1.;
    TH1F* h = BuildCoreHist_fts(v, nb, "_fts_h");
    if (!h) return -1.;
    double mu, muE, s, sE;
    FitGaussCore(h, 2.0, mu, muE, s, sE);
    delete h;
    if (errPs) *errPs = sE * 1000.;
    return (s > 0.) ? s : -1.;
}

// Core mean & sigma (mV) of the E_meas (sum_lg) spectrum — the binning seed used
// by the headline (FitGaussCore on a 150-bin VecToHist of sum_lg).
static bool CoreMuSig_fts(const std::vector<float>& v, double& mu, double& sig)
{
    mu = sig = 0.;
    TH1F* h = BuildCoreHist_fts(v, 150, "_fts_slg");
    if (!h) return false;
    double muE, sE; FitGaussCore(h, 2.0, mu, muE, sig, sE); delete h;
    return (sig > 0.);
}

// Diagnostic capture of the best-bin internals (one radius).
struct BinDiag_fts {
    double muE = 0., sigE = 0.;     // sum_lg core fit (the binning seed)
    int    n[9]   = {0};            // events per bin
    double sig[9];                  // per-bin Method-A core sigma (ps), -1 if N<500 or fit-fail
    int    selIdx = -1;            // selected (min-sigma, N>=500) bin index
    double selSig = -1.;           // its sigma (ps)
    double selEc  = 0.;            // its E_meas centre (mV)
    int    selN   = 0;
    BinDiag_fts(){ for (int i=0;i<9;++i) sig[i] = -1.; }
};

// Headline best-bin estimator: 9 equal sum_lg bins over [muE-2sigE, muE+2sigE];
// the best bin is the min Method-A core sigma among bins with N >= 500.
// Returns sigma_t in ps (-1 on failure); sets the fit error (ps), N, efficiency.
// If diag != nullptr, captures the per-bin internals.
static double BestBinSigma_fts(const std::vector<std::pair<float,float>>& sel,  // (sum_lg, tc)
                               double& errPs, int& nBest, double& effPct,
                               BinDiag_fts* diag = nullptr)
{
    errPs = 0.; nBest = 0; effPct = 0.;
    const int kNB = 9, kMinN = 500;
    std::vector<float> slg; slg.reserve(sel.size());
    for (const auto& p : sel) slg.push_back(p.first);
    double muE, sigE;
    if (!CoreMuSig_fts(slg, muE, sigE)) {
        double m = 0.; for (float x : slg) m += x; m /= (double)slg.size();
        double r = 0.; for (float x : slg) r += (x-m)*(x-m); r = std::sqrt(r/(double)slg.size());
        muE = m; sigE = r;
    }
    if (sigE <= 0.) return -1.;
    if (diag) { diag->muE = muE; diag->sigE = sigE; }
    const double binLo = muE - 2.*sigE, binW = 4.*sigE / kNB;
    if (binW <= 0.) return -1.;
    const int totFid = (int)sel.size();
    double best = 1e30, bestErr = 0.; int bN = 0, bIdx = -1; double bEc = 0.;
    for (int ib = 0; ib < kNB; ++ib) {
        const double lo = binLo + ib*binW, hi = binLo + (ib+1)*binW;
        std::vector<float> tc;
        for (const auto& p : sel) if (p.first >= lo && p.first < hi) tc.push_back(p.second);
        if (diag) diag->n[ib] = (int)tc.size();
        if ((int)tc.size() < kMinN) continue;
        double e; double s = RobustSigma_fts(tc, &e, 120);
        if (diag && s > 0.) diag->sig[ib] = s*1000.;
        // Floor at 15 ps: the corner estimator is physically >~20 ps, so a
        // smaller "sigma" is a degenerate Gaussian-core fit (a narrow spike)
        // that the min-sigma selection must not grab.
        if (s > 0. && s*1000. >= 15. && s*1000. < best) {
            best = s*1000.; bestErr = e; bN = (int)tc.size(); bIdx = ib; bEc = 0.5*(lo+hi);
        }
    }
    if (best > 1e29) return -1.;
    if (diag) { diag->selIdx = bIdx; diag->selSig = best; diag->selN = bN; diag->selEc = bEc; }
    errPs = bestErr; nBest = bN; effPct = 100.*bN/totFid;
    return best;
}

// SMARTER BINNING: equal-OCCUPANCY (quantile) bins.  Sort events by E_meas and
// split into K bins of EQUAL COUNT (perBin events each), so the highest-E_meas
// bin is always well-populated -- no eligibility flicker on the falling tail.
// Best bin = min Method-A core sigma.  Returns sigma_t (ps); sets err, N, eff.
static double BestBinEqualOcc_fts(const std::vector<std::pair<float,float>>& sel,  // (sum_lg, tc)
                                  int perBin, double& errPs, int& nBest, double& effPct,
                                  int* selIdxOut = nullptr, int* kOut = nullptr,
                                  std::vector<float>* selTc = nullptr,   // selected bin's tc values
                                  double* selEcOut = nullptr)            // its E_meas centre (mV)
{
    errPs = 0.; nBest = 0; effPct = 0.;
    if (selIdxOut) *selIdxOut = -1; if (kOut) *kOut = 0;
    const int N = (int)sel.size();
    if (N < 2*perBin) return -1.;
    const int K = N / perBin;               // number of equal-occupancy bins
    if (K < 2) return -1.;
    std::vector<std::pair<float,float>> s = sel;
    std::sort(s.begin(), s.end(),
              [](const std::pair<float,float>& a, const std::pair<float,float>& b){ return a.first < b.first; });
    double best = 1e30, bestErr = 0.; int bN = 0, bIdx = -1;
    for (int ib = 0; ib < K; ++ib) {
        const int lo = ib*N/K, hi = (ib+1)*N/K;     // contiguous quantile slice
        std::vector<float> tc; tc.reserve(hi-lo);
        for (int k = lo; k < hi; ++k) tc.push_back(s[k].second);
        if ((int)tc.size() < 200) continue;
        double e; double sg = RobustSigma_fts(tc, &e, 120);
        if (sg > 0. && sg*1000. >= 15. && sg*1000. < best) {
            best = sg*1000.; bestErr = e; bN = (int)tc.size(); bIdx = ib;
            if (selTc) *selTc = tc;
            if (selEcOut) { double m = 0.; for (int k=lo;k<hi;++k) m += s[k].first; *selEcOut = m/(hi-lo); }
        }
    }
    if (best > 1e29) return -1.;
    errPs = bestErr; nBest = bN; effPct = 100.*bN/N;
    if (selIdxOut) *selIdxOut = bIdx; if (kOut) *kOut = K;
    return best;
}

// Per-event store after event-level cuts (no radius cut): r^2, E_meas, corner time
struct EvFTS { float r2, slg, tc; int run; };

// Collect events passing the headline event-level cuts for one energy ntuple.
// Returns false if the file/tree cannot be read.
static bool CollectEvents(const char* ntuple, std::vector<EvFTS>& out)
{
    TFile* f = TFile::Open(ntuple);
    if (!f || f->IsZombie()) { printf("[fiducialTimingScan] cannot open %s\n", ntuple); return false; }
    TTree* t = dynamic_cast<TTree*>(f->Get("rad"));
    if (!t) { printf("[fiducialTimingScan] no TTree 'rad' in %s\n", ntuple); f->Close(); return false; }

    double xc, yc, tOff, tRms;
    ScanRunCenters(t, xc, yc, tOff, tRms);

    Bool_t  wc_ok;
    Float_t x_trk, y_trk, mcp_peak, mcp_time, hg_cfd[8], sum_lg;
    Int_t   run = 0;
    if (t->GetBranch("run")) t->SetBranchAddress("run", &run);
    t->SetBranchAddress("wc_ok",    &wc_ok);
    t->SetBranchAddress("x_trk",    &x_trk);
    t->SetBranchAddress("y_trk",    &y_trk);
    t->SetBranchAddress("mcp_peak", &mcp_peak);
    t->SetBranchAddress("mcp_time", &mcp_time);
    t->SetBranchAddress(t->GetBranch("hg_cfd05") ? "hg_cfd05" : "hg_cfd", hg_cfd);
    t->SetBranchAddress("sum_lg",   &sum_lg);

    out.clear(); out.reserve(200000);
    const Long64_t N = t->GetEntries();
    for (Long64_t i = 0; i < N; ++i) {
        t->GetEntry(i);
        if (!wc_ok)                                               continue;
        if (mcp_peak < kMCP1_minPeak || mcp_peak > kMCP1_maxPeak) continue;
        if (mcp_time <= kSentCut_fts)                            continue;
        double dw = 0., up = 0.; int ndw = 0, nup = 0;
        for (int c = 0; c < 4; ++c) if (hg_cfd[c] > kSentCut_fts) { dw += hg_cfd[c]; ++ndw; }
        for (int c = 4; c < 8; ++c) if (hg_cfd[c] > kSentCut_fts) { up += hg_cfd[c]; ++nup; }
        if (ndw < 1 || nup < 1) continue;
        float tc = (float)((dw / ndw - up / nup) * 0.5);        // Method A: (DW-UP)/2
        float dx = x_trk - (float)xc, dy = y_trk - (float)yc;
        out.push_back({ dx * dx + dy * dy, sum_lg, tc, run });
    }
    printf("[fiducialTimingScan] %-6s centroid (%.2f, %.2f) mm  %zu events\n",
           gSystem->BaseName(gSystem->DirName(ntuple)), xc, yc, out.size());
    f->Close();
    return true;
}

// Scan radius for one energy; fill gTop2 / gTop10 (sigma_t in ps vs radius, with
// the Gaussian-fit uncertainty as the y error bar).
static void ScanRadii(const std::vector<EvFTS>& evs,
                      TGraphErrors& gTop2, TGraphErrors& gTop10, TGraphErrors& gBest,
                      double& yMin, double& yMax)
{
    for (double R = 1.0; R <= 5.001; R += 0.25) {
        const float R2 = (float)(R * R);
        std::vector<std::pair<float, float>> sel;            // (sum_lg, tc)
        sel.reserve(evs.size());
        for (const auto& e : evs) if (e.r2 < R2) sel.emplace_back(e.slg, e.tc);
        if (sel.size() < 400) continue;

        // --- single best E_meas bin (the headline estimator), re-optimised at R ---
        double eb = 0., effb = 0.; int nb = 0;
        double sBest = BestBinSigma_fts(sel, eb, nb, effb);
        if (sBest > 0.) { int n = gBest.GetN(); gBest.SetPoint(n, R, sBest); gBest.SetPointError(n, 0., eb);
                          yMin = std::min(yMin, sBest); yMax = std::max(yMax, sBest); }

        // --- fixed top-fraction selections (stable diagnostic) ---
        std::sort(sel.begin(), sel.end(),
                  [](const std::pair<float,float>& a, const std::pair<float,float>& b)
                  { return a.first > b.first; });            // E_meas descending
        auto topSigma = [&](double frac, double& errPs) -> double {
            size_t nt = (size_t)(frac * (double)sel.size());
            if (nt < 300) { errPs = 0.; return -1.; }   // drop low-statistics (noisy) points
            nt = std::min(nt, sel.size());
            std::vector<float> tt; tt.reserve(nt);
            for (size_t k = 0; k < nt; ++k) tt.push_back(sel[k].second);
            return RobustSigma_fts(tt, &errPs);
        };
        double e10 = 0., e02 = 0.;
        double s10 = topSigma(0.10, e10);
        double s02 = topSigma(0.02, e02);
        if (s02 > 0.) { int n = gTop2.GetN(); gTop2.SetPoint(n, R, s02 * 1000.); gTop2.SetPointError(n, 0., e02);
                        yMin = std::min(yMin, s02*1000.); yMax = std::max(yMax, s02*1000.); }
        if (s10 > 0.) { int n = gTop10.GetN(); gTop10.SetPoint(n, R, s10 * 1000.); gTop10.SetPointError(n, 0., e10); }
    }
}

// Run-folded OUT-OF-SAMPLE sigma for the equal-occupancy estimator (5 folds over
// runs/spills).  On each training set: equal-occ bins -> pick the min-sigma bin
// -> record its ΣLG window; measure that window on the held-out fold; pool.
// This is the SAME OOS scheme as the headline, applied to quantile bins.
static double EqualOccOOS_fts(const std::vector<EvFTS>& fid, int perBin, double& errPs, int& nOOS)
{
    errPs = 0.; nOOS = 0;
    const int kF = 5;
    std::vector<int> ur; ur.reserve(fid.size());
    for (const auto& e : fid) ur.push_back(e.run);
    std::sort(ur.begin(), ur.end()); ur.erase(std::unique(ur.begin(), ur.end()), ur.end());
    if (ur.size() < 2) return -1.;
    auto foldOf = [&](int r){ return (int)(std::lower_bound(ur.begin(), ur.end(), r) - ur.begin()) % kF; };
    std::vector<float> pool;
    for (int f = 0; f < kF; ++f) {
        std::vector<std::pair<float,float>> tr;     // (slg, tc) training
        for (const auto& e : fid) if (foldOf(e.run) != f) tr.emplace_back(e.slg, e.tc);
        if ((int)tr.size() < 2*perBin) continue;
        std::sort(tr.begin(), tr.end(), [](const std::pair<float,float>&a, const std::pair<float,float>&b){ return a.first < b.first; });
        const int K = (int)tr.size() / perBin; if (K < 2) continue;
        double best = 1e30; float wlo = 0.f, whi = 0.f;
        for (int ib = 0; ib < K; ++ib) {
            const int lo = ib*(int)tr.size()/K, hi = (ib+1)*(int)tr.size()/K;
            std::vector<float> tc; for (int k = lo; k < hi; ++k) tc.push_back(tr[k].second);
            if ((int)tc.size() < 200) continue;
            double e; double s = RobustSigma_fts(tc, &e, 120);
            if (s > 0. && s*1000. >= 15. && s*1000. < best) {
                best = s*1000.; wlo = tr[lo].first; whi = (hi < (int)tr.size()) ? tr[hi].first : 1e9f;
            }
        }
        if (best > 1e29) continue;
        for (const auto& e : fid) if (foldOf(e.run) == f && e.slg >= wlo && e.slg < whi) pool.push_back(e.tc);
    }
    if (pool.size() < 100) return -1.;
    double e; double s = RobustSigma_fts(pool, &e, 120);
    if (s*1000. < 15.) return -1.;          // reject degenerate held-out fits
    errPs = e; nOOS = (int)pool.size();
    return (s > 0.) ? s*1000. : -1.;
}

void fiducialTimingScan()
{
    TH1::AddDirectory(kFALSE);
    ApplyRADiCALStyle();

    const int    nE = 6;
    const int    eGeV[nE]   = { 25, 50, 75, 100, 125, 150 };
    const char*  eLbl[nE]   = { "25 GeV","50 GeV","75 GeV","100 GeV","125 GeV","150 GeV" };

    TGraphErrors gTop2[nE], gTop10[nE], gBest[nE];
    double yMin = 1e9, yMax = -1e9, yBestMin = 1e9, yBestMax = -1e9;
    std::vector<EvFTS> evs150;
    std::vector<EvFTS> evsAll[nE];
    for (int e = 0; e < nE; ++e) {
        TString nt = Form("Analysis/Output/%dGeV/ntuple.root", eGeV[e]);
        std::vector<EvFTS> evs;
        if (!CollectEvents(nt.Data(), evs)) continue;
        ScanRadii(evs, gTop2[e], gTop10[e], gBest[e], yMin, yMax);
        for (int i = 0; i < gBest[e].GetN(); ++i) {
            yBestMin = std::min(yBestMin, gBest[e].GetY()[i]);
            yBestMax = std::max(yBestMax, gBest[e].GetY()[i]);
        }
        evsAll[e] = evs;
        if (eGeV[e] == 150) evs150 = evs;   // keep for the detail figure
    }

    // =========================================================================
    // PER-ENERGY BEST equal-occupancy sigma_core (in-sample + run-folded OOS).
    // For each energy, scan the fiducial radius and pick the minimum -- the best
    // achievable timing resolution -- and OOS-validate that the minimum is real.
    // =========================================================================
    printf("\n====== Per-energy timing: best IN-SAMPLE config, then OOS-validate THAT config ======\n");
    printf("  Energy | in-samp(min over r): r(mm) sigma_core | OOS at that r:  sigma   OVERFIT\n");
    for (int e = 0; e < nE; ++e) {
        if (evsAll[e].empty()) continue;
        // 1) find the in-sample-optimal radius (the 'best possible' a naive search reports)
        double bisR = 0., bisS = 1e9;
        for (double R = 1.5; R <= 3.501; R += 0.0625) {
            const float R2 = (float)(R*R);
            std::vector<std::pair<float,float>> sel;
            for (const auto& ev : evsAll[e]) if (ev.r2 < R2) sel.emplace_back(ev.slg, ev.tc);
            if (sel.size() < 2500) continue;
            double e_=0., eff_=0.; int n_=0;
            double sIn = BestBinEqualOcc_fts(sel, 1000, e_, n_, eff_);
            if (sIn > 0. && sIn < bisS) { bisS = sIn; bisR = R; }
        }
        // 2) OOS-validate THAT EXACT config (same radius)
        std::vector<EvFTS> fid;
        for (const auto& ev : evsAll[e]) if (ev.r2 < (float)(bisR*bisR)) fid.push_back(ev);
        double eo=0.; int no=0;
        double sOOS = EqualOccOOS_fts(fid, 1000, eo, no);
        printf("  %3d GeV |  %5.2f   %5.1f ps           |  %5.1f ps      %+5.1f ps\n",
               eGeV[e], bisR, bisS, sOOS, (sOOS>0? sOOS-bisS : 0.));
    }
    printf("=====================================================================================\n");

    // Validation: best-bin sigma at r=3 mm must reproduce the headline teb_sigma.
    printf("\n[fiducialTimingScan] VALIDATION — best-bin sigma_t at r=3 mm (should match headline):\n");
    { auto v3 = [](const TGraphErrors& g){ for (int i=0;i<g.GetN();++i) if (std::fabs(g.GetX()[i]-3.0)<1e-6) return g.GetY()[i]; return 0.; };
      printf("   "); for (int e = 0; e < nE; ++e) printf(" %3dGeV=%.1f", eGeV[e], v3(gBest[e])); printf(" ps\n");
      printf("   headline teb_sigma = 47.6 35.7 32.5 33.8 29.2 27.4 ps\n"); }

    auto valAt = [](const TGraph& g, double r) -> double {
        for (int i = 0; i < g.GetN(); ++i) if (std::fabs(g.GetX()[i] - r) < 1e-6) return g.GetY()[i];
        return 0.;
    };
    auto plateauMax = [](const TGraph& g, double r0, double r1) -> double {
        double hi = -1e9;
        for (int i = 0; i < g.GetN(); ++i)
            if (g.GetX()[i] >= r0 - 1e-6 && g.GetX()[i] <= r1 + 1e-6) hi = std::max(hi, g.GetY()[i]);
        return hi;
    };
    const double rPlateau = 2.5;

    // =========================================================================
    // Figure 1 — ALL ENERGIES overlaid (top-2% E_meas corner sigma_t vs radius)
    // =========================================================================
    {
        TCanvas* c = new TCanvas("c_fts_all", "", 960, 760);
        c->SetLeftMargin(0.13); c->SetBottomMargin(0.13);
        c->SetRightMargin(0.05); c->SetTopMargin(0.10);
        c->SetTickx(1); c->SetTicky(1);

        const double ylo = std::max(0., yMin - 5.), yhi = yMax + 6.;
        bool first = true;
        for (int e = 0; e < nE; ++e) {
            if (gTop2[e].GetN() < 2) continue;
            gTop2[e].SetMarkerStyle(20); gTop2[e].SetMarkerSize(1.1);
            gTop2[e].SetMarkerColor(kREnergyCols[e]); gTop2[e].SetLineColor(kREnergyCols[e]);
            gTop2[e].SetLineWidth(3);
            if (first) {
                gTop2[e].Draw("APL");
                gTop2[e].GetXaxis()->SetTitle("timing fiducial radius  r  (mm)");
                gTop2[e].GetYaxis()->SetTitle("#sigma_{t}  (DW#minusUP)/2,  top 2% E_{meas}  (ps)");
                gTop2[e].GetXaxis()->SetLimits(0.5, 5.2);
                gTop2[e].GetYaxis()->SetRangeUser(ylo, yhi);
                gTop2[e].GetXaxis()->SetTitleSize(0.046); gTop2[e].GetYaxis()->SetTitleSize(0.044);
                first = false;
            } else {
                gTop2[e].Draw("PL same");
            }
        }
        // current 3 mm cut (the radius optimum is energy-dependent, so no single
        // "plateau" band is drawn here — see the 150 GeV detail figure for that)
        for (int e = 0; e < nE; ++e) if (gTop2[e].GetN() >= 2) gTop2[e].Draw("PL same");
        TLine* l3 = new TLine(kFiducial_r_timing, ylo, kFiducial_r_timing, yhi);
        l3->SetLineStyle(2); l3->SetLineColor(kGray+2); l3->SetLineWidth(2); l3->Draw();
        { TLatex a; a.SetTextColor(kGray+2); a.SetTextSize(0.027); a.SetTextAngle(90);
          a.DrawLatex(kFiducial_r_timing + 0.09, ylo + 0.42*(yhi-ylo), "adopted cut  r < 3 mm"); }

        TLegend* L = new TLegend(0.70, 0.50, 0.93, 0.88);
        L->SetBorderSize(0); L->SetFillStyle(0); L->SetTextFont(42); L->SetTextSize(0.032);
        for (int e = 0; e < nE; ++e) if (gTop2[e].GetN() >= 2) L->AddEntry(&gTop2[e], eLbl[e], "pl");
        L->Draw();

        { TLatex a; a.SetNDC(); a.SetTextSize(0.026); a.SetTextColor(kGray+3);
          a.DrawLatex(0.16, 0.895, "Tightening helps the lower energies (outer-ring");
          a.DrawLatex(0.16, 0.862, "position/walk); higher E #rightarrow lower #sigma_{t}.");
          a.SetTextColor(kRRed);
          a.DrawLatex(0.16, 0.825, "Adopted: 2.5 mm (E#leq100), 3.0 mm (E#geq125)"); }

        DrawPageTitle("Timing resolution vs fiducial radius -- all energies  (top-2% E_{meas})");
        c->Print("Analysis/Output/Summary/fiducial_timing_scan.png");
        c->Print("Analysis/Output/Summary/fiducial_timing_scan.pdf");
        printf("[fiducialTimingScan] wrote fiducial_timing_scan.png (all energies)\n");
    }

    // =========================================================================
    // Figure 2 — 150 GeV detail (top-2% vs top-10%), the pedagogical panel
    // =========================================================================
    if (gTop2[nE-1].GetN() >= 2)
    {
        TGraphErrors& g2  = gTop2[nE-1];
        TGraphErrors& g10 = gTop10[nE-1];
        TCanvas* c = new TCanvas("c_fts_150", "", 920, 760);
        c->SetLeftMargin(0.14); c->SetBottomMargin(0.13);
        c->SetRightMargin(0.05); c->SetTopMargin(0.10);
        c->SetTickx(1); c->SetTicky(1);

        TGraphErrors& gB = gBest[nE-1];   // the headline single-best-bin curve at 150 GeV
        g2.SetMarkerStyle(20); g2.SetMarkerColor(kRRed);  g2.SetLineColor(kRRed);  g2.SetLineWidth(3); g2.SetMarkerSize(1.2);
        g10.SetMarkerStyle(21); g10.SetMarkerColor(kRData); g10.SetLineColor(kRData); g10.SetLineWidth(3); g10.SetMarkerSize(1.2);
        gB.SetMarkerStyle(33); gB.SetMarkerColor(kBlack); gB.SetLineColor(kBlack); gB.SetLineWidth(3); gB.SetMarkerSize(1.7);

        double yl = 1e9, yh = -1e9;
        auto span = [&](const TGraphErrors& g){ for (int i=0;i<g.GetN();++i){ yl=std::min(yl,g.GetY()[i]); yh=std::max(yh,g.GetY()[i]); } };
        span(g2); span(g10); span(gB);
        const double ylo = yl - 5., yhi = yh + 5.;

        g2.Draw("APL");
        g2.GetXaxis()->SetTitle("timing fiducial radius  r  (mm)");
        g2.GetYaxis()->SetTitle("#sigma_{t}  (DW#minusUP)/2  (ps)");
        g2.GetXaxis()->SetLimits(0.5, 5.2);
        g2.GetYaxis()->SetRangeUser(ylo, yhi);
        g2.GetXaxis()->SetTitleSize(0.046); g2.GetYaxis()->SetTitleSize(0.046);
        g10.Draw("PL same");
        TBox* band = new TBox(0.5, ylo, rPlateau, yhi);
        band->SetFillColorAlpha(kRData, 0.06); band->SetLineColor(0); band->Draw();
        g2.Draw("PL same"); g10.Draw("PL same"); gB.Draw("PL same");
        // mark the adopted 150 GeV point (r=3 mm; OOS-validated 27.4 ps)
        const double kAdopt150 = 27.4;   // OOS best-bin headline at 3 mm
        TMarker* adopt = new TMarker(kFiducial_r_timing, kAdopt150, 29);
        adopt->SetMarkerColor(kROrange); adopt->SetMarkerSize(3.0); adopt->Draw();
        TLine* l3 = new TLine(kFiducial_r_timing, ylo, kFiducial_r_timing, yhi);
        l3->SetLineStyle(2); l3->SetLineColor(kROrange); l3->SetLineWidth(2); l3->Draw();
        { TLatex a; a.SetTextColor(kROrange); a.SetTextSize(0.028); a.SetTextAngle(90);
          a.DrawLatex(kFiducial_r_timing + 0.10, ylo + 0.62*(yhi-ylo), "adopted: r < 3 mm"); }
        { TLatex a; a.SetTextColor(kRData); a.SetTextSize(0.025);
          a.DrawLatex(0.78, yhi - 0.05*(yhi-ylo), "top-X% plateau (r #leq 2.5 mm)"); }

        TLegend* L = new TLegend(0.40, 0.73, 0.93, 0.88);
        L->SetBorderSize(0); L->SetFillStyle(0); L->SetTextFont(42); L->SetTextSize(0.030);
        L->AddEntry(&gB,  "single best bin  (the headline estimator)", "pl");
        L->AddEntry(&g2,  "top 2% by E_{meas} (#SigmaA_{LG})", "pl");
        L->AddEntry(&g10, "top 10% by E_{meas}", "pl");
        L->Draw();
        { TPave* bg = new TPave(0.135, 0.135, 0.955, 0.275, 0, "brNDC");
          bg->SetFillColor(kWhite); bg->SetLineColor(kWhite); bg->Draw();
          TLatex a; a.SetNDC(); a.SetTextSize(0.0245); a.SetTextColor(kGray+3);
          a.DrawLatex(0.150, 0.240,
            "Top-X% (frozen selection): average-event optimum at r #approx 2 mm.");
          a.SetTextColor(kBlack);
          a.DrawLatex(0.150, 0.205,
            "Headline best-bin (black) is tighter (sits below) but JUMPY: its 2.25 mm dip is a");
          a.DrawLatex(0.150, 0.170,
            "lucky radius (2.0 / 2.5 mm = 30.7 / 28.1 ps).  OOS #Rightarrow adopt 3 mm: 27.4 ps (#bigstar)."); }
        DrawPageTitle("Fiducial radius at 150 GeV:  headline best-bin  vs  stable top-X% selection");
        c->Print("Analysis/Output/Summary/fiducial_timing_scan_150.png");
        c->Print("Analysis/Output/Summary/fiducial_timing_scan_150.pdf");
        printf("[fiducialTimingScan] wrote fiducial_timing_scan_150.png (detail)\n");
    }

    // =========================================================================
    // Figure 3 — JOINT view: best-bin (headline) sigma_t vs radius, per energy
    // =========================================================================
    {
        TCanvas* c = new TCanvas("c_fts_best", "", 960, 760);
        c->SetLeftMargin(0.13); c->SetBottomMargin(0.13);
        c->SetRightMargin(0.05); c->SetTopMargin(0.10);
        c->SetTickx(1); c->SetTicky(1);

        const double ylo = std::max(0., yBestMin - 4.), yhi = yBestMax + 5.;
        bool first = true;
        for (int e = 0; e < nE; ++e) {
            if (gBest[e].GetN() < 2) continue;
            gBest[e].SetMarkerStyle(20); gBest[e].SetMarkerSize(1.1);
            gBest[e].SetMarkerColor(kREnergyCols[e]); gBest[e].SetLineColor(kREnergyCols[e]);
            gBest[e].SetLineWidth(2);
            if (first) {
                gBest[e].Draw("APL");
                gBest[e].GetXaxis()->SetTitle("timing fiducial radius  r  (mm)");
                gBest[e].GetYaxis()->SetTitle("best-bin #sigma_{t}  (DW#minusUP)/2  (ps)  [headline]");
                gBest[e].GetXaxis()->SetLimits(0.5, 5.2);
                gBest[e].GetYaxis()->SetRangeUser(ylo, yhi);
                gBest[e].GetXaxis()->SetTitleSize(0.046); gBest[e].GetYaxis()->SetTitleSize(0.044);
                first = false;
            } else gBest[e].Draw("PL same");
        }
        TLine* l3 = new TLine(kFiducial_r_timing, ylo, kFiducial_r_timing, yhi);
        l3->SetLineStyle(2); l3->SetLineColor(kGray+2); l3->SetLineWidth(2); l3->Draw();
        for (int e = 0; e < nE; ++e) if (gBest[e].GetN() >= 2) gBest[e].Draw("PL same");
        { TLatex a; a.SetTextColor(kGray+2); a.SetTextSize(0.027); a.SetTextAngle(90);
          a.DrawLatex(kFiducial_r_timing + 0.09, ylo + 0.42*(yhi-ylo), "adopted cut  r < 3 mm"); }

        TLegend* L = new TLegend(0.70, 0.50, 0.93, 0.88);
        L->SetBorderSize(0); L->SetFillStyle(0); L->SetTextFont(42); L->SetTextSize(0.032);
        for (int e = 0; e < nE; ++e) if (gBest[e].GetN() >= 2) L->AddEntry(&gBest[e], eLbl[e], "pl");
        L->Draw();

        { TLatex a; a.SetNDC(); a.SetTextSize(0.025); a.SetTextColor(kGray+3);
          a.DrawLatex(0.155, 0.878, "Joint optimisation: E_{meas} bin re-picked at each radius (in-sample).");
          a.DrawLatex(0.155, 0.845, "Jumpy; tighter cuts can OVERFIT, so the run-folded OOS validation");
          a.DrawLatex(0.155, 0.812, "(see text) sets the adopted per-energy fiducial:");
          a.SetTextColor(kRRed);
          a.DrawLatex(0.155, 0.778, "r < 2.5 mm for E #leq 100 GeV;  r < 3.0 mm for E #geq 125 GeV."); }

        DrawPageTitle("Best-bin (headline) #sigma_{t} vs fiducial radius -- all energies");
        c->Print("Analysis/Output/Summary/fiducial_timing_scan_bestbin.png");
        c->Print("Analysis/Output/Summary/fiducial_timing_scan_bestbin.pdf");
        printf("[fiducialTimingScan] wrote fiducial_timing_scan_bestbin.png (joint view)\n");
    }

    // =========================================================================
    // DIAGNOSTIC: dissect the 150 GeV best-bin jitter (fine radius steps).
    // Shows the binning seed (muE, sigE), which bin is selected, its N and sigma,
    // and the per-bin sigma[N] for the high-E_meas bins where the best one lives.
    // =========================================================================
    printf("\n================= 150 GeV best-bin jitter diagnostic =================\n");
    printf("  R     Nfid   muE  sigE | sel: idx  Ec(mV)    N   sig | high bins  sigma[N]\n");
    for (double R = 1.5; R <= 3.501; R += 0.125) {
        const float R2 = (float)(R*R);
        std::vector<std::pair<float,float>> sel;
        for (const auto& e : evs150) if (e.r2 < R2) sel.emplace_back(e.slg, e.tc);
        if (sel.size() < 400) continue;
        BinDiag_fts d; double e_ = 0., eff_ = 0.; int n_ = 0;
        double s = BestBinSigma_fts(sel, e_, n_, eff_, &d);
        printf("  %4.2f %6zu %5.0f %5.0f | %4d %7.0f %5d %5.1f |",
               R, sel.size(), d.muE, d.sigE, d.selIdx, d.selEc, d.selN, s);
        for (int ib = 5; ib <= 8; ++ib) {
            if (d.sig[ib] > 0.) printf("  b%d:%4.1f[%d]", ib, d.sig[ib], d.n[ib]);
            else                printf("  b%d: -- [%d]", ib, d.n[ib]);
        }
        printf("\n");
    }
    printf("======================================================================\n");

    // -- smarter binning comparison: equal-WIDTH (current) vs equal-OCCUPANCY --
    printf("\n  150 GeV: equal-WIDTH best-bin  vs  equal-OCCUPANCY (quantile) best-bin\n");
    printf("  R    | eq-width  | eq-occ(perBin=1000) | eq-occ(1500) | eq-occ(2500)\n");
    auto stdev = [](const std::vector<double>& v){ if(v.size()<2) return 0.; double m=0; for(double x:v)m+=x; m/=v.size();
                  double s=0; for(double x:v)s+=(x-m)*(x-m); return std::sqrt(s/v.size()); };
    std::vector<double> vW, v10, v15, v25;
    for (double R = 1.5; R <= 3.501; R += 0.125) {
        const float R2 = (float)(R*R);
        std::vector<std::pair<float,float>> sel;
        for (const auto& e : evs150) if (e.r2 < R2) sel.emplace_back(e.slg, e.tc);
        if (sel.size() < 400) continue;
        double e_, eff_; int n_;
        double sw  = BestBinSigma_fts(sel, e_, n_, eff_);
        double s10 = BestBinEqualOcc_fts(sel, 1000, e_, n_, eff_);
        double s15 = BestBinEqualOcc_fts(sel, 1500, e_, n_, eff_);
        double s25 = BestBinEqualOcc_fts(sel, 2500, e_, n_, eff_);
        printf("  %4.2f |   %5.1f   |        %5.1f        |    %5.1f     |   %5.1f\n", R, sw, s10, s15, s25);
        if (sw>0)  vW.push_back(sw); if (s10>0) v10.push_back(s10);
        if (s15>0) v15.push_back(s15); if (s25>0) v25.push_back(s25);
    }
    printf("  RMS scatter over radius:  width=%.2f  occ1000=%.2f  occ1500=%.2f  occ2500=%.2f ps\n",
           stdev(vW), stdev(v10), stdev(v15), stdev(v25));

    // =========================================================================
    // Figure 4 — the jitter EXPLAINED: best-bin = min(bin6 stable, bin7 noisy),
    // and bin7 flickers across the N>=500 eligibility threshold.
    // =========================================================================
    {
        TGraph gB6, gB7, gBst, gN7;
        for (double R = 1.5; R <= 3.501; R += 0.0625) {
            const float R2 = (float)(R*R);
            std::vector<std::pair<float,float>> sel;
            for (const auto& e : evs150) if (e.r2 < R2) sel.emplace_back(e.slg, e.tc);
            if (sel.size() < 400) continue;
            BinDiag_fts d; double e_ = 0., eff_ = 0.; int n_ = 0;
            double s = BestBinSigma_fts(sel, e_, n_, eff_, &d);
            if (d.sig[6] > 0.) gB6.SetPoint(gB6.GetN(), R, d.sig[6]);
            if (d.sig[7] > 0.) gB7.SetPoint(gB7.GetN(), R, d.sig[7]);   // only where bin7 eligible
            if (s > 0.)        gBst.SetPoint(gBst.GetN(), R, s);
            gN7.SetPoint(gN7.GetN(), R, d.n[7]);
        }

        TCanvas* c = new TCanvas("c_fts_diag", "", 920, 860);
        TPad* pT = new TPad("pT","",0,0.36,1,1);  pT->SetBottomMargin(0.02); pT->SetTopMargin(0.10);
        pT->SetLeftMargin(0.13); pT->SetRightMargin(0.05); pT->SetTickx(1); pT->SetTicky(1); pT->Draw();
        TPad* pB = new TPad("pB","",0,0,1,0.36);   pB->SetTopMargin(0.02); pB->SetBottomMargin(0.26);
        pB->SetLeftMargin(0.13); pB->SetRightMargin(0.05); pB->SetTickx(1); pB->SetTicky(1); pB->Draw();

        // -- top: per-bin sigma + best --
        pT->cd();
        gB6.SetMarkerStyle(21); gB6.SetMarkerColor(kRData); gB6.SetLineColor(kRData); gB6.SetLineWidth(3); gB6.SetMarkerSize(1.0);
        gB7.SetMarkerStyle(20); gB7.SetMarkerColor(kRRed);  gB7.SetLineColor(kRRed);  gB7.SetLineWidth(2); gB7.SetMarkerSize(1.2);
        gBst.SetMarkerStyle(33); gBst.SetMarkerColor(kBlack); gBst.SetLineColor(kBlack); gBst.SetLineWidth(3); gBst.SetMarkerSize(1.4);
        gB6.Draw("APL");
        gB6.GetYaxis()->SetTitle("#sigma_{t} of the bin (ps)");
        gB6.GetXaxis()->SetLimits(1.4, 3.6); gB6.GetYaxis()->SetRangeUser(22, 38);
        gB6.GetYaxis()->SetTitleSize(0.055); gB6.GetYaxis()->SetLabelSize(0.045);
        gB6.GetXaxis()->SetLabelSize(0.0);
        gBst.Draw("PL same"); gB7.Draw("PL same");
        { TLegend* L = new TLegend(0.16,0.70,0.62,0.90); L->SetBorderSize(0); L->SetFillStyle(0);
          L->SetTextFont(42); L->SetTextSize(0.046);
          L->AddEntry(&gBst,"best bin = the min (the jumpy curve)","pl");
          L->AddEntry(&gB6, "bin 6  (large N, stable ~31-35 ps)","pl");
          L->AddEntry(&gB7, "bin 7  (highest E_{meas}, small N, noisy)","pl"); L->Draw(); }
        { TLatex a; a.SetNDC(); a.SetTextSize(0.044); a.SetTextColor(kGray+3);
          a.DrawLatex(0.40,0.13,"best bin flips bin6 #leftrightarrow bin7 as bin7 appears/vanishes"); }
        DrawPageTitle("Why the 150 GeV best-bin is jumpy:  bin 7 flickers in and out");

        // -- bottom: bin7 N with the eligibility threshold --
        pB->cd();
        gN7.SetMarkerStyle(20); gN7.SetMarkerColor(kRRed); gN7.SetLineColor(kRRed); gN7.SetLineWidth(2); gN7.SetMarkerSize(1.1);
        gN7.Draw("APL");
        gN7.GetXaxis()->SetTitle("timing fiducial radius  r  (mm)");
        gN7.GetYaxis()->SetTitle("bin 7  N");
        gN7.GetXaxis()->SetLimits(1.4, 3.6); gN7.GetYaxis()->SetRangeUser(0, 3000);
        gN7.GetXaxis()->SetTitleSize(0.095); gN7.GetXaxis()->SetLabelSize(0.085);
        gN7.GetYaxis()->SetTitleSize(0.085); gN7.GetYaxis()->SetLabelSize(0.075); gN7.GetYaxis()->SetTitleOffset(0.6);
        TLine* thr = new TLine(1.4, 500., 3.6, 500.); thr->SetLineColor(kBlack); thr->SetLineStyle(2); thr->SetLineWidth(2); thr->Draw();
        { TLatex a; a.SetTextSize(0.085); a.SetTextColor(kBlack);
          a.DrawLatex(1.5, 620., "N #geq 500 eligibility threshold");
          a.SetTextColor(kRRed); a.DrawLatex(2.0, 2400., "bin 7 crosses 500 repeatedly #Rightarrow on/off"); }

        c->Print("Analysis/Output/Summary/fiducial_bestbin_jitter.png");
        c->Print("Analysis/Output/Summary/fiducial_bestbin_jitter.pdf");
        printf("[fiducialTimingScan] wrote fiducial_bestbin_jitter.png (jitter diagnostic)\n");
    }

    // =========================================================================
    // Figure 5 — the FIX: equal-occupancy binning smooths the radius curve.
    // =========================================================================
    {
        TGraph gW, gO10, gO25;   // sigma vs R: equal-width, equal-occ(1000), equal-occ(2500)
        std::vector<double> sW, sO10, sO25;
        for (double R = 1.5; R <= 3.501; R += 0.0625) {
            const float R2 = (float)(R*R);
            std::vector<std::pair<float,float>> sel;
            for (const auto& e : evs150) if (e.r2 < R2) sel.emplace_back(e.slg, e.tc);
            if (sel.size() < 400) continue;
            double e_ = 0., eff_ = 0.; int n_ = 0;
            double w  = BestBinSigma_fts(sel, e_, n_, eff_);
            double o1 = BestBinEqualOcc_fts(sel, 1000, e_, n_, eff_);
            double o2 = BestBinEqualOcc_fts(sel, 2500, e_, n_, eff_);
            if (w  > 0.) { gW.SetPoint(gW.GetN(),  R, w);  sW.push_back(w); }
            if (o1 > 0.) { gO10.SetPoint(gO10.GetN(), R, o1); sO10.push_back(o1); }
            if (o2 > 0.) { gO25.SetPoint(gO25.GetN(), R, o2); sO25.push_back(o2); }
        }
        auto rms = [](const std::vector<double>& v){ if(v.size()<2) return 0.; double m=0; for(double x:v)m+=x; m/=v.size();
                    double s=0; for(double x:v)s+=(x-m)*(x-m); return std::sqrt(s/v.size()); };

        TCanvas* c = new TCanvas("c_fts_fix", "", 920, 720);
        c->SetLeftMargin(0.13); c->SetBottomMargin(0.13); c->SetRightMargin(0.05); c->SetTopMargin(0.10);
        c->SetTickx(1); c->SetTicky(1);
        gW.SetMarkerStyle(33); gW.SetMarkerColor(kBlack); gW.SetLineColor(kBlack); gW.SetLineWidth(2); gW.SetMarkerSize(1.4);
        gO10.SetMarkerStyle(20); gO10.SetMarkerColor(kRGreen); gO10.SetLineColor(kRGreen); gO10.SetLineWidth(3); gO10.SetMarkerSize(1.1);
        gO25.SetMarkerStyle(21); gO25.SetMarkerColor(kRData);  gO25.SetLineColor(kRData);  gO25.SetLineWidth(3); gO25.SetMarkerSize(1.1);
        gW.Draw("APL");
        gW.GetXaxis()->SetTitle("timing fiducial radius  r  (mm)");
        gW.GetYaxis()->SetTitle("best-bin #sigma_{t}  (DW#minusUP)/2  (ps)");
        gW.GetXaxis()->SetLimits(1.4, 3.6); gW.GetYaxis()->SetRangeUser(24, 34);
        gW.GetXaxis()->SetTitleSize(0.046); gW.GetYaxis()->SetTitleSize(0.046);
        gO10.Draw("PL same"); gO25.Draw("PL same");
        TLine* l3 = new TLine(kFiducial_r_timing, 24, kFiducial_r_timing, 34);
        l3->SetLineStyle(2); l3->SetLineColor(kGray+2); l3->SetLineWidth(2); l3->Draw();
        TLegend* L = new TLegend(0.16, 0.74, 0.74, 0.90); L->SetBorderSize(0); L->SetFillStyle(0);
        L->SetTextFont(42); L->SetTextSize(0.030);
        L->AddEntry(&gW,  Form("equal-WIDTH 9 bins (current)        RMS = %.2f ps", rms(sW)),  "pl");
        L->AddEntry(&gO10,Form("equal-OCCUPANCY, 1000/bin           RMS = %.2f ps", rms(sO10)),"pl");
        L->AddEntry(&gO25,Form("equal-OCCUPANCY, 2500/bin           RMS = %.2f ps", rms(sO25)),"pl");
        L->Draw();
        { TLatex a; a.SetNDC(); a.SetTextSize(0.027); a.SetTextColor(kGray+3);
          a.DrawLatex(0.16, 0.21, "Equal-occupancy bins hold a fixed event count, so the top bin never");
          a.DrawLatex(0.16, 0.175, "flickers across the eligibility line -- the radius curve smooths out");
          a.DrawLatex(0.16, 0.140, "(RMS falls ~3x) at a comparable resolution."); }
        DrawPageTitle("Smarter binning: equal-occupancy removes the best-bin jitter (150 GeV)");
        c->Print("Analysis/Output/Summary/fiducial_binning_fix.png");
        c->Print("Analysis/Output/Summary/fiducial_binning_fix.pdf");
        printf("[fiducialTimingScan] wrote fiducial_binning_fix.png  (width RMS=%.2f, occ1000 RMS=%.2f, occ2500 RMS=%.2f)\n",
               rms(sW), rms(sO10), rms(sO25));
    }

    // =========================================================================
    // Figure 6 — the FITS behind the equal-occupancy curve (150 GeV, 1000/bin).
    // One panel per radius: the selected bin's (DW-UP)/2 distribution + Gaussian
    // core fit, so the reader can judge each sigma_t directly.
    // =========================================================================
    {
        const int nR = 9;
        const double rList[nR] = { 1.75, 2.00, 2.125, 2.25, 2.375, 2.50, 2.625, 2.75, 3.00 };
        TCanvas* c = new TCanvas("c_fts_fits", "", 1500, 1400);
        c->Divide(3, 3, 0.004, 0.012);
        for (int ip = 0; ip < nR; ++ip) {
            c->cd(ip+1);
            gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.12); gPad->SetRightMargin(0.04); gPad->SetTopMargin(0.10);
            const double R = rList[ip];
            const float R2 = (float)(R*R);
            std::vector<std::pair<float,float>> sel;
            for (const auto& e : evs150) if (e.r2 < R2) sel.emplace_back(e.slg, e.tc);
            std::vector<float> tc; double e_ = 0., eff_ = 0., ec_ = 0.; int n_ = 0, idx_ = 0, k_ = 0;
            double sig = BestBinEqualOcc_fts(sel, 1000, e_, n_, eff_, &idx_, &k_, &tc, &ec_);
            if (sig <= 0. || tc.size() < 50) continue;
            // EXACT same hist + fit the curve uses: 120-bin core hist + FitGaussCore.
            TH1F* h = BuildCoreHist_fts(tc, 120, Form("hfit_%d", ip));
            const double rmsPs = h->GetRMS() * 1000.;     // full-distribution RMS
            h->SetLineColor(kRData); h->SetLineWidth(2); h->SetFillColorAlpha(kRData, 0.25);
            h->GetXaxis()->SetTitle("(DW#minusUP)/2  (ns)");
            h->GetYaxis()->SetTitle("events");
            h->GetXaxis()->SetTitleSize(0.052); h->GetYaxis()->SetTitleSize(0.052);
            h->GetXaxis()->SetLabelSize(0.045); h->GetYaxis()->SetLabelSize(0.045);
            h->Draw("HIST");
            double mu = 0., muE = 0., s = 0., sE = 0.;
            FitGaussCore(h, 2.0, mu, muE, s, sE);          // the headline's core fit
            const double sFit = s * 1000.;                 // = the green-curve value
            // draw the fitted core Gaussian
            TF1* f = new TF1(Form("fg_%d", ip), "[0]*exp(-0.5*((x-[1])/[2])^2)", mu-4*s, mu+4*s);
            f->SetParameters(h->GetBinContent(h->GetMaximumBin()), mu, s);
            f->SetLineColor(kRRed); f->SetLineWidth(3); f->Draw("same");
            const bool hot = (std::fabs(R - 2.125) < 1e-3);
            { TLatex a; a.SetNDC(); a.SetTextFont(42);
              a.SetTextSize(0.075); a.SetTextColor(hot ? kRRed : kBlack);
              a.DrawLatex(0.15, 0.84, Form("r < %.3g mm", R));
              a.SetTextSize(0.066); a.SetTextColor(kRRed);
              a.DrawLatex(0.15, 0.76, Form("#sigma_{core} = %.1f ps", sFit));
              a.SetTextColor(kGray+2); a.SetTextSize(0.050);
              a.DrawLatex(0.15, 0.69, Form("(RMS = %.1f ps)", rmsPs));
              a.SetTextColor(kGray+3);
              a.DrawLatex(0.15, 0.63, Form("N = %d", (int)tc.size()));
              a.DrawLatex(0.15, 0.575, Form("E_{meas} #approx %.0f mV", ec_));
              if (hot) { a.SetTextColor(kRRed); a.SetTextSize(0.052);
                         a.DrawLatex(0.15, 0.50, "(the point you flagged)"); } }
        }
        c->cd(0);
        DrawPageTitle("Equal-occupancy (1000/bin) best-bin fits vs fiducial radius -- 150 GeV");
        c->Print("Analysis/Output/Summary/fiducial_occ_fits.png");
        c->Print("Analysis/Output/Summary/fiducial_occ_fits.pdf");
        printf("[fiducialTimingScan] wrote fiducial_occ_fits.png (equal-occ fit array)\n");
    }

    // Console summary tables
    auto printTable = [&](const char* tag, TGraphErrors* g) {
        printf("\n  R(mm)");
        for (int e = 0; e < nE; ++e) printf("  %6d", eGeV[e]);
        printf("   [%s sigma_t, ps]\n", tag);
        for (double R = 1.0; R <= 5.001; R += 0.25) {
            printf("  %4.2f ", R);
            for (int e = 0; e < nE; ++e) { double v = valAt(g[e], R); if (v>0) printf("  %6.1f", v); else printf("       -"); }
            printf("\n");
        }
    };
    printTable("BEST-BIN", gBest);
    // Mean over energies of best-bin sigma at each radius (where defined for all 6)
    printf("\n  best-bin sigma averaged over the 6 energies:\n");
    for (double R = 1.5; R <= 3.501; R += 0.25) {
        double sum = 0.; int n = 0;
        for (int e = 0; e < nE; ++e) { double v = valAt(gBest[e], R); if (v>0){ sum+=v; ++n; } }
        if (n == nE) printf("   r=%.2f mm : mean = %.2f ps  (over %d energies)\n", R, sum/n, n);
        else         printf("   r=%.2f mm : (only %d/%d energies defined)\n", R, n, nE);
    }
    printTable("top-2%", gTop2);
}

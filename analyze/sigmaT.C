// ============================================================================
// sigmaT.C — build-agnostic timing resolution driver (uses the shared method)
// ----------------------------------------------------------------------------
// Thin driver over rad::timingBestBin (analyze/RadTiming.h), which IS the
// published headline (DW-UP)/2 Method A pipeline. The same rigorous method is
// applied to every build/run/year. Reproduces DSB1/150 = 27.4 ps.
//
//   ROOT_INCLUDE_PATH=lib/waveform:lib/io:lib/physics:lib/viz root -l -b -q \
//     'analyze/sigmaT.C+("data/2023/configs/DSB1.json", 150)'
// ============================================================================
#include "RadTiming.h"       // radcore: timingBestBin (the one method)
#include "DataPaths.h"       // Analysis: radReduced
#include "TFile.h"
#include "TTree.h"
#include <cstdio>

void sigmaT(const char* configPath, double energy = 150) {
    rad::BuildConfig cfg = rad::BuildConfig::Load(configPath);
    if (!cfg.valid()) { printf("config load failed: %s\n", cfg.error()); return; }

    TFile* fp = TFile::Open(radReduced(cfg.build.c_str(), energy));
    if (!fp || fp->IsZombie()) { printf("no reduced file for %s @ %.0f\n", cfg.build.c_str(), energy); return; }
    TTree* t = (TTree*)fp->Get("rad");
    rad::RadView v; v.attach(t, &cfg);

    int nT = 0; for (int i=0;i<cfg.nend;++i) if (v.is_timing(i)) ++nT;
    rad::TimingResult r = rad::timingBestBin(v, energy);

    printf("build %-8s @ %.0f GeV [%-9s]  %d timing ends; centroid (%.2f,%.2f) rFid %.1f mm\n",
           cfg.build.c_str(), energy, v.named?"canonical":"slots", nT, r.xc, r.yc, r.rFid);
    printf("  fiducial=%zu  muE=%.0f sigE=%.0f mV  best bin %d (sum_lg~%.0f)\n",
           r.nFid, r.muE, r.sigE, r.best_bin, r.bestE);
    printf("  (DW-UP)/2 best-bin sigma_t = %.1f ps\n", r.sigma_ps);
    fp->Close();
}

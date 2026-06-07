// ============================================================================
// timingLadder.C — cross-build sigma_t(E) via the ONE shared headline method.
// ----------------------------------------------------------------------------
// Runs rad::timingBestBin (the published (DW-UP)/2 Method A pipeline) for every
// build at every energy, so the whole table is produced by the identical
// world-class method. Paper-grade cross-build timing comparison.
//
//   ROOT_INCLUDE_PATH=lib/waveform:lib/io:lib/physics:lib/viz root -l -b -q analyze/timingLadder.C+
// ============================================================================
#include "RadTiming.h"
#include "DataPaths.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include <cstdio>

void timingLadder() {
    const char* builds[4] = {"DSB1","LUAG","MIXED","TENERGY"};
    double Es[6] = {25,50,75,100,125,150};

    printf("\n=== sigma_t(E) [ps] — headline (DW-UP)/2 method, all builds (energy caps excluded by role) ===\n");
    printf("%-9s", "build");
    for (double E : Es) printf("%8.0f", E);
    printf("   method = rad::timingBestBin (== timingEnergyBins.C Method A)\n");

    for (int b = 0; b < 4; ++b) {
        rad::BuildConfig cfg = rad::BuildConfig::Load(Form("data/2023/configs/%s.json", builds[b]));
        printf("%-9s", builds[b]);
        for (double E : Es) {
            TString f = radReduced(builds[b], E);
            TFile* fp = TFile::Open(f);
            if (!fp || fp->IsZombie()) { printf("%8s", "-"); continue; }
            TTree* t = (TTree*)fp->Get("rad");
            if (!t) { printf("%8s", "-"); fp->Close(); continue; }
            rad::RadView v; v.attach(t, &cfg);
            rad::TimingResult r = rad::timingBestBin(v, E);
            if (r.sigma_ps > 0) printf("%8.1f", r.sigma_ps); else printf("%8s", "-");
            fp->Close();
        }
        printf("\n");
    }
    printf("(DSB1 published headline ladder: 47.1 32.7 30.8 30.2 29.2 27.4 ps)\n");
}

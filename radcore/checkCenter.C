// checkCenter.C — does RadView.beamCenter match the real ScanRunCenters?
#include "RadView.h"
#include "PlotUtils.h"     // ScanRunCenters
#include "DataPaths.h"
#include "TFile.h"
#include "TTree.h"
#include <cstdio>
void checkCenter() {
    TFile* fp = TFile::Open(radReduced("DSB1",150));
    TTree* t = (TTree*)fp->Get("rad");
    double xc,yc,to,tr; ScanRunCenters(t, xc, yc, to, tr);
    printf("ScanRunCenters    : (%.3f, %.3f) mm   CFD off %.2f\n", xc, yc, to);
    rad::BuildConfig cfg = rad::BuildConfig::Load("datasets/2023/configs/DSB1.json");
    rad::RadView v; v.attach(t, &cfg);
    double mx,my; v.beamCenter(mx,my);
    printf("RadView.beamCenter: (%.3f, %.3f) mm\n", mx, my);
    fp->Close();
}

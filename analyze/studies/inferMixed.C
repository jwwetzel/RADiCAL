// inferMixed.C — infer which MIXED corners are DSB1 vs LuAG from the data.
// Per corner, compare MIXED timing-cap brightness to pure DSB1 and pure LUAG
// (same corner, so SiPM/geometry cancels); assign the closer material.
//   ROOT_INCLUDE_PATH=lib/waveform:lib/io:lib/physics:lib/viz root -l -b -q analyze/studies/inferMixed.C+
#include "BuildConfig.h"
#include "Schema.h"
#include "DataPaths.h"
#include "TFile.h"
#include "TTree.h"
#include <string>
#include <cmath>
#include <cstdio>

static void cornerBright(const char* build, double E, const rad::BuildConfig& map, double out[4]) {
    TFile* fp = TFile::Open(radReduced(build, E));
    if (!fp || fp->IsZombie()) { for (int c=0;c<4;++c) out[c]=-1; return; }
    TTree* t = (TTree*)fp->Get("rad");
    rad::RadEvent ev; ev.ConnectBranches(t);
    bool isDSB1 = (std::string(build) == "DSB1");      // DSB1 file = processRun format
    double sum[4]={0}; long n[4]={0}; long N=t->GetEntries();
    for (long i=0;i<N;++i){ t->GetEntry(i);
        if (!ev.wc_ok || ev.mcp1_peak<200 || ev.mcp1_peak>750) continue;
        for (int c=0;c<4;++c){
            double d = isDSB1 ? ev.hg_peak[c]   : ev.s_peak[map.end[c].hg/1024];
            double u = isDSB1 ? ev.hg_peak[c+4] : ev.s_peak[map.end[c+4].hg/1024];
            if (d>5 && u>5){ sum[c]+=0.5*(d+u); ++n[c]; }
        }
    }
    for (int c=0;c<4;++c) out[c] = n[c] ? sum[c]/n[c] : 0;
    fp->Close();
}

void inferMixed() {
    rad::BuildConfig map = rad::BuildConfig::Load("data/2023/configs/DSB1.json");
    double dsb1[4], luag[4], mix[4];
    cornerBright("DSB1", 150, map, dsb1);
    cornerBright("LUAG", 150, map, luag);
    cornerBright("MIXED",150, map, mix);
    const char* cn[4] = {"NW","NE","SE","SW"};
    printf("\nper-corner mean timing-cap peak [mV] @150 GeV\n");
    printf("%4s %10s %10s %10s    %s\n","corner","DSB1","LUAG","MIXED","-> MIXED is");
    int nD=0, nL=0;
    for (int c=0;c<4;++c){
        double dD=std::fabs(mix[c]-dsb1[c]), dL=std::fabs(mix[c]-luag[c]);
        bool isD = dD < dL;
        printf("%4s %10.1f %10.1f %10.1f    -> %s\n", cn[c], dsb1[c], luag[c], mix[c], isD?"DSB1":"LuAG");
        if (isD) ++nD; else ++nL;
    }
    printf("\ninferred MIXED fill: %d DSB1 + %d LuAG  (%s)\n", nD, nL,
           (nD==2 && nL==2) ? "consistent with '2xDSB1, 2xLuAG'" : "CHECK — not 2+2");
}

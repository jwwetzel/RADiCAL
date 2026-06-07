// gateDSB1.C — read the existing full-stat DSB1 file via the canonical reader
// and compute a sanity sigma_t (exercises the tolerant mcp_peak->mcp1_peak bind).
//   ROOT_INCLUDE_PATH=lib/waveform:lib/io:lib/physics:lib/viz root -l -b -q analyze/studies/gateDSB1.C+
#include "Schema.h"
#include "TFile.h"
#include "TTree.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>

void gateDSB1(const char* f = "data/2023/reduced/DSB1/150GeV.root") {
    TFile* fp = TFile::Open(f); if (!fp || fp->IsZombie()) { printf("no file %s\n", f); return; }
    TTree* t = (TTree*)fp->Get("rad");
    rad::RadEvent ev; ev.ConnectBranches(t);
    long N = t->GetEntries();
    double xs=0, ys=0; long nw=0;
    for (long i=0;i<N && nw<50000;++i){ t->GetEntry(i); if (ev.wc_ok && ev.x_trk>-100 && ev.x_trk<100){xs+=ev.x_trk;ys+=ev.y_trk;++nw;} }
    double xc=xs/nw, yc=ys/nw;
    std::vector<float> tv; double mcpsum=0; long mok=0;
    for (long i=0;i<N;++i){ t->GetEntry(i);
        if (!ev.wc_ok || ev.mcp1_peak<200 || ev.mcp1_peak>750) continue;   // mcp1 via tolerant reader
        double dx=ev.x_trk-xc, dy=ev.y_trk-yc; if (dx*dx+dy*dy>=9.0) continue;
        double ds=0,us=0; int dn=0,un=0;
        for (int c=0;c<4;++c) if (ev.hg_cfd05[c]>-1e5){ ds+=ev.hg_cfd05[c]; ++dn; }
        for (int c=4;c<8;++c) if (ev.hg_cfd05[c]>-1e5){ us+=ev.hg_cfd05[c]; ++un; }
        if (dn>=1 && un>=1) tv.push_back(0.5f*(float)(ds/dn-us/un));
        mcpsum+=ev.mcp1_peak; ++mok;
    }
    std::sort(tv.begin(),tv.end()); double mu=tv[tv.size()/2], s=0.1;
    for (int it=0;it<6;++it){ double a=0,a2=0; long n=0; for(float x:tv) if(std::fabs(x-mu)<2.5*s){a+=x;a2+=x*x;++n;} if(n<50)break; mu=a/n; double v=a2/n-mu*mu; s=v>0?std::sqrt(v):s; }
    printf("entries=%ld  fiducial-timing events=%zu\n", N, tv.size());
    printf("mcp1 via tolerant reader: mean=%0.f mV over %ld evts  => %s\n",
           mok?mcpsum/mok:0, mok, mok>0 ? "reader OK" : "READER FAILED");
    printf("all-fiducial (DW-UP)/2 core sigma_t = %.1f ps  (published headline = 27.4 ps, best-bin OOS)\n", s*1000);
    fp->Close();
}

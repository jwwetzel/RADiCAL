// ============================================================================
// verify.C — post-reduction cleanliness gate. For each merged reduced file,
// confirm the NEW branches are present (hg_lgcfd, tr0a_time, tr0b_time) and
// report entries per energy. Run this on Argon after the merge, BEFORE rsyncing
// home — a stale file (old reduction lingering) shows hg_lgcfd=N and is the
// signal to re-merge rather than ship a mixed bag.
//
//   ROOT_INCLUDE_PATH=lib/waveform:lib/io:lib/physics:lib/viz root -l -b -q 'reduce/verify.C+("data/2023/reduced")'
//   # or a single build:  ...'reduce/verify.C+("data/2023/reduced","DSB1")'
// ============================================================================
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TString.h"
#include <cstdio>
#include <vector>
#include <string>
#include <algorithm>

static bool hasBranch(TTree* t, const char* b){ return t && t->GetBranch(b) != nullptr; }

static void verifyBuild(const TString& root, const char* build){
    TString dir = Form("%s/%s", root.Data(), build);
    void* dp = gSystem->OpenDirectory(dir);
    if (!dp){ printf("  %-8s  (no directory %s)\n", build, dir.Data()); return; }
    std::vector<std::string> files;
    const char* e;
    while ((e = gSystem->GetDirEntry(dp))){ TString f(e); if (f.EndsWith("GeV.root")) files.push_back(e); }
    gSystem->FreeDirectory(dp);
    std::sort(files.begin(), files.end());
    if (files.empty()){ printf("  %-8s  (empty)\n", build); return; }
    printf("  %-8s\n", build);
    for (auto& fn : files){
        TString path = Form("%s/%s", dir.Data(), fn.c_str());
        TFile* f = TFile::Open(path);
        if (!f || f->IsZombie()){ printf("    %-12s  CANNOT OPEN\n", fn.c_str()); if(f) f->Close(); continue; }
        TTree* t = (TTree*)f->Get("rad");
        if (!t){ printf("    %-12s  no 'rad' tree\n", fn.c_str()); f->Close(); continue; }
        bool lg = hasBranch(t,"hg_lgcfd"), ta = hasBranch(t,"tr0a_time"), tb = hasBranch(t,"tr0b_time");
        bool ok = lg && ta && tb;
        printf("    %-12s  %9lld ev   hg_lgcfd=%s tr0a_time=%s tr0b_time=%s   %s\n",
               fn.c_str(), (long long)t->GetEntries(),
               lg?"Y":"N", ta?"Y":"N", tb?"Y":"N", ok?"OK":"*** STALE / RE-MERGE ***");
        f->Close();
    }
}

void verify(const char* reducedRoot = "data/2023/reduced", const char* build = ""){
    TString root(reducedRoot);
    printf("[verify] %s\n", root.Data());
    if (build && build[0]){ verifyBuild(root, build); return; }
    const char* builds[] = {"DSB1","LUAG","MIXED","TENERGY"};
    for (const char* b : builds) verifyBuild(root, b);
    printf("[verify] All builds must show hg_lgcfd=Y tr0a_time=Y tr0b_time=Y.\n");
    printf("         DSB1 should carry 25 GeV; LUAG/MIXED/TENERGY start at 50.\n");
}

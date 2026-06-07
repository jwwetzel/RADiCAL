// ============================================================================
// validateReduce.C — prove the unified Reducer reproduces processRun.C exactly
// ----------------------------------------------------------------------------
// Compares a NEW canonical reduced file (radcore/Reducer.C) against the OLD
// processRun.C DSB1 ntuple, event-by-event, on every shared branch. They must
// be identical (the only difference is the mcp_peak->mcp1_peak rename and the
// added generic s_peak[36] arrays).
//
//   ROOT_INCLUDE_PATH=radcore:Analysis root -l -b -q \
//     'radcore/validateReduce.C+("data/2023/reduced/DSB1/25GeV.root","data/2023/reduced/_validate/25GeV.root")'
// ============================================================================
#include "Schema.h"
#include "TFile.h"
#include "TTree.h"
#include <cstdio>
#include <cmath>
#include <algorithm>

void validateReduce(const char* oldF, const char* newF) {
    TFile* fo = TFile::Open(oldF); TFile* fn = TFile::Open(newF);
    if (!fo || fo->IsZombie() || !fn || fn->IsZombie()) { printf("cannot open inputs\n"); return; }
    TTree* to = (TTree*)fo->Get("rad"); TTree* tn = (TTree*)fn->Get("rad");
    Long64_t no = to->GetEntries(), nn = tn->GetEntries();
    printf("entries: old=%lld  new=%lld  %s\n", no, nn, no==nn ? "(match)" : "(MISMATCH!)");

    // OLD (legacy processRun) branch names
    Int_t orun, oevt; Bool_t owc; Float_t ox, oy, oslg, ospb, omcp, omcp2t;
    Float_t ohgcfd[8], ohgpk[8], olgpk[8], ohgq[8], ohg03[8], olgq[8];
    to->SetBranchAddress("run",&orun); to->SetBranchAddress("event",&oevt);
    to->SetBranchAddress("wc_ok",&owc); to->SetBranchAddress("x_trk",&ox); to->SetBranchAddress("y_trk",&oy);
    to->SetBranchAddress("sum_lg",&oslg); to->SetBranchAddress("sum_pb",&ospb);
    to->SetBranchAddress("mcp_peak",&omcp); to->SetBranchAddress("mcp2_time",&omcp2t);
    to->SetBranchAddress("hg_cfd05",ohgcfd); to->SetBranchAddress("hg_cfd03",ohg03);
    to->SetBranchAddress("hg_peak",ohgpk); to->SetBranchAddress("lg_peak",olgpk);
    to->SetBranchAddress("hg_charge",ohgq); to->SetBranchAddress("lg_charge",olgq);

    rad::RadEvent ev; ev.ConnectBranches(tn);

    Long64_t N = std::min(no, nn);
    long meta_bad = 0, wc_bad = 0;
    double d_mcp=0, d_slg=0, d_spb=0, d_xy=0, d_cfd=0, d_cfd03=0, d_pk=0, d_lg=0, d_q=0, d_lq=0;
    auto up = [](double& m, double v){ if (v > m) m = v; };

    for (Long64_t i = 0; i < N; ++i) {
        to->GetEntry(i); tn->GetEntry(i);
        if (orun != ev.run || oevt != ev.event) ++meta_bad;
        if (owc != ev.wc_ok) ++wc_bad;
        up(d_mcp, std::fabs(omcp - ev.mcp1_peak));
        up(d_slg, std::fabs(oslg - ev.sum_lg));
        up(d_spb, std::fabs(ospb - ev.sum_pb));
        up(d_xy,  std::max(std::fabs(ox - ev.x_trk), std::fabs(oy - ev.y_trk)));
        for (int c = 0; c < 8; ++c) {
            up(d_cfd,  std::fabs(ohgcfd[c] - ev.hg_cfd05[c]));
            up(d_cfd03,std::fabs(ohg03[c]  - ev.hg_cfd03[c]));
            up(d_pk,   std::fabs(ohgpk[c]  - ev.hg_peak[c]));
            up(d_lg,   std::fabs(olgpk[c]  - ev.lg_peak[c]));
            up(d_q,    std::fabs(ohgq[c]   - ev.hg_charge[c]));
            up(d_lq,   std::fabs(olgq[c]   - ev.lg_charge[c]));
        }
    }
    printf("\nmax |old-new| over %lld events:\n", N);
    printf("  meta(run,event) mismatches : %ld\n", meta_bad);
    printf("  wc_ok mismatches           : %ld\n", wc_bad);
    printf("  mcp peak                   : %.3g\n", d_mcp);
    printf("  sum_lg / sum_pb            : %.3g / %.3g\n", d_slg, d_spb);
    printf("  x_trk/y_trk                : %.3g\n", d_xy);
    printf("  hg_cfd05 / hg_cfd03 [ns]   : %.3g / %.3g\n", d_cfd, d_cfd03);
    printf("  hg_peak / lg_peak          : %.3g / %.3g\n", d_pk, d_lg);
    printf("  hg_charge / lg_charge      : %.3g / %.3g\n", d_q, d_lq);

    double worst = std::max({d_mcp,d_slg,d_spb,d_xy,d_cfd,d_cfd03,d_pk,d_lg,d_q,d_lq});
    bool ok = (no==nn) && meta_bad==0 && wc_bad==0 && worst < 1e-3;
    printf("\n==== %s (worst diff = %.3g) ====\n",
           ok ? "IDENTICAL — unified Reducer reproduces processRun.C" : "DIFFERENCES FOUND", worst);
    fo->Close(); fn->Close();
}

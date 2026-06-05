// ============================================================================
// Schema.h — the ONE canonical RADiCAL reduced-ntuple schema (TTree "rad")
// ----------------------------------------------------------------------------
// A strict SUPERSET of the two legacy formats, so a single file feeds every
// analysis:
//   * DSB1-style (processRun.C): hg_peak/hg_cfd03/05/10/20/30/50/led/tot/charge,
//     hg_ped_rms, hg_saturated/spike, lg_peak/charge, sum_lg/sum_pb, pb_peak,
//     wc_ok/x_trk/y_trk/wc_peak, mcp(1)_peak/time, mcp2_peak/time, stopcell.
//   * config-style (reduceRaw.C): generic per-slot s_peak[36]/s_cfd05[36]/
//     s_charge[36] for rediscovery & diagnostics.
//
// The per-capillary arrays [8] are stored in a fixed CANONICAL END ORDER
// (NW-D, NE-D, SE-D, SW-D, NW-U, NE-U, SE-U, SW-U) defined by the build config,
// so index i always means the same physical corner/end across builds and years.
//
// The single naming change vs the legacy DSB1 format: the MCP1 reference is
// named mcp1_peak / mcp1_time (legacy DSB1 called it mcp_peak / mcp_time).
//
// RadEvent holds all buffers. CreateBranches(tree) wires them for WRITING
// (the reducer); ConnectBranches(tree) wires them for READING (analyses).
// ============================================================================
#ifndef RADCORE_SCHEMA_H
#define RADCORE_SCHEMA_H

#include "TTree.h"
#include "TFile.h"
#include "TNamed.h"

namespace rad {

static const int NCAP  = 8;   // 4 capillaries x 2 ends (down/up) = 8 SiPM readouts
static const int NSLOT = 36;  // 2 DRS4 boards x 2 groups x 9 channels
static const int NWC   = 4;   // wire-chamber planes: R, L, D, U
static const int NPB   = 4;   // PbGlass reference blocks
static const int NSTOP = 4;   // DRS4 stop cells: D0G0, D0G1, D1G0, D1G1
static const float kNoTime = -1.0e6f;  // sentinel (matches WaveformUtils.h)

// canonical end order for the [8] arrays
static const char* const kEndName[NCAP] =
    { "NW-D","NE-D","SE-D","SW-D","NW-U","NE-U","SE-U","SW-U" };

// ---------------------------------------------------------------------------
struct RadEvent {
    // --- metadata ---
    Int_t   run = 0, event = 0;
    Float_t beam_energy = 0.f;

    // --- wire chamber (config-invariant) ---
    Bool_t  wc_ok = false, in_fiducial = false;
    Float_t x_trk = kNoTime, y_trk = kNoTime;
    Float_t wc_peak[NWC]   = {0};

    // --- MCP timing references (config-invariant) ---
    Float_t mcp1_peak = 0.f, mcp1_time = kNoTime;   // legacy DSB1 name: mcp_peak/mcp_time
    Float_t mcp2_peak = 0.f, mcp2_time = kNoTime;

    // --- DRS4 stop cells (D0G0,D0G1,D1G0,D1G1) ---
    Int_t   stopcell[NSTOP] = {0};

    // --- per capillary end [8] (role-resolved, canonical order) ---
    Float_t hg_peak[NCAP]    = {0};
    Float_t hg_ped_rms[NCAP] = {0};
    Bool_t  hg_saturated[NCAP] = {false};
    Bool_t  hg_spike[NCAP]   = {false};
    Float_t hg_cfd03[NCAP]   = {0};   // MCP-referenced HG CFD times [ns]
    Float_t hg_cfd05[NCAP]   = {0};
    Float_t hg_cfd[NCAP]     = {0};   // CFD-20% (legacy name kept)
    Float_t hg_cfd10[NCAP]   = {0};
    Float_t hg_cfd30[NCAP]   = {0};
    Float_t hg_cfd50[NCAP]   = {0};
    Float_t hg_led[NCAP]     = {0};
    Float_t hg_tot[NCAP]     = {0};
    Float_t hg_charge[NCAP]  = {0};
    Float_t lg_peak[NCAP]    = {0};
    Float_t lg_charge[NCAP]  = {0};

    // --- derived energy proxies ---
    Float_t sum_lg = 0.f, sum_pb = 0.f;
    Float_t pb_peak[NPB] = {0};

    // --- generic per-slot features (rediscovery / diagnostics) ---
    Float_t s_peak[NSLOT]   = {0};
    Float_t s_cfd05[NSLOT]  = {0};
    Float_t s_charge[NSLOT] = {0};

    // ----------------------------------------------------------------------
    // WRITE side — create all branches on a fresh tree
    void CreateBranches(TTree* t) {
        t->Branch("run",         &run,         "run/I");
        t->Branch("event",       &event,       "event/I");
        t->Branch("beam_energy", &beam_energy, "beam_energy/F");
        t->Branch("wc_ok",       &wc_ok,       "wc_ok/O");
        t->Branch("x_trk",       &x_trk,       "x_trk/F");
        t->Branch("y_trk",       &y_trk,       "y_trk/F");
        t->Branch("in_fiducial", &in_fiducial, "in_fiducial/O");
        t->Branch("wc_peak",      wc_peak,     "wc_peak[4]/F");
        t->Branch("mcp1_peak",   &mcp1_peak,   "mcp1_peak/F");
        t->Branch("mcp1_time",   &mcp1_time,   "mcp1_time/F");
        t->Branch("mcp2_peak",   &mcp2_peak,   "mcp2_peak/F");
        t->Branch("mcp2_time",   &mcp2_time,   "mcp2_time/F");
        t->Branch("stopcell",     stopcell,    "stopcell[4]/I");
        t->Branch("hg_peak",      hg_peak,     "hg_peak[8]/F");
        t->Branch("hg_ped_rms",   hg_ped_rms,  "hg_ped_rms[8]/F");
        t->Branch("hg_saturated", hg_saturated,"hg_saturated[8]/O");
        t->Branch("hg_spike",     hg_spike,    "hg_spike[8]/O");
        t->Branch("hg_cfd03",     hg_cfd03,    "hg_cfd03[8]/F");
        t->Branch("hg_cfd05",     hg_cfd05,    "hg_cfd05[8]/F");
        t->Branch("hg_cfd",       hg_cfd,      "hg_cfd[8]/F");
        t->Branch("hg_cfd10",     hg_cfd10,    "hg_cfd10[8]/F");
        t->Branch("hg_cfd30",     hg_cfd30,    "hg_cfd30[8]/F");
        t->Branch("hg_cfd50",     hg_cfd50,    "hg_cfd50[8]/F");
        t->Branch("hg_led",       hg_led,      "hg_led[8]/F");
        t->Branch("hg_tot",       hg_tot,      "hg_tot[8]/F");
        t->Branch("hg_charge",    hg_charge,   "hg_charge[8]/F");
        t->Branch("lg_peak",      lg_peak,     "lg_peak[8]/F");
        t->Branch("lg_charge",    lg_charge,   "lg_charge[8]/F");
        t->Branch("sum_lg",      &sum_lg,      "sum_lg/F");
        t->Branch("sum_pb",      &sum_pb,      "sum_pb/F");
        t->Branch("pb_peak",      pb_peak,     "pb_peak[4]/F");
        t->Branch("s_peak",       s_peak,      "s_peak[36]/F");
        t->Branch("s_cfd05",      s_cfd05,     "s_cfd05[36]/F");
        t->Branch("s_charge",     s_charge,    "s_charge[36]/F");
    }

    // ----------------------------------------------------------------------
    // READ side — connect to an existing tree (only branches that exist)
    void ConnectBranches(TTree* t) {
        auto SB = [&](const char* n, void* p){ if (t->GetBranch(n)) t->SetBranchAddress(n, p); };
        SB("run",&run); SB("event",&event); SB("beam_energy",&beam_energy);
        SB("wc_ok",&wc_ok); SB("x_trk",&x_trk); SB("y_trk",&y_trk);
        SB("in_fiducial",&in_fiducial); SB("wc_peak",wc_peak);
        SB("mcp1_peak",&mcp1_peak); SB("mcp1_time",&mcp1_time);
        SB("mcp2_peak",&mcp2_peak); SB("mcp2_time",&mcp2_time);
        // transitional: legacy DSB1 (processRun.C) files name the MCP1 reference
        // mcp_peak/mcp_time — bind those to mcp1_* so the canonical reader works
        // on both old and new files during the migration.
        if (!t->GetBranch("mcp1_peak") && t->GetBranch("mcp_peak")) SB("mcp_peak",&mcp1_peak);
        if (!t->GetBranch("mcp1_time") && t->GetBranch("mcp_time")) SB("mcp_time",&mcp1_time);
        SB("stopcell",stopcell);
        SB("hg_peak",hg_peak); SB("hg_ped_rms",hg_ped_rms);
        SB("hg_saturated",hg_saturated); SB("hg_spike",hg_spike);
        SB("hg_cfd03",hg_cfd03); SB("hg_cfd05",hg_cfd05); SB("hg_cfd",hg_cfd);
        SB("hg_cfd10",hg_cfd10); SB("hg_cfd30",hg_cfd30); SB("hg_cfd50",hg_cfd50);
        SB("hg_led",hg_led); SB("hg_tot",hg_tot); SB("hg_charge",hg_charge);
        SB("lg_peak",lg_peak); SB("lg_charge",lg_charge);
        SB("sum_lg",&sum_lg); SB("sum_pb",&sum_pb); SB("pb_peak",pb_peak);
        SB("s_peak",s_peak); SB("s_cfd05",s_cfd05); SB("s_charge",s_charge);
    }
};

// Stamp/read the build name on the output file (so a reduced file is
// self-describing about which configuration produced it).
inline void StampBuild(TFile* f, const char* build, const char* year) {
    f->cd();
    TNamed("rad_build", build).Write();
    TNamed("rad_year",  year ).Write();
}

} // namespace rad

#endif // RADCORE_SCHEMA_H

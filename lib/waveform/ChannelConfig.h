// ============================================================================
// ChannelConfig.h — RADiCAL CERN May 2023 Test Beam
// ============================================================================
//
// Two CAEN DT5742 digitizers (DRS0, DRS1).  Each module has 2 inner groups
// of 9 channels (8 signal + 1 trigger).  Data is stored flat per event:
//
//   channel_offset(drs, grp, ch) = (1024*9*2)*drs + (1024*9)*grp + 1024*ch
//   time_offset(drs, grp)        = (1024*2)*drs + 1024*grp
//
// The amplitude array has  2*2*9*1024 = 36864 floats per event.
// The timevalue array has  2*2*1024   =  4096 floats per event.
// ============================================================================

#ifndef CHANNELCONFIG_H
#define CHANNELCONFIG_H

#include "SelectionCuts.h"   // all numerical thresholds — single source of truth
#include "TString.h"

// ---------------------------------------------------------------------------
// Index helpers
// ---------------------------------------------------------------------------
inline int chanOff(int drs, int grp, int ch) {
    return (1024*9*2)*drs + (1024*9)*grp + 1024*ch;
}
inline int timeOff(int drs, int grp) {
    return (1024*2)*drs + 1024*grp;
}

// Time-axis offsets — one shared axis per DRS inner group
//   DRS0 Group0: channels  0-8  (high-gain timing capillaries + MCP1)
//   DRS0 Group1: channels  9-17 (SW-Up capillary, PbGlass, MCP2)
//   DRS1 Group0: channels  0-8  (low-gain energy capillaries)
//   DRS1 Group1: channels  9-17 (wire chambers, scintillators)
static const int kT_D0G0 = timeOff(0,0);  // = 0
static const int kT_D0G1 = timeOff(0,1);  // = 1024
static const int kT_D1G0 = timeOff(1,0);  // = 2048
static const int kT_D1G1 = timeOff(1,1);  // = 3072

// ---------------------------------------------------------------------------
// MCP (Micro Channel Plate — beam timing reference)
// ---------------------------------------------------------------------------
static const int kMCP1     = chanOff(0,0,7);  // DRS0 G0 Ch7   (ch 7 in 0-17 scheme)
static const int kMCP2     = chanOff(0,1,7);  // DRS0 G1 Ch7   (ch16 in 0-17 scheme)
static const int kMCP1_t   = kT_D0G0;
static const int kMCP2_t   = kT_D0G1;

// ---------------------------------------------------------------------------
// Wire chamber channels (DRS1 Group1)
// Position formula: x = kWC_Scale * (t_Right - t_Left)  [mm]
//                   y = kWC_Scale * (t_Down  - t_Up  )  [mm]
// ---------------------------------------------------------------------------
static const int    kWC_R    = chanOff(1,1,1);  // WC-X-Right
static const int    kWC_L    = chanOff(1,1,2);  // WC-X-Left
static const int    kWC_D    = chanOff(1,1,3);  // WC-Y-Down
static const int    kWC_U    = chanOff(1,1,5);  // WC-Y-Up
static const int    kWC_t    = kT_D1G1;
static const double kWC_Scale = 7.0/36.0;       // mm/ns

// kCalo_x0, kCalo_y0, kFiducial_r_energy, kFiducial_r_timing
// are defined in SelectionCuts.h (included above).

// ---------------------------------------------------------------------------
// Capillary pairs — 8 total
// Ordering: NW-D, NE-D, SE-D, SW-D, NW-U, NE-U, SE-U, SW-U
//   HG = high-gain DRS0 channel (WLS at shower max  → timing)
//   LG = low-gain  DRS1 channel (WLS full length    → energy)
// ---------------------------------------------------------------------------
struct CapCfg {
    int    hg;       // amplitude offset, HG channel
    int    hg_t;     // time offset for HG
    int    lg;       // amplitude offset, LG channel
    int    lg_t;     // time offset for LG
    int    mcp;      // amplitude offset, MCP reference for this capillary
    int    mcp_t;    // time offset for the MCP reference
    bool   use_mcp2; // true only for SW-U (shares DRS0 Group1 with MCP2)
    const char* name;
};

static const CapCfg kCap[8] = {
    // idx 0: NW-Down
    { chanOff(0,0,1), kT_D0G0, chanOff(1,0,1), kT_D1G0, kMCP1, kMCP1_t, false, "NW-D" },
    // idx 1: NE-Down
    { chanOff(0,0,2), kT_D0G0, chanOff(1,0,2), kT_D1G0, kMCP1, kMCP1_t, false, "NE-D" },
    // idx 2: SE-Down
    { chanOff(0,0,3), kT_D0G0, chanOff(1,0,3), kT_D1G0, kMCP1, kMCP1_t, false, "SE-D" },
    // idx 3: SW-Down
    { chanOff(0,0,0), kT_D0G0, chanOff(1,0,0), kT_D1G0, kMCP1, kMCP1_t, false, "SW-D" },
    // idx 4: NW-Up
    { chanOff(0,0,5), kT_D0G0, chanOff(1,0,5), kT_D1G0, kMCP1, kMCP1_t, false, "NW-U" },
    // idx 5: NE-Up
    { chanOff(0,0,4), kT_D0G0, chanOff(1,0,4), kT_D1G0, kMCP1, kMCP1_t, false, "NE-U" },
    // idx 6: SE-Up
    { chanOff(0,0,6), kT_D0G0, chanOff(1,0,6), kT_D1G0, kMCP1, kMCP1_t, false, "SE-U" },
    // idx 7: SW-Up  (DRS0 G1 Ch0 = ch9 in 0-17; shares time axis with MCP2)
    { chanOff(0,1,0), kT_D0G1, chanOff(1,0,7), kT_D1G0, kMCP2, kMCP2_t, true,  "SW-U" },
};
static const int kNCap = 8;

// ---------------------------------------------------------------------------
// PbGlass reference channels (DRS0 Group1, channels 1–4)
// ---------------------------------------------------------------------------
static const int kPbGlass[4] = {
    chanOff(0,1,1),  // PB2 SW  (ch10 in 0-17)
    chanOff(0,1,2),  // PB1 NW  (ch11)
    chanOff(0,1,3),  // PB3 NE  (ch12)
    chanOff(0,1,4),  // PB4 SE  (ch13)
};
static const int kPbGlass_t = kT_D0G1;

// ---------------------------------------------------------------------------
// Run list
//
// inFiles is a semicolon-separated list of file paths (or glob patterns).
// processRun.C splits on ';' and adds each entry to a TChain, so multiple
// runs at the same energy are transparently chained together.
//
// Example (single file):   "data/2023/raw/RUN1211_25_GeV.root"
// Example (four runs):     "data/2023/raw/RUN1258_150_GeV.root;data/2023/raw/RUN1259_150_GeV.root;..."
// ---------------------------------------------------------------------------
struct RunCfg {
    double  energy_GeV;
    TString inFiles;   // semicolon-separated list
    TString label;
};

static const RunCfg kRuns[] = {
    {  25., "data/2023/raw/RUN1211_25_GeV.root",  "25GeV"  },
    {  50., "data/2023/raw/RUN1148_50_GeV.root",  "50GeV"  },
    {  75., "data/2023/raw/RUN1112_75_GeV.root",  "75GeV"  },
    { 100., "data/2023/raw/RUN1075_100_GeV.root", "100GeV" },
    { 125., "data/2023/raw/RUN1034_125_GeV.root", "125GeV" },
    // 150 GeV: four runs chained (~8.5 GB total, ~120k events)
    { 150., "data/2023/raw/RUN1258_150_GeV.root"
            ";data/2023/raw/RUN1259_150_GeV.root"
            ";data/2023/raw/RUN1260_150_GeV.root"
            ";data/2023/raw/RUN1261_150_GeV.root", "150GeV" },
};
static const int kNRuns = 6;

#endif // CHANNELCONFIG_H

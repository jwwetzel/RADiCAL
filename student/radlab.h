// ===========================================================================
// radlab.h  —  tiny shared helpers for the RADiCAL student data lab.
//
// Everything you need to turn a raw DRS4 waveform into a number, written out
// the simple way so you can SEE what it does. The "real" analysis uses the more
// careful versions in Analysis/WaveformUtils.h and Analysis/ChannelConfig.h —
// read those once these make sense.
//
// A raw event holds two big flat arrays:
//    amplitude[36864]  =  36 channels ("slots") x 1024 samples   (millivolts)
//    timevalue[4096]   =   4 DRS groups x 1024 samples           (nanoseconds)
// Each slot's 1024 samples start at  slot*1024  inside amplitude.
// Pulses point DOWN (negative-going), so "peak height" = pedestal - minimum.
// ===========================================================================
#pragma once

// Where slot `s` (0..35) lives inside the two flat arrays.
inline int ampOffset (int s){ return s * 1024; }                 // into amplitude[]
inline int timeOffset(int s){ int drs = s/18, grp = (s/9)%2;     // into timevalue[]
                              return (drs*2 + grp) * 1024; }

// Baseline ("pedestal"): the average level before any signal arrives. We use
// samples 3..52 (skip the first few noisy DRS4 cells; stop before the pulse).
inline float findPedestal(const float* a){
    double ped = 0.0;
    for (int i = 3; i < 53; ++i) ped += a[i];
    return (float)(ped / 50.0);
}

// Peak HEIGHT in mV (positive). `a` points at this slot's 1024 samples:
//   const float* a = &amplitude[0] + ampOffset(slot);
inline float findPeak(const float* a){
    float ped  = findPedestal(a);
    float vmin = a[3];
    for (int i = 4; i < 1024; ++i) if (a[i] < vmin) vmin = a[i];   // most-negative sample
    return ped - vmin;                                             // flip to a positive height
}

// TIME of the peak sample in ns. `a` = slot samples, `t` = that group's time axis:
//   const float* t = &timevalue[0] + timeOffset(slot);
inline float findPeakTime(const float* a, const float* t){
    float vmin = a[3]; int imin = 3;
    for (int i = 4; i < 1024; ++i) if (a[i] < vmin){ vmin = a[i]; imin = i; }
    return t[imin];
}

// ---- The RADiCAL channel map (which slot is what) -------------------------
// 8 high-gain TIMING capillaries (fast, used for the 27 ps result):
static const int HG_SLOT[8] = { 1, 2, 3, 0, 5, 4, 6, 9 };
// 8 low-gain ENERGY capillaries (slow & wide, never clipped):
static const int LG_SLOT[8] = { 19,20,21,18,23,22,24,25 };
static const char* CAP_NAME[8] = { "NW-D","NE-D","SE-D","SW-D","NW-U","NE-U","SE-U","SW-U" };
// Single detectors:
static const int SLOT_1x1  = 15;   // the 1x1 cm trigger scintillator
static const int SLOT_MCP1 =  7;   // micro-channel-plate reference clock 1
static const int SLOT_MCP2 = 16;   // micro-channel-plate reference clock 2
// Wire-chamber planes (beam position), all in DRS1 group 1:
static const int WC_R = 28, WC_L = 29, WC_D = 30, WC_U = 32;
static const double WC_SCALE = 7.0/36.0;   // mm per ns  (left/right & up/down timing -> mm)
// Nominal calorimeter face centre in wire-chamber coordinates [mm]:
static const double CALO_X0 = 6.6, CALO_Y0 = 4.7;

// ============================================================================
// DataPaths.h — THE one place that resolves every RADiCAL data file.
//
// Canonical layout (what lives on CERNBox and what a newcomer downloads):
//     $RAD_DATA/data/<year>/raw/RUN<n>_<E>_GeV.root         (raw 'pulse' tree)
//     $RAD_DATA/data/<year>/reduced/<BUILD>/<E>GeV.root      (analysis 'rad' tree)
//     $RAD_DATA/data/<year>/configs/<BUILD>.json|.hglg       (build configs)
//     $RAD_DATA/data/<year>/metadata/                        (channel map, run table)
//
// $RAD_DATA defaults to the repo root ("."). One canonical path each — no legacy
// fallbacks. If a file is missing, the returned path points at where it SHOULD be,
// so error messages are actionable.
//
// Usage in a macro:
//     #include "DataPaths.h"
//     TFile* f = TFile::Open(radReduced("DSB1", 150));   // DSB1 @150 GeV
//     TFile* g = TFile::Open(radRaw("RUN1258_150_GeV.root"));
//     std::string cfg = radConfig("DSB1").Data();        // build JSON path
// ============================================================================
#ifndef DATAPATHS_H
#define DATAPATHS_H

#include "TString.h"
#include "TSystem.h"

// Default campaign under data/. Override per-shell with `export RAD_YEAR=2024`
// so a new dataset is analyzed without recompiling — the year flows as a
// parameter through every radReduced()/radConfig()/radRaw() call.
inline int radYear(){
    const char* e = gSystem->Getenv("RAD_YEAR");
    if (e && *e){ TString s(e); if (s.IsDigit()){ int y = s.Atoi(); if (y > 2000) return y; } }
    return 2023;
}

inline TString radDataBase(){
    const char* e = gSystem->Getenv("RAD_DATA");
    return (e && *e) ? TString(e) : TString(".");
}

// ---- raw waveform run, by basename e.g. "RUN1258_150_GeV.root" -------------
inline TString radRaw(const char* basename, int year){
    return radDataBase() + Form("/data/%d/raw/%s", year, basename);
}
inline TString radRaw(const char* basename){ return radRaw(basename, radYear()); }

// ---- reduced ntuple for a build + energy ----------------------------------
inline TString radReduced(const char* build, double E, int year){
    return radDataBase() + Form("/data/%d/reduced/%s/%.0fGeV.root", year, build, E);
}
inline TString radReduced(const char* build, double E){ return radReduced(build, E, radYear()); }

// ---- build config JSON + its HG/LG calibration sidecar --------------------
inline TString radConfig(const char* build, int year){
    return radDataBase() + Form("/data/%d/configs/%s.json", year, build);
}
inline TString radConfig(const char* build){ return radConfig(build, radYear()); }

inline TString radHglg(const char* build, int year){
    return radDataBase() + Form("/data/%d/configs/%s.hglg", year, build);
}
inline TString radHglg(const char* build){ return radHglg(build, radYear()); }

#endif // DATAPATHS_H

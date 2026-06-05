// ============================================================================
// DataPaths.h — ONE place that resolves every RADiCAL data file.
//
// Canonical layout (what lives on CERNBox and what a newcomer downloads):
//     $RAD_DATA/datasets/<year>/raw/RUN<n>_<E>_GeV.root        (raw 'pulse' tree)
//     $RAD_DATA/datasets/<year>/reduced/<BUILD>/<E>GeV.root    (analysis 'rad' tree)
//
// $RAD_DATA defaults to the repo root (".").  Each resolver tries the canonical
// path first and FALLS BACK to the legacy in-repo location (Data/, reduced/,
// Analysis/Output/) so existing checkouts keep working during the migration.
//
// Usage in a macro:
//     #include "DataPaths.h"
//     TFile* f = TFile::Open(radReduced("DSB1", 150));   // DSB1 @150 GeV
//     TFile* g = TFile::Open(radRaw("RUN1258_150_GeV.root"));
// ============================================================================
#ifndef DATAPATHS_H
#define DATAPATHS_H

#include "TString.h"
#include "TSystem.h"

static const char* kRadYear = "2023";   // campaign under datasets/

inline TString radDataBase(){
    const char* e = gSystem->Getenv("RAD_DATA");
    return (e && *e) ? TString(e) : TString(".");
}
inline bool radExists(const TString& p){ return gSystem->AccessPathName(p) == kFALSE; }

// return the first path that exists; otherwise the canonical one (so error
// messages point at the intended location).
inline TString radPick(const TString& canonical, const TString& legacy){
    if (radExists(canonical)) return canonical;
    if (radExists(legacy))    return legacy;
    return canonical;
}

// ---- raw waveform run, by basename e.g. "RUN1258_150_GeV.root" -------------
inline TString radRaw(const char* basename){
    TString base = radDataBase();
    return radPick(base + "/datasets/" + kRadYear + "/raw/" + basename,
                   base + "/Data/" + basename);
}

// ---- reduced ntuple for a build + energy ----------------------------------
//   canonical:  datasets/<year>/reduced/<BUILD>/<E>GeV.root
//   legacy:     DSB1 -> Analysis/Output/<E>GeV/ntuple.root
//               others -> reduced/<BUILD>/<E>GeV.root
inline TString radReduced(const char* build, double E){
    TString base = radDataBase();
    TString canonical = base + Form("/datasets/%s/reduced/%s/%.0fGeV.root", kRadYear, build, E);
    TString legacy = (TString(build) == "DSB1")
        ? base + Form("/Analysis/Output/%.0fGeV/ntuple.root", E)
        : base + Form("/reduced/%s/%.0fGeV.root", build, E);
    return radPick(canonical, legacy);
}

#endif // DATAPATHS_H

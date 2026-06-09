// ============================================================================
// FigPaths.h — campaign-namespaced figure paths, so a new dataset's figures
// never overwrite another's. Mirrors DataPaths.h: the year flows as a parameter
// (RAD_YEAR), so the same macro writes figures/2023/... or figures/2024/...
// without a source edit.
//
// Two ways to use it:
//   1. Drop-in wrapper for an existing relative "figures/<sub>/<name>.png" path —
//      just wrap the Print() argument, no other change:
//          c->Print(radFigP("figures/narrative/floor_model_DSB1.png"));
//      ->  figures/2023/narrative/floor_model_DSB1.png   (and mkdir's the dir)
//   2. Build a path from parts:
//          radFig("floor_model", build);   // figures/<year>/narrative/floor_model_<build>.png
//          radFig("optimization");         // figures/<year>/narrative/optimization.png
//          radFig("etype_vs_ttype", "", "capillary");  // figures/<year>/capillary/...
// ============================================================================
#ifndef FIGPATHS_H
#define FIGPATHS_H

#include "DataPaths.h"   // radYear()
#include "TString.h"
#include "TSystem.h"

// Insert the campaign year into a "figures/<sub>/..." path and ensure the dir
// exists. Idempotent: a path already namespaced (figures/<year>/...) is untouched.
// Non-"figures/" paths are returned unchanged (only the dir is created).
inline TString radFigP(const char* relpath){
    TString p(relpath);
    TString tag = Form("figures/%d/", radYear());
    if (p.BeginsWith("figures/") && !p.BeginsWith(tag)){
        p.Replace(0, 8, tag);              // "figures/" -> "figures/<year>/"
    }
    TString dir = gSystem->DirName(p);
    if (dir.Length()) gSystem->mkdir(dir, kTRUE);
    return p;
}

// Build a namespaced figure path from parts: figures/<year>/<sub>/<name>[_<build>].png
inline TString radFig(const char* name, const char* build=nullptr, const char* sub="narrative"){
    TString dir = Form("figures/%d/%s", radYear(), sub);
    gSystem->mkdir(dir, kTRUE);
    if (build && *build) return dir + Form("/%s_%s.png", name, build);
    return dir + Form("/%s.png", name);
}

#endif // FIGPATHS_H

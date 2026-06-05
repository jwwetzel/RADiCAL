// ============================================================================
// RadView.h — format-agnostic, role-resolved view over a reduced "rad" tree
// ----------------------------------------------------------------------------
// The analysis-layer API. Analyses ask for per-capillary-end quantities by
// CANONICAL INDEX (NW-D..SW-U) and RadView returns them whether the file is the
// canonical schema (role-resolved branches present) OR a legacy reduceRaw file
// (only generic s_peak[36]/s_cfd05[36] present), deriving the role-resolved
// values from the BuildConfig channel map in the latter case.
//
// This lets the same analysis run on every build today (DSB1 canonical-ish,
// LUAG/MIXED/TENERGY legacy slots) and unchanged after the cluster re-reduces
// everything to the canonical schema.
//
//   rad::RadView v; v.attach(tree, &cfg);
//   for (Long64_t i=0;i<v.entries();++i){ v.get(i);
//       if (!v.wc_ok() || v.mcp1_peak()<200) continue;
//       double t = v.cfd05(c);  double e = v.sum_lg(); ... }
// ============================================================================
#ifndef RADCORE_RADVIEW_H
#define RADCORE_RADVIEW_H

#include "Schema.h"
#include "BuildConfig.h"
#include "TTree.h"

namespace rad {

struct RadView {
    RadEvent ev;
    const BuildConfig* cfg = nullptr;
    TTree* t = nullptr;
    bool named = false;          // file has role-resolved branches (hg_cfd05[8] etc.)

    void attach(TTree* tree, const BuildConfig* c) {
        t = tree; cfg = c;
        ev.ConnectBranches(tree);
        named = (tree->GetBranch("hg_cfd05") != nullptr);
    }
    Long64_t entries() const { return t ? t->GetEntries() : 0; }
    void     get(Long64_t i) { t->GetEntry(i); }

    // config-invariant passthroughs
    bool   wc_ok()     const { return ev.wc_ok; }
    float  x_trk()     const { return ev.x_trk; }
    float  y_trk()     const { return ev.y_trk; }
    float  mcp1_peak() const { return ev.mcp1_peak; }
    float  mcp2_time() const { return ev.mcp2_time; }
    int    run()       const { return ev.run; }

    // role-resolved per-end accessors (canonical index 0..7)
    float cfd05(int i) const {
        if (named) return ev.hg_cfd05[i];
        float s = ev.s_cfd05[cfg->end[i].hg/1024];
        float r = cfg->end[i].use_mcp2 ? ev.mcp2_time : ev.mcp1_time;
        return (s > -1e5f && r > -1e5f) ? s - r : kNoTime;   // MCP-referenced
    }
    float hg_peak(int i) const { return named ? ev.hg_peak[i] : ev.s_peak[cfg->end[i].hg/1024]; }
    float lg_peak(int i) const { return named ? ev.lg_peak[i] : ev.s_peak[cfg->end[i].lg/1024]; }

    float sum_lg() const {
        if (named) return ev.sum_lg;
        float s = 0; for (int i = 0; i < cfg->nend; ++i) s += ev.s_peak[cfg->end[i].lg/1024]; return s;
    }

    // role of canonical end i (from per-corner config); default "timing"
    bool is_timing(int i) const {
        for (auto& c : cfg->caps) if (c.corner == cfg->end[i].corner) return c.role == "timing";
        return true;
    }
};

} // namespace rad

#endif // RADCORE_RADVIEW_H

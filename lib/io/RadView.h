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

    // selectable per-end timing source (all MCP-referenced in the schema)
    // kLGCFD = branch hg_lgcfd = the paper's "srCFD" (saturation-recovered CFD, PRIMARY
    // for bright/clipped builds); kLED = adopted source for the dim builds (LuAG, TENERGY)
    enum { kCFD03=0, kCFD05, kCFD10, kCFD20, kCFD30, kCFD50, kLED, kLGCFD, kNSrc };
    bool srcAvail[kNSrc] = {false};
    static const char* srcName(int s) {
        static const char* n[kNSrc] = {"cfd03","cfd05","cfd10","cfd20","cfd30","cfd50","led","lgcfd"};
        return (s>=0 && s<kNSrc) ? n[s] : "?";
    }
    bool hasSrc(int s) const { return s>=0 && s<kNSrc && srcAvail[s]; }

    void attach(TTree* tree, const BuildConfig* c) {
        t = tree; cfg = c;
        ev.ConnectBranches(tree);
        named = (tree->GetBranch("hg_cfd05") != nullptr);
        if (named) {
            srcAvail[kCFD03]=tree->GetBranch("hg_cfd03"); srcAvail[kCFD05]=tree->GetBranch("hg_cfd05");
            srcAvail[kCFD10]=tree->GetBranch("hg_cfd10"); srcAvail[kCFD20]=tree->GetBranch("hg_cfd");
            srcAvail[kCFD30]=tree->GetBranch("hg_cfd30"); srcAvail[kCFD50]=tree->GetBranch("hg_cfd50");
            srcAvail[kLED]  =tree->GetBranch("hg_led");   srcAvail[kLGCFD]=tree->GetBranch("hg_lgcfd");
        } else {
            srcAvail[kCFD05] = (tree->GetBranch("s_cfd05") != nullptr);   // legacy: only cfd05 derivable
        }
    }

    // generic MCP-referenced per-end time for source s (kNoTime if unavailable)
    float timeOf(int i, int s) const {
        if (!named) return (s==kCFD05) ? cfd05(i) : kNoTime;
        switch (s) {
            case kCFD03: return ev.hg_cfd03[i]; case kCFD05: return ev.hg_cfd05[i];
            case kCFD10: return ev.hg_cfd10[i]; case kCFD20: return ev.hg_cfd[i];
            case kCFD30: return ev.hg_cfd30[i]; case kCFD50: return ev.hg_cfd50[i];
            case kLED:   return ev.hg_led[i];   case kLGCFD: return ev.hg_lgcfd[i];
        }
        return kNoTime;
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

    // ScanRunCenters-equivalent beam center: ONE sum_lg-light-weighted centroid
    // over the whole file (matches lib/viz/PlotUtils.h ScanRunCenters). This is
    // the center of the timing fiducial circle — the world-class method applied
    // uniformly. Cuts from SelectionCuts.h: mcp>kMCP_minPeak_E(50), sum_lg>
    // kSumLG_centroid(300). Format-agnostic via the accessors, so EVERY build
    // gets it. Does a full pre-pass over the tree.
    void beamCenter(double& xc, double& yc) {
        double wx = 0., wy = 0., w = 0.;
        Long64_t N = entries();
        for (Long64_t i = 0; i < N; ++i) {
            get(i);
            if (!ev.wc_ok || ev.mcp1_peak < 50.0f) continue;   // kMCP_minPeak_E
            double slg = sum_lg();
            if (slg > 300.0) { wx += ev.x_trk * slg; wy += ev.y_trk * slg; w += slg; }  // kSumLG_centroid
        }
        xc = (w > 0.) ? wx / w : 0.;
        yc = (w > 0.) ? wy / w : 0.;
    }
};

} // namespace rad

#endif // RADCORE_RADVIEW_H

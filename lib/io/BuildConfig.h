// ============================================================================
// BuildConfig.h — load a RADiCAL build config (JSON) into resolved DRS4 offsets
// ----------------------------------------------------------------------------
// One JSON file fully describes a build: the DRS4 channel map (which board/
// group/channel is each capillary end, MCP, wire-chamber plane, PbGlass), the
// per-capillary material/role, the module geometry, and the run list. Dropping
// in new data = writing one of these files (no recompile).
//
// The JSON stores hardware coordinates as [drs, grp, ch]; this loader converts
// them to the flat amplitude/time offsets the reducer indexes into, using the
// SAME arithmetic as the legacy Analysis/ChannelConfig.h:
//     chanOff(drs,grp,ch) = (1024*9*2)*drs + (1024*9)*grp + 1024*ch
//     timeOff(drs,grp)    = (1024*2)*drs   + 1024*grp
//
//   rad::BuildConfig cfg = rad::BuildConfig::Load("data/2023/configs/DSB1.json");
//   if (!cfg.valid()) { ... cfg.error ... }
//   for (int i=0;i<cfg.nend;++i) extract(A + cfg.end[i].hg, ...);
// ============================================================================
#ifndef RADCORE_BUILDCONFIG_H
#define RADCORE_BUILDCONFIG_H

#include "MiniJson.h"
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>

namespace rad {

inline int chanOff(int drs, int grp, int ch) { return (1024*9*2)*drs + (1024*9)*grp + 1024*ch; }
inline int timeOff(int drs, int grp)         { return (1024*2)*drs + 1024*grp; }

// resolve a [drs,grp,ch] JSON triple to a flat amplitude offset
inline int offOf(const mj::Value& triple) {
    return chanOff(triple[0].asInt(), triple[1].asInt(), triple[2].asInt());
}
inline int timeOf(const mj::Value& triple) {
    return timeOff(triple[0].asInt(), triple[1].asInt());
}

struct EndMap {
    int  hg = 0, hg_t = 0, lg = 0, lg_t = 0, mcp = 0, mcp_t = 0;
    bool use_mcp2 = false;
    std::string corner, end, name;   // e.g. "NW","D","NW-D"
};

struct CapMap {
    std::string corner, material, role;   // role: "timing" | "energy"
};

struct BuildConfig {
    bool        loaded = false;
    std::string err;
    std::string year, build, description;

    // digitizer geometry
    int    boards = 2, groups_per_board = 2, chans_per_group = 9, samples = 1024;
    double sample_ns = 0.2;
    double lgcfd_frac = 0.15;   // hg_lgcfd (srCFD) threshold = frac * (LG-predicted TRUE HG peak).
                                // WHY 0.15: the steep-but-safe sweet spot — low enough that the
                                // threshold stays below the ~820 mV clip even for the brightest
                                // events (higher fracs re-hit the clip and blow up), high enough
                                // to sit on the steep recovered edge (~x3.7 the foot slope).
                                // Baked in AT REDUCTION, so it is not a post-hoc systematics knob;
                                // the fraction sweep lives in the App. A optimization scan.
    double hg_sat_mV = 950.0;   // HG positive clip/saturation level [mV]. 2023: ~800
                                // because the DT5742 window was DC-offset down to
                                // capture the negative afterpulse (default 950 = full rail)

    // channel map
    int    nend = 0;
    EndMap end[8];
    int    mcp1 = 0, mcp1_t = 0, mcp2 = 0, mcp2_t = 0;
    int    tr0a = 0, tr0a_t = 0, tr0b = 0, tr0b_t = 0;  // TR0 trigger scint, split into ch8 of each group
    double hg_lg_a[8] = {0,0,0,0,0,0,0,0};      // HG_true = hg_lg_a + hg_lg_b*LG_peak (per end)
    double hg_lg_b[8] = {5,5,5,5,5,5,5,5};      // from a <build>.hglg sidecar (calibHGLG.C)
    bool   has_lgcal  = false;                  // true if the sidecar was found -> reducer skips its pre-pass
    int    wc_r = 0, wc_l = 0, wc_d = 0, wc_u = 0, wc_t = 0;
    double wc_scale = 7.0/36.0;
    int    npb = 0, pb[4] = {0,0,0,0}, pb_t = 0;

    // capillaries (per corner)
    std::vector<CapMap> caps;

    // module geometry / physics metadata
    double X0_mm = 5.4, RM_mm = 13.7, depth_X0 = 25.0;
    double module_x0 = 0.0, module_y0 = 0.0;
    double mcp_ref_jitter_ps = 73.0;

    // runs: energy [GeV] -> list of raw file basenames (resolved via DataPaths radRaw)
    std::map<double, std::vector<std::string> > runs;

    bool        valid() const { return loaded && nend > 0; }
    const char* error() const { return err.c_str(); }

    // --------------------------------------------------------------------
    static BuildConfig Load(const std::string& path) {
        BuildConfig c;
        std::string e;
        mj::Value j = mj::parseFile(path, &e);
        if (!e.empty() || j.isNull()) { c.err = e.empty() ? ("empty/invalid config: " + path) : e; return c; }

        c.year        = j["year"].type == mj::Value::Num ? std::to_string(j["year"].asInt()) : j["year"].asStr();
        c.build       = j["build"].asStr();
        c.description = j["description"].asStr();

        const mj::Value& dg = j["digitizer"];
        if (dg.isObj()) {
            c.boards          = dg["boards"].asInt() ? dg["boards"].asInt() : c.boards;
            c.groups_per_board= dg["groups_per_board"].asInt() ? dg["groups_per_board"].asInt() : c.groups_per_board;
            c.chans_per_group = dg["chans_per_group"].asInt() ? dg["chans_per_group"].asInt() : c.chans_per_group;
            c.samples         = dg["samples"].asInt() ? dg["samples"].asInt() : c.samples;
            if (dg.has("sample_ns")) c.sample_ns = dg["sample_ns"].asDouble(c.sample_ns);
            if (dg.has("hg_sat_mV")) c.hg_sat_mV = dg["hg_sat_mV"].asDouble(c.hg_sat_mV);
            if (dg.has("lgcfd_frac")) c.lgcfd_frac = dg["lgcfd_frac"].asDouble(c.lgcfd_frac);
        }

        const mj::Value& cm = j["channel_map"];
        if (!cm.isObj()) { c.err = "missing channel_map in " + path; return c; }

        // MCP references (resolve first — needed by per-end mcp_ref)
        c.mcp1 = offOf(cm["mcp1"]); c.mcp1_t = timeOf(cm["mcp1"]);
        c.mcp2 = offOf(cm["mcp2"]); c.mcp2_t = timeOf(cm["mcp2"]);

        // TR0 trigger scintillator (optional): same PMT pulse split into ch8 of both
        // DRS0 groups -> stored ABSOLUTE for inter-group sync (tr0a - tr0b).
        if (cm.has("tr0a")) { c.tr0a = offOf(cm["tr0a"]); c.tr0a_t = timeOf(cm["tr0a"]); }
        if (cm.has("tr0b")) { c.tr0b = offOf(cm["tr0b"]); c.tr0b_t = timeOf(cm["tr0b"]); }

        // ends (capillary readouts) — canonical order as listed
        const mj::Value& ends = cm["ends"];
        int n = (int)ends.size(); if (n > 8) n = 8;
        for (int i = 0; i < n; ++i) {
            const mj::Value& ed = ends[i];
            EndMap& m = c.end[i];
            m.corner = ed["corner"].asStr();
            m.end    = ed["end"].asStr();
            m.name   = m.corner + "-" + m.end;
            m.hg = offOf(ed["hg"]); m.hg_t = timeOf(ed["hg"]);
            m.lg = offOf(ed["lg"]); m.lg_t = timeOf(ed["lg"]);
            int ref = ed["mcp_ref"].asInt(); if (ref == 0) ref = 1;
            m.use_mcp2 = (ref == 2);
            m.mcp   = m.use_mcp2 ? c.mcp2   : c.mcp1;
            m.mcp_t = m.use_mcp2 ? c.mcp2_t : c.mcp1_t;
        }
        // GUARD (code audit 2026-07): rad::eventDWUP hardcodes canonical index order
        // 0-3 = Down ends, 4-7 = Up ends. A config that lists ends in any other order
        // would silently corrupt the (DW-UP)/2 observable — warn loudly.
        if (n == 8) {
            bool d4u4 = true;
            for (int i = 0; i < 4; ++i) if (c.end[i].end   != "D") d4u4 = false;
            for (int i = 4; i < 8; ++i) if (c.end[i].end   != "U") d4u4 = false;
            if (!d4u4)
                fprintf(stderr, "WARNING [BuildConfig] %s: channel_map.ends is NOT in the "
                        "canonical D,D,D,D,U,U,U,U order that rad::eventDWUP assumes — "
                        "(DW-UP)/2 results from this config would be corrupted!\n", path.c_str());
        }
        c.nend = n;

        // wire chamber
        const mj::Value& wc = cm["wc"];
        c.wc_r = offOf(wc["xr"]); c.wc_l = offOf(wc["xl"]);
        c.wc_d = offOf(wc["yd"]); c.wc_u = offOf(wc["yu"]);
        c.wc_t = timeOf(wc["xr"]);
        if (wc.has("scale_mm_per_ns_frac")) {            // exact fraction [num, den]
            const mj::Value& fr = wc["scale_mm_per_ns_frac"];
            if (fr.isArr() && fr.size() >= 2 && fr[1].asDouble() != 0.0)
                c.wc_scale = fr[0].asDouble() / fr[1].asDouble();
        } else if (wc.has("scale_mm_per_ns")) {
            c.wc_scale = wc["scale_mm_per_ns"].asDouble(c.wc_scale);
        }

        // PbGlass
        const mj::Value& pg = cm["pbglass"];
        int np = (int)pg.size(); if (np > 4) np = 4;
        for (int i = 0; i < np; ++i) c.pb[i] = offOf(pg[i]);
        if (np > 0) c.pb_t = timeOf(pg[0]);
        c.npb = np;

        // capillaries
        const mj::Value& caps = j["capillaries"];
        for (size_t i = 0; i < caps.size(); ++i) {
            CapMap cap;
            cap.corner   = caps[i]["corner"].asStr();
            cap.material = caps[i]["material"].asStr();
            cap.role     = caps[i]["role"].asStr();
            c.caps.push_back(cap);
        }

        // geometry
        const mj::Value& g = j["geometry"];
        if (g.isObj()) {
            c.X0_mm    = g["X0_mm"].asDouble(c.X0_mm);
            c.RM_mm    = g["RM_mm"].asDouble(c.RM_mm);
            c.depth_X0 = g["depth_X0"].asDouble(c.depth_X0);
            const mj::Value& ctr = g["module_center_wc_mm"];
            if (ctr.isArr() && ctr.size() >= 2) { c.module_x0 = ctr[0].asDouble(); c.module_y0 = ctr[1].asDouble(); }
        }
        if (j.has("mcp_ref_jitter_ps")) c.mcp_ref_jitter_ps = j["mcp_ref_jitter_ps"].asDouble(c.mcp_ref_jitter_ps);

        // runs
        const mj::Value& rn = j["runs"];
        if (rn.isObj()) {
            for (auto& kv : rn.obj) {
                double E = std::strtod(kv.first.c_str(), nullptr);
                std::vector<std::string> files;
                for (size_t i = 0; i < kv.second.size(); ++i) files.push_back(kv.second[i].asStr());
                c.runs[E] = files;
            }
        }

        // HG-vs-LG calibration sidecar (<config>.hglg): lines "end_idx a b".
        // Energy-independent per-channel gain ratio used by hg_lgcfd; computed by
        // calibHGLG.C from clean low-energy data (robust vs over-clipped / bad runs).
        { std::string side = path;
          size_t dot = side.rfind(".json"); if (dot != std::string::npos) side = side.substr(0,dot);
          side += ".hglg";
          std::ifstream fin(side.c_str());
          if (fin) { std::string line; int got = 0;
              while (std::getline(fin, line)) {
                  if (line.empty() || line[0]=='#') continue;
                  std::istringstream iss(line); int idx; double a, b;
                  if ((iss >> idx >> a >> b) && idx>=0 && idx<8) { c.hg_lg_a[idx]=a; c.hg_lg_b[idx]=b; ++got; } }
              if (got > 0) c.has_lgcal = true; }
        }

        c.loaded = true;
        return c;
    }
};

} // namespace rad

#endif // RADCORE_BUILDCONFIG_H

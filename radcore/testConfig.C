// ============================================================================
// testConfig.C — Phase-1 gate: JSON config must reproduce the hardcoded map
// ----------------------------------------------------------------------------
// Loads data/2023/configs/DSB1.json via radcore/BuildConfig.h and checks
// that every resolved DRS4 offset EQUALS the legacy Analysis/ChannelConfig.h
// kCap[]/kMCP*/kWC_*/kPbGlass map. If this passes, the declarative config is a
// faithful drop-in replacement for the hardcoded channel map.
//
//   ROOT_INCLUDE_PATH=radcore:Analysis root -l -b -q radcore/testConfig.C+
// ============================================================================
#include "BuildConfig.h"        // radcore
#include "ChannelConfig.h"      // Analysis (legacy hardcoded map)
#include <cstdio>

static int g_fail = 0;
static void CK(const char* what, long got, long exp) {
    bool ok = (got == exp);
    if (!ok) ++g_fail;
    printf("  [%s] %-16s got=%-8ld exp=%-8ld\n", ok ? "OK" : "FAIL", what, got, exp);
}

void testConfig() {
    rad::BuildConfig cfg = rad::BuildConfig::Load("data/2023/configs/DSB1.json");
    if (!cfg.valid()) { printf("CONFIG LOAD FAILED: %s\n", cfg.error()); return; }

    printf("Loaded build '%s' (year %s): %s\n", cfg.build.c_str(), cfg.year.c_str(), cfg.description.c_str());
    printf("nend=%d  caps=%zu  energies=%zu\n\n", cfg.nend, cfg.caps.size(), cfg.runs.size());

    printf("--- per-end channel map vs kCap[] ---\n");
    for (int i = 0; i < kNCap; ++i) {
        printf(" end %d (%s):\n", i, cfg.end[i].name.c_str());
        CK("hg",       cfg.end[i].hg,       kCap[i].hg);
        CK("lg",       cfg.end[i].lg,       kCap[i].lg);
        CK("hg_t",     cfg.end[i].hg_t,     kCap[i].hg_t);
        CK("lg_t",     cfg.end[i].lg_t,     kCap[i].lg_t);
        CK("mcp",      cfg.end[i].mcp,      kCap[i].mcp);
        CK("mcp_t",    cfg.end[i].mcp_t,    kCap[i].mcp_t);
        CK("use_mcp2", cfg.end[i].use_mcp2, kCap[i].use_mcp2);
    }

    printf("--- references ---\n");
    CK("mcp1", cfg.mcp1, kMCP1);   CK("mcp2", cfg.mcp2, kMCP2);
    CK("wc_r", cfg.wc_r, kWC_R);   CK("wc_l", cfg.wc_l, kWC_L);
    CK("wc_d", cfg.wc_d, kWC_D);   CK("wc_u", cfg.wc_u, kWC_U);
    for (int i = 0; i < 4; ++i) CK(Form("pb%d", i), cfg.pb[i], kPbGlass[i]);
    CK("wc_scale_eq", (long)(cfg.wc_scale * 1e9 + 0.5), (long)(kWC_Scale * 1e9 + 0.5));

    printf("--- geometry / runs ---\n");
    printf("  X0=%.1f mm  RM=%.1f mm  depth=%.0f X0  center=(%.1f,%.1f)  jitter=%.0f ps\n",
           cfg.X0_mm, cfg.RM_mm, cfg.depth_X0, cfg.module_x0, cfg.module_y0, cfg.mcp_ref_jitter_ps);
    for (auto& kv : cfg.runs)
        printf("  %.0f GeV: %zu file(s)  (e.g. %s)\n", kv.first, kv.second.size(),
               kv.second.empty() ? "-" : kv.second[0].c_str());

    printf("\n==== %s ====\n", g_fail == 0 ? "ALL CHECKS PASSED — config == hardcoded map" :
                                              Form("%d CHECK(S) FAILED", g_fail));
}

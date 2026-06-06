// ============================================================================
// calibHGLG.C — HG-vs-LG gain calibration for hg_lgcfd (writes a <build>.hglg sidecar).
// ----------------------------------------------------------------------------
// The per-channel slope HG_peak = a + b*LG_peak is an energy-INDEPENDENT gain ratio.
// It must be fit where HG is UNCLIPPED -> use the LOWEST energies (25/50 GeV). High-E
// files are fully clipped (no lever arm) and some single runs are degenerate
// (e.g. RUN1261 slope=0), so a per-file fit is unreliable; this pools clean low-E
// runs with robust (3-sigma trimmed) fits and writes the result for the reducer.
//
// Raw files come from the config's "runs" block (DSB1) OR, for manifest-driven
// builds (LUAG/MIXED/TENERGY), from the SGE tasklist that submit_reduce.sh builds
// (tab-separated: run, label, energy, raw_path) passed as the 2nd argument.
//
//   # DSB1 (has a runs block):
//   root -l -b -q 'radcore/calibHGLG.C+("datasets/2023/configs/DSB1.json")'
//   # LUAG/MIXED/TENERGY (tasklist from your reduce run):
//   root -l -b -q 'radcore/calibHGLG.C+("datasets/2023/configs/LUAG.json","'$RAD_WORK'/tasks_LUAG.txt")'
//   -> writes datasets/2023/configs/<BUILD>.hglg  (read automatically by BuildConfig)
// ============================================================================
#include "BuildConfig.h"
#include "WaveformUtils.h"
#include "DataPaths.h"
#include "TChain.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>

void calibHGLG(const char* configPath, const char* tasklist = "", int nLowE = 2, long maxEvtPerE = 120000){
    rad::BuildConfig cfg = rad::BuildConfig::Load(configPath);
    if (!cfg.valid()) { printf("config load failed: %s\n", cfg.error()); return; }

    // energy -> list of FULL raw paths, from the tasklist (preferred) or cfg.runs
    std::map<double, std::vector<std::string>> runs;
    if (tasklist && tasklist[0]) {
        std::ifstream tf(tasklist);
        if (!tf) { printf("tasklist not found: %s\n", tasklist); return; }
        std::string line;
        while (std::getline(tf, line)) {
            if (line.empty() || line[0]=='#') continue;
            std::vector<std::string> f; std::string cur;
            for (char ch : line) { if (ch=='\t') { f.push_back(cur); cur.clear(); } else cur.push_back(ch); }
            f.push_back(cur);
            if (f.size() < 4) continue;                       // run, label, energy, raw_path
            double E = std::strtod(f[2].c_str(), nullptr);
            if (E > 0) runs[E].push_back(f[3]);
        }
    } else {
        for (auto& kv : cfg.runs) runs[kv.first] = std::vector<std::string>();
        for (auto& kv : cfg.runs) for (auto& b : kv.second) runs[kv.first].push_back(radRaw(b.c_str()).Data());
    }
    std::vector<double> Es; for (auto& kv : runs) Es.push_back(kv.first);
    std::sort(Es.begin(), Es.end());
    if (Es.empty()) { printf("no runs found (give a tasklist for manifest-driven builds)\n"); return; }
    if ((int)Es.size() > nLowE) Es.resize(nLowE);
    printf("[calibHGLG] %s : fitting HG=a+b*LG from energies", cfg.build.c_str());
    for (double E : Es) printf(" %.0f", E); printf(" GeV\n");

    // collect unclipped (LG,HG) pairs per end
    std::vector<std::vector<std::pair<float,float>>> pts(8);
    for (double E : Es) {
        TChain ch("pulse");
        for (auto& path : runs[E]) ch.Add(path.c_str());
        TTreeReader rd(&ch);
        TTreeReaderArray<float> amp(rd,"amplitude"), tim(rd,"timevalue");
        long cnt = 0;
        while (rd.Next() && cnt < maxEvtPerE) { ++cnt;
            const float* T=&tim[0]; const float* A=&amp[0];
            for (int i = 0; i < cfg.nend; ++i) { const rad::EndMap& c = cfg.end[i];
                Pulse hg = ExtractPulse(T + c.hg_t, A + c.hg, 0.20f, 5.f);
                Pulse lg = ExtractPulse(T + c.lg_t, A + c.lg, 0.20f, 5.f);
                if (hg.peak > 30.f && hg.peak < 700.f && lg.peak > 10.f)
                    pts[i].push_back({lg.peak, hg.peak}); }
        }
    }

    // robust per-end linear fit: fit, then 3-sigma trim on residuals, refit
    double A8[8]={0}, B8[8]={5,5,5,5,5,5,5,5};
    for (int i = 0; i < cfg.nend; ++i) {
        auto fit = [&](const std::vector<std::pair<float,float>>& v, double& a, double& b)->long{
            double sx=0,sy=0,sxx=0,sxy=0; long n=v.size(); if(n<50) return 0;
            for (auto& p : v){ sx+=p.first; sy+=p.second; sxx+=p.first*p.first; sxy+=p.first*p.second; }
            double d = n*sxx - sx*sx; if (std::fabs(d)<1e-9) return 0;
            b = (n*sxy - sx*sy)/d; a = (sy - b*sx)/n; return n; };
        double a=0,b=5; long n0=fit(pts[i],a,b);
        if (n0>0){ double s2=0; for(auto&p:pts[i]){double r=p.second-(a+b*p.first); s2+=r*r;} double rms=std::sqrt(s2/n0);
            std::vector<std::pair<float,float>> kept; for(auto&p:pts[i]) if(std::fabs(p.second-(a+b*p.first))<3*rms) kept.push_back(p);
            fit(kept,a,b); n0=kept.size(); }
        A8[i]=a; B8[i]=b;
        const char* flag = (B8[i]>1.0 && B8[i]<8.0 && n0>500) ? "" : "  <-- SUSPECT (check)";
        printf("  %-5s HG = %6.1f + %.3f*LG   (n=%ld)%s\n", cfg.end[i].name.c_str(), A8[i], B8[i], n0, flag);
    }

    // write sidecar  <config without .json>.hglg
    std::string side = configPath; size_t dot = side.rfind(".json");
    if (dot != std::string::npos) side = side.substr(0,dot); side += ".hglg";
    std::ofstream out(side.c_str());
    out << "# HG_true = a + b*LG_peak per end (calibHGLG.C). cols: end_idx a b\n";
    for (int i = 0; i < cfg.nend; ++i) out << i << " " << A8[i] << " " << B8[i] << "\n";
    out.close();
    printf("[calibHGLG] wrote %s\n", side.c_str());
}

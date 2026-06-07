// ============================================================================
// calibHGLG.C — HG-vs-LG gain calibration for hg_lgcfd (writes a <build>.hglg sidecar).
// ----------------------------------------------------------------------------
// The per-channel slope HG_peak = a + b*LG_peak is an energy-INDEPENDENT gain ratio,
// but it is only cleanly fit at energies where the HG-LG diagonal is well-populated and
// not over-clipped. Which energies are clean varies by BUILD: DSB1 (high-light LYSO) is
// clean at 25-125 but degenerate at 150 (94% clipped -> biased unclipped tail); MIXED is
// degenerate at 50/75 (faint/bad runs) but clean at 100-150. So: fit EACH energy
// separately (spike cut + 3-sigma trim), keep only SANE fits (slope 1.5-8, n>300), and
// take the per-channel MEDIAN. Auto-rejects the bad energies on either end; build-median
// fallback for any channel with no sane energy. Default nLowE=0 -> use ALL energies.
//
// Raw files come from the config's "runs" block (DSB1) OR, for manifest-driven
// builds (LUAG/MIXED/TENERGY), from the SGE tasklist that submit_reduce.sh builds
// (tab-separated: run, label, energy, raw_path) passed as the 2nd argument.
//
//   # DSB1 (has a runs block):
//   root -l -b -q 'reduce/calibHGLG.C+("data/2023/configs/DSB1.json")'
//   # LUAG/MIXED/TENERGY (tasklist from your reduce run):
//   root -l -b -q 'reduce/calibHGLG.C+("data/2023/configs/LUAG.json","'$RAD_WORK'/tasks_LUAG.txt")'
//   -> writes data/2023/configs/<BUILD>.hglg  (read automatically by BuildConfig)
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

void calibHGLG(const char* configPath, const char* tasklist = "", int nLowE = 0, long maxEvtPerE = 120000){
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
    if (nLowE > 0 && (int)Es.size() > nLowE) Es.resize(nLowE);   // nLowE<=0 -> use ALL energies
    printf("[calibHGLG] %s : per-energy HG=a+b*LG over", cfg.build.c_str());
    for (double E : Es) printf(" %.0f", E); printf(" GeV  (median of SANE slopes)\n");

    struct P { float lg, hg, shp; };
    auto fit = [](const std::vector<std::pair<float,float>>& w, double& a, double& b)->long{
        double sx=0,sy=0,sxx=0,sxy=0; long n=w.size(); if(n<30) return 0;
        for (auto& p : w){ sx+=p.first; sy+=p.second; sxx+=p.first*p.first; sxy+=p.first*p.second; }
        double d = n*sxx - sx*sx; if (std::fabs(d)<1e-9) return 0;
        b = (n*sxy - sx*sy)/d; a = (sy - b*sx)/n; return n; };
    // robust per-energy per-channel fit: spike cut (shape<0.6*median) + 3-sigma trim
    auto fitClean = [&](std::vector<P>& raw, double& a, double& b)->long{
        std::vector<float> sh; for (auto& p : raw) sh.push_back(p.shp);
        if (sh.empty()) return 0; std::sort(sh.begin(),sh.end()); double shCut=0.6*sh[sh.size()/2];
        std::vector<std::pair<float,float>> v; for (auto& p : raw) if (p.shp>shCut) v.push_back({p.lg,p.hg});
        long n=fit(v,a,b); if(n<=0) return 0;
        double s2=0; for(auto&p:v){double r=p.second-(a+b*p.first); s2+=r*r;} double rms=std::sqrt(s2/n);
        std::vector<std::pair<float,float>> k; for(auto&p:v) if(std::fabs(p.second-(a+b*p.first))<3*rms) k.push_back(p);
        return fit(k,a,b); };

    // The HG/LG slope is energy-INDEPENDENT (gain ratio). Fit each energy separately,
    // keep only SANE fits (slope 1.5-8, n>300), and take the per-channel MEDIAN. This
    // auto-rejects energies where the diagonal is corrupted/faint (e.g. MIXED's low-E)
    // and uses the clean ones (DSB1: low-E; MIXED: high-E).
    double A8[8]={0}, B8[8]={5,5,5,5,5,5,5,5};
    std::vector<std::vector<double>> slope(8), inter(8);
    for (double E : Es) {
        std::vector<std::vector<P>> pts(8);
        TChain ch("pulse"); for (auto& path : runs[E]) ch.Add(path.c_str());
        TTreeReader rd(&ch); TTreeReaderArray<float> amp(rd,"amplitude"), tim(rd,"timevalue");
        long cnt=0;
        while (rd.Next() && cnt<maxEvtPerE) { ++cnt; const float* T=&tim[0]; const float* A=&amp[0];
            for (int i=0;i<cfg.nend;++i){ const rad::EndMap& c=cfg.end[i];
                Pulse hg=ExtractPulse(T+c.hg_t,A+c.hg,0.20f,5.f), lg=ExtractPulse(T+c.lg_t,A+c.lg,0.20f,5.f);
                if (hg.peak>30.f&&hg.peak<700.f&&lg.peak>10.f) pts[i].push_back({lg.peak,hg.peak,hg.charge/hg.peak}); } }
        printf("  E=%-4.0f:", E);
        for (int i=0;i<cfg.nend;++i){ double a=0,b=0; long n=fitClean(pts[i],a,b);
            bool sane=(b>1.5&&b<8.0&&n>300);
            printf(" %s=%.2f%s", cfg.end[i].name.c_str(), b, sane?"":"*");
            if (sane){ slope[i].push_back(b); inter[i].push_back(a); } }
        printf("   (*=rejected)\n");
    }
    for (int i=0;i<cfg.nend;++i) if (!slope[i].empty()){
        std::sort(slope[i].begin(),slope[i].end()); B8[i]=slope[i][slope[i].size()/2];
        std::sort(inter[i].begin(),inter[i].end()); A8[i]=inter[i][inter[i].size()/2]; }
    // fallback: channels with NO sane energy -> build-median of the well-determined ones
    { std::vector<double> good; for (int i=0;i<cfg.nend;++i) if (!slope[i].empty()) good.push_back(B8[i]);
      if (!good.empty()){ std::sort(good.begin(),good.end()); double medB=good[good.size()/2];
        for (int i=0;i<cfg.nend;++i) if (slope[i].empty()){
            printf("  %-5s FALLBACK: no sane energy -> build-median %.3f\n", cfg.end[i].name.c_str(), medB);
            B8[i]=medB; A8[i]=0.0; } } }
    printf("  RESULT:");
    for (int i=0;i<cfg.nend;++i) printf(" %s=%.0f+%.2fLG", cfg.end[i].name.c_str(), A8[i], B8[i]);
    printf("\n");

    // write sidecar  <config without .json>.hglg
    std::string side = configPath; size_t dot = side.rfind(".json");
    if (dot != std::string::npos) side = side.substr(0,dot); side += ".hglg";
    std::ofstream out(side.c_str());
    out << "# HG_true = a + b*LG_peak per end (calibHGLG.C). cols: end_idx a b\n";
    for (int i = 0; i < cfg.nend; ++i) out << i << " " << A8[i] << " " << B8[i] << "\n";
    out.close();
    printf("[calibHGLG] wrote %s\n", side.c_str());
}

// ============================================================================
// Reducer.C — the ONE config-driven raw -> canonical reduced ntuple maker
// ----------------------------------------------------------------------------
// Replaces the forked processRun.C (DSB1-only, rich) and reduceRaw.C (generic
// slots). Driven entirely by a BuildConfig JSON, it writes the canonical
// rad::RadEvent schema: the full per-capillary timing/energy feature set
// (faithful to processRun.C, so DSB1 numbers reproduce exactly) PLUS the
// generic per-slot arrays (faithful to reduceRaw.C, for rediscovery).
//
//   ROOT_INCLUDE_PATH=radcore:Analysis \
//     root -l -b -q 'radcore/Reducer.C+("datasets/2023/configs/DSB1.json", 150)'
//
//   energy < 0  -> reduce every energy listed in the config's run list.
//   outDir ""   -> datasets/<year>/reduced/<BUILD>/   (canonical location)
// ============================================================================
#include "BuildConfig.h"       // radcore
#include "Schema.h"            // radcore
#include "WaveformUtils.h"     // Analysis: ExtractPulse, ExtractPulseMulti, Pulse, PulseMulti, kNoTime
#include "DRS4Calibration.h"   // Analysis: drs4::FindStopCell
#include "SelectionCuts.h"     // Analysis: kHG_LED_thresh, kWC_minPeak, kFiducial_r_energy

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TSystem.h"
#include <cmath>
#include <cstdio>

// resolve a raw file basename to a path (year-aware canonical, legacy fallback)
static TString radRawPath(const rad::BuildConfig& cfg, const std::string& base) {
    TString p = Form("datasets/%s/raw/%s", cfg.year.c_str(), base.c_str());
    if (!gSystem->AccessPathName(p)) return p;                 // canonical exists
    TString legacy = Form("Data/%s", base.c_str());            // legacy location
    if (!gSystem->AccessPathName(legacy)) return legacy;
    return p;                                                  // canonical (will error if missing)
}

// ---------------------------------------------------------------------------
// Core: reduce an already-built 'pulse' TChain -> outFile (canonical schema).
static long reduceChain(const rad::BuildConfig& cfg, TChain* chain, double energy, const char* outFile) {
    TFile* fout = TFile::Open(outFile, "RECREATE");
    fout->SetCompressionSettings(505);                 // ZSTD level 5 — small, copyable
    TTree* tree = new TTree("rad", Form("RADiCAL canonical reduced - %s %.0f GeV", cfg.build.c_str(), energy));
    tree->SetAutoSave(100000);
    rad::RadEvent ev;
    ev.CreateBranches(tree);
    rad::StampBuild(fout, cfg.build.c_str(), cfg.year.c_str());
    ev.beam_energy = (Float_t)energy;

    // time-axis offsets (config-invariant): D0G0,D0G1,D1G0,D1G1
    const int tD0G0 = rad::timeOff(0,0), tD0G1 = rad::timeOff(0,1);
    const int tD1G0 = rad::timeOff(1,0), tD1G1 = rad::timeOff(1,1);
    const int NS = 1024;

    TTreeReader reader(chain);
    TTreeReaderValue<int>   run_v (reader, "run");
    TTreeReaderValue<int>   evt_v (reader, "event0");
    TTreeReaderArray<float> time_v(reader, "timevalue");
    TTreeReaderArray<float> amp_v (reader, "amplitude");

    long nTot = 0, nWC = 0;
    while (reader.Next()) {
        const float* T = &time_v[0];
        const float* A = &amp_v[0];
        ev.run = *run_v; ev.event = *evt_v;

        // --- DRS4 stop cells ---
        ev.stopcell[0] = drs4::FindStopCell(T + tD0G0);
        ev.stopcell[1] = drs4::FindStopCell(T + tD0G1);
        ev.stopcell[2] = drs4::FindStopCell(T + tD1G0);
        ev.stopcell[3] = drs4::FindStopCell(T + tD1G1);

        // --- wire chamber ---
        Pulse wcR = ExtractPulse(T + cfg.wc_t, A + cfg.wc_r, 0.50f, kWC_minPeak);
        Pulse wcL = ExtractPulse(T + cfg.wc_t, A + cfg.wc_l, 0.50f, kWC_minPeak);
        Pulse wcD = ExtractPulse(T + cfg.wc_t, A + cfg.wc_d, 0.50f, kWC_minPeak);
        Pulse wcU = ExtractPulse(T + cfg.wc_t, A + cfg.wc_u, 0.50f, kWC_minPeak);
        ev.wc_peak[0]=wcR.peak; ev.wc_peak[1]=wcL.peak; ev.wc_peak[2]=wcD.peak; ev.wc_peak[3]=wcU.peak;
        ev.wc_ok = (wcR.valid && wcL.valid && wcD.valid && wcU.valid &&
                    wcR.peakTime>0.f && wcL.peakTime>0.f && wcD.peakTime>0.f && wcU.peakTime>0.f);
        if (ev.wc_ok) {
            ev.x_trk = (Float_t)(cfg.wc_scale * (wcR.peakTime - wcL.peakTime));
            ev.y_trk = (Float_t)(cfg.wc_scale * (wcD.peakTime - wcU.peakTime));
            ++nWC;
        } else { ev.x_trk = kNoTime; ev.y_trk = kNoTime; }
        {
            float dx = ev.x_trk - (Float_t)cfg.module_x0, dy = ev.y_trk - (Float_t)cfg.module_y0;
            ev.in_fiducial = ev.wc_ok && (std::sqrt(dx*dx+dy*dy) < (float)kFiducial_r_energy);
        }

        // --- MCP references ---
        Pulse mcp1 = ExtractPulse(T + cfg.mcp1_t, A + cfg.mcp1, 0.20f, 30.f);
        Pulse mcp2 = ExtractPulse(T + cfg.mcp2_t, A + cfg.mcp2, 0.20f, 30.f);
        ev.mcp1_peak = mcp1.peak; ev.mcp1_time = mcp1.crossingTime;
        ev.mcp2_peak = mcp2.peak; ev.mcp2_time = mcp2.crossingTime;

        // --- TR0 trigger scintillator (ch8 of each DRS0 group), ABSOLUTE leading-edge ---
        // Same PMT pulse split into both groups -> (tr0a_time - tr0b_time) gives the
        // inter-group (mezzanine) timing offset for syncing the two free-running DRS4 chips.
        ev.tr0a_time = kNoTime; ev.tr0a_peak = 0.f; ev.tr0b_time = kNoTime; ev.tr0b_peak = 0.f;
        if (cfg.tr0a) { PulseMulti t0a = ExtractPulseMulti(T + cfg.tr0a_t, A + cfg.tr0a, kHG_LED_thresh, 30.f, 2000.f);
                        ev.tr0a_peak = t0a.peak; ev.tr0a_time = t0a.ledTime; }
        if (cfg.tr0b) { PulseMulti t0b = ExtractPulseMulti(T + cfg.tr0b_t, A + cfg.tr0b, kHG_LED_thresh, 30.f, 2000.f);
                        ev.tr0b_peak = t0b.peak; ev.tr0b_time = t0b.ledTime; }

        // --- per capillary end ---
        ev.sum_lg = 0.f;
        for (int i = 0; i < cfg.nend; ++i) {
            const rad::EndMap& c = cfg.end[i];
            PulseMulti hg = ExtractPulseMulti(T + c.hg_t, A + c.hg, kHG_LED_thresh, 5.f, (float)cfg.hg_sat_mV);
            Pulse      lg = ExtractPulse(T + c.lg_t, A + c.lg, 0.20f, 5.f);

            // HG pedestal RMS (samples 3..52)
            double rms2 = 0.; const float* ha = A + c.hg; float ped = hg.pedestal;
            for (int s = 3; s < 53; ++s) { double d = ha[s] - ped; rms2 += d*d; }
            ev.hg_ped_rms[i] = (Float_t)std::sqrt(rms2 / 50.0);

            ev.hg_saturated[i] = hg.saturated;
            ev.hg_spike[i]     = hg.spike;
            float ref = c.use_mcp2 ? mcp2.crossingTime : mcp1.crossingTime;

            ev.hg_peak[i] = hg.peak; ev.lg_peak[i] = lg.peak; ev.sum_lg += lg.peak;
            ev.hg_cfd[i]   = (hg.cfd20  > -1e5f && ref > -1e5f) ? hg.cfd20  - ref : kNoTime;
            ev.hg_cfd10[i] = (hg.cfd10  > -1e5f && ref > -1e5f) ? hg.cfd10  - ref : kNoTime;
            ev.hg_cfd30[i] = (hg.cfd30  > -1e5f && ref > -1e5f) ? hg.cfd30  - ref : kNoTime;
            ev.hg_cfd50[i] = (hg.cfd50  > -1e5f && ref > -1e5f) ? hg.cfd50  - ref : kNoTime;
            ev.hg_cfd03[i] = (hg.cfd03  > -1e5f && ref > -1e5f) ? hg.cfd03  - ref : kNoTime;
            ev.hg_cfd05[i] = (hg.cfd05  > -1e5f && ref > -1e5f) ? hg.cfd05  - ref : kNoTime;
            ev.hg_led[i]   = (hg.ledTime> -1e5f && ref > -1e5f) ? hg.ledTime- ref : kNoTime;
            ev.hg_tot[i]   = hg.totTime;
            ev.hg_charge[i]= hg.charge;
            ev.lg_charge[i]= lg.charge;
        }

        // --- PbGlass reference ---
        ev.sum_pb = 0.f;
        for (int i = 0; i < cfg.npb; ++i) {
            Pulse pb = ExtractPulse(T + cfg.pb_t, A + cfg.pb[i], 0.20f, 5.f);
            ev.pb_peak[i] = pb.peak; ev.sum_pb += pb.peak;
        }

        // --- generic per-slot features (rediscovery / diagnostics) ---
        for (int s = 0; s < rad::NSLOT; ++s) {
            int drs = s/18, grp = (s/9)%2;
            const float* a  = A + s*NS;
            const float* tt = T + ((drs*2+grp)*NS);
            Pulse p = ExtractPulse(tt, a, 0.05f, 3.f);
            ev.s_peak[s]=p.peak; ev.s_cfd05[s]=p.crossingTime; ev.s_charge[s]=p.charge;
        }

        tree->Fill();
        if (++nTot % 5000 == 0) { printf("  %ld events\r", nTot); fflush(stdout); }
    }
    printf("\n[Reducer] %s %.0f GeV: %ld events (%ld good WC) -> %s\n",
           cfg.build.c_str(), energy, nTot, nWC, outFile);
    fout->cd(); tree->Write(); fout->Close();
    return nTot;
}

// Reduce one (build, energy) using the config's run list -> outFile.
static long reduceOne(const rad::BuildConfig& cfg, double energy, const char* outFile) {
    auto it = cfg.runs.find(energy);
    if (it == cfg.runs.end() || it->second.empty()) {
        printf("[Reducer] %s: no runs for %.0f GeV\n", cfg.build.c_str(), energy);
        return -1;
    }
    TChain* chain = new TChain("pulse");
    int nfiles = 0;
    for (const std::string& base : it->second) {
        TString rp = radRawPath(cfg, base);
        int n = chain->Add(rp);
        printf("[Reducer]   + %s  (%d file)\n", rp.Data(), n);
        nfiles += n;
    }
    if (nfiles == 0) { printf("[Reducer] no raw files for %.0f GeV\n", energy); delete chain; return -1; }
    long n = reduceChain(cfg, chain, energy, outFile);
    delete chain;
    return n;
}

// HPC / per-run entry: reduce a SINGLE raw file with a config -> outFile.
// (SGE array job calls this once per run; merge per energy with hadd afterwards)
void ReduceFile(const char* configPath, const char* rawFile, double energy, const char* outFile) {
    rad::BuildConfig cfg = rad::BuildConfig::Load(configPath);
    if (!cfg.valid()) { printf("[Reducer] config load failed: %s\n", cfg.error()); return; }
    TChain* chain = new TChain("pulse");
    if (chain->Add(rawFile) == 0) { printf("[Reducer] no file: %s\n", rawFile); delete chain; return; }
    reduceChain(cfg, chain, energy, outFile);
    delete chain;
}

// ---------------------------------------------------------------------------
void Reducer(const char* configPath, double energy = -1, const char* outDir = "") {
    rad::BuildConfig cfg = rad::BuildConfig::Load(configPath);
    if (!cfg.valid()) { printf("[Reducer] config load failed: %s\n", cfg.error()); return; }

    TString dir = outDir[0] ? TString(outDir)
                            : Form("datasets/%s/reduced/%s", cfg.year.c_str(), cfg.build.c_str());
    gSystem->mkdir(dir, kTRUE);

    if (energy >= 0) {
        reduceOne(cfg, energy, Form("%s/%.0fGeV.root", dir.Data(), energy));
    } else {
        for (auto& kv : cfg.runs)
            reduceOne(cfg, kv.first, Form("%s/%.0fGeV.root", dir.Data(), kv.first));
    }
}

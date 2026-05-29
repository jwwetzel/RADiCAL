// ============================================================================
// processRun.C — convert raw RADiCAL run file(s) into a compact analysis ntuple
// ============================================================================
//
// Reads flat waveform data from one or more runs (chained) and extracts
// per-event physics quantities into a TTree for fast downstream analysis.
//
// inFiles is a semicolon-separated list of file paths.  A single file or a
// glob pattern (e.g. "Data/RUN125?_150_GeV.root") also work via TChain::Add().
//
// Usage (from the repository root directory):
//   Single run:
//     root -l -b -q 'Analysis/processRun.C+("Data/RUN1034_125_GeV.root", 125., "125GeV")'
//   Chained runs:
//     root -l -b -q 'Analysis/processRun.C+("Data/RUN1258_150_GeV.root;Data/RUN1259_150_GeV.root;Data/RUN1260_150_GeV.root;Data/RUN1261_150_GeV.root", 150., "150GeV")'
//
//   The '+' suffix triggers ACLiC compilation — ~10x faster than interpreted mode.
//
// Output:
//   Analysis/Output/<label>/ntuple.root   containing TTree "rad"
//
// Branches in the output TTree:
//   run           Int_t    run number
//   event         Int_t    event number within run
//   beam_energy   Float_t  nominal beam energy [GeV] (from argument)
//   wc_ok         Bool_t   wire-chamber reconstruction succeeded (all 4 planes valid)
//   x_trk         Float_t  beam x position [mm]  (kNoTime if !wc_ok)
//   y_trk         Float_t  beam y position [mm]  (kNoTime if !wc_ok)
//   wc_peak[4]    Float_t  WC plane amplitudes [mV]: [0]=R [1]=L [2]=D [3]=U
//                          stored regardless of wc_ok — for plane-efficiency monitoring
//   in_fiducial   Bool_t   beam within kFiducial_r_energy of nominal kCalo_x0/y0
//                          NOTE: downstream macros recompute from data-derived centroid
//   mcp_peak      Float_t  MCP1 peak amplitude [mV]
//   mcp_time      Float_t  MCP1 CFD-20% crossing time [ns]  (kNoTime if not found)
//   mcp2_peak     Float_t  MCP2 peak amplitude [mV]
//   mcp2_time     Float_t  MCP2 CFD-20% crossing time [ns]  (kNoTime if not found)
//   hg_peak[8]    Float_t  HG capillary peak amplitudes [mV]
//   hg_ped_rms[8] Float_t  HG pedestal RMS (std-dev of samples 3-52) [mV]
//                          healthy channels: ~1-3 mV; noisy channels: >5 mV
//   hg_saturated[8] Bool_t   true if hg_peak >= 950 mV (DRS4 near-saturation)
//   hg_spike[8]     Bool_t   true if any pedestal sample (3-52) deviates > 5xRMS
//   stopcell[4]     Int_t    DRS4 stop (trigger) cell per group
//                            [0]=D0G0 [1]=D0G1 [2]=D1G0 [3]=D1G1.  Enables the
//                            data-driven stop-cell timing correction downstream
//                            (see DRS4Calibration.h / drs4TimeBase.C).
//   hg_cfd[8]     Float_t  HG CFD-20% time minus MCP reference [ns]
//   hg_cfd10[8]   Float_t  HG CFD-10% time minus MCP reference [ns]
//   hg_cfd30[8]   Float_t  HG CFD-30% time minus MCP reference [ns]
//   hg_cfd50[8]   Float_t  HG CFD-50% time minus MCP reference [ns]
//   hg_led[8]     Float_t  HG fixed-threshold LED time minus MCP reference [ns]
//                          LED absolute threshold = kHG_LED_thresh (20 mV) above pedestal
//   hg_tot[8]     Float_t  HG time over threshold [ns]  (same kHG_LED_thresh threshold)
//                          kNoTime if trailing edge not found in window
//   lg_peak[8]    Float_t  low-gain capillary peak amplitudes [mV]
//   pb_peak[4]    Float_t  PbGlass peak amplitudes [mV]
//   sum_lg        Float_t  sum of all 8 low-gain peaks [mV]  (energy proxy)
//   sum_pb        Float_t  sum of all 4 PbGlass peaks [mV]  (reference)
// ============================================================================

#include "ChannelConfig.h"
#include "WaveformUtils.h"
#include "DRS4Calibration.h"   // FindStopCell — DRS4 stop-cell recovery

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TSystem.h"
#include "TObjArray.h"
#include "TObjString.h"

#include <cmath>
#include <iostream>

// LED/TOT threshold: kHG_LED_thresh from SelectionCuts.h (via ChannelConfig.h).
// 20 mV is well above the ~5 mV HG noise floor and matches the per-channel
// quality cut, so every event used for timing also reaches the LED threshold.

void processRun(const char* inFiles, double beamEnergy_GeV, const char* label)
{
    // -------------------------------------------------------------------------
    // Build TChain — split inFiles on ';', add each segment (file or glob)
    // -------------------------------------------------------------------------
    TChain* chain = new TChain("pulse");
    TString patterns(inFiles);
    TObjArray* parts = patterns.Tokenize(";");
    int nFilesAdded = 0;
    for (int i = 0; i < parts->GetEntries(); ++i) {
        TString pat = ((TObjString*)parts->At(i))->GetString().Strip();
        int n = chain->Add(pat);
        std::cout << "[processRun]   + " << pat << "  (" << n << " file(s))\n";
        nFilesAdded += n;
    }
    delete parts;
    if (nFilesAdded == 0) {
        std::cerr << "[processRun] No files matched: " << inFiles << std::endl;
        delete chain;
        return;
    }
    std::cout << "[processRun] " << label << ": " << nFilesAdded
              << " file(s) chained, " << chain->GetEntries() << " total events\n";

    // -------------------------------------------------------------------------
    // Prepare output
    // -------------------------------------------------------------------------
    TString outDir = TString("Analysis/Output/") + label;
    gSystem->mkdir(outDir, kTRUE);
    TFile* fout = new TFile(outDir + "/ntuple.root", "RECREATE");
    TTree* tree = new TTree("rad", Form("RADiCAL analysis ntuple - %s", label));
    tree->SetAutoSave(100000);  // auto-save every 100k events

    // -------------------------------------------------------------------------
    // Declare output branches
    // -------------------------------------------------------------------------
    Int_t   br_run, br_event;
    Bool_t  br_wc_ok, br_in_fiducial;
    Float_t br_x_trk, br_y_trk;
    Float_t br_mcp_peak,  br_mcp_time;
    Float_t br_mcp2_peak, br_mcp2_time;
    Float_t br_hg_peak[8];
    Float_t br_hg_ped_rms[8]; // HG pedestal RMS (baseline noise) per channel [mV]
    Bool_t  br_hg_saturated[8];  // HG peak near DRS4 saturation range
    Bool_t  br_hg_spike[8];      // spike in HG pedestal window
    Int_t   br_stopcell[4];      // DRS4 stop cell per group (D0G0,D0G1,D1G0,D1G1)
    Float_t br_hg_cfd03[8];      // HG CFD-3% time minus MCP ref [ns] (Ledovskoy sim optimum)
    Float_t br_hg_cfd05[8];      // HG CFD-5% time minus MCP ref [ns]
    Float_t br_hg_charge[8];     // HG integrated charge (energy proxy) [mV*sample]
    Float_t br_lg_charge[8];     // LG integrated charge (energy proxy) [mV*sample]
    Float_t br_hg_cfd[8];    // CFD 20% — backwards-compatible
    Float_t br_hg_cfd10[8];  // CFD 10%
    Float_t br_hg_cfd30[8];  // CFD 30%
    Float_t br_hg_cfd50[8];  // CFD 50%
    Float_t br_hg_led[8];    // LED fixed-threshold (20 mV)
    Float_t br_hg_tot[8];    // time over threshold [ns]
    Float_t br_lg_peak[8];
    Float_t br_wc_peak[4];   // individual WC plane peak amplitudes [mV]
    Float_t br_pb_peak[4];
    Float_t br_sum_lg, br_sum_pb;
    Float_t br_beam_energy;

    tree->Branch("run",          &br_run,          "run/I");
    tree->Branch("event",        &br_event,        "event/I");
    tree->Branch("beam_energy",  &br_beam_energy,  "beam_energy/F");
    tree->Branch("wc_ok",        &br_wc_ok,        "wc_ok/O");
    tree->Branch("x_trk",        &br_x_trk,        "x_trk/F");
    tree->Branch("y_trk",        &br_y_trk,        "y_trk/F");
    tree->Branch("in_fiducial",  &br_in_fiducial,  "in_fiducial/O");
    tree->Branch("mcp_peak",     &br_mcp_peak,     "mcp_peak/F");
    tree->Branch("mcp_time",     &br_mcp_time,     "mcp_time/F");
    tree->Branch("mcp2_peak",    &br_mcp2_peak,    "mcp2_peak/F");
    tree->Branch("mcp2_time",    &br_mcp2_time,    "mcp2_time/F");
    tree->Branch("hg_peak",       br_hg_peak,      "hg_peak[8]/F");
    tree->Branch("hg_ped_rms",    br_hg_ped_rms,   "hg_ped_rms[8]/F");
    tree->Branch("hg_saturated",  br_hg_saturated,  "hg_saturated[8]/O");
    tree->Branch("hg_spike",       br_hg_spike,       "hg_spike[8]/O");
    tree->Branch("stopcell",       br_stopcell,       "stopcell[4]/I");
    tree->Branch("hg_cfd03",       br_hg_cfd03,      "hg_cfd03[8]/F");
    tree->Branch("hg_cfd05",       br_hg_cfd05,      "hg_cfd05[8]/F");
    tree->Branch("hg_charge",      br_hg_charge,     "hg_charge[8]/F");
    tree->Branch("lg_charge",      br_lg_charge,     "lg_charge[8]/F");
    tree->Branch("hg_cfd",        br_hg_cfd,       "hg_cfd[8]/F");
    tree->Branch("hg_cfd10",      br_hg_cfd10,     "hg_cfd10[8]/F");
    tree->Branch("hg_cfd30",      br_hg_cfd30,     "hg_cfd30[8]/F");
    tree->Branch("hg_cfd50",      br_hg_cfd50,     "hg_cfd50[8]/F");
    tree->Branch("hg_led",        br_hg_led,       "hg_led[8]/F");
    tree->Branch("hg_tot",        br_hg_tot,       "hg_tot[8]/F");
    tree->Branch("lg_peak",       br_lg_peak,      "lg_peak[8]/F");
    tree->Branch("wc_peak",       br_wc_peak,      "wc_peak[4]/F");
    tree->Branch("pb_peak",       br_pb_peak,      "pb_peak[4]/F");
    tree->Branch("sum_lg",       &br_sum_lg,       "sum_lg/F");
    tree->Branch("sum_pb",       &br_sum_pb,       "sum_pb/F");

    br_beam_energy = static_cast<Float_t>(beamEnergy_GeV);

    // -------------------------------------------------------------------------
    // Event loop
    // -------------------------------------------------------------------------
    TTreeReader reader(chain);
    TTreeReaderValue<int>   run_v  (reader, "run");
    TTreeReaderValue<int>   evt_v  (reader, "event0");
    TTreeReaderArray<float> time_v (reader, "timevalue");
    TTreeReaderArray<float> amp_v  (reader, "amplitude");

    int nTotal = 0, nGoodWC = 0, nFiducial = 0;

    while (reader.Next()) {
        const float* T = &time_v[0];
        const float* A = &amp_v[0];

        br_run   = *run_v;
        br_event = *evt_v;

        // --- DRS4 stop cell per group --------------------------------------
        // Recovered from the single zero-width step in each group's nominal
        // time axis.  Stored so downstream timing macros can apply the
        // data-driven stop-cell correction (see DRS4Calibration.h).
        br_stopcell[0] = drs4::FindStopCell(T + kT_D0G0);
        br_stopcell[1] = drs4::FindStopCell(T + kT_D0G1);
        br_stopcell[2] = drs4::FindStopCell(T + kT_D1G0);
        br_stopcell[3] = drs4::FindStopCell(T + kT_D1G1);

        // --- Wire chamber reconstruction -----------------------------------
        Pulse wcR = ExtractPulse(T + kWC_t, A + kWC_R, 0.50f, kWC_minPeak);
        Pulse wcL = ExtractPulse(T + kWC_t, A + kWC_L, 0.50f, kWC_minPeak);
        Pulse wcD = ExtractPulse(T + kWC_t, A + kWC_D, 0.50f, kWC_minPeak);
        Pulse wcU = ExtractPulse(T + kWC_t, A + kWC_U, 0.50f, kWC_minPeak);

        // Individual WC plane amplitudes — for plane-efficiency monitoring
        br_wc_peak[0] = wcR.peak;
        br_wc_peak[1] = wcL.peak;
        br_wc_peak[2] = wcD.peak;
        br_wc_peak[3] = wcU.peak;

        br_wc_ok = (wcR.valid && wcL.valid && wcD.valid && wcU.valid &&
                    wcR.peakTime > 0.f && wcL.peakTime > 0.f &&
                    wcD.peakTime > 0.f && wcU.peakTime > 0.f);

        if (br_wc_ok) {
            br_x_trk = static_cast<Float_t>(kWC_Scale * (wcR.peakTime - wcL.peakTime));
            br_y_trk = static_cast<Float_t>(kWC_Scale * (wcD.peakTime - wcU.peakTime));
            ++nGoodWC;
        } else {
            br_x_trk = kNoTime;
            br_y_trk = kNoTime;
        }

        // in_fiducial uses the nominal calorimeter center (kCalo_x0/y0) and the
        // energy fiducial radius.  Downstream timing macros recompute their own
        // fiducial from the data-derived centroid (ScanRunCenters), so this stored
        // branch is used only for a quick first-pass energy-analysis selection.
        float dx = br_x_trk - static_cast<Float_t>(kCalo_x0);
        float dy = br_y_trk - static_cast<Float_t>(kCalo_y0);
        br_in_fiducial = br_wc_ok &&
                         (std::sqrt(dx*dx + dy*dy) < static_cast<float>(kFiducial_r_energy));
        if (br_in_fiducial) ++nFiducial;

        // --- MCP reference timing ------------------------------------------
        // MCP1 is used for capillaries 0-6; MCP2 for capillary 7 (SW-U).
        // Store both MCPs so timingMethods.C can apply MCP amplitude walk correction.
        Pulse mcp1 = ExtractPulse(T + kMCP1_t, A + kMCP1, 0.20f, 30.f);
        Pulse mcp2 = ExtractPulse(T + kMCP2_t, A + kMCP2, 0.20f, 30.f);
        br_mcp_peak  = mcp1.peak;
        br_mcp_time  = mcp1.crossingTime;
        br_mcp2_peak = mcp2.peak;
        br_mcp2_time = mcp2.crossingTime;

        // --- 8 capillary pairs ---------------------------------------------
        // ExtractPulseMulti gives CFD at four fractions, LED, and TOT in one
        // scan — no extra overhead vs calling ExtractPulse four times.
        br_sum_lg = 0.f;
        for (int i = 0; i < kNCap; ++i) {
            const CapCfg& c = kCap[i];

            PulseMulti hg = ExtractPulseMulti(T + c.hg_t, A + c.hg,
                                               kHG_LED_thresh, 5.f, 950.0f);
            Pulse lg      = ExtractPulse(T + c.lg_t, A + c.lg, 0.20f, 5.f);

            // HG pedestal RMS: std-dev of samples 3..52 around the pedestal mean.
            // A healthy channel has RMS ~1-3 mV.  Values > ~5 mV indicate a noisy
            // DRS4 capacitor cell or ADC instability and should be flagged in
            // qualityPlots.C.  Computed after ExtractPulseMulti so hg.pedestal exists.
            {
                double ped_rms2 = 0.;
                const float* ha = A + c.hg;
                float         ped = hg.pedestal;
                for (int s = 3; s < 53; ++s) {
                    double dev = ha[s] - ped;
                    ped_rms2 += dev * dev;
                }
                br_hg_ped_rms[i] = static_cast<Float_t>(std::sqrt(ped_rms2 / 50.0));
            }

            br_hg_saturated[i] = hg.saturated;
            br_hg_spike[i]     = hg.spike;

            // Choose the correct MCP reference for this capillary
            float mcpRef_t20 = c.use_mcp2 ? mcp2.crossingTime : mcp1.crossingTime;

            br_hg_peak[i] = hg.peak;
            br_lg_peak[i] = lg.peak;
            br_sum_lg    += lg.peak;

            // Helper: compute relative time; store kNoTime if either is missing
            // CFD-20% (backwards compatible — hg.cfd20 == ExtractPulse with 0.20)
            br_hg_cfd[i]   = (hg.cfd20  > -1e5f && mcpRef_t20 > -1e5f)
                               ? hg.cfd20  - mcpRef_t20 : kNoTime;
            br_hg_cfd10[i] = (hg.cfd10  > -1e5f && mcpRef_t20 > -1e5f)
                               ? hg.cfd10  - mcpRef_t20 : kNoTime;
            br_hg_cfd30[i] = (hg.cfd30  > -1e5f && mcpRef_t20 > -1e5f)
                               ? hg.cfd30  - mcpRef_t20 : kNoTime;
            br_hg_cfd50[i] = (hg.cfd50  > -1e5f && mcpRef_t20 > -1e5f)
                               ? hg.cfd50  - mcpRef_t20 : kNoTime;
            br_hg_cfd03[i] = (hg.cfd03  > -1e5f && mcpRef_t20 > -1e5f)
                               ? hg.cfd03  - mcpRef_t20 : kNoTime;
            br_hg_cfd05[i] = (hg.cfd05  > -1e5f && mcpRef_t20 > -1e5f)
                               ? hg.cfd05  - mcpRef_t20 : kNoTime;
            br_hg_led[i]   = (hg.ledTime > -1e5f && mcpRef_t20 > -1e5f)
                               ? hg.ledTime - mcpRef_t20 : kNoTime;
            // TOT is a duration (not relative to MCP), stored directly
            br_hg_tot[i]   = hg.totTime;
            // Integrated charge (energy proxy)
            br_hg_charge[i] = hg.charge;
            br_lg_charge[i] = lg.charge;
        }

        // --- PbGlass reference ---------------------------------------------
        br_sum_pb = 0.f;
        for (int i = 0; i < 4; ++i) {
            Pulse pb      = ExtractPulse(T + kPbGlass_t, A + kPbGlass[i], 0.20f, 5.f);
            br_pb_peak[i] = pb.peak;
            br_sum_pb    += pb.peak;
        }

        tree->Fill();
        ++nTotal;

        if (nTotal % 5000 == 0)
            std::cout << "  " << nTotal << " events processed\r" << std::flush;
    }

    std::cout << "\n[processRun] " << label
              << ": " << nTotal    << " total events"
              << ", " << nGoodWC   << " with good WC"
              << ", " << nFiducial << " in fiducial\n";
    std::cout << "Output: " << outDir << "/ntuple.root\n";

    fout->Write();
    fout->Close();
    delete chain;
}

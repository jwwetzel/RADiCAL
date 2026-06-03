// ============================================================================
// reduceRaw.C — config-AGNOSTIC raw->reduced ntuple maker (for any capillary
// configuration: DSB1, LuAG, mixed, timing+energy).
// ----------------------------------------------------------------------------
// The capillary FILL changes between configs, but the wire chamber, MCPs and
// the DRS4 slot/time-axis scheme do NOT.  So this reducer does NOT assume which
// slot is which capillary.  For EVERY one of the 36 DRS slots it stores the
// pulse features (peak, peak-time, CFD-5% time, charge); it also precomputes the
// config-invariant WC track position and both MCP references for convenience.
//
// => One reduced file works for any config.  Role assignment (which slot is a
//    timing vs energy capillary, corner/depth geometry) and the timing/energy
//    resolution are done LOCALLY on these small files, per config, with
//    discoverChannels.C + the analysis chain.
//
// Per-event size ~0.6 kB; a merged per-energy file is tens of MB (copyable).
//
//   root -l 'Analysis/reduceRaw.C+("RUN1234.root", 150., "out.root")'
// ============================================================================
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "ChannelConfig.h"     // chanOff/timeOff, kWC_*, kMCP*, kT_* (config-invariant)
#include "WaveformUtils.h"     // ExtractPulse, Pulse, kNoTime
#include "DRS4Calibration.h"   // drs4::FindStopCell
#include "SelectionCuts.h"     // kWC_minPeak, kWC_Scale via ChannelConfig
#include <cstdio>

void reduceRaw(const char* rawFile, double beamEnergy_GeV, const char* outFile)
{
    TFile* fin = TFile::Open(rawFile);
    if(!fin || fin->IsZombie()){ printf("[reduceRaw] cannot open %s\n", rawFile); return; }
    TTree* tin = (TTree*)fin->Get("pulse");
    if(!tin){ printf("[reduceRaw] no 'pulse' tree in %s\n", rawFile); return; }

    const int NSLOT = 36, NS = 1024;

    TFile* fout = TFile::Open(outFile, "RECREATE");
    fout->SetCompressionSettings(505);   // ZSTD (algo 5) level 5 — strong+fast, small copies
    TTree* rad  = new TTree("rad", "config-agnostic reduced ntuple");

    Int_t   run, event; Float_t energy;
    Bool_t  wc_ok; Float_t x_trk, y_trk;
    Float_t mcp1_peak, mcp1_time, mcp2_peak, mcp2_time;
    Int_t   stopcell[4];
    Float_t s_peak[NSLOT], s_cfd05[NSLOT], s_charge[NSLOT];

    rad->Branch("run",&run); rad->Branch("event",&event); rad->Branch("beam_energy",&energy);
    rad->Branch("wc_ok",&wc_ok); rad->Branch("x_trk",&x_trk); rad->Branch("y_trk",&y_trk);
    rad->Branch("mcp1_peak",&mcp1_peak); rad->Branch("mcp1_time",&mcp1_time);
    rad->Branch("mcp2_peak",&mcp2_peak); rad->Branch("mcp2_time",&mcp2_time);
    rad->Branch("stopcell",stopcell,"stopcell[4]/I");
    // per-slot features (slot s: drs=s/18, grp=(s/9)%2, ch=s%9).  WC track (x_trk,
    // y_trk) is already reconstructed above, so per-slot peak-time isn't stored.
    rad->Branch("s_peak",  s_peak,  Form("s_peak[%d]/F",NSLOT));    // amplitude [mV]
    rad->Branch("s_cfd05", s_cfd05, Form("s_cfd05[%d]/F",NSLOT));   // CFD-5% crossing time [ns]
    rad->Branch("s_charge",s_charge,Form("s_charge[%d]/F",NSLOT));  // pulse integral

    TTreeReader reader(tin);
    TTreeReaderValue<int>   run_v (reader, "run");
    TTreeReaderValue<int>   evt_v (reader, "event0");
    TTreeReaderArray<float> time_v(reader, "timevalue");
    TTreeReaderArray<float> amp_v (reader, "amplitude");

    long n = 0;
    while(reader.Next()){
        const float* T = &time_v[0];
        const float* A = &amp_v[0];
        run = *run_v; event = *evt_v; energy = (Float_t)beamEnergy_GeV;

        stopcell[0]=drs4::FindStopCell(T+kT_D0G0); stopcell[1]=drs4::FindStopCell(T+kT_D0G1);
        stopcell[2]=drs4::FindStopCell(T+kT_D1G0); stopcell[3]=drs4::FindStopCell(T+kT_D1G1);

        // config-invariant WC + MCP
        Pulse wcR=ExtractPulse(T+kWC_t,A+kWC_R,0.50f,kWC_minPeak);
        Pulse wcL=ExtractPulse(T+kWC_t,A+kWC_L,0.50f,kWC_minPeak);
        Pulse wcD=ExtractPulse(T+kWC_t,A+kWC_D,0.50f,kWC_minPeak);
        Pulse wcU=ExtractPulse(T+kWC_t,A+kWC_U,0.50f,kWC_minPeak);
        wc_ok = (wcR.valid&&wcL.valid&&wcD.valid&&wcU.valid&&
                 wcR.peakTime>0&&wcL.peakTime>0&&wcD.peakTime>0&&wcU.peakTime>0);
        x_trk = wc_ok ? (Float_t)(kWC_Scale*(wcR.peakTime-wcL.peakTime)) : kNoTime;
        y_trk = wc_ok ? (Float_t)(kWC_Scale*(wcD.peakTime-wcU.peakTime)) : kNoTime;

        Pulse m1=ExtractPulse(T+kMCP1_t,A+kMCP1,0.20f,30.f);
        Pulse m2=ExtractPulse(T+kMCP2_t,A+kMCP2,0.20f,30.f);
        mcp1_peak=m1.peak; mcp1_time=m1.crossingTime; mcp2_peak=m2.peak; mcp2_time=m2.crossingTime;

        // every slot, generically
        for(int s=0;s<NSLOT;++s){
            int drs=s/18, grp=(s/9)%2;
            const float* a = A + s*NS;
            const float* tt= T + ((drs*2+grp)*NS);
            Pulse p = ExtractPulse(tt, a, 0.05f, 3.f);
            s_peak[s]=p.peak; s_cfd05[s]=p.crossingTime; s_charge[s]=p.charge;
        }
        rad->Fill(); ++n;
    }
    fin->Close();
    fout->cd(); rad->Write(); fout->Close();
    printf("[reduceRaw] %s -> %s : %ld events @ %.0f GeV\n", rawFile, outFile, n, beamEnergy_GeV);
}

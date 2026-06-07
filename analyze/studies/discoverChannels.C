// ============================================================================
// discoverChannels.C — infer the slot->role channel map from a RAW run file.
// ----------------------------------------------------------------------------
// The raw `pulse` tree stores amplitude[36*1024] (36 slots = 2 DRS x 2 groups
// x 9 ch) and timevalue[4*1024] (one time axis per DRS-group).  Nothing in the
// file says what each slot IS — that mapping is hand-coded in ChannelConfig.h
// for DSB1 only.  This tool characterises every slot (occupancy, amplitude,
// edge speed, polarity) so a NEW configuration's map (LuAG, ...) can be derived
// from data instead of guessed.
//
// Validation: run on a DSB1 file; the inferred roles must match ChannelConfig.h.
//
//   root -l 'Analysis/discoverChannels.C+("Data/RUN1075_100_GeV.root")'
// ============================================================================
#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "WaveformUtils.h"
#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>

void discoverChannels(const char* rawFile="Data/RUN1075_100_GeV.root",
                      long maxEvents=3000)
{
    TFile* f=TFile::Open(rawFile);
    if(!f||f->IsZombie()){ printf("[discover] cannot open %s\n",rawFile); return; }
    TTree* t=(TTree*)f->Get("pulse");
    if(!t){ printf("[discover] no 'pulse' tree\n"); return; }

    const int NSLOT=36, NS=1024;
    // chanOff(drs,grp,ch)=1024*(drs*18+grp*9+ch); timeOff(drs,grp)=1024*(drs*2+grp)
    auto ampOff =[&](int s){ return s*NS; };
    auto timeOffOf=[&](int drs,int grp){ return (drs*2+grp)*NS; };

    TTreeReader reader(t);
    TTreeReaderArray<float> amp (reader,"amplitude");
    TTreeReaderArray<float> tim (reader,"timevalue");

    std::vector<long>   nValid(NSLOT,0), nTot(NSLOT,0);
    std::vector<double> sumPeak(NSLOT,0), sumSpeed(NSLOT,0), sumCharge(NSLOT,0);
    long nev=0;
    while(reader.Next() && nev<maxEvents){
        ++nev;
        for(int s=0;s<NSLOT;++s){
            int drs=s/18, grp=(s/9)%2;
            const float* A=&amp[ampOff(s)];
            const float* T=&tim[timeOffOf(drs,grp)];
            Pulse p=ExtractPulse(T,A,0.50f,8.0f);
            ++nTot[s];
            if(p.valid){
                ++nValid[s]; sumPeak[s]+=p.peak; sumCharge[s]+=p.charge;
                if(p.crossingTime>kNoTime/2 && p.peakTime>p.crossingTime)
                    sumSpeed[s]+=(p.peakTime-p.crossingTime);  // 50%->peak rise [ns]
            }
        }
    }
    f->Close();

    // ---- classify ----
    printf("\n  Discovered %ld events from %s\n",nev,rawFile);
    printf("  %-9s %5s %8s %8s  %-14s\n","slot(D,G,c)","occ%","peak[mV]","rise[ns]","inferred role");
    printf("  %s\n",std::string(60,'-').c_str());
    const char* roleName[NSLOT];
    for(int s=0;s<NSLOT;++s){
        int drs=s/18, grp=(s/9)%2, ch=s%9;
        double occ = nTot[s]? 100.*nValid[s]/nTot[s] : 0;
        double pk  = nValid[s]? sumPeak[s]/nValid[s] : 0;
        double rt  = nValid[s]? sumSpeed[s]/nValid[s] : 0;
        // ch8 of every DRS group is the DT5742 trigger/sync input (hardware, all
        // configs): occ~100%, very fast, identical across groups.
        const char* role;
        if(ch==8)                              role="TR/sync (DAQ)";
        else if(occ<5)                         role="dead/unused";
        else if(drs==1 && grp==1)              role="WC / scint";
        else if(drs==1 && grp==0)              role="LG energy (slow)";
        else if(drs==0 && rt<3.4 && pk<260)    role="MCP / fast ref";
        else if(drs==0 && pk>=350)             role="HG timing";
        else if(drs==0 && grp==1)              role="PbGlass / aux";
        else                                   role="(weak/ambiguous)";
        roleName[s]=role;
        printf("  D%dG%dc%-3d %5.0f %8.0f %8.2f  %-14s\n",drs,grp,ch,occ,pk,rt,role);
    }

    // ---- validate against the hand-coded DSB1 map ----
    printf("\n  DSB1 cross-check (expected live slots from ChannelConfig.h):\n");
    printf("   HG timing  : D0G0 c0-6 + D0G1 c0   (8 capillaries)\n");
    printf("   MCP refs   : D0G0 c7 (MCP1), D0G1 c7 (MCP2)\n");
    printf("   LG energy  : D1G0 c0-7            (8 capillaries)\n");
    printf("   WC planes  : D1G1 c1,c2,c3,c5\n");
    printf("  -> compare the 'inferred role' column above; for a NEW config the\n");
    printf("     live-slot pattern reveals its capillary count, timing vs energy split.\n");
}

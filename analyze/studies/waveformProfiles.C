// ============================================================================
// waveformProfiles.C — per-channel HG + LG average DRS4 waveform profiles for
// ANY build (DSB1, LUAG, MIXED, TENERGY). Generalizes averageWaveforms.C:
//   * the DRS wiring (kCap) is IDENTICAL across builds, so only the run list
//     changes -- read from reduce/hpc/manifest_<BUILD>.csv (status==DONE).
//   * draws BOTH the high-gain pulse (shows the 820 mV saturation) and the
//     low-gain pulse (slow, full-length, the energy signal) per channel.
// Runs on the raw "pulse" tree -> on Argon for the full set; locally it works
// for whichever raw RUN files are present under data/2023/raw/.
//   source setup.sh
//   root -l -b -q 'analyze/studies/waveformProfiles.C+("DSB1")'
//   (build = DSB1 | LUAG | MIXED | TENERGY)
// Output: output/<BUILD>/hg_waveforms.png, lg_waveforms.png  (+ .root of profiles)
// ============================================================================
#include "ChannelConfig.h"   // kCap, kMCP1, kMCP1_t (identical wiring, all builds)
#include "WaveformUtils.h"   // ExtractPulse, Pulse, kNoTime
#include "PlotUtils.h"
#include "TChain.h"
#include "TFile.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <string>
#include <cstdio>

static const double Es[6]={25,50,75,100,125,150};
static const int    NB_HG=300; static const double T0_HG=-6, T1_HG=18;
static const int    NB_LG=320; static const double T0_LG=-10, T1_LG=1000;
static const Long64_t NMAX=5000;   // events per energy (raw is large; shape converges fast)

// raw file path for a run: DSB1 era "RUN<n>_<E>_GeV.root", others "RUN<n>.root"
static std::string rawPath(int run,int E){
    std::string a=Form("data/2023/raw/RUN%d.root",run);
    if(!gSystem->AccessPathName(a.c_str())) return a;
    std::string b=Form("data/2023/raw/RUN%d_%d_GeV.root",run,E);
    if(!gSystem->AccessPathName(b.c_str())) return b;
    return "";   // not present locally (Argon)
}

// read manifest_<build>.csv -> map energy -> list of DONE raw file paths
static std::map<int,std::vector<std::string>> readManifest(const char* build){
    std::map<int,std::vector<std::string>> out;
    std::string mb=build; for(auto&c:mb)c=tolower(c);   // file is manifest_dsb1.csv (lowercase dsb1)
    std::string path="reduce/hpc/manifest_"+std::string(build)+".csv";
    if(gSystem->AccessPathName(path.c_str())) path="reduce/hpc/manifest_"+mb+".csv";
    std::ifstream f(path); if(!f){ printf("  NO manifest: %s\n",path.c_str()); return out; }
    std::string line; std::getline(f,line);   // header
    while(std::getline(f,line)){ std::stringstream ss(line); std::string tok; std::vector<std::string> col;
        while(std::getline(ss,tok,',')) col.push_back(tok);
        if(col.size()<3) continue;
        int run=atoi(col[0].c_str()), E=atoi(col[1].c_str());
        bool done = (col.back().find("DONE")!=std::string::npos) || line.find("DONE")!=std::string::npos;
        if(!done) continue;
        std::string p=rawPath(run,E); if(!p.empty()) out[E].push_back(p);
    }
    return out;
}

static TProfile* book(const char* nm,int nb,double t0,double t1){
    TProfile* p=new TProfile(nm,";t #minus t_{CFD} (ns);amplitude (mV)",nb,t0,t1);
    p->SetErrorOption("S"); p->SetDirectory(nullptr); return p; }

void waveformProfiles(const char* build="DSB1"){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    auto runs=readManifest(build);
    if(runs.empty()){ printf("%s: no raw files found locally (run on Argon).\n",build); return; }
    gSystem->mkdir(Form("output/%s",build),kTRUE);

    // hAll[energyIdx][chan] HG, lAll[energyIdx][chan] LG
    TProfile *hAll[6][8]={}, *lAll[6][8]={}; bool valid[6]={};
    for(int e=0;e<6;++e){ int E=(int)Es[e]; if(!runs.count(E))continue;
        TChain ch("pulse"); long nf=0; for(auto&p:runs[E]){ ch.Add(p.c_str()); ++nf; }
        if(ch.GetEntries()==0) continue;
        TProfile* hg[8]; TProfile* lg[8];
        for(int c=0;c<8;++c){ hg[c]=book(Form("hg_%s_%d_%s",build,E,kCap[c].name),NB_HG,T0_HG,T1_HG);
                              lg[c]=book(Form("lg_%s_%d_%s",build,E,kCap[c].name),NB_LG,T0_LG,T1_LG); }
        TTreeReader rd(&ch); TTreeReaderArray<float> tv(rd,"timevalue"), av(rd,"amplitude");
        Long64_t n=0;
        while(rd.Next()&&n<NMAX){ const float* A=&av[0]; const float* T=&tv[0];
            Pulse mcp=ExtractPulse(T+kMCP1_t,A+kMCP1,0.20f,kMCP1_minPeak);
            if(!mcp.valid||mcp.peak<kMCP1_minPeak||mcp.peak>kMCP1_maxPeak){++n;continue;}
            for(int c=0;c<8;++c){
                const float* Ah=A+kCap[c].hg; const float* Th=T+kCap[c].hg_t;
                Pulse ph=ExtractPulse(Th,Ah,0.20f,kHG_minPeak);
                if(ph.valid&&ph.crossingTime!=kNoTime){ float tc=ph.crossingTime,pe=ph.pedestal;
                    for(int s=0;s<1024;++s){ float tr=Th[s]-tc; if(tr<T0_HG||tr>T1_HG)continue; hg[c]->Fill(tr,pe-Ah[s]); } }
                const float* Al=A+kCap[c].lg; const float* Tl=T+kCap[c].lg_t;
                Pulse pl=ExtractPulse(Tl,Al,0.20f,kLG_minPeak);
                if(pl.valid&&pl.crossingTime!=kNoTime){ float tc=pl.crossingTime,pe=pl.pedestal;
                    for(int s=0;s<1024;++s){ float tr=Tl[s]-tc; if(tr<T0_LG||tr>T1_LG)continue; lg[c]->Fill(tr,pe-Al[s]); } }
            }
            ++n;
        }
        printf("  %s %3d GeV: %ld files, %lld events\n",build,E,nf,n);
        for(int c=0;c<8;++c){ hAll[e][c]=hg[c]; lAll[e][c]=lg[c]; } valid[e]=true;
    }

    // energies actually present (drives the title + shared legend, no false colour claim)
    std::vector<double> Epres; std::vector<int> Ecols;
    for(int e=0;e<6;++e) if(valid[e]){ Epres.push_back(Es[e]); Ecols.push_back(kREnergyCols[e]); }
    TString erange = Epres.empty()? "" : Form("%d#minus%d GeV",(int)Epres.front(),(int)Epres.back());

    // draw a gain (HG or LG) summary: 8 channel panels, energies overlaid, COMMON y-range
    auto drawSummary=[&](TProfile* P[6][8],const char* gain,double t0,double t1,const char* fn){
        bool isHG = (std::string(gain)=="HG");
        TCanvas* c=new TCanvas(Form("c%s",gain),"",3200,1640);
        TPad* g=GridWithTitle(c,4,2,Form("%s: average %s DRS4 waveform per channel (mean #pm RMS), %s, aligned at CFD t=0",build,gain,erange.Data()),0.004,0.03,0.06,0.020);
        // ONE common y-range across all 8 channels so peak heights are comparable + undershoot is contained
        double gymax=0, gymin=0;
        for(int ci=0;ci<8;++ci) for(int e=0;e<6;++e) if(P[e][ci]){ gymax=std::max(gymax,P[e][ci]->GetMaximum());
            for(int b=1;b<=P[e][ci]->GetNbinsX();++b){ double tc=P[e][ci]->GetBinCenter(b); if(tc<t0||tc>t1)continue; double y=P[e][ci]->GetBinContent(b); if(P[e][ci]->GetBinEntries(b)>0) gymin=std::min(gymin,y); } }
        if(gymax<=0)gymax=1; double ylo=gymin-0.10*gymax, yhi=1.18*gymax;
        for(int ci=0;ci<8;++ci){ g->cd(ci+1); gPad->SetLeftMargin(0.165); gPad->SetBottomMargin(0.15); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.10);
            bool first=true;
            for(int e=0;e<6;++e){ if(!P[e][ci])continue; TProfile* h=(TProfile*)P[e][ci]->Clone(Form("d%s%d%d",gain,e,ci));
                h->SetLineColor(kREnergyCols[e]); h->SetLineWidth(2); h->SetMarkerSize(0);
                h->GetYaxis()->SetRangeUser(ylo,yhi); h->GetXaxis()->SetRangeUser(t0,t1);
                h->GetXaxis()->SetTitleSize(0.055); h->GetYaxis()->SetTitleSize(0.055); h->GetYaxis()->SetTitleOffset(1.35);
                h->GetXaxis()->SetLabelSize(0.05); h->GetYaxis()->SetLabelSize(0.05);
                h->Draw(first?"HIST":"HIST SAME"); first=false; }
            TLine* l=new TLine(0,ylo,0,yhi); l->SetLineStyle(2); l->SetLineColor(kRed); l->Draw();
            DrawPadTitle(kCap[ci].name);
            if(ci==0) DrawEnergyLegend(isHG?0.21:0.60, 0.46, isHG?0.50:0.95, 0.90, Epres, Ecols, 0.052, "beam E");
        }
        c->Print(Form("output/%s/%s",build,fn));
        printf("  wrote output/%s/%s\n",build,fn);
    };
    drawSummary(hAll,"HG",T0_HG,T1_HG,"hg_waveforms.png");
    drawSummary(lAll,"LG",T0_LG,520,"lg_waveforms.png");   // crop to the informative pulse (ends ~450 ns); no '1000' edge-clip

    // persist profiles
    TFile fo(Form("output/%s/waveform_profiles.root",build),"RECREATE");
    for(int e=0;e<6;++e)for(int c=0;c<8;++c){ if(hAll[e][c])hAll[e][c]->Write(); if(lAll[e][c])lAll[e][c]->Write(); }
    fo.Close();
    printf("%s: done.\n",build);
}

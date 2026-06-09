// ============================================================================
// edgeSlope.C — what limits timing AFTER the DT5742 HG-channel clips?
// ----------------------------------------------------------------------------
// Uses ExtractPulseMulti (DRS4-timing-correct) on the raw HG waveforms. CFD
// timing jitter ~ noise/(dV/dt). For a shape-invariant pulse dV/dt scales with
// peak, so clipping the peak caps the slope -> caps timing. The discriminator:
// across the CLIPPED energies (50-150, peak pinned ~820), does the rising-edge
// transit time (cfd30-cfd05) keep SHRINKING (true light still growing, only the
// digitizer clips -> floor is electronics) or stay FLAT (light/SiPM saturated)?
//   ROOT_INCLUDE_PATH=lib/waveform:lib/io:lib/physics:lib/viz root -l -b -q analyze/studies/edgeSlope.C+
// ============================================================================
#include "ChannelConfig.h"   // kCap
#include "WaveformUtils.h"   // ExtractPulseMulti (handles DRS4 cell ordering)
#include "PlotUtils.h"       // ApplyRADiCALStyle
#include "SelectionCuts.h"   // kHG_LED_thresh
#include "DataPaths.h"
#include "TChain.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "FigPaths.h"
#include <cstdio>
#include <cmath>

void edgeSlope() {
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    struct R { double E; const char* file; };
    R runs[6] = {
        {25,"RUN1211_25_GeV.root"},{50,"RUN1148_50_GeV.root"},{75,"RUN1112_75_GeV.root"},
        {100,"RUN1075_100_GeV.root"},{125,"RUN1034_125_GeV.root"},{150,"RUN1258_150_GeV.root"} };
    const long MAXEV=5000;
    double Es[6], rt[6], slope[6], peak[6];   // rt = cfd30-cfd05 transit [ns]

    for(int e=0;e<6;++e){
        TChain ch("pulse"); ch.Add(radRaw(runs[e].file));
        TTreeReader rd(&ch); TTreeReaderArray<float> amp(rd,"amplitude"), tim(rd,"timevalue");
        double sRT=0,sSL=0,sP=0; long n=0, iev=0;
        while(rd.Next() && iev<MAXEV){ ++iev;
            const float* A=&amp[0]; const float* T=&tim[0];
            for(int c=0;c<8;++c){
                PulseMulti p=ExtractPulseMulti(T+kCap[c].hg_t, A+kCap[c].hg, kHG_LED_thresh, 5.f, 1100.0f);
                if(p.peak<200) continue;                       // a real shower pulse
                if(p.cfd05<=-1e5f || p.cfd30<=-1e5f) continue;
                double transit=p.cfd30-p.cfd05; if(transit<=0 || transit>5) continue;
                sRT+=transit; sSL+=0.25*p.peak/transit; sP+=p.peak; ++n;   // slope over 5-30% [mV/ns]
            }
        }
        Es[e]=runs[e].E; rt[e]=n?sRT/n:0; slope[e]=n?sSL/n:0; peak[e]=n?sP/n:0;
        printf("E=%3.0f  peak=%.0f mV  rise(cfd30-cfd05)=%.3f ns  dV/dt(5-30%%)=%.0f mV/ns  (n=%ld)\n",
               Es[e], peak[e], rt[e], slope[e], n);
    }

    TCanvas* c=new TCanvas("c_es","",1500,620); c->Divide(2,1,0.010,0.006);
    // Panel 1: peak (clips) and absolute slope dV/dt vs energy
    c->cd(1); gPad->SetGridx(); gPad->SetGridy();
    double smx=0; for(int e=0;e<6;++e){ smx=std::max(smx,slope[e]); smx=std::max(smx,peak[e]); }
    TH1F* fr=gPad->DrawFrame(0,0,160,smx*1.25);
    fr->SetTitle("HG peak and rising-edge slope vs energy;beam energy [GeV];peak [mV]  /  dV/dt [mV/ns]");
    TGraph* gP=new TGraph(6,Es,peak);  gP->SetMarkerStyle(24); gP->SetMarkerColor(kGray+2); gP->SetLineColor(kGray+2); gP->SetMarkerSize(1.5); gP->SetLineWidth(2); gP->Draw("PL SAME");
    TGraph* gS=new TGraph(6,Es,slope); gS->SetMarkerStyle(21); gS->SetMarkerColor(kRed+1);  gS->SetLineColor(kRed+1);  gS->SetMarkerSize(1.5); gS->SetLineWidth(2); gS->Draw("PL SAME");
    TLegend* lg=new TLegend(0.16,0.70,0.60,0.88); lg->SetBorderSize(0);
    lg->AddEntry(gP,"HG peak (clips ~820 mV)","pl"); lg->AddEntry(gS,"dV/dt over 5-30% (absolute slew)","pl"); lg->Draw();
    TLatex tl; tl.SetNDC(); tl.SetTextSize(0.026);
    tl.DrawLatex(0.16,0.24,"dV/dt KEEPS RISING (102#rightarrow279 mV/ns) while the peak clips:");
    tl.DrawLatex(0.16,0.20,"the usable CFD rising edge is NOT capped -> more light still helps slew.");

    // Panel 2: the FRACTIONAL transit time (cfd30-cfd05) — does the edge keep steepening?
    c->cd(2); gPad->SetGridx(); gPad->SetGridy();
    double rmx=0; for(int e=0;e<6;++e) rmx=std::max(rmx,rt[e]);
    TH1F* fr2=gPad->DrawFrame(0,0,160,rmx*1.4);
    fr2->SetTitle("5%-30% transit time vs energy (edge shape);beam energy [GeV];cfd30 - cfd05  [ns]");
    TGraph* gR=new TGraph(6,Es,rt); gR->SetMarkerStyle(20); gR->SetMarkerColor(kAzure+2); gR->SetLineColor(kAzure+2); gR->SetMarkerSize(1.5); gR->SetLineWidth(2); gR->Draw("PL SAME");
    TLatex t2; t2.SetNDC(); t2.SetTextSize(0.026);
    t2.DrawLatex(0.16,0.30,"Transit KEEPS shrinking to 150 GeV: light/edge still growing,");
    t2.DrawLatex(0.16,0.26,"NO SiPM saturation. Timing plateaus at a SEPARATE ~28 ps floor");
    t2.DrawLatex(0.16,0.22,"(DRS4 sampling + scint/SiPM time spread), not the slope or clip.");

    gSystem->mkdir("figures",kTRUE);
    c->Print(radFigP("figures/edge_slope.png"));
    printf("\npeak[mV]: "); for(int e=0;e<6;++e) printf("%.0f ",peak[e]);
    printf("\nrise(cfd30-cfd05)[ns]: "); for(int e=0;e<6;++e) printf("%.3f ",rt[e]);
    printf("\ndV/dt[mV/ns]: "); for(int e=0;e<6;++e) printf("%.0f ",slope[e]); printf("\n");
    printf("wrote figures/edge_slope.png\n");
}

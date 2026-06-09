// ============================================================================
// slopeVsThresh.C — WHERE on the rising edge do we time, and does the slope
//                   there actually rise with energy?  (tests the user's point)
// ----------------------------------------------------------------------------
// We time at a LOW level: hg_led = fixed 20 mV, hg_cfd05 = 5% of the CLIPPED
// ~820 mV peak ~ 41 mV. Hypothesis: at the FOOT of the pulse dV/dt is low AND
// nearly energy-independent (the early rise looks the same regardless of total
// amplitude), so sigma_slew = noise/(dV/dt) is large and FLAT -> timing does not
// improve even though the STEEP part of the edge gets faster with energy.
// Measure, from raw, the local dV/dt at the crossing of several FIXED absolute
// thresholds vs energy, and the implied single-channel sigma_slew = ped_rms/(dV/dt).
//   ROOT_INCLUDE_PATH=lib/waveform:lib/io:lib/physics:lib/viz root -l -b -q analyze/studies/slopeVsThresh.C+
// ============================================================================
#include "ChannelConfig.h"
#include "WaveformUtils.h"
#include "PlotUtils.h"
#include "SelectionCuts.h"
#include "DataPaths.h"
#include "TChain.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "FigPaths.h"
#include <vector>
#include <cmath>
#include <cstdio>

static const double LV[6]={20,50,100,200,300,400};   // fixed absolute thresholds [mV]

void slopeVsThresh(){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    struct R{double E; const char* f;};
    R runs[6]={{25,"RUN1211_25_GeV.root"},{50,"RUN1148_50_GeV.root"},{75,"RUN1112_75_GeV.root"},
               {100,"RUN1075_100_GeV.root"},{125,"RUN1034_125_GeV.root"},{150,"RUN1258_150_GeV.root"}};
    double Es[6]; double dvdt[6][6], pedrms[6];   // [energy][level]
    const int G0[7]={0,1,2,3,4,5,6};

    for(int e=0;e<6;++e){ Es[e]=runs[e].E;
        TChain ch("pulse"); ch.Add(radRaw(runs[e].f)); TTreeReader rd(&ch);
        TTreeReaderArray<float> amp(rd,"amplitude"), tim(rd,"timevalue");
        double sdvdt[6]={0}; long ndvdt[6]={0}; double sped=0; long nped=0; long iev=0;
        while(rd.Next() && iev<20000){ ++iev; const float* A=&amp[0]; const float* T=&tim[0];
            Pulse m1=ExtractPulse(T+kMCP1_t,A+kMCP1,0.2f,30.f); if(m1.peak<kMCP1_minPeak||m1.peak>kMCP1_maxPeak) continue;
            for(int j=0;j<7;++j){ int c=G0[j]; const float* a=A+kCap[c].hg; const float* tt=T+kCap[c].hg_t;
                // pedestal + peak (negative pulse: sig = ped - a)
                double ped=0; for(int s=3;s<53;++s) ped+=a[s]; ped/=50.0;
                double rms2=0; for(int s=3;s<53;++s){double d=a[s]-ped; rms2+=d*d;} double prms=std::sqrt(rms2/50.0);
                int imin=3; float vmin=a[3]; for(int s=3;s<1000;++s) if(a[s]<vmin){vmin=a[s];imin=s;}
                double peak=ped-vmin; if(peak<60) continue;          // need a real pulse
                sped+=prms; ++nped;
                for(int L=0;L<6;++L){ double tgt=LV[L]; if(tgt>=0.9*peak) continue;   // unclipped level only
                    for(int s=3;s<imin;++s){ double si=ped-a[s], si1=ped-a[s+1];
                        if(si<tgt && si1>=tgt){ double dt=tt[s+1]-tt[s]; if(dt>0){ double slope=(si1-si)/dt; // mV/ns
                            if(slope>0){ sdvdt[L]+=slope; ++ndvdt[L]; } } break; } }
                }
            }
        }
        pedrms[e]= nped? sped/nped : 0;
        printf("E=%3.0f  ped_rms=%.2f mV |", Es[e], pedrms[e]);
        for(int L=0;L<6;++L){ dvdt[e][L]= ndvdt[L]? sdvdt[L]/ndvdt[L] : 0;
            printf("  %3.0fmV:%.0f", LV[L], dvdt[e][L]); }
        printf("  (mV/ns)\n");
    }
    // sigma_slew = ped_rms / dvdt  at each level, vs energy
    printf("\nimplied single-channel sigma_slew = ped_rms/(dV/dt)  [ps]:\n");
    for(int e=0;e<6;++e){ printf("E=%3.0f |", Es[e]);
        for(int L=0;L<6;++L){ double s = dvdt[e][L]>0 ? 1000.0*pedrms[e]/dvdt[e][L] : 0; printf("  %3.0fmV:%4.1f", LV[L], s); }
        printf("\n"); }

    // plot dV/dt vs energy for low (20 mV) vs steep (300 mV) thresholds
    TCanvas* c=new TCanvas("c_sv","",900,680); c->SetGridx(); c->SetGridy();
    TH1F* fr=c->DrawFrame(0,0,160,700);
    fr->SetTitle("Rising-edge dV/dt at the crossing: foot vs steep part;beam energy (GeV);dV/dt at threshold (mV/ns)");
    int Ls[3]={0,3,5}; int col[3]={(int)(kGray+2),(int)(kAzure+1),(int)(kRed+1)};
    for(int k=0;k<3;++k){ int L=Ls[k]; double y[6]; for(int e=0;e<6;++e) y[e]=dvdt[e][L];
        TGraph* g=new TGraph(6,Es,y); g->SetMarkerStyle(20+k); g->SetMarkerColor(col[k]); g->SetLineColor(col[k]);
        g->SetMarkerSize(1.6); g->SetLineWidth(2); g->Draw("PL SAME"); }
    TLegend* lg=new TLegend(0.16,0.70,0.60,0.88); lg->SetBorderSize(0);
    lg->AddEntry((TObject*)0,"dV/dt at crossing of:","");
    for(int k=0;k<3;++k){ TGraph* gg=new TGraph(); gg->SetMarkerStyle(20+k); gg->SetMarkerColor(col[k]); gg->SetLineColor(col[k]);
        lg->AddEntry(gg,Form("%.0f mV %s",LV[Ls[k]], k==0?"(where we time now)":""),"pl"); }
    lg->Draw();
    gSystem->mkdir("figures",kTRUE); c->Print(radFigP("figures/slope_vs_thresh.png"));
    printf("wrote figures/slope_vs_thresh.png\n");
}

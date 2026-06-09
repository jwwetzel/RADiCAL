// ============================================================================
// lgCFD.C — CFD on the LG-predicted TRUE peak (the user's method), from raw.
// ----------------------------------------------------------------------------
// The HG clips ~820 mV so its peak is unknown -> CFD-of-clipped-peak sits at the
// foot (cfd05 ~ 41 mV) where dV/dt is low. Instead:
//   1. Fit HG_peak = a + b*LG_peak per channel on UNCLIPPED events (HG in [30,700])
//      -- the linear part of the "hockey stick".
//   2. Per event: HG_true = a + b*LG_peak  (LG never clips -> stable).
//   3. Threshold = frac * HG_true, on the STEEP part of the edge (below the clip).
//   4. Interpolate the exact crossing time at that ABSOLUTE threshold.
// Compare best-bin (bright-slice) sigma_t vs energy: baseline cfd05 vs lgCFD@frac,
// for DEPTH (DW-UP)/2 and ABSOLUTE mean(hg)-MCP.
//   ROOT_INCLUDE_PATH=lib/waveform:lib/io:lib/physics:lib/viz root -l -b -q analyze/studies/lgCFD.C+
// ============================================================================
#include "ChannelConfig.h"
#include "WaveformUtils.h"
#include "PlotUtils.h"
#include "RadTiming.h"
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
#include <algorithm>
#include <cmath>
#include <cstdio>

// crossing time at an ABSOLUTE threshold thr (mV) on the rising edge (negative pulse)
static double crossAbs(const float* t, const float* a, double ped, int imin, double thr){
    if(thr<=0) return -1e9;
    for(int i=3;i<imin;++i){ double si=ped-a[i], si1=ped-a[i+1];
        if(si<thr && si1>=thr){ double dt=t[i+1]-t[i]; if(dt<=0) return -1e9;
            double slope=(si1-si)/dt; if(slope<=0) return -1e9;
            return t[i] + (thr-si)/slope; } }
    return -1e9;
}

struct Ev { float slg,x,y; float dN,uN,absN; int dn,un,an; float dB,uB; int dnb,unb; float mcp; };

static double medOf(std::vector<float> v){ if(v.empty())return 0; size_t k=v.size()/2; std::nth_element(v.begin(),v.begin()+k,v.end()); return v[k]; }
static double bc(std::vector<std::pair<float,float>>& v){
    if(v.size()<500) return -1; std::sort(v.begin(),v.end());
    std::vector<float> t; for(size_t i=0.55*v.size();i<v.size();++i) t.push_back(v[i].second);
    double m=medOf(t); std::vector<float> c; for(float x:t) if(std::fabs(x-m)<1.0) c.push_back(x);
    if(c.size()<250) return -1; return rad::tebSigma(c);
}

void lgCFD(){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    struct R{double E; const char* f;};
    R runs[6]={{25,"RUN1211_25_GeV.root"},{50,"RUN1148_50_GeV.root"},{75,"RUN1112_75_GeV.root"},
               {100,"RUN1075_100_GeV.root"},{125,"RUN1034_125_GeV.root"},{150,"RUN1258_150_GeV.root"}};
    const int G0[7]={0,1,2,3,4,5,6};
    const double FR[3]={0.15,0.20,0.30};
    double Es[6], sBase_d[6], sBase_a[6], sLG_d[3][6], sLG_a[3][6];

    // (1) global per-channel HG=a+b*LG fit from the 150 GeV run (best stats, full dynamic range)
    double aa[8]={0}, bb[8]={5};
    { TChain ch("pulse"); ch.Add(radRaw(runs[5].f)); TTreeReader rd(&ch);
      TTreeReaderArray<float> amp(rd,"amplitude"), tim(rd,"timevalue");
      double sx[8]={0},sy[8]={0},sxx[8]={0},sxy[8]={0}; long n[8]={0};
      while(rd.Next()){ const float*A=&amp[0]; const float*T=&tim[0];
        for(int c=0;c<8;++c){ Pulse hg=ExtractPulse(T+kCap[c].hg_t, A+kCap[c].hg,0.2f,5.f);
            Pulse lg=ExtractPulse(T+kCap[c].lg_t, A+kCap[c].lg,0.2f,5.f);
            if(hg.peak>30&&hg.peak<700&&lg.peak>10){ double X=lg.peak,Y=hg.peak; sx[c]+=X;sy[c]+=Y;sxx[c]+=X*X;sxy[c]+=X*Y;++n[c]; } } }
      for(int c=0;c<8;++c){ if(n[c]>100){ double d=n[c]*sxx[c]-sx[c]*sx[c]; bb[c]=(n[c]*sxy[c]-sx[c]*sy[c])/d; aa[c]=(sy[c]-bb[c]*sx[c])/n[c]; } }
      printf("HG=a+b*LG fit (150 GeV): "); for(int c=0;c<8;++c) printf("c%d:%.0f+%.2fLG ",c,aa[c],bb[c]); printf("\n");
    }

    for(int e=0;e<6;++e){ Es[e]=runs[e].E;
        TChain ch("pulse"); ch.Add(radRaw(runs[e].f)); TTreeReader rd(&ch);
        TTreeReaderArray<float> amp(rd,"amplitude"), tim(rd,"timevalue");
        std::vector<std::pair<float,float>> base_d, base_a, lg_d[3], lg_a[3];
        std::vector<float> X,Y,SLG;  // for fiducial
        struct Row{float x,y,slg; float bd,bu; int bdn,bun; float Ld[3],Lu[3],La[3]; int Ldn[3],Lun[3],Lan[3]; float mcp;};
        std::vector<Row> rows;
        while(rd.Next()){ const float*A=&amp[0]; const float*T=&tim[0];
            Pulse wcR=ExtractPulse(T+kWC_t,A+kWC_R,0.5f,kWC_minPeak), wcL=ExtractPulse(T+kWC_t,A+kWC_L,0.5f,kWC_minPeak),
                  wcD=ExtractPulse(T+kWC_t,A+kWC_D,0.5f,kWC_minPeak), wcU=ExtractPulse(T+kWC_t,A+kWC_U,0.5f,kWC_minPeak);
            if(!(wcR.valid&&wcL.valid&&wcD.valid&&wcU.valid)) continue;
            Pulse m1=ExtractPulse(T+kMCP1_t,A+kMCP1,0.2f,30.f); if(m1.peak<kMCP1_minPeak||m1.peak>kMCP1_maxPeak) continue;
            Row r; r.x=kWC_Scale*(wcR.peakTime-wcL.peakTime); r.y=kWC_Scale*(wcD.peakTime-wcU.peakTime); r.mcp=m1.crossingTime;
            r.bd=r.bu=0; r.bdn=r.bun=0; for(int k=0;k<3;++k){r.Ld[k]=r.Lu[k]=r.La[k]=0; r.Ldn[k]=r.Lun[k]=r.Lan[k]=0;}
            double slg=0;
            for(int j=0;j<7;++j){ int c=G0[j]; const float* a=A+kCap[c].hg; const float* tt=T+kCap[c].hg_t;
                double ped=0; for(int s=3;s<53;++s) ped+=a[s]; ped/=50.0;
                int imin=3; float vmin=a[3]; for(int s=3;s<1000;++s) if(a[s]<vmin){vmin=a[s];imin=s;}
                double hgpk=ped-vmin; if(hgpk<60) continue;
                Pulse lg=ExtractPulse(T+kCap[c].lg_t,A+kCap[c].lg,0.2f,5.f); slg+=lg.peak;
                double HGtrue=aa[c]+bb[c]*lg.peak;
                // baseline cfd05 (5% of MEASURED/clipped peak)
                double t05=crossAbs(tt,a,ped,imin,0.05*hgpk);
                if(t05>-1e5){ double tt05=t05-m1.crossingTime; if(c<4){r.bd+=tt05;++r.bdn;} else {r.bu+=tt05;++r.bun;} }
                // lgCFD at frac of TRUE peak
                for(int k=0;k<3;++k){ double thr=FR[k]*HGtrue;
                    if(thr<20||thr>780) continue;                 // clean steep edge, below clip
                    double tc=crossAbs(tt,a,ped,imin,thr);
                    if(tc>-1e5){ double tcm=tc-m1.crossingTime; if(c<4){r.Ld[k]+=tcm;++r.Ldn[k];} else {r.Lu[k]+=tcm;++r.Lun[k];}
                                 r.La[k]+=tcm; ++r.Lan[k]; } }
            }
            r.slg=slg; rows.push_back(r);
        }
        // fiducial
        double wx=0,wy=0,w=0; for(auto&r:rows){wx+=r.x*r.slg;wy+=r.y*r.slg;w+=r.slg;}
        double xc=w>0?wx/w:0,yc=w>0?wy/w:0,rF=TimingFiducialR(Es[e]),r2=rF*rF;
        for(auto&r:rows){ double dx=r.x-xc,dy=r.y-yc; if(dx*dx+dy*dy>=r2) continue;
            if(r.bdn>=1&&r.bun>=1){ base_d.push_back({r.slg,0.5f*(r.bd/r.bdn-r.bu/r.bun)});
                // baseline absolute: mean of all baseline hg (down+up) /count
                base_a.push_back({r.slg,(r.bd+r.bu)/(r.bdn+r.bun)}); }
            for(int k=0;k<3;++k){ if(r.Ldn[k]>=1&&r.Lun[k]>=1) lg_d[k].push_back({r.slg,0.5f*(r.Ld[k]/r.Ldn[k]-r.Lu[k]/r.Lun[k])});
                if(r.Lan[k]>=4) lg_a[k].push_back({r.slg,r.La[k]/r.Lan[k]}); }
        }
        sBase_d[e]=bc(base_d); sBase_a[e]=bc(base_a);
        printf("E=%3.0f  DEPTH: cfd05=%.1f", Es[e], sBase_d[e]);
        for(int k=0;k<3;++k){ sLG_d[k][e]=bc(lg_d[k]); printf("  lg%.0f%%=%.1f", FR[k]*100, sLG_d[k][e]); }
        printf(" ps | ABS: cfd05=%.1f", sBase_a[e]);
        for(int k=0;k<3;++k){ sLG_a[k][e]=bc(lg_a[k]); printf("  lg%.0f%%=%.1f", FR[k]*100, sLG_a[k][e]); }
        printf(" ps\n");
    }
    // plot DEPTH: cfd05 vs best lgCFD (20%)
    TCanvas* c=new TCanvas("c_lg","",900,680); c->SetGridx(); c->SetGridy();
    TH1F* fr=c->DrawFrame(0,0,160,80);
    fr->SetTitle("CFD on LG-predicted TRUE peak vs cfd05-of-clipped;beam energy (GeV);best-bin #sigma_{t} (ps)");
    TGraph* gB=new TGraph(6,Es,sBase_d); gB->SetMarkerStyle(24); gB->SetMarkerColor(kGray+2); gB->SetLineColor(kGray+2); gB->SetMarkerSize(1.5); gB->SetLineWidth(2); gB->Draw("PL SAME");
    TGraph* gL=new TGraph(6,Es,sLG_d[1]); gL->SetMarkerStyle(20); gL->SetMarkerColor(kRed+1); gL->SetLineColor(kRed+1); gL->SetMarkerSize(1.6); gL->SetLineWidth(2); gL->Draw("PL SAME");
    TGraph* gBa=new TGraph(6,Es,sBase_a); gBa->SetMarkerStyle(26); gBa->SetMarkerColor(kAzure+1); gBa->SetLineColor(kAzure+1); gBa->SetMarkerSize(1.5); gBa->SetLineWidth(2); gBa->SetLineStyle(2); gBa->Draw("PL SAME");
    TGraph* gLa=new TGraph(6,Es,sLG_a[1]); gLa->SetMarkerStyle(22); gLa->SetMarkerColor(kGreen+2); gLa->SetLineColor(kGreen+2); gLa->SetMarkerSize(1.6); gLa->SetLineWidth(2); gLa->SetLineStyle(2); gLa->Draw("PL SAME");
    TLegend* lg=new TLegend(0.36,0.66,0.88,0.88); lg->SetBorderSize(0);
    lg->AddEntry(gB,"DEPTH cfd05 (clipped)","pl"); lg->AddEntry(gL,"DEPTH lgCFD 20% (true peak)","pl");
    lg->AddEntry(gBa,"ABS cfd05 (clipped)","pl"); lg->AddEntry(gLa,"ABS lgCFD 20% (true peak)","pl"); lg->Draw();
    gSystem->mkdir("figures",kTRUE); c->Print(radFigP("figures/lg_cfd.png"));
    printf("wrote figures/lg_cfd.png\n");
}

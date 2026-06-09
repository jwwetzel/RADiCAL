// ============================================================================
// narrativeFigs.C — the two crux figures for the timing narrative:
//   (1) narrative_clip.png  — HG peak spectrum piling up at the 820 mV clip (the wall)
//   (2) narrative_pulse.png — ONE clipped HG pulse: cfd05 stuck at the low-slope
//       foot vs hg_lgcfd on the steep edge, with the LG-predicted true peak.
//   source setup.sh
//   root -l -b -q 'analyze/studies/narrativeFigs.C+("DSB1","data/2023/raw/RUN1034_125_GeV.root",125,1)'
// ============================================================================
#include "BuildConfig.h"
#include "WaveformUtils.h"
#include "PlotUtils.h"
#include "DataPaths.h"
#include "TChain.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLatex.h"
#include "TMarker.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>
using namespace rad;

void narrativeFigs(const char* build="DSB1", const char* rawPath="data/2023/raw/RUN1034_125_GeV.root",
                   double E=125, int ch=1){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    BuildConfig cfg = BuildConfig::Load(radConfig(build).Data());
    if(!cfg.valid()){ printf("config load failed\n"); return; }
    const EndMap& c = cfg.end[ch];
    double a = cfg.hg_lg_a[ch], b = cfg.hg_lg_b[ch];
    printf("[narrativeFigs] %s ch%d=%s  HG_true = %.1f + %.2f*LG  (frac=%.2f)\n",
           build, ch, c.name.c_str(), a, b, cfg.lgcfd_frac);

    TChain tr("pulse"); tr.Add(rawPath); TTreeReader rd(&tr);
    TTreeReaderArray<float> amp(rd,"amplitude"), tim(rd,"timevalue");

    TH1F* hClip=new TH1F("hClip","",170,0,950); hClip->SetDirectory(0);
    // pick a representative clipped event: flat-topped HG near the clip, clean LG
    std::vector<double> wt, wv; double clipPk=0, hgTrue=0, lgPk=0, ped=0, pol=1; bool got=false;
    long cnt=0;
    while(rd.Next() && cnt<60000){ ++cnt;
        const float* A=&amp[0]; const float* T=&tim[0];
        Pulse hg=ExtractPulse(T+c.hg_t,A+c.hg,0.20f,5.f), lg=ExtractPulse(T+c.lg_t,A+c.lg,0.20f,5.f);
        if(hg.peak>40&&hg.peak<950) hClip->Fill(hg.peak);
        if(got) continue;
        if(!(hg.peak>805&&hg.peak<835&&lg.peak>20)) continue;     // clipped + good LG
        // build the displayed (pedestal-subtracted, positive-going) HG waveform
        const float* ah=A+c.hg; const float* th=T+c.hg_t;
        double p=0; for(int s=2;s<52;++s) p+=ah[s]; p/=50.0;
        int xs=2; for(int s=2;s<1000;++s) if(std::fabs(ah[s]-p)>std::fabs(ah[xs]-p)) xs=s;
        double pl=(ah[xs]-p)>=0?1.0:-1.0;
        std::vector<double> sig(1024); double mx=-1e9; int flat=0;
        for(int s=0;s<1024;++s){ sig[s]=pl*(ah[s]-p); if(sig[s]>mx) mx=sig[s]; }
        for(int s=0;s<1024;++s) if(mx-sig[s]<6) ++flat;
        if(mx<790||flat<4) continue;                               // require a real clip plateau
        wt.assign(th,th+1024); wv=sig; ped=p; pol=pl; clipPk=mx; lgPk=lg.peak;
        hgTrue=a+b*lg.peak; got=true;
    }
    if(!got){ printf("  no clean clipped event found (try another run/channel)\n"); }

    // ---- Figure 1: the clip spectrum -----------------------------------------
    { TCanvas* c1=new TCanvas("c_clip","",820,560); c1->SetLeftMargin(0.13); c1->SetRightMargin(0.05);
      hClip->SetLineColor(kAzure+2); hClip->SetLineWidth(2); hClip->SetFillColorAlpha(kAzure+2,0.25);
      hClip->GetXaxis()->SetTitle("HG pulse peak amplitude (mV)");
      hClip->GetYaxis()->SetTitle("events"); hClip->Draw("HIST");
      double ymax=hClip->GetMaximum();
      TLine* lc=new TLine(820,0,820,ymax*1.02); lc->SetLineColor(kRed+1); lc->SetLineStyle(2); lc->SetLineWidth(2); lc->Draw();
      TLatex t; t.SetNDC(); t.SetTextFont(62); t.SetTextSize(0.045);
      t.DrawLatex(0.16,0.93,Form("%s  %s  @ %.0f GeV : HG saturates at the 820 mV clip", build, c.name.c_str(), E));
      t.SetTextFont(42); t.SetTextSize(0.038); t.SetTextColor(kRed+1);
      t.DrawLatex(0.60,0.78,"clip @ 820 mV"); t.SetTextColor(kGray+3);
      t.DrawLatex(0.60,0.72,"#rightarrow true peak unknown");
      gSystem->mkdir("figures/narrative",kTRUE); c1->Print("figures/narrative/narrative_clip.png");
      printf("  wrote figures/narrative/narrative_clip.png\n"); }

    if(!got) return;
    // ---- Figure 2: one clipped pulse, cfd05-foot vs lgcfd-edge ---------------
    // crossing time on the rising edge at absolute threshold thr
    auto cross=[&](double thr)->double{ for(size_t s=1;s<wv.size();++s)
        if(wv[s-1]<thr && wv[s]>=thr && s>3){ double f=(thr-wv[s-1])/(wv[s]-wv[s-1]); return wt[s-1]+f*(wt[s]-wt[s-1]); } return -1; };
    double thr05=0.05*clipPk, thrLG=cfg.lgcfd_frac*hgTrue;
    double t05=cross(thr05), tLG=cross(thrLG);
    // local slope dV/dt around each crossing (mV/ns)
    auto slope=[&](double tc)->double{ int s=1; while(s<(int)wt.size()&&wt[s]<tc) ++s; if(s<2||s>=(int)wt.size()) return 0;
        return (wv[s]-wv[s-2])/(wt[s]-wt[s-2]); };

    // zoom window around the rising edge
    double tlo=(t05>0?t05:wt[0])-3, thi=(tLG>0?tLG:wt[0])+8;
    TGraph* g=new TGraph(); for(size_t s=0;s<wt.size();++s) if(wt[s]>tlo-2&&wt[s]<thi+6) g->SetPoint(g->GetN(),wt[s],wv[s]);
    TCanvas* c2=new TCanvas("c_pulse","",900,620); c2->SetLeftMargin(0.12); c2->SetRightMargin(0.05);
    g->SetTitle(""); g->SetLineColor(kBlack); g->SetLineWidth(3);
    g->GetXaxis()->SetTitle("time (ns)"); g->GetYaxis()->SetTitle("HG signal, pedestal-subtracted (mV)");
    g->GetXaxis()->SetLimits(tlo,thi); g->GetYaxis()->SetRangeUser(-40, std::max(hgTrue,clipPk)*1.12);
    g->Draw("AL");
    auto hline=[&](double y,int col,int sty){ TLine* L=new TLine(tlo,y,thi,y); L->SetLineColor(col); L->SetLineStyle(sty); L->SetLineWidth(2); L->Draw(); return L; };
    hline(clipPk,kGray+1,2); hline(hgTrue,kGreen+2,2); hline(thr05,kRed+1,1); hline(thrLG,kAzure+2,1);
    auto vmark=[&](double tc,double thr,int col){ if(tc<0)return; TMarker* m=new TMarker(tc,thr,20); m->SetMarkerColor(col); m->SetMarkerSize(1.6); m->Draw();
        TLine* L=new TLine(tc,-40,tc,thr); L->SetLineColor(col); L->SetLineStyle(3); L->Draw(); };
    vmark(t05,thr05,kRed+1); vmark(tLG,thrLG,kAzure+2);
    TLatex t; t.SetTextFont(42); t.SetTextSize(0.034);
    double xa=tlo+0.46*(thi-tlo), xb=tlo+0.30*(thi-tlo);
    t.SetTextColor(kGray+2);  t.DrawLatex(xa,clipPk+60,"clipped peak (820 mV)");
    t.SetTextColor(kGreen+3); t.DrawLatex(xb,hgTrue-120,Form("LG-predicted true peak = %.0f mV",hgTrue));
    t.SetTextColor(kRed+1);   t.DrawLatex(tlo+0.2,thr05+14,Form("cfd05 @ %.0f mV  (foot: %.0f mV/ns)",thr05,std::fabs(slope(t05))));
    t.SetTextColor(kAzure+2); t.DrawLatex(tlo+0.2,thrLG+14,Form("hg_lgcfd @ %.0f mV  (edge: %.0f mV/ns)",thrLG,std::fabs(slope(tLG))));
    TLatex tt; tt.SetNDC(); tt.SetTextFont(62); tt.SetTextSize(0.040);
    tt.DrawLatex(0.13,0.94,Form("%s %s @ %.0f GeV: time the steep edge, not the clipped foot", build, c.name.c_str(), E));
    c2->Print("figures/narrative/narrative_pulse.png");
    printf("  wrote figures/narrative/narrative_pulse.png  (cfd05 t=%.3f slope=%.0f | lgcfd t=%.3f slope=%.0f mV/ns)\n",
           t05,std::fabs(slope(t05)),tLG,std::fabs(slope(tLG)));
}

// ============================================================================
// trSync.C — TR0 (the 2x2 trigger scintillator+PMT, recorded in ch8 of BOTH DRS4
//            groups) as (1) inter-group sync and (2) absolute timing reference.
// ----------------------------------------------------------------------------
// TR0 is split into ch8 of group 0 (chanOff(0,0,8)) and ch8 of group 1
// (chanOff(0,1,8)) -- the SAME PMT pulse in both groups (like the split MCP),
// ~315 mV, clean/fast. Two questions, from RAW waveforms, vs energy:
//   (1) SYNC: RMS(tr0_G0 - tr0_G1) = the inter-group (mezzanine) timing offset.
//             Per-event it can align group 1 to group 0 (fixes SW-U / cross-group).
//   (2) REFERENCE: is the absolute shower time referenced to TR0 much tighter
//       than to the MCP?  Recall SiPM-MCP ~ 3 ns is dominated by a 2.88 ns MCP
//       term. If SiPM-TR0 is tight (and improves with light), TR0 is the correct
//       absolute reference and the depth-independent shower time is unlocked.
//   ROOT_INCLUDE_PATH=lib/waveform:lib/io:lib/physics:lib/viz root -l -b -q analyze/studies/trSync.C+
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
#include "TH1F.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "FigPaths.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>

static const int kTR0a = 8192;          // chanOff(0,0,8)  (group 0 ch8)
static const int kTR0b = 1024*9+1024*8; // chanOff(0,1,8)  (group 1 ch8) = 17408

struct Ev { float slg,x,y; float sip,m1,tr0a,tr0b; float pm1,pa,pb; };

static double medOf(std::vector<float> v){ if(v.empty())return 0; size_t k=v.size()/2;
    std::nth_element(v.begin(),v.begin()+k,v.end()); return v[k]; }
struct Acc{double s=0,ss=0;long n=0;void add(double x){s+=x;ss+=x*x;++n;}
    double rms(){return n?std::sqrt(ss/n-(s/n)*(s/n)):0;} double mean(){return n?s/n:0;}};
static double corr(std::vector<Ev>&ev, float Ev::*a, float Ev::*b){
    double sa=0,sb=0,saa=0,sbb=0,sab=0; long n=0;
    for(auto&e:ev){ double X=e.*a,Y=e.*b; sa+=X;sb+=Y;saa+=X*X;sbb+=Y*Y;sab+=X*Y;++n; }
    double va=saa/n-(sa/n)*(sa/n), vb=sbb/n-(sb/n)*(sb/n), cov=sab/n-(sa/n)*(sb/n);
    return (va>0&&vb>0)? cov/std::sqrt(va*vb):0;
}
// bright-slice core sigma of (e.*p - e.*q) [ps]: top 45% by sum_lg, tight median
// core then tebSigma. Robust at the modest single-run statistics here.
static double brightCore(std::vector<Ev>&ev, float Ev::*p, float Ev::*q){
    if(ev.size()<800) return -1;
    std::vector<Ev> s=ev; std::sort(s.begin(),s.end(),[](const Ev&A,const Ev&B){return A.slg<B.slg;});
    size_t lo=(size_t)(0.55*s.size()); std::vector<float> t;
    for(size_t i=lo;i<s.size();++i) t.push_back(s[i].*p - s[i].*q);
    if(t.size()<300) return -1;
    double m=medOf(t); std::vector<float> c; for(float x:t) if(std::fabs(x-m)<1.5) c.push_back(x);
    if(c.size()<200) return -1; return rad::tebSigma(c);
}

static void processE(double E, const char* file, std::vector<Ev>& out){
    TChain ch("pulse"); ch.Add(radRaw(file)); TTreeReader rd(&ch);
    TTreeReaderArray<float> amp(rd,"amplitude"), tim(rd,"timevalue");
    while(rd.Next()){ const float* A=&amp[0]; const float* T=&tim[0];
        Pulse wcR=ExtractPulse(T+kWC_t,A+kWC_R,0.5f,kWC_minPeak), wcL=ExtractPulse(T+kWC_t,A+kWC_L,0.5f,kWC_minPeak),
              wcD=ExtractPulse(T+kWC_t,A+kWC_D,0.5f,kWC_minPeak), wcU=ExtractPulse(T+kWC_t,A+kWC_U,0.5f,kWC_minPeak);
        if(!(wcR.valid&&wcL.valid&&wcD.valid&&wcU.valid)) continue;
        PulseMulti m1=ExtractPulseMulti(T+kMCP1_t,A+kMCP1,20.f,30.f,2000.f);
        PulseMulti ta=ExtractPulseMulti(T+kT_D0G0,A+kTR0a,20.f,30.f,2000.f);
        PulseMulti tb=ExtractPulseMulti(T+kT_D0G1,A+kTR0b,20.f,30.f,2000.f);
        if(m1.peak<kMCP1_minPeak||m1.peak>kMCP1_maxPeak) continue;
        if(ta.peak<60.f||tb.peak<60.f) continue;
        if(m1.ledTime<-1e4||ta.ledTime<-1e4||tb.ledTime<-1e4) continue;
        // mean SiPM leading edge over the 6 pure-G0 corners (NW/NE/SE D+U) + sum_lg
        double sip=0; int ns=0,slgN=0; double slg=0;
        const int G0[6]={0,1,2,4,5,6};
        for(int j=0;j<6;++j){ int c=G0[j];
            PulseMulti hg=ExtractPulseMulti(T+kCap[c].hg_t,A+kCap[c].hg,20.f,5.f,(float)820.0);
            if(hg.peak>=kHG_minPeak&&hg.ledTime>-1e4){ sip+=hg.ledTime; ++ns; } }
        for(int c=0;c<8;++c){ Pulse lg=ExtractPulse(T+kCap[c].lg_t,A+kCap[c].lg,0.2f,5.f); slg+=lg.peak; }
        if(ns<4) continue;
        Ev e; e.slg=(float)slg; e.x=(float)(kWC_Scale*(wcR.peakTime-wcL.peakTime)); e.y=(float)(kWC_Scale*(wcD.peakTime-wcU.peakTime));
        e.sip=(float)(sip/ns); e.m1=m1.ledTime; e.tr0a=ta.ledTime; e.tr0b=tb.ledTime;
        e.pm1=m1.peak; e.pa=ta.peak; e.pb=tb.peak; out.push_back(e);
    }
    // fiducial
    double wx=0,wy=0,w=0; for(auto&e:out){wx+=e.x*e.slg;wy+=e.y*e.slg;w+=e.slg;}
    double xc=w>0?wx/w:0,yc=w>0?wy/w:0,rFid=TimingFiducialR(E),r2=rFid*rFid;
    std::vector<Ev> fid; for(auto&e:out){double dx=e.x-xc,dy=e.y-yc; if(dx*dx+dy*dy<r2) fid.push_back(e);} out.swap(fid);
}

void trSync(){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    double Es[6]={25,50,75,100,125,150}; double sMCP[6], sTR0[6], sInter[6];
    for(int e=0;e<6;++e){ std::vector<Ev> ev; processE(Es[e],Form("*_%.0f_GeV.root",Es[e]),ev);
        if(ev.size()<800){ printf("E=%3.0f  too few (%zu)\n",Es[e],ev.size()); sMCP[e]=sTR0[e]=sInter[e]=-1; continue; }
        Acc pa,pb,pm; for(auto&v:ev){pa.add(v.pa);pb.add(v.pb);pm.add(v.pm1);}
        double cAB=corr(ev,&Ev::pa,&Ev::pb);                      // TR0 same pulse in both groups?
        double cS_t=corr(ev,&Ev::sip,&Ev::tr0a), cS_m=corr(ev,&Ev::sip,&Ev::m1), cT_m=corr(ev,&Ev::tr0a,&Ev::m1);
        sInter[e]=brightCore(ev,&Ev::tr0a,&Ev::tr0b);                 // inter-group sync via TR0
        sMCP[e]  =brightCore(ev,&Ev::sip,&Ev::m1);                    // absolute SiPM time vs MCP
        sTR0[e]  =brightCore(ev,&Ev::sip,&Ev::tr0a);                  // absolute SiPM time vs TR0
        printf("E=%3.0f N=%6zu | TR0amp G0/G1=%.0f/%.0f corr=%.3f | corr(SiPM,TR0)=%.3f corr(SiPM,MCP)=%.3f corr(TR0,MCP)=%.3f\n",
               Es[e], ev.size(), pa.mean(), pb.mean(), cAB, cS_t, cS_m, cT_m);
        printf("        best-bin: inter-group(TR0a-TR0b)=%.1f ps | SiPM-MCP=%.1f ps | SiPM-TR0=%.1f ps\n",
               sInter[e], sMCP[e], sTR0[e]);
    }
    TCanvas* c=new TCanvas("c_tr","",900,680);
    TH1F* fr=c->DrawFrame(0,0,160,300);
    fr->SetTitle("Absolute shower time: TR0 reference vs MCP reference;beam energy (GeV);best-bin #sigma_{t} (ps)");
    TGraph* gM=new TGraph(6,Es,sMCP); gM->SetMarkerStyle(24); gM->SetMarkerColor(kGray+2); gM->SetLineColor(kGray+2); gM->SetMarkerSize(1.5); gM->SetLineWidth(2); gM->Draw("PL SAME");
    TGraph* gT=new TGraph(6,Es,sTR0); gT->SetMarkerStyle(20); gT->SetMarkerColor(kRed+1);  gT->SetLineColor(kRed+1);  gT->SetMarkerSize(1.7); gT->SetLineWidth(2); gT->Draw("PL SAME");
    TGraph* gI=new TGraph(6,Es,sInter);gI->SetMarkerStyle(21);gI->SetMarkerColor(kAzure+1);gI->SetLineColor(kAzure+1);gI->SetMarkerSize(1.5); gI->SetLineWidth(2); gI->Draw("PL SAME");
    TLegend* lg=new TLegend(0.40,0.70,0.88,0.88); lg->SetBorderSize(0);
    lg->AddEntry(gM,"SiPM #minus MCP (current reference)","pl");
    lg->AddEntry(gT,"SiPM #minus TR0 (scintillator reference)","pl");
    lg->AddEntry(gI,"inter-group sync: TR0(G0) #minus TR0(G1)","pl"); lg->Draw();
    gSystem->mkdir("figures",kTRUE); c->Print(radFigP("figures/tr0_sync.png"));
    printf("wrote figures/tr0_sync.png\n");
}

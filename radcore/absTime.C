// ============================================================================
// absTime.C — ABSOLUTE shower time vs the MCP, done right (leading-edge + walk).
// ----------------------------------------------------------------------------
// Diagnosis (see depthCorr.C / mcpJitterCanon.C): SiPM<->SiPM timing is coherent
// even across DRS4 groups (cross-group RMS = same-group RMS ~0.6 ns), so there is
// NO global-timebase problem.  The incoherence is SiPM-vs-MCP, and it collapses
// from 26 ns (CFD-5% of the clipped peak) to ~3 ns (leading edge).  It was the
// time EXTRACTION on slow, clipped LYSO pulses -- not the DRS4 domino phase.
//
// Here: build the depth-independent absolute shower time
//     B = mean_G0( t_led ) - mcp_time           (z and the G0 domino phase cancel)
// with a per-channel amplitude-walk correction t_led -> t_led - w(amp), and
// measure best-bin core sigma_t(E).  If B IMPROVES with energy it is light-limited
// -- the absolute-timing path the (DW-UP)/2 depth estimator can't give.
// Compare: A=(DW-UP)/2 (depth, led), B (no walk), Bw (walk-corrected).
//   ROOT_INCLUDE_PATH=radcore:Analysis root -l -b -q radcore/absTime.C+
// ============================================================================
#include "RadTiming.h"
#include "DataPaths.h"
#include "SelectionCuts.h"
#include "PlotUtils.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>

static const int G0[6]={0,1,2,4,5,6};   // NW/NE/SE Down+Up, all DRS0-G0, MCP1 ref

struct Ev { float slg, A, B, Bw, mt; float led[8], amp[8], lg[8]; };

static double medOf(std::vector<float> v){ if(v.empty())return 0; size_t k=v.size()/2;
    std::nth_element(v.begin(),v.begin()+k,v.end()); return v[k]; }

// best-bin core sigma of a per-event value; tight median window unwraps any tail
static double bestBin(std::vector<Ev>& ev, int which){  // 0=A,1=B,2=Bw
    std::vector<float> slg,t;
    for(auto&e:ev){ float v= which==0?e.A : (which==1?e.B:e.Bw); slg.push_back(e.slg); t.push_back(v);}
    if(slg.size()<2000) return -1;
    double smin=*std::min_element(slg.begin(),slg.end()), smax=*std::max_element(slg.begin(),slg.end());
    TH1F hS("hS","",150,smin,smax); hS.SetDirectory(0); for(float x:slg) hS.Fill(x);
    double muE,muEe,sigE,sigEe; FitGaussCore(&hS,2.0,muE,muEe,sigE,sigEe); if(sigE<=0){muE=hS.GetMean();sigE=hS.GetRMS();}
    double lo=muE-2*sigE, bw=4*sigE/9.0, best=1e9;
    for(int b=0;b<9;++b){ double blo=lo+b*bw,bhi=blo+bw; std::vector<float> vt;
        for(size_t i=0;i<slg.size();++i) if(slg[i]>=blo&&slg[i]<bhi) vt.push_back(t[i]);
        if(vt.size()<300) continue;
        if(which>0){ double m=medOf(vt); std::vector<float> c; for(float x:vt) if(std::fabs(x-m)<1.5) c.push_back(x); vt.swap(c); }
        if(vt.size()<200) continue;
        double s=rad::tebSigma(vt); if(s>3&&s<best) best=s; }
    return best;
}

void absTime(){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    const double Es[6]={25,50,75,100,125,150};
    double sA[6],sB[6],sBw[6];
    for(int e=0;e<6;++e){
        TFile* fp=TFile::Open(radReduced("DSB1",Es[e])); TTree* t=(TTree*)fp->Get("rad");
        Bool_t wc; Float_t x,y,mcp,mt,slg,led[8],hgp[8],lgp[8];
        t->SetBranchAddress("wc_ok",&wc); t->SetBranchAddress("x_trk",&x); t->SetBranchAddress("y_trk",&y);
        t->SetBranchAddress("mcp_peak",&mcp); t->SetBranchAddress("mcp_time",&mt); t->SetBranchAddress("sum_lg",&slg);
        t->SetBranchAddress("hg_led",led); t->SetBranchAddress("hg_peak",hgp); t->SetBranchAddress("lg_peak",lgp);
        long N=t->GetEntries(); double wx=0,wy=0,w=0;
        for(long i=0;i<N;++i){ t->GetEntry(i); if(wc&&mcp>kMCP_minPeak_E&&slg>kSumLG_centroid){wx+=x*slg;wy+=y*slg;w+=slg;} }
        double xc=w>0?wx/w:0,yc=w>0?wy/w:0,rFid=TimingFiducialR(Es[e]),r2=rFid*rFid;

        std::vector<Ev> ev; ev.reserve(N);
        for(long i=0;i<N;++i){ t->GetEntry(i);
            if(!wc||mcp<kMCP1_minPeak||mcp>kMCP1_maxPeak||mt<-1e4) continue;
            double dx=x-xc,dy=y-yc; if(dx*dx+dy*dy>=r2) continue;
            Ev e; e.slg=slg; e.mt=mt; bool ok=true;
            for(int c=0;c<8;++c){ e.led[c]=led[c]; e.amp[c]=hgp[c]; e.lg[c]=lgp[c]; }
            double ds=0,us=0; int dn=0,un=0;
            for(int j=0;j<3;++j){int c=G0[j]; if(hgp[c]>=kHG_minPeak&&led[c]>-1e4){ds+=led[c];++dn;} else if(hgp[c]<kHG_minPeak){} }
            for(int j=3;j<6;++j){int c=G0[j]; if(hgp[c]>=kHG_minPeak&&led[c]>-1e4){us+=led[c];++un;} }
            if(dn<1||un<1){ok=false;}
            if(!ok) continue;
            double td=ds/dn, tu=us/un;
            e.A=0.5f*(float)(td-tu); e.B=0.5f*(float)(td+tu)-mt; e.Bw=e.B;
            ev.push_back(e);
        }
        // per-channel time-walk correction of the 20 mV leading edge: walk ~ thresh/slope ~ 1/amplitude.
        // Use the UNSATURATED LG amplitude (hg_peak is clipped at high E -> useless as a walk variable).
        for(int c=0;c<8;++c){
            double sx=0,sy=0,sxx=0,sxy=0; long n=0;
            for(auto&e:ev){ if(e.amp[c]<kHG_minPeak||e.led[c]<-1e4||e.lg[c]<10.f) continue; double X=1.0/e.lg[c]; sx+=X;sy+=e.led[c];sxx+=X*X;sxy+=X*e.led[c];++n; }
            if(n<200) continue; double d=n*sxx-sx*sx; if(std::fabs(d)<1e-12) continue;
            double b=(n*sxy-sx*sy)/d, a=(sy-b*sx)/n;
            for(auto&e:ev){ if(e.amp[c]<kHG_minPeak||e.led[c]<-1e4||e.lg[c]<10.f) continue; double X=1.0/e.lg[c]; e.led[c]-= (float)(a+b*X); }
        }
        // recompute Bw with walk-corrected led
        for(auto&e:ev){ double ds=0,us=0; int dn=0,un=0;
            for(int j=0;j<3;++j){int c=G0[j]; if(e.amp[c]>=kHG_minPeak&&e.led[c]>-1e4){ds+=e.led[c];++dn;} }
            for(int j=3;j<6;++j){int c=G0[j]; if(e.amp[c]>=kHG_minPeak&&e.led[c]>-1e4){us+=e.led[c];++un;} }
            if(dn>=1&&un>=1) e.Bw=0.5f*(float)((ds/dn+us/un))-e.mt; }

        sA[e]=bestBin(ev,0); sB[e]=bestBin(ev,1); sBw[e]=bestBin(ev,2);
        printf("E=%3.0f  A=(DW-UP)/2 led depth=%.1f ps | B=(DW+UP)/2-MCP led=%.1f ps | Bw walk-corrected=%.1f ps  (N=%zu)\n",
               Es[e], sA[e], sB[e], sBw[e], ev.size());
        fp->Close();
    }
    TCanvas* c=new TCanvas("c_at","",900,680);
    TH1F* fr=c->DrawFrame(0,0,160,200);
    fr->SetTitle("Absolute shower time vs MCP (leading-edge + walk) -- does it improve with light?;beam energy (GeV);best-bin #sigma_{t} (ps)");
    TGraph* gA=new TGraph(6,Es,sA); gA->SetMarkerStyle(20); gA->SetMarkerColor(kAzure+1); gA->SetLineColor(kAzure+1); gA->SetMarkerSize(1.6); gA->SetLineWidth(2); gA->Draw("PL SAME");
    TGraph* gB=new TGraph(6,Es,sB); gB->SetMarkerStyle(24); gB->SetMarkerColor(kGray+2); gB->SetLineColor(kGray+2); gB->SetMarkerSize(1.5); gB->SetLineWidth(2); gB->Draw("PL SAME");
    TGraph* gW=new TGraph(6,Es,sBw);gW->SetMarkerStyle(21); gW->SetMarkerColor(kRed+1); gW->SetLineColor(kRed+1); gW->SetMarkerSize(1.6); gW->SetLineWidth(2); gW->Draw("PL SAME");
    TLegend* lg=new TLegend(0.40,0.70,0.88,0.88); lg->SetBorderSize(0);
    lg->AddEntry(gA,"(DW#minusUP)/2 depth (led)","pl");
    lg->AddEntry(gB,"(DW+UP)/2 #minus MCP (led, raw)","pl");
    lg->AddEntry(gW,"(DW+UP)/2 #minus MCP (led + walk)","pl"); lg->Draw();
    gSystem->mkdir("radcore/figs",kTRUE); c->Print("radcore/figs/abs_time.png");
    printf("wrote radcore/figs/abs_time.png\n");
}

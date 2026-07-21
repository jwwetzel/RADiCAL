// ============================================================================
// cfdInterp.C — FULL-STATISTICS test (no re-reduction) of whether a better, FIXED
// absolute timing threshold recovers sigma_t improvement with energy.
// ----------------------------------------------------------------------------
// The reduced data stores 6 CFD crossing times (3/5/10/20/30/50% of the clipped
// peak). The crossing-vs-fraction curve is smooth, so we interpolate the crossing
// at any FIXED absolute level L (mV) -> f = L/hg_peak in [5%,50%] -> bracketed
// interpolation. This uses a fixed per-channel level (NOT the per-event LG, which
// would inject ~tens of ps of LG photostatistics noise). Compare best-bin
// (DW-UP)/2 sigma_t(E) for standard CFD-5% vs the best fixed-L, on full stats.
//   ROOT_INCLUDE_PATH=lib/waveform:lib/io:lib/physics:lib/viz root -l -b -q analyze/studies/cfdInterp.C+
// ============================================================================
#include "RadTiming.h"       // rad::tebSigma
#include "SelectionCuts.h"
#include "DataPaths.h"
#include "PlotUtils.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "FigPaths.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>

static const double FR[6]={0.03,0.05,0.10,0.20,0.30,0.50};
// interpolate the crossing time at fraction f from the 6 stored crossings
static double interp(const float c[6], double f){
    if(f<FR[0]||f>FR[5]) return -1e9;
    for(int k=0;k<5;++k) if(f>=FR[k] && f<=FR[k+1]){ if(c[k]<-1e5||c[k+1]<-1e5) return -1e9;
        return c[k]+(c[k+1]-c[k])*(f-FR[k])/(FR[k+1]-FR[k]); }
    return -1e9;
}
struct Ev { float slg; float hgp[8]; float cfd[8][6]; };

static double bestBin(const std::vector<Ev>& ev, double L /*mV, <=0 = plain CFD-5%*/) {
    // build (sum_lg, (DW-UP)/2) then 9-bin best-bin core sigma
    std::vector<float> slg, tt;
    for(const auto&e:ev){ double d=0,u=0; int nd=0,nu=0;
        for(int c=0;c<8;++c){ if(e.hgp[c]<kHG_minPeak) continue;
            double tc; if(L<=0) tc=e.cfd[c][1];          // stored CFD-5%
                       else { double f=L/e.hgp[c]; tc=interp(e.cfd[c],f); }
            if(tc<-1e5) continue; if(c<4){d+=tc;++nd;}else{u+=tc;++nu;} }
        if(nd>=1&&nu>=1){ slg.push_back(e.slg); tt.push_back(0.5f*(float)(d/nd-u/nu)); } }
    if(slg.size()<2000) return -1;
    double smin=*std::min_element(slg.begin(),slg.end()), smax=*std::max_element(slg.begin(),slg.end());
    TH1F hS("hS","",150,smin,smax); hS.SetDirectory(0); for(float x:slg) hS.Fill(x);
    double muE,muEe,sigE,sigEe; FitGaussCore(&hS,2.0,muE,muEe,sigE,sigEe); if(sigE<=0){muE=hS.GetMean();sigE=hS.GetRMS();}
    double lo=muE-2*sigE, bw=4*sigE/9.0, best=1e9;
    for(int b=0;b<9;++b){ double blo=lo+b*bw,bhi=blo+bw; std::vector<float> vt;
        for(size_t i=0;i<slg.size();++i) if(slg[i]>=blo&&slg[i]<bhi) vt.push_back(tt[i]);
        if(vt.size()<400) continue; double s=rad::tebSigma(vt); if(s>5&&s<best) best=s; }
    return best;
}

static void loadE(double E, std::vector<Ev>& out) {
    TFile* fp=TFile::Open(radReduced("DSB1",E)); TTree* t=(TTree*)fp->Get("rad");
    Bool_t wc; Float_t x,y,mcp,slg,hgp[8],c03[8],c05[8],c10[8],c20[8],c30[8],c50[8];
    t->SetBranchAddress("wc_ok",&wc); t->SetBranchAddress("x_trk",&x); t->SetBranchAddress("y_trk",&y);
    t->SetBranchAddress(t->GetBranch("mcp1_peak")?"mcp1_peak":"mcp_peak",&mcp); t->SetBranchAddress("sum_lg",&slg); t->SetBranchAddress("hg_peak",hgp);
    t->SetBranchAddress("hg_cfd03",c03); t->SetBranchAddress("hg_cfd05",c05); t->SetBranchAddress("hg_cfd10",c10);
    t->SetBranchAddress("hg_cfd",c20); t->SetBranchAddress("hg_cfd30",c30); t->SetBranchAddress("hg_cfd50",c50);
    long N=t->GetEntries(); double wx=0,wy=0,w=0;
    for(long i=0;i<N;++i){ t->GetEntry(i); if(wc&&mcp>kMCP_minPeak_E&&slg>kSumLG_centroid){wx+=x*slg;wy+=y*slg;w+=slg;} }
    double xc=w>0?wx/w:0,yc=w>0?wy/w:0,rFid=TimingFiducialR(E),r2=rFid*rFid;
    for(long i=0;i<N;++i){ t->GetEntry(i);
        if(!wc||mcp<kMCP1_minPeak||mcp>kMCP1_maxPeak) continue; double dx=x-xc,dy=y-yc; if(dx*dx+dy*dy>=r2) continue;
        Ev e; e.slg=slg; for(int c=0;c<8;++c){ e.hgp[c]=hgp[c];
            e.cfd[c][0]=c03[c];e.cfd[c][1]=c05[c];e.cfd[c][2]=c10[c];e.cfd[c][3]=c20[c];e.cfd[c][4]=c30[c];e.cfd[c][5]=c50[c]; }
        out.push_back(e); }
    fp->Close();
}

void cfdInterp() {
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    double Es[6]={25,50,75,100,125,150};
    double Ls[5]={120,160,200,240,300};
    double s5[6], sBest[6]; double Lbest[6];
    for(int e=0;e<6;++e){ std::vector<Ev> ev; loadE(Es[e],ev);
        s5[e]=bestBin(ev,-1);                              // standard CFD-5%
        sBest[e]=1e9; Lbest[e]=0;
        printf("E=%3.0f  CFD-5%%=%.1f  fixed-L:", Es[e], s5[e]);
        for(double L:Ls){ double s=bestBin(ev,L); printf(" %.0f=%.1f", L, s);
            if(s>0&&s<sBest[e]){sBest[e]=s;Lbest[e]=L;} }
        printf("  -> best L=%.0f mV: %.1f ps\n", Lbest[e], sBest[e]); }

    TCanvas* c=new TCanvas("c_ci","",900,680); c->SetGridx(); c->SetGridy();
    TH1F* fr=c->DrawFrame(0,15,160,55);
    fr->SetTitle("Full-stat timing: best fixed-absolute threshold vs CFD-5%;beam energy [GeV];best-bin (DW#minusUP)/2 #sigma_{t} [ps]");
    TGraph* g5=new TGraph(6,Es,s5);   g5->SetMarkerStyle(24); g5->SetMarkerColor(kGray+2); g5->SetLineColor(kGray+2); g5->SetMarkerSize(1.5); g5->SetLineWidth(2); g5->Draw("PL SAME");
    TGraph* gB=new TGraph(6,Es,sBest);gB->SetMarkerStyle(20); gB->SetMarkerColor(kRed+1);  gB->SetLineColor(kRed+1);  gB->SetMarkerSize(1.6); gB->SetLineWidth(2); gB->Draw("PL SAME");
    TLegend* lg=new TLegend(0.40,0.74,0.88,0.88); lg->SetBorderSize(0);
    lg->AddEntry(g5,"standard CFD-5% (of clipped peak)","pl");
    lg->AddEntry(gB,"best fixed-absolute threshold (per E)","pl"); lg->Draw();
    gSystem->mkdir("figures",kTRUE); c->Print(radFigP("figures/cfd_interp.png"));
    printf("wrote figures/cfd_interp.png\n");
}

// ============================================================================
// threshScan.C — does timing HIGHER on the edge improve sigma_t, and make it
//                track energy?  (reduced DSB1; hg_cfd* are MCP-referenced.)
// ----------------------------------------------------------------------------
// Two observables vs threshold (cfd fraction of the clipped peak) vs energy:
//   * SLEW PROXY: sigma(hg[NW-D] - hg[NE-D])/sqrt2  -- two G0 down channels, same
//     shower; MCP and shower-arrival cancel, leaving slew (+) electronics. Its
//     THRESHOLD dependence isolates slew.
//   * ABSOLUTE: sigma(mean hg over G0)  -- the absolute shower time.
// If higher fractions give smaller sigma AND a steeper improvement with energy,
// the foot (cfd05) is the limitation, as argued.
//   ROOT_INCLUDE_PATH=radcore:Analysis root -l -b -q radcore/threshScan.C+
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

static double medOf(std::vector<float> v){ if(v.empty())return 0; size_t k=v.size()/2;
    std::nth_element(v.begin(),v.begin()+k,v.end()); return v[k]; }
// bright-slice (top 45% sum_lg) tight-core sigma of a value list
static double bc(std::vector<std::pair<float,float>>& v){
    if(v.size()<600) return -1; std::sort(v.begin(),v.end());
    std::vector<float> t; for(size_t i=0.55*v.size();i<v.size();++i) t.push_back(v[i].second);
    double m=medOf(t); std::vector<float> c; for(float x:t) if(std::fabs(x-m)<1.0) c.push_back(x);
    if(c.size()<300) return -1; return rad::tebSigma(c);
}

void threshScan(){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    const double Es[6]={25,50,75,100,125,150};
    const char* frac[6]={"hg_cfd05","hg_cfd10","hg_cfd","hg_cfd30","hg_cfd50","hg_led"};
    const char* flab[6]={"cfd05(foot)","cfd10","cfd20","cfd30","cfd50(steep)","led20mV"};
    double slew[6][6], absT[6][6];   // [frac][energy]
    for(int e=0;e<6;++e){
        TFile* fp=TFile::Open(radReduced("DSB1",Es[e])); TTree* t=(TTree*)fp->Get("rad");
        Bool_t wc; Float_t x,y,mcp,slg,hgp[8]; Float_t v[6][8];
        t->SetBranchAddress("wc_ok",&wc); t->SetBranchAddress("x_trk",&x); t->SetBranchAddress("y_trk",&y);
        t->SetBranchAddress("mcp_peak",&mcp); t->SetBranchAddress("sum_lg",&slg); t->SetBranchAddress("hg_peak",hgp);
        for(int k=0;k<6;++k) t->SetBranchAddress(frac[k], v[k]);
        long N=t->GetEntries(); double wx=0,wy=0,w=0;
        for(long i=0;i<N;++i){ t->GetEntry(i); if(wc&&mcp>kMCP_minPeak_E&&slg>kSumLG_centroid){wx+=x*slg;wy+=y*slg;w+=slg;} }
        double xc=w>0?wx/w:0,yc=w>0?wy/w:0,rFid=TimingFiducialR(Es[e]),r2=rFid*rFid;
        std::vector<std::pair<float,float>> sl[6], ab[6];
        for(long i=0;i<N;++i){ t->GetEntry(i);
            if(!wc||mcp<kMCP1_minPeak||mcp>kMCP1_maxPeak) continue;
            double dx=x-xc,dy=y-yc; if(dx*dx+dy*dy>=r2) continue;
            for(int k=0;k<6;++k){
                if(hgp[0]>=kHG_minPeak&&hgp[1]>=kHG_minPeak&&v[k][0]>-1e4&&v[k][1]>-1e4)
                    sl[k].push_back({slg,0.70710678f*(v[k][0]-v[k][1])});   // /sqrt2
                double s=0;int n=0; for(int c=0;c<7;++c) if(hgp[c]>=kHG_minPeak&&v[k][c]>-1e4){s+=v[k][c];++n;}
                if(n>=4) ab[k].push_back({slg,(float)(s/n)});
            }
        }
        printf("E=%3.0f  slew-proxy sigma(NW-D - NE-D)/sqrt2 [ps]:", Es[e]);
        for(int k=0;k<6;++k){ slew[k][e]=bc(sl[k]); printf("  %s=%.0f", flab[k], slew[k][e]); }
        printf("\n        absolute mean(hg) [ps]:                   ");
        for(int k=0;k<6;++k){ absT[k][e]=bc(ab[k]); printf("  %s=%.0f", flab[k], absT[k][e]); }
        printf("\n");
        fp->Close();
    }
    // plot slew proxy vs energy for cfd05 (foot) vs cfd50 (steep)
    TCanvas* c=new TCanvas("c_ts","",900,680); c->SetGridx(); c->SetGridy();
    TH1F* fr=c->DrawFrame(0,0,160,60);
    fr->SetTitle("Single-channel slew vs energy: foot (cfd05) vs steep (cfd50);beam energy (GeV);#sigma(NW-D #minus NE-D)/#sqrt{2} (ps)");
    double y05[6],y50[6]; for(int e=0;e<6;++e){y05[e]=slew[0][e];y50[e]=slew[4][e];}
    TGraph* g0=new TGraph(6,Es,y05); g0->SetMarkerStyle(24); g0->SetMarkerColor(kGray+2); g0->SetLineColor(kGray+2); g0->SetMarkerSize(1.6); g0->SetLineWidth(2); g0->Draw("PL SAME");
    TGraph* g5=new TGraph(6,Es,y50); g5->SetMarkerStyle(20); g5->SetMarkerColor(kRed+1);  g5->SetLineColor(kRed+1);  g5->SetMarkerSize(1.6); g5->SetLineWidth(2); g5->Draw("PL SAME");
    TLegend* lg=new TLegend(0.40,0.72,0.88,0.88); lg->SetBorderSize(0);
    lg->AddEntry(g0,"cfd05 (foot, where we time now)","pl");
    lg->AddEntry(g5,"cfd50 (steep part of the edge)","pl"); lg->Draw();
    gSystem->mkdir("radcore/figs",kTRUE); c->Print("radcore/figs/thresh_scan.png");
    printf("wrote radcore/figs/thresh_scan.png\n");
}

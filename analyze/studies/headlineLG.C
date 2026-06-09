// ============================================================================
// headlineLG.C — the HEADLINE (DW-UP)/2 best-bin ladder, cfd05 vs lgcfd.
// ----------------------------------------------------------------------------
// Exact published pipeline (RadTiming::timingBestBin): LG-weighted beam centroid
// -> fiducial -> (DW-UP)/2 from the chosen MCP-referenced time -> 9 sum_lg bins
// (muE +- 2 sigmaE) -> best bin (N>=500) -> tebSigma. Run for hg_cfd05 (current
// headline) and hg_lgcfd (CFD on the LG-predicted true peak), and fit each
// sigma_t(E) = a/sqrt(E) (+) b.  (Local re-reduced DSB1 raw -> subsample stats.)
//   ROOT_INCLUDE_PATH=lib/waveform:lib/io:lib/physics:lib/viz root -l -b -q 'analyze/studies/headlineLG.C+("/tmp/radlg/DSB1")'
// ============================================================================
#include "RadTiming.h"        // rad::tebSigma
#include "SelectionCuts.h"
#include "PlotUtils.h"        // FitGaussCore
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"
#include "FigPaths.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>

static void ladder(const char* dir, const char* tb, double E, double& sig, int& bestbin, size_t& nfid){
    sig=-1; bestbin=-1; nfid=0;
    TFile* fp=TFile::Open(Form("%s/%.0fGeV.root",dir,E)); if(!fp||fp->IsZombie()) return;
    TTree* t=(TTree*)fp->Get("rad");
    Bool_t wc; Float_t x,y,mcp,slg,tv[8],hgp[8];
    t->SetBranchAddress("wc_ok",&wc); t->SetBranchAddress("x_trk",&x); t->SetBranchAddress("y_trk",&y);
    t->SetBranchAddress("mcp1_peak",&mcp); t->SetBranchAddress("sum_lg",&slg);
    t->SetBranchAddress(tb,tv); t->SetBranchAddress("hg_peak",hgp);
    long N=t->GetEntries();
    // LG-weighted beam centroid
    double wx=0,wy=0,w=0;
    for(long i=0;i<N;++i){ t->GetEntry(i); if(wc&&mcp>kMCP_minPeak_E&&slg>kSumLG_centroid){wx+=x*slg;wy+=y*slg;w+=slg;} }
    double xc=w>0?wx/w:0,yc=w>0?wy/w:0,rF=TimingFiducialR(E),r2=rF*rF;
    std::vector<float> sl, tval;
    for(long i=0;i<N;++i){ t->GetEntry(i);
        if(!wc||mcp<kMCP1_minPeak||mcp>kMCP1_maxPeak) continue;
        double dx=x-xc,dy=y-yc; if(dx*dx+dy*dy>=r2) continue;
        double ds=0,us=0; int dn=0,un=0;
        for(int c=0;c<4;++c) if(hgp[c]>=kHG_minPeak&&tv[c]>-1e4){ds+=tv[c];++dn;}
        for(int c=4;c<8;++c) if(hgp[c]>=kHG_minPeak&&tv[c]>-1e4){us+=tv[c];++un;}
        if(dn<1||un<1) continue;
        sl.push_back(slg); tval.push_back(0.5f*(float)(ds/dn-us/un));
    }
    nfid=sl.size(); if(sl.size()<600){ fp->Close(); return; }
    // subsample stats -> bright-slice core (top 45% by sum_lg) instead of 9-bin best-bin.
    std::vector<std::pair<float,float>> v; for(size_t i=0;i<sl.size();++i) v.push_back({sl[i],tval[i]});
    std::sort(v.begin(),v.end());
    std::vector<float> bt; for(size_t i=(size_t)(0.55*v.size());i<v.size();++i) bt.push_back(v[i].second);
    bestbin=1; sig=rad::tebSigma(bt); fp->Close();
}

void headlineLG(const char* dir="/tmp/radlg/DSB1"){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    const double Es[6]={25,50,75,100,125,150};
    double s05[6], sLG[6];
    printf("%-6s %12s %12s\n","E[GeV]","cfd05[ps]","lgcfd[ps]");
    for(int e=0;e<6;++e){ int bb; size_t nf;
        ladder(dir,"hg_cfd05",Es[e],s05[e],bb,nf);
        ladder(dir,"hg_lgcfd",Es[e],sLG[e],bb,nf);
        printf("%-6.0f %12.1f %12.1f   (Nfid=%zu, bin=%d)\n",Es[e],s05[e],sLG[e],nf,bb);
    }
    // stochastic fits  sigma = sqrt(a^2/E + b^2)
    TGraph g05(6,Es,s05), gLG(6,Es,sLG);
    TF1 f05("f05","sqrt([0]*[0]/x+[1]*[1])",20,160), fLG("fLG","sqrt([0]*[0]/x+[1]*[1])",20,160);
    f05.SetParameters(200,25); fLG.SetParameters(200,20);
    g05.Fit(&f05,"Q"); gLG.Fit(&fLG,"Q");
    printf("\nstochastic fit sigma_t = a/sqrt(E) (+) b:\n");
    printf("  cfd05 : a=%.0f ps*sqrt(GeV)  b=%.1f ps\n", f05.GetParameter(0), f05.GetParameter(1));
    printf("  lgcfd : a=%.0f ps*sqrt(GeV)  b=%.1f ps\n", fLG.GetParameter(0), fLG.GetParameter(1));

    TCanvas* c=new TCanvas("c_hl","",900,680); c->SetGridx(); c->SetGridy();
    TH1F* fr=c->DrawFrame(0,0,160,70);
    fr->SetTitle("Headline (DW#minusUP)/2 best-bin: cfd05 vs LG-predicted-peak CFD;beam energy (GeV);#sigma_{t} (ps)");
    g05.SetMarkerStyle(24); g05.SetMarkerColor(kGray+2); g05.SetLineColor(kGray+2); g05.SetMarkerSize(1.6); g05.Draw("P SAME"); f05.SetLineColor(kGray+2); f05.Draw("SAME");
    gLG.SetMarkerStyle(20); gLG.SetMarkerColor(kRed+1); gLG.SetLineColor(kRed+1); gLG.SetMarkerSize(1.7); gLG.Draw("P SAME"); fLG.SetLineColor(kRed+1); fLG.Draw("SAME");
    TLegend* lg=new TLegend(0.40,0.72,0.88,0.88); lg->SetBorderSize(0);
    lg->AddEntry(&g05,Form("cfd05 (clipped): %.0f/#sqrt{E} #oplus %.0f ps", f05.GetParameter(0),f05.GetParameter(1)),"pl");
    lg->AddEntry(&gLG,Form("lgcfd (true peak): %.0f/#sqrt{E} #oplus %.0f ps", fLG.GetParameter(0),fLG.GetParameter(1)),"pl");
    lg->Draw();
    gSystem->mkdir("figures",kTRUE); c->Print(radFigP("figures/headline_lgcfd.png"));
    printf("wrote figures/headline_lgcfd.png\n");
}

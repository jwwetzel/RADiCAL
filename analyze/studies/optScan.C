// ============================================================================
// optScan.C — expose EVERY value the best-bin estimator selects, at EVERY energy.
// timingBestBin slices the fiducial sample into 9 sum_lg bins (muE +- 2 sigmaE)
// and keeps the lowest-sigma bin with N>=500. This dumps, per energy, the full
// landscape: muE, sigmaE, and for each of the 9 bins its sum_lg center, N, and
// (DW-UP)/2 sigma_t for cfd05 AND lgcfd — so the trend (and where 150 GeV
// deviates) is plain. 6 panels (one per energy) + a console table.
//   source setup.sh; root -l -b -q 'analyze/studies/optScan.C+("DSB1")'
// ============================================================================
#include "RadView.h"
#include "RadTiming.h"
#include "SelectionCuts.h"
#include "DataPaths.h"
#include "PlotUtils.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TMarker.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>
using namespace rad;

static bool depth(RadView& v,double xc,double yc,double r2,int src,float& out){
    if(!v.wc_ok()||v.mcp1_peak()<kMCP1_minPeak||v.mcp1_peak()>kMCP1_maxPeak) return false;
    double dx=v.x_trk()-xc, dy=v.y_trk()-yc; if(dx*dx+dy*dy>=r2) return false;
    double ds=0,us=0;int dn=0,un=0;
    for(int c=0;c<4;++c) if(v.is_timing(c)&&v.hg_peak(c)>=kHG_minPeak){float tc=v.timeOf(c,src);if(tc>-1e5){ds+=tc;++dn;}}
    for(int c=4;c<8;++c) if(v.is_timing(c)&&v.hg_peak(c)>=kHG_minPeak){float tc=v.timeOf(c,src);if(tc>-1e5){us+=tc;++un;}}
    if(dn<1||un<1)return false; out=0.5f*(float)(ds/dn-us/un); return true;
}

void optScan(const char* build="DSB1"){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    BuildConfig cfg=BuildConfig::Load(radConfig(build).Data());
    const double Es[]={25,50,75,100,125,150}; int nE=6;
    TCanvas* c=new TCanvas("c_opt","",1500,900); c->Divide(3,2,0.008,0.03);
    printf("\n%s  best-bin optimization landscape (9 sum_lg bins, mu_E +- 2 sigma_E):\n",build);
    for(int e=0;e<nE;++e){ double E=Es[e];
        TFile* fp=TFile::Open(radReduced(build,E)); if(!fp||fp->IsZombie()){continue;}
        TTree* t=(TTree*)fp->Get("rad"); RadView v; v.attach(t,&cfg);
        double xc,yc; v.beamCenter(xc,yc); double r2=9.0; Long64_t N=v.entries();
        // gather fiducial sum_lg + both depth series
        std::vector<float> slg, d5, dL; std::vector<char> ok5, okL;
        for(Long64_t i=0;i<N;++i){ v.get(i);
            if(!v.wc_ok()||v.mcp1_peak()<kMCP1_minPeak||v.mcp1_peak()>kMCP1_maxPeak) continue;
            double dx=v.x_trk()-xc, dy=v.y_trk()-yc; if(dx*dx+dy*dy>=r2) continue;
            float a5,aL; bool b5=depth(v,xc,yc,r2,RadView::kCFD05,a5), bL=depth(v,xc,yc,r2,RadView::kLGCFD,aL);
            slg.push_back(v.sum_lg()); d5.push_back(b5?a5:-1e9); dL.push_back(bL?aL:-1e9); }
        double smin=*std::min_element(slg.begin(),slg.end()), smax=*std::max_element(slg.begin(),slg.end());
        TH1F h("_h","",150,smin,smax); for(float x:slg)h.Fill(x); double mu,sig,me,se; FitGaussCore(&h,2.0,mu,me,sig,se);
        if(sig<=0){mu=h.GetMean();sig=h.GetRMS();}
        double binLo=mu-2*sig, binW=4*sig/9.0;
        std::vector<double> bc, s5, sL; double best5=1e9,bestL=1e9; int ib5=-1,ibL=-1;
        printf("  E=%3.0f  muE=%.0f sigE=%.0f | bin: sumlg   N    cfd05  lgcfd\n",E,mu,sig);
        for(int b=0;b<9;++b){ double blo=binLo+b*binW,bhi=blo+binW; std::vector<float> v5,vL;
            for(size_t i=0;i<slg.size();++i) if(slg[i]>=blo&&slg[i]<bhi){ if(d5[i]>-1e5)v5.push_back(d5[i]); if(dL[i]>-1e5)vL.push_back(dL[i]); }
            double c5=(v5.size()>=500)?tebSigma(v5):-1, cL=(vL.size()>=500)?tebSigma(vL):-1;
            double xc2=0.5*(blo+bhi); bc.push_back(xc2); s5.push_back(c5); sL.push_back(cL);
            if(c5>15&&c5<best5){best5=c5;ib5=b;} if(cL>15&&cL<bestL){bestL=cL;ibL=b;}
            printf("        %2d: %5.0f %6zu  %5.1f  %5.1f\n",b,xc2,vL.size(),c5,cL); }
        printf("    -> best cfd05=%.1f (bin %d) ; best lgcfd=%.1f (bin %d)\n",best5,ib5,bestL,ibL);
        // panel
        c->cd(e+1); gPad->SetGridy(); gPad->SetLeftMargin(0.13); gPad->SetRightMargin(0.04); gPad->SetTopMargin(0.10);
        std::vector<double> bx5,by5,bxL,byL; for(int b=0;b<9;++b){ if(s5[b]>0){bx5.push_back(bc[b]);by5.push_back(s5[b]);} if(sL[b]>0){bxL.push_back(bc[b]);byL.push_back(sL[b]);} }
        TGraph* g5=new TGraph(bx5.size(),&bx5[0],&by5[0]); g5->SetMarkerStyle(24); g5->SetMarkerColor(kRed+1); g5->SetLineColor(kRed+1); g5->SetMarkerSize(1.3);
        g5->SetTitle(Form("%.0f GeV;#Sigma LG (bin center);(DW#minusUP)/2 #sigma_{t} (ps)",E));
        g5->GetYaxis()->SetRangeUser(18,46); g5->Draw("ALP");
        TGraph* gL=new TGraph(bxL.size(),&bxL[0],&byL[0]); gL->SetMarkerStyle(20); gL->SetMarkerColor(kAzure+2); gL->SetLineColor(kAzure+2); gL->SetMarkerSize(1.3); gL->SetLineWidth(2); gL->Draw("LP SAME");
        if(ibL>=0){ TMarker* m=new TMarker(bc[ibL],sL[ibL],29); m->SetMarkerColor(kAzure+3); m->SetMarkerSize(2.4); m->Draw(); }
        TLatex tl; tl.SetNDC(); tl.SetTextSize(0.058); tl.SetTextColor(kAzure+3);
        tl.DrawLatex(0.16,0.16,Form("best lgcfd = %.1f ps",bestL));
        fp->Close();
    }
    c->cd(0); TLatex t0; t0.SetNDC(); t0.SetTextFont(62); t0.SetTextSize(0.022);
    t0.DrawLatex(0.18,0.985,Form("%s  (DW#minusUP)/2 #sigma_{t} per sum_lg bin, every energy   (#circ cfd05  #bullet lgcfd, #star = selected best bin)",build));
    gSystem->mkdir("figures/narrative",kTRUE); c->Print("figures/narrative/opt_scan.png");
    printf("  wrote figures/narrative/opt_scan.png\n");
}

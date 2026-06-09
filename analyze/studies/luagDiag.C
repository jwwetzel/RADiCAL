// ============================================================================
// luagDiag.C — why is LUAG non-monotonic? Diagnose how each (E) point is built.
// Per energy: Nfid, sum_lg mu/sigma, HG amplitude + clip fraction, and the FULL
// sigma_t-vs-(quantile-bin) landscape for cfd05 AND lgcfd with per-bin N. Tells
// us whether 75<100 is a real per-bin effect (e.g. clipping onset) or a binning
// /selection artifact, and whether a source switch is needed.
//   source setup.sh; root -l -b -q 'analyze/studies/luagDiag.C+("LUAG")'
// ============================================================================
#include "RadView.h"
#include "RadTiming.h"
#include "SelectionCuts.h"
#include "DataPaths.h"
#include "PlotUtils.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TMarker.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>
using namespace rad;

struct E4 { float slg, d5, dL, hgmax; };

void luagDiag(const char* build="LUAG"){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    BuildConfig cfg=BuildConfig::Load(radConfig(build).Data());
    const double Es[]={50,75,100,125,150}; int nE=5;
    TCanvas* c=new TCanvas("cl","",1500,560); c->Divide(5,1,0.006,0.02);
    printf("\n=== %s diagnostic: per-energy build of each sigma_t point ===\n",build);
    for(int e=0;e<nE;++e){ double E=Es[e];
        TFile* fp=TFile::Open(radReduced(build,E)); if(!fp||fp->IsZombie()){printf("  %.0f: no file\n",E);continue;}
        TTree* t=(TTree*)fp->Get("rad"); RadView v; v.attach(t,&cfg);
        double xc,yc; v.beamCenter(xc,yc); double r2=9.0; Long64_t N=v.entries();
        std::vector<E4> ev; long nclip=0,nhg=0;
        for(Long64_t i=0;i<N;++i){ v.get(i);
            if(!v.wc_ok()||v.mcp1_peak()<kMCP1_minPeak||v.mcp1_peak()>kMCP1_maxPeak) continue;
            double dx=v.x_trk()-xc,dy=v.y_trk()-yc; if(dx*dx+dy*dy>=r2) continue;
            double d5d=0,d5u=0,dLd=0,dLu=0;int e5d=0,e5u=0,eLd=0,eLu=0; float hgmax=0;
            for(int ch=0;ch<4;++ch) if(v.is_timing(ch)&&v.hg_peak(ch)>=kHG_minPeak){float t5=v.timeOf(ch,RadView::kCFD05),tl=v.timeOf(ch,RadView::kLGCFD);if(t5>-1e5){d5d+=t5;++e5d;}if(tl>-1e5){dLd+=tl;++eLd;} if(v.hg_peak(ch)>hgmax)hgmax=v.hg_peak(ch);}
            for(int ch=4;ch<8;++ch) if(v.is_timing(ch)&&v.hg_peak(ch)>=kHG_minPeak){float t5=v.timeOf(ch,RadView::kCFD05),tl=v.timeOf(ch,RadView::kLGCFD);if(t5>-1e5){d5u+=t5;++e5u;}if(tl>-1e5){dLu+=tl;++eLu;} if(v.hg_peak(ch)>hgmax)hgmax=v.hg_peak(ch);}
            if(e5d&&e5u){ E4 q{(float)v.sum_lg(),0.5f*(float)(d5d/e5d-d5u/e5u),(eLd&&eLu)?0.5f*(float)(dLd/eLd-dLu/eLu):-1e9f,hgmax};
                ev.push_back(q); ++nhg; if(hgmax>800) ++nclip; } }
        // sum_lg stats
        double smu=0,srms=0; for(auto&q:ev) smu+=q.slg; smu/=ev.size();
        for(auto&q:ev) srms+=(q.slg-smu)*(q.slg-smu); srms=std::sqrt(srms/ev.size());
        double hgm=0; for(auto&q:ev) hgm+=q.hgmax; hgm/=ev.size();
        printf("  E=%3.0f  Nfid=%6zu  sum_lg=%.0f+-%.0f  <HGmax>=%.0f mV  clip(>800)=%.1f%%\n",
               E,ev.size(),smu,srms,hgm,100.0*nclip/std::max(1L,nhg));
        // quantile bins
        std::sort(ev.begin(),ev.end(),[](const E4&a,const E4&b){return a.slg<b.slg;});
        int NB=9; size_t per=ev.size()/NB; std::vector<double> bx,b5,bL; double best5=1e9,bestL=1e9; int ib5=-1;
        printf("        bin: sumlg     N    cfd05  lgcfd\n");
        for(int b=0;b<NB;++b){ size_t lo=(size_t)b*per,hi=(b==NB-1)?ev.size():lo+per;
            std::vector<float> v5,vL; double sc=0; for(size_t k=lo;k<hi;++k){ v5.push_back(ev[k].d5); if(ev[k].dL>-1e5)vL.push_back(ev[k].dL); sc+=ev[k].slg; }
            double s5=(v5.size()>=300)?tebSigma(v5):-1, sL=(vL.size()>=300)?tebSigma(vL):-1; double xc2=sc/(hi-lo);
            bx.push_back(xc2); b5.push_back(s5); bL.push_back(sL);
            if(s5>15&&s5<best5){best5=s5;ib5=b;} if(sL>15&&sL<bestL)bestL=sL;
            printf("        %2d: %6.0f %6zu  %5.1f  %5.1f\n",b,xc2,hi-lo,s5,sL); }
        printf("    -> best cfd05=%.1f (bin %d)  best lgcfd=%.1f\n",best5,ib5,bestL);
        // panel
        c->cd(e+1); gPad->SetGridy(); gPad->SetLeftMargin(0.16); gPad->SetTopMargin(0.10);
        std::vector<double> x5,y5,xL,yL; for(int b=0;b<NB;++b){ if(b5[b]>0){x5.push_back(bx[b]);y5.push_back(b5[b]);} if(bL[b]>0){xL.push_back(bx[b]);yL.push_back(bL[b]);} }
        TGraph* g5=new TGraph(x5.size(),&x5[0],&y5[0]); g5->SetMarkerStyle(24); g5->SetMarkerColor(kRed+1); g5->SetLineColor(kRed+1); g5->SetMarkerSize(1.2);
        g5->SetTitle(Form("%.0f GeV;#Sigma LG;#sigma_{t} (ps)",E)); g5->GetYaxis()->SetRangeUser(30,70); g5->Draw("ALP");
        TGraph* gL=new TGraph(xL.size(),&xL[0],&yL[0]); gL->SetMarkerStyle(20); gL->SetMarkerColor(kAzure+2); gL->SetLineColor(kAzure+2); gL->SetMarkerSize(1.2); gL->SetLineWidth(2); gL->Draw("LP SAME");
        if(ib5>=0){TMarker* m=new TMarker(bx[ib5],b5[ib5],29); m->SetMarkerColor(kRed+2); m->SetMarkerSize(2.0); m->Draw();}
        fp->Close();
    }
    c->cd(0); TLatex t0; t0.SetNDC(); t0.SetTextFont(62); t0.SetTextSize(0.030);
    t0.DrawLatex(0.20,0.965,Form("%s: #sigma_{t} vs sum_lg quantile bin, per energy  (#circ cfd05  #bullet lgcfd, #star=best cfd05)",build));
    gSystem->mkdir("figures/narrative",kTRUE); c->Print("figures/narrative/luag_diag.png");
    printf("wrote figures/narrative/luag_diag.png\n");
}

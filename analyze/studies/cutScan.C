// ============================================================================
// cutScan.C — the "real experiment" view: time resolution vs AMPLITUDE THRESHOLD.
// For each candidate cut X, sigma_t of ALL events with sum_lg > X (CUMULATIVE), per
// beam energy. Raise the cut -> keep only brighter, better-timed showers -> sigma_t
// falls and flattens toward the floor. The "brightest 1000" is just ONE point on this
// curve (the threshold where N_above = 1000), shown for the top energy. Unlike a fixed
// count, a threshold is dataset-size-INDEPENDENT -- this is how you'd quote it for real.
//   source setup.sh; root -l -b -q 'analyze/studies/cutScan.C+("DSB1")'
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
#include "TGraphErrors.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
using namespace rad;

struct EV{float slg,dmu;};
static int robustSrc(const char* b){ std::string s=b; if(s=="LUAG"||s=="TENERGY") return RadView::kLED; return RadView::kLGCFD; }
static const char* srcName(int s){ return s==RadView::kLED?"led":(s==RadView::kLGCFD?"lgcfd":"cfd05"); }

void cutScan(const char* build="DSB1"){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    BuildConfig cfg=BuildConfig::Load(radConfig(build).Data());
    int SRC=robustSrc(build);
    const double Es[]={25,50,75,100,125,150}; int nE=6;
    const double cutFrac=0.40;
    std::vector<std::vector<EV>> byE(nE);
    for(int e=0;e<nE;++e){ double E=Es[e];
        TFile* fp=TFile::Open(radReduced(build,E)); if(!fp||fp->IsZombie())continue;
        TTree* t=(TTree*)fp->Get("rad"); RadView v; v.attach(t,&cfg);
        double xc,yc; v.beamCenter(xc,yc); double r2=9.0; Long64_t N=v.entries();
        std::vector<EV> tmp;
        for(Long64_t i=0;i<N;++i){ v.get(i);
            if(!v.wc_ok()||v.mcp1_peak()<kMCP1_minPeak||v.mcp1_peak()>kMCP1_maxPeak) continue;
            double dx=v.x_trk()-xc,dy=v.y_trk()-yc; if(dx*dx+dy*dy>=r2) continue;
            double ds=0,us=0;int dn=0,un=0;
            for(int c=0;c<4;++c) if(v.is_timing(c)&&v.hg_peak(c)>=kHG_minPeak){float tc=v.timeOf(c,SRC);if(tc>-1e5){ds+=tc;++dn;}}
            for(int c=4;c<8;++c) if(v.is_timing(c)&&v.hg_peak(c)>=kHG_minPeak){float tc=v.timeOf(c,SRC);if(tc>-1e5){us+=tc;++un;}}
            if(dn<2||un<2) continue;
            tmp.push_back({(float)v.sum_lg(),0.5f*(float)(ds/dn-us/un)});
        }
        fp->Close();
        if(tmp.empty())continue;
        std::vector<float> sl; for(auto&q:tmp)sl.push_back(q.slg);
        std::nth_element(sl.begin(),sl.begin()+sl.size()/2,sl.end()); double cut=cutFrac*sl[sl.size()/2];
        for(auto&q:tmp) if(q.slg>cut) byE[e].push_back(q);
    }

    // floor (brightest-1000 of each energy, slew fit, weighted -- identical to pub_res)
    std::vector<double> fX,fY,fE;
    for(int e=0;e<nE;++e){ auto V=byE[e]; if((int)V.size()<1000)continue;
        std::sort(V.begin(),V.end(),[](const EV&a,const EV&b){return a.slg<b.slg;});
        std::vector<float> d; double sa=0; for(size_t k=V.size()-1000;k<V.size();++k){d.push_back(V[k].dmu);sa+=V[k].slg;}
        double sg=tebSigma(d); if(sg>0){fX.push_back(sa/1000.0);fY.push_back(sg);fE.push_back(sg/std::sqrt(2.0*1000.0));} }
    double bF=0;
    if(fX.size()>=3){ TGraphErrors g(fX.size(),&fX[0],&fY[0],0,&fE[0]);
        double xlo=*std::min_element(fX.begin(),fX.end()),xhi=*std::max_element(fX.begin(),fX.end());
        TF1 ff("ff","sqrt([0]*[0]/(x*x)+[1]*[1])",xlo,xhi); ff.SetParameters(60000,18); g.Fit(&ff,"Q0"); bF=std::fabs(ff.GetParameter(1)); }

    // CUMULATIVE threshold scan for the highest energies
    int showE[]={5,4,3}; int nShow=3;   // 150, 125, 100
    const int ecol[6]={kViolet+1,kAzure+2,kTeal+2,kSpring-6,kOrange+7,kRed+1};
    std::vector<TGraph*> gC(nE,nullptr);
    double ymax=0,xmin=1e9,xmax=-1e9, thr1000=0,sig1000=0;
    printf("\n=== %s (%s): resolution vs amplitude threshold ===\n",build,srcName(SRC));
    for(int si=0;si<nShow;++si){ int e=showE[si]; auto& V=byE[e]; if((int)V.size()<1200)continue;
        std::sort(V.begin(),V.end(),[](const EV&a,const EV&b){return a.slg<b.slg;});
        int Ntot=V.size(); int step=std::max(1,Ntot/140);
        std::vector<double> tx,ty;
        for(int i=0;i+600<Ntot;i+=step){
            std::vector<float> d; d.reserve(Ntot-i); for(int k=i;k<Ntot;++k) d.push_back(V[k].dmu);
            double sg=tebSigma(d); if(sg>0){ tx.push_back(V[i].slg); ty.push_back(sg);
                ymax=std::max(ymax,sg); xmin=std::min(xmin,(double)V[i].slg); xmax=std::max(xmax,(double)V[i].slg); } }
        if(tx.size()<2)continue;
        TGraph* g=new TGraph(tx.size(),&tx[0],&ty[0]); g->SetLineColor(ecol[e]); g->SetLineWidth(si==0?4:2); gC[e]=g;
        if(e==5){ int i0=Ntot-1000; thr1000=V[i0].slg;
            std::vector<float> d; for(int k=i0;k<Ntot;++k)d.push_back(V[k].dmu); sig1000=tebSigma(d);
            printf("  150 GeV: all-contained=%.1f ps, brightest-1000 cut SLG>%.0f -> %.1f ps\n",ty.front(),thr1000,sig1000); }
    }
    if(ymax>110)ymax=110;
    printf("  slew floor b=%.1f ps\n",bF);

    // draw
    TCanvas* c=new TCanvas("ts","",1040,720);
    c->SetLeftMargin(0.10); c->SetRightMargin(0.035); c->SetTopMargin(0.085); c->SetBottomMargin(0.135); c->SetGridy();
    TH1F* fr=c->DrawFrame(xmin*0.9,0,xmax*1.05,ymax*1.14);
    fr->SetTitle(";amplitude threshold  #SigmaLG > X  (a.u.);#sigma_{t} of all events above X  (ps)");
    fr->GetYaxis()->SetTitleSize(0.044); fr->GetYaxis()->SetTitleOffset(1.02); fr->GetYaxis()->SetLabelSize(0.038);
    fr->GetXaxis()->SetTitleSize(0.044); fr->GetXaxis()->SetTitleOffset(1.32); fr->GetXaxis()->SetLabelSize(0.038);
    for(int si=nShow-1;si>=0;--si){ int e=showE[si]; if(gC[e]) gC[e]->Draw("L SAME"); }
    if(bF>0){ TLine* fl=new TLine(xmin*0.9,bF,xmax*1.05,bF); fl->SetLineColor(kBlack); fl->SetLineStyle(2); fl->SetLineWidth(2); fl->Draw();
        TLatex tb; tb.SetTextColor(kBlack); tb.SetTextSize(0.036); tb.DrawLatex(xmin*0.96,bF+ymax*0.025,Form("slew floor  b = %.1f ps  (#SigmaLG #rightarrow #infty)",bF)); }
    if(thr1000>0&&sig1000>0){
        TLine* vl=new TLine(thr1000,0,thr1000,sig1000); vl->SetLineColor(kRed+2); vl->SetLineStyle(3); vl->SetLineWidth(2); vl->Draw();
        TGraph* gm=new TGraph(1,&thr1000,&sig1000); gm->SetMarkerStyle(29); gm->SetMarkerColor(kRed+2); gm->SetMarkerSize(2.7); gm->Draw("P SAME");
        TLatex tm; tm.SetTextColor(kRed+2); tm.SetTextSize(0.033); tm.SetTextAlign(31);
        tm.DrawLatex(thr1000-0.02*(xmax-xmin),sig1000+ymax*0.11,Form("\"brightest 1000\" = cut at #SigmaLG > %.0f  #rightarrow  %.1f ps",thr1000,sig1000)); }
    TLegend* lg=new TLegend(0.74,0.70,0.95,0.90); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.044);
    for(int si=0;si<nShow;++si){int e=showE[si]; if(gC[e]) lg->AddEntry(gC[e],Form("%.0f GeV",Es[e]),"l");} lg->Draw();
    TLatex tt; tt.SetNDC(); tt.SetTextFont(62); tt.SetTextSize(0.038);
    tt.DrawLatex(0.10,0.945,Form("%s (%s): resolution vs amplitude THRESHOLD  #SigmaLG > X",build,srcName(SRC)));
    gSystem->mkdir("figures/narrative",kTRUE); c->Print(Form("figures/narrative/pub_thresh_%s.png",build));
    printf("  wrote figures/narrative/pub_thresh_%s.png\n",build);
}

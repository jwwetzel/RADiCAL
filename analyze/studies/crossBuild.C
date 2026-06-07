// ============================================================================
// crossBuild.C — the full spectrum: (DW-UP)/2 timing for ALL four material
// builds with the de-biased estimator. One tree pass per build/energy collects
// every candidate source; for each source it computes the quantile best-bin
// (OOS-robust) AND the brightest-1% slice (best case, OOS-stable). Picks the
// best source per build, fits sigma_t = a/sqrtE (+) b, and draws all four
// ladders on one canvas + a console table.
//   source setup.sh; root -l -b -q 'analyze/studies/crossBuild.C+'
// ============================================================================
#include "RadView.h"
#include "RadTiming.h"
#include "SelectionCuts.h"
#include "DataPaths.h"
#include "PlotUtils.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TCanvas.h"
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

// quantile best-bin sigma on in-memory (sum_lg, depth) pairs; sets N=winning-bin count
static double quantBest(std::vector<std::pair<float,float>> v, long& N){
    N=0; if(v.size()<2000) return -1; std::sort(v.begin(),v.end());
    int NB=9; size_t per=v.size()/NB; double best=1e9;
    for(int b=0;b<NB;++b){ size_t lo=(size_t)b*per,hi=(b==NB-1)?v.size():lo+per; if(hi-lo<600)continue;
        std::vector<float> t; for(size_t k=lo;k<hi;++k) t.push_back(v[k].second); double s=tebSigma(t);
        if(s>22&&s<best){best=s;N=hi-lo;} } return best<1e8?best:-1;   // >22: no real RADiCAL quantile (DW-UP)/2 bin is lower
}
// brightest-fraction f sigma (deterministic, OOS-stable)
static double brightFrac(std::vector<std::pair<float,float>> v,double f){
    int K=(int)(f*v.size()); if(K<200) return -1;
    std::nth_element(v.begin(),v.begin()+K,v.end(),[](auto&a,auto&b){return a.first>b.first;});
    std::vector<float> t; for(int i=0;i<K;++i) t.push_back(v[i].second); return tebSigma(t);
}

void crossBuild(){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    const char* builds[]={"DSB1","LUAG","MIXED","TENERGY"};
    const char* mat[]={"LYSO:Ce (high light)","LuAG:Ce (low light)","LYSO+LuAG mix","3xLYSO + Energy"};
    const char* matS[]={"LYSO:Ce","LuAG:Ce","LYSO/LuAG mix","3xLYSO+Energy"};
    int col[4]={kAzure+2,kGreen+3,kOrange+8,kViolet+1};
    const int SRC[]={RadView::kCFD05,RadView::kLED,RadView::kLGCFD};
    const char* sname[]={"cfd05","led","lgcfd"}; int nS=3;
    const double Es[]={25,50,75,100,125,150};

    TCanvas* c=new TCanvas("cx","",960,680); c->SetLeftMargin(0.12); c->SetRightMargin(0.04); c->SetGridy();
    TLegend* lg=new TLegend(0.50,0.70,0.95,0.90); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.032);
    TH1F* fr=c->DrawFrame(0,20,160,74); fr->SetTitle("RADiCAL: (DW#minusUP)/2 best-bin #sigma_{t} vs energy, four builds (de-biased);beam energy (GeV);#sigma_{t} (ps)");
    printf("\n=== CROSS-BUILD full spectrum (de-biased) ===\n");
    printf("build    material              best src  sigma@150[quant]  [bright1%%]  floor b\n");
    for(int bi=0;bi<4;++bi){ const char* B=builds[bi];
        BuildConfig cfg=BuildConfig::Load(Form("data/2023/configs/%s.json",B)); if(!cfg.valid())continue;
        // per source: quantile ladder + brightest ladder
        std::vector<std::vector<double>> q(nS),br(nS),Ev(nS),eq(nS);
        for(double E:Es){ TFile* fp=TFile::Open(radReduced(B,E)); if(!fp||fp->IsZombie())continue;
            TTree* t=(TTree*)fp->Get("rad"); RadView v; v.attach(t,&cfg);
            double xc,yc; v.beamCenter(xc,yc); double r2=9.0; Long64_t N=v.entries();
            std::vector<std::vector<std::pair<float,float>>> P(nS);
            for(Long64_t i=0;i<N;++i){ v.get(i);
                if(!v.wc_ok()||v.mcp1_peak()<kMCP1_minPeak||v.mcp1_peak()>kMCP1_maxPeak) continue;
                double dx=v.x_trk()-xc,dy=v.y_trk()-yc; if(dx*dx+dy*dy>=r2) continue;
                for(int s=0;s<nS;++s){ double ds=0,us=0;int dn=0,un=0;
                    for(int ch=0;ch<4;++ch) if(v.is_timing(ch)&&v.hg_peak(ch)>=kHG_minPeak){float tc=v.timeOf(ch,SRC[s]);if(tc>-1e5){ds+=tc;++dn;}}
                    for(int ch=4;ch<8;++ch) if(v.is_timing(ch)&&v.hg_peak(ch)>=kHG_minPeak){float tc=v.timeOf(ch,SRC[s]);if(tc>-1e5){us+=tc;++un;}}
                    if(dn>=1&&un>=1) P[s].push_back({(float)v.sum_lg(),0.5f*(float)(ds/dn-us/un)}); } }
            for(int s=0;s<nS;++s){ long nn; double sq=quantBest(P[s],nn), sb=brightFrac(P[s],0.02);
                if(sq>0){q[s].push_back(sq);br[s].push_back(sb);Ev[s].push_back(E);eq[s].push_back(sq/std::sqrt(2.0*nn));} }
            fp->Close();
        }
        // pick best source by sigma@150 (quantile)
        int best=-1; double bestS=1e9;
        for(int s=0;s<nS;++s){ if(q[s].empty())continue; double s150=q[s].back();
            if(Ev[s].back()==150 && s150>0 && s150<bestS){bestS=s150;best=s;} }
        if(best<0) continue;
        TF1 f("f","sqrt([0]*[0]/x+[1]*[1])",20,160); f.SetParameters(200,25);
        std::vector<double> ze(q[best].size(),0.0);
        TGraphErrors* g=new TGraphErrors(q[best].size(),&Ev[best][0],&q[best][0],&ze[0],&eq[best][0]); g->Fit(&f,"Q");
        // inflate per-point errors by sqrt(chi2/ndf): the bare sigma/sqrt(2N) misses the real
        // point-to-point systematic scatter (run-to-run, selection) -- worst for low-light LuAG.
        double cn=f.GetChisquare()/std::max(1,f.GetNDF()), esc=(cn>1?std::sqrt(cn):1.0);
        for(int p=0;p<g->GetN();++p) g->SetPointError(p,0,g->GetErrorY(p)*esc);
        double floor=std::fabs(f.GetParameter(1)), s150q=q[best].back(), s150b=br[best].back();
        char fb[16]; if(floor>=5) snprintf(fb,16,"%.0f",floor); else snprintf(fb,16,"--");
        char bb[16]; if(s150b>0) snprintf(bb,16,"%.1f",s150b); else snprintf(bb,16,"--");
        printf("%-7s  %-20s  %-7s   %5.1f            %5s       %5s\n",B,mat[bi],sname[best],s150q,bb,fb);
        g->SetMarkerStyle(20+bi); g->SetMarkerColor(col[bi]); g->SetLineColor(col[bi]); g->SetMarkerSize(1.5); g->SetLineWidth(2);
        g->Draw("LP SAME"); TF1* ff=(TF1*)f.Clone(); ff->SetLineColor(col[bi]); ff->SetLineStyle(2); ff->SetLineWidth(1); ff->Draw("SAME");
        lg->AddEntry(g,Form("%s  #color[920]{%s}  (%s)  %.0f ps",B,matS[bi],sname[best],s150q),"lp");
    }
    lg->Draw();
    gSystem->mkdir("figures/narrative",kTRUE); c->Print("figures/narrative/cross_build.png");
    printf("wrote figures/narrative/cross_build.png\n");
}

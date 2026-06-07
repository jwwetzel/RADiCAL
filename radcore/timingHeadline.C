// ============================================================================
// timingHeadline.C — world-class (DW-UP)/2 best-bin ladder for EVERY timing source.
// ----------------------------------------------------------------------------
// One tree pass per energy: LG-weighted beam centroid -> timing fiducial ->
// (DW-UP)/2 for every available per-end time source (cfd03/05/10/20/30/50, led,
// lgcfd, all MCP-referenced) -> sum_lg Gaussian -> 9 energy bins -> best bin ->
// tebSigma (falls back to a bright-slice core when stats are thin). Each source is
// fit sigma_t = a/sqrt(E) (+) b and ranked. cfd05 = the published headline; lgcfd =
// CFD on the LG-predicted true peak.
//   ROOT_INCLUDE_PATH=radcore:Analysis root -l -b -q \
//     'radcore/timingHeadline.C+("DSB1","datasets/2023/reduced/DSB1")'
// ============================================================================
#include "RadView.h"
#include "RadTiming.h"          // rad::tebSigma
#include "SelectionCuts.h"
#include "PlotUtils.h"          // FitGaussCore
#include "DataPaths.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>

using namespace rad;

static double medOf(std::vector<float> v){ if(v.empty())return 0; size_t k=v.size()/2; std::nth_element(v.begin(),v.begin()+k,v.end()); return v[k]; }

// best-bin (>=500) tebSigma; bright-slice core fallback for thin samples
static double bestBinOrBright(std::vector<float>& slg, std::vector<float>& t){
    // t aligned to slg; t<-1e4 = invalid for this source
    std::vector<std::pair<float,float>> v;
    for(size_t i=0;i<slg.size();++i) if(t[i]>-1e4) v.push_back({slg[i],t[i]});
    if(v.size()<800) return -1;
    // 9-bin best-bin
    std::vector<float> sl; for(auto&p:v) sl.push_back(p.first);
    double smin=*std::min_element(sl.begin(),sl.end()), smax=*std::max_element(sl.begin(),sl.end());
    TH1F hS("hS","",150,smin,smax); hS.SetDirectory(0); for(float x:sl) hS.Fill(x);
    double muE,muEe,sigE,sigEe; FitGaussCore(&hS,2.0,muE,muEe,sigE,sigEe); if(sigE<=0){muE=hS.GetMean();sigE=hS.GetRMS();}
    double lo=muE-2*sigE, bw=4*sigE/9.0, best=1e9;
    for(int b=0;b<9;++b){ double blo=lo+b*bw,bhi=blo+bw; std::vector<float> vt;
        for(auto&p:v) if(p.first>=blo&&p.first<bhi) vt.push_back(p.second);
        if(vt.size()<500) continue; double s=rad::tebSigma(vt); if(s>10&&s<best) best=s; }
    if(best<1e8) return best;
    // fallback: bright slice (top 45%) tight core
    std::sort(v.begin(),v.end()); std::vector<float> bt;
    for(size_t i=(size_t)(0.55*v.size());i<v.size();++i) bt.push_back(v[i].second);
    double m=medOf(bt); std::vector<float> c; for(float x:bt) if(std::fabs(x-m)<1.0) c.push_back(x);
    return c.size()>200 ? rad::tebSigma(c) : -1;
}

void timingHeadline(const char* build, const char* dir){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    BuildConfig cfg = BuildConfig::Load(Form("datasets/2023/configs/%s.json",build));
    if(!cfg.valid()){ printf("config load failed: %s\n", cfg.error()); return; }
    const double allE[6]={25,50,75,100,125,150};
    std::vector<double> Es; for(double E:allE) if(!gSystem->AccessPathName(Form("%s/%.0fGeV.root",dir,E))) Es.push_back(E);
    int nE=Es.size(); if(nE<2){ printf("need >=2 energies in %s\n",dir); return; }

    double sig[RadView::kNSrc][6]; bool avail[RadView::kNSrc]={false};
    for(int s=0;s<RadView::kNSrc;++s) for(int e=0;e<6;++e) sig[s][e]=-1;

    for(int e=0;e<nE;++e){ double E=Es[e];
        TFile* fp=TFile::Open(Form("%s/%.0fGeV.root",dir,E)); TTree* t=(TTree*)fp->Get("rad");
        RadView v; v.attach(t,&cfg);
        double xc,yc; v.beamCenter(xc,yc); double rF=TimingFiducialR(E),r2=rF*rF;
        std::vector<float> slg; std::vector<float> tval[RadView::kNSrc];
        Long64_t N=v.entries();
        for(Long64_t i=0;i<N;++i){ v.get(i);
            if(!v.wc_ok()||v.mcp1_peak()<kMCP1_minPeak||v.mcp1_peak()>kMCP1_maxPeak) continue;
            double dx=v.x_trk()-xc, dy=v.y_trk()-yc; if(dx*dx+dy*dy>=r2) continue;
            float lg=v.sum_lg(); slg.push_back(lg);
            for(int s=0;s<RadView::kNSrc;++s){ float tv=kNoTime;
                if(v.hasSrc(s)){ double ds=0,us=0; int dn=0,un=0;
                    for(int c=0;c<4;++c) if(v.is_timing(c)&&v.hg_peak(c)>=kHG_minPeak){float tc=v.timeOf(c,s); if(tc>-1e4){ds+=tc;++dn;}}
                    for(int c=4;c<8;++c) if(v.is_timing(c)&&v.hg_peak(c)>=kHG_minPeak){float tc=v.timeOf(c,s); if(tc>-1e4){us+=tc;++un;}}
                    if(dn>=1&&un>=1) tv=0.5f*(float)(ds/dn-us/un); }
                tval[s].push_back(tv); }
        }
        for(int s=0;s<RadView::kNSrc;++s){ if(!v.hasSrc(s)) continue; avail[s]=true;
            sig[s][e]=bestBinOrBright(slg,tval[s]); }
        printf("E=%3.0f  Nfid=%zu", E, slg.size());
        for(int s=0;s<RadView::kNSrc;++s) if(avail[s]) printf("  %s=%.1f", RadView::srcName(s), sig[s][e]);
        printf(" ps\n");
        fp->Close();
    }

    // stochastic fit + ranking
    printf("\n=== %s : sigma_t = a/sqrt(E) (+) b  (ranked by b) ===\n", build);
    struct Row{int s; double a,b,s150;};
    std::vector<Row> rows; double Ea[6]; for(int e=0;e<nE;++e) Ea[e]=Es[e];
    for(int s=0;s<RadView::kNSrc;++s){ if(!avail[s]) continue;
        double y[6]; int n=0; for(int e=0;e<nE;++e) if(sig[s][e]>0){y[n]=sig[s][e];Ea[n]=Es[e];++n;}
        if(n<2) continue;
        TGraph g(n,Ea,y); TF1 f("f","sqrt([0]*[0]/x+[1]*[1])",20,160); f.SetParameters(200,25); g.Fit(&f,"Q");
        double s150=-1; for(int e=0;e<nE;++e) if(std::fabs(Es[e]-150)<1) s150=sig[s][e];
        rows.push_back({s,f.GetParameter(0),std::fabs(f.GetParameter(1)),s150}); }
    std::sort(rows.begin(),rows.end(),[](const Row&A,const Row&B){return A.b<B.b;});
    for(auto&r:rows) printf("  %-6s : a=%4.0f ps*sqrt(GeV)  b=%5.1f ps   (sigma@150=%.1f ps)\n",
                            RadView::srcName(r.s), r.a, r.b, r.s150);

    // plot: published cfd05 vs best (lgcfd if present)
    int sB=RadView::kCFD05, sBest=rows.empty()?RadView::kCFD05:rows[0].s;
    TCanvas* c=new TCanvas("c_th","",900,680); c->SetGridx(); c->SetGridy();
    TH1F* fr=c->DrawFrame(0,0,160,std::max(60.0,sig[sB][0]*1.2));
    fr->SetTitle(Form("%s headline (DW#minusUP)/2 best-bin: every timing source;beam energy (GeV);#sigma_{t} (ps)",build));
    int cols[3]={(int)(kGray+2),(int)(kRed+1),(int)(kAzure+1)};
    auto draw=[&](int s,int col,int mk){ double y[6];int n=0; for(int e=0;e<nE;++e){ if(sig[s][e]>0){y[n]=sig[s][e];Ea[n]=Es[e];++n;} }
        TGraph* g=new TGraph(n,Ea,y); g->SetMarkerStyle(mk); g->SetMarkerColor(col); g->SetLineColor(col); g->SetMarkerSize(1.6); g->SetLineWidth(2); g->Draw("PL SAME"); return g; };
    TGraph* g0=draw(sB,cols[0],24);
    TGraph* g1=(sBest!=sB)?draw(sBest,cols[1],20):nullptr;
    TLegend* lg=new TLegend(0.40,0.74,0.88,0.88); lg->SetBorderSize(0);
    lg->AddEntry(g0,"cfd05 (published headline)","pl");
    if(g1) lg->AddEntry(g1,Form("%s (best)",RadView::srcName(sBest)),"pl"); lg->Draw();
    gSystem->mkdir("radcore/figs",kTRUE); c->Print(Form("radcore/figs/timing_headline_%s.png",build));
    printf("wrote radcore/figs/timing_headline_%s.png\n", build);
}

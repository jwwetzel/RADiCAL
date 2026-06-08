// ============================================================================
// optAppendix.C — Paper 1 selection/optimization appendix, all four builds in one
// figure. The four choices that define the headline brightest-1000 (DW-UP)/2
// estimator, each justified by a scan at the headline energy (150 GeV):
//
//   (a) timing source : sigma_t vs CFD fraction on the CLIPPED peak {3,5,10,30,50}%
//                       (solid) and the recovered-edge source we ADOPT (lgcfd/led,
//                       dashed horizontal) -- which beats the whole clipped family.
//   (b) brightest-N   : sigma_t vs K (200..10000) -- the stat/selection trade-off,
//                       why K=1000.
//   (c) fiducial r    : sigma_t vs r_fid (1.5..4 mm) -- the plateau, why r=3.0.
//   (d) HG threshold  : sigma_t vs per-channel HG min (10..50 mV) -- robustness to
//                       the channel noise floor, why 20 mV.
//
// One event pass per build at 150 GeV caches every per-event quantity (position,
// MCP, per-channel HG amplitude, all seven candidate source times); every scan
// point is then computed in-memory. Reduced files = complete statistics; identical
// on Argon. Complements optScan.C (the per-build per-energy best-bin landscape).
//
//   source setup.sh
//   root -l -b -q 'analyze/studies/optAppendix.C+'
// Output: figures/narrative/optimization.png
// ============================================================================
#include "RadView.h"
#include "RadTiming.h"        // rad::tebSigma
#include "DataPaths.h"
#include "BuildConfig.h"
#include "SelectionCuts.h"
#include "PlotUtils.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TGraph.h"
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

static const int    SRC[7] = { RadView::kCFD03, RadView::kCFD05, RadView::kCFD10,
                               RadView::kCFD30, RadView::kCFD50, RadView::kLED, RadView::kLGCFD };
static const double CFDFRAC[5] = {3,5,10,30,50};
static const int    IDX_LED=5, IDX_LGCFD=6;

static int  robustIdx(const char* b){ std::string s=b; return (s=="LUAG"||s=="TENERGY")?IDX_LED:IDX_LGCFD; }

struct Ev { float slg, x, y, mcp; float hg[8]; float t[8][7]; };

static void loadFile(const char* path, BuildConfig& cfg, std::vector<Ev>& out, double& xc, double& yc){
    out.clear(); xc=yc=0;
    TFile* fp=TFile::Open(path); if(!fp||fp->IsZombie()){ if(fp)fp->Close(); return; }
    TTree* t=(TTree*)fp->Get("rad"); if(!t){ fp->Close(); return; }
    RadView v; v.attach(t,&cfg);
    Long64_t N=v.entries(); out.reserve(N);
    double wx=0,wy=0,w=0;
    for(Long64_t i=0;i<N;++i){ v.get(i);
        if(!v.wc_ok()) continue;
        double slg=v.sum_lg();
        if(v.mcp1_peak()>=50.0f && slg>300.0){ wx+=v.x_trk()*slg; wy+=v.y_trk()*slg; w+=slg; }
        Ev e; e.slg=(float)slg; e.x=v.x_trk(); e.y=v.y_trk(); e.mcp=v.mcp1_peak();
        for(int c=0;c<8;++c){ bool tim=v.is_timing(c);
            e.hg[c]= tim? v.hg_peak(c) : -1.f;
            for(int s=0;s<7;++s) e.t[c][s]= tim? v.timeOf(c,SRC[s]) : kNoTime; }
        out.push_back(e);
    }
    xc=(w>0)?wx/w:0; yc=(w>0)?wy/w:0;
    fp->Close();
}

static double sig(const std::vector<Ev>& ev, double xc, double yc,
                  double rFid, int K, int sidx, float mcpHi, float hgMin){
    double r2=rFid*rFid;
    std::vector<std::pair<float,float>> sd; sd.reserve(ev.size()/4+16);
    for(const Ev& e : ev){
        if(e.mcp<kMCP1_minPeak || e.mcp>mcpHi) continue;
        double dx=e.x-xc, dy=e.y-yc; if(dx*dx+dy*dy>=r2) continue;
        double ds=0,us=0; int dn=0,un=0;
        for(int c=0;c<4;++c) if(e.hg[c]>=hgMin){ float tc=e.t[c][sidx]; if(tc>-1e5f){ds+=tc;++dn;} }
        for(int c=4;c<8;++c) if(e.hg[c]>=hgMin){ float tc=e.t[c][sidx]; if(tc>-1e5f){us+=tc;++un;} }
        if(dn<1||un<1) continue;
        sd.push_back({ e.slg, 0.5f*(float)(ds/dn-us/un) });
    }
    if((int)sd.size()<K) return -1;
    std::nth_element(sd.begin(), sd.begin()+K, sd.end(),
                     [](const std::pair<float,float>&a,const std::pair<float,float>&b){ return a.first>b.first; });
    std::vector<float> vt; vt.reserve(K); for(int i=0;i<K;++i) vt.push_back(sd[i].second);
    return tebSigma(vt);
}

void optAppendix(){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    const char* builds[]={"DSB1","TENERGY","MIXED","LUAG"};
    int col[4]={kAzure+2,kViolet+1,kOrange+8,kGreen+3};
    int mk [4]={20,21,22,23};
    const double E=150;

    std::vector<Ev> cache[4]; double XC[4],YC[4]; double rNom=TimingFiducialR(E);
    for(int bi=0;bi<4;++bi){ BuildConfig cfg=BuildConfig::Load(Form("data/2023/configs/%s.json",builds[bi]));
        if(!cfg.valid()){printf("%s: no config\n",builds[bi]);continue;}
        TString p=radReduced(builds[bi],E); if(gSystem->AccessPathName(p.Data())){printf("%s: no 150 file\n",builds[bi]);continue;}
        loadFile(p.Data(),cfg,cache[bi],XC[bi],YC[bi]);
        printf("  %s: cached %zu events @150 GeV\n",builds[bi],cache[bi].size());
    }

    const int nK=6;  double Kg[nK]={200,500,1000,2000,5000,10000};
    const int nR=6;  double Rg[nR]={1.5,2.0,2.5,3.0,3.5,4.0};
    const int nH=5;  double Hg[nH]={10,20,30,40,50};
    const double YLO=15, YHI=78;

    TCanvas* c=new TCanvas("opt","",1180,920); c->Divide(2,2,0.006,0.006);
    TLegend* leg=new TLegend(0.17,0.66,0.52,0.90); leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.042);

    // (a) CFD-fraction family + adopted recovered-edge reference
    c->cd(1); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14); gPad->SetLogx();
    { TH1F* fr=gPad->DrawFrame(2,YLO,70,YHI); fr->SetTitle("(a) timing source: CFD fraction on the clipped peak;CFD fraction (%);#sigma_{t}(150) (ps)");
      fr->GetXaxis()->SetTitleSize(0.05);fr->GetYaxis()->SetTitleSize(0.05);fr->GetYaxis()->SetTitleOffset(1.2);
      fr->GetXaxis()->SetMoreLogLabels(); fr->GetXaxis()->SetNoExponent();
      for(int bi=0;bi<4;++bi){ if(cache[bi].empty())continue;
        TGraph* g=new TGraph(); int np=0;
        for(int s=0;s<5;++s){ double y=sig(cache[bi],XC[bi],YC[bi],rNom,1000,s,750.f,20.f); if(y>0)g->SetPoint(np++,CFDFRAC[s],y); }
        g->SetLineColor(col[bi]);g->SetMarkerColor(col[bi]);g->SetMarkerStyle(mk[bi]);g->SetMarkerSize(1.3);g->SetLineWidth(2);g->Draw("PL SAME");
        int ridx=robustIdx(builds[bi]); double yr=sig(cache[bi],XC[bi],YC[bi],rNom,1000,ridx,750.f,20.f);
        if(yr>0){ TLine* l=new TLine(2,yr,70,yr); l->SetLineColor(col[bi]); l->SetLineStyle(2); l->SetLineWidth(2); l->Draw(); }
        leg->AddEntry(g,builds[bi],"pl");
      }
      TLatex t;t.SetNDC();t.SetTextSize(0.033);t.SetTextColor(kGray+3);
      t.DrawLatex(0.40,0.32,"solid: clipped-peak CFD family");
      t.DrawLatex(0.40,0.265,"dashed: adopted recovered edge");
      t.DrawLatex(0.40,0.21,"(lgcfd / led) #minus below the family");
      leg->Draw();
    }
    // (b) brightest-N ladder
    c->cd(2); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14); gPad->SetLogx();
    { TH1F* fr=gPad->DrawFrame(150,YLO,13000,YHI); fr->SetTitle("(b) selection tightness: brightest-N;K (brightest showers by #Sigma LG);#sigma_{t}(150) (ps)");
      fr->GetXaxis()->SetTitleSize(0.05);fr->GetYaxis()->SetTitleSize(0.05);fr->GetYaxis()->SetTitleOffset(1.2);
      for(int bi=0;bi<4;++bi){ if(cache[bi].empty())continue; int ridx=robustIdx(builds[bi]);
        TGraph* g=new TGraph(); int np=0;
        for(int k=0;k<nK;++k){ double y=sig(cache[bi],XC[bi],YC[bi],rNom,(int)Kg[k],ridx,750.f,20.f); if(y>0)g->SetPoint(np++,Kg[k],y); }
        g->SetLineColor(col[bi]);g->SetMarkerColor(col[bi]);g->SetMarkerStyle(mk[bi]);g->SetMarkerSize(1.3);g->SetLineWidth(2);g->Draw("PL SAME"); }
      TLine* l=new TLine(1000,YLO,1000,YHI); l->SetLineColor(kGray+1); l->SetLineStyle(2); l->Draw();
      TLatex t;t.SetNDC();t.SetTextSize(0.036);t.SetTextColor(kGray+3); t.DrawLatex(0.42,0.84,"adopted K=1000");
    }
    // (c) fiducial radius
    c->cd(3); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14);
    { TH1F* fr=gPad->DrawFrame(1.2,YLO,4.3,YHI); fr->SetTitle("(c) fiducial radius;timing fiducial r (mm);#sigma_{t}(150) (ps)");
      fr->GetXaxis()->SetTitleSize(0.05);fr->GetYaxis()->SetTitleSize(0.05);fr->GetYaxis()->SetTitleOffset(1.2);
      for(int bi=0;bi<4;++bi){ if(cache[bi].empty())continue; int ridx=robustIdx(builds[bi]);
        TGraph* g=new TGraph(); int np=0;
        for(int r=0;r<nR;++r){ double y=sig(cache[bi],XC[bi],YC[bi],Rg[r],1000,ridx,750.f,20.f); if(y>0)g->SetPoint(np++,Rg[r],y); }
        g->SetLineColor(col[bi]);g->SetMarkerColor(col[bi]);g->SetMarkerStyle(mk[bi]);g->SetMarkerSize(1.3);g->SetLineWidth(2);g->Draw("PL SAME"); }
      TLine* l=new TLine(3.0,YLO,3.0,YHI); l->SetLineColor(kGray+1); l->SetLineStyle(2); l->Draw();
      TLatex t;t.SetNDC();t.SetTextSize(0.036);t.SetTextColor(kGray+3); t.DrawLatex(0.54,0.84,"adopted r=3.0 mm");
    }
    // (d) HG threshold
    c->cd(4); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14);
    { TH1F* fr=gPad->DrawFrame(5,YLO,55,YHI); fr->SetTitle("(d) per-channel HG threshold;HG min amplitude (mV);#sigma_{t}(150) (ps)");
      fr->GetXaxis()->SetTitleSize(0.05);fr->GetYaxis()->SetTitleSize(0.05);fr->GetYaxis()->SetTitleOffset(1.2);
      for(int bi=0;bi<4;++bi){ if(cache[bi].empty())continue; int ridx=robustIdx(builds[bi]);
        TGraph* g=new TGraph(); int np=0;
        for(int h=0;h<nH;++h){ double y=sig(cache[bi],XC[bi],YC[bi],rNom,1000,ridx,750.f,(float)Hg[h]); if(y>0)g->SetPoint(np++,Hg[h],y); }
        g->SetLineColor(col[bi]);g->SetMarkerColor(col[bi]);g->SetMarkerStyle(mk[bi]);g->SetMarkerSize(1.3);g->SetLineWidth(2);g->Draw("PL SAME"); }
      TLine* l=new TLine(20,YLO,20,YHI); l->SetLineColor(kGray+1); l->SetLineStyle(2); l->Draw();
      TLatex t;t.SetNDC();t.SetTextSize(0.036);t.SetTextColor(kGray+3); t.DrawLatex(0.52,0.84,"adopted 20 mV");
    }

    c->cd(0); TLatex tt; tt.SetNDC(); tt.SetTextFont(62); tt.SetTextSize(0.026);
    tt.DrawLatex(0.06,0.965,"Selection optimization at 150 GeV: the four choices behind the headline estimator");
    gSystem->mkdir("figures/narrative",kTRUE);
    c->Print("figures/narrative/optimization.png");
    printf("  wrote figures/narrative/optimization.png\n");
}

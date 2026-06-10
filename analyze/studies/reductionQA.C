// reductionQA.C — 16-panel reduction QA per build: what the reducer extracted from
// the raw waveforms, energies overlaid (house palette, shared legend). Panels 1-8:
// per-channel HG peak distributions, LOG-y (pedestal + ~820 mV clip pile-up both
// visible), red line at 820. 9: sum_lg shower spectrum + energy legend. 10: MCP1
// reference (+200-750 mV window). 11-12: beam x/y. 13: beam x:y 2D. 14: mean HG peak
// per channel (markers). 15: cfd05 (DW-UP)/2. 16: fiducial events per energy.
//   source setup.sh; root -l -b -q 'analyze/studies/reductionQA.C+("DSB1")'
// Output: figures/<year>/narrative/reduction_<BUILD>.png
#include "RadView.h"
#include "DataPaths.h"
#include "FigPaths.h"
#include "PlotUtils.h"       // ApplyRADiCALStyle, DrawPadTitle, DrawSuperTitle, DrawEnergyLegend, kREnergyCols
#include "SelectionCuts.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include <vector>
#include <cstdio>
using namespace rad;

static const char* CN[8]={"NW-D","NE-D","SE-D","SW-D","NW-U","NE-U","SE-U","SW-U"};

void reductionQA(const char* build="DSB1"){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    BuildConfig cfg=BuildConfig::Load(radConfig(build).Data());
    const double Es[6]={25,50,75,100,125,150};
    TH1F* hpk[6][8]={}; TH1F* hslg[6]={}; TH1F* hmcp[6]={}; TH1F* hx[6]={}; TH1F* hy[6]={}; TH1F* hdw[6]={};
    TH2F* hbeam=new TH2F(Form("beam_%s",build),";x_{trk} (mm);y_{trk} (mm)",80,-15,15,80,-15,15); hbeam->SetDirectory(nullptr);
    double meanpk[6][8]={}; long npk[6][8]={}; long nfid[6]={}; bool valid[6]={};
    std::vector<double> Epresent; std::vector<int> Ecols;
    int topE=-1;
    for(int e=0;e<6;++e){ int E=(int)Es[e]; TString p=radReduced(build,E); if(gSystem->AccessPathName(p.Data()))continue;
        TFile* f=TFile::Open(p.Data()); if(!f||f->IsZombie()){if(f)f->Close();continue;} TTree* t=(TTree*)f->Get("rad"); if(!t){f->Close();continue;}
        RadView v; v.attach(t,&cfg); double xc,yc; v.beamCenter(xc,yc); double r2=TimingFiducialR(E)*TimingFiducialR(E);
        for(int c=0;c<8;++c){ hpk[e][c]=new TH1F(Form("pk%d_%d_%s",e,c,build),"",95,0,950); hpk[e][c]->SetDirectory(nullptr); }
        hslg[e]=new TH1F(Form("slg%d_%s",e,build),"",140,0,7000); hslg[e]->SetDirectory(nullptr);
        hmcp[e]=new TH1F(Form("mcp%d_%s",e,build),"",100,0,800); hmcp[e]->SetDirectory(nullptr);
        hx[e]=new TH1F(Form("x%d_%s",e,build),"",80,-15,15); hx[e]->SetDirectory(nullptr);
        hy[e]=new TH1F(Form("y%d_%s",e,build),"",80,-15,15); hy[e]->SetDirectory(nullptr);
        hdw[e]=new TH1F(Form("dw%d_%s",e,build),"",120,-0.6,0.4); hdw[e]->SetDirectory(nullptr);
        Long64_t N=v.entries();
        for(Long64_t i=0;i<N;++i){ v.get(i);
            hmcp[e]->Fill(v.mcp1_peak());
            for(int c=0;c<8;++c){ float pk=v.hg_peak(c); hpk[e][c]->Fill(pk); if(pk>20){meanpk[e][c]+=pk;++npk[e][c];} }
            if(!v.wc_ok()||v.mcp1_peak()<kMCP1_minPeak||v.mcp1_peak()>kMCP1_maxPeak) continue;
            double dx=v.x_trk()-xc,dy=v.y_trk()-yc; if(dx*dx+dy*dy>=r2) continue;
            ++nfid[e]; hslg[e]->Fill(v.sum_lg()); hx[e]->Fill(v.x_trk()); hy[e]->Fill(v.y_trk());
            if(topE<0||e>topE) hbeam->Fill(v.x_trk(),v.y_trk());
            double ds=0,us=0;int dn=0,un=0;
            for(int c=0;c<4;++c) if(v.is_timing(c)&&v.hg_peak(c)>=kHG_minPeak){float tc=v.timeOf(c,RadView::kCFD05);if(tc>-1e5f){ds+=tc;++dn;}}
            for(int c=4;c<8;++c) if(v.is_timing(c)&&v.hg_peak(c)>=kHG_minPeak){float tc=v.timeOf(c,RadView::kCFD05);if(tc>-1e5f){us+=tc;++un;}}
            if(dn>=1&&un>=1) hdw[e]->Fill(0.5*(ds/dn-us/un));
        }
        for(int c=0;c<8;++c) if(npk[e][c]) meanpk[e][c]/=npk[e][c];
        valid[e]=true; topE=e; Epresent.push_back(Es[e]); Ecols.push_back(kREnergyCols[e]); f->Close();
    }
    // highest filled sum_lg bin across energies -> tight common x-range
    double slgHi=1000; for(int e=0;e<6;++e) if(hslg[e]) for(int b=hslg[e]->GetNbinsX();b>=1;--b) if(hslg[e]->GetBinContent(b)>3){ slgHi=std::max(slgHi,hslg[e]->GetXaxis()->GetBinUpEdge(b)); break; }

    auto overlay=[&](TVirtualPad* pad,TH1F* h[6],const char* title,double clipLine,bool logy,double xhi){
        pad->cd(); gPad->SetLeftMargin(0.135);gPad->SetBottomMargin(0.135);gPad->SetRightMargin(0.035);gPad->SetTopMargin(0.10); if(logy)gPad->SetLogy();
        double ymax=0; for(int e=0;e<6;++e) if(h[e])ymax=std::max(ymax,h[e]->GetMaximum()); if(ymax<=0)ymax=1; bool first=true;
        for(int e=0;e<6;++e){ if(!h[e]||h[e]->GetEntries()<5)continue; h[e]->SetLineColor(kREnergyCols[e]);h[e]->SetLineWidth(2);
            h[e]->SetMaximum(logy?ymax*3:1.30*ymax); if(logy)h[e]->SetMinimum(0.5);
            if(xhi>0)h[e]->GetXaxis()->SetRangeUser(0,xhi);
            h[e]->SetTitle(""); h[e]->GetXaxis()->SetLabelSize(0.055);h[e]->GetYaxis()->SetLabelSize(0.05);
            h[e]->GetYaxis()->SetNoExponent(); h[e]->GetXaxis()->SetNdivisions(508);   // wide 0.135 margin holds 4-digit counts; no clipped x10^n
            h[e]->Draw(first?"HIST":"HIST SAME"); first=false; }
        if(clipLine>0){ TLine* l=new TLine(clipLine,logy?0.5:0,clipLine,logy?ymax*3:1.30*ymax); l->SetLineStyle(2);l->SetLineColor(kRed);l->SetLineWidth(2);l->Draw(); }
        DrawPadTitle(title);
    };
    TCanvas* c=new TCanvas(Form("rq_%s",build),"",1860,1860);
    TPad* g=GridWithTitle(c,4,4,Form("%s reduction QA: extracted features (energies in legend); red dash = 820 mV HG clip, gray = 200 mV MCP gate",build),0.006,0.014,0.045,0.019);
    // 1-8: per-channel HG peaks, LOG-y, clip at 820, x to 950 (so '900' is inside the frame)
    for(int ci=0;ci<8;++ci){ TH1F* col[6]; for(int e=0;e<6;++e) col[e]=hpk[e][ci];
        overlay(g->cd(ci+1),col,Form("HG peak %s (mV)",CN[ci]),820,true,0); }
    overlay(g->cd(9), hslg,"sum_lg (shower energy, mV)",0,false,slgHi*1.08);
    // energy legend lives in the (empty upper-right of the) sum_lg panel
    g->cd(9); DrawEnergyLegend(0.62,0.50,0.97,0.90,Epresent,Ecols,0.060,"beam E");
    overlay(g->cd(10),hmcp,"MCP1 peak (mV)",0,false,0); { g->cd(10); TLine*a=new TLine(200,0,200,1e7);a->SetLineStyle(2);a->SetLineColor(kGray+2);a->Draw(); }
    overlay(g->cd(11),hx,"beam x_{trk} (mm)",0,false,0);
    overlay(g->cd(12),hy,"beam y_{trk} (mm)",0,false,0);
    // 13: beam 2D
    g->cd(13); gPad->SetLeftMargin(0.155);gPad->SetBottomMargin(0.135);gPad->SetRightMargin(0.16);gPad->SetTopMargin(0.10);
    hbeam->GetXaxis()->SetTitle("x_{trk} (mm)"); hbeam->GetYaxis()->SetTitle("y_{trk} (mm)");
    hbeam->GetXaxis()->SetTitleSize(0.05);hbeam->GetYaxis()->SetTitleSize(0.05);hbeam->GetYaxis()->SetTitleOffset(1.35);
    hbeam->GetXaxis()->SetLabelSize(0.05);hbeam->GetYaxis()->SetLabelSize(0.05);hbeam->GetZaxis()->SetLabelSize(0.04);
    hbeam->Draw("COLZ"); DrawPadTitle("beam profile x:y (top E)");
    // 14: mean HG peak per channel — MARKERS ONLY (no connecting sawtooth)
    g->cd(14); gPad->SetLeftMargin(0.135);gPad->SetBottomMargin(0.135);gPad->SetRightMargin(0.035);gPad->SetTopMargin(0.10);
    TH1F* fr=(TH1F*)gPad->DrawFrame(-0.5,0,7.5,950); fr->SetTitle(";channel index;mean HG peak (mV)");
    fr->GetXaxis()->SetTitleSize(0.05);fr->GetYaxis()->SetTitleSize(0.05);fr->GetYaxis()->SetMaxDigits(3);fr->GetXaxis()->SetNdivisions(509);
    for(int e=0;e<6;++e){ if(!valid[e])continue; TGraph* gr=new TGraph(); for(int ci=0;ci<8;++ci) gr->SetPoint(ci,ci,meanpk[e][ci]);
        gr->SetMarkerStyle(20);gr->SetMarkerColor(kREnergyCols[e]);gr->SetMarkerSize(1.3);gr->Draw("P SAME"); }
    { TLine* l=new TLine(-0.5,820,7.5,820);l->SetLineStyle(2);l->SetLineColor(kRed);l->SetLineWidth(2);l->Draw(); } DrawPadTitle("mean HG peak per channel");
    overlay(g->cd(15),hdw,"(DW#minusUP)/2 cfd05 (ns)",0,false,0);
    // 16: fiducial events per energy
    g->cd(16); gPad->SetLeftMargin(0.155);gPad->SetBottomMargin(0.135);gPad->SetRightMargin(0.035);gPad->SetTopMargin(0.10);
    TH1F* hbar=new TH1F(Form("nfid_%s",build),";beam energy (GeV);fiducial events (10^{3})",6,0,6); hbar->SetDirectory(nullptr);
    for(int e=0;e<6;++e){ hbar->GetXaxis()->SetBinLabel(e+1,Form("%d",(int)Es[e])); hbar->SetBinContent(e+1,nfid[e]/1000.0); }  // 10^3 units -> no clipped exponent
    hbar->SetMinimum(0); hbar->GetYaxis()->SetNoExponent(); hbar->GetXaxis()->SetLabelSize(0.06); hbar->GetYaxis()->SetTitleOffset(1.4);
    hbar->GetXaxis()->SetTitleSize(0.05);hbar->GetYaxis()->SetTitleSize(0.05);
    hbar->SetFillColorAlpha(kAzure+1,0.6);hbar->SetLineColor(kAzure+2);hbar->SetBarWidth(0.85);hbar->Draw("HIST"); DrawPadTitle("fiducial events per energy");
    c->Print(radFigP(Form("figures/narrative/reduction_%s.png",build)));
    printf("wrote figures/narrative/reduction_%s.png  (fid/E:",build); for(int e=0;e<6;++e) if(valid[e]) printf(" %d:%ld",(int)Es[e],nfid[e]); printf(")\n");
}

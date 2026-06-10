// ============================================================================
// mixedSeparate.C — measure the MIXED module's two materials SEPARATELY, by
// forming (DW-UP)/2 only over the same-material diagonal corners.
//   MIXED corner map (verified vs clip fraction): NE,SW = DSB1 (bright, clip
//   33-50%); NW,SE = LuAG (dim, clip 0-8%).
//   kCap idx: 0NW-D 1NE-D 2SE-D 3SW-D 4NW-U 5NE-U 6SE-U 7SW-U.
//   DSB1 diagonal: down{NE-D=1, SW-D=3}, up{NE-U=5, SW-U=7}  -> hg_lgcfd (they clip)
//   LuAG diagonal: down{NW-D=0, SE-D=2}, up{NW-U=4, SE-U=6}  -> hg_led   (no clip)
// Each material's (DW-UP)/2 uses ONE consistent clipping regime + ONE method, so
// there is no mixed-method offset and no heterogeneous-corner averaging. If the two
// come out clean/monotonic with DSB1 < LuAG, the corrupted full-MIXED sigma(E) was
// the heterogeneous average -- and light yield, not the module, sets the timing.
// Reference: the pure DSB1 and LuAG BUILDS (4 corners) for context.
// All crossings are already in the reduced tree -> no re-reduction.
//   source setup.sh; root -l -b -q 'analyze/studies/mixedSeparate.C+'
// Output: figures/<year>/narrative/mixed_separate.png
// ============================================================================
#include "RadView.h"
#include "RadTiming.h"       // rad::tebSigma
#include "DataPaths.h"
#include "FigPaths.h"
#include "PlotUtils.h"
#include "SelectionCuts.h"
#include "TFile.h"
#include "TTree.h"
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
using namespace rad;

// (DW-UP)/2 over a subset of corners (down[],up[]), src, brightest-K by TOTAL sum_lg
static double sigmaSubset(RadView& v, double E, const int* dn, const int* up, int nc, int src, int K=1000){
    double xc,yc; v.beamCenter(xc,yc); double r2=TimingFiducialR(E)*TimingFiducialR(E);
    std::vector<std::pair<float,float>> sd; Long64_t N=v.entries();
    for(Long64_t i=0;i<N;++i){ v.get(i);
        if(!v.wc_ok()||v.mcp1_peak()<kMCP1_minPeak||v.mcp1_peak()>kMCP1_maxPeak) continue;
        double dx=v.x_trk()-xc,dy=v.y_trk()-yc; if(dx*dx+dy*dy>=r2) continue;
        double ds=0,us=0;int nd=0,nu=0;
        for(int k=0;k<nc;++k){ int c=dn[k]; if(v.is_timing(c)&&v.hg_peak(c)>=kHG_minPeak){float tc=v.timeOf(c,src); if(tc>-1e5f){ds+=tc;++nd;}} }
        for(int k=0;k<nc;++k){ int c=up[k]; if(v.is_timing(c)&&v.hg_peak(c)>=kHG_minPeak){float tc=v.timeOf(c,src); if(tc>-1e5f){us+=tc;++nu;}} }
        if(nd<1||nu<1) continue;
        sd.push_back({(float)v.sum_lg(), 0.5f*(float)(ds/nd-us/nu)});
    }
    if((int)sd.size()<K) return -1;
    std::nth_element(sd.begin(),sd.begin()+K,sd.end(),[](const std::pair<float,float>&a,const std::pair<float,float>&b){return a.first>b.first;});
    std::vector<float> vt; for(int i=0;i<K;++i) vt.push_back(sd[i].second);
    return tebSigma(vt);
}

struct Ld { std::vector<double> E,S,Se,ze; double a=0,b=0,be=0,chi=0; };
static void fitAB(Ld& L){ if(L.E.size()<3) return;
    TGraphErrors g(L.E.size(),&L.E[0],&L.S[0],&L.ze[0],&L.Se[0]);
    TF1 f("f","sqrt([0]*[0]/x+[1]*[1])",20,160); f.SetParameters(300,18); g.Fit(&f,"QN");
    L.a=std::fabs(f.GetParameter(0)); L.b=std::fabs(f.GetParameter(1));
    L.chi=f.GetChisquare()/std::max(1,f.GetNDF()); L.be=f.GetParError(1); if(L.chi>1) L.be*=std::sqrt(L.chi); }
static void add(Ld& L,double e,double s){ if(s>0){ L.E.push_back(e); L.S.push_back(s); L.Se.push_back(s/std::sqrt(2000.0)); L.ze.push_back(0); } }

void mixedSeparate(){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    const double Es[6]={25,50,75,100,125,150};
    // MIXED diagonals
    int dsbD[2]={1,3}, dsbU[2]={5,7};   // DSB1 corners NE,SW
    int luaD[2]={0,2}, luaU[2]={4,6};   // LuAG corners NW,SE
    int allD[4]={0,1,2,3}, allU[4]={4,5,6,7};

    Ld mD, mL, bD, bL, bD2, bL2;   // MIXED-DSB1, MIXED-LuAG, builds 4-corner, builds 2-corner-diagonal CONTROL
    // --- MIXED, split by material ---
    { BuildConfig cfg=BuildConfig::Load(radConfig("MIXED").Data());
      for(double e:Es){ TString p=radReduced("MIXED",e); if(gSystem->AccessPathName(p.Data()))continue;
        TFile* f=TFile::Open(p.Data()); TTree* t=(TTree*)f->Get("rad"); if(!t){f->Close();continue;}
        RadView v; v.attach(t,&cfg);
        add(mD,e,sigmaSubset(v,e,dsbD,dsbU,2,RadView::kLGCFD));   // DSB1 corners -> lgcfd
        add(mL,e,sigmaSubset(v,e,luaD,luaU,2,RadView::kLED));     // LuAG corners -> led
        f->Close(); } }
    // --- reference builds: full 4 corners AND the SAME 2-corner diagonal (geometry control) ---
    { BuildConfig cfg=BuildConfig::Load(radConfig("DSB1").Data());
      for(double e:Es){ TString p=radReduced("DSB1",e); if(gSystem->AccessPathName(p.Data()))continue;
        TFile* f=TFile::Open(p.Data()); TTree* t=(TTree*)f->Get("rad"); if(!t){f->Close();continue;}
        RadView v; v.attach(t,&cfg);
        add(bD ,e,sigmaSubset(v,e,allD,allU,4,RadView::kLGCFD));
        add(bD2,e,sigmaSubset(v,e,dsbD,dsbU,2,RadView::kLGCFD));   // NE,SW diagonal control
        f->Close(); } }
    { BuildConfig cfg=BuildConfig::Load(radConfig("LUAG").Data());
      for(double e:Es){ TString p=radReduced("LUAG",e); if(gSystem->AccessPathName(p.Data()))continue;
        TFile* f=TFile::Open(p.Data()); TTree* t=(TTree*)f->Get("rad"); if(!t){f->Close();continue;}
        RadView v; v.attach(t,&cfg);
        add(bL ,e,sigmaSubset(v,e,allD,allU,4,RadView::kLED));
        add(bL2,e,sigmaSubset(v,e,luaD,luaU,2,RadView::kLED));     // NW,SE diagonal control
        f->Close(); } }
    fitAB(mD); fitAB(mL); fitAB(bD); fitAB(bL); fitAB(bD2); fitAB(bL2);

    printf("\n========== MIXED by material vs pure-build 2-corner controls (same geometry) ==========\n");
    printf("  E    mixDSB1  pureDSB1-2c  | mixLuAG  pureLuAG-2c  ||  pureDSB1-4c  pureLuAG-4c\n");
    for(double e:Es){ auto g=[&](Ld&L){ for(size_t j=0;j<L.E.size();++j) if(std::fabs(L.E[j]-e)<1) return L.S[j]; return -1.0; };
        double a=g(mD),a2=g(bD2),b=g(mL),b2=g(bL2),c=g(bD),d=g(bL);
        printf("  %3.0f  %7.1f  %10.1f  | %7.1f  %10.1f  ||  %10.1f  %10.1f\n",e,a,a2,b,b2,c,d); }
    printf("  MIXED-DSB1 (NE,SW lgcfd):       b=%.1f+-%.1f chi2/ndf=%.1f\n",mD.b,mD.be,mD.chi);
    printf("  pure DSB1  (NE,SW lgcfd, ctrl): b=%.1f+-%.1f chi2/ndf=%.1f   <- 2-corner geometry penalty\n",bD2.b,bD2.be,bD2.chi);
    printf("  MIXED-LuAG (NW,SE led):         b=%.1f+-%.1f chi2/ndf=%.1f\n",mL.b,mL.be,mL.chi);
    printf("  pure LuAG  (NW,SE led, ctrl):   b=%.1f+-%.1f chi2/ndf=%.1f   <- 2-corner geometry penalty\n",bL2.b,bL2.be,bL2.chi);
    printf("  pure DSB1 (4-corner):           b=%.1f+-%.1f chi2/ndf=%.1f\n",bD.b,bD.be,bD.chi);
    printf("  pure LuAG (4-corner):           b=%.1f+-%.1f chi2/ndf=%.1f\n",bL.b,bL.be,bL.chi);

    TCanvas* c=new TCanvas("ms","",920,680); c->SetLeftMargin(0.12); c->SetRightMargin(0.04);
    c->SetTopMargin(0.06); c->SetBottomMargin(0.13); c->SetGridy();
    TH1F* fr=c->DrawFrame(0,0,165,80);
    fr->SetTitle(";beam energy E (GeV);brightest-1000  (DW#minusUP)/2  #sigma_{t} (ps)");
    fr->GetYaxis()->SetTitleSize(0.044); fr->GetYaxis()->SetTitleOffset(1.25); fr->GetXaxis()->SetTitleSize(0.044);
    auto draw=[&](Ld&L,int col,int mk,int ls,double sz){ if(L.E.empty())return (TGraphErrors*)nullptr;
        TGraphErrors* g=new TGraphErrors(L.E.size(),&L.E[0],&L.S[0],&L.ze[0],&L.Se[0]);
        g->SetMarkerStyle(mk); g->SetMarkerColor(col); g->SetLineColor(col); g->SetMarkerSize(sz);
        if(L.a>0&&L.chi<8){TF1* f=new TF1(Form("f%d%d",col,mk),"sqrt([0]*[0]/x+[1]*[1])",20,160);f->SetParameters(L.a,L.b);f->SetLineColor(col);f->SetLineStyle(ls);f->SetLineWidth(ls==1?3:2);f->Draw("SAME");}
        g->Draw("P SAME"); return g; };
    TGraphErrors* gMD =draw(mD ,kAzure+2,20,1,1.8);   // MIXED DSB1 corners (filled)
    TGraphErrors* gBD2=draw(bD2,kAzure+2,24,7,1.4);   // pure DSB1 same diagonal (open, dotted)
    TGraphErrors* gML =draw(mL ,kGreen+3,21,1,1.8);   // MIXED LuAG corners (filled)
    TGraphErrors* gBL2=draw(bL2,kGreen+3,25,7,1.4);   // pure LuAG same diagonal (open, dotted)
    TLegend* lg=new TLegend(0.34,0.63,0.96,0.91); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.032);
    if(gMD)  lg->AddEntry(gMD ,Form("MIXED DSB1 corners NE,SW (lgcfd):  b=%.0f#pm%.0f ps",mD.b,mD.be),"lp");
    if(gBD2) lg->AddEntry(gBD2,Form("pure DSB1, same NE,SW diagonal:  b=%.0f ps",bD2.b),"lp");
    if(gML)  lg->AddEntry(gML ,Form("MIXED LuAG corners NW,SE (led):  b=%.0f#pm%.0f ps",mL.b,mL.be),"lp");
    if(gBL2) lg->AddEntry(gBL2,Form("pure LuAG, same NW,SE diagonal:  b=%.0f ps",bL2.b),"lp");
    lg->Draw();
    TLatex t; t.SetNDC(); t.SetTextSize(0.028); t.SetTextColor(kGray+3);
    t.DrawLatex(0.14,0.21,"Filled = MIXED material diagonal;  open = pure build, SAME 2-corner diagonal (geometry control).");
    t.DrawLatex(0.14,0.165,"Gap between filled & open at fixed colour = MIXED-specific run-merge / beam-drift penalty.");
    c->Print(radFigP("figures/narrative/mixed_separate.png"));
    printf("  wrote figures/narrative/mixed_separate.png\n");
}

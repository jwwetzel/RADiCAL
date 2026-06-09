// ============================================================================
// floorModel.C — the timing floor is MODEL-DEPENDENT. Fit the brightest-1000
// sigma_t(E) ladder with BOTH:
//   PHOTOSTAT: sigma_t = a/sqrt(E) (+) b   (the published NIM A 1068 form -> ~17.5 ps)
//   SLEW:      sigma_t = a/E       (+) b   (sigma ~ 1/light -> ~25 ps)
// The within-energy amplitude data (sigma ~ 1/amplitude, pub_res_*) favors SLEW, so
// the ~25 ps floor is the physical one; 17.5 ps is the 1/sqrt(E) artifact.
//   source setup.sh; root -l -b -q 'analyze/studies/floorModel.C+("DSB1")'
// ============================================================================
#include "RadView.h"
#include "RadTiming.h"
#include "SelectionCuts.h"
#include "DataPaths.h"
#include "PlotUtils.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include <vector>
#include <cmath>
#include <cstdio>
#include <string>
using namespace rad;

static int robustSrc(const char* b){ std::string s=b; if(s=="LUAG"||s=="TENERGY") return RadView::kLED; return RadView::kLGCFD; }
static const char* srcName(int s){ return s==RadView::kLED?"led":(s==RadView::kLGCFD?"lgcfd":"cfd05"); }

struct FitOut { double b, be, chi; };
static FitOut fitForm(TGraphErrors* g, const char* form, double b0){
    TF1 f("f",form,20,160); f.SetParameters(200,b0);
    g->Fit(&f,"QN");
    double cn=f.GetChisquare()/std::max(1,f.GetNDF());
    double be=f.GetParError(1); if(cn>1) be*=std::sqrt(cn);
    return {std::fabs(f.GetParameter(1)), be, cn};
}

void floorModel(const char* build="DSB1"){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    BuildConfig cfg=BuildConfig::Load(radConfig(build).Data());
    int SRC=robustSrc(build);
    const double Es[]={25,50,75,100,125,150};
    std::vector<double> E,S,Se,ze;
    printf("\n=== %s (%s): brightest-1000 sigma_t(E) ===\n",build,srcName(SRC));
    for(double e:Es){ TFile* fp=TFile::Open(radReduced(build,e)); if(!fp||fp->IsZombie())continue;
        TTree* t=(TTree*)fp->Get("rad"); RadView v; v.attach(t,&cfg);
        TimingResult r=timingBrightestK(v,e,SRC,1000);
        if(r.sigma_ps>0){ E.push_back(e); S.push_back(r.sigma_ps); Se.push_back(r.sigma_ps/std::sqrt(2.0*1000.0)); ze.push_back(0);
            printf("  %3.0f GeV : %.1f ps\n",e,r.sigma_ps); }
        fp->Close();
    }
    if(E.size()<3){ printf("too few points\n"); return; }
    TGraphErrors* g=new TGraphErrors(E.size(),&E[0],&S[0],&ze[0],&Se[0]);
    FitOut ph=fitForm(g,"sqrt([0]*[0]/x+[1]*[1])",18);       // photostat 1/sqrt(E)
    FitOut sl=fitForm(g,"sqrt([0]*[0]/(x*x)+[1]*[1])",24);   // slew 1/E (1/light)
    printf("  PHOTOSTAT  a/sqrt(E)(+)b : floor b = %.1f +- %.1f ps  (chi2/ndf=%.2f)\n",ph.b,ph.be,ph.chi);
    printf("  SLEW       a/E      (+)b : floor b = %.1f +- %.1f ps  (chi2/ndf=%.2f)\n",sl.b,sl.be,sl.chi);

    // figure
    TCanvas* c=new TCanvas("fm","",900,650);
    c->SetLeftMargin(0.12); c->SetRightMargin(0.04); c->SetTopMargin(0.09); c->SetBottomMargin(0.13); c->SetGridy();
    g->SetMarkerStyle(20); g->SetMarkerColor(kBlack); g->SetMarkerSize(1.5); g->SetLineColor(kBlack); g->SetLineWidth(2);
    g->SetTitle(";beam energy E (GeV);brightest-1000  (DW#minusUP)/2  #sigma_{t} (ps)");
    g->GetXaxis()->SetLimits(0,165); g->GetYaxis()->SetRangeUser(0,std::max(50.0,S[0]*1.15));
    g->GetYaxis()->SetTitleSize(0.045); g->GetYaxis()->SetTitleOffset(1.25); g->GetXaxis()->SetTitleSize(0.045);
    g->Draw("AP");
    TF1* fph=new TF1("fph","sqrt([0]*[0]/x+[1]*[1])",20,160); fph->SetParameters(200,ph.b);
    g->Fit(fph,"QN"); fph->SetLineColor(kAzure+1); fph->SetLineWidth(3); fph->SetLineStyle(1); fph->Draw("SAME");
    TF1* fsl=new TF1("fsl","sqrt([0]*[0]/(x*x)+[1]*[1])",20,160); g->Fit(fsl,"QN"); fsl->SetLineColor(kRed+1); fsl->SetLineWidth(3); fsl->SetLineStyle(1); fsl->Draw("SAME");
    TLine* lp=new TLine(0,ph.b,165,ph.b); lp->SetLineColor(kAzure+1); lp->SetLineStyle(2); lp->SetLineWidth(2); lp->Draw();
    TLine* ls=new TLine(0,sl.b,165,sl.b); ls->SetLineColor(kRed+1); ls->SetLineStyle(2); ls->SetLineWidth(2); ls->Draw();
    g->Draw("P SAME");
    TLegend* lg=new TLegend(0.40,0.66,0.95,0.88); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.036);
    lg->AddEntry(fph,Form("photostat  a/#sqrt{E} #oplus b   #rightarrow  floor %.1f #pm %.1f ps",ph.b,ph.be),"l");
    lg->AddEntry(fsl,Form("slew  a/E #oplus b   #rightarrow  floor %.1f #pm %.1f ps",sl.b,sl.be),"l");
    lg->Draw();
    TLatex tt; tt.SetNDC(); tt.SetTextFont(62); tt.SetTextSize(0.036);
    tt.DrawLatex(0.12,0.945,Form("%s (%s): the floor is the FIT FORM -- 1/#sqrt{E} gives ~%.0f ps, 1/light gives ~%.0f ps",build,srcName(SRC),ph.b,sl.b));
    gSystem->mkdir("figures/narrative",kTRUE); c->Print(Form("figures/narrative/floor_model_%s.png",build));
    printf("  wrote figures/narrative/floor_model_%s.png\n",build);
}

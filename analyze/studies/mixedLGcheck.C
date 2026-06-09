// ============================================================================
// mixedLGcheck.C — is MIXED's low-gain readout "toast", or just dim for LuAG?
// The decisive test: per channel, does HG_peak still track LG_peak linearly?
// A WORKING LG channel shows a tight HG = a + b*LG correlation (this IS the
// calibration lgcfd needs); a dead/decorrelated channel shows a low correlation
// factor and a flat/garbage scatter. We overlay DSB1 (known-good) for reference.
// Canonical end order: NW-D NE-D SE-D SW-D NW-U NE-U SE-U SW-U.
// MIXED materials: NE/SW = DSB1 (bright), NW/SE = LuAG (dim).
//   source setup.sh; root -l -b -q 'analyze/studies/mixedLGcheck.C+("MIXED",150)'
// ============================================================================
#include "RadView.h"
#include "RadTiming.h"
#include "SelectionCuts.h"
#include "DataPaths.h"
#include "PlotUtils.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TF1.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include "FigPaths.h"
#include <vector>
#include <cmath>
#include <cstdio>
#include <string>
using namespace rad;

void mixedLGcheck(const char* build="MIXED", double E=150){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    BuildConfig cfg=BuildConfig::Load(radConfig(build).Data());
    const char* end[8]={"NW-D","NE-D","SE-D","SW-D","NW-U","NE-U","SE-U","SW-U"};
    // MIXED material per canonical end (NW/SE=LuAG, NE/SW=DSB1); DSB1 = all LYSO
    std::string b=build; bool isMix=(b=="MIXED");
    auto mat=[&](int c)->const char*{ if(!isMix) return "LYSO";
        // corner = end[c] prefix
        std::string e=end[c]; std::string cor=e.substr(0,2);
        return (cor=="NE"||cor=="SW")?"DSB1":"LuAG"; };

    TFile* fp=TFile::Open(radReduced(build,E)); if(!fp||fp->IsZombie()){printf("no file\n");return;}
    TTree* t=(TTree*)fp->Get("rad"); RadView v; v.attach(t,&cfg);
    double xc,yc; v.beamCenter(xc,yc); double r2=9.0; Long64_t N=v.entries();
    std::vector<TH2F*> h(8);
    for(int c=0;c<8;++c){ h[c]=new TH2F(Form("h%d",c),Form("%s (%s);LG peak (mV);HG peak (mV)",end[c],mat(c)),
        80,0,400,80,0,900); h[c]->SetDirectory(nullptr); }
    for(Long64_t i=0;i<N;++i){ v.get(i);
        if(!v.wc_ok()||v.mcp1_peak()<kMCP1_minPeak||v.mcp1_peak()>kMCP1_maxPeak) continue;
        double dx=v.x_trk()-xc,dy=v.y_trk()-yc; if(dx*dx+dy*dy>=r2) continue;
        for(int c=0;c<8;++c){ float hg=v.hg_peak(c), lg=v.lg_peak(c); if(lg>3&&hg>10) h[c]->Fill(lg,hg); }
    }
    fp->Close();

    TCanvas* cv=new TCanvas("lg","",1600,820); cv->Divide(4,2,0.004,0.03);
    printf("\n=== %s @%.0f GeV: per-channel HG-vs-LG correlation (LG health) ===\n",build,E);
    printf("  end    material   N      corr    slope b   (working if corr>~0.8)\n");
    for(int c=0;c<8;++c){ cv->cd(c+1); gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.14); gPad->SetRightMargin(0.04);
        double corr=h[c]->GetCorrelationFactor(); long n=(long)h[c]->GetEntries();
        h[c]->SetMarkerStyle(1);
        const char* m=mat(c); int col=(std::string(m)=="LuAG")?(kGreen+3):(kAzure+2);
        h[c]->SetMarkerColor(col); h[c]->GetXaxis()->SetTitleSize(0.06); h[c]->GetYaxis()->SetTitleSize(0.06);
        h[c]->GetXaxis()->SetLabelSize(0.05); h[c]->GetYaxis()->SetLabelSize(0.05);
        h[c]->Draw("");
        TProfile* pr=h[c]->ProfileX(Form("p%d",c)); pr->SetLineColor(kRed+1); pr->SetLineWidth(2); pr->SetMarkerColor(kRed+1); pr->SetMarkerStyle(20); pr->SetMarkerSize(0.5); pr->Draw("SAME");
        double slope=0; TF1 f("f","[0]+[1]*x",10,300); if(n>50){ pr->Fit(&f,"QNR"); slope=f.GetParameter(1); }
        TLatex tx; tx.SetNDC(); tx.SetTextSize(0.062); tx.SetTextColor(kBlack);
        tx.DrawLatex(0.20,0.86,Form("#rho=%.2f",corr));
        tx.DrawLatex(0.20,0.79,Form("b=%.1f",slope));
        printf("  %-5s  %-7s  %7ld  %5.2f   %6.1f\n",end[c],m,n,corr,slope);
    }
    cv->cd(0); TLatex tt; tt.SetNDC(); tt.SetTextFont(62); tt.SetTextSize(0.028);
    tt.DrawLatex(0.04,0.965,Form("%s @%.0f GeV: HG peak vs LG peak per channel  (#rho = correlation; working LG #Rightarrow tight linear band)  [blue=DSB1/LYSO, green=LuAG]",build,E));
    gSystem->mkdir(Form("figures/%d/narrative",radYear()),kTRUE); cv->Print(radFigP(Form("figures/narrative/lgcheck_%s.png",build)));
    printf("  wrote figures/narrative/lgcheck_%s.png\n",build);
}

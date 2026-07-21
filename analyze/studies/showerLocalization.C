// showerLocalization.C — Paper 2 (energy+position): reconstruct the transverse
//   EM-shower impact point from the FOUR corner-capillary amplitudes and measure
//   its resolution against the wire-chamber truth. The spatial half of the
//   shower-max measurement (published 2401.01747 Section 7 item 1).
//
//   Rather than assume the module->WC orientation, we fit the best LINEAR
//   estimator of the WC position from the four corner light fractions (energy per
//   corner = up+down low-gain peak), and take the residual as the position
//   resolution. This is geometry-agnostic and gives the optimal linear localizer.
//   root -l 'Analysis/showerLocalization.C+(150)'
#include "PlotUtils.h"
#include "TFile.h"
#include "DataPaths.h"
#include "TTree.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TLinearFitter.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "FigPaths.h"
#include <cstdio>
#include <cmath>

void showerLocalization(double E=150){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    TFile* fp=TFile::Open(radReduced("DSB1",E));
    if(!fp||fp->IsZombie()){ printf("no DSB1 ntuple at %.0f GeV\n",E); return; }
    TTree* t=(TTree*)fp->Get("rad");
    Bool_t wc; Float_t x,y,mp,lg[8];
    t->SetBranchAddress("wc_ok",&wc); t->SetBranchAddress("x_trk",&x); t->SetBranchAddress("y_trk",&y);
    t->SetBranchAddress(t->GetBranch("mcp1_peak")?"mcp1_peak":"mcp_peak",&mp); t->SetBranchAddress("lg_peak",lg);
    const double cx0=6.6, cy0=4.7;                 // module centre in WC coords

    // best linear estimator of WC (x,y) from 3 independent corner fractions (the
    // 4th is fixed by sum=1). hyp3 = const + 3 linear terms.
    TLinearFitter lx(3,"hyp3"), ly(3,"hyp3");
    TH2F* hx=new TH2F("hx","x: capillary estimator vs wire chamber;x_{WC} (mm);x_{capillary} (mm)",70,-7,7,70,-7,7);
    TH2F* hy=new TH2F("hy","y: capillary estimator vs wire chamber;y_{WC} (mm);y_{capillary} (mm)",70,-7,7,70,-7,7);
    hx->SetDirectory(0); hy->SetDirectory(0);
    std::vector<float> f0v,f1v,f2v,xwv,ywv;
    long N=t->GetEntries(); double sxw=0,sxw2=0; long n=0;
    for(long i=0;i<N;++i){ t->GetEntry(i);
        if(!wc||mp<200||mp>750) continue;
        double xw=x-cx0, yw=y-cy0; if(std::fabs(xw)>6||std::fabs(yw)>6) continue;
        double A[4]={ (double)lg[0]+lg[4], (double)lg[1]+lg[5], (double)lg[2]+lg[6], (double)lg[3]+lg[7] };
        double s=A[0]+A[1]+A[2]+A[3]; if(s<500) continue;
        double f[3]={A[0]/s,A[1]/s,A[2]/s};
        lx.AddPoint(f,xw,1.0); ly.AddPoint(f,yw,1.0);
        f0v.push_back(f[0]); f1v.push_back(f[1]); f2v.push_back(f[2]); xwv.push_back(xw); ywv.push_back(yw);
        sxw+=xw; sxw2+=xw*xw; ++n;
    }
    fp->Close();
    if(n<2000){ printf("low stats (%ld)\n",n); return; }
    lx.Eval(); ly.Eval();
    double resx=std::sqrt(lx.GetChisquare()/n), resy=std::sqrt(ly.GetChisquare()/n);  // residual vs WC, mm
    double cx[4]={lx.GetParameter(0),lx.GetParameter(1),lx.GetParameter(2),lx.GetParameter(3)};
    double cy[4]={ly.GetParameter(0),ly.GetParameter(1),ly.GetParameter(2),ly.GetParameter(3)};
    double spread=std::sqrt(sxw2/n-(sxw/n)*(sxw/n));   // beam-spot spread (no-info baseline)
    // fill reco-vs-WC maps
    for(size_t i=0;i<xwv.size();++i){
        double xr=cx[0]+cx[1]*f0v[i]+cx[2]*f1v[i]+cx[3]*f2v[i];
        double yr=cy[0]+cy[1]*f0v[i]+cy[2]*f1v[i]+cy[3]*f2v[i];
        hx->Fill(xwv[i],xr); hy->Fill(ywv[i],yr);
    }
    printf("\n=== shower localization from 4-capillary light (DSB1, %.0f GeV, N=%ld) ===\n",E,n);
    printf("  beam-spot spread, no-info baseline (sigma of WC positions in-cut): %.2f mm\n",spread);
    printf("  x: residual (capillary estimator - WC) = %.2f mm\n",resx);
    printf("  y: residual (capillary estimator - WC) = %.2f mm\n",resy);
    printf("  => the 4-capillary light localizes the shower to ~%.1f mm, well below the %.1f mm spot\n",0.5*(resx+resy),spread);
    printf("     spread (real position info). This residual = capillary (+) WC truth in quadrature, so\n");
    printf("     both the capillary and the WC resolution are <= ~%.1f mm here.\n",0.5*(resx+resy));

    TCanvas* c=new TCanvas("c_loc","",1300,560); c->Divide(2,1,0.006,0.006);
    c->cd(1); gPad->SetRightMargin(0.13); hx->Draw("COLZ"); { TProfile* p=hx->ProfileX("px"); p->SetMarkerStyle(20); p->SetLineColor(kBlack); p->Draw("SAME"); }
    c->cd(2); gPad->SetRightMargin(0.13); hy->Draw("COLZ"); { TProfile* p=hy->ProfileX("py"); p->SetMarkerStyle(20); p->SetLineColor(kBlack); p->Draw("SAME"); }
    c->cd(0); TLatex tl; tl.SetNDC(); tl.SetTextSize(0.026);
    tl.DrawLatex(0.06,0.955,Form("Shower position from 4-capillary light vs wire chamber (DSB1, %.0f GeV): residual %.1f / %.1f mm (x/y); spot spread %.1f mm",E,resx,resy,spread));
    c->Print(radFigP("figures/narrative/shower_localization.png"));
    printf("wrote figures/narrative/shower_localization.png\n");
}

// compareConfigsPlot.C — cross-config timing & energy resolution vs energy.
// First-look numbers (all-fiducial, (DW-UP)/2 CFD-5%, truncated-core), produced by
// configResolution.C (LuAG/MIXED/TENERGY, reduced) + configResolutionDSB1.C (DSB1,
// processRun).  NOT best-bin / OOS — a same-method comparison, not the headline.
#include "PlotUtils.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TLatex.h"

void compareConfigsPlot(){
  ApplyRADiCALStyle();
  // E grids
  double Ed[6]={25,50,75,100,125,150};               // DSB1 (has 25)
  double Eo[5]={50,75,100,125,150};                  // others (no 25)
  double stD[6]={70.2,60.6,59.1,55.4,56.2,58.4},  seD[6]={19.5,17.5,16.6,16.7,16.8,16.2};
  double stL[5]={72.4,63.4,61.8,55.5,49.8},        seL[5]={14.0,13.5,13.6,14.0,14.5};
  double stM[5]={66.0,60.8,58.8,56.6,54.5},        seM[5]={15.2,14.6,14.7,14.7,14.8};
  double stT[5]={66.6,62.3,63.8,59.4,57.9},        seT[5]={14.9,14.3,14.4,14.2,13.5};
  int cD=kRData, cL=kRGreen+1, cM=kROrange, cT=kRRed;

  auto mk=[&](int n,double*E,double*Y,int col,int mst){TGraph*g=new TGraph(n,E,Y);
    g->SetLineColor(col);g->SetMarkerColor(col);g->SetMarkerStyle(mst);g->SetMarkerSize(1.5);g->SetLineWidth(3);return g;};

  TCanvas*c=new TCanvas("c_cmp","",1340,560); c->Divide(2,1,0.015,0.02);

  // LEFT: timing
  c->cd(1); gPad->SetLeftMargin(0.13);gPad->SetBottomMargin(0.13);gPad->SetTopMargin(0.10);
  TH1F*f1=gPad->DrawFrame(0,40,165,80);
  f1->GetXaxis()->SetTitle("beam energy (GeV)");f1->GetYaxis()->SetTitle("#sigma_{t} (DW#minusUP)/2 (ps)");
  f1->GetXaxis()->SetTitleSize(0.047);f1->GetYaxis()->SetTitleSize(0.047);
  mk(6,Ed,stD,cD,20)->Draw("PL"); mk(5,Eo,stL,cL,21)->Draw("PL");
  mk(5,Eo,stM,cM,22)->Draw("PL"); mk(5,Eo,stT,cT,23)->Draw("PL");
  {TLegend*L=new TLegend(0.52,0.66,0.95,0.89);L->SetBorderSize(0);L->SetFillStyle(0);L->SetTextSize(0.036);
   L->AddEntry(mk(6,Ed,stD,cD,20),"DSB1 (saturates)","pl");
   L->AddEntry(mk(5,Eo,stL,cL,21),"LuAG","pl");
   L->AddEntry(mk(5,Eo,stM,cM,22),"2xDSB1+2xLuAG","pl");
   L->AddEntry(mk(5,Eo,stT,cT,23),"3xDSB1+1xEnergy","pl");L->Draw();}
  DrawPadTitle("Timing: LuAG wins at high E (no saturation)",0.05);

  // RIGHT: energy
  c->cd(2); gPad->SetLeftMargin(0.13);gPad->SetBottomMargin(0.13);gPad->SetTopMargin(0.10);
  TH1F*f2=gPad->DrawFrame(0,10,165,22);
  f2->GetXaxis()->SetTitle("beam energy (GeV)");f2->GetYaxis()->SetTitle("#sigma_{E}/E (%)");
  f2->GetXaxis()->SetTitleSize(0.047);f2->GetYaxis()->SetTitleSize(0.047);
  mk(6,Ed,seD,cD,20)->Draw("PL"); mk(5,Eo,seL,cL,21)->Draw("PL");
  mk(5,Eo,seM,cM,22)->Draw("PL"); mk(5,Eo,seT,cT,23)->Draw("PL");
  {TLatex t;t.SetTextSize(0.030);t.SetTextColor(kGray+3);
   t.DrawLatex(40,12.2,"leakage-dominated (short 4 X_{0} prototype):");
   t.DrawLatex(40,11.3,"flat constant term ~14#minus17% for all configs");}
  DrawPadTitle("Energy: all leakage-limited (~14#minus17%)",0.05);

  c->cd(0);
  DrawPageTitle("Capillary configs: timing & energy resolution (first-look, all-fiducial, same method)");
  c->Print("/tmp/compare_configs.png");
  printf("wrote /tmp/compare_configs.png\n");
}

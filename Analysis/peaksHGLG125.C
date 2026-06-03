// peaksHGLG125.C — reproduce the legacy "125 GeV peaks" plot (HG peak vs LG peak,
// 8 capillaries) and measure the scatter of the UNSATURATED linear branch at high
// energy (where de-saturation operates), to compare against the noisy low-E fit.
//   root -l 'Analysis/peaksHGLG125.C+'
#include "PlotUtils.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TF1.h"
#include "TLatex.h"
#include <cstdio>
#include <vector>
#include <cmath>
void peaksHGLG125(){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    const char* nm[8]={"NW-D","NE-D","SE-D","SW-D","NW-U","NE-U","SE-U","SW-U"};
    TFile* fp=TFile::Open("Analysis/Output/125GeV/ntuple.root");
    if(!fp||fp->IsZombie()){printf("no 125 GeV ntuple\n");return;}
    TTree* t=(TTree*)fp->Get("rad");
    Bool_t wc,fid;Float_t hp[8],lp[8];
    t->SetBranchAddress("wc_ok",&wc);t->SetBranchAddress("in_fiducial",&fid);
    t->SetBranchAddress("hg_peak",hp);t->SetBranchAddress("lg_peak",lp);
    TH2F* h[8];for(int i=0;i<8;++i){h[i]=new TH2F(Form("h%d",i),Form("%s;LG peak (mV);HG peak (mV)",nm[i]),100,0,500,100,0,1000);h[i]->SetDirectory(0);}
    std::vector<float> ux[8],uy[8];   // unsaturated branch (HG<780) for scatter
    long N=t->GetEntries();
    for(long e=0;e<N;++e){t->GetEntry(e);if(!wc||!fid)continue;   // <-- CONTAINED showers only
        for(int i=0;i<8;++i){ if(lp[i]>2&&hp[i]>2){h[i]->Fill(lp[i],hp[i]);
            if(hp[i]<780&&hp[i]>40&&lp[i]>5){ux[i].push_back(lp[i]);uy[i].push_back(hp[i]);}}}}
    fp->Close();
    printf("\n=== 125 GeV HG-vs-LG peak: UNSATURATED-branch linear scatter ===\n");
    for(int i=0;i<8;++i){ long n=ux[i].size(); if(n<300){printf("%-6s low stats\n",nm[i]);continue;}
        double sx=0,sy=0,sxx=0,sxy=0,syy=0;for(long k=0;k<n;++k){sx+=ux[i][k];sy+=uy[i][k];sxx+=ux[i][k]*ux[i][k];sxy+=ux[i][k]*uy[i][k];syy+=uy[i][k]*uy[i][k];}
        double mx=sx/n,my=sy/n,cxy=sxy/n-mx*my,vx=sxx/n-mx*mx,vy=syy/n-my*my,b=cxy/vx,r=cxy/std::sqrt(vx*vy);
        double sc=100*std::sqrt((vy-b*cxy)>0?(vy-b*cxy):0)/my;
        printf("%-6s  r=%.4f  slope=%.2f  linear-branch scatter=%.1f%%\n",nm[i],r,b,sc);}
    TCanvas* c=new TCanvas("c_pk","",1500,760);c->Divide(4,2,0.004,0.03);
    for(int i=0;i<8;++i){c->cd(i+1);gPad->SetRightMargin(0.02);h[i]->Draw("COL");
        DrawPadTitle(Form("%s  HG vs LG peak (125 GeV)",nm[i]),0.07);}
    c->Print("Analysis/capillary_figs/peaks_hglg_125.png");
    printf("\nwrote Analysis/capillary_figs/peaks_hglg_125.png\n");
}

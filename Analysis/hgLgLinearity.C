// hgLgLinearity.C — the HG-vs-LG dual-readout linearity (DSB1), the basis for an
// LG-referenced de-saturation of the HG CFD threshold. Fit a line to UNSATURATED
// HG events; report correlation + fractional scatter for BOTH peak-vs-peak and
// charge-vs-charge (peak ratio can scatter from HG/LG bandwidth differences even
// for the same signal; charge is bandwidth-invariant).
//   root -l 'Analysis/hgLgLinearity.C+'
#include "PlotUtils.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TF1.h"
#include "TH2F.h"
#include "TLatex.h"
#include <cstdio>
#include <vector>
#include <cmath>
static void fitline(std::vector<float>&X,std::vector<float>&Y,double&b,double&r,double&relscat,double&my){
    double sx=0,sy=0,sxx=0,sxy=0,syy=0;long n=X.size();
    for(long i=0;i<n;++i){sx+=X[i];sy+=Y[i];sxx+=X[i]*X[i];sxy+=X[i]*Y[i];syy+=Y[i]*Y[i];}
    double mx=sx/n;my=sy/n;double cxy=sxy/n-mx*my,vx=sxx/n-mx*mx,vy=syy/n-my*my;
    b=cxy/vx;r=cxy/std::sqrt(vx*vy);double resid=vy-b*cxy;relscat=100*std::sqrt(resid>0?resid:0)/my;
}
void hgLgLinearity(){
    ApplyRADiCALStyle();
    int hgS[8],lgS[8];const char* nm[8]={"NW-D","NE-D","SE-D","SW-D","NW-U","NE-U","SE-U","SW-U"};
    // slot indices via the kCap order from ChannelConfig is not included here; use the known map:
    // HG timing slots (DRS0): 1,2,3,0,5,4,6,9 ; LG energy slots (DRS1 G0): 19,20,21,18,23,22,24,25
    int hgmap[8]={1,2,3,0,5,4,6,9}, lgmap[8]={19,20,21,18,23,22,24,25};
    for(int i=0;i<8;++i){hgS[i]=hgmap[i];lgS[i]=lgmap[i];}
    // collect unsaturated (lg_peak, hg_peak) and (lg_charge, hg_charge) per channel from 25+50 GeV
    std::vector<float> px[8],py[8],cx[8],cy[8];
    for(double E:{25.,50.}){ TFile* fp=TFile::Open(Form("Analysis/Output/%.0fGeV/ntuple.root",E));
        if(!fp||fp->IsZombie())continue; TTree* t=(TTree*)fp->Get("rad");if(!t){fp->Close();continue;}
        Bool_t wc;Float_t hp[8],lp[8],hc[8],lc[8];
        t->SetBranchAddress("wc_ok",&wc);t->SetBranchAddress("hg_peak",hp);t->SetBranchAddress("lg_peak",lp);
        t->SetBranchAddress("hg_charge",hc);t->SetBranchAddress("lg_charge",lc);
        long N=t->GetEntries();
        for(long e=0;e<N;++e){t->GetEntry(e);if(!wc)continue;
            for(int i=0;i<8;++i){ if(hp[i]>30&&hp[i]<680&&lp[i]>10){px[i].push_back(lp[i]);py[i].push_back(hp[i]);}
                                  if(hc[i]>0&&lc[i]>0){cx[i].push_back(lc[i]);cy[i].push_back(hc[i]);} }}
        fp->Close(); }
    printf("\n=== DSB1 HG-vs-LG dual-readout linearity (unsaturated events) ===\n");
    printf("%-6s %22s %22s\n","chan","PEAK: r / slope / scatter","CHARGE: r / slope / scatter");
    for(int i=0;i<8;++i){ double bp,rp,sp,myp,bc,rc,sc,myc;
        if(px[i].size()<500||cx[i].size()<500){printf("%-6s (low stats)\n",nm[i]);continue;}
        fitline(px[i],py[i],bp,rp,sp,myp); fitline(cx[i],cy[i],bc,rc,sc,myc);
        printf("%-6s   r=%.3f b=%5.2f sc=%4.1f%%      r=%.3f b=%5.2f sc=%4.1f%%\n",nm[i],rp,bp,sp,rc,bc,sc);}
    // scatter figure for NE-D (peak) and its charge
    TCanvas* c=new TCanvas("c_lin","",1200,540);c->Divide(2,1);
    c->cd(1);gPad->SetLeftMargin(0.14);gPad->SetBottomMargin(0.13);
    TH2F* h1=new TH2F("h1","NE-D  HG peak vs LG peak (unsat);LG peak (mV);HG peak (mV)",80,0,160,80,0,700);
    for(size_t k=0;k<px[1].size();++k)h1->Fill(px[1][k],py[1][k]); h1->Draw("COLZ");
    {double b,r,s,my;fitline(px[1],py[1],b,r,s,my);TF1*f=new TF1("f","[0]*x",0,160);f->SetParameter(0,b);f->SetLineColor(kRRed);f->SetLineWidth(2);f->Draw("same");
     TLatex t;t.SetNDC();t.SetTextSize(0.04);t.SetTextColor(kRRed);t.DrawLatex(0.18,0.84,Form("r=%.3f  scatter=%.1f%%",r,s));}
    c->cd(2);gPad->SetLeftMargin(0.14);gPad->SetBottomMargin(0.13);
    TH2F* h2=new TH2F("h2","NE-D  HG charge vs LG charge (unsat);LG charge;HG charge",80,0,0,80,0,0);
    h2->SetBins(80,0,*std::max_element(cx[1].begin(),cx[1].end())*1.05,80,0,*std::max_element(cy[1].begin(),cy[1].end())*1.05);
    for(size_t k=0;k<cx[1].size();++k)h2->Fill(cx[1][k],cy[1][k]); h2->Draw("COLZ");
    {double b,r,s,my;fitline(cx[1],cy[1],b,r,s,my);TLatex t;t.SetNDC();t.SetTextSize(0.04);t.SetTextColor(kRRed);t.DrawLatex(0.18,0.84,Form("r=%.3f  scatter=%.1f%%",r,s));}
    DrawPageTitle("DSB1 HG-vs-LG dual-readout linearity (NE-D, unsaturated)");
    c->Print("/tmp/hglg_linearity.png");printf("\nwrote /tmp/hglg_linearity.png\n");
}

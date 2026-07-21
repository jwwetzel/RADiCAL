// ============================================================================
// pbLeakage.C — is the LG sub-linearity at high E longitudinal leakage?
// ----------------------------------------------------------------------------
// The PbGlass blocks sit behind the 25 X0 module and catch longitudinal leakage.
// If the LG (sum_lg) droop is leakage, then (a) sum_pb grows with energy and the
// pb/lg fraction rises, and (b) the leakage-corrected energy sum_lg + alpha*sum_pb
// becomes LINEAR again. We fit the single alpha that flattens the response.
//   ROOT_INCLUDE_PATH=lib/waveform:lib/io:lib/physics:lib/viz root -l -b -q analyze/studies/pbLeakage.C+
// ============================================================================
#include "PlotUtils.h"
#include "SelectionCuts.h"   // kMCP1_minPeak, kMCP_minPeak_E, kSumLG_centroid, kFiducial_r_energy
#include "DataPaths.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "FigPaths.h"
#include <cstdio>
#include <cmath>

void pbLeakage() {
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    double Es[6]={25,50,75,100,125,150};
    double mLG[6], mPB[6];

    for(int e=0;e<6;++e){
        TFile* fp=TFile::Open(radReduced("DSB1",Es[e])); TTree* t=(TTree*)fp->Get("rad");
        Bool_t wc; Float_t x,y,mcp,slg,spb;
        t->SetBranchAddress("wc_ok",&wc); t->SetBranchAddress("x_trk",&x); t->SetBranchAddress("y_trk",&y);
        t->SetBranchAddress(t->GetBranch("mcp1_peak")?"mcp1_peak":"mcp_peak",&mcp); t->SetBranchAddress("sum_lg",&slg); t->SetBranchAddress("sum_pb",&spb);
        long N=t->GetEntries();
        double wx=0,wy=0,w=0; for(long i=0;i<N;++i){ t->GetEntry(i); if(wc&&mcp>kMCP_minPeak_E&&slg>kSumLG_centroid){wx+=x*slg;wy+=y*slg;w+=slg;} }
        double xc=w>0?wx/w:0, yc=w>0?wy/w:0, r2=kFiducial_r_energy*kFiducial_r_energy;
        double sLG=0,sPB=0; long n=0;
        for(long i=0;i<N;++i){ t->GetEntry(i);
            if(!wc||mcp<kMCP1_minPeak||mcp>kMCP1_maxPeak) continue;
            double dx=x-xc,dy=y-yc; if(dx*dx+dy*dy>=r2) continue;
            sLG+=slg; sPB+=spb; ++n;
        }
        mLG[e]=n?sLG/n:0; mPB[e]=n?sPB/n:0;
        printf("E=%3.0f  sum_lg=%.0f  sum_pb=%.0f mV  pb/lg=%.3f\n", Es[e], mLG[e], mPB[e], mPB[e]/mLG[e]);
        fp->Close();
    }

    // PbGlass baseline: sum_pb at 25 GeV is dominated by pedestal/noise (a 25 GeV
    // shower is ~contained in 25 X0), so the LEAKAGE signal is sum_pb - sum_pb[25].
    double pb0=mPB[0];
    // scan alpha to flatten the corrected response (mLG + alpha*leakage_pb)/E
    auto cv=[&](double a)->double{ double r[6],m=0; for(int e=0;e<6;++e){ r[e]=(mLG[e]+a*(mPB[e]-pb0))/Es[e]; m+=r[e]; } m/=6;
        double v=0; for(int e=0;e<6;++e) v+=(r[e]-m)*(r[e]-m); return std::sqrt(v/6)/m; };
    double bestA=0, bestCV=cv(0); for(double a=0; a<=30; a+=0.05){ double q=cv(a); if(q<bestCV){bestCV=q;bestA=a;} }
    double cv0=cv(0);
    printf("\nPbGlass leakage (sum_pb - %.0f baseline): ", pb0); for(int e=0;e<6;++e) printf("%.0f ",mPB[e]-pb0);
    printf("\nresponse non-linearity (CV of response/GeV): uncorrected=%.1f%%  ->  +%.2f*(sum_pb-base) = %.1f%%\n",
           100*cv0, bestA, 100*bestCV);

    // normalized responses (to 25 GeV) for the plot
    double rU[6], rC[6], frac[6];
    double u0=mLG[0]/Es[0], c0=mLG[0]/Es[0];
    for(int e=0;e<6;++e){ rU[e]=(mLG[e]/Es[e])/u0; rC[e]=((mLG[e]+bestA*(mPB[e]-pb0))/Es[e])/c0; frac[e]=mPB[e]/mLG[e]; }

    TCanvas* c=new TCanvas("c_pb","",1500,620); c->Divide(2,1,0.010,0.006);
    // Panel 1: raw sum_lg and sum_pb vs E + pb/lg fraction
    c->cd(1); gPad->SetGridx(); gPad->SetGridy();
    double ymx=0; for(int e=0;e<6;++e) ymx=std::max(ymx,mLG[e]);
    TH1F* fr=gPad->DrawFrame(0,0,160,ymx*1.2);
    fr->SetTitle("Module (LG) energy and PbGlass leakage vs energy;beam energy [GeV];mean amplitude [mV]");
    TGraph* gLG=new TGraph(6,Es,mLG); gLG->SetMarkerStyle(20); gLG->SetMarkerColor(kAzure+2); gLG->SetLineColor(kAzure+2); gLG->SetMarkerSize(1.5); gLG->SetLineWidth(2); gLG->Draw("PL SAME");
    TGraph* gPB=new TGraph(6,Es,mPB); gPB->SetMarkerStyle(21); gPB->SetMarkerColor(kRed+1); gPB->SetLineColor(kRed+1); gPB->SetMarkerSize(1.5); gPB->SetLineWidth(2); gPB->Draw("PL SAME");
    TLegend* lg=new TLegend(0.16,0.70,0.55,0.88); lg->SetBorderSize(0);
    lg->AddEntry(gLG,"sum_lg (module energy)","pl"); lg->AddEntry(gPB,"sum_pb (PbGlass leakage)","pl"); lg->Draw();
    TLatex tl; tl.SetNDC(); tl.SetTextSize(0.027);
    tl.DrawLatex(0.16,0.30,"PbGlass leakage grows with energy as the shower leaks out the back");
    tl.DrawLatex(0.16,0.26,Form("of the 25 X_{0} module (pb/lg: %.2f at 25 -> %.2f at 150 GeV).",frac[0],frac[5]));

    // Panel 2: response/GeV normalized — uncorrected droops, leakage-corrected flat
    c->cd(2); gPad->SetGridx(); gPad->SetGridy();
    TH1F* fr2=gPad->DrawFrame(0,0.7,160,1.15);
    fr2->SetTitle("Energy response per GeV: leakage correction restores linearity;beam energy [GeV];response/GeV (norm. to 25 GeV)");
    TGraph* one=new TGraph(2); one->SetPoint(0,0,1); one->SetPoint(1,160,1); one->SetLineColor(kGray+2); one->SetLineStyle(2); one->Draw("L SAME");
    TGraph* gU=new TGraph(6,Es,rU); gU->SetMarkerStyle(20); gU->SetMarkerColor(kAzure+2); gU->SetLineColor(kAzure+2); gU->SetMarkerSize(1.5); gU->SetLineWidth(2); gU->Draw("PL SAME");
    TGraph* gC=new TGraph(6,Es,rC); gC->SetMarkerStyle(21); gC->SetMarkerColor(kGreen+3); gC->SetLineColor(kGreen+3); gC->SetMarkerSize(1.5); gC->SetLineWidth(2); gC->Draw("PL SAME");
    TLegend* lg2=new TLegend(0.16,0.18,0.62,0.40); lg2->SetBorderSize(0);
    lg2->AddEntry(gU,Form("sum_lg only (droops to %.0f%% at 150)",100*rU[5]),"pl");
    lg2->AddEntry(gC,Form("sum_lg + %.1f*(sum_pb#minusbase): %.0f%%->%.0f%% nonlin",bestA,100*cv0,100*bestCV),"pl");
    lg2->AddEntry(one,"linear","l"); lg2->Draw();

    gSystem->mkdir("figures",kTRUE);
    c->Print(radFigP("figures/pb_leakage.png"));
    printf("wrote figures/pb_leakage.png\n");
}

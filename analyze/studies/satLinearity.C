// ============================================================================
// satLinearity.C — HG (timing) channel clips at the offset-window rail while
// LG (energy) stays LINEAR with energy => the ~800 mV ceiling is the HG READOUT
// window, not SiPM pixel saturation (the SiPM light, witnessed by LG, is fine).
//   ROOT_INCLUDE_PATH=lib/waveform:lib/io:lib/physics:lib/viz root -l -b -q analyze/studies/satLinearity.C+
// ============================================================================
#include "BuildConfig.h"     // radcore: cfg.hg_sat_mV
#include "PlotUtils.h"       // ApplyRADiCALStyle
#include "SelectionCuts.h"
#include "DataPaths.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include <cstdio>

void satLinearity() {
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    rad::BuildConfig cfg = rad::BuildConfig::Load(radConfig("DSB1").Data());
    double Es[6]={25,50,75,100,125,150};
    double mHG[6]={0}, mLG[6]={0};
    long   nCapHit[6]={0}, nCapSat[6]={0};        // HG timing-cap hits + those at/above the clip
    const int CAP=1;                              // NE-D, a bright timing cap
    TH2F* h2=new TH2F("h2","",120,0,1500,120,0,1050); h2->SetDirectory(0);

    for(int e=0;e<6;++e){
        TFile* fp=TFile::Open(radReduced("DSB1",Es[e])); if(!fp||fp->IsZombie()) continue;
        TTree* t=(TTree*)fp->Get("rad");
        Bool_t wc; Float_t x,y,mcp,slg,hgp[8],lgp[8];
        t->SetBranchAddress("wc_ok",&wc); t->SetBranchAddress("x_trk",&x); t->SetBranchAddress("y_trk",&y);
        t->SetBranchAddress("mcp_peak",&mcp); t->SetBranchAddress("sum_lg",&slg);
        t->SetBranchAddress("hg_peak",hgp); t->SetBranchAddress("lg_peak",lgp);
        long N=t->GetEntries();
        double wx=0,wy=0,w=0; for(long i=0;i<N;++i){ t->GetEntry(i); if(wc&&mcp>kMCP_minPeak_E&&slg>kSumLG_centroid){wx+=x*slg;wy+=y*slg;w+=slg;} }
        double xc=w>0?wx/w:0, yc=w>0?wy/w:0, rFid=TimingFiducialR(Es[e]), r2=rFid*rFid;
        double sHG=0,sLG=0; long n=0;
        for(long i=0;i<N;++i){ t->GetEntry(i);
            if(!wc||mcp<kMCP1_minPeak||mcp>kMCP1_maxPeak) continue;
            double dx=x-xc,dy=y-yc; if(dx*dx+dy*dy>=r2) continue;
            double aHG=0,aLG=0; int na=0;
            for(int c=0;c<8;++c) if(hgp[c]>=kHG_minPeak){ aHG+=hgp[c]; aLG+=lgp[c]; ++na;
                ++nCapHit[e]; if(hgp[c]>=cfg.hg_sat_mV) ++nCapSat[e]; }
            if(na>0){ sHG+=aHG/na; sLG+=aLG/na; ++n; }
            if(lgp[CAP]>5 && hgp[CAP]>5) h2->Fill(lgp[CAP],hgp[CAP]);
        }
        mHG[e]=n?sHG/n:0; mLG[e]=n?sLG/n:0;
        printf("E=%3.0f  meanHG=%.0f  meanLG=%.0f mV  HG-cap saturated(>=%.0f)=%.0f%%\n",
               Es[e], mHG[e], mLG[e], cfg.hg_sat_mV, nCapHit[e]?100.0*nCapSat[e]/nCapHit[e]:0);
        fp->Close();
    }

    double nx[6],nhg[6],nlg[6],ideal[6];
    for(int e=0;e<6;++e){ nx[e]=Es[e]; nhg[e]=mHG[e]/mHG[0]; nlg[e]=mLG[e]/mLG[0]; ideal[e]=Es[e]/Es[0]; }

    TCanvas* c=new TCanvas("c_sl","",1500,640); c->Divide(2,1,0.010,0.006);
    // Panel 1: normalized HG vs LG vs energy
    c->cd(1); gPad->SetGridx(); gPad->SetGridy();
    TH1F* fr=gPad->DrawFrame(0,0,160,7);
    fr->SetTitle("HG saturates, LG stays linear (each normalized to its 25 GeV value);beam energy [GeV];response / response(25 GeV)");
    TGraph* gid=new TGraph(6,nx,ideal); gid->SetLineColor(kGray+2); gid->SetLineStyle(2); gid->SetLineWidth(2); gid->Draw("L SAME");
    TGraph* gLG=new TGraph(6,nx,nlg); gLG->SetMarkerStyle(20); gLG->SetMarkerColor(kAzure+2); gLG->SetLineColor(kAzure+2); gLG->SetMarkerSize(1.6); gLG->SetLineWidth(2); gLG->Draw("PL SAME");
    TGraph* gHG=new TGraph(6,nx,nhg); gHG->SetMarkerStyle(21); gHG->SetMarkerColor(kRed+1);  gHG->SetLineColor(kRed+1);  gHG->SetMarkerSize(1.6); gHG->SetLineWidth(2); gHG->Draw("PL SAME");
    TLegend* lg=new TLegend(0.16,0.64,0.58,0.88); lg->SetBorderSize(0);
    lg->AddEntry(gLG,"LG (energy channel) - linear","pl");
    lg->AddEntry(gHG,"HG (timing channel) - saturates","pl");
    lg->AddEntry(gid,"ideal linear (#propto E)","l"); lg->Draw();
    TLatex tl; tl.SetNDC(); tl.SetTextSize(0.026);
    tl.DrawLatex(0.16,0.34,"Same showers: LG tracks energy (#propto E); HG flat-tops by ~50 GeV.");
    tl.DrawLatex(0.16,0.30,"=> DT5742 HG-CHANNEL saturation at the offset window (~800 mV),");
    tl.DrawLatex(0.16,0.26,"   not SiPM saturation (LG, at lower gain, stays unsaturated).");
    tl.DrawLatex(0.16,0.20,Form("HG-cap saturated fraction: 25 GeV %.0f%% (mostly unsaturated) #rightarrow 150 GeV %.0f%%",
                                100.0*nCapSat[0]/nCapHit[0], 100.0*nCapSat[5]/nCapHit[5]));

    // Panel 2: per-event HG vs LG (one bright cap, all energies)
    c->cd(2); gPad->SetGridx(); gPad->SetGridy(); gPad->SetRightMargin(0.13); gPad->SetLogz();
    h2->SetTitle("Per-event HG vs LG, cap NE-D (all energies): HG clips, LG continues;LG peak (energy) [mV];HG peak (timing) [mV]");
    h2->Draw("COLZ");
    TProfile* pf=h2->ProfileX("pf"); pf->SetLineColor(kBlack); pf->SetLineWidth(3); pf->SetMarkerStyle(20); pf->SetMarkerSize(0.8); pf->Draw("SAME");
    TLine* sl=new TLine(0,cfg.hg_sat_mV,1500,cfg.hg_sat_mV); sl->SetLineColor(kRed); sl->SetLineWidth(2); sl->SetLineStyle(2); sl->Draw();
    TLatex ts; ts.SetTextColor(kRed); ts.SetTextSize(0.030); ts.DrawLatex(120,cfg.hg_sat_mV+22,Form("DT5742 HG-channel clip ~%.0f mV",cfg.hg_sat_mV));

    gSystem->mkdir("figures",kTRUE);
    c->Print("figures/hg_lg_saturation.png");
    printf("\nnormalized @150: HG=%.2fx, LG=%.2fx (ideal linear=%.1fx)  -> HG saturates, LG linear\n", nhg[5], nlg[5], ideal[5]);
    printf("wrote figures/hg_lg_saturation.png\n");
}

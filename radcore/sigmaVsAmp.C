// ============================================================================
// sigmaVsAmp.C — is DSB1 timing slew-limited (1/amplitude) + a floor, with the
// energies collapsing onto ONE sigma_t(amplitude) curve, and NO DT5742 saturation?
// ----------------------------------------------------------------------------
// For each beam energy: same fiducial/quality selection as the headline method
// (ScanRunCenters centroid, r<TimingFiducialR, MCP[200,750], per-channel
// hg_peak>=20). Per event: t=(DW-UP)/2 of the timing caps, and A = mean HG peak
// amplitude of the channels used. Then sigma_t in bins of A, overlaid across
// energies (collapse test) + the per-event HG-amplitude distributions with the
// DRS4 saturation rail marked.
//   ROOT_INCLUDE_PATH=radcore:Analysis root -l -b -q radcore/sigmaVsAmp.C+
// ============================================================================
#include "PlotUtils.h"       // FitGaussCore, ApplyRADiCALStyle
#include "SelectionCuts.h"   // kMCP1_minPeak.., kHG_minPeak, kMCP_minPeak_E, kSumLG_centroid, TimingFiducialR
#include "DataPaths.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>

// VecToHist_teb + FitGaussCore (same as the headline sigma estimator)
static double tebSigma(std::vector<float>& v) {
    if (v.size() < 300) return -1;
    double mu1=0; for(float x:v) mu1+=x; mu1/=v.size();
    double ms1=0; for(float x:v) ms1+=(x-mu1)*(x-mu1); ms1=std::sqrt(ms1/v.size()); if(ms1<0.008) ms1=0.1;
    double mu2=0; int n2=0; for(float x:v) if(std::fabs(x-mu1)<5*ms1){mu2+=x;++n2;}
    double ms2=ms1; if(n2>0){ mu2/=n2; ms2=0; for(float x:v) if(std::fabs(x-mu1)<5*ms1) ms2+=(x-mu2)*(x-mu2); ms2=std::sqrt(ms2/n2); if(ms2<0.008) ms2=0.1; } else mu2=mu1;
    TH1F h("_sva","",120,mu2-4*ms2,mu2+4*ms2); h.SetDirectory(nullptr); for(float x:v) h.Fill(x);
    double mu,muE,s,sE; FitGaussCore(&h,2.0,mu,muE,s,sE); return s>0? s*1000.0:-1;
}

void sigmaVsAmp() {
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    double Es[6]={25,50,75,100,125,150};
    int    col[6]={kAzure+2,kCyan+2,kGreen+2,kOrange+7,kRed+1,kMagenta+2};
    const int NA=24; double Alo=0, Ahi=1200, dA=(Ahi-Alo)/NA;   // HG amplitude bins [mV]
    const double kSat=950.0;                                    // hg_saturated threshold (DT5742 rail)

    std::vector<double> allA, allS, allSe;                      // pooled sigma_t(A) points for the fit
    TGraphErrors* g[6]; TH1F* hA[6];
    double meanA[6]={0}; double satf[6]={0};                    // per-energy mean amplitude + saturated fraction
    printf("\n%-6s %10s %10s %10s\n","E","meanA[mV]","sat-frac","Nfid");

    for (int e=0;e<6;++e) {
        TFile* fp=TFile::Open(radReduced("DSB1",Es[e]));
        if(!fp||fp->IsZombie()){ g[e]=nullptr; hA[e]=nullptr; continue; }
        TTree* t=(TTree*)fp->Get("rad");
        Bool_t wc; Float_t x,y,mcp,slg,hgc[8],hgp[8]; Bool_t hsat[8];
        t->SetBranchAddress("wc_ok",&wc); t->SetBranchAddress("x_trk",&x); t->SetBranchAddress("y_trk",&y);
        t->SetBranchAddress("mcp_peak",&mcp); t->SetBranchAddress("sum_lg",&slg);
        t->SetBranchAddress("hg_cfd05",hgc); t->SetBranchAddress("hg_peak",hgp); t->SetBranchAddress("hg_saturated",hsat);
        long N=t->GetEntries();
        // pass 1: LG-weighted centroid (ScanRunCenters)
        double wx=0,wy=0,w=0; for(long i=0;i<N;++i){ t->GetEntry(i); if(wc&&mcp>kMCP_minPeak_E&&slg>kSumLG_centroid){wx+=x*slg;wy+=y*slg;w+=slg;} }
        double xc=w>0?wx/w:0, yc=w>0?wy/w:0; double rFid=TimingFiducialR(Es[e]), r2=rFid*rFid;
        // pass 2: per-A-bin (DW-UP)/2 + amplitude distribution + saturation count
        std::vector<float> Tbin[NA]; hA[e]=new TH1F(Form("hA%d",e),"",60,0,1200); hA[e]->SetDirectory(nullptr);
        long nfid=0, nsat=0; double sumA=0;
        for(long i=0;i<N;++i){ t->GetEntry(i);
            if(!wc||mcp<kMCP1_minPeak||mcp>kMCP1_maxPeak) continue;
            double dx=x-xc,dy=y-yc; if(dx*dx+dy*dy>=r2) continue;
            double ds=0,us=0,asum=0; int dn=0,un=0,na=0; bool anySat=false;
            for(int c=0;c<4;++c) if(hgp[c]>=kHG_minPeak && hgc[c]>-1e5){ ds+=hgc[c];++dn; asum+=hgp[c];++na; if(hsat[c])anySat=true; }
            for(int c=4;c<8;++c) if(hgp[c]>=kHG_minPeak && hgc[c]>-1e5){ us+=hgc[c];++un; asum+=hgp[c];++na; if(hsat[c])anySat=true; }
            if(dn<1||un<1) continue;
            float tt=0.5f*(float)(ds/dn-us/un), A=asum/na;
            hA[e]->Fill(A); sumA+=A; ++nfid; if(anySat)++nsat;
            int b=(int)((A-Alo)/dA); if(b>=0&&b<NA) Tbin[b].push_back(tt);
        }
        // sigma_t per A-bin
        std::vector<double> gx,gy,gxe,gye;
        for(int b=0;b<NA;++b){ if(Tbin[b].size()<500) continue; double s=tebSigma(Tbin[b]); if(s>10){
            gx.push_back(Alo+(b+0.5)*dA); gy.push_back(s); gxe.push_back(0); gye.push_back(s/std::sqrt(2.0*Tbin[b].size()));
            allA.push_back(Alo+(b+0.5)*dA); allS.push_back(s); allSe.push_back(s/std::sqrt(2.0*Tbin[b].size())); } }
        g[e]=new TGraphErrors(gx.size(), gx.data(), gy.data(), gxe.data(), gye.data());
        g[e]->SetMarkerColor(col[e]); g[e]->SetLineColor(col[e]); g[e]->SetMarkerStyle(20+e%4); g[e]->SetMarkerSize(1.4);
        if(hA[e]->Integral()>0) hA[e]->Scale(1.0/hA[e]->Integral());
        hA[e]->SetLineColor(col[e]); hA[e]->SetLineWidth(2);
        meanA[e]=nfid?sumA/nfid:0; satf[e]=nfid?(double)nsat/nfid:0;
        printf("%-6.0f %10.0f %10.3f %10ld\n", Es[e], meanA[e], satf[e], nfid);
        fp->Close();
    }

    // ---- fit the pooled sigma_t(A) = sqrt((k/A)^2 + c^2) ----
    TGraphErrors* gAll=new TGraphErrors(allA.size(), allA.data(), allS.data(), nullptr, allSe.data());
    TF1* f=new TF1("f","sqrt([0]*[0]/(x*x)+[1]*[1])", 200, 1000); f->SetParameters(20000,35); f->SetParLimits(1,0,60);
    gAll->Fit(f,"RQN");
    double k=f->GetParameter(0), c=f->GetParameter(1);

    TCanvas* cv=new TCanvas("c_sva","",1860,560); cv->Divide(3,1,0.010,0.006);

    // Panel 1: mean HG amplitude vs beam energy — the SOFT CEILING
    cv->cd(1); gPad->SetGridx(); gPad->SetGridy();
    TH1F* fr1=gPad->DrawFrame(0,0,160,1050);
    fr1->SetTitle("HG timing amplitude vs energy: soft ceiling (no hard saturation);beam energy [GeV];mean HG peak amplitude  A  [mV]");
    TGraph* gMA=new TGraph(6,Es,meanA); gMA->SetMarkerStyle(20); gMA->SetMarkerSize(1.6); gMA->SetMarkerColor(kAzure+2); gMA->SetLineColor(kAzure+2); gMA->SetLineWidth(2);
    gMA->Draw("PL SAME");
    TLine* rail=new TLine(0,kSat,160,kSat); rail->SetLineColor(kRed); rail->SetLineWidth(2); rail->Draw();
    double ceil=0; for(int e=1;e<6;++e) ceil+=meanA[e]; ceil/=5;
    TLine* cl=new TLine(0,ceil,160,ceil); cl->SetLineColor(kGray+2); cl->SetLineStyle(2); cl->SetLineWidth(2); cl->Draw();
    TLatex t1; t1.SetTextSize(0.030);
    t1.SetTextColor(kRed);    t1.DrawLatex(8,kSat+12,"DT5742 hard rail ~950 mV  (sat-frac = 0)");
    t1.SetTextColor(kGray+3); t1.DrawLatex(8,ceil-45,Form("soft ceiling ~%.0f mV (SiPM/amp compression)",ceil));
    t1.SetTextColor(kAzure+3);t1.DrawLatex(30,300,"amplitude plateaus by 50 GeV");
    t1.DrawLatex(30,250,"=> slew rate fixed => #sigma_{t} flat 50-150");

    // Panel 2: HG amplitude distributions per energy + DRS4 rail
    cv->cd(2); gPad->SetGridx(); gPad->SetGridy();
    double ymx=0; for(int e=0;e<6;++e) if(hA[e]) ymx=std::max(ymx,hA[e]->GetMaximum());
    TH1F* fr2=gPad->DrawFrame(0,0,1050,ymx*1.30);
    fr2->SetTitle("Per-event HG amplitude by energy (spikes pile below the rail);mean HG peak amplitude  A  [mV];fraction of events");
    TLegend* lg2=new TLegend(0.16,0.58,0.42,0.90); lg2->SetBorderSize(0);
    for(int e=0;e<6;++e) if(hA[e]){ hA[e]->Draw("HIST SAME"); lg2->AddEntry(hA[e],Form("%.0f GeV",Es[e]),"l"); }
    TLine* sat=new TLine(kSat,0,kSat,ymx*1.30); sat->SetLineColor(kRed); sat->SetLineWidth(2); sat->Draw();
    TLatex ts; ts.SetTextColor(kRed); ts.SetTextSize(0.030); ts.SetTextAngle(90); ts.DrawLatex(kSat-22,ymx*0.35,"DT5742 rail ~950 mV");
    lg2->Draw();

    // Panel 3: sigma_t vs amplitude — the collapse + floor
    cv->cd(3); gPad->SetGridx(); gPad->SetGridy();
    TH1F* fr3=gPad->DrawFrame(200,0,1000,90);
    fr3->SetTitle("#sigma_{t} vs HG amplitude: one curve, all energies;mean HG peak amplitude  A  [mV];#sigma_{t}  (DW#minusUP)/2  [ps]");
    f->SetLineColor(kBlack); f->SetLineWidth(2); f->SetLineStyle(2); f->Draw("SAME");
    TLine* fl=new TLine(200,c,1000,c); fl->SetLineColor(kGray+2); fl->SetLineStyle(3); fl->SetLineWidth(2); fl->Draw();
    TLegend* lg=new TLegend(0.16,0.15,0.45,0.48); lg->SetBorderSize(0);
    for(int e=0;e<6;++e) if(g[e]&&g[e]->GetN()>0){ g[e]->Draw("P SAME"); lg->AddEntry(g[e],Form("%.0f GeV",Es[e]),"p"); }
    lg->AddEntry(f,Form("#sqrt{(k/A)^{2}+c^{2}}",k),"l"); lg->AddEntry(fl,Form("floor c = %.1f ps",c),"l"); lg->Draw();
    TLatex tl; tl.SetNDC(); tl.SetTextSize(0.030);
    tl.DrawLatex(0.50,0.82,"25 GeV: low A (on the ramp)");
    tl.DrawLatex(0.50,0.77,"50-150 GeV: pinned at the ceiling A");
    tl.DrawLatex(0.50,0.72,"=> energy enters ONLY via amplitude");

    gSystem->mkdir("radcore/figs",kTRUE);
    cv->Print("radcore/figs/sigma_vs_amplitude.png");
    printf("\nmeanA[mV] by E: "); for(int e=0;e<6;++e) printf("%.0f ",meanA[e]); printf(" (rail 950, sat-frac=0)\n");
    printf("POOLED FIT  sigma_t(A) = sqrt((k/A)^2 + c^2):  k = %.0f mV*ps,  all-event FLOOR c = %.1f ps\n", k, c);
    printf("wrote radcore/figs/sigma_vs_amplitude.png\n");
}

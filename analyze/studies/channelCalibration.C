// ============================================================================
// channelCalibration.C — inter-channel relative gain calibration  (Gap #4)
// ============================================================================
//
// Ledovskoy (Nov–Dec 2021) equalised the capillary channels using amplitude
// correlations (A_j = C * A_i) and the MIP peak.  Our electron-only data has no
// clean MIP peak (verified) so absolute deposited-energy calibration via MIPs is
// not available — but the INTER-CHANNEL RELATIVE calibration is, and is what
// makes the channel-combination and charge-sharing observables well-behaved.
//
// METHOD
//   - Use the LG (energy) amplitudes of EM-shower events.
//   - Relative gain g_i = <A_i> / <A_ref> over a clean shower sample, with the
//     reference taken as the geometric-mean channel (so g's straddle 1).
//   - Equalised amplitude A_i' = A_i / g_i.
//   - VALIDATION (Ledovskoy's "equal within 10%"): spread of the 8 equalised
//     channel amplitudes within an event, before vs after equalisation.
//   - Report g_i vs energy (stability) and the per-energy uniformity.
//
// Output:
//   output/Summary/channel_calibration.pdf  (3 pages)
//   output/Summary/channel_calibration.root (g_i graphs, constants)
//
// Usage:  root -l -b -q 'Analysis/channelCalibration.C+'
// ============================================================================

#include "ChannelConfig.h"
#include "PlotUtils.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TSystem.h"
#include "TParameter.h"

#include <cmath>
#include <iostream>
#include <vector>

void channelCalibration()
{
    TH1::AddDirectory(kFALSE);
    ApplyRADiCALStyle();
    gSystem->mkdir("output/Summary", kTRUE);

    // Per-energy: mean LG amplitude per channel (shower events) and the
    // per-event channel-spread before/after equalisation.
    std::vector<double> vE;
    std::vector<double> vGain[8];          // relative gain g_i vs energy
    std::vector<double> vSpreadBefore, vSpreadAfter;  // RMS/mean of 8 channels per event [%]

    for (int iRun = 0; iRun < kNRuns; ++iRun) {
        const RunCfg& rc = kRuns[iRun];
        TFile* f = TFile::Open(TString("output/") + rc.label + "/ntuple.root");
        if (!f || f->IsZombie()) { std::cout << "  skip " << rc.label << "\n"; continue; }
        TTree* t = (TTree*)f->Get("rad");
        if (!t || t->GetEntries() == 0) { f->Close(); continue; }

        double xc, yc, tcfd, trms; ScanRunCenters(t, xc, yc, tcfd, trms);

        Bool_t  wc_ok; Float_t x_trk, y_trk, mcp_peak, sum_lg, sum_pb, lg_peak[8];
        t->SetBranchAddress("wc_ok",   &wc_ok);
        t->SetBranchAddress("x_trk",   &x_trk);
        t->SetBranchAddress("y_trk",   &y_trk);
        t->SetBranchAddress(t->GetBranch("mcp1_peak")?"mcp1_peak":"mcp_peak",&mcp_peak);
        t->SetBranchAddress("sum_lg",  &sum_lg);
        t->SetBranchAddress("sum_pb",  &sum_pb);
        t->SetBranchAddress("lg_peak",  lg_peak);

        double sumA[8] = {}; long nA[8] = {};
        // collect events for the spread computation (after we know the gains we
        // need a second loop; cache the per-event amplitudes)
        std::vector<std::vector<float>> cache;  cache.reserve(20000);

        Long64_t N = t->GetEntries();
        for (Long64_t e = 0; e < N; ++e) {
            t->GetEntry(e);
            if (!wc_ok || mcp_peak < kMCP1_minPeak || mcp_peak > kMCP1_maxPeak) continue;
            float dx = x_trk - (float)xc, dy = y_trk - (float)yc;
            if (std::sqrt(dx*dx + dy*dy) >= (float)TimingFiducialR(rc.energy_GeV)) continue;
            if (sum_lg <= kSumLG_centroid) continue;                       // need a real shower
            if (sum_pb >= kPb_maxRatio * sum_lg) continue;                 // contained

            bool allOk = true;
            for (int i = 0; i < 8; ++i) if (lg_peak[i] < kLG_minPeak) { allOk = false; break; }
            if (!allOk) continue;   // require all 8 channels to have signal (uniform sample)

            std::vector<float> a(8);
            for (int i = 0; i < 8; ++i) { a[i] = lg_peak[i]; sumA[i] += lg_peak[i]; ++nA[i]; }
            cache.push_back(a);
        }
        t->ResetBranchAddresses(); f->Close();

        if (cache.size() < 200) { std::cout << "  " << rc.label << ": too few\n"; continue; }

        // Relative gains: g_i = <A_i> / geomean(<A>)
        double mean[8], logsum = 0.;
        for (int i = 0; i < 8; ++i) { mean[i] = (nA[i] > 0) ? sumA[i]/nA[i] : 0.; logsum += std::log(std::max(mean[i],1e-6)); }
        double geo = std::exp(logsum / 8.);
        double g[8]; for (int i = 0; i < 8; ++i) g[i] = (geo > 0.) ? mean[i]/geo : 1.;

        // Per-event channel spread (RMS/mean of the 8 channels), before vs after
        double sB = 0., sA = 0.; long nS = 0;
        for (auto& a : cache) {
            double mb = 0., ma = 0.;
            for (int i = 0; i < 8; ++i) { mb += a[i]; ma += a[i]/g[i]; }
            mb /= 8.; ma /= 8.;
            double vb = 0., va = 0.;
            for (int i = 0; i < 8; ++i) { vb += (a[i]-mb)*(a[i]-mb); va += (a[i]/g[i]-ma)*(a[i]/g[i]-ma); }
            if (mb > 0.) sB += std::sqrt(vb/8.)/mb;
            if (ma > 0.) sA += std::sqrt(va/8.)/ma;
            ++nS;
        }
        vE.push_back(rc.energy_GeV);
        for (int i = 0; i < 8; ++i) vGain[i].push_back(g[i]);
        vSpreadBefore.push_back(100.*sB/nS);
        vSpreadAfter .push_back(100.*sA/nS);

        std::cout << Form("  %-7s  channel spread (RMS/mean): %.1f%% -> %.1f%% after equalisation\n",
                          rc.label.Data(), 100.*sB/nS, 100.*sA/nS);
        if (iRun == kNRuns - 1) {
            std::cout << "  150 GeV relative gains g_i: ";
            for (int i = 0; i < 8; ++i) std::cout << Form("%s=%.2f ", kCap[i].name, g[i]);
            std::cout << "\n";
        }
    }

    if (vE.empty()) { std::cout << "channelCalibration: no data.\n"; return; }

    // =========================================================================
    // Output
    // =========================================================================
    TString outPDF = "output/Summary/channel_calibration.pdf";
    TString outROOT= "output/Summary/channel_calibration.root";
    TCanvas c("c_cc", "", 1000, 720);

    // -- Page 1: relative gain g_i vs energy (stability) -----------------------
    c.Clear(); c.cd(); StylePad(false, true);
    TH1F* fr = (TH1F*)c.DrawFrame(15., 0.4, 165., 1.8, ";Beam energy (GeV);relative gain g_{i}");
    fr->GetXaxis()->SetTitleSize(0.048); fr->GetYaxis()->SetTitleSize(0.048);
    fr->GetYaxis()->SetTitleOffset(1.30);
    TLine* l1 = new TLine(15., 1., 165., 1.); l1->SetLineStyle(2); l1->SetLineColor(kGray+1); l1->Draw();
    TGraph* gG[8]; TLegend* L = MakeLegend(8);
    for (int i = 0; i < 8; ++i) {
        gG[i] = new TGraph();
        for (size_t k = 0; k < vE.size(); ++k) gG[i]->SetPoint(gG[i]->GetN(), vE[k], vGain[i][k]);
        gG[i]->SetLineColor(kRChannelCols[i]); gG[i]->SetMarkerColor(kRChannelCols[i]);
        gG[i]->SetMarkerStyle(20+i); gG[i]->SetMarkerSize(1.0); gG[i]->SetLineWidth(2);
        gG[i]->Draw("PL SAME");
        L->AddEntry(gG[i], kCap[i].name, "lp");
    }
    L->Draw();
    DrawPageTitle("Inter-channel relative gain vs energy (LG; stable g #Rightarrow good calibration)");
    c.Print(outPDF + "(");

    // -- Page 2: per-event channel uniformity before/after equalisation --------
    c.Clear(); c.cd(); StylePad();
    double yMax = 0.; for (double s : vSpreadBefore) yMax = std::max(yMax, s);
    TH1F* fr2 = (TH1F*)c.DrawFrame(15., 0., 165., 1.30*yMax,
                                   ";Beam energy (GeV);8-channel spread RMS/mean (%)");
    fr2->GetXaxis()->SetTitleSize(0.048); fr2->GetYaxis()->SetTitleSize(0.048);
    fr2->GetYaxis()->SetTitleOffset(1.30);
    auto mk=[&](std::vector<double>& y,int col,int mk2){ TGraph* g=new TGraph();
        for(size_t k=0;k<vE.size();++k) g->SetPoint(g->GetN(),vE[k],y[k]);
        g->SetLineColor(col);g->SetMarkerColor(col);g->SetMarkerStyle(mk2);g->SetMarkerSize(1.4);g->SetLineWidth(2);return g;};
    TGraph* gB = mk(vSpreadBefore, kGray+2, 20);
    TGraph* gA = mk(vSpreadAfter,  kRGreen,  21);
    gB->Draw("PL SAME"); gA->Draw("PL SAME");
    TLegend leg(0.55,0.74,0.93,0.88); leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.038);
    leg.AddEntry(gB, "raw", "lp"); leg.AddEntry(gA, "after equalisation", "lp");
    leg.Draw();
    {
        TLatex a; a.SetNDC(); a.SetTextSize(0.034); a.SetTextColor(kGray+2);
        a.DrawLatex(0.16, 0.30, "Spread of the 8 channel amplitudes within an event.");
        a.DrawLatex(0.16, 0.25, "Equalisation removes the fixed per-channel gain differences;");
        a.DrawLatex(0.16, 0.20, "the residual is genuine shower-position charge sharing.");
    }
    DrawPageTitle("Per-event 8-channel uniformity: raw vs equalised");
    c.Print(outPDF + ")");

    // -- ROOT output -----------------------------------------------------------
    TFile* fo = new TFile(outROOT, "RECREATE");
    for (int i = 0; i < 8; ++i) gG[i]->Write(Form("gGain_%s", kCap[i].name));
    gB->Write("gSpreadRaw"); gA->Write("gSpreadEqualised");
    // store the 150 GeV (last) gains as parameters
    for (int i = 0; i < 8; ++i)
        if (!vGain[i].empty())
            TParameter<double>(Form("gain_%s", kCap[i].name), vGain[i].back()).Write();
    fo->Close(); delete fo;

    std::cout << "channelCalibration: done -> " << outPDF << "\n";
}

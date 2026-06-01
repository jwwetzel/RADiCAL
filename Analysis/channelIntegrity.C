// ============================================================================
// channelIntegrity.C -- Layer 1 Hardware Integrity Report for RADiCAL
// ============================================================================
//
// Produces a 5-page PDF assessing inter-channel consistency and signal
// integrity across all six beam energies.  Focus is on hardware health
// rather than physics performance.
//
//   Page 1  Channel active fraction vs beam energy (8 TGraphs)
//   Page 2  Mean HG/LG amplitude ratio per channel vs energy (8 TGraphs)
//   Page 3  Per-channel timing spread -- sigma_t vs energy (from summary.root)
//   Page 4  HG pedestal RMS distributions, all energies overlaid
//   Page 5  Cross-channel amplitude correlation matrix at 150 GeV (heat map)
//
// Prerequisites:
//   processRun.C   -- must have produced Analysis/Output/<label>/ntuple.root
//   analyzeResolution.C -- must have produced Analysis/Output/Summary/summary.root
//                          (Page 3 gracefully skipped if summary.root is absent)
//
// Output:
//   Analysis/Output/Summary/channel_integrity.pdf
//
// Usage (from repository root):
//   root -l -b -q 'Analysis/channelIntegrity.C+'
// ============================================================================

#include "ChannelConfig.h"
#include "PlotUtils.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TBox.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TAxis.h"
#include "TPad.h"

#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <numeric>

// ---------------------------------------------------------------------------
// Marker styles: one per energy run (matches kEMark used in compareEnergies.C)
// ---------------------------------------------------------------------------
static const int kEMark[6] = { 20, 21, 22, 23, 24, 25 };

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------
void channelIntegrity()
{
    TH1::AddDirectory(kFALSE);
    ApplyRADiCALStyle();

    gSystem->mkdir("Analysis/Output/Summary", kTRUE);

    const TString pdfPath = "Analysis/Output/Summary/channel_integrity.pdf";

    // =========================================================================
    // Accumulators filled during the per-energy loop
    //   [iRun][iCh]
    // =========================================================================
    double energies[kNRuns]  = {};        // beam energy in GeV
    bool   runDone[kNRuns]   = {};        // true if ntuple was read

    // Page 1: active fraction
    double activeFrac[kNRuns][8]   = {};

    // Page 2: mean HG/LG ratio
    double hglgRatio[kNRuns][8]    = {};

    // Page 4: per-energy pedestal RMS histograms (booked after loop)
    // We store them as TH1F* for later drawing.
    TH1F* hPedRMS[kNRuns] = {};

    // Page 5: correlation matrix data at 150 GeV (iRun == kNRuns-1)
    // Accumulators: mean[i], mean[i]*mean[j], var[i][j] (Pearson)
    double corrSum[8]    = {};   // sum of hg_peak[i]
    double corrSum2[8][8]= {};   // sum of hg_peak[i]*hg_peak[j]
    long   corrN         = 0;

    // =========================================================================
    // Per-energy loop
    // =========================================================================
    for (int iRun = 0; iRun < kNRuns; ++iRun)
    {
        const RunCfg& rc = kRuns[iRun];
        TString ntFile   = TString("Analysis/Output/") + rc.label + "/ntuple.root";

        TFile* fin = TFile::Open(ntFile);
        if (!fin || fin->IsZombie()) {
            std::cout << "[channelIntegrity] Skipping " << rc.label
                      << " -- ntuple not found.\n";
            continue;
        }
        TTree* tree = (TTree*)fin->Get("rad");
        if (!tree || tree->GetEntries() == 0) { fin->Close(); continue; }

        std::cout << "[channelIntegrity] " << rc.label
                  << " -- " << tree->GetEntries() << " events\n";

        energies[iRun] = rc.energy_GeV;

        // Derive beam centroid for fiducial radius
        double xc, yc, tOff, tRms;
        ScanRunCenters(tree, xc, yc, tOff, tRms);
        const float xcf = static_cast<float>(xc);
        const float ycf = static_cast<float>(yc);

        // ── Branch addresses ──────────────────────────────────────────────────
        Bool_t  wc_ok;
        Float_t x_trk, y_trk, mcp_peak, mcp2_peak;
        Float_t hg_peak[8], lg_peak[8], hg_ped_rms_v[8];

        tree->SetBranchAddress("wc_ok",    &wc_ok);
        tree->SetBranchAddress("x_trk",    &x_trk);
        tree->SetBranchAddress("y_trk",    &y_trk);
        tree->SetBranchAddress("mcp_peak", &mcp_peak);
        tree->SetBranchAddress("hg_peak",  hg_peak);
        tree->SetBranchAddress("lg_peak",  lg_peak);

        bool hasMCP2   = (tree->GetBranch("mcp2_peak") != nullptr);
        bool hasPedRMS = (tree->GetBranch("hg_ped_rms") != nullptr);
        mcp2_peak = 0.f;
        for (int i = 0; i < 8; ++i) hg_ped_rms_v[i] = 0.f;

        if (hasMCP2)   tree->SetBranchAddress("mcp2_peak",  &mcp2_peak);
        if (hasPedRMS) tree->SetBranchAddress("hg_ped_rms",  hg_ped_rms_v);

        // ── Book the pedestal RMS histogram for this energy ───────────────────
        TH1F* hPed = new TH1F(Form("hPedRMS_%s", rc.label.Data()),
            ";HG ped RMS (mV);Normalised", 100, 0., 10.);
        hPedRMS[iRun] = hPed;

        // ── Accumulators for active fraction and HG/LG ratio ──────────────────
        long nDenom       = 0;             // events passing base quality cuts
        long nActive[8]   = {};            // active channel count
        double sumHG[8]   = {};            // sum of HG peak for active events
        double sumLG[8]   = {};            // sum of LG peak for active events
        long   nHGLG[8]   = {};            // count of dual-active events

        const bool is150GeV = (iRun == kNRuns - 1);

        Long64_t nEntries = tree->GetEntries();
        for (Long64_t ev = 0; ev < nEntries; ++ev)
        {
            tree->GetEntry(ev);

            // Base quality: wc_ok + MCP quality + timing fiducial
            if (!wc_ok) continue;

            // MCP quality -- use the appropriate MCP per channel convention.
            // For the denominator count we use MCP1 quality (majority of channels).
            // Channel 7 (SW-U) uses MCP2 but its active fraction denominator
            // is still defined by MCP1 quality for inter-channel comparability.
            bool mcpOk = (mcp_peak >= kMCP1_minPeak && mcp_peak <= kMCP1_maxPeak);
            if (!mcpOk) continue;

            float dx = x_trk - xcf, dy = y_trk - ycf;
            float r  = std::sqrt(dx*dx + dy*dy);
            if (r >= static_cast<float>(kFiducial_r_timing)) continue;

            ++nDenom;

            // ── Pedestal RMS (fill for all qualifying events) ─────────────────
            if (hasPedRMS) {
                for (int i = 0; i < 8; ++i)
                    hPed->Fill(hg_ped_rms_v[i]);
            }

            // ── Per-channel quantities ─────────────────────────────────────────
            for (int i = 0; i < 8; ++i)
            {
                bool active = (hg_peak[i] > kHG_minPeak);
                if (active) {
                    ++nActive[i];
                    // HG/LG ratio accumulation
                    if (lg_peak[i] > kLG_minPeak) {
                        sumHG[i] += hg_peak[i];
                        sumLG[i] += lg_peak[i];
                        ++nHGLG[i];
                    }
                }
            }

            // ── Page 5: correlation matrix at 150 GeV ─────────────────────────
            if (is150GeV) {
                // Check that at least one channel fired (very loose)
                bool anyActive = false;
                for (int i = 0; i < 8; ++i)
                    if (hg_peak[i] > kHG_minPeak) { anyActive = true; break; }
                if (!anyActive) continue;

                ++corrN;
                for (int i = 0; i < 8; ++i)
                    corrSum[i] += hg_peak[i];
                for (int i = 0; i < 8; ++i)
                    for (int j = 0; j < 8; ++j)
                        corrSum2[i][j] += hg_peak[i] * hg_peak[j];
            }

        } // event loop

        fin->Close();   // branch addresses no longer needed; tree is gone

        runDone[iRun] = true;

        // Store results
        if (nDenom > 0) {
            for (int i = 0; i < 8; ++i) {
                activeFrac[iRun][i] = 100.0 * nActive[i] / nDenom;
                if (nHGLG[i] > 0)
                    hglgRatio[iRun][i] = sumHG[i] / sumLG[i];
            }
        }

        std::cout << "  Denominator events: " << nDenom << "\n";
        for (int i = 0; i < 8; ++i)
            std::cout << "    " << kCap[i].name
                      << "  active=" << Form("%.1f%%", activeFrac[iRun][i])
                      << "  HG/LG=" << Form("%.2f", hglgRatio[iRun][i]) << "\n";

    } // run loop

    // =========================================================================
    // PAGE 1 -- Channel active fraction vs beam energy
    // =========================================================================
    {
        TCanvas c("cCI1", "Channel Active Fraction", 1200, 700);
        c.cd();
        StylePad(false, true);   // sidebar for legend

        // Build TMultiGraph: one TGraph per channel
        TMultiGraph* mg = new TMultiGraph();
        TLegend* leg = MakeLegend(8);

        TGraph* gAct[8];
        for (int i = 0; i < 8; ++i)
            gAct[i] = new TGraph();

        for (int iRun = 0; iRun < kNRuns; ++iRun) {
            if (!runDone[iRun]) continue;
            for (int i = 0; i < 8; ++i) {
                int n = gAct[i]->GetN();
                gAct[i]->SetPoint(n, energies[iRun], activeFrac[iRun][i]);
            }
        }

        for (int i = 0; i < 8; ++i) {
            gAct[i]->SetMarkerStyle(kEMark[i % 6]);
            gAct[i]->SetMarkerSize(1.2);
            gAct[i]->SetMarkerColor(kRChannelCols[i]);
            gAct[i]->SetLineColor(kRChannelCols[i]);
            gAct[i]->SetLineWidth(2);
            mg->Add(gAct[i], "PL");
            leg->AddEntry(gAct[i], kCap[i].name, "lp");
        }

        mg->Draw("A");
        mg->GetXaxis()->SetTitle("Beam energy (GeV)");
        mg->GetYaxis()->SetTitle("Active fraction (%)");
        mg->GetXaxis()->SetRangeUser(10., 165.);
        // Reasonable y range: channels should be 70-100%
        mg->GetYaxis()->SetRangeUser(0., 110.);

        DrawPadTitle("Channel Active Fraction vs Beam Energy");
        leg->Draw();

        // Annotation: NW-U (ch 4) expected lowest
        {
            TLatex ann; ann.SetNDC(); ann.SetTextSize(0.036);
            ann.SetTextColor(kGray+2);
            ann.DrawLatex(0.14, 0.14, "NW-U (ch 4) expected lowest");
        }

        c.cd(0); DrawPageTitle("RADiCAL Channel Integrity -- Page 1: Active Fraction");
        c.Print(pdfPath + "(");
        delete mg;
    }

    // =========================================================================
    // PAGE 2 -- Mean HG/LG amplitude ratio per channel vs energy
    // =========================================================================
    {
        TCanvas c("cCI2", "HG/LG Ratio", 1200, 700);
        c.cd();
        StylePad(false, true);

        TMultiGraph* mg = new TMultiGraph();
        TLegend* leg = MakeLegend(8);

        TGraph* gRat[8];
        for (int i = 0; i < 8; ++i)
            gRat[i] = new TGraph();

        for (int iRun = 0; iRun < kNRuns; ++iRun) {
            if (!runDone[iRun]) continue;
            for (int i = 0; i < 8; ++i) {
                if (hglgRatio[iRun][i] <= 0.) continue;
                int n = gRat[i]->GetN();
                gRat[i]->SetPoint(n, energies[iRun], hglgRatio[iRun][i]);
            }
        }

        double yMax = 0.;
        for (int i = 0; i < 8; ++i) {
            gRat[i]->SetMarkerStyle(kEMark[i % 6]);
            gRat[i]->SetMarkerSize(1.2);
            gRat[i]->SetMarkerColor(kRChannelCols[i]);
            gRat[i]->SetLineColor(kRChannelCols[i]);
            gRat[i]->SetLineWidth(2);
            if (gRat[i]->GetN() > 0) mg->Add(gRat[i], "PL");
            leg->AddEntry(gRat[i], kCap[i].name, "lp");
            for (int p = 0; p < gRat[i]->GetN(); ++p)
                yMax = std::max(yMax, gRat[i]->GetY()[p]);
        }

        mg->Draw("A");
        mg->GetXaxis()->SetTitle("Beam energy (GeV)");
        mg->GetYaxis()->SetTitle("Mean HG/LG ratio");
        mg->GetXaxis()->SetRangeUser(10., 165.);
        mg->GetYaxis()->SetRangeUser(0., yMax * 1.35);

        DrawPadTitle("Mean HG/LG Amplitude Ratio vs Beam Energy");
        leg->Draw();

        {
            TLatex ann; ann.SetNDC(); ann.SetTextSize(0.034);
            ann.SetTextColor(kGray+2);
            ann.DrawLatex(0.14, 0.17, "Stable ratio: both channels healthy");
            ann.DrawLatex(0.14, 0.11, "Rising ratio at high energy: HG saturation");
        }

        c.cd(0); DrawPageTitle("RADiCAL Channel Integrity -- Page 2: HG/LG Ratio");
        c.Print(pdfPath);
        delete mg;
    }

    // =========================================================================
    // PAGE 3 -- Per-channel timing spread from summary.root
    // =========================================================================
    {
        TCanvas c("cCI3", "Timing Spread", 1200, 700);
        c.cd();
        StylePad(false, true);

        const TString sumFile = "Analysis/Output/Summary/summary.root";
        TFile* fSum = TFile::Open(sumFile);

        if (!fSum || fSum->IsZombie()) {
            // Graceful fallback: print a notice page
            TLatex msg; msg.SetNDC(); msg.SetTextSize(0.052);
            msg.SetTextAlign(22);
            msg.DrawLatex(0.50, 0.58,
                "Page 3 requires Analysis/Output/Summary/summary.root");
            msg.SetTextSize(0.040); msg.SetTextColor(kGray+2);
            msg.DrawLatex(0.50, 0.48, "Run analyzeResolution.C first, then re-run this macro.");
            c.cd(0); DrawPageTitle("RADiCAL Channel Integrity -- Page 3: Timing Spread (MISSING DATA)");
            c.Print(pdfPath);
            std::cout << "[channelIntegrity] WARNING: summary.root not found -- "
                      << "Page 3 shows placeholder.\n";
        } else {
            // Read per-channel timing graphs
            TGraphErrors* gTResCap[8];
            TGraphErrors* gTResAvg = (TGraphErrors*)fSum->Get("gTimingResolution");
            bool anyGraph = false;
            double yMax = 0.;

            for (int i = 0; i < 8; ++i) {
                gTResCap[i] = (TGraphErrors*)fSum->Get(
                    Form("gTRes_%s", kCap[i].name));
                if (gTResCap[i]) {
                    // Restyle: use channel colors
                    gTResCap[i]->SetMarkerStyle(20 + i);
                    gTResCap[i]->SetMarkerSize(1.1);
                    gTResCap[i]->SetMarkerColor(kRChannelCols[i]);
                    gTResCap[i]->SetLineColor(kRChannelCols[i]);
                    gTResCap[i]->SetLineWidth(2);
                    for (int p = 0; p < gTResCap[i]->GetN(); ++p)
                        yMax = std::max(yMax, gTResCap[i]->GetY()[p]);
                    anyGraph = true;
                }
            }

            if (gTResAvg) {
                gTResAvg->SetMarkerStyle(29);  // full star
                gTResAvg->SetMarkerSize(1.8);
                gTResAvg->SetMarkerColor(kBlack);
                gTResAvg->SetLineColor(kBlack);
                gTResAvg->SetLineWidth(3);
                for (int p = 0; p < gTResAvg->GetN(); ++p)
                    yMax = std::max(yMax, gTResAvg->GetY()[p]);
            }

            if (!anyGraph) {
                TLatex msg; msg.SetNDC(); msg.SetTextSize(0.048);
                msg.SetTextAlign(22);
                msg.DrawLatex(0.50, 0.55,
                    "No gTRes_* graphs found in summary.root");
            } else {
                // Reference band: +-30% around 150 GeV best (published 27 ps)
                const double tBest150 = 27.;  // ps, published result
                const double refLo    = tBest150 * 0.70;
                const double refHi    = tBest150 * 1.30;

                if (yMax < refHi * 1.5) yMax = refHi * 2.0;
                yMax *= 1.20;

                // Draw the reference band first (behind graphs)
                TBox refBand(10., refLo, 165., refHi);
                refBand.SetFillColor(kGray);
                refBand.SetFillStyle(3005);
                refBand.SetLineColor(kGray+1);
                refBand.SetLineWidth(1);

                // Need a frame before Box (Box needs axes to exist)
                TH1F* fr = gPad->DrawFrame(10., 0., 165., yMax,
                    ";Beam energy (GeV);#sigma_{t} (ps)");
                (void)fr;

                refBand.Draw("SAME");

                for (int i = 0; i < 8; ++i) {
                    if (gTResCap[i] && gTResCap[i]->GetN() > 0)
                        gTResCap[i]->Draw("P SAME");
                }
                if (gTResAvg && gTResAvg->GetN() > 0)
                    gTResAvg->Draw("P SAME");

                // Legend
                TLegend* leg = MakeLegend(10);
                for (int i = 0; i < 8; ++i)
                    if (gTResCap[i] && gTResCap[i]->GetN() > 0)
                        leg->AddEntry(gTResCap[i], kCap[i].name, "lp");
                if (gTResAvg && gTResAvg->GetN() > 0)
                    leg->AddEntry(gTResAvg, "8-ch avg", "lp");
                leg->Draw();

                // Band annotation
                {
                    TLatex ann; ann.SetNDC(); ann.SetTextSize(0.030);
                    ann.SetTextColor(kGray+2);
                    // The band is the COMBINED (DW-UP)/2 headline (~27 ps).  Single
                    // channels (~250 ps) sit far above it BY DESIGN: each is
                    // MCP-reference-limited; combining the 8 channels reaches the band.
                    ann.DrawLatex(0.14, 0.90,
                        Form("Gray band: combined (DW#minusUP)/2 headline 27 ps #pm30%% "
                             "(%.0f-%.0f ps)", refLo, refHi));
                    ann.DrawLatex(0.14, 0.85,
                        "Single-channel #sigma_{t} (~250 ps) is MCP-limited; "
                        "the 8-channel combination reaches the band.");
                }
            }

            DrawPadTitle("Per-channel Timing Resolution vs Beam Energy");
            c.cd(0); DrawPageTitle("RADiCAL Channel Integrity -- Page 3: Timing Spread");
            c.Print(pdfPath);

            fSum->Close();
        }
    }

    // =========================================================================
    // PAGE 4 -- HG pedestal RMS distributions (all energies overlaid)
    // =========================================================================
    {
        TCanvas c("cCI4", "Pedestal RMS", 1200, 700);
        c.cd();
        StylePad(false, true);

        // Normalise to unit area and overlay
        TLegend* leg = MakeLegend(7);   // 6 energies + threshold line

        bool first = true;
        double yMax = 0.;

        for (int iRun = 0; iRun < kNRuns; ++iRun) {
            if (!runDone[iRun] || !hPedRMS[iRun]) continue;
            if (hPedRMS[iRun]->GetEntries() < 10) continue;

            double integral = hPedRMS[iRun]->Integral();
            if (integral > 0.) hPedRMS[iRun]->Scale(1.0 / integral);

            hPedRMS[iRun]->SetLineColor(kREnergyCols[iRun]);
            hPedRMS[iRun]->SetLineWidth(2);

            yMax = std::max(yMax, hPedRMS[iRun]->GetMaximum());

            hPedRMS[iRun]->Draw(first ? "HIST" : "HIST SAME");
            first = false;
            leg->AddEntry(hPedRMS[iRun],
                Form("%.0f GeV", energies[iRun]), "l");
        }

        if (!first) {
            // Re-set y-axis range on the first histogram (drawn first)
            for (int iRun = 0; iRun < kNRuns; ++iRun) {
                if (!runDone[iRun] || !hPedRMS[iRun]) continue;
                if (hPedRMS[iRun]->GetEntries() < 10) continue;
                hPedRMS[iRun]->GetYaxis()->SetRangeUser(0., yMax * 1.30);
                break;
            }

            // Threshold line at 5 mV
            const double thresh = 5.0;
            TLine* lThr = new TLine(thresh, 0., thresh, yMax * 1.25);
            lThr->SetLineColor(kRed+1);
            lThr->SetLineStyle(2);
            lThr->SetLineWidth(2);
            lThr->Draw("SAME");
            leg->AddEntry(lThr, "5 mV threshold", "l");

            leg->Draw();

            {
                TLatex ann; ann.SetNDC(); ann.SetTextSize(0.034);
                ann.SetTextColor(kGray+2);
                ann.DrawLatex(0.14, 0.87,
                    "Narrow dist. = uniform DRS4 noise");
                ann.DrawLatex(0.14, 0.81,
                    "Broad dist. = noisy cells or pickup");
            }
        } else {
            TLatex msg; msg.SetNDC(); msg.SetTextSize(0.048);
            msg.SetTextAlign(22);
            msg.DrawLatex(0.50, 0.55, "No hg_ped_rms branch found in ntuples.");
            msg.SetTextSize(0.038); msg.SetTextColor(kGray+2);
            msg.DrawLatex(0.50, 0.47, "Re-run processRun.C to add this branch.");
        }

        DrawPadTitle("HG Pedestal RMS -- All Energies (normalised)");
        c.cd(0); DrawPageTitle("RADiCAL Channel Integrity -- Page 4: Pedestal RMS");
        c.Print(pdfPath);
    }

    // =========================================================================
    // PAGE 5 -- Cross-channel amplitude correlation matrix at 150 GeV
    // =========================================================================
    {
        TCanvas c("cCI5", "Correlation Matrix", 900, 850);
        c.cd();
        StylePad(true);   // COLZ -- wide right margin for palette, no grid

        // Compute Pearson correlation coefficients from the accumulators
        // rho(i,j) = (E[XiXj] - E[Xi]*E[Xj]) / (sigma_i * sigma_j)

        // Correlation matrix TH2F -- 8x8
        TH2F* h2Corr = new TH2F("hCorr150",
            ";Channel;Channel;Pearson #rho",
            8, -0.5, 7.5, 8, -0.5, 7.5);

        bool haveCorr = (corrN > 10);
        if (haveCorr) {
            // Mean per channel
            double mean[8] = {};
            for (int i = 0; i < 8; ++i)
                mean[i] = corrSum[i] / corrN;

            // Variance per channel (for normalisation)
            double var[8] = {};
            for (int i = 0; i < 8; ++i) {
                double e2 = corrSum2[i][i] / corrN;
                var[i]    = e2 - mean[i] * mean[i];
                if (var[i] < 0.) var[i] = 0.;
            }

            for (int i = 0; i < 8; ++i) {
                for (int j = 0; j < 8; ++j) {
                    double eij  = corrSum2[i][j] / corrN;
                    double cov  = eij - mean[i] * mean[j];
                    double denom = std::sqrt(var[i] * var[j]);
                    double rho  = (denom > 0.) ? cov / denom : 0.;
                    // Clamp to [-1, 1] for numerical safety
                    if (rho >  1.) rho =  1.;
                    if (rho < -1.) rho = -1.;
                    // TH2F: bin (ix, iy), ROOT bins start at 1
                    h2Corr->SetBinContent(i+1, j+1, rho);
                }
            }
        } else {
            // No 150 GeV data -- fill diagonal with 1s for a placeholder
            for (int i = 0; i < 8; ++i)
                h2Corr->SetBinContent(i+1, i+1, 1.0);
            std::cout << "[channelIntegrity] WARNING: 150 GeV ntuple not found or "
                      << "too few events for correlation matrix.\n";
        }

        // Axis labels: channel names
        for (int i = 0; i < 8; ++i) {
            h2Corr->GetXaxis()->SetBinLabel(i+1, kCap[i].name);
            h2Corr->GetYaxis()->SetBinLabel(i+1, kCap[i].name);
        }
        h2Corr->GetXaxis()->SetLabelSize(0.052);
        h2Corr->GetYaxis()->SetLabelSize(0.052);
        // Clear the redundant "Channel" x-axis title (cells are already labelled
        // NW-D, NE-D, ...) so it doesn't collide with the bottom annotation.
        h2Corr->GetXaxis()->SetTitle("");

        // Zoom the colour scale to the populated band.  All correlations are
        // 0.88-1.00 (channels see the same shower), so a 0-1 scale rendered
        // every cell the same bright colour and hid the structure.  Min 0.85
        // clips nothing (observed minimum is 0.88).
        h2Corr->SetMinimum(0.85);
        h2Corr->SetMaximum(1.00);

        // Palette (inverted kCherry) already set by ApplyRADiCALStyle
        h2Corr->Draw("COLZ");

        // Print rho values as text in each cell
        if (haveCorr) {
            TLatex cell; cell.SetNDC(kFALSE);
            cell.SetTextSize(0.028);
            cell.SetTextAlign(22);
            for (int i = 0; i < 8; ++i) {
                for (int j = 0; j < 8; ++j) {
                    double rho = h2Corr->GetBinContent(i+1, j+1);
                    // Bin center in user coordinates
                    double bx = h2Corr->GetXaxis()->GetBinCenter(i+1);
                    double by = h2Corr->GetYaxis()->GetBinCenter(j+1);
                    // On the zoomed 0.85-1.0 scale: bright (high-rho) cells get
                    // dark text, darker (low-rho) cells get white text.
                    cell.SetTextColor(rho > 0.93 ? kBlack : kWhite);
                    cell.DrawLatex(bx, by, Form("%.2f", rho));
                }
            }
        }

        DrawPadTitle("Amplitude Correlation Matrix -- 150 GeV");
        {
            TLatex ann; ann.SetNDC(); ann.SetTextSize(0.026); ann.SetTextAlign(22);
            ann.SetTextColor(kGray+2);
            // The high off-diagonal rho (~0.9) is COMMON-MODE shower sharing (all
            // channels see the same shower), NOT electronic cross-talk.
            // Centred + shorter so it no longer clips the right frame edge.
            ann.DrawLatex(0.52, 0.03,
                "High #rho: common-mode shower sharing, not electronic cross-talk");
        }
        c.cd(0); DrawPageTitle("RADiCAL Channel Integrity -- Page 5: Correlation Matrix (150 GeV)");
        c.Print(pdfPath + ")");
    }

    // ── Persist per-channel vitals for the Layer-1 hero plot ────────────────
    //   gActiveFrac_<chan>, gHGLG_<chan> : value vs beam energy.  layer1Summary.C
    //   reads these for the "channel vitals" figure.
    {
        TFile fout("Analysis/Output/Summary/channel_integrity.root", "RECREATE");
        for (int i = 0; i < 8; ++i) {
            TGraph gAF, gRR;
            for (int iRun = 0; iRun < kNRuns; ++iRun) {
                if (activeFrac[iRun][i] > 0.)
                    gAF.SetPoint(gAF.GetN(), energies[iRun], activeFrac[iRun][i]);
                if (hglgRatio[iRun][i] > 0.)
                    gRR.SetPoint(gRR.GetN(), energies[iRun], hglgRatio[iRun][i]);
            }
            gAF.Write(Form("gActiveFrac_%s", kCap[i].name));
            gRR.Write(Form("gHGLG_%s",       kCap[i].name));
        }
        fout.Close();
        std::cout << "[channelIntegrity] Wrote "
                  << "Analysis/Output/Summary/channel_integrity.root\n";
    }

    std::cout << "\n[channelIntegrity] Done.\n"
              << "  " << pdfPath << "\n";
}

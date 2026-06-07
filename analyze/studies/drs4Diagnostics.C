// ============================================================================
// drs4Diagnostics.C -- DRS4/DT5742 hardware integrity diagnostics
// ============================================================================
//
// Layer 1 of the analysis hierarchy: verify DRS4 time calibration quality
// and flag per-channel hardware issues.  Reads the processed ntuples
// (not raw waveform data) and diagnoses DRS4/DT5742 behavior from the
// branches already present in the ntuple.
//
// Pages 3 and 4 require the hg_saturated[8] and hg_spike[8] branches added
// by processRun.C upgrade.  If those branches are absent, those pages are
// gracefully skipped with a printed instruction.
//
// ── Output ───────────────────────────────────────────────────────────────────
//
//   output/Summary/drs4_diagnostics.pdf  (6 pages)
//   output/Summary/drs4_diagnostics.root
//
//   Page 1: HG pedestal RMS (noise floor) per channel vs beam energy
//             8 channels overlaid; reference line at 5 mV noise threshold
//
//   Page 2: HG peak amplitude distributions -- 8 panels (4x2)
//             One histogram per channel; 6 energies overlaid; normalised
//
//   Page 3: HG saturation fraction per channel vs beam energy
//             Uses hg_saturated[8] branch; warns if any channel > 5%
//
//   Page 4: HG spike rate per channel vs beam energy
//             Uses hg_spike[8] branch; note on spike definition
//
//   Page 5: Mean pedestal RMS vs energy -- DRS4 temperature drift diagnostic
//             Mean + envelope (max/min band) over all 8 channels
//
//   Page 6: HG/LG amplitude correlation at 150 GeV -- cross-talk check
//             2x4 grid of TH2F scatter plots (COLZ, log-z)
//
// Usage (from repository root):
//   root -l -b -q 'Analysis/drs4Diagnostics.C+'
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
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TBox.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TAxis.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

// ---------------------------------------------------------------------------
// OpenNtuple
// ---------------------------------------------------------------------------
static TFile* OpenNtuple(int r, TTree*& tree)
{
    TString path = Form("output/%s/ntuple.root",
                        kRuns[r].label.Data());
    TFile* f = TFile::Open(path);
    if (!f || f->IsZombie()) {
        std::cerr << "drs4Diagnostics: cannot open " << path << "\n";
        tree = nullptr; return nullptr;
    }
    tree = static_cast<TTree*>(f->Get("rad"));
    if (!tree) {
        std::cerr << "drs4Diagnostics: tree 'rad' missing in " << path << "\n";
        f->Close(); delete f; tree = nullptr; return nullptr;
    }
    return f;
}

// ---------------------------------------------------------------------------
// DrawNoDataNote -- used for graceful skip of pages 3 and 4
// ---------------------------------------------------------------------------
static void DrawNoDataNote(const char* msg)
{
    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.042);
    lat.SetTextAlign(22);
    lat.SetTextColor(kGray + 1);
    lat.DrawLatex(0.50, 0.55, msg);
    lat.SetTextSize(0.034);
    lat.DrawLatex(0.50, 0.45,
        "Run processRun.C to add hg_saturated / hg_spike branches");
}

// ===========================================================================
// Main
// ===========================================================================
void drs4Diagnostics()
{
    TH1::AddDirectory(kFALSE);
    ApplyRADiCALStyle();

    gSystem->mkdir("output/Summary", kTRUE);

    std::cout << "drs4Diagnostics: building DRS4 hardware diagnostics\n";

    // =========================================================================
    // Step 1 -- accumulators
    // =========================================================================

    // Per-run, per-channel mean pedestal RMS
    double pedRms[kNRuns][8]   = {};
    int    pedRmsN[kNRuns][8]  = {};

    // Saturation and spike fractions (require new branches)
    double satFrac[kNRuns][8]  = {};
    double spkFrac[kNRuns][8]  = {};
    bool   hasSatBranch        = false;
    bool   hasSpkBranch        = false;

    // Peak amplitude histograms: [channel][run], 6 energies overlaid per channel
    TH1F* hPeak[8][kNRuns] = {};
    for (int i = 0; i < 8; ++i)
        for (int r = 0; r < kNRuns; ++r) {
            hPeak[i][r] = new TH1F(
                Form("hPeak_ch%d_r%d", i, r), "",
                100, 0., 1000.);
            hPeak[i][r]->SetDirectory(nullptr);
        }

    // HG vs LG scatter at 150 GeV (run index 5)
    TH2F* hCorr[8] = {};
    for (int i = 0; i < 8; ++i) {
        hCorr[i] = new TH2F(
            Form("hCorr_ch%d", i), "",
            80, 0., 1000.,   // LG x-axis
            80, 0., 1000.);  // HG y-axis
        hCorr[i]->SetDirectory(nullptr);
    }

    // =========================================================================
    // Step 2 -- event loop over all runs
    // =========================================================================
    for (int r = 0; r < kNRuns; ++r) {
        TTree* t = nullptr;
        TFile* fin = OpenNtuple(r, t);
        if (!fin) continue;

        // Check for new branches (test once on first valid file)
        if (r == 0) {
            hasSatBranch = (t->GetBranch("hg_saturated") != nullptr);
            hasSpkBranch = (t->GetBranch("hg_spike")     != nullptr);
            if (!hasSatBranch)
                std::cout << "  [INFO] hg_saturated branch absent -- "
                             "pages 3/4 will be skipped\n";
        }

        Float_t hg_peak[8], hg_ped_rms[8], lg_peak[8];
        Bool_t  hg_saturated[8], hg_spike[8];
        Bool_t  wc_ok;
        Float_t mcp_peak;

        t->SetBranchAddress("wc_ok",       &wc_ok);
        t->SetBranchAddress("mcp_peak",    &mcp_peak);
        t->SetBranchAddress("hg_peak",      hg_peak);
        t->SetBranchAddress("hg_ped_rms",   hg_ped_rms);
        t->SetBranchAddress("lg_peak",      lg_peak);

        if (hasSatBranch) t->SetBranchAddress("hg_saturated", hg_saturated);
        if (hasSpkBranch) t->SetBranchAddress("hg_spike",     hg_spike);

        // Initialize per-run fraction counters
        long long nTotal[8]  = {};
        long long nSat[8]    = {};
        long long nSpk[8]    = {};

        Long64_t nEv = t->GetEntries();
        for (Long64_t ev = 0; ev < nEv; ++ev) {
            t->GetEntry(ev);
            if (!wc_ok) continue;
            if (mcp_peak < kMCP_minPeak_E) continue;

            for (int i = 0; i < 8; ++i) {
                // Pedestal RMS -- accumulate for all valid channels
                if (hg_ped_rms[i] > 0.f) {
                    pedRms[r][i]  += hg_ped_rms[i];
                    pedRmsN[r][i] += 1;
                }

                // Peak amplitude histogram
                hPeak[i][r]->Fill(hg_peak[i]);

                // Saturation and spike fractions
                ++nTotal[i];
                if (hasSatBranch && hg_saturated[i]) ++nSat[i];
                if (hasSpkBranch && hg_spike[i])     ++nSpk[i];

                // HG/LG correlation: 150 GeV only (last run index kNRuns-1)
                if (r == kNRuns - 1) {
                    hCorr[i]->Fill(lg_peak[i], hg_peak[i]);
                }
            }
        }

        t->ResetBranchAddresses();
        fin->Close();
        delete fin;

        // Finalise means and fractions for this run
        for (int i = 0; i < 8; ++i) {
            if (pedRmsN[r][i] > 0)
                pedRms[r][i] /= pedRmsN[r][i];
            if (nTotal[i] > 0) {
                satFrac[r][i] = 100.0 * static_cast<double>(nSat[i])
                                      / static_cast<double>(nTotal[i]);
                spkFrac[r][i] = 100.0 * static_cast<double>(nSpk[i])
                                      / static_cast<double>(nTotal[i]);
            }
        }

        std::cout << "  " << kRuns[r].label
                  << Form(": %.0f GeV -- %lld events processed\n",
                          kRuns[r].energy_GeV, nEv);
    }

    // =========================================================================
    // Step 3 -- Build TGraphs for per-channel vs energy plots
    // =========================================================================

    // Page 1: pedestal RMS per channel
    TGraph* gPedRms[8]  = {};
    // Page 3: saturation fraction per channel
    TGraph* gSatFrac[8] = {};
    // Page 4: spike fraction per channel
    TGraph* gSpkFrac[8] = {};
    // Page 5: mean pedestal RMS over channels (+ envelope)
    TGraph* gPedMean    = new TGraph();
    TGraph* gPedMax     = new TGraph();
    TGraph* gPedMin     = new TGraph();

    for (int i = 0; i < 8; ++i) {
        gPedRms[i]  = new TGraph();
        gSatFrac[i] = new TGraph();
        gSpkFrac[i] = new TGraph();
    }

    for (int r = 0; r < kNRuns; ++r) {
        double E = kRuns[r].energy_GeV;

        // Per-channel points
        for (int i = 0; i < 8; ++i) {
            if (pedRmsN[r][i] > 0) {
                int np = gPedRms[i]->GetN();
                gPedRms[i]->SetPoint(np, E, pedRms[r][i]);
            }
            {
                int np = gSatFrac[i]->GetN();
                gSatFrac[i]->SetPoint(np, E, satFrac[r][i]);
            }
            {
                int np = gSpkFrac[i]->GetN();
                gSpkFrac[i]->SetPoint(np, E, spkFrac[r][i]);
            }
        }

        // Mean/max/min over channels for page 5
        double sumR = 0., maxR = -1e9, minR = 1e9;
        int    nValid = 0;
        for (int i = 0; i < 8; ++i) {
            if (pedRmsN[r][i] > 0) {
                sumR += pedRms[r][i];
                if (pedRms[r][i] > maxR) maxR = pedRms[r][i];
                if (pedRms[r][i] < minR) minR = pedRms[r][i];
                ++nValid;
            }
        }
        if (nValid > 0) {
            int np = gPedMean->GetN();
            gPedMean->SetPoint(np, E, sumR / nValid);
            gPedMax ->SetPoint(np, E, maxR);
            gPedMin ->SetPoint(np, E, minR);
        }
    }

    // Style helper lambda: apply per-channel color to a TGraph
    auto styleGraph = [](TGraph* g, int col) {
        g->SetLineColor(col);
        g->SetMarkerColor(col);
        g->SetMarkerStyle(20);
        g->SetMarkerSize(1.2);
        g->SetLineWidth(2);
    };

    for (int i = 0; i < 8; ++i) {
        styleGraph(gPedRms[i],  kRChannelCols[i]);
        styleGraph(gSatFrac[i], kRChannelCols[i]);
        styleGraph(gSpkFrac[i], kRChannelCols[i]);
    }

    // =========================================================================
    // Step 4 -- Open output PDF and ROOT file
    // =========================================================================
    TString outPDF  = "output/Summary/drs4_diagnostics.pdf";
    TString outROOT = "output/Summary/drs4_diagnostics.root";

    // Use a 4:3 canvas; multi-panel pages use Divide()
    TCanvas c("c_drs4", "", 960, 720);

    // =========================================================================
    // Page 1: HG pedestal RMS per channel vs beam energy
    // =========================================================================
    c.Clear(); c.cd();
    StylePad(false, true);   // sidebar for legend

    // Determine y-axis range
    double yMaxP1 = 0.;
    for (int i = 0; i < 8; ++i)
        for (int r = 0; r < kNRuns; ++r)
            if (pedRmsN[r][i] > 0)
                yMaxP1 = std::max(yMaxP1, pedRms[r][i]);
    yMaxP1 = std::max(yMaxP1 * 1.30, 7.0);

    // Draw frame first
    TH1F* frame1 = static_cast<TH1F*>(
        c.DrawFrame(15., 0., 165., yMaxP1,
                    ";Beam energy (GeV);Pedestal RMS (mV)"));
    frame1->GetXaxis()->SetTitleSize(0.048);
    frame1->GetYaxis()->SetTitleSize(0.048);
    frame1->GetYaxis()->SetTitleOffset(1.35);

    // Reference line at 5 mV
    TLine* lNoise = new TLine(15., 5., 165., 5.);
    lNoise->SetLineColor(kGray + 1);
    lNoise->SetLineStyle(2);
    lNoise->SetLineWidth(2);
    lNoise->Draw("SAME");

    {
        TLatex ann; ann.SetNDC(); ann.SetTextSize(0.036);
        ann.SetTextColor(kGray + 1); ann.SetTextAlign(12);
        // Place label near the right end of the reference line (NDC)
        // 5 mV is at data y=5; frame goes 0 to yMaxP1
        ann.DrawLatex(0.74, 0.12 + 5.0 / yMaxP1 * 0.78,
                      "noise threshold");
    }

    TLegend* leg1 = MakeLegend(8);
    for (int i = 0; i < 8; ++i) {
        if (gPedRms[i]->GetN() > 0) {
            gPedRms[i]->Draw("PL SAME");
            leg1->AddEntry(gPedRms[i], kCap[i].name, "lp");
        }
    }
    leg1->Draw();

    DrawPadTitle("HG pedestal RMS per channel vs beam energy");
    c.Print(outPDF + "(");

    // =========================================================================
    // Page 2: HG peak amplitude distributions -- 4x2 canvas
    // =========================================================================
    c.Clear();
    c.Divide(4, 2, 0.002f, 0.002f);

    // Normalise histograms to unit area
    for (int i = 0; i < 8; ++i)
        for (int r = 0; r < kNRuns; ++r) {
            double integral = hPeak[i][r]->Integral();
            if (integral > 0.) hPeak[i][r]->Scale(1.0 / integral);
        }

    for (int i = 0; i < 8; ++i) {
        c.cd(i + 1);
        StylePad();

        double yMaxP2 = 0.;
        for (int r = 0; r < kNRuns; ++r)
            yMaxP2 = std::max(yMaxP2, hPeak[i][r]->GetMaximum());
        yMaxP2 *= 1.25;
        if (yMaxP2 <= 0.) yMaxP2 = 0.10;

        bool firstDraw = true;
        for (int r = 0; r < kNRuns; ++r) {
            hPeak[i][r]->SetLineColor(kREnergyCols[r]);
            hPeak[i][r]->SetLineWidth(2);
            hPeak[i][r]->GetXaxis()->SetTitle("HG peak amplitude (mV)");
            hPeak[i][r]->GetYaxis()->SetTitle("Events (normalised)");
            hPeak[i][r]->GetXaxis()->SetRangeUser(0., 1000.);
            hPeak[i][r]->GetYaxis()->SetRangeUser(0., yMaxP2);
            if (firstDraw) {
                hPeak[i][r]->Draw("HIST");
                firstDraw = false;
            } else {
                hPeak[i][r]->Draw("HIST SAME");
            }
        }

        // Lower cut line
        TLine* lLow = new TLine(kHG_minPeak, 0., kHG_minPeak, yMaxP2);
        lLow->SetLineColor(kGray + 1);
        lLow->SetLineStyle(2);
        lLow->SetLineWidth(1);
        lLow->Draw("SAME");

        // Saturation cut line
        TLine* lSat = new TLine(kHG_maxPeak, 0., kHG_maxPeak, yMaxP2);
        lSat->SetLineColor(kRed + 1);
        lSat->SetLineStyle(2);
        lSat->SetLineWidth(1);
        lSat->Draw("SAME");

        // Labels for cut lines (only on first panel to avoid clutter)
        if (i == 0) {
            TLatex lab; lab.SetNDC(); lab.SetTextSize(0.065);
            lab.SetTextColor(kGray + 1); lab.SetTextAlign(12);
            // kHG_minPeak is at x=20 out of [0,1000]; frame width (NDC) ~0.83
            lab.DrawLatex(0.17, 0.80, "kHG_minPeak");
            lab.SetTextColor(kRed + 1);
            lab.DrawLatex(0.87, 0.80, "sat");
        }

        DrawPadTitle(kCap[i].name);
    }

    // Energy legend on the canvas frame (cd(0))
    c.cd(0);
    {
        TLegend legE(0.01f, 0.00f, 0.99f, 0.035f);
        legE.SetNColumns(kNRuns);
        legE.SetBorderSize(0);
        legE.SetFillStyle(0);
        legE.SetTextSize(0.022f);
        for (int r = 0; r < kNRuns; ++r) {
            // dummy line entry
            TLine* dummy = new TLine();
            dummy->SetLineColor(kREnergyCols[r]);
            dummy->SetLineWidth(2);
            legE.AddEntry(dummy,
                          Form("%.0f GeV", kRuns[r].energy_GeV), "l");
        }
        legE.Draw();
        DrawPageTitle("HG peak amplitude distributions -- all channels (normalised)");
    }
    c.Print(outPDF);

    // =========================================================================
    // Page 3: HG saturation fraction per channel vs beam energy
    // =========================================================================
    c.Clear(); c.cd();
    StylePad(false, true);

    if (!hasSatBranch) {
        DrawNoDataNote("hg_saturated branch not found in ntuples");
        DrawPadTitle("HG saturation fraction -- data unavailable");
    } else {
        double yMaxP3 = 0.;
        for (int i = 0; i < 8; ++i)
            for (int r = 0; r < kNRuns; ++r)
                yMaxP3 = std::max(yMaxP3, satFrac[r][i]);
        yMaxP3 = std::max(yMaxP3 * 1.35, 2.0);

        TH1F* frame3 = static_cast<TH1F*>(
            c.DrawFrame(15., 0., 165., yMaxP3,
                        ";Beam energy (GeV);Fraction saturated (%)"));
        frame3->GetXaxis()->SetTitleSize(0.048);
        frame3->GetYaxis()->SetTitleSize(0.048);
        frame3->GetYaxis()->SetTitleOffset(1.35);
        (void)frame3;

        // Warning threshold line at 5%
        if (yMaxP3 >= 5.) {
            TLine* lWarn = new TLine(15., 5., 165., 5.);
            lWarn->SetLineColor(kOrange + 2);
            lWarn->SetLineStyle(2);
            lWarn->SetLineWidth(2);
            lWarn->Draw("SAME");
            TLatex ann; ann.SetNDC(); ann.SetTextSize(0.036);
            ann.SetTextColor(kOrange + 2); ann.SetTextAlign(12);
            ann.DrawLatex(0.74, 0.12 + 5.0 / yMaxP3 * 0.78, "5% warning");
        }

        TLegend* leg3 = MakeLegend(8);
        bool anyWarning = false;
        for (int i = 0; i < 8; ++i) {
            if (gSatFrac[i]->GetN() > 0) {
                gSatFrac[i]->Draw("PL SAME");
                leg3->AddEntry(gSatFrac[i], kCap[i].name, "lp");
            }
            for (int r = 0; r < kNRuns; ++r)
                if (satFrac[r][i] > 5.) anyWarning = true;
        }
        leg3->Draw();

        if (anyWarning) {
            TLatex warn; warn.SetNDC(); warn.SetTextSize(0.038);
            warn.SetTextColor(kRed + 1); warn.SetTextAlign(12);
            warn.DrawLatex(0.14, 0.88,
                "WARNING: one or more channels exceed 5% saturation -- "
                "consider adjusting gain");
        } else {
            // All channels read ~0 — make the empty plot informative.
            TLatex ok; ok.SetNDC(); ok.SetTextSize(0.044); ok.SetTextAlign(22);
            ok.SetTextColor(kRGreen + 1);
            ok.DrawLatex(0.52, 0.55,
                "All channels < 0.1% saturated at every energy");
            ok.SetTextSize(0.034); ok.SetTextColor(kGray + 2);
            ok.DrawLatex(0.52, 0.48,
                "(no DRS4 saturation; HG stays below the 950 mV rail)");
        }

        // Smaller, shorter title so it does not clip the frame edges
        DrawPadTitle("HG saturation fraction vs beam energy", 0.050);
    }
    c.Print(outPDF);

    // =========================================================================
    // Page 4: HG spike rate per channel vs beam energy
    // =========================================================================
    c.Clear(); c.cd();
    StylePad(false, true);

    if (!hasSpkBranch) {
        DrawNoDataNote("hg_spike branch not found in ntuples");
        DrawPadTitle("HG spike rate -- data unavailable");
    } else {
        double yMaxP4 = 0.;
        for (int i = 0; i < 8; ++i)
            for (int r = 0; r < kNRuns; ++r)
                yMaxP4 = std::max(yMaxP4, spkFrac[r][i]);
        yMaxP4 = std::max(yMaxP4 * 1.35, 2.0);

        TH1F* frame4 = static_cast<TH1F*>(
            c.DrawFrame(15., 0., 165., yMaxP4,
                        ";Beam energy (GeV);Spike fraction (%) [hg_spike flag]"));
        frame4->GetXaxis()->SetTitleSize(0.048);
        frame4->GetYaxis()->SetTitleSize(0.048);
        frame4->GetYaxis()->SetTitleOffset(1.35);
        (void)frame4;

        TLegend* leg4 = MakeLegend(8);
        for (int i = 0; i < 8; ++i) {
            if (gSpkFrac[i]->GetN() > 0) {
                gSpkFrac[i]->Draw("PL SAME");
                leg4->AddEntry(gSpkFrac[i], kCap[i].name, "lp");
            }
        }
        leg4->Draw();

        {
            TLatex ann; ann.SetNDC(); ann.SetTextSize(0.036);
            ann.SetTextColor(kGray + 1); ann.SetTextAlign(12);
            ann.DrawLatex(0.14, 0.14,
                "spike = pedestal sample > 5#sigma_{RMS} from mean");
        }

        DrawPadTitle("HG spike rate per channel vs beam energy");
    }
    c.Print(outPDF);

    // =========================================================================
    // Page 5: Mean pedestal RMS vs energy (temperature drift diagnostic)
    // =========================================================================
    c.Clear(); c.cd();
    StylePad();

    {
        double yMaxP5 = 0., yMinP5 = 1e9;
        for (int p = 0; p < gPedMax->GetN(); ++p) {
            yMaxP5 = std::max(yMaxP5, gPedMax->GetY()[p]);
            yMinP5 = std::min(yMinP5, gPedMin->GetY()[p]);
        }
        yMaxP5 = std::max(yMaxP5 * 1.35, 7.0);
        yMinP5 = std::max(0., yMinP5 * 0.70);

        TH1F* frame5 = static_cast<TH1F*>(
            c.DrawFrame(15., yMinP5, 165., yMaxP5,
                        ";Beam energy (GeV);Mean HG pedestal RMS (mV)"));
        frame5->GetXaxis()->SetTitleSize(0.048);
        frame5->GetYaxis()->SetTitleSize(0.048);
        frame5->GetYaxis()->SetTitleOffset(1.35);
        (void)frame5;

        // Shaded envelope between min and max channels
        // Draw as filled polygon using TGraph
        int nPts = gPedMax->GetN();
        if (nPts > 0) {
            int nPoly = 2 * nPts;
            std::vector<double> xPoly(nPoly), yPoly(nPoly);
            for (int p = 0; p < nPts; ++p) {
                xPoly[p]            = gPedMax->GetX()[p];
                yPoly[p]            = gPedMax->GetY()[p];
                xPoly[nPoly-1-p]    = gPedMin->GetX()[p];
                yPoly[nPoly-1-p]    = gPedMin->GetY()[p];
            }
            TGraph* gBand = new TGraph(nPoly, xPoly.data(), yPoly.data());
            gBand->SetFillColorAlpha(kGray, 0.35f);
            gBand->SetLineColor(0);
            gBand->SetFillStyle(1001);
            gBand->Draw("F SAME");

            // Also draw the max and min edges as thin dashed lines
            gPedMax->SetLineColor(kGray + 1);
            gPedMax->SetLineStyle(2);
            gPedMax->SetLineWidth(1);
            gPedMax->SetMarkerSize(0.);
            gPedMax->Draw("L SAME");

            gPedMin->SetLineColor(kGray + 1);
            gPedMin->SetLineStyle(2);
            gPedMin->SetLineWidth(1);
            gPedMin->SetMarkerSize(0.);
            gPedMin->Draw("L SAME");
        }

        // Mean line on top
        gPedMean->SetLineColor(kRData);
        gPedMean->SetMarkerColor(kRData);
        gPedMean->SetMarkerStyle(20);
        gPedMean->SetMarkerSize(1.3);
        gPedMean->SetLineWidth(2);
        gPedMean->Draw("PL SAME");

        // Reference line at 5 mV
        TLine* lRef = new TLine(15., 5., 165., 5.);
        lRef->SetLineColor(kOrange + 2);
        lRef->SetLineStyle(2);
        lRef->SetLineWidth(2);
        lRef->Draw("SAME");

        {
            TLatex ann; ann.SetNDC(); ann.SetTextSize(0.036);
            ann.SetTextColor(kGray + 1); ann.SetTextAlign(12);
            ann.DrawLatex(0.14, 0.85,
                "Flat trend: stable DRS4 temperature");
            ann.DrawLatex(0.14, 0.79,
                "Rising/drifting trend: sampling rate or temperature shift");
            ann.SetTextColor(kGray + 2); ann.SetTextAlign(12);
            ann.DrawLatex(0.14, 0.73,
                "Shaded band: min/max envelope over all 8 channels");
        }

        TLegend legP5(0.65, 0.65, 0.93, 0.82);
        legP5.SetBorderSize(0); legP5.SetFillStyle(0); legP5.SetTextSize(0.040);
        legP5.AddEntry(gPedMean, "Mean over 8 channels", "lp");

        TGraph* gBandLeg = new TGraph();
        gBandLeg->SetFillColorAlpha(kGray, 0.35f);
        gBandLeg->SetLineColor(0);
        gBandLeg->SetFillStyle(1001);
        legP5.AddEntry(gBandLeg, "Channel envelope", "f");
        legP5.Draw();

        DrawPadTitle("Mean HG pedestal RMS vs beam energy -- DRS4 drift diagnostic");
    }
    c.Print(outPDF);

    // =========================================================================
    // Page 6: HG/LG amplitude correlation at 150 GeV (2x4 grid)
    // =========================================================================
    c.Clear();
    c.Divide(4, 2, 0.002f, 0.002f);

    for (int i = 0; i < 8; ++i) {
        c.cd(i + 1);
        StylePad(true);   // COLZ -- disable grid, wide right margin
        gPad->SetLogz(1);

        hCorr[i]->GetXaxis()->SetTitle("LG peak (mV)");
        hCorr[i]->GetYaxis()->SetTitle("HG peak (mV)");
        hCorr[i]->GetYaxis()->SetTitleOffset(1.50);
        hCorr[i]->Draw("COLZ");

        DrawPadTitle(kCap[i].name);
    }

    c.cd(0);
    DrawPageTitle(
        "HG vs LG amplitude correlation -- 150 GeV only"
        " (linear = expected; deviations indicate cross-talk)");
    c.Print(outPDF + ")");

    // =========================================================================
    // Step 5 -- Write summary ROOT file
    // =========================================================================
    TFile* fOut = new TFile(outROOT, "RECREATE");

    for (int i = 0; i < 8; ++i) {
        gPedRms[i]->Write(Form("gPedRms_ch%d", i));
        if (hasSatBranch) gSatFrac[i]->Write(Form("gSatFrac_ch%d", i));
        if (hasSpkBranch) gSpkFrac[i]->Write(Form("gSpkFrac_ch%d", i));
    }
    gPedMean->Write("gPedMean");
    gPedMax ->Write("gPedMax");
    gPedMin ->Write("gPedMin");

    for (int i = 0; i < 8; ++i)
        hCorr[i]->Write(Form("hCorr_ch%d", i));

    fOut->Close();
    delete fOut;

    // =========================================================================
    // Step 6 -- Cleanup
    // =========================================================================
    for (int i = 0; i < 8; ++i) {
        for (int r = 0; r < kNRuns; ++r) delete hPeak[i][r];
        delete hCorr[i];
        delete gPedRms[i];
        delete gSatFrac[i];
        delete gSpkFrac[i];
    }
    delete gPedMean;
    delete gPedMax;
    delete gPedMin;

    std::cout << "drs4Diagnostics: done -> " << outPDF  << "\n";
    std::cout << "                 root -> " << outROOT << "\n";
}

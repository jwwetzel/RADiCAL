// ============================================================================
// chargeProfiles.C — HG amplitude maps for all 8 capillary channels
// ============================================================================
//
// Produces TProfile2D maps of the HG (high-gain) peak amplitude for each
// capillary versus beam impact position (x_trk, y_trk).  HG channels carry
// the timing signal; their spatial response map is needed to diagnose
// per-channel efficiency differences (e.g. NW-Up underperformance).
//
// Output per energy:
//   Analysis/Output/<label>/hg_charge_profiles.pdf
//       — 4×2 grid of TProfile2D maps (one per capillary, COLZ, inverted-kRust palette)
//
// Summary outputs:
//   Analysis/Output/Summary/hg_amplitude_vs_energy.pdf
//       — 4×2 TGraph overlay: mean HG amplitude vs beam energy per capillary
//   Analysis/Output/Summary/hg_lg_ratio_maps.pdf
//       — 4×2 TProfile2D maps of hg_peak[i]/lg_peak[i] at 150 GeV
//
// Usage (from repository root):
//   root -l -b -q 'Analysis/chargeProfiles.C+'
//
// Prerequisites: run processRun.C for each energy first (see runAll.sh).
// ============================================================================

#include "ChannelConfig.h"
#include "PlotUtils.h"     // StylePad, ScanRunCenters

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TProfile2D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLatex.h"
#include "TPad.h"
#include "TMath.h"

#include <cmath>
#include <iostream>
#include <vector>

// DrawPadTitle and DrawPageTitle are provided by RADiCALStyle.h (via PlotUtils.h).

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------
void chargeProfiles()
{
    ApplyRADiCALStyle();

    gSystem->mkdir("Analysis/Output/Summary", kTRUE);

    // -------------------------------------------------------------------------
    // Storage for multi-energy summary (one mean HG amplitude per cap per run)
    // -------------------------------------------------------------------------
    std::vector<double> vE;                     // beam energies [GeV]
    std::vector<double> vMeanHG[8];             // mean HG amplitude per cap [mV]

    // Index of the 150 GeV run (for HG/LG ratio maps) — set during loop
    int iRun150 = -1;

    // =========================================================================
    // Per-energy loop
    // =========================================================================
    for (int iRun = 0; iRun < kNRuns; ++iRun)
    {
        const RunCfg& rc   = kRuns[iRun];
        TString ntupleFile = TString("Analysis/Output/") + rc.label + "/ntuple.root";

        TFile* f = TFile::Open(ntupleFile);
        if (!f || f->IsZombie()) {
            std::cout << "[chargeProfiles] Skipping " << rc.label
                      << " — ntuple not found (" << ntupleFile << ")\n"
                      << "   Run processRun.C for this energy first.\n";
            continue;
        }
        TTree* t = (TTree*)f->Get("rad");
        if (!t || t->GetEntries() == 0) {
            std::cout << "[chargeProfiles] Skipping " << rc.label
                      << " — empty tree.\n";
            f->Close();
            continue;
        }

        Long64_t nEntries = t->GetEntries();
        std::cout << "[chargeProfiles] " << rc.label
                  << " — " << nEntries << " events\n";

        // ---------------------------------------------------------------------
        // Pre-scan: derive beam centroid from signal-weighted mean position
        // ---------------------------------------------------------------------
        double x_center, y_center, t_cfd_offset, t_cfd_rms;
        ScanRunCenters(t, x_center, y_center, t_cfd_offset, t_cfd_rms);

        // ---------------------------------------------------------------------
        // Book histograms
        // ---------------------------------------------------------------------
        const double cm_hw = 6.;  // half-width of per-cap map [mm] (tight frame: ~4 mm spot)

        int nBinsCM = static_cast<int>(std::round(2.*cm_hw / kWC_resBin)); // 1 mm/bin
        TProfile2D* hHGMap[8];
        for (int i = 0; i < kNCap; ++i) {
            hHGMap[i] = new TProfile2D(
                Form("hHGMap_%s_%d", rc.label.Data(), i),
                ";x_{WC} (mm);y_{WC} (mm);A_{HG} (mV)",
                nBinsCM, x_center - cm_hw, x_center + cm_hw,
                nBinsCM, y_center - cm_hw, y_center + cm_hw);
            hHGMap[i]->SetDirectory(nullptr);  // detach from TFile — prevents
        }                                       // double-delete on f->Close()

        // HG/LG ratio maps — only needed for the highest energy run, but book
        // for every run so we can update iRun150 at the end.
        TProfile2D* hRatioMap[8];
        for (int i = 0; i < kNCap; ++i) {
            hRatioMap[i] = new TProfile2D(
                Form("hRatioMap_%s_%d", rc.label.Data(), i),
                ";x_{WC} (mm);y_{WC} (mm);A_{HG}/A_{LG}",
                nBinsCM, x_center - cm_hw, x_center + cm_hw,
                nBinsCM, y_center - cm_hw, y_center + cm_hw);
            hRatioMap[i]->SetDirectory(nullptr);
        }

        // ---------------------------------------------------------------------
        // Set branch addresses
        // ---------------------------------------------------------------------
        Bool_t  wc_ok;
        Float_t x_trk, y_trk, mcp_peak, sum_lg, sum_pb;
        Float_t hg_peak[8], lg_peak[8];

        t->SetBranchAddress("wc_ok",    &wc_ok);
        t->SetBranchAddress("x_trk",    &x_trk);
        t->SetBranchAddress("y_trk",    &y_trk);
        t->SetBranchAddress("mcp_peak", &mcp_peak);
        t->SetBranchAddress("sum_lg",   &sum_lg);
        t->SetBranchAddress("sum_pb",   &sum_pb);
        t->SetBranchAddress("hg_peak",   hg_peak);
        t->SetBranchAddress("lg_peak",   lg_peak);

        // ---------------------------------------------------------------------
        // Fill
        // ---------------------------------------------------------------------
        long nInFid = 0;
        double hgSum[8] = {};
        long   hgN  [8] = {};

        for (Long64_t ev = 0; ev < nEntries; ++ev) {
            t->GetEntry(ev);

            if (!wc_ok)                        continue;
            if (mcp_peak < kMCP_minPeak_E)     continue;

            // Dynamic fiducial cut centred on data-derived beam centroid
            float dx = x_trk - static_cast<float>(x_center);
            float dy = y_trk - static_cast<float>(y_center);
            if (std::sqrt(dx*dx + dy*dy) >= static_cast<float>(kFiducial_r_energy)) continue;

            // Shower containment cut
            if (sum_lg > kSumLG_centroid &&
                sum_pb / sum_lg > kPb_maxRatio) continue;

            ++nInFid;

            for (int i = 0; i < kNCap; ++i) {
                if (hg_peak[i] < kHG_minPeak) continue;  // below noise floor

                hHGMap[i]->Fill(x_trk, y_trk, hg_peak[i]);
                hgSum[i] += hg_peak[i];
                ++hgN[i];

                // HG/LG ratio: require valid LG signal
                if (lg_peak[i] > kLG_minPeak)
                    hRatioMap[i]->Fill(x_trk, y_trk,
                                       hg_peak[i] / lg_peak[i]);
            }
        }

        std::cout << "  Fiducial events: " << nInFid
                  << " / " << nEntries
                  << Form(" (%.1f%%)\n", 100. * nInFid / nEntries);

        for (int i = 0; i < kNCap; ++i)
            std::cout << "    " << kCap[i].name
                      << ": HG entries = " << hgN[i]
                      << "  mean = "
                      << (hgN[i] > 0 ? hgSum[i]/hgN[i] : 0.)
                      << " mV\n";

        // Store per-cap mean HG amplitude for summary graph
        vE.push_back(rc.energy_GeV);
        for (int i = 0; i < kNCap; ++i)
            vMeanHG[i].push_back(hgN[i] > 0 ? hgSum[i] / hgN[i] : 0.);

        // Track which run index is closest to 150 GeV (for ratio maps)
        if (std::fabs(rc.energy_GeV - 150.) < 1.)
            iRun150 = static_cast<int>(vE.size()) - 1;

        // Tighten z-axis range: set minimum to ~80% of the peak value so
        // intra-spot amplitude variation is visible (not washed out by noise floor)
        for (int i = 0; i < kNCap; ++i) {
            double zMax = hHGMap[i]->GetMaximum();
            if (zMax > 0.)
                hHGMap[i]->SetMinimum(0.80 * zMax);
        }

        // ---------------------------------------------------------------------
        // Per-energy PDF: 4×2 HG amplitude maps
        // ---------------------------------------------------------------------
        TString outDir = TString("Analysis/Output/") + rc.label;
        TString pdfPath = outDir + "/hg_charge_profiles.pdf";

        {
            TCanvas c("c_hg", "", 1600, 800);
            c.Divide(4, 2, 0.003, 0.035);

            for (int i = 0; i < kNCap; ++i) {
                c.cd(i + 1);
                gPad->SetRightMargin(0.18);
                StylePad(true);
                gPad->SetRightMargin(0.18);   // restore after StylePad overrides
                hHGMap[i]->Draw("COLZ");
                DrawPadTitle(Form("%.0f GeV  %s", rc.energy_GeV, kCap[i].name));
            }

            // Draw page title on the canvas (cd(0) = canvas level)
            c.cd(0);
            DrawPageTitle(Form("HG Amplitude Maps  --  %.0f GeV", rc.energy_GeV));

            c.Print(pdfPath + "(");   // open + first (only) page
            c.Print(pdfPath + ")");   // close
        }

        std::cout << "  -> " << pdfPath << "\n";

        // Keep ratio maps alive only for the 150 GeV run; delete for others
        // to save memory (they are on the heap as TProfile2D*)
        if (std::fabs(rc.energy_GeV - 150.) < 1.) {
            // Write to a temp ROOT file so we can re-draw them below.
            // (ROOT objects owned by this scope will be deleted when f->Close()
            //  is called; instead we save them to a dedicated file.)
            TString ratioFile = "Analysis/Output/Summary/_ratio_tmp.root";
            TFile* fTmp = new TFile(ratioFile, "RECREATE");
            for (int i = 0; i < kNCap; ++i)
                hRatioMap[i]->Write();
            fTmp->Close();
            delete fTmp;
        }

        f->Close();

        // Clean up TProfile2D objects (owned by the heap, not by f)
        for (int i = 0; i < kNCap; ++i) {
            delete hHGMap[i];
            delete hRatioMap[i];
        }

    } // end per-energy loop

    if (vE.empty()) {
        std::cout << "[chargeProfiles] No data processed — run processRun.C first.\n";
        return;
    }

    const int N = static_cast<int>(vE.size());

    // =========================================================================
    // Summary page: mean HG amplitude vs beam energy (one panel per channel)
    // =========================================================================
    {
        const Int_t* colors = kRChannelCols;   // Apple iOS palette from RADiCALStyle.h

        TCanvas cAmp("cAmpVsE", "", 1600, 800);
        cAmp.Divide(4, 2, 0.003, 0.035);

        for (int i = 0; i < kNCap; ++i) {
            cAmp.cd(i + 1);
            StylePad();

            TGraph* g = new TGraph(N);
            g->SetName(Form("gHGAmp_%s", kCap[i].name));
            for (int p = 0; p < N; ++p)
                g->SetPoint(p, vE[p], vMeanHG[i][p]);

            g->SetMarkerStyle(20);
            g->SetMarkerSize(1.2);
            g->SetMarkerColor(colors[i]);
            g->SetLineColor(colors[i]);
            g->SetLineWidth(2);
            g->GetXaxis()->SetTitle("Beam Energy (GeV)");
            g->GetYaxis()->SetTitle("Mean HG Amplitude (mV)");
            g->GetXaxis()->SetRangeUser(0, 180);
            g->Draw("APL");

            DrawPadTitle(kCap[i].name);
        }

        cAmp.cd(0);
        DrawPageTitle("Mean HG Amplitude vs Beam Energy  --  All Capillaries");

        TString ampPdf = "Analysis/Output/Summary/hg_amplitude_vs_energy.pdf";
        cAmp.Print(ampPdf);   // single page — "(" + ")" emitted two identical pages
        std::cout << "  -> " << ampPdf << "\n";

        // ── Persist mean-HG-amplitude-vs-energy graphs for the Layer-1 hero plot ─
        //   gHGAmp_<chan> : mean HG amplitude [mV] vs beam energy [GeV].
        TFile fAmp("Analysis/Output/Summary/hg_amplitude_vs_energy.root", "RECREATE");
        for (int i = 0; i < kNCap; ++i) {
            TGraph g(N);
            for (int p = 0; p < N; ++p) g.SetPoint(p, vE[p], vMeanHG[i][p]);
            g.Write(Form("gHGAmp_%s", kCap[i].name));
        }
        fAmp.Close();
        std::cout << "  -> Analysis/Output/Summary/hg_amplitude_vs_energy.root\n";
    }

    // =========================================================================
    // Bonus page: HG/LG ratio maps at 150 GeV
    // =========================================================================
    {
        TString ratioFile = "Analysis/Output/Summary/_ratio_tmp.root";
        TFile* fTmp = TFile::Open(ratioFile);
        if (!fTmp || fTmp->IsZombie()) {
            std::cout << "[chargeProfiles] HG/LG ratio maps skipped"
                      << " (150 GeV ntuple not available).\n";
        } else {
            TCanvas cRatio("cRatio", "", 1600, 800);
            cRatio.Divide(4, 2, 0.003, 0.035);

            bool anyMap = false;
            for (int i = 0; i < kNCap; ++i) {
                TProfile2D* hR = (TProfile2D*)fTmp->Get(
                    Form("hRatioMap_150GeV_%d", i));
                if (!hR) continue;
                anyMap = true;

                cRatio.cd(i + 1);
                gPad->SetRightMargin(0.18);
                StylePad(true);
                gPad->SetRightMargin(0.18);
                hR->Draw("COLZ");
                DrawPadTitle(Form("150 GeV  %s  HG/LG", kCap[i].name));
            }

            if (anyMap) {
                cRatio.cd(0);
                DrawPageTitle("HG/LG Amplitude Ratio Maps  --  150 GeV");

                TString ratioPdf = "Analysis/Output/Summary/hg_lg_ratio_maps.pdf";
                cRatio.Print(ratioPdf);   // single page — "(" + ")" emitted two identical pages
                std::cout << "  -> " << ratioPdf << "\n";
            }
            fTmp->Close();

            // Remove the temp file
            gSystem->Unlink(ratioFile);
        }
    }

    std::cout << "\n[chargeProfiles] Done.\n"
              << "  Per-energy PDFs : Analysis/Output/<label>/hg_charge_profiles.pdf\n"
              << "  Summary         : Analysis/Output/Summary/hg_amplitude_vs_energy.pdf\n"
              << "  Ratio maps      : Analysis/Output/Summary/hg_lg_ratio_maps.pdf\n";
}

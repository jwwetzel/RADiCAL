// ============================================================================
// uniformityScan.C -- spatial uniformity of A^2-weighted combo timing
// ============================================================================
//
// Measures the A^2-weighted combination timing resolution sigma_t as a
// function of beam impact position (x, y) within the timing fiducial.
//
// Motivation: CMS qualification requires that sigma_t does not degrade at
// crystal edges.  The NW-Up underperformance may be spatial (beam hits near
// the fiber boundary) rather than intrinsic to the channel.
//
// Timing estimator: A^2-weighted combination of all 8 HG CFD-5% channels.
//   t_combo = sum(hg_peak[i]^2 * hg_cfd[i]) / sum(hg_peak[i]^2)
//   Requires >= 2 valid channels (hg_peak[i] >= kHG_minPeak AND
//   hg_cfd[i] > -1e5f).
//
// Spatial grid: 5x5 cells of 1.2 mm x 1.2 mm spanning +-3 mm around the
// data-derived beam centroid (ScanRunCenters).  Cells with < 50 events show
// as sigma_t = 0.  Fit method: FitGaussCore with nsig=2.
//
// Event selection (applied at all energies):
//   wc_ok + MCP quality (kMCP1_minPeak..kMCP1_maxPeak) + timing fiducial
//   (r < kFiducial_r_timing) + containment cut (sum_pb < kPb_maxRatio *
//   sum_lg, only when sum_lg > kSumLG_centroid).
//
// -- Output -----------------------------------------------------------------
//   output/Summary/uniformity_scan.pdf   (3 pages)
//   output/Summary/uniformity_scan.root
//
//   Page 1: sigma_t(x,y) heat map at 150 GeV -- TProfile2D COLZ, timing
//           fiducial circle (orange) and energy fiducial circle (red dashed)
//   Page 2: sigma_t(x,y) at all 6 energies, 3x2 panel grid, shared colour scale
//   Page 3: sigma_t vs r from centroid -- 6 TGraphErrors, one per energy
//
// Usage (from repository root):
//   root -l -b -q 'Analysis/uniformityScan.C+'
// ============================================================================

#include "ChannelConfig.h"
#include "PlotUtils.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile2D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TEllipse.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TAxis.h"
#include "TColor.h"
#include "TPaletteAxis.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

// ---------------------------------------------------------------------------
// Grid parameters
// Cell size = kWC_resBin (1 mm) — matched to WC delay-line wire pitch.
// 6x6 grid spanning ±3 mm = full timing fiducial radius.
// ---------------------------------------------------------------------------
static const int    kNCell    = 6;                      // bins per axis (6×6 grid)
static const double kCellSize = kWC_resBin;             // mm per cell (= 1.0 mm)
static const double kGridHalf = kNCell * kCellSize * 0.5;  // = 3.0 mm
static const int    kMinEvt   = 50;        // minimum events per cell for fit

// Per-energy marker styles (matches the convention in the other macros)
static const int    kEMark[6] = { 20, 21, 22, 23, 24, 25 };

// Number of radial bins for Page 3
static const int    kNRBins   = 5;

// Sentinel for missing timing (must match timingContainmentScan.C)
static const float kNoTime = -1e9f;

// ---------------------------------------------------------------------------
// OpenNtuple -- open per-energy ntuple, return the "rad" tree
// ---------------------------------------------------------------------------
static TFile* OpenNtuple(int r, TTree*& tree)
{
    TString path = Form("output/%s/ntuple.root", kRuns[r].label.Data());
    TFile* f = TFile::Open(path);
    if (!f || f->IsZombie()) {
        std::cerr << "uniformityScan: cannot open " << path << "\n";
        tree = nullptr; return nullptr;
    }
    tree = static_cast<TTree*>(f->Get("rad"));
    if (!tree) {
        std::cerr << "uniformityScan: tree 'rad' missing in " << path << "\n";
        f->Close(); delete f; tree = nullptr; return nullptr;
    }
    return f;
}

// ---------------------------------------------------------------------------
// A2WeightedCombo -- A^2-weighted mean of hg_cfd[8]
//
// Returns kNoTime if fewer than 2 channels pass the amplitude + CFD threshold.
// ---------------------------------------------------------------------------
static float A2WeightedCombo(const Float_t hg_peak[8], const Float_t hg_cfd[8])
{
    double sw = 0., swt = 0.;
    int    nValid = 0;
    for (int ic = 0; ic < kNCap; ++ic) {
        if (hg_peak[ic] < kHG_minPeak) continue;
        if (hg_cfd[ic]  < -1e5f)       continue;  // CFD-failed sentinel
        double w  = static_cast<double>(hg_peak[ic]) * hg_peak[ic];
        sw  += w;
        swt += w * hg_cfd[ic];
        ++nValid;
    }
    if (nValid < 2 || sw <= 0.) return kNoTime;
    return static_cast<float>(swt / sw);
}

// ---------------------------------------------------------------------------
// FitCellSigma -- fit Gaussian to vector of t_combo values; return sigma [ps]
//
// Returns 0 if too few events or fit fails.
// ---------------------------------------------------------------------------
static double FitCellSigma(const std::vector<float>& v,
                            double histLo, double histHi)
{
    if (static_cast<int>(v.size()) < kMinEvt) return 0.;

    TH1F hc("_us_cell", "", 200, histLo, histHi);
    hc.SetDirectory(nullptr);
    for (float x : v) hc.Fill(x);

    double mu, muE, sig, sigE;
    FitGaussCore(&hc, 2.0, mu, muE, sig, sigE);
    if (sig <= 0.) return 0.;
    return sig * 1000.;   // ns -> ps
}

// ===========================================================================
// Main
// ===========================================================================
void uniformityScan()
{
    TH1::AddDirectory(kFALSE);
    ApplyRADiCALStyle();

    gSystem->mkdir("output/Summary", kTRUE);
    std::cout << "uniformityScan: running spatial uniformity study\n";

    // =========================================================================
    // Per-run outputs
    //   prof2D[r]  : TProfile2D of sigma_t vs (x,y), filled with sigma_t per cell
    //   gRadial[r] : TGraphErrors of sigma_t vs r_beam
    // =========================================================================
    TProfile2D* prof2D[kNRuns] = {};
    for (int r = 0; r < kNRuns; ++r) {
        prof2D[r] = new TProfile2D(
            Form("hSigT2D_%d", r),
            Form("#sigma_{t}(x,y)  --  %.0f GeV;x_{beam} (mm);y_{beam} (mm)",
                 kRuns[r].energy_GeV),
            kNCell, -kGridHalf, kGridHalf,
            kNCell, -kGridHalf, kGridHalf);
        prof2D[r]->SetDirectory(nullptr);
    }

    // Radial sigma_t results indexed [run][rbin]: sigma [ps] and counts
    double sigR   [kNRuns][kNRBins] = {};
    double sigRErr[kNRuns][kNRBins] = {};
    double rCenter[kNRBins]         = {};
    {
        double rbinW = kFiducial_r_timing / kNRBins;
        for (int b = 0; b < kNRBins; ++b)
            rCenter[b] = (b + 0.5) * rbinW;
    }

    // =========================================================================
    // Event loop -- one pass per run
    // =========================================================================
    for (int r = 0; r < kNRuns; ++r) {
        TTree* t = nullptr;
        TFile* fin = OpenNtuple(r, t);
        if (!fin) continue;

        // Pre-scan for beam centroid and CFD offset
        double xc, yc, tcfd, trms;
        ScanRunCenters(t, xc, yc, tcfd, trms);
        // ScanRunCenters calls t->ResetBranchAddresses() before returning

        // The histogram window is set AFTER the event loop from the combo's OWN
        // RMS — NOT from the pooled per-channel trms.  ScanRunCenters inflates
        // trms to ~6.5 ns via the inter-channel cable-delay spread; using 4×trms
        // gave a ±26 ns range (258 ps/bin) that cannot resolve the ~80-130 ps
        // combination core and biased the fitted sigma_t upward (same artefact
        // fixed in timingResolution.C).
        (void) tcfd; (void) trms;

        // Declare branches
        Float_t x_trk, y_trk, mcp_peak, sum_lg, sum_pb;
        Float_t hg_peak[8], hg_cfd[8];
        Bool_t  wc_ok;
        t->SetBranchAddress("x_trk",    &x_trk);
        t->SetBranchAddress("y_trk",    &y_trk);
        t->SetBranchAddress("wc_ok",    &wc_ok);
        t->SetBranchAddress(t->GetBranch("mcp1_peak")?"mcp1_peak":"mcp_peak", &mcp_peak);
        t->SetBranchAddress("sum_lg",   &sum_lg);
        t->SetBranchAddress("sum_pb",   &sum_pb);
        t->SetBranchAddress("hg_peak",   hg_peak);
        // CFD-5% (adopted headline fraction) when present; guarded fallback to
        // CFD-20% for pre-reprocess ntuples.  Consistent with the headline and the
        // per-channel diagnostics (removes the Down-capillary CFD-20% shoulder).
        t->SetBranchAddress(t->GetBranch("hg_cfd05") ? "hg_cfd05" : "hg_cfd", hg_cfd);

        // Accumulators: tvecGrid[ix][iy] and tvecRad[rb]
        std::vector<float> tvecGrid[kNCell][kNCell];
        std::vector<float> tvecRad[kNRBins];
        for (int ix = 0; ix < kNCell; ++ix)
            for (int iy = 0; iy < kNCell; ++iy)
                tvecGrid[ix][iy].reserve(2000);
        for (int b = 0; b < kNRBins; ++b)
            tvecRad[b].reserve(10000);
        std::vector<float> tvecAll;          // all valid combos → sets the fit window
        tvecAll.reserve(40000);

        Long64_t nEv = t->GetEntries();
        for (Long64_t ev = 0; ev < nEv; ++ev) {
            t->GetEntry(ev);

            // Standard event quality cuts
            if (!wc_ok)                                                continue;
            if (mcp_peak < kMCP1_minPeak || mcp_peak > kMCP1_maxPeak) continue;

            float dx   = x_trk - static_cast<float>(xc);
            float dy   = y_trk - static_cast<float>(yc);
            float rbeam = std::sqrt(dx*dx + dy*dy);

            // Timing fiducial
            if (rbeam >= static_cast<float>(kFiducial_r_timing)) continue;

            // Containment cut (skip events with too little LG signal)
            if (sum_lg > kSumLG_centroid &&
                sum_pb >= kPb_maxRatio * sum_lg)               continue;

            // A^2-weighted combo timing
            float tcombo = A2WeightedCombo(hg_peak, hg_cfd);
            if (tcombo <= kNoTime + 1e5f) continue;
            tvecAll.push_back(tcombo);

            // -- Grid cell --
            // Map (dx, dy) into cell indices
            // Cell (0,0) = bottom-left corner: dx in [-3, -1.8), dy in [-3, -1.8)
            int ix = static_cast<int>((dx + kGridHalf) / kCellSize);
            int iy = static_cast<int>((dy + kGridHalf) / kCellSize);
            if (ix >= 0 && ix < kNCell && iy >= 0 && iy < kNCell)
                tvecGrid[ix][iy].push_back(tcombo);

            // -- Radial bin --
            int rb = static_cast<int>(rbeam / (kFiducial_r_timing / kNRBins));
            if (rb >= 0 && rb < kNRBins)
                tvecRad[rb].push_back(tcombo);
        }

        t->ResetBranchAddresses();
        fin->Close(); delete fin;

        // Histogram window from the combo's OWN distribution (robust, matches
        // timingContainmentScan): ±4 σ_combo with a 200 ps floor → ~12 ps/bin.
        double cMean = 0., cRms = 0.;
        if (!tvecAll.empty()) {
            for (float v : tvecAll) cMean += v;
            cMean /= tvecAll.size();
            for (float v : tvecAll) cRms += (v - cMean) * (v - cMean);
            cRms = std::sqrt(cRms / tvecAll.size());
        }
        if (cRms < 0.050) cRms = 0.300;
        double histLo = cMean - 4.0 * cRms;
        double histHi = cMean + 4.0 * cRms;

        // -- Fit each grid cell and fill TProfile2D --
        for (int ix = 0; ix < kNCell; ++ix) {
            for (int iy = 0; iy < kNCell; ++iy) {
                double sig = FitCellSigma(tvecGrid[ix][iy], histLo, histHi);
                if (sig <= 0.) continue;

                // Center of this cell in (dx,dy) space
                double cx = -kGridHalf + (ix + 0.5) * kCellSize;
                double cy = -kGridHalf + (iy + 0.5) * kCellSize;
                // Fill the profile with the fitted sigma_t at the cell center
                // (profile stores running-mean -- filling once with the fitted
                //  value is the correct approach for a derived quantity)
                prof2D[r]->Fill(cx, cy, sig);
            }
        }

        // -- Fit each radial bin --
        for (int b = 0; b < kNRBins; ++b) {
            double sig = FitCellSigma(tvecRad[b], histLo, histHi);
            sigR[r][b] = sig;
            // Propagate statistical uncertainty from Gaussian fit
            int nv = static_cast<int>(tvecRad[b].size());
            sigRErr[r][b] = (nv > 2 && sig > 0.) ? sig / std::sqrt(2.0 * nv) : 0.;
        }

        std::cout << "  " << kRuns[r].label
                  << Form(": centroid (%.1f, %.1f) mm  CFD offset %.2f ns\n",
                           xc, yc, tcfd);
    }

    // =========================================================================
    // Build TGraphErrors for Page 3 (sigma_t vs r_beam)
    // =========================================================================
    TGraphErrors* gRadial[kNRuns] = {};
    double rbinW = kFiducial_r_timing / kNRBins;
    for (int r = 0; r < kNRuns; ++r) {
        gRadial[r] = new TGraphErrors();
        for (int b = 0; b < kNRBins; ++b) {
            if (sigR[r][b] <= 0.) continue;
            int n = gRadial[r]->GetN();
            gRadial[r]->SetPoint(n, rCenter[b], sigR[r][b]);
            gRadial[r]->SetPointError(n, rbinW * 0.5, sigRErr[r][b]);
        }
        gRadial[r]->SetLineColor(kREnergyCols[r]);
        gRadial[r]->SetMarkerColor(kREnergyCols[r]);
        gRadial[r]->SetMarkerStyle(kEMark[r]);
        gRadial[r]->SetMarkerSize(1.2);
        gRadial[r]->SetLineWidth(2);
    }

    // Determine sigma_t range from 150 GeV (last run) for colour scale
    double zMax150 = 0.;
    for (int ix = 1; ix <= kNCell; ++ix)
        for (int iy = 1; iy <= kNCell; ++iy) {
            double v = prof2D[kNRuns - 1]->GetBinContent(ix, iy);
            if (v > zMax150) zMax150 = v;
        }
    if (zMax150 < 10.) zMax150 = 200.;   // safety default

    // Determine sigma_t at 25 GeV for the shared Page 2 scale
    double zMax25 = 0.;
    for (int ix = 1; ix <= kNCell; ++ix)
        for (int iy = 1; iy <= kNCell; ++iy) {
            double v = prof2D[0]->GetBinContent(ix, iy);
            if (v > zMax25) zMax25 = v;
        }
    if (zMax25 < 10.) zMax25 = 400.;
    double zMaxShared = 1.3 * zMax25;

    // =========================================================================
    // Output PDF
    // =========================================================================
    TString outPDF  = "output/Summary/uniformity_scan.pdf";
    TString outROOT = "output/Summary/uniformity_scan.root";

    TCanvas c("c_us", "", 900, 800);

    // =========================================================================
    // Page 1: sigma_t(x,y) at 150 GeV -- single large COLZ map
    // =========================================================================
    c.Clear(); c.cd();
    StylePad(true);   // COLZ right palette, no grid

    TProfile2D* p150 = prof2D[kNRuns - 1];
    p150->SetMinimum(0.);
    p150->SetMaximum(200.);
    p150->GetXaxis()->SetTitle("x_{beam}  (mm)");
    p150->GetYaxis()->SetTitle("y_{beam}  (mm)");
    p150->GetZaxis()->SetTitle("#sigma_{t}  (ps)");
    p150->GetYaxis()->SetTitleOffset(1.35);
    p150->GetZaxis()->SetTitleOffset(1.55);
    p150->Draw("COLZ");

    // Timing fiducial circle (orange)
    TEllipse* elTim = new TEllipse(0., 0.,
                                   kFiducial_r_timing, kFiducial_r_timing);
    elTim->SetLineColor(kOrange+1);
    elTim->SetLineWidth(2);
    elTim->SetLineStyle(1);
    elTim->SetFillStyle(0);
    elTim->Draw("SAME");

    // Energy fiducial circle (red dashed)
    TEllipse* elEng = new TEllipse(0., 0.,
                                   kFiducial_r_energy, kFiducial_r_energy);
    elEng->SetLineColor(kRed+1);
    elEng->SetLineWidth(2);
    elEng->SetLineStyle(2);
    elEng->SetFillStyle(0);
    elEng->Draw("SAME");

    // Annotations
    {
        TLatex ann; ann.SetNDC(); ann.SetTextSize(0.034); ann.SetTextFont(42);
        ann.SetTextColor(kOrange+1);
        ann.DrawLatex(0.15, 0.86, Form("r = %.1f mm  (timing fiducial)",
                                        kFiducial_r_timing));
        ann.SetTextColor(kRed+1);
        ann.DrawLatex(0.15, 0.80, Form("r = %.1f mm  (energy fiducial)",
                                        kFiducial_r_energy));
    }
    DrawPadTitle("A^{2}-weighted combo #sigma_{t} vs beam position  --  150 GeV");
    c.Print(outPDF + "(");

    // =========================================================================
    // Page 2: sigma_t(x,y) at all 6 energies -- 3x2 panel grid
    // =========================================================================
    c.Clear();
    c.Divide(3, 2, 0.005f, 0.005f);

    for (int r = 0; r < kNRuns; ++r) {
        c.cd(r + 1);
        StylePad(true);   // COLZ palette, no grid
        // Smaller margins for multi-panel
        gPad->SetLeftMargin  (0.14f);
        gPad->SetBottomMargin(0.14f);
        gPad->SetTopMargin   (0.12f);
        gPad->SetRightMargin (0.20f);

        prof2D[r]->SetMinimum(0.);
        prof2D[r]->SetMaximum(zMaxShared);

        TProfile2D* pr = static_cast<TProfile2D*>(prof2D[r]->Clone(
                             Form("hClone2D_%d", r)));
        pr->SetDirectory(nullptr);
        pr->GetXaxis()->SetTitle("x_{beam}  (mm)");
        pr->GetYaxis()->SetTitle("y_{beam}  (mm)");
        pr->GetZaxis()->SetTitle("#sigma_{t}  (ps)");
        pr->GetXaxis()->SetTitleSize(0.060f);
        pr->GetYaxis()->SetTitleSize(0.060f);
        pr->GetZaxis()->SetTitleSize(0.058f);
        pr->GetXaxis()->SetLabelSize(0.050f);
        pr->GetYaxis()->SetLabelSize(0.050f);
        pr->GetZaxis()->SetLabelSize(0.050f);
        pr->GetYaxis()->SetTitleOffset(1.10f);
        pr->GetZaxis()->SetTitleOffset(1.25f);
        pr->Draw("COLZ");

        // Timing fiducial circle
        TEllipse* et = new TEllipse(0., 0.,
                                    kFiducial_r_timing, kFiducial_r_timing);
        et->SetLineColor(kOrange+1); et->SetLineWidth(2);
        et->SetLineStyle(1); et->SetFillStyle(0);
        et->Draw("SAME");

        // Panel title
        {
            TLatex tit; tit.SetNDC(); tit.SetTextFont(42);
            tit.SetTextSize(0.072f); tit.SetTextAlign(22);
            float lm = static_cast<float>(gPad->GetLeftMargin());
            float rm = static_cast<float>(gPad->GetRightMargin());
            float xc_pad = lm + (1.0f - lm - rm) * 0.5f;
            float yc_pad = 1.0f - static_cast<float>(gPad->GetTopMargin()) * 0.40f;
            tit.DrawLatex(xc_pad, yc_pad,
                          Form("%.0f GeV", kRuns[r].energy_GeV));
        }
    }

    // Page-level title on canvas frame
    c.cd(0);
    {
        TLatex ptit; ptit.SetNDC(); ptit.SetTextFont(42);
        ptit.SetTextSize(0.025f); ptit.SetTextAlign(22);
        ptit.DrawLatex(0.50f, 0.990f,
            "A^{2}-weighted combo #sigma_{t}(x,y)  --  all energies"
            "  (z_{max} = 1.3#times25 GeV value)");
    }
    c.Print(outPDF);

    // =========================================================================
    // Page 3: sigma_t vs r_beam -- 6 TGraphErrors overlaid
    // =========================================================================
    c.Clear(); c.cd(); StylePad(false, true);   // sidebar legend

    // Y range
    double yMaxR = 0.;
    for (int r = 0; r < kNRuns; ++r)
        for (int b = 0; b < kNRBins; ++b)
            if (sigR[r][b] > yMaxR) yMaxR = sigR[r][b];
    if (yMaxR < 10.) yMaxR = 400.;

    TH1F* frame3 = static_cast<TH1F*>(
        c.DrawFrame(0., 0., kFiducial_r_timing * 1.08, 1.25 * yMaxR,
                    ";r_{beam}  (mm);A^{2}-weighted combo #sigma_{t}  (ps)"));
    frame3->GetXaxis()->SetTitleSize(0.050f);
    frame3->GetYaxis()->SetTitleSize(0.050f);
    frame3->GetYaxis()->SetTitleOffset(1.30f);

    // Dashed vertical line at timing fiducial radius
    TLine* lFid = new TLine(kFiducial_r_timing, 0.,
                             kFiducial_r_timing, 1.25 * yMaxR);
    lFid->SetLineColor(kOrange+1);
    lFid->SetLineStyle(2);
    lFid->SetLineWidth(2);
    lFid->Draw("SAME");
    {
        TLatex ann; ann.SetNDC(); ann.SetTextSize(0.032f);
        ann.SetTextColor(kOrange+1);
        // Annotate just to the left of the fiducial line
        double lm = gPad->GetLeftMargin();
        double rm = gPad->GetRightMargin();
        double xFrac = lm + (kFiducial_r_timing / (kFiducial_r_timing * 1.08))
                       * (1.0 - lm - rm);
        ann.DrawLatex(static_cast<float>(xFrac - 0.12), 0.92f,
                      "timing fiducial");
    }

    TLegend* leg3 = MakeLegend(kNRuns);
    for (int r = 0; r < kNRuns; ++r) {
        if (!gRadial[r] || gRadial[r]->GetN() == 0) continue;
        gRadial[r]->Draw("PL SAME");
        leg3->AddEntry(gRadial[r],
                       Form("%.0f GeV", kRuns[r].energy_GeV), "lp");
    }
    leg3->Draw();

    DrawPadTitle("#sigma_{t} vs beam radius  --  spatial uniformity");
    c.Print(outPDF + ")");

    // =========================================================================
    // Write ROOT file
    // =========================================================================
    TFile* fOut = new TFile(outROOT, "RECREATE");
    for (int r = 0; r < kNRuns; ++r) {
        prof2D[r]->Write(Form("hSigT2D_%s", kRuns[r].label.Data()));
        if (gRadial[r] && gRadial[r]->GetN() > 0)
            gRadial[r]->Write(Form("gSigT_radius_%s",
                                   kRuns[r].label.Data()));
    }
    fOut->Close(); delete fOut;

    // Cleanup
    for (int r = 0; r < kNRuns; ++r) {
        delete prof2D[r];
        delete gRadial[r];
    }
    delete lFid;
    delete leg3;
    delete elTim;
    delete elEng;

    std::cout << "uniformityScan: done -> " << outPDF  << "\n"
              << "                root -> " << outROOT << "\n";
}

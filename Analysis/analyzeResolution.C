// ============================================================================
// analyzeResolution.C — multi-energy resolution analysis for RADiCAL
// ============================================================================
//
// Reads the analysis ntuples produced by processRun.C (one per beam energy)
// and produces the key physics plots:
//
//   Per energy (in Analysis/Output/<label>/):
//     energy_distribution.pdf  — LG amplitude sum + Gaussian fit
//     timing_distributions.pdf — per-capillary Δt = t_HG − t_MCP distributions
//     hit_map.pdf              — 2D mean signal vs beam impact position
//
//   Summary (in Analysis/Output/Summary/):
//     resolution_summary.pdf   — energy resolution, timing resolution, linearity
//     summary.root             — TGraphErrors for the three summary plots
//
// Usage (from repository root):
//   root -l -b -q 'Analysis/analyzeResolution.C+'
//
// Prerequisites: run processRun.C for each energy first (see runAll.sh).
// ============================================================================

#include "ChannelConfig.h"
#include "PlotUtils.h"     // StylePad, FitGaussCore, DrawFitOverlay, ScanRunCenters

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile2D.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLatex.h"
#include "TPad.h"
#include "TMath.h"
#include "TArc.h"
#include "TBox.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

// ScanRunCenters, FitGaussCore, StylePad — provided by PlotUtils.h

// ---------------------------------------------------------------------------
// Beam momentum spread at the CERN SPS H2 beam line (sigma_p/p, fractional).
//
// The H2 momentum bite is set by the collimator slits as a roughly
// ENERGY-INDEPENDENT fraction of the central momentum.  Published H2
// characterisation: maximum Delta p/p = +-2.0%, with ~1% single-particle
// momentum resolution typical for a resolution-optimised secondary beam
// (CERN North Area H2 beam-line handbook, nahandbook.web.cern.ch/h2).
// We therefore use a single literature-typical value of 1.0% at all energies.
//
// CAVEAT: this is a LITERATURE-TYPICAL value, NOT the run-specific collimator
// setting for the May-2023 period (which would come from the PS/SPS tune
// sheets).  IMPACT IS NEGLIGIBLE HERE: our measured sigma_E/E is 11-19%
// (leakage-dominated in the compact 14 mm module), so subtracting 1% in
// quadrature changes the constant term by < 0.1% -- the detector-only and
// measured fits are essentially identical.  The framework is in place should
// run-specific numbers (or a larger-module measurement) ever require it.
// ---------------------------------------------------------------------------
static const double kBeamSpread[6] = {
    0.010,  // 25 GeV
    0.010,  // 50 GeV
    0.010,  // 75 GeV
    0.010,  // 100 GeV
    0.010,  // 125 GeV
    0.010,  // 150 GeV
};

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------
void analyzeResolution()
{
    ApplyRADiCALStyle();   // RADiCALStyle.h via PlotUtils.h

    gSystem->mkdir("Analysis/Output/Summary", kTRUE);

    // Storage for multi-energy summary
    std::vector<double> vE, vEErr;
    std::vector<double> vERes, vEResErr;
    std::vector<double> vMeanLG, vMeanLGErr;
    std::vector<double> vTRes[8], vTResErr[8];   // per capillary

    // =========================================================================
    // Per-energy loop
    // =========================================================================
    for (int iRun = 0; iRun < kNRuns; ++iRun)
    {
        const RunCfg& rc   = kRuns[iRun];
        TString ntupleFile = TString("Analysis/Output/") + rc.label + "/ntuple.root";

        TFile* f = TFile::Open(ntupleFile);
        if (!f || f->IsZombie()) {
            std::cout << "[analyzeResolution] Skipping " << rc.label
                      << " — ntuple not found (" << ntupleFile << ")\n"
                      << "   Run processRun.C for this energy first.\n";
            continue;
        }
        TTree* t = (TTree*)f->Get("rad");
        if (!t || t->GetEntries() == 0) { f->Close(); continue; }

        Long64_t nEntries = t->GetEntries();
        std::cout << "[analyzeResolution] " << rc.label
                  << " — " << nEntries << " events\n";

        // ---------------------------------------------------------------------
        // Pre-scan: derive beam centroid and CFD offset from this run's data.
        // This avoids hardcoding position or timing constants that differ
        // between runs / beam optics settings.
        // ---------------------------------------------------------------------
        double x_center, y_center, t_cfd_offset, t_cfd_rms;
        ScanRunCenters(t, x_center, y_center, t_cfd_offset, t_cfd_rms);
        // Timing histogram: ±5σ around mean, minimum ±2 ns window
        double t_win = std::max(5.0 * t_cfd_rms, 2.0);
        double t_lo  = t_cfd_offset - t_win;
        double t_hi  = t_cfd_offset + t_win;

        // ---------------------------------------------------------------------
        // Book histograms (using data-derived ranges)
        // ---------------------------------------------------------------------

        // Energy: sum of 8 LG amplitudes
        TH1F* hSumLG = new TH1F(
            Form("hSumLG_%s", rc.label.Data()),
            ";Sum_{LG} (mV);Events",
            300, 0, 7000);

        // PbGlass reference sum
        TH1F* hSumPb = new TH1F(
            Form("hSumPb_%s", rc.label.Data()),
            ";Sum_{PbGlass} (mV);Events",
            300, 0, 7000);

        // Timing: Δt = t_HG[i] − t_MCP, centred on measured offset
        TH1F* hTiming[8];
        for (int i = 0; i < 8; ++i)
            hTiming[i] = new TH1F(
                Form("hTiming_%s_%d", rc.label.Data(), i),
                ";t_{HG}#minust_{MCP} (ns);Events",
                400, t_lo, t_hi);

        // 2D hit map: mean LG sum vs beam position (centred on beam)
        double hm_hw  = 20.;  // half-width of the hit map [mm]
        int    nBinsHM = static_cast<int>(std::round(2.*hm_hw / kWC_resBin)); // 1 mm/bin
        TProfile2D* hHitMap = new TProfile2D(
            Form("hHitMap_%s", rc.label.Data()),
            ";x_{WC} (mm);y_{WC} (mm);Sum_{LG} (mV)",
            nBinsHM, x_center-hm_hw, x_center+hm_hw,
            nBinsHM, y_center-hm_hw, y_center+hm_hw);

        // Per-capillary amplitude maps (centred on beam)
        double cm_hw   = 10.;  // half-width of the per-cap map [mm]
        int    nBinsCM = static_cast<int>(std::round(2.*cm_hw / kWC_resBin)); // 1 mm/bin
        TProfile2D* hCapMap[8];
        for (int i = 0; i < 8; ++i)
            hCapMap[i] = new TProfile2D(
                Form("hCapMap_%s_%d", rc.label.Data(), i),
                ";x_{WC} (mm);y_{WC} (mm);A_{LG} (mV)",
                nBinsCM, x_center-cm_hw, x_center+cm_hw,
                nBinsCM, y_center-cm_hw, y_center+cm_hw);

        // ---------------------------------------------------------------------
        // Fill histograms
        // ---------------------------------------------------------------------
        Bool_t  wc_ok;
        Float_t x_trk, y_trk, mcp_peak, sum_lg, sum_pb;
        Float_t hg_cfd[8], hg_peak[8], lg_peak[8];

        t->SetBranchAddress("wc_ok",    &wc_ok);
        t->SetBranchAddress("x_trk",    &x_trk);
        t->SetBranchAddress("y_trk",    &y_trk);
        t->SetBranchAddress("mcp_peak", &mcp_peak);
        t->SetBranchAddress("sum_lg",   &sum_lg);
        t->SetBranchAddress("sum_pb",   &sum_pb);
        // Per-channel timing uses CFD-5% (the adopted headline fraction), not
        // CFD-20%.  On the Down capillaries the leading-edge SHAPE jitters more
        // pulse-to-pulse high on the edge (not a mean-slope effect — the mean edge
        // is steeper at 20%), producing a broad, shouldered peak (~2x wider);
        // CFD-5% times low on the edge where it is most reproducible and removes it
        // (see timingMethods.C page 3, edgeMechanism.C, elbowInvestigation.C).
        // Falls back to hg_cfd (20%) for pre-reprocess ntuples without hg_cfd05.
        t->SetBranchAddress(t->GetBranch("hg_cfd05") ? "hg_cfd05" : "hg_cfd", hg_cfd);
        t->SetBranchAddress("hg_peak",   hg_peak);
        t->SetBranchAddress("lg_peak",   lg_peak);

        long nInFid = 0;
        for (Long64_t ev = 0; ev < nEntries; ++ev) {
            t->GetEntry(ev);

            if (mcp_peak < kMCP_minPeak_E) continue;   // require MCP signal

            // Hit map — all wire-chamber-OK events
            if (wc_ok) hHitMap->Fill(x_trk, y_trk, sum_lg);

            // Dynamic fiducial cut: within kFiducial_r of data-derived centroid
            if (!wc_ok) continue;
            float dx = x_trk - static_cast<float>(x_center);
            float dy = y_trk - static_cast<float>(y_center);
            if (std::sqrt(dx*dx + dy*dy) >= static_cast<float>(kFiducial_r_energy)) continue;
            ++nInFid;

            hSumLG->Fill(sum_lg);
            hSumPb->Fill(sum_pb);

            for (int i = 0; i < 8; ++i) {
                hCapMap[i]->Fill(x_trk, y_trk, lg_peak[i]);
                // Require hg_peak >= kHG_minPeak: sub-threshold pulses that
                // accidentally produce a CFD crossing add noise-dominated timing
                // entries and artificially degrade the reported sigma.
                if (hg_cfd[i] > -1e5f && hg_peak[i] >= kHG_minPeak)
                    hTiming[i]->Fill(hg_cfd[i]);
            }
        }

        std::cout << "  Fiducial events: " << nInFid
                  << " / " << nEntries
                  << Form(" (%.1f%%)\n", 100.*nInFid/nEntries);

        // ---------------------------------------------------------------------
        // Fit and extract physics quantities
        // ---------------------------------------------------------------------
        double mu, muErr, sig, sigErr;
        // nsig = 2.0: wider window captures the full Gaussian core without
        // being overly sensitive to non-Gaussian tails from shower leakage or
        // hadron contamination.  1.5 was too tight and could clip on a tail bin.
        FitGaussCore(hSumLG, 2.0, mu, muErr, sig, sigErr);
        double res    = (mu > 0) ? sig / mu : 0;
        double resErr = (mu > 0 && sig > 0)
                       ? res * std::sqrt(sigErr*sigErr/(sig*sig) + muErr*muErr/(mu*mu))
                       : 0;

        if (mu > 0) {
            vE.push_back(rc.energy_GeV);
            vEErr.push_back(0.);
            vERes.push_back(res * 100.);
            vEResErr.push_back(resErr * 100.);
            vMeanLG.push_back(mu);
            vMeanLGErr.push_back(muErr);
        }

        for (int i = 0; i < 8; ++i) {
            double tm, tmE, ts, tsE;
            FitGaussCore(hTiming[i], 2.0, tm, tmE, ts, tsE);
            if (ts > 0) {
                vTRes[i].push_back(ts * 1000.);     // ns → ps
                vTResErr[i].push_back(tsE * 1000.);
            } else {
                vTRes[i].push_back(0.);
                vTResErr[i].push_back(0.);
            }
        }

        // ---------------------------------------------------------------------
        // Per-energy plots
        // ---------------------------------------------------------------------
        TString outDir = TString("Analysis/Output/") + rc.label;

        // — Energy distribution —
        {
            TCanvas c("c_en", "", 700, 550);
            StylePad();
            hSumLG->SetLineColor(kBlue+1);
            hSumLG->SetLineWidth(2);
            hSumLG->Draw("HIST");

            TLatex ptit;
            ptit.SetNDC(); ptit.SetTextSize(0.046); ptit.SetTextAlign(22);
            ptit.DrawLatex(0.55, 0.94,
                Form("%.0f GeV  #minus  LG energy sum", rc.energy_GeV));

            if (mu > 0 && sig > 0) {
                TF1 fDraw("fDraw", "gaus",
                          mu - 2.*sig, mu + 2.*sig);
                fDraw.SetParameters(hSumLG->GetMaximum(), mu, sig);
                fDraw.SetLineColor(kRed+1);
                fDraw.SetLineWidth(2);
                fDraw.Draw("SAME");

                TLatex lat;
                lat.SetNDC();
                lat.SetTextSize(0.042);
                lat.DrawLatex(0.60, 0.77,
                    Form("#sigma/E = %.1f%%", res*100.));
                lat.DrawLatex(0.60, 0.70,
                    Form("#mu = %.0f mV", mu));
            }
            c.Print(outDir + "/energy_distribution.pdf");
        }

        // — Timing distributions —
        {
            TCanvas c("c_tim", "", 1400, 700);
            c.Divide(4, 2, 0.01, 0.01);
            for (int i = 0; i < 8; ++i) {
                c.cd(i+1);
                StylePad();
                hTiming[i]->SetLineColor(kBlue+1);
                hTiming[i]->SetLineWidth(2);

                // Fit before drawing so we can zoom x-axis to ±5.5σ around peak
                double tm2, tmE2, ts2, tsE2;
                FitGaussCore(hTiming[i], 2.0, tm2, tmE2, ts2, tsE2);
                if (ts2 > 0)
                    hTiming[i]->GetXaxis()->SetRangeUser(tm2 - 5.5*ts2,
                                                          tm2 + 5.5*ts2);

                hTiming[i]->Draw("HIST");
                DrawPadTitle(Form("%.0f GeV  %s", rc.energy_GeV, kCap[i].name));

                if (ts2 > 0) {
                    TF1 ftDraw("ftDraw", "gaus",
                               tm2 - 2.5*ts2, tm2 + 2.5*ts2);
                    ftDraw.SetParameters(hTiming[i]->GetMaximum(), tm2, ts2);
                    ftDraw.SetLineColor(kRed+1);
                    ftDraw.SetLineWidth(2);
                    ftDraw.Draw("SAME");
                    TLatex lat;
                    lat.SetNDC();
                    lat.SetTextSize(0.065);
                    lat.DrawLatex(0.62, 0.77,
                        Form("#sigma = %.0f ps", ts2*1000.));
                }
            }
            c.Print(outDir + "/timing_distributions.pdf");
        }

        // — Hit map —
        {
            TCanvas c("c_hm", "", 700, 600);
            c.SetRightMargin(0.16);
            StylePad(true);            // rightPalette=true → grid off, palette bar margin
            gPad->SetRightMargin(0.16);
            hHitMap->Draw("COLZ");

            TLatex ptit;
            ptit.SetNDC(); ptit.SetTextSize(0.046); ptit.SetTextAlign(22);
            ptit.DrawLatex(0.47, 0.94,
                Form("%.0f GeV  beam hit map", rc.energy_GeV));
            // Mark the DATA-DERIVED beam centroid (solid red = energy fiducial,
            // dashed red = timing fiducial — both drawn for reference)
            TArc* fidE = new TArc(x_center, y_center, kFiducial_r_energy);
            fidE->SetLineColor(kRed); fidE->SetLineWidth(2);
            fidE->SetLineStyle(1);   fidE->SetFillStyle(0);
            fidE->Draw("SAME");
            TArc* fidT = new TArc(x_center, y_center, kFiducial_r_timing);
            fidT->SetLineColor(kRed); fidT->SetLineWidth(1);
            fidT->SetLineStyle(2);   fidT->SetFillStyle(0);
            fidT->Draw("SAME");
            // 14x14 mm calorimeter face
            TBox* calFace = new TBox(kCalo_x0-7., kCalo_y0-7.,
                                     kCalo_x0+7., kCalo_y0+7.);
            calFace->SetLineColor(kWhite);
            calFace->SetLineWidth(2);
            calFace->SetFillStyle(0);
            calFace->Draw("SAME");
            c.Print(outDir + "/hit_map.pdf");
        }

        // — Per-capillary LG amplitude maps —
        {
            TCanvas c("c_caps", "", 1400, 700);
            c.Divide(4, 2, 0.02, 0.02);
            for (int i = 0; i < 8; ++i) {
                c.cd(i+1);
                gPad->SetRightMargin(0.16);
                StylePad(true);        // rightPalette=true → grid off, palette bar margin
                gPad->SetRightMargin(0.16);
                hCapMap[i]->Draw("COLZ");

                TLatex ptit;
                ptit.SetNDC(); ptit.SetTextSize(0.072); ptit.SetTextAlign(22);
                ptit.DrawLatex(0.46, 0.94,
                    Form("%.0f GeV  %s", rc.energy_GeV, kCap[i].name));
            }
            c.Print(outDir + "/capillary_maps.pdf");
        }

        f->Close();
        std::cout << "  -> sigma/E = " << res*100. << "% at "
                  << rc.energy_GeV << " GeV\n";
    } // end per-energy loop

    if (vE.empty()) {
        std::cout << "[analyzeResolution] No data — run processRun.C first.\n";
        return;
    }

    // =========================================================================
    // Multi-energy summary plots
    // =========================================================================
    const int N = static_cast<int>(vE.size());

    // --- Energy resolution ---
    TGraphErrors* gERes = new TGraphErrors(N);
    gERes->SetName("gEnergyResolution");
    gERes->SetTitle("RADiCAL Energy Resolution;"
                    "Beam Energy (GeV);#sigma_{E}/E (%)");
    for (int i = 0; i < N; ++i) {
        gERes->SetPoint(i,      vE[i],    vERes[i]);
        gERes->SetPointError(i, vEErr[i], vEResErr[i]);
    }
    gERes->SetMarkerStyle(20);
    gERes->SetMarkerSize(1.4);
    gERes->SetMarkerColor(kBlue+1);
    gERes->SetLineColor(kBlue+1);
    gERes->SetLineWidth(2);

    // Fit: sigma/E = sqrt( (a/sqrt(E))^2 + b^2 )  (a and b in %)
    TF1* fERes = new TF1("fERes", "sqrt(([0]*[0])/x + [1]*[1])", 10, 200);
    fERes->SetParameters(10., 1.);
    fERes->SetParNames("a (%#sqrt{GeV})", "b (%)");
    fERes->SetLineColor(kRed+1);
    fERes->SetLineWidth(2);

    // --- Linearity ---
    TGraphErrors* gLin = new TGraphErrors(N);
    gLin->SetName("gLinearity");
    gLin->SetTitle("RADiCAL Linearity;"
                   "Beam Energy (GeV);Mean Sum_{LG} (mV)");
    for (int i = 0; i < N; ++i) {
        gLin->SetPoint(i,      vE[i],    vMeanLG[i]);
        gLin->SetPointError(i, vEErr[i], vMeanLGErr[i]);
    }
    gLin->SetMarkerStyle(21);
    gLin->SetMarkerSize(1.4);
    gLin->SetMarkerColor(kBlue+1);
    gLin->SetLineColor(kBlue+1);
    gLin->SetLineWidth(2);

    // Use a+b*E (not forced through origin): a non-zero intercept reveals a
    // constant pedestal offset or non-linearity in the lowest-energy response.
    TF1* fLin = new TF1("fLin", "[0] + [1]*x", 0, 200);
    fLin->SetParameter(0, 0.);
    fLin->SetParameter(1, N > 0 ? vMeanLG[0]/vE[0] : 20.);
    fLin->SetParNames("offset (mV)", "slope (mV/GeV)");
    fLin->SetLineColor(kRed+1);
    fLin->SetLineWidth(2);

    // --- Average timing resolution (mean over all 8 capillaries) ---
    TGraphErrors* gTRes = new TGraphErrors(N);
    gTRes->SetName("gTimingResolution");
    gTRes->SetTitle("RADiCAL Timing Resolution (avg. 8 capillaries);"
                    "Beam Energy (GeV);#sigma_{t} (ps)");
    for (int i = 0; i < N; ++i) {
        double avg = 0, errSq = 0;
        int ng = 0;
        for (int j = 0; j < 8; ++j) {
            if (static_cast<int>(vTRes[j].size()) > i && vTRes[j][i] > 0) {
                avg   += vTRes[j][i];
                errSq += vTResErr[j][i] * vTResErr[j][i];
                ++ng;
            }
        }
        if (ng > 0) {
            gTRes->SetPoint(i,      vE[i],    avg/ng);
            gTRes->SetPointError(i, vEErr[i], std::sqrt(errSq)/ng);
        }
    }
    gTRes->SetMarkerStyle(22);
    gTRes->SetMarkerSize(1.4);
    gTRes->SetMarkerColor(kBlue+1);
    gTRes->SetLineColor(kBlue+1);
    gTRes->SetLineWidth(2);

    TF1* fTRes = new TF1("fTRes", "sqrt(([0]*[0])/x + [1]*[1])", 10, 200);
    fTRes->SetParameters(1500., 100.);
    fTRes->SetParNames("a (ps#sqrt{GeV})", "b (ps)");
    fTRes->SetLineColor(kRed+1);
    fTRes->SetLineWidth(2);

    // Per-capillary timing resolution (overlay at each energy)
    TGraphErrors* gTResCap[8];
    const int colors[8] = {kBlue+1, kRed+1, kGreen+2, kMagenta+1,
                            kOrange+2, kCyan+1, kViolet+1, kYellow+2};
    for (int j = 0; j < 8; ++j) {
        gTResCap[j] = new TGraphErrors();
        gTResCap[j]->SetName(Form("gTRes_%s", kCap[j].name));
        gTResCap[j]->SetTitle(kCap[j].name);
        gTResCap[j]->SetMarkerStyle(20+j);
        gTResCap[j]->SetMarkerSize(1.1);
        gTResCap[j]->SetMarkerColor(colors[j]);
        gTResCap[j]->SetLineColor(colors[j]);
        for (int i = 0; i < N; ++i) {
            if (static_cast<int>(vTRes[j].size()) > i && vTRes[j][i] > 0)
                gTResCap[j]->SetPoint(gTResCap[j]->GetN(),
                                      vE[i], vTRes[j][i]);
        }
    }

    // =========================================================================
    // Print summary canvas
    // =========================================================================
    TCanvas* cSum = new TCanvas("cSummary", "RADiCAL Summary", 1800, 600);
    cSum->Divide(3, 1, 0.02, 0.02);

    // Pad 1: energy resolution
    cSum->cd(1);
    StylePad();
    gERes->GetXaxis()->SetRangeUser(0, 180);
    gERes->Draw("AP");
    {
        TLatex tit1; tit1.SetNDC(); tit1.SetTextSize(0.052); tit1.SetTextAlign(22);
        tit1.DrawLatex(0.55, 0.94, "Energy Resolution");
    }
    if (N >= 3) {
        gERes->Fit(fERes, "RQ");
        fERes->Draw("SAME");
        // Fit parameters in a tidy lower-left block (clear of the descending curve)
        TLatex lat;
        lat.SetNDC();
        lat.SetTextSize(0.038);
        lat.SetTextColor(kBlue+1);
        lat.DrawLatex(0.20, 0.40,
            Form("measured: a = %.1f%%#sqrt{GeV}, b = %.2f%%",
                 fERes->GetParameter(0), fERes->GetParameter(1)));
    }

    // Beam-spread-subtracted fit
    TGraphErrors* gERes_det = new TGraphErrors(N);
    gERes_det->SetName("gEnergyResolution_detector");
    gERes_det->SetTitle("RADiCAL Energy Resolution (detector only)");
    for (int i = 0; i < N; ++i) {
        // Find the beam spread for this energy
        int ie = -1;
        for (int r = 0; r < kNRuns; ++r)
            if (std::fabs(kRuns[r].energy_GeV - vE[i]) < 1.) { ie = r; break; }
        double beam_frac = (ie >= 0) ? kBeamSpread[ie] * 100. : 0.;  // convert to %
        double res_meas  = vERes[i];  // measured sigma/E in %
        // Subtract beam spread in quadrature
        double res_det2  = res_meas*res_meas - beam_frac*beam_frac;
        double res_det   = (res_det2 > 0.) ? std::sqrt(res_det2) : 0.;
        // Propagate error conservatively (treat beam spread as exact)
        double res_det_err = (res_det > 0.) ? vEResErr[i] * res_meas / res_det : 0.;
        gERes_det->SetPoint(i, vE[i], res_det);
        gERes_det->SetPointError(i, vEErr[i], res_det_err);
    }
    gERes_det->SetMarkerStyle(24);  // open circle = corrected
    gERes_det->SetMarkerSize(1.4);
    gERes_det->SetMarkerColor(kRed+1);
    gERes_det->SetLineColor(kRed+1);
    gERes_det->SetLineWidth(2);
    gERes_det->Draw("P SAME");

    TF1* fERes_det = new TF1("fERes_det", "sqrt(([0]*[0])/x + [1]*[1])", 10, 200);
    fERes_det->SetParameters(10., 1.);
    fERes_det->SetLineColor(kRed+1); fERes_det->SetLineStyle(2); fERes_det->SetLineWidth(2);
    if (N >= 3) {
        gERes_det->Fit(fERes_det, "RQ");
        fERes_det->Draw("SAME");
        TLatex latDet; latDet.SetNDC(); latDet.SetTextSize(0.038);
        latDet.SetTextColor(kRed+1);
        latDet.DrawLatex(0.20, 0.33,
            Form("detector-only: a = %.1f%%#sqrt{GeV}, b = %.2f%%",
                 fERes_det->GetParameter(0), fERes_det->GetParameter(1)));
    }
    // Legend distinguishing measured vs detector-only
    TLegend* legERes = new TLegend(0.50, 0.66, 0.95, 0.80);
    legERes->SetBorderSize(0); legERes->SetFillStyle(0); legERes->SetTextSize(0.034);
    legERes->AddEntry(gERes,     "measured #sigma_{E}/E", "p");
    legERes->AddEntry(gERes_det, "detector-only (prelim.)", "p");
    legERes->Draw();

    // Pad 2: linearity
    cSum->cd(2);
    StylePad();
    gLin->GetXaxis()->SetRangeUser(0, 180);
    gLin->Draw("AP");
    {
        TLatex tit2; tit2.SetNDC(); tit2.SetTextSize(0.052); tit2.SetTextAlign(22);
        tit2.DrawLatex(0.55, 0.94, "Linearity");
    }
    gLin->Fit(fLin, "RQ");
    fLin->Draw("SAME");
    {
        TLatex lat2; lat2.SetNDC(); lat2.SetTextSize(0.042);
        lat2.DrawLatex(0.19, 0.38,
            Form("slope  = %.1f mV/GeV", fLin->GetParameter(1)));
        lat2.DrawLatex(0.19, 0.30,
            Form("offset = %.0f mV", fLin->GetParameter(0)));
    }

    // Pad 3: timing resolution (single-capillary average)
    cSum->cd(3);
    StylePad();
    gTRes->GetXaxis()->SetRangeUser(0, 180);
    gTRes->Draw("AP");
    {
        TLatex tit3; tit3.SetNDC(); tit3.SetTextSize(0.052); tit3.SetTextAlign(22);
        tit3.DrawLatex(0.55, 0.94, "Timing Resolution (avg. 8 cap.)");
    }
    if (N >= 3) {
        gTRes->Fit(fTRes, "RQ");
        fTRes->Draw("SAME");
        TLatex lat3; lat3.SetNDC(); lat3.SetTextSize(0.042);
        lat3.DrawLatex(0.19, 0.38,
            Form("a = %.0fps#sqrt{GeV}", fTRes->GetParameter(0)));
        lat3.DrawLatex(0.19, 0.30,
            Form("b = %.0fps",           fTRes->GetParameter(1)));
    }

    cSum->Print("Analysis/Output/Summary/resolution_summary.pdf");

    // Per-capillary timing overlay
    // Determine a sensible y-max from the actual data (add 20% headroom)
    double tCapYMax = 100.;
    for (int j = 0; j < 8; ++j)
        for (int i = 0; i < gTResCap[j]->GetN(); ++i)
            tCapYMax = std::max(tCapYMax, gTResCap[j]->GetY()[i]);
    tCapYMax *= 1.25;

    TCanvas* cTCap = new TCanvas("cTimingCap", "Timing per capillary", 800, 600);
    StylePad();
    bool first = true;
    for (int j = 0; j < 8; ++j) {
        if (gTResCap[j]->GetN() == 0) continue;
        gTResCap[j]->GetXaxis()->SetRangeUser(0, 180);
        gTResCap[j]->GetYaxis()->SetRangeUser(0, tCapYMax);
        gTResCap[j]->GetXaxis()->SetTitle("Beam Energy (GeV)");
        gTResCap[j]->GetYaxis()->SetTitle("#sigma_{t} (ps)");
        gTResCap[j]->Draw(first ? "AP" : "P SAME");
        first = false;
    }
    // Draw title AFTER graphs so Draw("AP") doesn't overwrite it
    {
        TLatex tit; tit.SetNDC(); tit.SetTextSize(0.046); tit.SetTextAlign(22);
        tit.DrawLatex(0.55, 0.95, "Timing Resolution per Capillary");
    }
    TLegend* leg = new TLegend(0.65, 0.55, 0.92, 0.92);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.032);
    for (int j = 0; j < 8; ++j)
        if (gTResCap[j]->GetN() > 0)
            leg->AddEntry(gTResCap[j], kCap[j].name, "lp");
    leg->Draw();
    cTCap->Print("Analysis/Output/Summary/timing_per_capillary.pdf");

    // =========================================================================
    // Write summary ROOT file
    // =========================================================================
    TFile* fSumOut = new TFile("Analysis/Output/Summary/summary.root", "RECREATE");
    gERes->Write();
    gERes_det->Write();
    gLin->Write();
    gTRes->Write();
    for (int j = 0; j < 8; ++j) gTResCap[j]->Write();
    fSumOut->Close();

    std::cout << "\n[analyzeResolution] Summary written to Analysis/Output/Summary/\n";
    std::cout << "  resolution_summary.pdf\n";
    std::cout << "  timing_per_capillary.pdf\n";
    std::cout << "  summary.root\n";
}

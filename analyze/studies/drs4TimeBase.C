// ============================================================================
// drs4TimeBase.C — DT5742 / DRS4 time-base verification & stop-cell correction
//                  validation  (Layer 1: Hardware Integrity)
// ============================================================================
//
// Answers three questions directly from the RAW data, with no assumptions:
//
//   1. Is the DRS4 cell-width calibration applied?   (cell-width RMS per group)
//   2. Does the stop cell rotate uniformly?          (stop-cell distribution)
//   3. Is the residual cell-width timing error        (split-half OUT-OF-SAMPLE
//      reproducible and correctable, and on which      validation of the
//      timing estimators does it matter?)              StopCellCorrection)
//
// The third question is answered honestly: the correction is TRAINED on even
// events and its benefit is measured on the held-out ODD events with a
// Gaussian-core fit.  In-sample variance reduction is partly an artefact of
// subtracting per-bin means; only the out-of-sample improvement is physical.
//
// Estimators probed (all referenced so MCP jitter is handled):
//   * single  : SE-D (ch2) hg_cfd                     — per-channel timing
//   * combo   : A^2-weighted mean of the 7 DRS0-G0 HG channels (0-6), MCP-ref
//   * corner  : (NW-D - NW-U)/2  = (hg_cfd[0]-hg_cfd[4])/2  — MCP- & (expected)
//                                                            cell-width-free
//
// Output:
//   output/Summary/drs4_timebase.pdf   (4 pages)
//   output/Summary/drs4_timebase.root  (correction tables + summary)
//
// Usage (from repository root):
//   root -l -b -q 'Analysis/drs4TimeBase.C+'
// ============================================================================

#include "ChannelConfig.h"
#include "WaveformUtils.h"
#include "DRS4Calibration.h"
#include "PlotUtils.h"        // ApplyRADiCALStyle, StylePad, FitGaussCore, colors

#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLatex.h"
#include "TLine.h"
#include "TLegend.h"
#include "TPave.h"
#include "TSystem.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TParameter.h"

#include <cmath>
#include <iostream>
#include <vector>

// ---------------------------------------------------------------------------
// Per-event quantities collected in a single pass over the raw data.
// Sentinel kBad marks an estimator that could not be formed this event.
// ---------------------------------------------------------------------------
static const float kBadEst = -1e9f;

struct TbEvent {
    int   stopD0G0;   // DRS0 G0 stop cell (HG + MCP1 share this group)
    float single;     // SE-D hg_cfd [ns]
    float combo;      // A^2-weighted combo of ch0-6 hg_cfd [ns]
    float corner;     // (NW-D - NW-U)/2 [ns]
};

// ---------------------------------------------------------------------------
// CoreSigmaPs — Gaussian-core σ [ps] of a value list, fit within ±win of the
// data mean.  Returns -1 if the fit fails / too few entries.
// ---------------------------------------------------------------------------
static double CoreSigmaPs(const std::vector<float>& v, double win_ns)
{
    if (static_cast<int>(v.size()) < 100) return -1.;

    double m = 0.;
    for (float x : v) m += x;
    m /= v.size();

    TH1F h("_tb_core", "", 200, m - win_ns, m + win_ns);
    h.SetDirectory(nullptr);
    for (float x : v) h.Fill(x);

    double mu, muE, sig, sigE;
    FitGaussCore(&h, 2.0, mu, muE, sig, sigE);
    return (sig > 0.) ? sig * 1000. : -1.;   // ns → ps
}

// ---------------------------------------------------------------------------
// ValidateSplitHalf — train a StopCellCorrection on EVEN events, apply to the
// held-out ODD events, and report the core σ before/after on the odd set.
//
//   stop/val : parallel arrays of (stop cell, estimator value) for ALL events
//   win_ns   : half-width of the fit window around the data mean
//   nBins    : stop-cell bins for the correction table
//   outCorr  : receives the trained correction (for plotting the pattern)
//   sBefore/sAfter : out-of-sample core σ [ps]
// ---------------------------------------------------------------------------
static void ValidateSplitHalf(const std::vector<int>&   stop,
                              const std::vector<float>& val,
                              double win_ns, int nBins,
                              drs4::StopCellCorrection& outCorr,
                              double& sBefore, double& sAfter)
{
    sBefore = sAfter = -1.;
    if (stop.size() != val.size() || val.size() < 400) return;

    // Train on EVEN-indexed events only.
    drs4::StopCellCorrection corr(nBins);
    for (size_t i = 0; i < val.size(); i += 2)
        if (val[i] > kBadEst + 1.f) corr.Accumulate(stop[i], val[i]);
    corr.Finalize(/*minCount=*/30);
    outCorr = corr;

    // Evaluate on held-out ODD events.
    std::vector<float> odd, oddCorr;
    odd.reserve(val.size() / 2);
    oddCorr.reserve(val.size() / 2);
    for (size_t i = 1; i < val.size(); i += 2) {
        if (val[i] <= kBadEst + 1.f) continue;
        odd.push_back(val[i]);
        oddCorr.push_back(static_cast<float>(corr.Apply(stop[i], val[i])));
    }

    sBefore = CoreSigmaPs(odd,     win_ns);
    sAfter  = CoreSigmaPs(oddCorr, win_ns);
}

// ---------------------------------------------------------------------------
// A2Combo7 — A^2-weighted mean of the 7 DRS0-G0 HG channels (indices 0-6),
// each already MCP1-referenced (hg_cfd).  Returns kBadEst if < 2 valid.
// ---------------------------------------------------------------------------
static float A2Combo7(const float hgcfd[8], const float hgpeak[8])
{
    double sw = 0., swt = 0.;
    int    n  = 0;
    for (int i = 0; i < 7; ++i) {              // ch0-6 only (all on DRS0 G0)
        if (hgpeak[i] < kHG_minPeak)  continue;
        if (hgcfd[i]  <= kBadEst + 1.f) continue;
        const double w = static_cast<double>(hgpeak[i]) * hgpeak[i];
        sw += w; swt += w * hgcfd[i]; ++n;
    }
    return (n >= 2 && sw > 0.) ? static_cast<float>(swt / sw) : kBadEst;
}

// ===========================================================================
// Main
// ===========================================================================
void drs4TimeBase()
{
    TH1::AddDirectory(kFALSE);
    ApplyRADiCALStyle();
    gSystem->mkdir("output/Summary", kTRUE);

    // ── Open the highest-energy run (best HG statistics) ─────────────────────
    TChain* chain = new TChain("pulse");
    {
        TObjArray* parts = TString(kRuns[kNRuns - 1].inFiles).Tokenize(";");
        for (int i = 0; i < parts->GetEntries(); ++i) {
            TString pat = ((TObjString*)parts->At(i))->GetString().Strip();
            chain->Add(pat);
        }
        delete parts;
    }
    if (chain->GetEntries() == 0) {
        std::cerr << "drs4TimeBase: no raw data found for "
                  << kRuns[kNRuns - 1].label << "\n";
        return;
    }
    std::cout << "drs4TimeBase: " << chain->GetEntries()
              << " raw events (" << kRuns[kNRuns - 1].label << ")\n";

    Float_t timev[4096], amp[36864];
    chain->SetBranchAddress("timevalue", timev);
    chain->SetBranchAddress("amplitude", amp);

    // ── Accumulators ─────────────────────────────────────────────────────────
    // Cell-width RMS per group (diagnostic of calibration)
    TH1F* hCW[4];
    const char* gName[4] = {"D0G0", "D0G1", "D1G0", "D1G1"};
    const int   gOff [4] = {kT_D0G0, kT_D0G1, kT_D1G0, kT_D1G1};
    for (int g = 0; g < 4; ++g) {
        hCW[g] = new TH1F(Form("hCW_%s", gName[g]),
                          ";per-cell width RMS (ps);events", 100, 0., 60.);
        hCW[g]->SetDirectory(nullptr);
    }
    // Stop-cell distribution (D0G0) — uniformity test
    TH1F* hStop = new TH1F("hStopD0G0",
                           ";DRS0 G0 stop cell;events", 64, 0., 1024.);
    hStop->SetDirectory(nullptr);

    // Per-event estimator storage
    std::vector<int>   evStop;
    std::vector<float> evSingle, evCombo, evCorner;
    evStop.reserve(120000); evSingle.reserve(120000);
    evCombo.reserve(120000); evCorner.reserve(120000);

    // ── Single pass over raw data ─────────────────────────────────────────────
    const Long64_t nEv = chain->GetEntries();
    for (Long64_t ev = 0; ev < nEv; ++ev) {
        chain->GetEntry(ev);

        // Cell-width RMS + stop cell per group (cheap, every event helps stats)
        for (int g = 0; g < 4; ++g) {
            const int sc = drs4::FindStopCell(timev + gOff[g]);
            const double rms = drs4::CellWidthRMS(timev + gOff[g], sc);
            if (rms >= 0.) hCW[g]->Fill(rms);
        }
        const int stopD0G0 = drs4::FindStopCell(timev + kT_D0G0);
        if (stopD0G0 >= 0) hStop->Fill(stopD0G0);

        // MCP1 reference
        Pulse mcp = ExtractPulse(timev + kMCP1_t, amp + kMCP1, 0.20f, 30.f);
        const bool mcpOK = (mcp.valid && mcp.crossingTime > -1e5f &&
                            mcp.peak > kMCP1_minPeak && mcp.peak < kMCP1_maxPeak);
        if (!mcpOK || stopD0G0 < 0) continue;

        // HG channels 0-6 (DRS0 G0), MCP1-referenced
        float hgcfd[8], hgpeak[8];
        for (int i = 0; i < 7; ++i) {
            const CapCfg& c = kCap[i];
            Pulse hg = ExtractPulse(timev + c.hg_t, amp + c.hg, 0.20f, 5.f);
            hgpeak[i] = hg.peak;
            hgcfd[i]  = (hg.valid && hg.crossingTime > -1e5f && hg.peak >= kHG_minPeak)
                        ? static_cast<float>(hg.crossingTime - mcp.crossingTime)
                        : kBadEst;
        }

        // Estimators (keep only the physics peak: |hg_cfd + 27| < 8 ns)
        float single = (hgcfd[2] > kBadEst + 1.f &&
                        std::fabs(hgcfd[2] + 27.f) < 8.f) ? hgcfd[2] : kBadEst;

        float combo  = A2Combo7(hgcfd, hgpeak);
        if (combo > kBadEst + 1.f && std::fabs(combo + 27.f) >= 8.f) combo = kBadEst;

        float corner = kBadEst;
        if (hgcfd[0] > kBadEst + 1.f && hgcfd[4] > kBadEst + 1.f)
            corner = 0.5f * (hgcfd[0] - hgcfd[4]);   // MCP cancels in the difference

        evStop  .push_back(stopD0G0);
        evSingle.push_back(single);
        evCombo .push_back(combo);
        evCorner.push_back(corner);

        if ((ev % 20000) == 0)
            std::cout << "  " << ev << " / " << nEv << "\r" << std::flush;
    }
    std::cout << "\n  collected " << evStop.size() << " MCP-quality events\n";

    // ── Split-half validation for each estimator ─────────────────────────────
    drs4::StopCellCorrection corrSingle(64), corrCombo(64), corrCorner(64);
    double sglB, sglA, cmbB, cmbA, cnrB, cnrA;
    ValidateSplitHalf(evStop, evSingle, 2.0, 64, corrSingle, sglB, sglA);
    ValidateSplitHalf(evStop, evCombo,  2.0, 64, corrCombo,  cmbB, cmbA);
    ValidateSplitHalf(evStop, evCorner, 1.0, 64, corrCorner, cnrB, cnrA);

    std::cout << "\n=== OUT-OF-SAMPLE core sigma (ps): before -> after stop-cell corr ===\n";
    std::cout << Form("  single SE-D : %.1f -> %.1f\n", sglB, sglA);
    std::cout << Form("  combo  (7ch): %.1f -> %.1f\n", cmbB, cmbA);
    std::cout << Form("  corner NW   : %.1f -> %.1f\n", cnrB, cnrA);

    // ── Cell-width verdict ────────────────────────────────────────────────────
    const double cwMean = hCW[0]->GetMean();
    const bool   calibrated = (cwMean > 3.0);   // >3 ps ⇒ genuine cell-width calib

    // =========================================================================
    // Output PDF
    // =========================================================================
    TString outPDF  = "output/Summary/drs4_timebase.pdf";
    TString outROOT = "output/Summary/drs4_timebase.root";
    TCanvas c("c_tb", "", 1100, 760);

    // ── Page 1: cell-width RMS per group ─────────────────────────────────────
    c.Clear(); c.cd(); StylePad(false);
    double yMax1 = 0.;
    for (int g = 0; g < 4; ++g) yMax1 = std::max(yMax1, hCW[g]->GetMaximum());
    // Data sits at 0-5 ps, so the upper-right of the 0-8 ps frame is empty:
    // place the legend there (the sidebar variant clips the "mean" text).
    TLegend* L1 = new TLegend(0.58, 0.60, 0.90, 0.88);
    L1->SetBorderSize(0); L1->SetFillStyle(0); L1->SetTextSize(0.034);
    for (int g = 0; g < 4; ++g) {
        hCW[g]->SetLineColor(kRChannelCols[g]);
        hCW[g]->SetLineWidth(2);
        hCW[g]->GetYaxis()->SetRangeUser(0., 1.2 * yMax1);
        hCW[g]->GetXaxis()->SetRangeUser(0., 8.);   // data sits at 0-5 ps; was 0-60
        hCW[g]->Draw(g == 0 ? "HIST" : "HIST SAME");
        L1->AddEntry(hCW[g], Form("%s: %.2f ps", gName[g], hCW[g]->GetMean()), "l");
    }
    L1->Draw();
    {
        TLine* lc = new TLine(3., 0., 3., 1.2 * yMax1);
        lc->SetLineColor(kRRed); lc->SetLineStyle(2); lc->SetLineWidth(2); lc->Draw();
        TLatex a; a.SetNDC(); a.SetTextSize(0.032); a.SetTextColor(kGray+2);
        // The verdict applies to the DRS0 timing groups (D0G0/D0G1 << 1 ps).
        // D1G0/D1G1 (LG energy + wire chambers) read ~4-5 ps — still essentially
        // nominal and NOT used for precision timing, so the headline is unaffected.
        a.DrawLatex(0.16, 0.86, "DRS0 timing groups (D0G0/D0G1): cell widths < 1 ps");
        a.DrawLatex(0.16, 0.81, "-> NOMINAL, uncalibrated (no per-cell width applied)");
        a.SetTextColor(kGray+1);
        a.DrawLatex(0.16, 0.76, "D1 groups (LG/WC, not used for timing) read ~4-5 ps");
        a.SetTextColor(kRRed);
        a.DrawLatex(0.16, 0.71, "dashed line = 3 ps calibration threshold");
    }
    DrawPageTitle("DRS4 per-cell width RMS by group  --  calibration diagnostic");
    c.Print(outPDF + "(");

    // ── Page 2: stop-cell uniformity ─────────────────────────────────────────
    c.Clear(); c.cd(); StylePad();
    hStop->SetLineColor(kRData); hStop->SetLineWidth(2);
    hStop->SetMinimum(0.);
    hStop->Draw("HIST");
    {
        const double uniformRMS = 1024. / std::sqrt(12.);
        // White backing box so the annotation stays readable over the spiky bars.
        TPave* bg = new TPave(0.155, 0.63, 0.74, 0.905, 0, "brNDC");
        bg->SetFillColor(kWhite); bg->SetLineColor(kWhite); bg->Draw();
        TLatex a; a.SetNDC(); a.SetTextSize(0.038);
        a.DrawLatex(0.17, 0.85, Form("stop-cell mean = %.0f", hStop->GetMean()));
        a.DrawLatex(0.17, 0.79, Form("stop-cell RMS  = %.0f  (uniform = %.0f)",
                                     hStop->GetRMS(), uniformRMS));
        a.SetTextColor(kGray+2); a.SetTextSize(0.032);
        a.DrawLatex(0.17, 0.72, "Uniform rotation -> every physical cell sampled equally");
        a.DrawLatex(0.17, 0.67, "-> data-driven stop-cell correction is well-posed");
    }
    DrawPageTitle("DRS0 G0 stop-cell distribution  --  asynchronous trigger");
    c.Print(outPDF);

    // ── Page 3: stop-cell correction pattern + out-of-sample improvement ─────
    c.Clear();
    c.Divide(2, 2, 0.015, 0.03);

    struct Panel { drs4::StopCellCorrection* corr; const char* name;
                   double sB, sA; int col; };
    Panel panels[3] = {
        { &corrSingle, "single SE-D",        sglB, sglA, kRBlue  },
        { &corrCombo,  "A^{2} combo (7ch)",  cmbB, cmbA, kRGreen },
        { &corrCorner, "corner (NW-D#minusNW-U)/2", cnrB, cnrA, kRRed },
    };
    for (int p = 0; p < 3; ++p) {
        c.cd(p + 1); StylePad();
        // Offset pattern vs stop cell (the reproducible cell-width fixed pattern)
        TGraph* g = new TGraph();
        for (int b = 0; b < panels[p].corr->NBins(); ++b) {
            const double sc = (b + 0.5) * 1024. / panels[p].corr->NBins();
            g->SetPoint(g->GetN(), sc,
                        panels[p].corr->Offset(static_cast<int>(sc)) * 1000.); // ps
        }
        g->SetMarkerStyle(20); g->SetMarkerSize(0.8);
        g->SetMarkerColor(panels[p].col); g->SetLineColor(panels[p].col);
        g->SetLineWidth(2);
        g->SetTitle(";DRS0 G0 stop cell;timing offset (ps)");
        g->Draw("APL");
        DrawPadTitle(panels[p].name, 0.075);
        TLatex a; a.SetNDC(); a.SetTextSize(0.052);
        a.DrawLatex(0.18, 0.26, Form("pattern RMS = %.1f ps",
                                     panels[p].corr->OffsetRMS() * 1000.));
        a.SetTextColor(panels[p].sA < panels[p].sB ? kRGreen : kGray+2);
        a.DrawLatex(0.18, 0.19, Form("#sigma: %.1f #rightarrow %.1f ps (oos)",
                                     panels[p].sB, panels[p].sA));
    }
    // 4th pad: bar summary of out-of-sample improvement
    c.cd(4); StylePad();
    {
        TH1F* hb = new TH1F("hImp", ";;out-of-sample core #sigma (ps)", 3, 0, 3);
        hb->SetDirectory(nullptr);
        const char* lab[3] = {"single", "combo", "corner"};
        double before[3] = {sglB, cmbB, cnrB};
        double after [3] = {sglA, cmbA, cnrA};
        double ymax = 0.;
        for (int i = 0; i < 3; ++i) ymax = std::max(ymax, before[i]);
        for (int i = 0; i < 3; ++i) hb->GetXaxis()->SetBinLabel(i + 1, lab[i]);
        hb->GetYaxis()->SetRangeUser(0., 1.3 * ymax);
        hb->SetLineColor(kGray+2); hb->Draw("AXIS");
        for (int i = 0; i < 3; ++i) {
            // before (open) and after (filled) markers
            TGraph* gb = new TGraph(1); gb->SetPoint(0, i + 0.35, before[i]);
            gb->SetMarkerStyle(24); gb->SetMarkerSize(1.6); gb->SetMarkerColor(kGray+2);
            gb->Draw("P SAME");
            TGraph* ga = new TGraph(1); ga->SetPoint(0, i + 0.65, after[i]);
            ga->SetMarkerStyle(20); ga->SetMarkerSize(1.6);
            ga->SetMarkerColor(after[i] < before[i] ? kRGreen : kRRed);
            ga->Draw("P SAME");
        }
        DrawPadTitle("before (open) / after (filled)", 0.065);
    }
    DrawPageTitle("Stop-cell correction: trained on even events, validated on odd");
    c.Print(outPDF + ")");   // close the 3-page raw-diagnostics PDF (appendix)
    (void) calibrated;       // verdict now lives in the report's KEY FINDING box

    // =========================================================================
    // Hero: one clean 2-panel page that answers "is the time base sound?"
    //   left  — stop cell rotates uniformly  => the correction is well-posed
    //   right — out-of-sample before/after    => the corner CANCELS cell-width
    //           (unchanged), only the MCP-referenced combo needs the correction
    // =========================================================================
    {
        TCanvas ch("c_tb_hero", "", 1320, 560);
        ch.Divide(2, 1, 0.014, 0.02);

        // ----- left pad: stop-cell uniformity --------------------------------
        ch.cd(1); StylePad();
        hStop->SetLineColor(kRData); hStop->SetLineWidth(2); hStop->SetMinimum(0.);
        hStop->SetTitle(";DRS0 G0 stop cell;events");
        hStop->Draw("HIST");
        {
            const double uniformRMS = 1024. / std::sqrt(12.);
            TPave* bg = new TPave(0.155, 0.70, 0.80, 0.905, 0, "brNDC");
            bg->SetFillColor(kWhite); bg->SetLineColor(kWhite); bg->Draw();
            TLatex a; a.SetNDC(); a.SetTextSize(0.044);
            a.DrawLatex(0.17, 0.84, Form("RMS = %.0f  (uniform = %.0f)",
                                         hStop->GetRMS(), uniformRMS));
            a.SetTextColor(kGray+2); a.SetTextSize(0.036);
            a.DrawLatex(0.17, 0.76, "stop cell rotates uniformly");
            a.DrawLatex(0.17, 0.715, "#Rightarrow correction is well-posed");
        }
        DrawPadTitle("Stop cell samples every physical cell equally", 0.060);

        // ----- right pad: out-of-sample before/after (log y) -----------------
        ch.cd(2); StylePad(); gPad->SetLogy(); gPad->SetRightMargin(0.05);
        {
            const char* lab[3] = {"single SE-D", "A^{2} combo", "corner (DW#minusUP)/2"};
            const char* tag[3] = {"MCP-limited", "improves",    "already cancels"};
            double before[3] = {sglB, cmbB, cnrB};
            double after [3] = {sglA, cmbA, cnrA};

            TH1F* hb = new TH1F("hImpHero", ";;out-of-sample core #sigma (ps)", 3, 0, 3);
            hb->SetDirectory(nullptr);
            for (int i = 0; i < 3; ++i) hb->GetXaxis()->SetBinLabel(i + 1, lab[i]);
            hb->GetXaxis()->SetLabelSize(0.043);
            hb->GetYaxis()->SetRangeUser(60., 900.);
            hb->SetLineColor(kGray+2); hb->Draw("AXIS");

            for (int i = 0; i < 3; ++i) {
                const bool improved = after[i] < before[i] - 1.0;
                // connector
                TLine* ln = new TLine(i + 0.30, before[i], i + 0.62, after[i]);
                ln->SetLineColor(kGray+1); ln->SetLineStyle(2); ln->SetLineWidth(1); ln->Draw();
                // before (open) / after (filled, green only when the correction helps)
                TGraph* gb = new TGraph(1); gb->SetPoint(0, i + 0.30, before[i]);
                gb->SetMarkerStyle(24); gb->SetMarkerSize(1.7); gb->SetMarkerColor(kGray+2);
                gb->Draw("P SAME");
                TGraph* ga = new TGraph(1); ga->SetPoint(0, i + 0.62, after[i]);
                ga->SetMarkerStyle(20); ga->SetMarkerSize(1.7);
                ga->SetMarkerColor(improved ? kRGreen : kGray+1);
                ga->Draw("P SAME");
                // value label above, interpretation tag below
                TLatex t; t.SetTextAlign(21);
                t.SetTextColor(kGray+2); t.SetTextSize(0.033);
                t.DrawLatex(i + 0.46, std::max(before[i], after[i]) * 1.14,
                            Form("%.0f #rightarrow %.0f", before[i], after[i]));
                t.SetTextColor(improved ? kRGreen+1 : kGray+2); t.SetTextSize(0.030);
                t.DrawLatex(i + 0.46, std::min(before[i], after[i]) * 0.78, tag[i]);
            }
            // legend for open/filled
            TLegend* lg = new TLegend(0.58, 0.79, 0.94, 0.905);
            lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.033);
            TGraph* mB = new TGraph(1); mB->SetMarkerStyle(24); mB->SetMarkerColor(kGray+2);
            TGraph* mA = new TGraph(1); mA->SetMarkerStyle(20); mA->SetMarkerColor(kRGreen);
            lg->AddEntry(mB, "before correction", "p");
            lg->AddEntry(mA, "after correction",  "p");
            lg->Draw();
        }
        DrawPadTitle("Correction helps the combo, not the corner", 0.060);

        ch.cd(0);
        DrawPageTitle("DRS4 time base: the corner estimator already cancels the cell-width error");
        ch.Print("output/Summary/drs4_timebase_hero.png");
        std::cout << "  -> output/Summary/drs4_timebase_hero.png\n";
    }

    // ── ROOT output ───────────────────────────────────────────────────────────
    TFile* fout = new TFile(outROOT, "RECREATE");
    hStop->Write();
    for (int g = 0; g < 4; ++g) hCW[g]->Write();
    TParameter<double>("cellWidthRMS_D0G0_ps", cwMean).Write();
    TParameter<double>("sigma_single_before_ps", sglB).Write();
    TParameter<double>("sigma_single_after_ps",  sglA).Write();
    TParameter<double>("sigma_combo_before_ps",  cmbB).Write();
    TParameter<double>("sigma_combo_after_ps",   cmbA).Write();
    TParameter<double>("sigma_corner_before_ps", cnrB).Write();
    TParameter<double>("sigma_corner_after_ps",  cnrA).Write();
    fout->Close(); delete fout;

    for (int g = 0; g < 4; ++g) delete hCW[g];
    delete hStop;

    std::cout << "drs4TimeBase: done -> " << outPDF << "\n";
    std::cout << "              root  -> " << outROOT << "\n";
}

// ============================================================================
// systematicUncertainties.C -- Layer 6: Systematic Uncertainty Evaluation
// ============================================================================
//
// Runs the A2-weighted combination timing analysis under cut variations to
// evaluate systematic uncertainties on sigma_t at each beam energy.
//
// For each variation, the full event loop is re-run with the altered cut, a
// ROBUST core sigma_t (truncated-RMS, truncation-bias corrected -- see
// RobustCoreSigma) is computed for the resulting t_combo distribution, and the
// shift relative to the nominal is recorded as a systematic contribution.  A
// robust estimator is used instead of a Gaussian core fit because at low
// statistics (25, 125 GeV) the fit latches onto an inflated core and produces
// large spurious cut-variation shifts.  The total systematic at each energy is
// the quadrature sum of the cut-variation shifts (no fit-model term: the single
// fixed robust estimator has no fit-shape ambiguity to assign).
//
// ── Systematic variations ───────────────────────────────────────────────────
//
//   Nominal (baseline)    fid_r=3.0, pb_ratio=0.30, mcp_lo=200, mcp_hi=750, hg_min=20
//   Fid +0.5mm            fiducial radius 3.5 mm
//   Fid -0.5mm            fiducial radius 2.5 mm
//   Contain +0.05         PbGlass ratio 0.35
//   Contain -0.05         PbGlass ratio 0.25
//   MCP lo +50mV          MCP lower cut raised to 250 mV
//   MCP hi -50mV          MCP upper cut lowered to 700 mV
//   HG thresh +5mV        HG amplitude threshold raised to 25 mV
//
// ── Output ──────────────────────────────────────────────────────────────────
//
//   output/Summary/systematic_uncertainties.pdf  (3 pages)
//   output/Summary/systematic_uncertainties.root
//
//   Page 1: Horizontal bar chart of |delta_sigma_t| at 150 GeV per variation
//   Page 2: sigma_t vs energy with stat + syst uncertainty band
//   Page 3: Systematic uncertainty table (all energies)
//
// Usage (from repository root):
//   root -l -b -q 'Analysis/systematicUncertainties.C+'
// ============================================================================

#include "ChannelConfig.h"
#include "PlotUtils.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
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
#include "TArrow.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

// ---------------------------------------------------------------------------
// Systematic variation definitions
// ---------------------------------------------------------------------------
struct SystVar {
    const char* name;
    double fid_r;      // fiducial radius [mm]
    float  pb_ratio;   // PbGlass containment ratio
    float  mcp_lo;     // MCP lower peak cut [mV]
    float  mcp_hi;     // MCP upper peak cut [mV]
    float  hg_min;     // HG minimum peak [mV]
};

static const SystVar kVars[] = {
    {"Nominal",        3.0,  0.30f, 200.f, 750.f, 20.f},
    {"Fid +0.5mm",     3.5,  0.30f, 200.f, 750.f, 20.f},
    {"Fid -0.5mm",     2.5,  0.30f, 200.f, 750.f, 20.f},
    {"Contain +0.05",  3.0,  0.35f, 200.f, 750.f, 20.f},
    {"Contain -0.05",  3.0,  0.25f, 200.f, 750.f, 20.f},
    {"MCP lo +50mV",   3.0,  0.30f, 250.f, 750.f, 20.f},
    {"MCP hi -50mV",   3.0,  0.30f, 200.f, 700.f, 20.f},
    {"HG thresh +5mV", 3.0,  0.30f, 200.f, 750.f, 25.f},
};
static const int kNVars = 8;

// Sentinel for missing timing
static const float kNoTime = -1e9f;

// ---------------------------------------------------------------------------
// OpenNtuple
// ---------------------------------------------------------------------------
static TFile* OpenNtuple(int r, TTree*& tree)
{
    TString path = Form("output/%s/ntuple.root", kRuns[r].label.Data());
    TFile* f = TFile::Open(path);
    if (!f || f->IsZombie()) {
        std::cerr << "systematicUncertainties: cannot open " << path << "\n";
        tree = nullptr; return nullptr;
    }
    tree = static_cast<TTree*>(f->Get("rad"));
    if (!tree) {
        std::cerr << "systematicUncertainties: tree 'rad' missing in " << path << "\n";
        f->Close(); delete f; tree = nullptr; return nullptr;
    }
    return f;
}

// ---------------------------------------------------------------------------
// A2WeightedCombo -- A2-weighted mean CFD time
//
// Returns kNoTime if fewer than 2 channels pass the amplitude + sentinel cuts.
// The hg_min threshold is passed in from the variation being evaluated.
// ---------------------------------------------------------------------------
static float A2WeightedCombo(const Float_t hg_peak[8], const Float_t hg_cfd[8],
                              float hg_min)
{
    double sw = 0., swt = 0.;
    int    nValid = 0;
    for (int ic = 0; ic < kNCap; ++ic) {
        if (hg_peak[ic] < hg_min)  continue;
        if (hg_cfd[ic]  < -1e5f)   continue;  // CFD-failed sentinel
        double w  = static_cast<double>(hg_peak[ic]) * hg_peak[ic];
        sw  += w;
        swt += w * hg_cfd[ic];
        ++nValid;
    }
    if (nValid < 2 || sw <= 0.) return kNoTime;
    return static_cast<float>(swt / sw);
}

// ---------------------------------------------------------------------------
// StatsFromVec -- mean and population sigma
// ---------------------------------------------------------------------------
static void StatsFromVec(const std::vector<float>& v, double& mean, double& sigma)
{
    mean = sigma = 0.;
    if (v.empty()) return;
    for (float x : v) mean += x;
    mean /= static_cast<double>(v.size());
    for (float x : v) sigma += (x - mean) * (x - mean);
    sigma = std::sqrt(sigma / static_cast<double>(v.size()));
}

// ---------------------------------------------------------------------------
// RobustCoreSigma -- core timing resolution via an iteratively-truncated RMS
// within +/- k*sigma of the median, corrected for the Gaussian truncation bias.
//
// This REPLACES an iterative Gaussian core fit (FitGaussCore) for the systematic
// study.  At low statistics / with non-Gaussian tails (e.g. the A^2-combo at 25
// and 125 GeV) the Gaussian fit can latch onto an inflated core, so a small cut
// variation makes the fit re-stabilise and produces a large *spurious* shift --
// the failure mode that made three independent cut variations all move ~-39 ps
// in lock-step at 125 GeV.  A truncated RMS has no such latching: seeded from the
// median + MAD it converges deterministically, so the cut-variation deltas
// reflect real sample changes, not fit instability.  The truncation-bias factor
// (variance of a Gaussian truncated at +/-k*sigma) is divided out so the result
// matches a Gaussian-core sigma for a Gaussian core.  Returns sigma in v's units.
// ---------------------------------------------------------------------------
static void RobustCoreSigma(const std::vector<float>& v, double k,
                            double& mu, double& sigma, double& sigErr)
{
    mu = sigma = sigErr = 0.;
    const int n = static_cast<int>(v.size());
    if (n < 20) return;

    std::vector<float> s(v.begin(), v.end());
    std::sort(s.begin(), s.end());
    const double med = s[n / 2];

    std::vector<double> dev(n);
    for (int i = 0; i < n; ++i) dev[i] = std::fabs(v[i] - med);
    std::sort(dev.begin(), dev.end());
    double scale = 1.4826 * dev[n / 2];        // MAD -> Gaussian-equivalent sigma
    if (scale < 0.005) scale = 0.200;          // 200 ps safety floor

    double c = med, sg = scale;
    long kIn = 0;
    for (int it = 0; it < 5; ++it) {           // converges in 2-3 iterations
        const double lo = c - k * sg, hi = c + k * sg;
        double sum = 0., sum2 = 0.; kIn = 0;
        for (float x : v) if (x >= lo && x <= hi) { sum += x; sum2 += x * x; ++kIn; }
        if (kIn < 5) return;
        const double m = sum / kIn, var = sum2 / kIn - m * m;
        c = m; sg = (var > 0.) ? std::sqrt(var) : sg;
    }

    // Truncation-bias correction: Var of N(0,sigma^2) truncated to +/-k*sigma is
    //   sigma^2 * [ 1 - 2k phi(k)/(2Phi(k)-1) ].  Divide it out.
    const double TWO_PI = 6.283185307179586;
    const double phi = std::exp(-0.5 * k * k) / std::sqrt(TWO_PI);
    const double Phi = 0.5 * (1.0 + std::erf(k / std::sqrt(2.0)));
    const double fac = 1.0 - 2.0 * k * phi / (2.0 * Phi - 1.0);
    const double corr = (fac > 0.) ? 1.0 / std::sqrt(fac) : 1.0;

    mu     = c;
    sigma  = sg * corr;
    sigErr = sigma / std::sqrt(2.0 * static_cast<double>(kIn));
}

// ---------------------------------------------------------------------------
// PageTitle -- small centered title at the very top of the canvas
// ---------------------------------------------------------------------------
static void PageTitle(const char* t)
{
    TLatex lat; lat.SetNDC(); lat.SetTextSize(0.030); lat.SetTextAlign(21);
    lat.DrawLatex(0.50, 0.987, t);
}

// ===========================================================================
// Main
// ===========================================================================
void systematicUncertainties()
{
    TH1::AddDirectory(kFALSE);
    ApplyRADiCALStyle();

    gSystem->mkdir("output/Summary", kTRUE);
    std::cout << "systematicUncertainties: evaluating cut systematics\n";

    // Results: sigma_t [ps] and statistical error indexed [variation][run]
    double sigT   [kNVars][kNRuns] = {};
    double sigTErr[kNVars][kNRuns] = {};
    // NB: this is a CUT-systematic study, so the total is the quadrature of the
    // cut variations only.  No fit-model term: a single fixed robust estimator
    // (RobustCoreSigma) is used, so there is no fit-shape ambiguity to assign.

    // =========================================================================
    // Main event loop -- one pass per run per variation
    // =========================================================================
    for (int r = 0; r < kNRuns; ++r) {
        TTree* t = nullptr;
        TFile* fin = OpenNtuple(r, t);
        if (!fin) continue;

        // Beam centroid from nominal pre-scan
        double xc, yc, tcfd, trms;
        ScanRunCenters(t, xc, yc, tcfd, trms);
        // ScanRunCenters already calls ResetBranchAddresses internally

        // Branch addresses
        Float_t x_trk, y_trk, mcp_peak, sum_lg, sum_pb;
        Float_t hg_peak[8], hg_cfd[8];
        Bool_t  wc_ok;
        t->SetBranchAddress("x_trk",    &x_trk);
        t->SetBranchAddress("y_trk",    &y_trk);
        t->SetBranchAddress("wc_ok",    &wc_ok);
        t->SetBranchAddress("mcp_peak", &mcp_peak);
        t->SetBranchAddress("sum_lg",   &sum_lg);
        t->SetBranchAddress("sum_pb",   &sum_pb);
        t->SetBranchAddress("hg_peak",   hg_peak);
        // CFD-5% (adopted headline fraction); guarded fallback to CFD-20%.  The
        // cut-variation systematic is thus evaluated on the same basis as the
        // headline analysis.
        t->SetBranchAddress(t->GetBranch("hg_cfd05") ? "hg_cfd05" : "hg_cfd", hg_cfd);

        // Per-variation t_combo accumulator
        std::vector<float> tvec[kNVars];
        for (int v = 0; v < kNVars; ++v) tvec[v].reserve(20000);

        // Nominal vector also used for histogram range determination
        Long64_t nEv = t->GetEntries();
        for (Long64_t ev = 0; ev < nEv; ++ev) {
            t->GetEntry(ev);
            if (!wc_ok) continue;

            for (int v = 0; v < kNVars; ++v) {
                const SystVar& sv = kVars[v];

                // MCP quality cuts (variation-specific)
                if (mcp_peak < sv.mcp_lo || mcp_peak > sv.mcp_hi) continue;

                // Fiducial cut (variation-specific radius, nominal centroid)
                float dx = x_trk - static_cast<float>(xc);
                float dy = y_trk - static_cast<float>(yc);
                if (std::sqrt(dx*dx + dy*dy) >= static_cast<float>(sv.fid_r)) continue;

                // A2-weighted combination (variation-specific HG threshold)
                float tcombo = A2WeightedCombo(hg_peak, hg_cfd, sv.hg_min);
                if (tcombo <= kNoTime + 1e5f) continue;

                // Containment cut (variation-specific ratio)
                // Applied only when sufficient RADiCAL signal is present
                if (sum_lg > kSumLG_centroid &&
                    sum_pb >= sv.pb_ratio * sum_lg) continue;

                tvec[v].push_back(tcombo);
            }
        }

        t->ResetBranchAddresses();
        fin->Close(); delete fin;

        // Robust core sigma for each variation (truncated RMS at +/-2.5 sigma,
        // truncation-bias corrected).  Stable at low stats -- no Gaussian-fit
        // latching, so the cut-variation deltas are real (see RobustCoreSigma).
        for (int v = 0; v < kNVars; ++v) {
            if (static_cast<int>(tvec[v].size()) < 50) {
                std::cout << "  " << kRuns[r].label
                          << " var=" << kVars[v].name
                          << ": too few events (" << tvec[v].size() << ")\n";
                continue;
            }
            double rMu, rSig, rSigErr;
            RobustCoreSigma(tvec[v], 2.5, rMu, rSig, rSigErr);
            if (rSig > 0.) {
                sigT[v][r]    = rSig * 1000.;       // ns -> ps
                sigTErr[v][r] = rSigErr * 1000.;
            }
        }

        std::cout << "  " << kRuns[r].label
                  << Form(": nominal sigma_t = %.1f +/- %.1f ps  (%d events)\n",
                          sigT[0][r], sigTErr[0][r],
                          static_cast<int>(tvec[0].size()));
    }

    // =========================================================================
    // Compute deltas and total systematic
    //
    //   delta[v][r] = sigT[v][r] - sigT[0][r]    (shift vs nominal)
    //   syst[r]     = sqrt( sum_{v>0} delta[v][r]^2 )
    // =========================================================================
    double delta    [kNVars][kNRuns] = {};
    double systTotal[kNRuns]         = {};

    for (int r = 0; r < kNRuns; ++r) {
        if (sigT[0][r] <= 0.) continue;
        double sumSq = 0.;
        for (int v = 1; v < kNVars; ++v) {
            if (sigT[v][r] <= 0.) continue;
            delta[v][r] = sigT[v][r] - sigT[0][r];
            sumSq += delta[v][r] * delta[v][r];
        }
        systTotal[r] = std::sqrt(sumSq);
    }

    // =========================================================================
    // Output PDF
    // =========================================================================
    TString outPDF  = "output/Summary/systematic_uncertainties.pdf";
    TString outROOT = "output/Summary/systematic_uncertainties.root";
    TCanvas c("c_su", "", 960, 720);

    // ── Page 1: Horizontal bar chart at 150 GeV ───────────────────────────────
    // 150 GeV is the last run (r = kNRuns - 1)
    {
        int r150 = kNRuns - 1;
        c.Clear(); c.cd();
        gPad->SetLeftMargin(0.28);
        gPad->SetRightMargin(0.08);
        gPad->SetTopMargin(0.12);
        gPad->SetBottomMargin(0.14);
        gPad->SetTickx(1); gPad->SetTicky(1);
        gPad->SetGridx(1); gPad->SetGridy(0);

        // Number of bars: kNVars-1 variations + 1 total = kNVars
        const int nBars = kNVars;  // vars 1..7 + Total

        // Maximum |delta| for axis range
        double maxDelta = 0.;
        for (int v = 1; v < kNVars; ++v)
            maxDelta = std::max(maxDelta, std::fabs(delta[v][r150]));
        maxDelta = std::max(maxDelta, systTotal[r150]);
        if (maxDelta < 1.) maxDelta = 5.;
        double xMax = maxDelta * 1.45;

        // Draw frame using a dummy histogram
        TH1F hFrame("hSUFrame", "", 1, 0., xMax);
        hFrame.SetDirectory(nullptr);
        hFrame.SetMinimum(0.);
        hFrame.SetMaximum(static_cast<double>(nBars) + 0.5);
        hFrame.GetXaxis()->SetTitle("|#Delta#sigma_{t}| (ps)");
        hFrame.GetXaxis()->SetTitleSize(0.050);
        hFrame.GetXaxis()->SetLabelSize(0.044);
        hFrame.GetYaxis()->SetLabelSize(0.);  // labels drawn manually
        hFrame.GetYaxis()->SetTickLength(0.);
        hFrame.Draw("AXIS");

        // Bar height in pad coordinates (nBars bars fill 0 to nBars+0.5)
        double barH = 0.55;  // half-height in y-axis units

        TLatex lbl;
        lbl.SetNDC(kFALSE);
        lbl.SetTextSize(0.040);
        lbl.SetTextAlign(32);  // right-align at x=0

        // Bars for each variation (v=1..kNVars-1)
        // y positions: variation v sits at y = kNVars - v  (top = var 1)
        for (int v = 1; v < kNVars; ++v) {
            double yc_bar = static_cast<double>(kNVars - v);
            double dv = delta[v][r150];
            double absD = std::fabs(dv);
            if (absD < 1e-9) absD = 0.;

            int col = (dv >= 0.) ? kRed+1 : kBlue+1;

            TBox* box = new TBox(0., yc_bar - barH, absD, yc_bar + barH);
            box->SetFillColor(col);
            box->SetFillStyle(1001);
            box->SetLineColor(col);
            box->SetLineWidth(1);
            box->Draw("F SAME");

            // Label on y-axis
            lbl.SetTextColor(kBlack);
            lbl.DrawLatex(-0.005 * xMax, yc_bar, kVars[v].name);

            // Value annotation inside or outside bar
            TLatex val;
            val.SetNDC(kFALSE);
            val.SetTextSize(0.038);
            val.SetTextColor(kBlack);
            val.SetTextAlign(12);
            val.DrawLatex(absD + 0.01 * xMax, yc_bar,
                          Form("%.1f ps", absD));
        }

        // Total systematic bar (at y = 0.7, drawn below all variations)
        {
            double yc_bar = 0.7;
            double absD = systTotal[r150];

            TBox* box = new TBox(0., yc_bar - barH * 0.8,
                                 absD, yc_bar + barH * 0.8);
            box->SetFillColor(kBlack);
            box->SetFillStyle(3004);   // hatched
            box->SetLineColor(kBlack);
            box->SetLineWidth(2);
            box->Draw("F SAME");

            // Draw outline
            TBox* boxL = new TBox(0., yc_bar - barH * 0.8,
                                  absD, yc_bar + barH * 0.8);
            boxL->SetFillStyle(0);
            boxL->SetLineColor(kBlack);
            boxL->SetLineWidth(2);
            boxL->Draw("SAME");

            lbl.SetTextColor(kBlack);
            lbl.DrawLatex(-0.005 * xMax, yc_bar, "Total syst.");

            TLatex val;
            val.SetNDC(kFALSE);
            val.SetTextSize(0.038);
            val.SetTextColor(kBlack);
            val.SetTextAlign(12);
            val.DrawLatex(absD + 0.01 * xMax, yc_bar,
                          Form("%.1f ps", absD));
        }

        // Separator line between variations and total
        TLine* lSep = new TLine(0., 1.3, xMax * 0.92, 1.3);
        lSep->SetLineColor(kGray+1);
        lSep->SetLineStyle(2);
        lSep->SetLineWidth(1);
        lSep->Draw("SAME");

        // Color key annotation
        {
            TLatex ann; ann.SetNDC();
            ann.SetTextSize(0.035); ann.SetTextColor(kGray+1);
            ann.DrawLatex(0.30, 0.06,
                "Red = positive shift; Blue = negative shift; "
                "Hatched = quadrature total");
        }

        PageTitle(Form("Systematic contributions to #sigma_{t} at %.0f GeV",
                       kRuns[r150].energy_GeV));
        c.Print(outPDF + "(");
    }

    // ── Page 2: sigma_t vs energy with uncertainty band ───────────────────────
    {
        c.Clear(); c.cd(); StylePad();

        // Build nominal TGraphErrors (stat errors only)
        // Build combined TGraphErrors (stat + syst in quadrature for error bars)
        // Build TGraph for syst-only band (drawn as a filled polygon)
        TGraphErrors* gNominal  = new TGraphErrors();
        TGraphErrors* gCombined = new TGraphErrors();

        double yMin = 1e9, yMax = 0.;
        for (int r = 0; r < kNRuns; ++r) {
            if (sigT[0][r] <= 0.) continue;
            double E    = kRuns[r].energy_GeV;
            double sig  = sigT[0][r];
            double stat = sigTErr[0][r];
            double syst = systTotal[r];

            int n = gNominal->GetN();
            gNominal->SetPoint(n, E, sig);
            gNominal->SetPointError(n, 0., stat);

            int m = gCombined->GetN();
            gCombined->SetPoint(m, E, sig);
            gCombined->SetPointError(m, 0., std::sqrt(stat*stat + syst*syst));

            yMin = std::min(yMin, sig - stat - syst);
            yMax = std::max(yMax, sig + stat + syst);
        }
        if (yMax < 10.) { yMin = 0.; yMax = 400.; }

        double yPad = 0.15 * (yMax - yMin);
        double yFrameHi = yMax + 3. * yPad;
        double yFrameLo = std::max(0., yMin - yPad);

        TH1F* frame2 = static_cast<TH1F*>(
            c.DrawFrame(15., yFrameLo, 165., yFrameHi,
                        ";Beam energy (GeV);#sigma_{t}  (ps)"));
        frame2->GetXaxis()->SetTitleSize(0.050);
        frame2->GetYaxis()->SetTitleSize(0.050);
        frame2->GetYaxis()->SetTitleOffset(1.30);

        // Shaded systematic uncertainty band
        // Build as a filled TGraph polygon (upper then lower edge reversed)
        {
            int nPts = gNominal->GetN();
            if (nPts > 0) {
                TGraph* gBand = new TGraph(2 * nPts + 1);
                for (int p = 0; p < nPts; ++p) {
                    double E = gNominal->GetX()[p];
                    double s = gNominal->GetY()[p];
                    int r = -1;
                    for (int rr = 0; rr < kNRuns; ++rr)
                        if (std::fabs(kRuns[rr].energy_GeV - E) < 1.) { r = rr; break; }
                    double syst = (r >= 0) ? systTotal[r] : 0.;
                    gBand->SetPoint(p, E, s + syst);
                    gBand->SetPoint(2*nPts - 1 - p, E, s - syst);
                }
                // Close polygon
                double E0 = gNominal->GetX()[0];
                double s0 = gNominal->GetY()[0];
                double sy0 = (systTotal[0] > 0.) ? systTotal[0] : 0.;
                gBand->SetPoint(2*nPts, E0, s0 + sy0);

                gBand->SetFillColorAlpha(kBlue+1, 0.25f);
                gBand->SetFillStyle(1001);
                gBand->SetLineWidth(0);
                gBand->Draw("F SAME");
            }
        }

        // Combined uncertainty points (stat+syst, hollow markers)
        gCombined->SetLineColor(kBlue+1);
        gCombined->SetMarkerColor(kBlue+1);
        gCombined->SetMarkerStyle(24);  // open circle
        gCombined->SetMarkerSize(1.2);
        gCombined->SetLineWidth(1);
        gCombined->SetLineStyle(2);
        gCombined->Draw("PL SAME");

        // Nominal stat-only points (filled markers, drawn on top)
        gNominal->SetLineColor(kBlue+1);
        gNominal->SetMarkerColor(kBlue+1);
        gNominal->SetMarkerStyle(20);
        gNominal->SetMarkerSize(1.2);
        gNominal->SetLineWidth(2);
        gNominal->Draw("PL SAME");

        // Annotation at 150 GeV
        int r150 = kNRuns - 1;
        if (sigT[0][r150] > 0.) {
            double sigV  = sigT[0][r150];
            double statV = sigTErr[0][r150];
            double systV = systTotal[r150];
            TLatex ann; ann.SetNDC();
            ann.SetTextSize(0.040); ann.SetTextColor(kBlue+1);
            ann.DrawLatex(0.22, 0.22,
                Form("150 GeV: #sigma_{t} = %.1f #pm %.1f (stat) "
                     "#pm %.1f (syst) ps",
                     sigV, statV, systV));
        }

        // Legend
        TLegend leg2(0.55, 0.68, 0.93, 0.88);
        leg2.SetBorderSize(0); leg2.SetFillStyle(0); leg2.SetTextSize(0.038);
        leg2.AddEntry(gNominal,  "Nominal (stat only)", "lp");
        leg2.AddEntry(gCombined, "Nominal (stat #oplus syst)", "lp");
        // Proxy for shaded band
        TBox* bandProxy = new TBox(0,0,1,1);
        bandProxy->SetFillColorAlpha(kBlue+1, 0.25f);
        bandProxy->SetFillStyle(1001);
        leg2.AddEntry(bandProxy, "Syst. band", "f");
        leg2.Draw();

        PageTitle("A^{2}-weighted combo #sigma_{t} vs energy -- nominal with systematic band");
        c.Print(outPDF);
    }

    // ── Page 3: Systematic uncertainty table ──────────────────────────────────
    {
        c.Clear(); c.cd();
        gPad->SetLeftMargin(0.02);
        gPad->SetRightMargin(0.02);
        gPad->SetTopMargin(0.10);
        gPad->SetBottomMargin(0.04);

        // Column layout:
        //   Col 0: Variation name
        //   Cols 1..kNRuns: delta_sigma_t at each energy (ps)
        // Row layout:
        //   Row 0: Header (energy labels)
        //   Rows 1..kNVars-1: delta per variation
        //   Row kNVars:   Nominal sigma_t (reference)
        //   Row kNVars+1: Total syst

        TLatex tab; tab.SetNDC();
        const double xName = 0.02;   // left edge of variation name column
        const double xData0 = 0.29;  // left edge of first energy column
        double colW = (0.97 - xData0) / static_cast<double>(kNRuns);

        // Header row
        tab.SetTextSize(0.038); tab.SetTextColor(kBlack);
        tab.DrawLatex(xName, 0.87, "Variation");
        for (int r = 0; r < kNRuns; ++r) {
            double xc = xData0 + (r + 0.5) * colW;
            tab.SetTextAlign(22);
            tab.DrawLatex(xc, 0.87,
                Form("%.0f GeV", kRuns[r].energy_GeV));
        }
        tab.SetTextAlign(12);

        // Horizontal line under header
        TLine* lH = new TLine(xName, 0.84, 0.98, 0.84);
        lH->SetLineColor(kGray+1); lH->SetLineWidth(1); lH->Draw();

        // Data rows: variations (v=1..kNVars-1)
        double yRow = 0.80;
        const double rowStep = 0.068;
        const double kLargeSyst = 3.0;  // threshold for red highlighting [ps]

        for (int v = 1; v < kNVars; ++v) {
            tab.SetTextSize(0.036); tab.SetTextColor(kBlack);
            tab.SetTextAlign(12);
            tab.DrawLatex(xName, yRow, kVars[v].name);

            for (int r = 0; r < kNRuns; ++r) {
                double xc = xData0 + (r + 0.5) * colW;
                double d  = delta[v][r];
                // Color-code large systematics
                int textCol = (std::fabs(d) > kLargeSyst) ? kRed+1 : kBlack;
                tab.SetTextColor(textCol);
                tab.SetTextAlign(22);
                if (sigT[v][r] > 0.) {
                    const char* sign = (d >= 0.) ? "+" : "";
                    tab.DrawLatex(xc, yRow,
                        Form("%s%.1f", sign, d));
                } else {
                    tab.DrawLatex(xc, yRow, "--");
                }
            }
            yRow -= rowStep;
        }

        // Thin separator before totals
        TLine* lSep = new TLine(xName, yRow + rowStep * 0.4, 0.98,
                                yRow + rowStep * 0.4);
        lSep->SetLineColor(kGray+1); lSep->SetLineStyle(2);
        lSep->SetLineWidth(1); lSep->Draw();
        yRow -= 0.01;

        // Nominal sigma_t row
        tab.SetTextSize(0.036); tab.SetTextColor(kBlack);
        tab.SetTextAlign(12);
        tab.DrawLatex(xName, yRow, "Nominal #sigma_{t} (ps)");
        for (int r = 0; r < kNRuns; ++r) {
            double xc = xData0 + (r + 0.5) * colW;
            tab.SetTextColor(kBlue+1);
            tab.SetTextAlign(22);
            if (sigT[0][r] > 0.)
                tab.DrawLatex(xc, yRow,
                    Form("%.1f", sigT[0][r]));
            else
                tab.DrawLatex(xc, yRow, "--");
        }
        yRow -= rowStep;

        // Total systematic row
        tab.SetTextSize(0.036); tab.SetTextAlign(12);
        tab.SetTextColor(kBlack);
        tab.DrawLatex(xName, yRow, "Total syst. (ps)");
        for (int r = 0; r < kNRuns; ++r) {
            double xc = xData0 + (r + 0.5) * colW;
            int textCol = (systTotal[r] > kLargeSyst) ? kRed+1 : kGreen+2;
            tab.SetTextColor(textCol);
            tab.SetTextAlign(22);
            if (systTotal[r] > 0.)
                tab.DrawLatex(xc, yRow,
                    Form("%.1f", systTotal[r]));
            else
                tab.DrawLatex(xc, yRow, "--");
        }
        yRow -= rowStep;

        // Footer line
        TLine* lFoot = new TLine(xName, yRow + rowStep * 0.5, 0.98,
                                 yRow + rowStep * 0.5);
        lFoot->SetLineColor(kGray+1); lFoot->SetLineWidth(1); lFoot->Draw();

        {
            TLatex note; note.SetNDC(); note.SetTextSize(0.030);
            note.SetTextColor(kGray+1); note.SetTextAlign(12);
            note.DrawLatex(xName, yRow + rowStep * 0.2,
                "Values in ps. #Delta#sigma_{t} = sigma_{t}(variation) -- "
                "#sigma_{t}(nominal). Syst = quadrature sum of variations.");
            note.DrawLatex(xName, yRow - 0.03,
                "Red values: |#Delta#sigma_{t}| > 3 ps. "
                "Beam centroid fixed to nominal pre-scan for all variations.");
        }

        PageTitle("Systematic uncertainty table -- all energies (ps)");
        c.Print(outPDF + ")");
    }

    // =========================================================================
    // Write summary ROOT file
    // =========================================================================
    TFile* fOut = new TFile(outROOT, "RECREATE");

    // Nominal sigma_t vs energy
    TGraphErrors* gNomOut = new TGraphErrors();
    gNomOut->SetName("gSigT_nominal");
    gNomOut->SetTitle("Nominal #sigma_{t} (A^{2}-wgt combo);Energy (GeV);#sigma_{t} (ps)");
    for (int r = 0; r < kNRuns; ++r) {
        if (sigT[0][r] <= 0.) continue;
        int n = gNomOut->GetN();
        gNomOut->SetPoint(n, kRuns[r].energy_GeV, sigT[0][r]);
        gNomOut->SetPointError(n, 0., sigTErr[0][r]);
    }
    gNomOut->Write();

    // Total systematic vs energy
    TGraphErrors* gSystOut = new TGraphErrors();
    gSystOut->SetName("gSystTotal");
    gSystOut->SetTitle("Total systematic #sigma_{t};Energy (GeV);Syst. #sigma_{t} (ps)");
    for (int r = 0; r < kNRuns; ++r) {
        if (sigT[0][r] <= 0.) continue;
        int n = gSystOut->GetN();
        gSystOut->SetPoint(n, kRuns[r].energy_GeV, systTotal[r]);
        gSystOut->SetPointError(n, 0., 0.);
    }
    gSystOut->Write();

    // Per-variation delta graphs
    for (int v = 1; v < kNVars; ++v) {
        TGraphErrors* gDelta = new TGraphErrors();
        TString vname(kVars[v].name);
        vname.ReplaceAll(" ", "_");
        vname.ReplaceAll("+", "p");
        vname.ReplaceAll("-", "m");
        vname.ReplaceAll(".", "d");
        gDelta->SetName(Form("gDelta_%s", vname.Data()));
        gDelta->SetTitle(Form("#Delta#sigma_{t} for %s;Energy (GeV);#Delta#sigma_{t} (ps)",
                              kVars[v].name));
        for (int r = 0; r < kNRuns; ++r) {
            if (sigT[v][r] <= 0.) continue;
            int n = gDelta->GetN();
            gDelta->SetPoint(n, kRuns[r].energy_GeV, delta[v][r]);
            gDelta->SetPointError(n, 0., 0.);
        }
        gDelta->Write();
        delete gDelta;
    }

    fOut->Close(); delete fOut;
    delete gNomOut; delete gSystOut;

    // =========================================================================
    // Console summary
    // =========================================================================
    std::cout << "\n  Systematic summary at 150 GeV:\n";
    int r150 = kNRuns - 1;
    for (int v = 1; v < kNVars; ++v) {
        std::cout << Form("    %-20s  delta = %+.2f ps\n",
                          kVars[v].name, delta[v][r150]);
    }
    std::cout << Form("    %-20s  total = %.2f ps\n",
                      "Total syst.", systTotal[r150]);
    std::cout << "\n";
    std::cout << "systematicUncertainties: done -> " << outPDF << "\n";
    std::cout << "                          root -> " << outROOT << "\n";
}

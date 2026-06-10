// ============================================================================
// PlotUtils.h — shared ROOT plotting utilities for the RADiCAL analysis chain
// ============================================================================
//
// All functions are declared static so this header can be safely included in
// multiple ACLiC-compiled translation units without link-time symbol conflicts.
// Each macro gets its own compiled copy; no ODR violations occur.
//
// Functions:
//   StylePad         — defined in RADiCALStyle.h (re-exported via this header)
//   FitGaussCore     — iterative two-pass Gaussian fit with error output
//   DrawFitOverlay   — draw fitted Gaussian + σ annotation on current pad
//   ScanRunCenters   — pre-pass on a TTree to derive beam centroid + CFD offset
// ============================================================================

#ifndef PLOTUTILS_H
#define PLOTUTILS_H

#include "SelectionCuts.h"
#include "RADiCALStyle.h"   // ApplyRADiCALStyle, StylePad, color constants

#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include "TLatex.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraph.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

// StylePad is defined in RADiCALStyle.h (included above).

// ---------------------------------------------------------------------------
// DrawSuperTitle — canvas super-title that AUTO-FITS its width so it can never
// run off the right edge (the recurring "...HG cli"/"...adopted m" truncation).
// Call after canvas->cd(0). Shrinks the text size until it fits [x0, 1-x0].
// ---------------------------------------------------------------------------
static void DrawSuperTitle(const char* text, float maxSize = 0.020f)
{
    const float x0 = 0.012f, avail = 1.0f - 2.0f * x0;
    TLatex lab; lab.SetNDC(); lab.SetTextFont(62); lab.SetTextAlign(12);
    float s = maxSize;
    for (int it = 0; it < 22; ++it) {
        lab.SetTextSize(s); lab.SetText(x0, 0.985f, text);
        if (lab.GetXsize() <= avail) break;
        s *= 0.93f; if (s < 0.009f) { s = 0.009f; break; }
    }
    lab.SetTextSize(s); lab.SetTextColor(kBlack);
    lab.DrawLatex(x0, 0.985f, text);
}

// ---------------------------------------------------------------------------
// GridWithTitle — divide a canvas into nx*ny pads BELOW a reserved top band that
// holds the (auto-fitting) super-title, so the title can never overprint the
// row-1 panel titles. Returns the grid pad; call grid->cd(i) for panel i.
// ---------------------------------------------------------------------------
static TPad* GridWithTitle(TCanvas* c, int nx, int ny, const char* title,
                           double gx = 0.004, double gy = 0.012,
                           double bandFrac = 0.05, double titleSize = 0.020)
{
    c->cd();
    TPad* g = new TPad(Form("grid_%s", c->GetName()), "", 0.0, 0.0, 1.0, 1.0 - bandFrac);
    g->SetFillStyle(0); g->SetBorderSize(0); g->Draw();
    g->Divide(nx, ny, gx, gy);
    c->cd(0); DrawSuperTitle(title, titleSize);
    return g;
}

// ---------------------------------------------------------------------------
// DrawEnergyLegend — shared "colour -> GeV" key for the multi-energy overlays,
// driven off the energies ACTUALLY drawn (so it can never drift from the data).
// E and cols are parallel; draws short colour-line swatches with "<E> GeV".
// ---------------------------------------------------------------------------
static void DrawEnergyLegend(double x1, double y1, double x2, double y2,
                             const std::vector<double>& E, const std::vector<int>& cols,
                             double txt = 0.030, const char* header = nullptr)
{
    TLegend* lg = new TLegend(x1, y1, x2, y2);
    lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(txt); lg->SetMargin(0.30);
    if (header) lg->SetHeader(header);
    for (size_t i = 0; i < E.size(); ++i) {
        TGraph* g = new TGraph(1); g->SetPoint(0, 0, 0);
        g->SetLineColor(cols[i]); g->SetLineWidth(5); g->SetMarkerColor(cols[i]);
        lg->AddEntry(g, Form("%d GeV", (int)E[i]), "l");
    }
    lg->Draw();
}

// ---------------------------------------------------------------------------
// FitGaussCore — iterative two-pass Gaussian fit to the core of a histogram
//
// First pass starts from the histogram RMS; the second pass re-fits within
// ±nsig of the first-pass result, giving sub-sample precision on σ.
// Returns zeros in all output parameters on failure (too few entries, etc.).
//
// Typical call:
//   double mu, muErr, sig, sigErr;
//   FitGaussCore(h, 2.0, mu, muErr, sig, sigErr);
// ---------------------------------------------------------------------------
static void FitGaussCore(TH1F* h, double nsig,
                          double& mean,  double& meanErr,
                          double& sigma, double& sigmaErr)
{
    mean = meanErr = sigma = sigmaErr = 0.;
    if (!h || h->GetEntries() < 50) return;

    // Use the histogram mean as the seed peak estimate rather than the
    // maximum bin.  GetMaximumBin() can be fooled by a narrow background
    // spike (e.g. MIP muons in an energy histogram) that has a higher
    // count-per-bin than the broad signal peak.  The histogram mean is
    // pulled toward the signal because signal events typically dominate
    // by number, making it a more robust starting point.
    double mu0  = h->GetMean();
    double sig0 = h->GetRMS();
    if (sig0 <= 0.) return;

    // Clamp sig0 so the initial fit range stays within the histogram
    // axis.  A very inflated RMS (from a bimodal distribution) can push
    // the fit range beyond the axis edges, causing ROOT's minimiser to
    // diverge to unphysical values.
    const double xlo = h->GetXaxis()->GetXmin();
    const double xhi = h->GetXaxis()->GetXmax();
    sig0 = std::min(sig0, (xhi - xlo) / (2.0 * nsig + 1.0));

    TF1 f("_pu_fg", "gaus", mu0 - nsig*sig0, mu0 + nsig*sig0);
    f.SetParameters(h->GetMaximum(), mu0, sig0);
    h->Fit(&f, "RQN");

    mu0  = f.GetParameter(1);
    sig0 = std::fabs(f.GetParameter(2));
    if (sig0 <= 0.) return;

    f.SetRange(mu0 - nsig*sig0, mu0 + nsig*sig0);
    h->Fit(&f, "RQN");

    mean     = f.GetParameter(1);
    meanErr  = f.GetParError(1);
    sigma    = std::fabs(f.GetParameter(2));
    sigmaErr = f.GetParError(2);
}

// ---------------------------------------------------------------------------
// DrawFitOverlay — fit a Gaussian to h, draw it on the current pad, and
//                  annotate σ in picoseconds at the given NDC position.
//
// Returns σ [ns], or -1 on failure.  The histogram must already be drawn
// (Draw("HIST")) before calling this function.
//
// Default annotation position (0.62, 0.78) sits below the top frame line
// for pads with SetTopMargin(0.12).
// ---------------------------------------------------------------------------
static double DrawFitOverlay(TH1F* h,
                              int    color = kRed+1,
                              double xNDC  = 0.62,
                              double yNDC  = 0.78)
{
    double mu, muErr, sig, sigErr;
    FitGaussCore(h, 2.0, mu, muErr, sig, sigErr);
    if (sig <= 0.) return -1.;

    TF1 fg("_pu_fo", "gaus", mu - 2.5*sig, mu + 2.5*sig);
    fg.SetParameters(h->GetMaximum(), mu, sig);
    fg.SetLineColor(color);
    fg.SetLineWidth(2);
    fg.DrawCopy("SAME");

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.055);
    lat.DrawLatex(xNDC, yNDC, Form("#sigma = %.0f ps", sig * 1000.));

    return sig;   // [ns]
}

// ---------------------------------------------------------------------------
// ScanRunCenters — derive per-run calibrations via a quick pre-pass
//
// Computes:
//   x_center, y_center   LG-signal-weighted beam centroid [mm] — used as
//                         the centre of the fiducial circle for this run
//   t_cfd_offset          mean of all valid hg_cfd values [ns] — used to
//                         centre the timing histogram window
//   t_cfd_rms             RMS of the timing distribution [ns]
//
// Event selection in this pre-scan (independent of analysis-level cuts):
//   wc_ok AND mcp_peak > kMCP_minPeak_E (loose MCP cut — we just want beam)
//   For centroid:  sum_lg > kSumLG_centroid
//   For timing:    hg_peak[i] > kHG_minPeak_prescan AND hg_cfd[i] valid
//
// Using the stricter kHG_minPeak_prescan in the timing offset estimate
// avoids biasing the mean with noisy low-amplitude crossings.
// ---------------------------------------------------------------------------
static void ScanRunCenters(TTree* t,
                            double& x_center,    double& y_center,
                            double& t_cfd_offset, double& t_cfd_rms)
{
    // Sensible defaults (nominal geometry + typical offset) if scan fails
    x_center     = kCalo_x0;
    y_center     = kCalo_y0;
    t_cfd_offset = -27.;
    t_cfd_rms    =   1.;

    Float_t x_trk, y_trk, sum_lg, mcp_peak;
    Float_t hg_cfd[8], hg_peak[8];
    Bool_t  wc_ok;
    t->SetBranchAddress("x_trk",    &x_trk);
    t->SetBranchAddress("y_trk",    &y_trk);
    t->SetBranchAddress("wc_ok",    &wc_ok);
    t->SetBranchAddress("sum_lg",   &sum_lg);
    t->SetBranchAddress("mcp_peak", &mcp_peak);
    t->SetBranchAddress("hg_cfd",    hg_cfd);
    t->SetBranchAddress("hg_peak",   hg_peak);

    double wx = 0., wy = 0., w = 0.;
    double t_sum = 0., t_sum2 = 0.;
    int    nt = 0;

    Long64_t nEntries = t->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        t->GetEntry(i);
        if (!wc_ok || mcp_peak < kMCP_minPeak_E) continue;

        if (sum_lg > kSumLG_centroid) {
            wx += x_trk * sum_lg;
            wy += y_trk * sum_lg;
            w  += sum_lg;
        }
        for (int j = 0; j < 8; ++j) {
            if (hg_peak[j] > kHG_minPeak_prescan && hg_cfd[j] > -1e5f) {
                t_sum  += hg_cfd[j];
                t_sum2 += hg_cfd[j] * hg_cfd[j];
                ++nt;
            }
        }
    }

    if (w  > 0.) { x_center = wx / w;  y_center = wy / w; }
    if (nt > 0) {
        t_cfd_offset = t_sum / nt;
        double var   = t_sum2 / nt - t_cfd_offset * t_cfd_offset;
        t_cfd_rms    = (var > 0.) ? std::sqrt(var) : 0.1;
    }

    std::cout << "  Pre-scan: centroid ("
              << std::fixed << std::setprecision(1)
              << x_center << ", " << y_center << ") mm"
              << "  CFD offset "
              << std::setprecision(2) << t_cfd_offset
              << " +/- " << t_cfd_rms << " ns\n"
              << std::defaultfloat;

    // Clear all branch addresses set during the pre-scan so that any
    // subsequent SetBranchAddress calls in the calling macro see a clean state.
    // Without this, branch addresses pointing at our now-dead local variables
    // would remain active and cause a segfault on the next GetEntry() call.
    t->ResetBranchAddresses();
}

// ---------------------------------------------------------------------------
// CalibrateCableDelays — per-channel cable delay from Gaussian fit to hg_cfd
//
// For each channel i, fits a Gaussian to the hg_cfd[i] distribution and
// stores the peak position as delay[i] in nanoseconds (each channel fitted
// independently, unlike ScanRunCenters which pools all 8).
//
// MEASURED FINDING (probe, May 2026 — do NOT naively wire into combinations):
//   The inter-channel delay spread is only ~1.1 ns (a ~0.8 ns Down-vs-Up
//   offset + ~0.3 ns within-group), NOT the ~5 ns the pooled t_cfd_rms (6.5 ns)
//   suggested — that pooled RMS is dominated by per-event jitter, not delay.
//   Because each delay is CONSTANT per channel, it only broadens a combination
//   through event-to-event multiplicity/weight variation, which is small.
//   Split-half out-of-sample tests at 150 and 25 GeV:
//       M2 mean-all : 74.5->73.2 ps (150), 104.8->104.6 ps (25)  — negligible
//       M3 A2-wgt   : 72.2->76.6 ps (150), 135.9->135.1 ps (25)  — WORSE at 150
//   Subtracting the UNWEIGHTED mean from an AMPLITUDE-WEIGHTED combo mismatches
//   the weighted centre (time-walk interaction) and injects ~3.5 ps/channel of
//   delay-estimation noise.  The residual amplitude-dependence is the proper
//   job of time-walk correction (timingMethods.C M5/M6/M7), not delay subtraction.
//   => This function is retained for diagnostics / histogram centring only; it
//      is deliberately NOT applied as a combination-alignment correction.
//
// Requirements: hg_peak[i] > kHG_minPeak_prescan AND hg_cfd[i] valid.
// Uses MCP quality cut kMCP1_minPeak (loose, same as ScanRunCenters).
// Needs <TH1F.h>, <TF1.h> (already included in PlotUtils.h).
// ---------------------------------------------------------------------------
static void CalibrateCableDelays(TTree* t, double delay[8])
{
    // Initialise to a sensible default (typical DRS4 absolute time)
    for (int i = 0; i < 8; ++i) delay[i] = -27.0;

    // Per-channel histograms: wide window ±10 ns, 500 bins → 40 ps/bin
    TH1F* hCFD[8];
    for (int i = 0; i < 8; ++i) {
        hCFD[i] = new TH1F(Form("_ccd_h%d", i), "", 500, -37., -17.);
        hCFD[i]->SetDirectory(nullptr);
    }

    Float_t mcp_peak, hg_cfd[8], hg_peak[8];
    Bool_t  wc_ok;
    t->SetBranchAddress("wc_ok",    &wc_ok);
    t->SetBranchAddress("mcp_peak", &mcp_peak);
    t->SetBranchAddress("hg_cfd",    hg_cfd);
    t->SetBranchAddress("hg_peak",   hg_peak);

    Long64_t nEv = t->GetEntries();
    for (Long64_t ev = 0; ev < nEv; ++ev) {
        t->GetEntry(ev);
        if (!wc_ok || mcp_peak < kMCP1_minPeak || mcp_peak > kMCP1_maxPeak) continue;
        for (int i = 0; i < 8; ++i) {
            if (hg_peak[i] > kHG_minPeak_prescan && hg_cfd[i] > -1e5f)
                hCFD[i]->Fill(hg_cfd[i]);
        }
    }
    t->ResetBranchAddresses();

    for (int i = 0; i < 8; ++i) {
        double mu, muE, sig, sigE;
        FitGaussCore(hCFD[i], 2.0, mu, muE, sig, sigE);
        if (sig > 0. && std::fabs(mu + 27.) < 10.) delay[i] = mu;
        delete hCFD[i];
    }
}

// ---------------------------------------------------------------------------
// FitCrystalBall — Crystal Ball fit to a timing or energy distribution
//
// The Crystal Ball function is a Gaussian core with a power-law tail on the
// low side.  It is more robust than a pure Gaussian when hadronic events
// produce a low-side tail (earlier timing / lower energy).
//
// Parameters: [0]=norm, [1]=mean (mu), [2]=sigma, [3]=alpha, [4]=n
// Initial seeds: alpha=1.5, n=3.0 (typical for EM calorimeter timing)
//
// Returns the Gaussian-core sigma.  Returns zeros on failure.
// Interface is identical to FitGaussCore for drop-in replacement.
// ---------------------------------------------------------------------------
static void FitCrystalBall(TH1F* h, double nsig,
                            double& mean,  double& meanErr,
                            double& sigma, double& sigmaErr)
{
    mean = meanErr = sigma = sigmaErr = 0.;
    if (!h || h->GetEntries() < 50) return;

    // Gaussian pre-fit to seed parameters
    double mu0, muE0, sig0, sigE0;
    FitGaussCore(h, nsig, mu0, muE0, sig0, sigE0);
    if (sig0 <= 0.) return;

    // Crystal Ball: Gaussian core + power-law tail on low side
    // ROOT built-in "crystalball": norm * CrystalBall(x; mean, sigma, alpha, n)
    TF1 fcb("_pu_cb", "crystalball",
             mu0 - nsig*sig0, mu0 + nsig*sig0);
    fcb.SetParameters(h->GetMaximum(), mu0, sig0, 1.5, 3.0);
    fcb.SetParLimits(2, 0., 10.*sig0);  // sigma > 0
    fcb.SetParLimits(3, 0.1, 10.);      // alpha > 0
    fcb.SetParLimits(4, 1.01, 50.);     // n > 1 for normalisable tail

    h->Fit(&fcb, "RQN");

    // Second pass: re-seed from first result, tighter window
    fcb.SetRange(mu0 - nsig * std::fabs(fcb.GetParameter(2)),
                 mu0 + nsig * std::fabs(fcb.GetParameter(2)));
    h->Fit(&fcb, "RQN");

    mean     = fcb.GetParameter(1);
    meanErr  = fcb.GetParError(1);
    sigma    = std::fabs(fcb.GetParameter(2));
    sigmaErr = fcb.GetParError(2);
}

#endif // PLOTUTILS_H

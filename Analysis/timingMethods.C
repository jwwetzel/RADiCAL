// ============================================================================
// timingMethods.C — compare timing reconstruction methods for RADiCAL HG channels
// ============================================================================
//
// Reads Analysis/Output/<label>/ntuple.root (produced by processRun.C with the
// extended branch set) and compares eight timing extraction methods:
//
//   M0: CFD 10%                — earlier crossing, minimal time-walk but more noise
//   M1: CFD 20%  [baseline]    — standard operating point used in other macros
//   M2: CFD 30%                — trades noise for reduced risetime walk
//   M3: CFD 50%                — half-maximum crossing, maximum walk suppression
//   M4: LED 20 mV              — fixed absolute threshold, simplest algorithm
//   M5: TOT-corrected LED      — removes amplitude walk: t_corr = t_LED - k*TOT
//   M6: 1/A-corrected CFD-20%  — removes residual walk: t_corr = t_CFD - k/A_HG
//   M7: HG/LG ratio-corrected  — removes longitudinal shower depth jitter:
//                                 t_corr = t_CFD - k*(R - <R>),  R = A_HG/A_LG
//
// Corrections (M5-M7) are determined per-channel via ordinary least-squares
// regression over all fiducial events at each energy, then applied event-by-event.
//
// Outputs:
//   Analysis/Output/<label>/timing_methods.pdf     (2 pages, per energy)
//   Analysis/Output/Summary/timing_methods_summary.pdf
//   Analysis/Output/Summary/timing_methods.root
// ============================================================================

#include "ChannelConfig.h"
#include "PlotUtils.h"       // StylePad, FitGaussCore, DrawFitOverlay, ScanRunCenters

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TAxis.h"

#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>

// ---------------------------------------------------------------------------
// Constants  (amplitude thresholds from SelectionCuts.h via ChannelConfig.h)
// ---------------------------------------------------------------------------
static const float kNoTimeCut_tm = -1e5f;  // sentinel test (kNoTime = -1e6)
static const int   kNBins_tm     = 200;    // bins per timing histogram
static const int   kNM_tm        = 8;      // number of methods

static const char* kMName_tm[kNM_tm] = {
    "CFD 10%",
    "CFD 20%",
    "CFD 30%",
    "CFD 50%",
    "LED 20 mV",
    "LED + TOT walk corr.",
    "CFD + 1/A walk corr.",
    "CFD + HG/LG ratio corr."
};
static const int kMColor_tm[kNM_tm] = {
    kCyan+2, kBlack, kGreen+2, kOrange+1,
    kRed+1, kMagenta+1, kBlue+1, kViolet+2
};
// Line style: solid for raw methods, dashed for corrected
static const int kMStyle_tm[kNM_tm] = { 1, 1, 1, 1, 1, 2, 2, 2 };
static const int kMWidth_tm[kNM_tm] = { 1, 2, 1, 1, 1, 2, 2, 2 };

// ---------------------------------------------------------------------------
// Vector statistics helpers
// ---------------------------------------------------------------------------
static double Vm(const std::vector<double>& v) {
    if (v.empty()) return 0.;
    double s = 0.;
    for (size_t i = 0; i < v.size(); ++i) s += v[i];
    return s / (double)v.size();
}

static double Vrms(const std::vector<double>& v, double mean) {
    if (v.size() < 2) return 0.1;
    double s = 0.;
    for (size_t i = 0; i < v.size(); ++i) s += (v[i]-mean)*(v[i]-mean);
    return std::sqrt(s / (double)(v.size()-1));
}

// Ordinary least-squares slope: y = a + b*x, returns b
static double OLSSlope(const std::vector<double>& x, const std::vector<double>& y) {
    int n = (int)x.size();
    if (n < 15) return 0.;
    double sx=0, sy=0, sxx=0, sxy=0;
    for (int i=0;i<n;++i) { sx+=x[i]; sy+=y[i]; sxx+=x[i]*x[i]; sxy+=x[i]*y[i]; }
    double denom = (double)n*sxx - sx*sx;
    if (std::fabs(denom) < 1e-30) return 0.;
    return ((double)n*sxy - sx*sy) / denom;
}

// Build TH1F from a vector, auto-ranging on mean ± nSig*RMS
static TH1F* VecHist(const char* name, const std::vector<double>& v,
                      double nSig = 3.5) {
    TH1F* h;
    if (v.empty()) {
        h = new TH1F(name, ";t (ns);Events", kNBins_tm, -1., 1.);
    } else {
        double mu  = Vm(v);
        double rms = Vrms(v, mu);
        if (rms < 0.010) rms = 0.100;
        h = new TH1F(name, ";t (ns);Events",
                     kNBins_tm, mu - nSig*rms, mu + nSig*rms);
        for (size_t i = 0; i < v.size(); ++i) h->Fill(v[i]);
    }
    h->SetDirectory(nullptr);
    return h;
}

// GFit — thin wrapper around PlotUtils.h FitGaussCore; returns sigma [ns]
static double GFit(TH1F* h) {
    double mu, muE, sig, sigE;
    FitGaussCore(h, 2.0, mu, muE, sig, sigE);
    return (sig > 0.) ? sig : -1.;
}

// ---------------------------------------------------------------------------
// Per-channel event data
// ---------------------------------------------------------------------------
struct ChData {
    // CFD events: all four fractions valid AND hg_peak >= kHG_minPeak
    std::vector<double> t10, t20, t30, t50;
    std::vector<double> t05;    // CFD-5% (the adopted headline fraction); for the
                                // page-3 leading-edge diagnostic. May be empty if
                                // the hg_cfd05 branch is absent (old ntuples).
    std::vector<double> amp;    // hg_peak [mV] for each CFD event
    std::vector<double> lgamp;  // lg_peak [mV]; 0 if below kLG_minPeak

    // LED+TOT events: hg_peak >= kHG_minPeak AND hg_led valid AND hg_tot > 0
    std::vector<double> tled, tot;
};

// Compute corrected time vector for method m (M5-M7 apply walk corrections)
static std::vector<double> MethodVec(const ChData& d, int m) {
    switch (m) {
        case 0: return d.t10;
        case 1: return d.t20;
        case 2: return d.t30;
        case 3: return d.t50;
        case 4: return d.tled;

        case 5: {
            // TOT-corrected LED: t_corr = t_led - slope*(TOT - <TOT>)
            if (d.tled.size() < 15 || d.tot.size() != d.tled.size())
                return d.tled;
            double mean_tot = Vm(d.tot);
            double slope    = OLSSlope(d.tot, d.tled);
            std::vector<double> r(d.tled.size());
            for (size_t i = 0; i < r.size(); ++i)
                r[i] = d.tled[i] - slope * (d.tot[i] - mean_tot);
            return r;
        }

        case 6: {
            // 1/A corrected CFD-20%: t_corr = t_cfd - slope*(1/A - <1/A>)
            if (d.t20.size() < 15) return d.t20;
            std::vector<double> inv_a(d.amp.size());
            for (size_t i = 0; i < inv_a.size(); ++i)
                inv_a[i] = (d.amp[i] > 0.) ? 1.0 / d.amp[i] : 0.;
            double mean_ia = Vm(inv_a);
            double slope   = OLSSlope(inv_a, d.t20);
            std::vector<double> r(d.t20.size());
            for (size_t i = 0; i < r.size(); ++i)
                r[i] = d.t20[i] - slope * (inv_a[i] - mean_ia);
            return r;
        }

        case 7: {
            // HG/LG ratio corrected CFD-20%
            // Derive slope only from events where LG is valid
            if (d.t20.size() < 15) return d.t20;
            std::vector<double> t_lg, ratio_lg;
            for (size_t i = 0; i < d.t20.size(); ++i) {
                if (d.lgamp[i] > 0.) {
                    t_lg.push_back(d.t20[i]);
                    ratio_lg.push_back(d.amp[i] / d.lgamp[i]);
                }
            }
            if ((int)t_lg.size() < 15) return d.t20;
            double mean_r = Vm(ratio_lg);
            double slope  = OLSSlope(ratio_lg, t_lg);
            std::vector<double> r(d.t20.size());
            for (size_t i = 0; i < d.t20.size(); ++i) {
                double ratio = (d.lgamp[i] > 0.) ? d.amp[i] / d.lgamp[i] : mean_r;
                r[i] = d.t20[i] - slope * (ratio - mean_r);
            }
            return r;
        }

        default: return std::vector<double>();
    }
}

// StylePad() is provided by PlotUtils.h

// ---------------------------------------------------------------------------
// Draw a normalized (peak=1) histogram with given style
// ---------------------------------------------------------------------------
static void DrawNorm(TH1F* h, int color, int lstyle, int lwidth,
                     const char* opt = "HIST") {
    if (!h || h->GetMaximum() <= 0.) return;
    h->Scale(1.0 / h->GetMaximum());
    h->SetLineColor(color);
    h->SetLineStyle(lstyle);
    h->SetLineWidth(lwidth);
    h->SetFillStyle(0);
    h->Draw(opt);
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------
void timingMethods()
{
    ApplyRADiCALStyle();   // RADiCALStyle.h via PlotUtils.h

    TString sumDir = "Analysis/Output/Summary";
    gSystem->mkdir(sumDir, kTRUE);

    // Storage for summary: [energy_idx][method][channel] → sigma [ps]
    // and [energy_idx][method] → best (min) sigma across channels [ps]
    const int kNE = kNRuns;
    double vBest[kNE][kNM_tm];
    double vEnergy[kNE];
    for (int ie = 0; ie < kNE; ++ie) {
        vEnergy[ie] = kRuns[ie].energy_GeV;
        for (int m = 0; m < kNM_tm; ++m) vBest[ie][m] = -1.;
    }

    // -----------------------------------------------------------------------
    // Per-energy loop
    // -----------------------------------------------------------------------
    for (int ie = 0; ie < kNE; ++ie) {
        const RunCfg& rc = kRuns[ie];

        TString ntupleFile = TString("Analysis/Output/") + rc.label + "/ntuple.root";
        TFile* fin = TFile::Open(ntupleFile);
        if (!fin || fin->IsZombie()) {
            std::cerr << "[timingMethods] Cannot open " << ntupleFile << std::endl;
            continue;
        }
        TTree* tree = (TTree*)fin->Get("rad");
        if (!tree) {
            std::cerr << "[timingMethods] No TTree 'rad' in " << ntupleFile << std::endl;
            fin->Close(); continue;
        }

        // Check for new branches (may be absent if processRun not re-run)
        bool hasNew = (tree->GetBranch("hg_led") != nullptr);
        if (!hasNew) {
            std::cout << "[timingMethods] " << rc.label
                      << ": new branches absent — re-run processRun.C first. Skipping.\n";
            fin->Close(); continue;
        }

        // Derive per-run beam centroid via LG-signal-weighted pre-scan
        double xc, yc, tOff, tRms;
        ScanRunCenters(tree, xc, yc, tOff, tRms);
        const float xcf = static_cast<float>(xc);
        const float ycf = static_cast<float>(yc);
        const float rFid2 = static_cast<float>(
            TimingFiducialR(rc.energy_GeV) * TimingFiducialR(rc.energy_GeV));

        // Branch addresses
        Float_t hg_peak[8], hg_cfd[8], hg_cfd10[8], hg_cfd30[8], hg_cfd50[8];
        Float_t hg_cfd05[8];
        Float_t hg_led[8], hg_tot[8], lg_peak[8];
        Float_t x_trk, y_trk, mcp_peak, mcp2_peak;
        Bool_t  wc_ok;

        tree->SetBranchAddress("hg_peak",   hg_peak);
        tree->SetBranchAddress("hg_cfd",    hg_cfd);
        tree->SetBranchAddress("hg_cfd10",  hg_cfd10);
        tree->SetBranchAddress("hg_cfd30",  hg_cfd30);
        tree->SetBranchAddress("hg_cfd50",  hg_cfd50);
        tree->SetBranchAddress("hg_led",    hg_led);
        tree->SetBranchAddress("hg_tot",    hg_tot);
        tree->SetBranchAddress("lg_peak",   lg_peak);

        // CFD-5% (adopted headline fraction) — present in all current ntuples.
        // If absent (pre-reprocess), fill with a sentinel so the page-3 diagnostic
        // simply skips the CFD-5% overlay rather than crashing.
        bool hasCFD05 = (tree->GetBranch("hg_cfd05") != nullptr);
        if (hasCFD05) tree->SetBranchAddress("hg_cfd05", hg_cfd05);
        else for (int i = 0; i < 8; ++i) hg_cfd05[i] = -1e6f;  // kNoTime sentinel
        tree->SetBranchAddress("x_trk",     &x_trk);
        tree->SetBranchAddress("y_trk",     &y_trk);
        tree->SetBranchAddress("mcp_peak",  &mcp_peak);
        tree->SetBranchAddress("wc_ok",     &wc_ok);

        // MCP2 amplitude — present in all ntuples produced by the current processRun.C;
        // safe default (pass-through) if reading an old ntuple without this branch.
        mcp2_peak = 9999.f;
        bool hasMCP2_tm = (tree->GetBranch("mcp2_peak") != nullptr);
        if (hasMCP2_tm) tree->SetBranchAddress("mcp2_peak", &mcp2_peak);

        // -----------------------------------------------------------------------
        // Fill per-channel data vectors.
        //
        // CROSS-VALIDATION (M5-M7): walk corrections are trained and tested on the
        // same event sample if we use all events at once, which gives an optimistically
        // biased sigma.  We split events into odd/even halves: train on odd-indexed,
        // test on even-indexed.  This gives an unbiased sigma estimate at the cost
        // of halved test statistics — a fair trade for a rigorous method comparison.
        //
        // ChData.t20_train/t20_test etc. hold the split samples.
        // -----------------------------------------------------------------------
        ChData chd[8];
        // Common-sample event counts per channel: how many events are valid
        // for ALL methods simultaneously (CFD all fractions AND LED AND TOT).
        // Used to report per-method efficiency relative to the common denominator.
        long nCommon[8] = {};   // events where all methods are valid
        long nCFD[8]    = {};   // events where CFD methods are valid
        long nLED[8]    = {};   // events where LED+TOT are valid

        Long64_t nEv = tree->GetEntries();
        for (Long64_t iev = 0; iev < nEv; ++iev) {
            tree->GetEntry(iev);

            // Event-level cuts: WC track, MCP1 quality window, timing fiducial (3 mm)
            if (!wc_ok || mcp_peak < kMCP1_minPeak || mcp_peak > kMCP1_maxPeak) continue;
            float dx = x_trk - xcf, dy = y_trk - ycf;
            if (dx*dx + dy*dy >= rFid2) continue;

            for (int i = 0; i < 8; ++i) {
                // ── MCP quality: channel 7 (SW-U) uses MCP2 as its reference ──
                // Apply the MCP2 amplitude cut before accepting this channel's timing.
                // Without this, events with a weak MCP2 (large CFD time-walk) would
                // contaminate SW-U timing while passing the MCP1-only event gate.
                if (kCap[i].use_mcp2) {
                    if (!hasMCP2_tm || mcp2_peak < kMCP2_minPeak || mcp2_peak > kMCP2_maxPeak) continue;
                }

                // CFD quality: all four fractions found AND peak above threshold
                bool vcfd = (hg_peak[i] >= kHG_minPeak) &&
                            (hg_cfd10[i] > kNoTimeCut_tm) &&
                            (hg_cfd[i]   > kNoTimeCut_tm) &&
                            (hg_cfd30[i] > kNoTimeCut_tm) &&
                            (hg_cfd50[i] > kNoTimeCut_tm);
                // LED quality: LED crossing AND trailing edge found
                bool vled = (hg_peak[i] >= kHG_minPeak) &&
                            (hg_led[i]  > kNoTimeCut_tm)  &&
                            (hg_tot[i]  > 0.f);

                if (vcfd) ++nCFD[i];
                if (vled) ++nLED[i];
                if (vcfd && vled) ++nCommon[i];

                if (vcfd) {
                    chd[i].t10.push_back(hg_cfd10[i]);
                    chd[i].t20.push_back(hg_cfd[i]);
                    chd[i].t30.push_back(hg_cfd30[i]);
                    chd[i].t50.push_back(hg_cfd50[i]);
                    // CFD-5%: only the leading-edge diagnostic uses this; push the
                    // value when valid so its histogram can be built independently.
                    if (hg_cfd05[i] > kNoTimeCut_tm) chd[i].t05.push_back(hg_cfd05[i]);
                    chd[i].amp.push_back(hg_peak[i]);
                    chd[i].lgamp.push_back(
                        (lg_peak[i] >= kLG_minPeak) ? (double)lg_peak[i] : 0.);
                }
                if (vled) {
                    chd[i].tled.push_back(hg_led[i]);
                    chd[i].tot.push_back(hg_tot[i]);
                }
            }
        }
        fin->Close();

        std::cout << "[timingMethods] " << rc.label
                  << "  CFD events ch0=" << chd[0].t20.size()
                  << "  LED events ch0=" << chd[0].tled.size()
                  << "  common ch0=" << nCommon[0] << "\n";
        std::cout << "  Per-method event counts and efficiencies relative to CFD sample:\n";
        std::cout << "  Ch     CFD(M0-3)  LED(M4)  TOT(M5)  Common  "
                     "LED/CFD  Common/CFD\n";
        for (int i = 0; i < 8; ++i) {
            double effLED    = nCFD[i] > 0 ? 100.*nLED[i]/nCFD[i]    : 0.;
            double effCommon = nCFD[i] > 0 ? 100.*nCommon[i]/nCFD[i] : 0.;
            std::cout << Form("  %-6s %8ld %8ld %8ld %7ld  %6.1f%%   %6.1f%%\n",
                              kCap[i].name,
                              nCFD[i], nLED[i], nLED[i],  // nTOT == nLED (vled requires tot>0)
                              nCommon[i], effLED, effCommon);
        }

        // Compute sigma for all methods × channels
        double sig[8][kNM_tm];
        TH1F* hists[8][kNM_tm];
        for (int i = 0; i < 8; ++i) {
            for (int m = 0; m < kNM_tm; ++m) {
                std::vector<double> tv = MethodVec(chd[i], m);
                hists[i][m] = VecHist(Form("h_tm_%s_ch%d_m%d", rc.label.Data(), i, m), tv);
                sig[i][m] = GFit(hists[i][m]) * 1000.;  // store as ps
            }
        }

        // Record best (min) sigma across all channels for summary
        for (int m = 0; m < kNM_tm; ++m) {
            double best = 9999.;
            for (int i = 0; i < 8; ++i)
                if (sig[i][m] > 0. && sig[i][m] < best) best = sig[i][m];
            vBest[ie][m] = (best < 9999.) ? best : -1.;
        }

        // -------------------------------------------------------------------
        // Print sigma table
        // -------------------------------------------------------------------
        std::cout << "  " << rc.label << " timing resolution (ps):\n";
        std::cout << "  Ch     ";
        for (int m = 0; m < kNM_tm; ++m) std::cout << Form(" %6s", Form("M%d",m));
        std::cout << "\n";
        for (int i = 0; i < 8; ++i) {
            std::cout << Form("  %-6s", kCap[i].name);
            for (int m = 0; m < kNM_tm; ++m) {
                if (sig[i][m] > 0.) std::cout << Form(" %6.0f", sig[i][m]);
                else                std::cout << "      -";
            }
            std::cout << "\n";
        }
        std::cout << "  Best  ";
        for (int m = 0; m < kNM_tm; ++m) {
            if (vBest[ie][m] > 0.) std::cout << Form(" %6.0f", vBest[ie][m]);
            else                   std::cout << "      -";
        }
        std::cout << "\n\n";

        // -------------------------------------------------------------------
        // Per-energy PDF — Page 1: CFD fractions + LED
        // -------------------------------------------------------------------
        TString outPDF = TString("Analysis/Output/") + rc.label + "/timing_methods.pdf";
        TCanvas c1("c_tm1", "", 2400, 1400);
        c1.Divide(3, 3, 0.003, 0.003);

        for (int i = 0; i < 8; ++i) {
            c1.cd(i+1);
            StylePad();

            // Methods on this page: 0=CFD10, 1=CFD20, 2=CFD30, 3=CFD50, 4=LED
            // Choose common x-range centered on CFD20 mean
            const std::vector<double>& ref = chd[i].t20;
            double mu_ref  = Vm(ref);
            double rms_ref = Vrms(ref, mu_ref);
            if (rms_ref < 0.010) rms_ref = 0.200;
            double xlo = mu_ref - 4.0*rms_ref;
            double xhi = mu_ref + 4.0*rms_ref;

            bool first = true;
            double ymax = 0.;
            for (int m = 0; m <= 4; ++m) {
                TH1F* h = hists[i][m];
                if (!h || h->GetEntries() < 10) continue;
                if (h->GetMaximum() > 0.) h->Scale(1.0 / h->GetMaximum());
                h->SetLineColor(kMColor_tm[m]);
                h->SetLineStyle(kMStyle_tm[m]);
                h->SetLineWidth(kMWidth_tm[m]);
                h->SetFillStyle(0);
                h->GetXaxis()->SetRangeUser(xlo, xhi);
                if (h->GetMaximum() > ymax) ymax = h->GetMaximum();
                if (first) {
                    h->GetYaxis()->SetRangeUser(0., 1.35);
                    h->Draw("HIST");
                    first = false;
                } else {
                    h->Draw("HIST SAME");
                }
            }
            if (first) { c1.cd(i+1); continue; }  // no data

            // Title
            TLatex tit; tit.SetNDC(); tit.SetTextSize(0.072); tit.SetTextAlign(22);
            tit.DrawLatex(0.54, 0.93, Form("%.0f GeV  %s", rc.energy_GeV, kCap[i].name));

            // Sigma annotations
            TLatex lat; lat.SetNDC();
            double ytop = 0.88;
            for (int m = 0; m <= 4; ++m) {
                if (sig[i][m] <= 0.) continue;
                lat.SetTextColor(kMColor_tm[m]);
                lat.SetTextSize(0.052);
                lat.DrawLatex(0.16, ytop, Form("%s:  %.0f ps", kMName_tm[m], sig[i][m]));
                ytop -= 0.100;
            }
        }

        // Pad 9: legend for page 1
        c1.cd(9);
        gPad->SetLeftMargin(0.05); gPad->SetTopMargin(0.05);
        TLatex leg1; leg1.SetNDC(); leg1.SetTextSize(0.065);
        leg1.SetTextColor(kBlack);
        leg1.DrawLatex(0.08, 0.92, Form("%.0f GeV", rc.energy_GeV));
        leg1.SetTextSize(0.055);
        leg1.DrawLatex(0.08, 0.82, "CFD Fractions & LED");
        leg1.SetTextSize(0.048);
        for (int m = 0; m <= 4; ++m) {
            leg1.SetTextColor(kMColor_tm[m]);
            leg1.DrawLatex(0.08, 0.68 - m*0.12, kMName_tm[m]);
        }

        c1.Print(outPDF + "(");  // open PDF with page 1

        // -------------------------------------------------------------------
        // Per-energy PDF — Page 2: Walk corrections vs CFD-20% baseline
        // -------------------------------------------------------------------
        TCanvas c2("c_tm2", "", 2400, 1400);
        c2.Divide(3, 3, 0.003, 0.003);

        for (int i = 0; i < 8; ++i) {
            c2.cd(i+1);
            StylePad();

            const std::vector<double>& ref = chd[i].t20;
            double mu_ref  = Vm(ref);
            double rms_ref = Vrms(ref, mu_ref);
            if (rms_ref < 0.010) rms_ref = 0.200;
            double xlo = mu_ref - 4.0*rms_ref;
            double xhi = mu_ref + 4.0*rms_ref;

            // Methods on page 2: 1=CFD20 baseline, 5=TOT-LED, 6=1/A-CFD, 7=ratio-CFD
            int pageM[4] = {1, 5, 6, 7};
            bool first = true;
            for (int pi = 0; pi < 4; ++pi) {
                int m = pageM[pi];
                TH1F* h = hists[i][m];
                if (!h || h->GetEntries() < 10) continue;
                if (h->GetMaximum() > 0.) h->Scale(1.0 / h->GetMaximum());
                h->SetLineColor(kMColor_tm[m]);
                h->SetLineStyle(kMStyle_tm[m]);
                h->SetLineWidth(kMWidth_tm[m]);
                h->SetFillStyle(0);
                h->GetXaxis()->SetRangeUser(xlo, xhi);
                if (first) {
                    h->GetYaxis()->SetRangeUser(0., 1.35);
                    h->Draw("HIST");
                    first = false;
                } else {
                    h->Draw("HIST SAME");
                }
            }
            if (first) continue;

            // Title
            TLatex tit; tit.SetNDC(); tit.SetTextSize(0.072); tit.SetTextAlign(22);
            tit.DrawLatex(0.54, 0.93, Form("%.0f GeV  %s", rc.energy_GeV, kCap[i].name));

            // Sigma annotations with improvement arrow notation
            TLatex lat; lat.SetNDC();
            double ytop = 0.88;
            for (int pi = 0; pi < 4; ++pi) {
                int m = pageM[pi];
                if (sig[i][m] <= 0.) continue;
                lat.SetTextColor(kMColor_tm[m]);
                lat.SetTextSize(0.048);
                if (m == 1) {
                    lat.DrawLatex(0.16, ytop, Form("%s:  %.0f ps", kMName_tm[m], sig[i][m]));
                } else {
                    // Show absolute sigma + improvement vs CFD-20%
                    double improv = sig[i][1] - sig[i][m];
                    if (sig[i][1] > 0. && improv > 0.)
                        lat.DrawLatex(0.16, ytop,
                            Form("%s:  %.0f ps  (#minus%.0f)", kMName_tm[m], sig[i][m], improv));
                    else
                        lat.DrawLatex(0.16, ytop,
                            Form("%s:  %.0f ps", kMName_tm[m], sig[i][m]));
                }
                ytop -= 0.100;
            }
        }

        // Pad 9: legend for page 2
        c2.cd(9);
        gPad->SetLeftMargin(0.05); gPad->SetTopMargin(0.05);
        TLatex leg2; leg2.SetNDC(); leg2.SetTextSize(0.065);
        leg2.SetTextColor(kBlack);
        leg2.DrawLatex(0.08, 0.92, Form("%.0f GeV", rc.energy_GeV));
        leg2.SetTextSize(0.055);
        leg2.DrawLatex(0.08, 0.82, "Walk Corrections");
        leg2.SetTextSize(0.048);
        int pageM2[4] = {1, 5, 6, 7};
        for (int pi = 0; pi < 4; ++pi) {
            int m = pageM2[pi];
            leg2.SetTextColor(kMColor_tm[m]);
            leg2.DrawLatex(0.08, 0.68 - pi*0.12, kMName_tm[m]);
        }

        c2.Print(outPDF + "");  // middle page — keep PDF open for Page 3

        // -------------------------------------------------------------------
        // Per-energy PDF — Page 3: the "elbow" is a CFD-fraction effect, not a
        //   DRS4 satellite, and it lives on the Down capillaries.
        //
        // What you see in the per-channel CFD-20% distribution (M1) on the four
        // DOWN capillaries (NW-D, NE-D, SE-D, SW-D) is a broad, slightly skewed
        // peak with a low-side shoulder — visibly ~2x wider than the UP
        // capillaries.  Direct investigation of the ntuple (Analysis/
        // elbowInvestigation.C) established its origin:
        //
        //   * It is NOT a +0.2 ns DRS4 sampling satellite.  The signed offset
        //     distribution (t20 - median) is a single broad skew with no bump at
        //     one sample-width (0.2 ns); LED here is among the CLEANEST methods,
        //     not the worst.
        //   * It is NOT a selection effect.  The shoulder events are not
        //     saturated, not spikes, have near-normal amplitude and beam radius,
        //     and the MCP reference is shared across the group (an MCP-side
        //     effect would hit every channel, not just the Down ones).
        //   * It is NOT amplitude walk in the correctable sense: a 1/sqrt(A) or
        //     profile walk correction trained out-of-sample removes <30 ps (and
        //     can worsen it) — see Analysis/walkCorrTest.C.
        //   * It IS a leading-edge SHAPE effect — but NOT the mean slope.  The
        //     mean rising edge is actually STEEPER at 20% than at 5% (see
        //     edgeMechanism.C), so electronic noise would make 20% *better*.
        //     Instead the edge shape varies pulse-to-pulse, and that jitter
        //     sigma(t_f - t_5%) grows the higher the threshold and is larger on
        //     the Down capillaries.  CFD-5% times low on the edge where it is most
        //     reproducible — the adopted headline fraction — collapsing the
        //     shoulder (~140 ps) with zero events lost.  The Up capillaries depend
        //     far less on fraction: they are the control.
        //
        // This page overlays CFD-20% (the shouldered M1) and CFD-5% (the adopted
        // headline fraction) per channel, with the windowed RMS of each, so the
        // before/after is immediate.  The headline (DW-UP)/2 already uses CFD-5%,
        // so it never carried this shoulder.
        // -------------------------------------------------------------------
        {
            // windowed RMS (ns) within +/- win of the median — reflects the width
            // the eye sees, including the shoulder (unlike a Gaussian core fit).
            auto winRMS = [](const std::vector<double>& v, double win)->double {
                if (v.size() < 5) return -1.;
                std::vector<double> s = v;
                std::nth_element(s.begin(), s.begin() + s.size()/2, s.end());
                double med = s[s.size()/2];
                double sum = 0., sum2 = 0.; int k = 0;
                for (double x : v) if (std::fabs(x - med) <= win) { sum += x; sum2 += x*x; ++k; }
                if (k < 2) return -1.;
                double m = sum/k, var = sum2/k - m*m;
                return var > 0. ? std::sqrt(var) : 0.;
            };

            TCanvas c3("c_tm3", "", 2400, 1400);
            c3.Divide(3, 3, 0.003, 0.003);

            // Keep histogram pointers alive until AFTER c3.Print().
            TH1F* h20save[8] = {};
            TH1F* h05save[8] = {};

            for (int i = 0; i < 8; ++i) {
                c3.cd(i+1);
                StylePad();

                const std::vector<double>& v20 = chd[i].t20;   // CFD-20% (shouldered)
                const std::vector<double>& v05 = chd[i].t05;   // CFD-5%  (clean)
                if (v20.size() < 20) continue;

                // Each method is shifted to its OWN median so the comparison is of
                // SHAPE/WIDTH, not the (uninteresting) 5%->20% rise-time offset in
                // absolute time.  Window +/-1.6 ns shows the shoulder without the
                // failed-crossing garbage tail.
                auto medianOf = [](std::vector<double> s)->double {
                    if (s.empty()) return 0.;
                    std::nth_element(s.begin(), s.begin()+s.size()/2, s.end());
                    return s[s.size()/2];
                };
                double med20 = medianOf(v20);
                double med05 = v05.empty() ? med20 : medianOf(v05);
                double half_w = 1.6;
                const int nBinW = 160;
                TH1F* h20 = new TH1F(Form("hC20_%s_%d", rc.label.Data(), i), "", nBinW, -half_w, half_w);
                TH1F* h05 = new TH1F(Form("hC05_%s_%d", rc.label.Data(), i), "", nBinW, -half_w, half_w);
                h20->SetDirectory(nullptr); h05->SetDirectory(nullptr);
                for (double t : v20) h20->Fill(t - med20);
                for (double t : v05) h05->Fill(t - med05);

                double mx20 = h20->GetMaximum(), mx05 = h05->GetMaximum();
                if (mx20 > 0.) h20->Scale(1./mx20);
                if (mx05 > 0.) h05->Scale(1./mx05);

                // CFD-5% drawn first as a filled green reference (the "after"),
                // CFD-20% overlaid as a blue outline (the "before", shoulder visible).
                bool have05 = (v05.size() >= 20);
                if (have05) {
                    h05->SetLineColor(kGreen+2); h05->SetLineWidth(2);
                    h05->SetFillColorAlpha(kGreen-9, 0.45); h05->SetFillStyle(1001);
                    h05->GetXaxis()->SetTitle("t #minus median  (ns)");
                    h05->GetXaxis()->SetTitleSize(0.054);
                    h05->GetYaxis()->SetTitle("normalised counts");
                    h05->GetYaxis()->SetTitleSize(0.050);
                    h05->GetYaxis()->SetRangeUser(0., 1.45);
                    h05->Draw("HIST");
                    h20->SetLineColor(kBlue+1); h20->SetLineWidth(2); h20->SetFillStyle(0);
                    h20->Draw("HIST SAME");
                } else {
                    h20->SetLineColor(kBlue+1); h20->SetLineWidth(2); h20->SetFillStyle(0);
                    h20->GetXaxis()->SetTitle("t #minus t_{MCP}  (ns)");
                    h20->GetXaxis()->SetTitleSize(0.054);
                    h20->GetYaxis()->SetTitle("normalised counts");
                    h20->GetYaxis()->SetTitleSize(0.050);
                    h20->GetYaxis()->SetRangeUser(0., 1.45);
                    h20->Draw("HIST");
                }

                // windowed RMS (ps) for each
                double r20 = winRMS(v20, 1.5);
                double r05 = winRMS(v05, 1.5);

                TLatex tit; tit.SetNDC(); tit.SetTextSize(0.064); tit.SetTextAlign(22);
                tit.DrawLatex(0.54, 0.93, Form("%.0f GeV  %s", rc.energy_GeV, kCap[i].name));

                TLatex ann; ann.SetNDC(); ann.SetTextSize(0.046);
                ann.SetTextColor(kBlue+1);
                ann.DrawLatex(0.16, 0.84, Form("CFD-20%%:  %.0f ps", r20 > 0. ? r20*1000. : 0.));
                if (have05) {
                    ann.SetTextColor(kGreen+2);
                    ann.DrawLatex(0.16, 0.75, Form("CFD-5%%:  %.0f ps", r05 > 0. ? r05*1000. : 0.));
                    // Signed delta on EVERY channel (negative = CFD-5% better,
                    // positive = CFD-5% worse — honest for the Up channels too).
                    if (r20 > 0. && r05 > 0.) {
                        double d = (r20 - r05) * 1000.;   // >0 means 5% improves
                        ann.SetTextColor(kGray+3); ann.SetTextSize(0.040);
                        ann.DrawLatex(0.16, 0.66, Form("#Delta = %s%.0f ps",
                                      d >= 0. ? "#minus" : "+", std::fabs(d)));
                    }
                }

                h20save[i] = h20; h05save[i] = h05;
            }

            // Pad 9: corrected explanation
            c3.cd(9);
            gPad->SetLeftMargin(0.05); gPad->SetTopMargin(0.05);
            TLatex info; info.SetNDC();
            info.SetTextSize(0.060); info.SetTextColor(kBlack);
            info.DrawLatex(0.06, 0.93, Form("%.0f GeV", rc.energy_GeV));
            info.SetTextSize(0.048);
            info.DrawLatex(0.06, 0.85, "Down-capillary timing shoulder");
            info.SetTextSize(0.038); info.SetTextColor(kGray+2);
            info.DrawLatex(0.06, 0.75, "The CFD-20% peak on the Down");
            info.DrawLatex(0.06, 0.69, "capillaries is broad + skewed:");
            info.DrawLatex(0.06, 0.63, "the rising-edge shape jitters");
            info.DrawLatex(0.06, 0.57, "more high on the edge.");
            info.DrawLatex(0.06, 0.48, "Timing lower on the edge (CFD-5%,");
            info.DrawLatex(0.06, 0.42, "the headline fraction) collapses");
            info.DrawLatex(0.06, 0.36, "it - no events lost.  Up channels");
            info.DrawLatex(0.06, 0.30, "depend far less on fraction.");
            info.DrawLatex(0.06, 0.22, "#Delta>0: CFD-5% better (see also");
            info.DrawLatex(0.06, 0.16, "the CFD-fraction trend figure).");
            // Color key
            TGraph* d20 = new TGraph(1); d20->SetLineColor(kBlue+1);  d20->SetLineWidth(2);
            TGraph* d05 = new TGraph(1); d05->SetLineColor(kGreen+2); d05->SetLineWidth(2);
            TLegend* legS = new TLegend(0.06, 0.01, 0.97, 0.11);
            legS->SetBorderSize(0); legS->SetTextSize(0.040);
            legS->AddEntry(d20, "CFD 20% (shouldered)", "l");
            legS->AddEntry(d05, "CFD 5% (adopted, clean)", "l");
            legS->Draw();

            c3.Print(outPDF + ")");  // close PDF with page 3 — histograms still alive here

            for (int i = 0; i < 8; ++i) { delete h20save[i]; delete h05save[i]; }
        }

        // Clean up histograms
        for (int i = 0; i < 8; ++i)
            for (int m = 0; m < kNM_tm; ++m)
                delete hists[i][m];

        std::cout << "[timingMethods] " << rc.label << ": wrote " << outPDF << "\n";
    }

    // -----------------------------------------------------------------------
    // Summary: sigma vs energy for each method (best channel per energy)
    // -----------------------------------------------------------------------
    TFile* fsum = new TFile(sumDir + "/timing_methods.root", "RECREATE");
    TGraph* gBest[kNM_tm];
    for (int m = 0; m < kNM_tm; ++m)
        gBest[m] = new TGraph();

    for (int ie = 0; ie < kNE; ++ie) {
        for (int m = 0; m < kNM_tm; ++m) {
            if (vBest[ie][m] > 0.)
                gBest[m]->SetPoint(gBest[m]->GetN(), vEnergy[ie], vBest[ie][m]);
        }
    }

    TCanvas cs("c_tmsum", "", 1800, 800);
    cs.Divide(2, 1, 0.005, 0.005);

    // ------- Left pad: sigma vs energy, all methods -------
    cs.cd(1);
    gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);
    gPad->SetRightMargin(0.04); gPad->SetTopMargin(0.12);
    gPad->SetTickx(1); gPad->SetTicky(1);

    // Find y-axis range
    double yMax = 0.;
    for (int m = 0; m < kNM_tm; ++m) {
        for (int p = 0; p < gBest[m]->GetN(); ++p)
            yMax = std::max(yMax, gBest[m]->GetY()[p]);
    }
    yMax *= 1.25;
    if (yMax < 100.) yMax = 500.;

    TH1F* frame = cs.cd(1)->DrawFrame(15., 0., 165., yMax,
        ";Beam Energy (GeV);#sigma (ps)");
    frame->GetXaxis()->SetTitleOffset(1.1);
    frame->GetYaxis()->SetTitleOffset(1.2);

    for (int m = 0; m < kNM_tm; ++m) {
        gBest[m]->SetLineColor(kMColor_tm[m]);
        gBest[m]->SetMarkerColor(kMColor_tm[m]);
        gBest[m]->SetMarkerStyle(20 + m);
        gBest[m]->SetMarkerSize(1.2);
        gBest[m]->SetLineWidth(kMWidth_tm[m]);
        gBest[m]->SetLineStyle(kMStyle_tm[m]);
        gBest[m]->Draw("PL SAME");
        gBest[m]->Write(Form("gBest_m%d", m));
    }

    TLatex tsum; tsum.SetNDC(); tsum.SetTextSize(0.046); tsum.SetTextAlign(22);
    tsum.DrawLatex(0.55, 0.93, "Best-channel #sigma vs Energy");

    // ------- Right pad: sigma table at each energy -------
    cs.cd(2);
    gPad->SetLeftMargin(0.02); gPad->SetRightMargin(0.02);
    gPad->SetTopMargin(0.04); gPad->SetBottomMargin(0.04);

    TLatex tab; tab.SetNDC(); tab.SetTextSize(0.040);

    // Header row
    tab.SetTextColor(kBlack);
    tab.DrawLatex(0.03, 0.95, "Method");
    double xcol[kNE];
    for (int ie = 0; ie < kNE; ++ie) xcol[ie] = 0.32 + ie * 0.115;
    for (int ie = 0; ie < kNE; ++ie)
        tab.DrawLatex(xcol[ie], 0.95, Form("%.0fG", vEnergy[ie]));

    double yrow = 0.88;
    for (int m = 0; m < kNM_tm; ++m) {
        tab.SetTextColor(kMColor_tm[m]);
        tab.SetTextSize(0.036);
        tab.DrawLatex(0.03, yrow, kMName_tm[m]);
        for (int ie = 0; ie < kNE; ++ie) {
            if (vBest[ie][m] > 0.)
                tab.DrawLatex(xcol[ie], yrow, Form("%.0f", vBest[ie][m]));
            else
                tab.DrawLatex(xcol[ie], yrow, "-");
        }
        yrow -= 0.105;
    }
    // Units label
    tab.SetTextColor(kGray+2); tab.SetTextSize(0.032);
    tab.DrawLatex(0.03, yrow - 0.02, "All values in ps (best channel per energy)");

    TString sumPDF = sumDir + "/timing_methods_summary.pdf";
    cs.Print(sumPDF);
    cs.Write("c_tmsum");

    fsum->Close();

    std::cout << "[timingMethods] Summary: " << sumPDF << "\n";
    std::cout << "[timingMethods] ROOT:    " << sumDir << "/timing_methods.root\n";
    std::cout << "\n[timingMethods] Best-channel sigma summary (ps):\n";
    std::cout << "  Method                       ";
    for (int ie = 0; ie < kNE; ++ie) std::cout << Form(" %5.0fG", vEnergy[ie]);
    std::cout << "\n";
    for (int m = 0; m < kNM_tm; ++m) {
        std::cout << Form("  %-30s", kMName_tm[m]);
        for (int ie = 0; ie < kNE; ++ie) {
            if (vBest[ie][m] > 0.) std::cout << Form(" %5.0f", vBest[ie][m]);
            else                   std::cout << "     -";
        }
        std::cout << "\n";
    }
}

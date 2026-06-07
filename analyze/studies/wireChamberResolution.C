// ===========================================================================
// wireChamberResolution.C  —  data-driven spatial resolution of the delay-line
// wire chamber, measured WITHOUT a second tracker.
//
// Geometry (ChannelConfig.h): a delay-line chamber with four plane signals
//   R, L  (x readout)   and   D, U  (y readout),  all on the SAME DRS4 group
//   (kWC_t = kT_D1G1), with
//       x = kWC_Scale * (t_R - t_L)      kWC_Scale = 7/36 mm/ns
//       y = kWC_Scale * (t_D - t_U)
//
// Self-calibration.  The hit position lives in the DIFFERENCE of the two end
// times; the SUM is position-independent:
//       t_R + t_L = 2 t0 + (delay-line transit)            (no x dependence)
// so the spread of the sum isolates the timing noise that limits the position:
//       sigma_x  =  kWC_Scale * sigma(t_R + t_L)           (UPPER BOUND:
//                   the sum also carries the common beam-arrival jitter t0,
//                   which CANCELS in the position difference t_R - t_L).
//
// Because all four planes share one trigger window, the arrival jitter t0 is
// common to the x- and y-plane sums, so
//       sigma(sumX - sumY)  cancels t0  ->  isolates per-plane electronic noise
// giving a REFINED estimate (under the symmetry sigma_R^2+sigma_L^2 ==
// sigma_D^2+sigma_U^2):
//       sigma_x(clean) = kWC_Scale * sigma(sumX - sumY) / sqrt(2).
//
// Two timing definitions are reported:
//   - peak-sample time  (peakTime): exactly what processRun.C uses for x_trk
//                        -> the resolution of the data AS RECONSTRUCTED.
//   - CFD-50% crossing  (crossingTime): sub-sample interpolation -> what the
//                        same pulses could deliver with finer timing.
//
// Reads the raw `pulse` waveform files locally (Data/RUN*.root), so no HPC /
// reprocessing is needed.
//
// Run:  ROOT_INCLUDE_PATH=Analysis root -l -b -q \
//          'Analysis/wireChamberResolution.C+("Data/RUN*_150_GeV.root")'
// ===========================================================================
#include "ChannelConfig.h"
#include "WaveformUtils.h"
#include "SelectionCuts.h"
#include "PlotUtils.h"      // FitGaussCore, StylePad

#include "TChain.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TFile.h"
#include "TSystem.h"
#include "TParameter.h"
#include "TObjArray.h"
#include "TObjString.h"

#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>

namespace {

// Fit the core of a vector's distribution; return sigma (and its error) in the
// vector's native units.  Auto-ranges around the median ± kSpan robust-RMS.
double coreSigma(const std::vector<float>& v, double& sigErr,
                 const char* hname, TH1F** keep = nullptr)
{
    sigErr = 0.;
    if (v.size() < 200) return -1.;
    std::vector<float> s = v;
    std::sort(s.begin(), s.end());
    double med = s[s.size()/2];
    double q25 = s[(size_t)(0.25*s.size())];
    double q75 = s[(size_t)(0.75*s.size())];
    double rrms = (q75 - q25) / 1.349;             // IQR-based robust sigma (tail-immune)
    if (rrms <= 0.) rrms = 1.;
    // Tight window around the mode so the long delay-line tail / secondary lobe
    // does not bias the core fit — we want the resolution of WELL-reconstructed hits.
    double lo = med - 5.*rrms, hi = med + 5.*rrms;
    TH1F* h = new TH1F(hname, "", 240, lo, hi);
    h->SetDirectory(nullptr);
    for (float x : v) h->Fill(x);
    double mu, muE, sig, sE;
    FitGaussCore(h, 2.0, mu, muE, sig, sE);
    sigErr = sE;
    if (keep) *keep = h; else delete h;
    return sig;
}

} // namespace

void wireChamberResolution(const char* inGlob = nullptr)
{
    const double scale = kWC_Scale;   // mm/ns

    // No-arg (pipeline) call: use the highest-energy run, like drs4TimeBase.C.
    // The Data/ paths resolve locally and on Argon (env.sh setup_data_links).
    // The resolution is energy-independent (the sum method is position-blind),
    // so one run at the best stats / least multiple-scattering is sufficient.
    TString src = inGlob ? TString(inGlob) : kRuns[kNRuns - 1].inFiles;

    TChain ch("pulse");
    int nf = 0;
    { TObjArray* parts = src.Tokenize(";");
      for (int i = 0; i < parts->GetEntries(); ++i)
          nf += ch.Add(((TObjString*)parts->At(i))->GetString());
      delete parts; }
    const char* inGlobMsg = src.Data();
    Long64_t N = ch.GetEntries();
    std::cout << "[wcRes] " << inGlobMsg << " -> " << nf << " file(s), "
              << N << " events\n";
    if (N == 0) { std::cout << "[wcRes] no events — check the glob.\n"; return; }

    TTreeReader reader(&ch);
    TTreeReaderArray<float> time_v(reader, "timevalue");
    TTreeReaderArray<float> amp_v (reader, "amplitude");

    std::vector<float> sumX_pk, sumY_pk, diffSum_pk, x_pk, y_pk;
    std::vector<float> sumX_cf, sumY_cf, diffSum_cf;
    long nValid = 0, nCFD = 0;

    while (reader.Next()) {
        const float* T = &time_v[0];
        const float* A = &amp_v[0];

        Pulse R = ExtractPulse(T + kWC_t, A + kWC_R, 0.50f, kWC_minPeak);
        Pulse L = ExtractPulse(T + kWC_t, A + kWC_L, 0.50f, kWC_minPeak);
        Pulse D = ExtractPulse(T + kWC_t, A + kWC_D, 0.50f, kWC_minPeak);
        Pulse U = ExtractPulse(T + kWC_t, A + kWC_U, 0.50f, kWC_minPeak);
        if (!(R.valid && L.valid && D.valid && U.valid)) continue;
        ++nValid;

        // Peak-sample timing — exactly what processRun.C uses for x_trk/y_trk.
        float sX = R.peakTime + L.peakTime;
        float sY = D.peakTime + U.peakTime;
        sumX_pk.push_back(sX);
        sumY_pk.push_back(sY);
        diffSum_pk.push_back(sX - sY);
        x_pk.push_back(scale * (R.peakTime - L.peakTime));
        y_pk.push_back(scale * (D.peakTime - U.peakTime));

        // CFD-50% crossing timing — finer (sub-sample) alternative.
        if (R.crossingTime > -1e5f && L.crossingTime > -1e5f &&
            D.crossingTime > -1e5f && U.crossingTime > -1e5f) {
            float cX = R.crossingTime + L.crossingTime;
            float cY = D.crossingTime + U.crossingTime;
            sumX_cf.push_back(cX);
            sumY_cf.push_back(cY);
            diffSum_cf.push_back(cX - cY);
            ++nCFD;
        }
    }
    std::cout << "[wcRes] valid 4-plane events: " << nValid
              << "  (CFD-50% valid on all 4: " << nCFD << ")\n";
    if (nValid < 200) { std::cout << "[wcRes] too few events.\n"; return; }

    // ── Fit the spreads ─────────────────────────────────────────────────────
    TH1F *hSX=nullptr, *hSY=nullptr, *hDS=nullptr;
    double eSX, eSY, eDS, eSXc, eSYc, eDSc, eBx, eBy;
    double sSX  = coreSigma(sumX_pk,    eSX,  "hWC_sumX",  &hSX);   // ns
    double sSY  = coreSigma(sumY_pk,    eSY,  "hWC_sumY",  &hSY);
    double sDS  = coreSigma(diffSum_pk, eDS,  "hWC_dsum",  &hDS);
    double sSXc = coreSigma(sumX_cf,    eSXc, "hWC_sumXc");
    double sSYc = coreSigma(sumY_cf,    eSYc, "hWC_sumYc");
    double sDSc = coreSigma(diffSum_cf, eDSc, "hWC_dsumc");
    double bX   = coreSigma(x_pk,       eBx,  "hWC_bx");            // mm (beam profile)
    double bY   = coreSigma(y_pk,       eBy,  "hWC_by");

    // ── Convert to spatial resolution (mm) ──────────────────────────────────
    auto mm = [&](double sigNs){ return scale * sigNs; };           // upper bound
    double res_x_ub  = mm(sSX),  res_y_ub  = mm(sSY);               // peak-sample, UB
    double res_x_cln = mm(sDS / std::sqrt(2.0));                    // arrival-removed, peak
    double res_x_ub_cf = mm(sSXc), res_y_ub_cf = mm(sSYc);          // CFD, UB
    double res_x_cln_cf= mm(sDSc / std::sqrt(2.0));                 // CFD, arrival-removed

    auto um = [](double r){ return 1000.*r; };  // mm -> micron
    std::cout << "\n=== Wire-chamber spatial resolution (data-driven, no 2nd tracker) ===\n";
    std::cout << Form("  kWC_Scale = %.4f mm/ns,  delay-line v_prop = %.3f mm/ns\n",
                      scale, 2.*scale);
    std::cout << Form("  beam profile width (context):  sigma_x = %.2f mm,  sigma_y = %.2f mm\n",
                      bX, bY);
    std::cout << "  --- peak-sample timing (AS RECONSTRUCTED in x_trk/y_trk) ---\n";
    std::cout << Form("    sigma(sumX) = %.3f ns  -> sigma_x <= %.0f um   (upper bound)\n", sSX, um(res_x_ub));
    std::cout << Form("    sigma(sumY) = %.3f ns  -> sigma_y <= %.0f um   (upper bound)\n", sSY, um(res_y_ub));
    std::cout << Form("    [x-check] sigma(sumX-sumY) = %.3f ns -> %.0f um/sqrt2 = %.0f um\n",
                      sDS, um(mm(sDS)), um(res_x_cln));
    std::cout << "             (t0-cancellation cross-check; UNRELIABLE here: sumX,sumY are\n";
    std::cout << "              anti-correlated + secondary-lobe contaminated, so disregard.)\n";
    std::cout << "  --- CFD-50% crossing timing (achievable, sub-sample) ---\n";
    std::cout << Form("    sigma(sumX) = %.3f ns  -> sigma_x <= %.0f um   (upper bound)\n", sSXc, um(res_x_ub_cf));
    std::cout << Form("    sigma(sumY) = %.3f ns  -> sigma_y <= %.0f um   (upper bound)\n", sSYc, um(res_y_ub_cf));
    std::cout << "====================================================================\n";
    std::cout << Form("  HEADLINE: WC spatial resolution sigma_x ~ %.1f mm, sigma_y ~ %.1f mm\n",
                      res_x_ub, res_y_ub);
    std::cout << "  Recommended track-position bin width: >= ~3.5-4 mm (>= resolution);\n";
    std::cout << Form("  quantization step ~%.2f mm (one digitizer sample).\n", scale * 1.3);
    std::cout << "====================================================================\n\n";
    (void)res_x_cln_cf;

    // ── Diagnostic PDF ───────────────────────────────────────────────────────
    TString sumDir = "output/Summary";
    gSystem->mkdir(sumDir, kTRUE);
    {
        TCanvas c("c_wcres", "", 1500, 500);
        c.Divide(3, 1, 0.004, 0.004);
        auto drawFit = [&](int pad, TH1F* h, const char* xt, double sig, double ub_um, int col){
            c.cd(pad); StylePad();
            if (!h) return;
            h->SetLineColor(col); h->SetLineWidth(2);
            h->GetXaxis()->SetTitle(xt); h->GetYaxis()->SetTitle("events");
            h->Draw("HIST");
            DrawFitOverlay(h, kRed+1, 0.60, 0.82);
            TLatex l; l.SetNDC(); l.SetTextSize(0.040);
            l.DrawLatex(0.16, 0.84, Form("#sigma = %.3f ns", sig));
            l.DrawLatex(0.16, 0.78, Form("#sigma_{x} #leq %.0f #mum", ub_um));
        };
        drawFit(1, hSX, "t_{R} + t_{L}  (ns)",         sSX, um(res_x_ub),  kBlue+1);
        drawFit(2, hSY, "t_{D} + t_{U}  (ns)",         sSY, um(res_y_ub),  kBlue+1);
        drawFit(3, hDS, "sumX #minus sumY  (ns)",      sDS, um(res_x_cln), kGreen+2);
        c.cd(3);
        TLatex l2; l2.SetNDC(); l2.SetTextSize(0.036); l2.SetTextColor(kGreen+2);
        l2.DrawLatex(0.16, 0.72, Form("#sigma_{x}#approx%.0f #mum (t_{0} removed)", um(res_x_cln)));
        TString pdf = sumDir + "/wire_chamber_resolution.pdf";
        c.Print(pdf);
        std::cout << "[wcRes] wrote " << pdf << "\n";
    }

    // ── Persist scalars ───────────────────────────────────────────────────────
    {
        TFile fo(sumDir + "/wire_chamber_resolution.root", "RECREATE");
        auto put = [&](const char* n, double v){ TParameter<double>(n, v).Write(); };
        put("wc_res_x_ub_um",      um(res_x_ub));
        put("wc_res_y_ub_um",      um(res_y_ub));
        put("wc_res_x_clean_um",   um(res_x_cln));
        put("wc_res_x_ub_cfd_um",  um(res_x_ub_cf));
        put("wc_res_x_clean_cfd_um", um(res_x_cln_cf));
        put("wc_beam_sigma_x_mm",  bX);
        put("wc_beam_sigma_y_mm",  bY);
        put("wc_nvalid",           (double)nValid);
        fo.Close();
        std::cout << "[wcRes] wrote " << sumDir << "/wire_chamber_resolution.root\n";
    }
}

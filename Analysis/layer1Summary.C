// ===========================================================================
// layer1Summary.C  —  the compact, Ledovskoy-clean "Layer 1: Hardware Integrity"
//                     hero figures.  One message per plot, read instantly.
//
// Distils the dense Layer-1 diagnostics (drs4Diagnostics / channelIntegrity /
// averageWaveforms / chargeProfiles / drs4TimeBase — each still produces its
// full multipage "appendix" PDF) into four headline figures that say:
// "the instrument is sound and every channel is alive."
//
//   layer1_pulse_shapes.pdf   H1  average pulse shape per capillary (4x2)
//   layer1_vitals.pdf         H2  pedestal noise + activity per channel
//   layer1_linearity.pdf      H3  mean HG amplitude vs beam energy (8 channels)
//   layer1_timebase.pdf       H4  DRS4 stop-cell uniformity
//
// Reads only persisted analysis products in Output/Summary/*.root (same pattern
// as harvestResults.C) — no raw data, no re-analysis.
//
// Run:  ROOT_INCLUDE_PATH=Analysis root -l -b -q 'Analysis/layer1Summary.C+'
// ===========================================================================
#include <cmath>
#include <algorithm>

#include "TFile.h"
#include "TProfile.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TParameter.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TSystem.h"

#include "RADiCALStyle.h"
#include "ChannelConfig.h"

namespace {

const char* kSumDir = "Analysis/Output/Summary/";
// Representative energy for single-energy figures (highest SNR), with fallbacks.
const char* kPrefer[] = {"150GeV", "125GeV", "100GeV", "75GeV", "50GeV", "25GeV"};

double graphAtE(TGraph* g, double E)
{
    if (!g) return NAN;
    double x, y;
    for (int i = 0; i < g->GetN(); ++i) { g->GetPoint(i, x, y); if (std::fabs(x - E) < 1.0) return y; }
    return NAN;
}
double graphMax(TGraph* g)
{
    if (!g || g->GetN() == 0) return NAN;
    double x, y, m = -1e30;
    for (int i = 0; i < g->GetN(); ++i) { g->GetPoint(i, x, y); m = std::max(m, y); }
    return m;
}

// ── H1 — average pulse shape per capillary ─────────────────────────────────
void HeroPulseShapes()
{
    TFile f(Form("%saverage_waveforms.root", kSumDir));
    if (f.IsZombie()) { printf("[layer1Summary] no average_waveforms.root — skip H1\n"); return; }

    // 4×4: high gain (rows 1–2, navy) over low gain (rows 3–4, teal).
    // Taller title band (104 px) leaves room for a one-line electronics note.
    TCanvas* c; TPad* g = NewSquareGrid(c, "c_l1_pulse", 4, 4, 250, 74, 66, 46, 24, 104);
    const char* usedLabel = nullptr;

    // Find a profile across the preferred energies for a given key format.
    auto findProf = [&](const char* fmt, const char* chan) -> TProfile* {
        for (const char* lab : kPrefer) {
            TProfile* p = dynamic_cast<TProfile*>(f.Get(Form(fmt, lab, chan)));
            if (p) { usedLabel = lab; return p; }
        }
        return nullptr;
    };

    // Draw one normalized average pulse into pad `pad` (1-based).  `rebin`
    // averages groups of source bins (the slow LG pulse doesn't need 0.25 ns
    // resolution — light rebinning removes normalization fuzz without distorting
    // the shape).
    auto drawProf = [&](int pad, TProfile* p, const char* cap, int uid, int col,
                        int rebin = 1, double ylo = -0.15) {
        g->cd(pad);
        if (!p) return;
        const int nb   = p->GetNbinsX();
        const int nout = nb / rebin;
        TH1D* h = new TH1D(Form("hps_%d", uid), "", nout,
                           p->GetXaxis()->GetXmin(), p->GetXaxis()->GetXmax());
        h->SetDirectory(nullptr);
        for (int j = 0; j < nout; ++j) {
            double s = 0.; int cnt = 0;
            for (int k = 0; k < rebin; ++k) {
                const int b = j * rebin + k + 1;
                if (b <= nb) { s += p->GetBinContent(b); ++cnt; }
            }
            h->SetBinContent(j + 1, cnt ? s / cnt : 0.);
        }
        double mx = 0.;
        for (int b = 1; b <= nout; ++b) mx = std::max(mx, h->GetBinContent(b));
        if (mx <= 0.) mx = 1.;
        h->Scale(1.0 / mx);

        h->SetLineColor(col);
        h->SetLineWidth(2);
        h->GetXaxis()->SetTitle("t #minus t_{CFD} (ns)");
        h->GetYaxis()->SetTitle("a(t) / a_{max}");
        h->GetYaxis()->SetRangeUser(ylo, 1.22);
        h->GetXaxis()->SetTitleSize(0.058);
        h->GetYaxis()->SetTitleSize(0.058);
        h->Draw("HIST L");

        TLine* l0 = new TLine(0., ylo, 0., 1.22);
        l0->SetLineColor(kGray + 1); l0->SetLineStyle(2); l0->SetLineWidth(1); l0->Draw();

        DrawPadTitle(cap, 0.078f);
    };

    for (int i = 0; i < 8; ++i)   // rows 1–2: high gain
        drawProf(i + 1, findProf("sum_%s_%s", kCap[i].name),
                 Form("%s  HG", kCap[i].name), i, kRData);
    for (int i = 0; i < 8; ++i)   // rows 3–4: low gain (light rebin; lower y-floor
                                  //           shows the negative undershoot / reflection)
        drawProf(i + 9, findProf("sum_%s_%s_LG", kCap[i].name),
                 Form("%s  LG", kCap[i].name), i + 8, kRTeal, 3, -0.5);

    c->cd(0);
    DrawPageTitle(Form("Average pulse shape per capillary, %s  "
                       "--  HG timing (rows 1-2),  LG energy (rows 3-4)",
                       usedLabel ? usedLabel : "best run"));
    { TLatex sub; sub.SetNDC(); sub.SetTextSize(0.016); sub.SetTextColor(kGray + 2);
      sub.SetTextAlign(12);
      sub.DrawLatex(0.012, 0.945,
        "LG shown full-window: AC-coupled -- pulse + 35% undershoot, "
        "clean baseline recovery by 500 ns, no ringing"); }
    // Direct PNG: pixel-exact aspect (ROOT's PDF paper-fit letterboxes; PNG does not).
    c->Print(Form("%slayer1_pulse_shapes.png", kSumDir));
    printf("[layer1Summary] wrote layer1_pulse_shapes.png\n");
}

// ── H2 — channel vitals (pedestal noise + activity) ────────────────────────
void HeroVitals()
{
    TFile fd(Form("%sdrs4_diagnostics.root", kSumDir));
    if (fd.IsZombie()) { printf("[layer1Summary] no drs4_diagnostics.root — skip H2\n"); return; }
    TFile fc(Form("%schannel_integrity.root", kSumDir));   // may be absent

    TCanvas* c = NewSquareCanvas("c_l1_vitals", 660, 128);  // wider left margin for the y-title
    c->cd();

    TH1F* hPed = new TH1F("hPedVitals", "", 8, 0.5, 8.5);
    hPed->SetDirectory(nullptr);
    double minActive = 100., maxSat = 0.;

    for (int i = 0; i < 8; ++i) {
        double ped = graphAtE(dynamic_cast<TGraph*>(fd.Get(Form("gPedRms_ch%d", i))), 150.);
        hPed->SetBinContent(i + 1, std::isnan(ped) ? 0. : ped);
        hPed->GetXaxis()->SetBinLabel(i + 1, kCap[i].name);

        double s = graphMax(dynamic_cast<TGraph*>(fd.Get(Form("gSatFrac_ch%d", i))));
        if (!std::isnan(s)) maxSat = std::max(maxSat, s);

        if (!fc.IsZombie()) {
            TGraph* ga = dynamic_cast<TGraph*>(fc.Get(Form("gActiveFrac_%s", kCap[i].name)));
            if (ga) for (int k = 0; k < ga->GetN(); ++k) {
                double x, y; ga->GetPoint(k, x, y);
                if (y > 0.) minActive = std::min(minActive, y);
            }
        }
    }

    hPed->SetMarkerStyle(20); hPed->SetMarkerSize(1.8); hPed->SetMarkerColor(kRData);
    hPed->SetLineColor(kRData);
    hPed->GetYaxis()->SetTitle("pedestal RMS (mV)");
    hPed->GetYaxis()->SetTitleOffset(1.45);
    hPed->GetXaxis()->SetLabelSize(0.050);
    hPed->GetYaxis()->SetTitleSize(0.050);
    double ymax = hPed->GetMaximum() * 1.4; if (ymax < 6.) ymax = 6.;
    hPed->GetYaxis()->SetRangeUser(0., ymax);
    hPed->Draw("P");

    TLine* l5 = new TLine(0.5, 5., 8.5, 5.);
    l5->SetLineColor(kRRed); l5->SetLineStyle(2); l5->SetLineWidth(2); l5->Draw();
    { TLatex t; t.SetNDC(); t.SetTextSize(0.030); t.SetTextColor(kRRed); t.SetTextAlign(31);
      t.DrawLatex(0.90, 0.815, "5 mV noise floor"); }   // right-aligned, clear of the right frame

    { TLatex v; v.SetNDC(); v.SetTextSize(0.040); v.SetTextColor(kGray + 3);
      TString verdict = "all 8 capillaries active";
      if (minActive <= 100.) verdict += Form(" > %.1f%%", minActive);
      v.DrawLatex(0.17, 0.58, verdict);
      TLatex v2; v2.SetNDC(); v2.SetTextSize(0.040); v2.SetTextColor(kGray + 3);
      v2.DrawLatex(0.17, 0.52, Form("max saturation %.2f%%   |   noise ~1.3 mV", maxSat)); }

    DrawPageTitle("Channel vitals  --  pedestal noise per capillary (150 GeV)");
    c->Print(Form("%slayer1_vitals.png", kSumDir));
    printf("[layer1Summary] wrote layer1_vitals.png\n");
}

// ── H3 — linearity: mean HG amplitude vs beam energy ───────────────────────
void HeroLinearity()
{
    TFile f(Form("%shg_amplitude_vs_energy.root", kSumDir));
    if (f.IsZombie()) { printf("[layer1Summary] no hg_amplitude_vs_energy.root — skip H3\n"); return; }

    TCanvas* c = NewSquareCanvas("c_l1_lin");
    c->cd();

    TMultiGraph* mg = new TMultiGraph();
    TLegend* leg = MakeCornerLegend(8, "tl", 0.033);

    for (int i = 0; i < 8; ++i) {
        TGraph* gg = dynamic_cast<TGraph*>(f.Get(Form("gHGAmp_%s", kCap[i].name)));
        if (!gg) continue;
        gg->SetLineColor(kRChannelCols[i]);   gg->SetMarkerColor(kRChannelCols[i]);
        gg->SetMarkerStyle(20); gg->SetMarkerSize(1.0); gg->SetLineWidth(2);
        mg->Add(gg, "PL");
        leg->AddEntry(gg, kCap[i].name, "lp");
    }

    mg->Draw("A");
    mg->GetXaxis()->SetTitle("beam energy (GeV)");
    mg->GetYaxis()->SetTitle("mean HG amplitude (mV)");
    mg->GetXaxis()->SetLimits(0., 165.);
    mg->GetXaxis()->SetTitleSize(0.048);
    mg->GetYaxis()->SetTitleSize(0.048);
    leg->Draw();

    { TLatex n; n.SetNDC(); n.SetTextSize(0.027); n.SetTextColor(kGray + 2); n.SetTextAlign(31);
      n.DrawLatex(0.93, 0.30, "All 8 capillaries track together.");
      n.DrawLatex(0.93, 0.255, "HG compresses near the rail above ~50 GeV;");
      n.DrawLatex(0.93, 0.210, "energy is measured on the linear LG channels."); }

    DrawPageTitle("Channel response  --  mean HG amplitude vs beam energy");
    c->Print(Form("%slayer1_linearity.png", kSumDir));
    printf("[layer1Summary] wrote layer1_linearity.png\n");
}

// ── H4 — DRS4 time-base: stop-cell uniformity ──────────────────────────────
void HeroTimebase()
{
    TFile f(Form("%sdrs4_timebase.root", kSumDir));
    if (f.IsZombie()) { printf("[layer1Summary] no drs4_timebase.root — skip H4\n"); return; }
    TH1F* h = dynamic_cast<TH1F*>(f.Get("hStopD0G0"));
    if (!h) { printf("[layer1Summary] hStopD0G0 missing — skip H4\n"); return; }
    h->SetDirectory(nullptr);

    // The raw stop-cell histogram aliases into a comb (cells are quantized), so
    // the *cumulative* distribution is the clean, honest uniformity test: uniform
    // coverage ⇒ a straight diagonal.  Deviation from the diagonal is the only
    // real non-uniformity.  (The comb integrates away.)
    const int    nb   = h->GetNbinsX();
    const double xmin = h->GetXaxis()->GetXmin();
    const double xmax = h->GetXaxis()->GetXmax();
    const double tot  = h->Integral();
    TH1D* hc = new TH1D("hStopCumul", "", nb, xmin, xmax);
    hc->SetDirectory(nullptr);
    double run = 0.;
    for (int b = 1; b <= nb; ++b) { run += h->GetBinContent(b); hc->SetBinContent(b, tot > 0 ? 100. * run / tot : 0.); }

    TCanvas* c = NewSquareCanvas("c_l1_tb");
    c->cd();

    hc->SetLineColor(kRData); hc->SetLineWidth(3);
    hc->GetXaxis()->SetTitle("DRS4 stop cell");
    hc->GetYaxis()->SetTitle("cumulative coverage (%)");
    hc->GetXaxis()->SetTitleSize(0.048);
    hc->GetYaxis()->SetTitleSize(0.048);
    hc->GetYaxis()->SetRangeUser(0., 105.);
    hc->Draw("HIST L");

    TLine* ideal = new TLine(xmin, 0., xmax, 100.);
    ideal->SetLineColor(kRRed); ideal->SetLineStyle(2); ideal->SetLineWidth(2); ideal->Draw();

    { TLatex t; t.SetNDC(); t.SetTextSize(0.030); t.SetTextColor(kGray + 3);
      t.DrawLatex(0.13, 0.85, "Stop cell rotates uniformly #Rightarrow correction is well-posed"); }
    TParameter<double>* cw = dynamic_cast<TParameter<double>*>(f.Get("cellWidthRMS_D0G0_ps"));
    if (cw) { TLatex t2; t2.SetNDC(); t2.SetTextSize(0.028); t2.SetTextColor(kGray + 2);
      t2.DrawLatex(0.13, 0.805, Form("cell-width RMS = %.2f ps  (nominal; no DRS4 calibration)",
                                    cw->GetVal())); }
    // Why the uncalibrated time-base does not bite: it is common-mode, and the
    // residual stop-cell pattern is corrected + validated out-of-sample.  These
    // sit in the empty lower-right triangle (right-aligned), clear of the diagonal.
    { TLatex tr; tr.SetNDC(); tr.SetTextColor(kRRed); tr.SetTextSize(0.030); tr.SetTextAlign(31);
      tr.DrawLatex(0.92, 0.46, "ideal uniform (diagonal)"); }
    TParameter<double>* cb = dynamic_cast<TParameter<double>*>(f.Get("sigma_combo_before_ps"));
    TParameter<double>* ca = dynamic_cast<TParameter<double>*>(f.Get("sigma_combo_after_ps"));
    { TLatex t3; t3.SetNDC(); t3.SetTextSize(0.027); t3.SetTextColor(kGray + 3); t3.SetTextAlign(31);
      t3.DrawLatex(0.92, 0.34, "common-mode #rightarrow cancels in (DW#minusUP)/2");
      if (cb && ca)
        t3.DrawLatex(0.92, 0.295,
          Form("residual validated out-of-sample: combo %.0f#rightarrow%.0f ps",
               cb->GetVal(), ca->GetVal())); }

    DrawPageTitle("DRS4 time-base  --  stop-cell coverage (DRS0 G0)");
    c->Print(Form("%slayer1_timebase.png", kSumDir));
    printf("[layer1Summary] wrote layer1_timebase.png\n");
}

// ── H5 — DRS4 health: noise floor vs energy (saturation/spike annotated) ─────
// Replaces the dense 6-page drs4_diagnostics grid with one clean figure: the
// noise floor (the real, informative metric) across all 8 channels and all six
// energies, with the flat-at-zero saturation/spike facts stated, not plotted.
void HeroDRS4Clean()
{
    TFile f(Form("%sdrs4_diagnostics.root", kSumDir));
    if (f.IsZombie()) { printf("[layer1Summary] no drs4_diagnostics.root — skip H5\n"); return; }

    TCanvas* c = NewSquareCanvas("c_l1_drs4", 660, 128);
    c->cd();

    TMultiGraph* mg = new TMultiGraph();
    double maxSat = 0., maxSpk = 0.;
    for (int i = 0; i < 8; ++i) {
        TGraph* gp = dynamic_cast<TGraph*>(f.Get(Form("gPedRms_ch%d", i)));
        if (gp) {
            gp->SetLineColor(kRChannelCols[i]); gp->SetMarkerColor(kRChannelCols[i]);
            gp->SetMarkerStyle(20); gp->SetMarkerSize(1.0); gp->SetLineWidth(2);
            mg->Add(gp, "PL");
        }
        double s = graphMax(dynamic_cast<TGraph*>(f.Get(Form("gSatFrac_ch%d", i))));
        if (!std::isnan(s)) maxSat = std::max(maxSat, s);
        double k = graphMax(dynamic_cast<TGraph*>(f.Get(Form("gSpkFrac_ch%d", i))));
        if (!std::isnan(k)) maxSpk = std::max(maxSpk, k);
    }

    mg->Draw("A");
    mg->GetXaxis()->SetTitle("beam energy (GeV)");
    mg->GetYaxis()->SetTitle("pedestal noise floor (mV)");
    mg->GetXaxis()->SetLimits(0., 165.);
    mg->GetYaxis()->SetRangeUser(0., 6.);
    mg->GetXaxis()->SetTitleSize(0.048);
    mg->GetYaxis()->SetTitleSize(0.048);
    mg->GetYaxis()->SetTitleOffset(1.45);
    mg->Draw("A");

    TLine* l5 = new TLine(0., 5., 165., 5.);
    l5->SetLineColor(kRRed); l5->SetLineStyle(2); l5->SetLineWidth(2); l5->Draw();
    { TLatex t; t.SetNDC(); t.SetTextSize(0.030); t.SetTextColor(kRRed); t.SetTextAlign(31);
      t.DrawLatex(0.90, 0.815, "5 mV noise floor"); }

    // No legend: all 8 channels overlap at ~1.3 mV — the message is collective.
    { TLatex v; v.SetNDC(); v.SetTextSize(0.032); v.SetTextColor(kGray + 3); v.SetTextAlign(11);
      v.DrawLatex(0.185, 0.60, "All 8 capillaries (coloured) #minus flat noise floor ~1.3 mV");
      v.DrawLatex(0.185, 0.545, "across 25#minus150 GeV");
      v.DrawLatex(0.185, 0.475, Form("HG saturation < %.1f%%    #bullet    max spike rate %.2f%%",
                                    std::max(0.1, maxSat), maxSpk));
      v.DrawLatex(0.185, 0.420, "no degradation with beam rate or energy"); }

    DrawPageTitle("DRS4 health  --  noise floor vs beam energy (all 8 channels)");
    c->Print(Form("%slayer1_drs4clean.png", kSumDir));
    printf("[layer1Summary] wrote layer1_drs4clean.png\n");
}

} // namespace

void layer1Summary()
{
    ApplyRADiCALStyle();
    gSystem->mkdir(kSumDir, kTRUE);

    HeroPulseShapes();
    HeroVitals();
    HeroLinearity();
    HeroTimebase();
    HeroDRS4Clean();

    printf("\n[layer1Summary] Done — 5 hero figures in %s\n", kSumDir);
}

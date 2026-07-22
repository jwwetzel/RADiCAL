// ============================================================================
// lgWindowCheck.C — INVESTIGATION (not a pipeline product):
// Where does the LG pulse sit in the RAW DRS1 acquisition window?
//
// The average-waveform figure is aligned on each channel's own CFD-20%
// crossing (t - t_CFD), which pins every pulse to the left edge BY
// CONSTRUCTION and hides the raw window position.  This macro answers the
// actual question with UNALIGNED data:
//   (1) unaligned average LG waveform vs raw time (0..~1024 ns),
//   (2) distribution of the LG CFD-20% crossing time in the raw window,
//   (3) same for HG (DRS0 @5 GS/s, ~205 ns window) for contrast,
//   (4) per-channel table: median crossing, pedestal-window overlap
//       (pedestal = samples 3..52 => first ~53 ns), edge-truncation fractions.
//
// Run from repo root:  root -l -b -q '<scratch>/lgWindowCheck.C+'
// ============================================================================
#include "ChannelConfig.h"
#include "WaveformUtils.h"
#include "SelectionCuts.h"

#include "TChain.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TProfile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TROOT.h"
#include <algorithm>
#include <cstdio>

static const Long64_t kNMaxEv = 20000;

void lgWindowCheck()
{
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);

    TChain* chain = new TChain("pulse");
    chain->Add("data/2023/raw/RUN1258_150_GeV.root");
    chain->Add("data/2023/raw/RUN1259_150_GeV.root");
    chain->Add("data/2023/raw/RUN1260_150_GeV.root");
    chain->Add("data/2023/raw/RUN1261_150_GeV.root");
    printf("[lgWindowCheck] chain entries: %lld\n", chain->GetEntries());

    TProfile* pLG[8];  TH1F* hLGc[8];  TH1F* hHGc[8];
    for (int i = 0; i < 8; ++i) {
        pLG[i]  = new TProfile(Form("pLG_%d", i), ";raw time in DRS1 window (ns);mean ped#minusA (mV)", 256, 0., 1024.);
        pLG[i]->SetDirectory(nullptr);
        hLGc[i] = new TH1F(Form("hLGc_%d", i), ";LG CFD-20%% crossing, raw window time (ns);events", 256, 0., 1024.);
        hLGc[i]->SetDirectory(nullptr);
        hHGc[i] = new TH1F(Form("hHGc_%d", i), ";HG CFD-20%% crossing, raw window time (ns);events", 210, 0., 210.);
        hHGc[i]->SetDirectory(nullptr);
    }
    double tLGmin = 1e9, tLGmax = -1e9;   // actual DRS1 time-axis span seen

    TTreeReader reader(chain);
    TTreeReaderArray<float> time_v(reader, "timevalue");
    TTreeReaderArray<float> amp_v (reader, "amplitude");

    Long64_t nProc = 0, nUsed = 0;
    while (reader.Next() && nProc < kNMaxEv) {
        ++nProc;
        const float* A = &amp_v[0];
        const float* T = &time_v[0];

        // same event gate as averageWaveforms.C
        Pulse mcp = ExtractPulse(T + kMCP1_t, A + kMCP1, 0.20f, kMCP1_minPeak);
        if (!mcp.valid || mcp.peak < kMCP1_minPeak || mcp.peak > kMCP1_maxPeak) continue;
        Pulse wcR = ExtractPulse(T + kWC_t, A + kWC_R, 0.50f, kWC_minPeak);
        Pulse wcL = ExtractPulse(T + kWC_t, A + kWC_L, 0.50f, kWC_minPeak);
        Pulse wcD = ExtractPulse(T + kWC_t, A + kWC_D, 0.50f, kWC_minPeak);
        Pulse wcU = ExtractPulse(T + kWC_t, A + kWC_U, 0.50f, kWC_minPeak);
        if (!(wcR.valid && wcL.valid && wcD.valid && wcU.valid)) continue;
        ++nUsed;

        for (int i = 0; i < 8; ++i) {
            // LG — unaligned
            const float* Aa = A + kCap[i].lg;
            const float* Ta = T + kCap[i].lg_t;
            Pulse p = ExtractPulse(Ta, Aa, 0.20f, kLG_minPeak);
            if (p.valid && p.crossingTime != kNoTime) {
                hLGc[i]->Fill(p.crossingTime);
                for (int s = 0; s < 1024; ++s) pLG[i]->Fill(Ta[s], p.pedestal - Aa[s]);
                if (Ta[0]    < tLGmin) tLGmin = Ta[0];
                if (Ta[1023] > tLGmax) tLGmax = Ta[1023];
            }
            // HG — crossing only
            Pulse q = ExtractPulse(T + kCap[i].hg_t, A + kCap[i].hg, 0.20f, kHG_minPeak);
            if (q.valid && q.crossingTime != kNoTime) hHGc[i]->Fill(q.crossingTime);
        }
    }
    printf("[lgWindowCheck] processed %lld, passed gate %lld\n", nProc, nUsed);
    printf("[lgWindowCheck] DRS1 time-axis span seen: %.1f .. %.1f ns\n", tLGmin, tLGmax);

    // ── per-channel table ───────────────────────────────────────────────────
    printf("\n%-6s %8s %8s %8s %10s %10s %10s\n",
           "chan", "N", "median", "p1..p99", "f(<53ns)", "f(<10ns)", "f(>900ns)");
    double medLG[8];
    for (int i = 0; i < 8; ++i) {
        double q[3], pr[3] = {0.01, 0.50, 0.99};
        hLGc[i]->GetQuantiles(3, q, pr);
        medLG[i] = q[1];
        double n     = hLGc[i]->GetEntries();
        double f53   = hLGc[i]->Integral(1, hLGc[i]->FindBin(53.)) / std::max(1., n);
        double f10   = hLGc[i]->Integral(1, hLGc[i]->FindBin(10.)) / std::max(1., n);
        double f900  = hLGc[i]->Integral(hLGc[i]->FindBin(900.), 257) / std::max(1., n);
        printf("%-6s %8.0f %7.1f %5.0f..%-5.0f %9.4f %9.4f %9.4f\n",
               kCap[i].name, n, q[1], q[0], q[2], f53, f10, f900);
    }

    // ── canvas ──────────────────────────────────────────────────────────────
    const int cols[8] = {kBlue+1, kAzure+7, kTeal+3, kGreen+2,
                         kOrange+7, kRed+1, kMagenta+2, kGray+2};
    TCanvas* c = new TCanvas("c_lgwin", "", 1500, 1000);
    c->Divide(2, 2, 0.01, 0.02);

    c->cd(1);  gPad->SetLeftMargin(0.11); gPad->SetBottomMargin(0.12);
    for (int i = 0; i < 8; ++i) {
        TH1D* h = pLG[i]->ProjectionX(Form("pLGpx_%d", i));
        double m = h->GetMaximum(); if (m > 0) h->Scale(1./m);
        h->SetLineColor(cols[i]); h->SetLineWidth(2);
        h->GetYaxis()->SetRangeUser(-0.55, 1.35);
        h->SetTitle(";raw time in DRS1 window (ns);a(t)/a_{max}  (unaligned mean)");
        h->Draw(i ? "HIST same" : "HIST");
    }
    { TLegend* L = new TLegend(0.55, 0.55, 0.93, 0.92); L->SetBorderSize(0); L->SetFillStyle(0);
      for (int i = 0; i < 8; ++i) L->AddEntry(pLG[i], kCap[i].name, "l"); L->Draw();
      TLatex t; t.SetNDC(); t.SetTextSize(0.045);
      t.DrawLatex(0.13, 0.94, "UNALIGNED average LG waveforms (raw window position)"); }

    c->cd(2);  gPad->SetLeftMargin(0.11); gPad->SetBottomMargin(0.12); gPad->SetLogy();
    for (int i = 0; i < 8; ++i) {
        hLGc[i]->SetLineColor(cols[i]); hLGc[i]->SetLineWidth(2);
        hLGc[i]->Draw(i ? "HIST same" : "HIST");
    }
    { TLine* l = new TLine(53., 0., 53., hLGc[0]->GetMaximum()*1.1);
      l->SetLineColor(kGray+2); l->SetLineStyle(2); l->Draw();
      TLatex t; t.SetNDC(); t.SetTextSize(0.045);
      t.DrawLatex(0.13, 0.94, "LG CFD-20%% crossing time in the raw window");
      t.SetTextSize(0.032); t.SetTextColor(kGray+3);
      t.DrawLatex(0.30, 0.85, "dashed: pedestal window edge (samples 3-52 #approx 53 ns)"); }

    c->cd(3);  gPad->SetLeftMargin(0.11); gPad->SetBottomMargin(0.12); gPad->SetLogy();
    for (int i = 0; i < 8; ++i) {
        hHGc[i]->SetLineColor(cols[i]); hHGc[i]->SetLineWidth(2);
        hHGc[i]->Draw(i ? "HIST same" : "HIST");
    }
    { TLatex t; t.SetNDC(); t.SetTextSize(0.045);
      t.DrawLatex(0.13, 0.94, "HG CFD-20%% crossing time (DRS0 window, ~205 ns)"); }

    c->cd(4);
    { TLatex t; t.SetNDC(); t.SetTextSize(0.038); t.SetTextFont(42);
      t.DrawLatex(0.06, 0.92, "Per-LG-channel raw-window placement (150 GeV):");
      for (int i = 0; i < 8; ++i)
          t.DrawLatex(0.08, 0.84 - 0.065*i,
              Form("%s:  median crossing = %.0f ns", kCap[i].name, medLG[i]));
      t.SetTextColor(kGray+3); t.SetTextSize(0.032);
      t.DrawLatex(0.06, 0.26, "Pedestal uses samples 3-52 (#approx 3-53 ns at 1 GS/s).");
      t.DrawLatex(0.06, 0.20, "If median crossing #lessgtr 53 ns the pedestal window");
      t.DrawLatex(0.06, 0.14, "does / does not overlap the rising edge."); }

    c->Print("output/investigations/lg_window_check.png");
    printf("[lgWindowCheck] wrote output/investigations/lg_window_check.png\n");
}

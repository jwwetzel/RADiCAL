// ===========================================================================
// layer2Summary.C  —  compact, Ledovskoy-clean "Layer 2: Reference
//                     Characterization" hero figures (square frames).
//
// The story: we understand our timing reference (the MCP) and it does NOT
// limit the headline result.
//
//   layer2_mcp_jitter.png   H1  sigma_MCP,single vs energy  (flat ~71 ps)
//   layer2_sub_mcp.png      H2  the MCP floor vs the (DW-UP)/2 headline:
//                               the headline sits BELOW the reference floor
//                               because the corner difference cancels the MCP jitter
//   layer2_per_channel.png  H3  single-channel sigma_t per capillary vs energy
//
// Reads persisted analysis products in Output/Summary/*.root (no raw data).
//
// Run:  ROOT_INCLUDE_PATH=Analysis root -l -b -q 'Analysis/layer2Summary.C+'
// ===========================================================================
#include <cmath>
#include <algorithm>

#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TSystem.h"

#include "RADiCALStyle.h"
#include "ChannelConfig.h"

namespace {

const char* kSumDir = "output/Summary/";

double graphMean(TGraph* g)
{
    if (!g || g->GetN() == 0) return NAN;
    double x, y, s = 0.; int n = g->GetN();
    for (int i = 0; i < n; ++i) { g->GetPoint(i, x, y); s += y; }
    return s / n;
}

// ── H1 — per-group (inter-mezzanine) reference jitter vs beam energy ───────
//   MCP1 and MCP2 are ONE MCP split into the two DT5742 groups (corr=1.000), so
//   the MCP's own jitter cancels in MCP1−MCP2; what remains is the group-to-group
//   DRS4 timing jitter.  gMCP_jitter = σ(MCP1−MCP2)/√2 = per-group reference jitter.
void HeroMcpJitter()
{
    TFile f(Form("%smcp_jitter.root", kSumDir));
    if (f.IsZombie()) { printf("[layer2Summary] no mcp_jitter.root — skip H1\n"); return; }
    TGraphErrors* g = dynamic_cast<TGraphErrors*>(f.Get("gMCP_jitter"));
    if (!g) { printf("[layer2Summary] gMCP_jitter missing — skip H1\n"); return; }

    TCanvas* c = NewSquareCanvas("c_l2_mcp");
    c->cd();

    const double mean = graphMean(g);
    g->SetMarkerStyle(20); g->SetMarkerSize(1.5); g->SetMarkerColor(kRData);
    g->SetLineColor(kRData); g->SetLineWidth(2);
    g->Draw("APL");
    g->GetXaxis()->SetTitle("beam energy (GeV)");
    g->GetYaxis()->SetTitle("#sigma(MCP1#minusMCP2) / #sqrt{2}   (ps)");
    g->GetXaxis()->SetLimits(0., 165.);
    g->GetYaxis()->SetRangeUser(0., 100.);
    g->GetXaxis()->SetTitleSize(0.046);
    g->GetYaxis()->SetTitleSize(0.046);

    TLine* lm = new TLine(0., mean, 165., mean);
    lm->SetLineColor(kGray + 1); lm->SetLineStyle(2); lm->SetLineWidth(1); lm->Draw();

    { TLatex t; t.SetNDC(); t.SetTextSize(0.033); t.SetTextColor(kGray + 3);
      t.DrawLatex(0.18, 0.34, Form("per-group reference jitter: %.0f ps", mean));
      t.SetTextSize(0.027); t.SetTextColor(kGray + 2);
      t.DrawLatex(0.18, 0.295, "one MCP split into both DT5742 groups #Rightarrow");
      t.DrawLatex(0.18, 0.255, "the MCP's own jitter cancels in MCP1#minusMCP2;");
      t.DrawLatex(0.18, 0.215, "this is the inter-group (DRS4) timing jitter,");
      t.DrawLatex(0.18, 0.175, "flat with energy (purely electronic)."); }

    DrawPageTitle("Per-group reference jitter vs beam energy");
    c->Print(Form("%slayer2_mcp_jitter.png", kSumDir));
    printf("[layer2Summary] wrote layer2_mcp_jitter.png\n");
}

// ── H2 — the headline sits below the per-group floor (cancellation) ────────
void HeroSubMcp()
{
    TFile fM(Form("%smcp_jitter.root", kSumDir));
    TFile fT(Form("%stiming_energy_bins.root", kSumDir));
    if (fM.IsZombie() || fT.IsZombie()) { printf("[layer2Summary] missing input — skip H2\n"); return; }
    TGraphErrors* gM = dynamic_cast<TGraphErrors*>(fM.Get("gMCP_jitter"));
    TGraphErrors* gT = dynamic_cast<TGraphErrors*>(fT.Get("gBestSigma_teb_m0"));
    if (!gM || !gT) { printf("[layer2Summary] graphs missing — skip H2\n"); return; }

    TCanvas* c = NewSquareCanvas("c_l2_submcp");
    c->cd();

    gM->SetMarkerStyle(24); gM->SetMarkerSize(1.4); gM->SetMarkerColor(kROrange);
    gM->SetLineColor(kROrange); gM->SetLineWidth(2); gM->SetLineStyle(2);
    gT->SetMarkerStyle(20); gT->SetMarkerSize(1.4); gT->SetMarkerColor(kRData);
    gT->SetLineColor(kRData); gT->SetLineWidth(2);

    TMultiGraph* mg = new TMultiGraph();
    mg->Add(gM, "PL");
    mg->Add(gT, "PL");
    mg->Draw("A");
    mg->GetXaxis()->SetTitle("beam energy (GeV)");
    mg->GetYaxis()->SetTitle("#sigma_{t} (ps)");
    mg->GetXaxis()->SetLimits(0., 165.);
    mg->GetYaxis()->SetRangeUser(0., 100.);
    mg->GetXaxis()->SetTitleSize(0.048);
    mg->GetYaxis()->SetTitleSize(0.048);

    TLegend* leg = new TLegend(0.40, 0.795, 0.95, 0.90);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextFont(42); leg->SetTextSize(0.033);
    leg->AddEntry(gM, "per-group reference floor", "lp");
    leg->AddEntry(gT, "(DW#minusUP)/2 headline", "lp");
    leg->Draw();

    { TLatex t; t.SetNDC(); t.SetTextSize(0.029); t.SetTextColor(kGray + 3);
      t.DrawLatex(0.17, 0.32, "Each channel is referenced to its OWN group's MCP;");
      t.DrawLatex(0.17, 0.280, "the corner difference cancels it (and the");
      t.DrawLatex(0.17, 0.240, "inter-group jitter), so the headline sits");
      t.DrawLatex(0.17, 0.200, "well below the per-group floor."); }

    DrawPageTitle("Why the headline beats the per-group floor");
    c->Print(Form("%slayer2_sub_mcp.png", kSumDir));
    printf("[layer2Summary] wrote layer2_sub_mcp.png\n");
}

// ── H3 — single-channel timing resolution per capillary vs energy ──────────
void HeroPerChannel()
{
    TFile f(Form("%ssummary.root", kSumDir));
    if (f.IsZombie()) { printf("[layer2Summary] no summary.root — skip H3\n"); return; }

    TCanvas* c = NewSquareCanvas("c_l2_chan");
    c->cd();

    TMultiGraph* mg = new TMultiGraph();
    TLegend* leg = MakeCornerLegend(8, "tr", 0.033);
    bool any = false;
    for (int i = 0; i < 8; ++i) {
        TGraphErrors* g = dynamic_cast<TGraphErrors*>(f.Get(Form("gTRes_%s", kCap[i].name)));
        if (!g) continue;
        any = true;
        g->SetLineColor(kRChannelCols[i]);   g->SetMarkerColor(kRChannelCols[i]);
        g->SetMarkerStyle(20); g->SetMarkerSize(1.0); g->SetLineWidth(2);
        mg->Add(g, "PL");
        leg->AddEntry(g, kCap[i].name, "lp");
    }
    if (!any) { printf("[layer2Summary] no gTRes_* — skip H3\n"); return; }

    mg->Draw("A");
    mg->GetXaxis()->SetTitle("beam energy (GeV)");
    mg->GetYaxis()->SetTitle("single-channel #sigma_{t} (ps)");
    mg->GetXaxis()->SetLimits(0., 165.);
    mg->GetXaxis()->SetTitleSize(0.048);
    mg->GetYaxis()->SetTitleSize(0.048);
    leg->Draw();

    DrawPageTitle("Single-channel timing resolution per capillary");
    c->Print(Form("%slayer2_per_channel.png", kSumDir));
    printf("[layer2Summary] wrote layer2_per_channel.png\n");
}

} // namespace

void layer2Summary()
{
    ApplyRADiCALStyle();
    gSystem->mkdir(kSumDir, kTRUE);

    HeroMcpJitter();
    HeroSubMcp();
    HeroPerChannel();

    printf("\n[layer2Summary] Done — hero figures in %s\n", kSumDir);
}

// ===========================================================================
// layer3Summary.C  —  compact, Ledovskoy-clean "Layer 3: Beam Characterization"
//                     hero figures (square frames).
//
// The story: we ran a clean, well-aimed electron beam, we SEE and remove the
// hadronic contamination, and the showers are well contained.
//
//   layer3_beam_map.png    H1  beam illumination (x,y) with the fiducial circles
//   layer3_containment.png H2  Sigma_PbGlass vs Sigma_LG with the 30% cut line
//                              (EM band below, hadronic punch-through above)
//   layer3_quality.png     H3  sample quality vs energy: shower containment (high)
//                              and hadronic punch-through (low), both vs beam energy
//
// H1/H2 are computed from the 150 GeV ntuple (tree 'rad'); H3 reads the
// persisted gContainedFrac / gPunchThroughFrac graphs.
//
// Run:  ROOT_INCLUDE_PATH=Analysis root -l -b -q 'Analysis/layer3Summary.C+'
// ===========================================================================
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TAxis.h"
#include "TEllipse.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TSystem.h"

#include "RADiCALStyle.h"
#include "ChannelConfig.h"
#include "SelectionCuts.h"

namespace {

const char* kSumDir = "Analysis/Output/Summary/";
const char* kRefNtuple = "Analysis/Output/150GeV/ntuple.root";   // representative

// ── H1 + H2 share one pass over the 150 GeV ntuple ─────────────────────────
void HeroBeamAndContainment()
{
    TFile* f = TFile::Open(kRefNtuple);
    if (!f || f->IsZombie()) { printf("[layer3Summary] no 150 GeV ntuple — skip H1/H2\n"); if (f) delete f; return; }
    TTree* t = dynamic_cast<TTree*>(f->Get("rad"));
    if (!t) { printf("[layer3Summary] no 'rad' tree — skip H1/H2\n"); delete f; return; }

    Float_t x_trk = 0, y_trk = 0, mcp_peak = 0, sum_lg = 0, sum_pb = 0;
    Bool_t  wc_ok = false, in_fid = false;
    t->SetBranchAddress("x_trk", &x_trk);
    t->SetBranchAddress("y_trk", &y_trk);
    t->SetBranchAddress("mcp_peak", &mcp_peak);
    t->SetBranchAddress("sum_lg", &sum_lg);
    t->SetBranchAddress("sum_pb", &sum_pb);
    t->SetBranchAddress("wc_ok", &wc_ok);
    t->SetBranchAddress("in_fiducial", &in_fid);

    const Long64_t N = t->GetEntries();

    // Pass 1: per-run fiducial centre = mean (x,y) of in-fiducial events.
    double sx = 0, sy = 0; long nfid = 0;
    for (Long64_t i = 0; i < N; ++i) { t->GetEntry(i); if (in_fid) { sx += x_trk; sy += y_trk; ++nfid; } }
    const double xc = nfid ? sx / nfid : 0., yc = nfid ? sy / nfid : 0.;

    // Bin at the canonical WC resolution (kWC_resBin = 1 mm, SelectionCuts.h) —
    // the same grid used by the qualityPlots and compareEnergies hit maps.
    // The wire-chamber study found the track position quantises in ~0.25 mm
    // steps (x) with a ~0.5 mm peak-sample comb (y); hit/illumination maps
    // belong at 1-2 mm.  The previous hardcoded 0.33 mm binning resolved that
    // comb -> cross-hatch + aliased dead-wire stripes, not the true beam profile.
    const int nBeam = static_cast<int>(std::round(20. / kWC_resBin));  // 20 mm span
    TH2F* hBeam = new TH2F("hBeam", "", nBeam, xc - 10., xc + 10., nBeam, yc - 10., yc + 10.);
    hBeam->SetDirectory(nullptr);
    TH2F* hPb = new TH2F("hPb", "", 120, 0., 6500., 120, 0., 3000.);
    hPb->SetDirectory(nullptr);

    // Pass 2: fill both maps.
    for (Long64_t i = 0; i < N; ++i) {
        t->GetEntry(i);
        if (wc_ok && mcp_peak > kMCP1_minPeak && mcp_peak < kMCP1_maxPeak)
            hBeam->Fill(x_trk, y_trk);
        if (wc_ok && in_fid)
            hPb->Fill(sum_lg, sum_pb);
    }

    gStyle->SetPalette(kRust); TColor::InvertPalette();

    // ── H1: beam illumination + fiducial circles ──────────────────────────
    {
        TCanvas* c = NewSquareCanvas("c_l3_beam", 600, 104, 92, 54, 138);
        c->cd();
        hBeam->GetXaxis()->SetTitle("x_{track} (mm)");
        hBeam->GetYaxis()->SetTitle("y_{track} (mm)");
        hBeam->GetZaxis()->SetTitle("events");
        hBeam->GetXaxis()->SetTitleSize(0.046);
        hBeam->GetYaxis()->SetTitleSize(0.046);
        hBeam->Draw("COLZ");

        TEllipse* eT = new TEllipse(xc, yc, kFiducial_r_timing, kFiducial_r_timing);
        eT->SetFillStyle(0); eT->SetLineColor(kROrange); eT->SetLineWidth(3); eT->Draw();
        TEllipse* eE = new TEllipse(xc, yc, kFiducial_r_energy, kFiducial_r_energy);
        eE->SetFillStyle(0); eE->SetLineColor(kRRed); eE->SetLineWidth(3); eE->SetLineStyle(2); eE->Draw();

        { TLatex l; l.SetNDC(); l.SetTextSize(0.030); l.SetTextColor(kROrange);
          l.DrawLatex(0.17, 0.885, Form("timing fiducial  r < %.0f mm", kFiducial_r_timing)); }
        { TLatex l; l.SetNDC(); l.SetTextSize(0.030); l.SetTextColor(kRRed);
          l.DrawLatex(0.17, 0.845, Form("energy fiducial  r < %.0f mm", kFiducial_r_energy)); }

        DrawPageTitle("Beam illumination & fiducial  (150 GeV, wire chamber)");
        c->Print(Form("%slayer3_beam_map.png", kSumDir));
        printf("[layer3Summary] wrote layer3_beam_map.png\n");
    }

    // ── H2: containment cut — Sigma_PbGlass vs Sigma_LG ────────────────────
    {
        TCanvas* c = NewSquareCanvas("c_l3_cont", 600, 104, 92, 54, 138);
        c->cd();
        gPad->SetLogz(1);
        hPb->GetXaxis()->SetTitle("#Sigma_{LG}  (mV)");
        hPb->GetYaxis()->SetTitle("#Sigma_{PbGlass}  (mV)");
        hPb->GetZaxis()->SetTitle("events");
        hPb->GetXaxis()->SetTitleSize(0.046);
        hPb->GetYaxis()->SetTitleSize(0.046);
        hPb->Draw("COLZ");

        TLine* cut = new TLine(0., 0., 6500., kPb_maxRatio * 6500.);
        cut->SetLineColor(kRRed); cut->SetLineStyle(2); cut->SetLineWidth(3); cut->Draw();

        { TLatex l; l.SetNDC(); l.SetTextSize(0.030); l.SetTextColor(kRRed); l.SetTextAlign(31);
          l.DrawLatex(0.62, 0.86, Form("#Sigma_{PbGlass} = %.0f%% #times #Sigma_{LG}", 100.*kPb_maxRatio)); }
        { TLatex l; l.SetNDC(); l.SetTextSize(0.029); l.SetTextColor(kGray + 3);
          l.DrawLatex(0.40, 0.30, "contained EM showers");
          l.DrawLatex(0.40, 0.255, "(kept, below the cut)");
          l.SetTextColor(kRRed);
          l.DrawLatex(0.18, 0.74, "hadronic punch-through"); }

        DrawPageTitle("Shower containment cut  (150 GeV, in-fiducial)");
        c->Print(Form("%slayer3_containment.png", kSumDir));
        printf("[layer3Summary] wrote layer3_containment.png\n");
    }

    delete hBeam; delete hPb; delete f;
}

// ── H3 — sample quality vs beam energy (containment + punch-through) ───────
void HeroQualityVsEnergy()
{
    TFile fC(Form("%scross_energy.root", kSumDir));
    TFile fP(Form("%spbglass_investigation.root", kSumDir));
    if (fC.IsZombie() || fP.IsZombie()) { printf("[layer3Summary] missing graph input — skip H3\n"); return; }
    TGraph* gCont  = dynamic_cast<TGraph*>(fC.Get("gContainedFrac"));
    TGraph* gPunch = dynamic_cast<TGraph*>(fP.Get("gPunchThroughFrac"));
    if (!gCont || !gPunch) { printf("[layer3Summary] graphs missing — skip H3\n"); return; }

    TCanvas* c = NewSquareCanvas("c_l3_quality");
    c->cd();

    gCont->SetMarkerStyle(20);  gCont->SetMarkerSize(1.4);  gCont->SetMarkerColor(kRGreen);
    gCont->SetLineColor(kRGreen);  gCont->SetLineWidth(2);
    gPunch->SetMarkerStyle(21); gPunch->SetMarkerSize(1.4); gPunch->SetMarkerColor(kRRed);
    gPunch->SetLineColor(kRRed); gPunch->SetLineWidth(2);

    TMultiGraph* mg = new TMultiGraph();
    mg->Add(gCont, "PL");
    mg->Add(gPunch, "PL");
    mg->Draw("A");
    mg->GetXaxis()->SetTitle("beam energy (GeV)");
    mg->GetYaxis()->SetTitle("fraction of in-fiducial events (%)");
    mg->GetXaxis()->SetLimits(0., 165.);
    mg->GetYaxis()->SetRangeUser(0., 105.);
    mg->GetXaxis()->SetTitleSize(0.046);
    mg->GetYaxis()->SetTitleSize(0.046);

    TLegend* leg = new TLegend(0.40, 0.45, 0.93, 0.58);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextFont(42); leg->SetTextSize(0.034);
    leg->AddEntry(gCont,  "shower containment", "lp");
    leg->AddEntry(gPunch, "hadronic punch-through", "lp");
    leg->Draw();

    { TLatex l; l.SetNDC(); l.SetTextSize(0.029); l.SetTextColor(kGray + 3);
      l.DrawLatex(0.18, 0.30, "Showers are >90% contained; hadronic");
      l.DrawLatex(0.18, 0.255, "punch-through is a few % (rising to 13.6%");
      l.DrawLatex(0.18, 0.210, "at 150 GeV from the SPS pion fraction)."); }

    DrawPageTitle("Sample quality vs beam energy");
    c->Print(Form("%slayer3_quality.png", kSumDir));
    printf("[layer3Summary] wrote layer3_quality.png\n");
}

} // namespace

void layer3Summary()
{
    ApplyRADiCALStyle();
    gSystem->mkdir(kSumDir, kTRUE);

    HeroBeamAndContainment();
    HeroQualityVsEnergy();

    printf("\n[layer3Summary] Done — hero figures in %s\n", kSumDir);
}

// ===========================================================================
// moduleCenter.C  —  data-driven shashlik centre from the calorimeter EDGES.
//
// The module centre is found geometrically (beam-independent): project the mean
// summed high-gain amplitude of the 8 RADiCAL capillaries (Sigma HG) onto x (for
// tracks in a central y-band) and onto y (central x-band).  Each projection is a
// flat-topped plateau with sharp shoulders at the module edges.  The low-gain
// energy sum (Sigma LG) is carried as an independent cross-check and gives the
// same centre to <0.05 mm.  The half-maximum crossings give the two edges;
// their midpoint is the centre on that axis:
//        x0 = (x_edgeL + x_edgeR) / 2 ,   y0 = (y_edgeL + y_edgeR) / 2 .
// Edges are located by scanning OUTWARD from the plateau peak to the first
// half-max crossing, so off-module structure (Pb-glass / leakage backsplash
// seen at large |y|) does not bias the result.
//
// Unlike an amplitude-weighted centroid, the edge midpoint is independent of the
// beam profile within the module — the module does not move with energy, so the
// per-energy results agree and are combined for the headline value.
//
// Reads the reduced per-energy ntuples (Output/<E>/ntuple.root); no raw data.
//
// Run:  ROOT_INCLUDE_PATH=Analysis root -l -b -q 'Analysis/moduleCenter.C+'
// ===========================================================================
#include "ChannelConfig.h"
#include "PlotUtils.h"        // ApplyRADiCALStyle, StylePad, colours
#include "RADiCALStyle.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLatex.h"
#include "TParameter.h"
#include "TSystem.h"
#include "TStyle.h"

#include <cstdio>

namespace {

// Linear-interpolated half-max crossing between two adjacent bins.
double crossing(TProfile* h, int b1, int b2, double half)
{
    double x1 = h->GetBinCenter(b1), x2 = h->GetBinCenter(b2);
    double v1 = h->GetBinContent(b1), v2 = h->GetBinContent(b2);
    return (v2 == v1) ? x1 : x1 + (half - v1) * (x2 - x1) / (v2 - v1);
}

// Find the two module edges (half-max, scanning outward from the plateau peak)
// and return the centre.  minN = minimum entries per bin to trust it.
double edgeCentre(TProfile* h, int minN, double& eL, double& eR, double& width)
{
    const int nb = h->GetNbinsX();
    int pb = 0; double plat = 0.;
    for (int b = 1; b <= nb; ++b)
        if (h->GetBinEntries(b) > minN && h->GetBinContent(b) > plat) { plat = h->GetBinContent(b); pb = b; }
    const double half = 0.5 * plat;
    int bl = pb; while (bl > 1)  { if (h->GetBinEntries(bl) > minN && h->GetBinContent(bl) < half) break; --bl; }
    int br = pb; while (br < nb) { if (h->GetBinEntries(br) > minN && h->GetBinContent(br) < half) break; ++br; }
    eL = crossing(h, bl, bl + 1, half);
    eR = crossing(h, br - 1, br, half);
    width = eR - eL;
    return 0.5 * (eL + eR);
}

// Fill the x- and y-edge profiles from a tree (appends; reuse across a chain).
// Primary signal = summed high-gain amplitude of the 8 RADiCAL capillaries
// (Sigma HG).  If the LG profiles are supplied, also fill them with the low-gain
// energy sum (sum_lg) as an independent cross-check.
void fillProfiles(TTree* t, TProfile* hx, TProfile* hy,
                  TProfile* hxLG, TProfile* hyLG,
                  double bandHalf, double seedX, double seedY)
{
    Float_t x, y, slg, hg[8]; Bool_t ok;
    t->SetBranchAddress("x_trk", &x); t->SetBranchAddress("y_trk", &y);
    t->SetBranchAddress("hg_peak", hg); t->SetBranchAddress("sum_lg", &slg);
    t->SetBranchAddress("wc_ok", &ok);
    const Long64_t N = t->GetEntries();
    for (Long64_t i = 0; i < N; ++i) {
        t->GetEntry(i);
        if (!ok) continue;
        double shg = 0.; for (int k = 0; k < 8; ++k) shg += hg[k];
        if (std::fabs(y - seedY) < bandHalf) { hx->Fill(x, shg); if (hxLG) hxLG->Fill(x, slg); }
        if (std::fabs(x - seedX) < bandHalf) { hy->Fill(y, shg); if (hyLG) hyLG->Fill(y, slg); }
    }
}

} // namespace

void moduleCenter()
{
    ApplyRADiCALStyle();
    gStyle->SetOptStat(0);
    const double seedX = kCalo_x0, seedY = kCalo_y0;   // nominal seed for the bands
    const double bandHalf = 4.0;                        // mm, central band half-width
    const int    minN = 50;

    // ── Per-energy centres (consistency check — module is energy-independent) ──
    std::printf("\n  Module centre per energy (edge midpoints):\n");
    std::printf("    E (GeV)    x0 (mm)   y0 (mm)   wx (mm)  wy (mm)\n");
    TGraph gX, gY; gX.SetName("gCenterX"); gY.SetName("gCenterY");
    for (int e = 0; e < kNRuns; ++e) {
        TString fn = Form("Analysis/Output/%s/ntuple.root", kRuns[e].label.Data());
        TFile* f = TFile::Open(fn);
        if (!f || f->IsZombie()) { if (f) delete f; continue; }
        TTree* t = (TTree*)f->Get("rad");
        if (!t) { delete f; continue; }
        TProfile hx("hx_e", "", 90, -15, 30), hy("hy_e", "", 90, -20, 25);
        fillProfiles(t, &hx, &hy, nullptr, nullptr, bandHalf, seedX, seedY);
        double exL, exR, wx, eyL, eyR, wy;
        double cx = edgeCentre(&hx, minN, exL, exR, wx);
        double cy = edgeCentre(&hy, minN, eyL, eyR, wy);
        std::printf("    %5.0f      %6.2f    %6.2f    %5.1f    %5.1f\n",
                    kRuns[e].energy_GeV, cx, cy, wx, wy);
        gX.SetPoint(gX.GetN(), kRuns[e].energy_GeV, cx);
        gY.SetPoint(gY.GetN(), kRuns[e].energy_GeV, cy);
        delete f;
    }

    // ── Combined (all energies chained) → headline centre + plot ───────────────
    TChain ch("rad");
    for (int e = 0; e < kNRuns; ++e)
        ch.Add(Form("Analysis/Output/%s/ntuple.root", kRuns[e].label.Data()));
    if (ch.GetEntries() == 0) { std::printf("[moduleCenter] no ntuples found.\n"); return; }

    TProfile* hX = new TProfile("hModEdgeX",
        "Shashlik edges (X);x Track (mm);mean #Sigma HG, 8 capillaries (mV)", 90, -15, 30);
    TProfile* hY = new TProfile("hModEdgeY",
        "Shashlik edges (Y);y Track (mm);mean #Sigma HG, 8 capillaries (mV)", 90, -20, 25);
    TProfile hXlg("hX_lg", "", 90, -15, 30), hYlg("hY_lg", "", 90, -20, 25);  // LG cross-check
    fillProfiles(&ch, hX, hY, &hXlg, &hYlg, bandHalf, seedX, seedY);

    double xL, xR, wx, yL, yR, wy;
    const double x0 = edgeCentre(hX, minN, xL, xR, wx);
    const double y0 = edgeCentre(hY, minN, yL, yR, wy);
    double lxL, lxR, lwx, lyL, lyR, lwy;
    const double lx0 = edgeCentre(&hXlg, minN, lxL, lxR, lwx);
    const double ly0 = edgeCentre(&hYlg, minN, lyL, lyR, lwy);
    std::printf("\n  >>> MODULE CENTRE (combined, #Sigma HG of 8 capillaries): "
                "(%.2f, %.2f) mm   widths (%.1f, %.1f) mm\n",
                x0, y0, wx, wy);
    std::printf("      cross-check (#Sigma LG energy):            "
                "(%.2f, %.2f) mm   widths (%.1f, %.1f) mm\n",
                lx0, ly0, lwx, lwy);
    std::printf("      nominal kCalo_x0/y0 = (%.2f, %.2f) mm\n\n", kCalo_x0, kCalo_y0);

    // ── Plot: two edge profiles with edge + centre markers ─────────────────────
    TString sumDir = "Analysis/Output/Summary";
    gSystem->mkdir(sumDir, kTRUE);
    TCanvas c("c_modctr", "", 900, 1000);   // portrait: two stacked panels
    c.Divide(1, 2, 0.004, 0.020);
    auto drawOne = [&](TProfile* h, double cc, double eL, double eR, double w){
        h->SetLineColor(kRBlue); h->SetLineWidth(2);
        h->SetMinimum(0.); h->SetMaximum(h->GetMaximum() * 1.12);
        h->Draw("HIST L");
        const double ymax = h->GetMaximum();
        for (double xx : {eL, eR}) {
            TLine* l = new TLine(xx, 0, xx, ymax);
            l->SetLineColor(kGray + 2); l->SetLineStyle(2); l->SetLineWidth(2); l->Draw();
        }
        TLine* lc = new TLine(cc, 0, cc, ymax);
        lc->SetLineColor(kRRed); lc->SetLineWidth(3); lc->Draw();
        TLatex tt; tt.SetNDC(); tt.SetTextSize(0.050); tt.SetTextColor(kRRed);
        tt.DrawLatex(0.40, 0.80, Form("centre = %.2f mm", cc));
        TLatex te; te.SetNDC(); te.SetTextSize(0.036); te.SetTextColor(kGray + 2);
        te.DrawLatex(0.40, 0.73, Form("edges %.1f / %.1f  (w %.1f)", eL, eR, w));
    };
    c.cd(1); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13); drawOne(hX, x0, xL, xR, wx);
    c.cd(2); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13); drawOne(hY, y0, yL, yR, wy);
    TString pdf = sumDir + "/module_center.pdf";
    PrintClean(&c, pdf);                       // PDF (archival)
    c.Print(sumDir + "/module_center.png");    // PNG = exact canvas pixels (report uses this)
    std::printf("[moduleCenter] wrote %s (+ .png)\n", pdf.Data());

    // ── Persist ────────────────────────────────────────────────────────────────
    TFile fo(sumDir + "/module_center.root", "RECREATE");
    TParameter<double>("module_center_x", x0).Write();   // Sigma HG (primary)
    TParameter<double>("module_center_y", y0).Write();
    TParameter<double>("module_width_x",  wx).Write();
    TParameter<double>("module_width_y",  wy).Write();
    TParameter<double>("module_center_x_lg", lx0).Write();  // Sigma LG cross-check
    TParameter<double>("module_center_y_lg", ly0).Write();
    hX->Write(); hY->Write(); gX.Write(); gY.Write();
    fo.Close();
    std::printf("[moduleCenter] wrote %s/module_center.root\n", sumDir.Data());
}

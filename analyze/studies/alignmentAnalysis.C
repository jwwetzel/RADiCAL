// ===========================================================================
// alignmentAnalysis.C  —  data-driven transverse alignment of the detectors in
// the wire-chamber track frame: RADiCAL shashlik, MCP, and the beam.
//
// Each detector's transverse centre is found with the method appropriate to its
// response shape:
//   RADiCAL : sum of the 8 high-gain capillaries (Sigma HG) — a flat-topped
//             plateau with sharp edges; centre = half-max edge midpoint (x & y).
//   MCP     : a smooth round acceptance blob — centre = FWHM midpoint (same
//             edge algorithm; the connector dip in the middle does not bias it).
//   Beam    : centroid (and RMS) of the track positions.
//
// The Pb-glass is shown for context only: it sits behind/beside the RADiCAL and
// responds to leakage (a low square at the RADiCAL footprint plus an asymmetric
// high lobe), so it has no clean geometric centre — none is claimed.
//
// Output: a 2x2 figure (RADiCAL / MCP / Pb-glass response maps + an alignment
// schematic overlaying the footprints and centres) and the centres + offsets,
// persisted to alignment.root.  Reads the reduced ntuples; no raw data.
//
// Run:  ROOT_INCLUDE_PATH=Analysis root -l -b -q 'Analysis/alignmentAnalysis.C+'
// ===========================================================================
#include "ChannelConfig.h"
#include "PlotUtils.h"
#include "RADiCALStyle.h"

#include "TChain.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TCanvas.h"
#include "TBox.h"
#include "TEllipse.h"
#include "TMarker.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TParameter.h"
#include "TFile.h"
#include "TSystem.h"
#include "TStyle.h"

#include <cmath>
#include <cstdio>

namespace {

double crossing(TProfile* h, int b1, int b2, double half)
{
    double x1 = h->GetBinCenter(b1), x2 = h->GetBinCenter(b2);
    double v1 = h->GetBinContent(b1), v2 = h->GetBinContent(b2);
    return (v2 == v1) ? x1 : x1 + (half - v1) * (x2 - x1) / (v2 - v1);
}

// Half-max edge midpoint (scan outward from the plateau/peak), returns centre
// and sets eL, eR.  Works for the RADiCAL plateau and the MCP blob alike.
double edgeCentre(TProfile* h, int minN, double& eL, double& eR)
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
    return 0.5 * (eL + eR);
}

} // namespace

void alignmentAnalysis()
{
    ApplyRADiCALStyle();
    gStyle->SetOptStat(0);
    const double seedX = kCalo_x0, seedY = kCalo_y0, band = 4.0;
    const int minN = 50;

    TChain ch("rad");
    for (int e = 0; e < kNRuns; ++e)
        ch.Add(Form("output/%s/ntuple.root", kRuns[e].label.Data()));
    if (ch.GetEntries() == 0) { std::printf("[alignment] no ntuples found.\n"); return; }

    Float_t x, y, spb, mcp, hg[8]; Bool_t ok;
    ch.SetBranchAddress("x_trk", &x); ch.SetBranchAddress("y_trk", &y);
    ch.SetBranchAddress("sum_pb", &spb); ch.SetBranchAddress(ch.GetBranch("mcp1_peak")?"mcp1_peak":"mcp_peak", &mcp);
    ch.SetBranchAddress("hg_peak", hg); ch.SetBranchAddress("wc_ok", &ok);

    // response maps (context) + projection profiles (centres)
    auto* mHG = new TProfile2D("aHG", "RADiCAL  #SigmaHG;x track (mm);y track (mm)", 80, -20, 35, 80, -25, 30);
    auto* mMC = new TProfile2D("aMC", "MCP;x track (mm);y track (mm)",              80, -20, 35, 80, -25, 30);
    auto* mPB = new TProfile2D("aPB", "Pb-glass  #Sigmapb (context);x track (mm);y track (mm)", 80, -20, 35, 80, -25, 30);
    TProfile hHGx("aHGx", "", 90, -15, 30), hHGy("aHGy", "", 90, -20, 25);
    TProfile hMCx("aMCx", "", 90, -15, 30), hMCy("aMCy", "", 90, -20, 25);
    for (int i=0;i<3;i++) { TProfile2D* m = (i==0?mHG:i==1?mMC:mPB); m->SetDirectory(nullptr);
                            m->GetZaxis()->SetTitle("mean amplitude (mV)"); }

    double bx=0, by=0, bxx=0, byy=0, bn=0;
    const Long64_t N = ch.GetEntries();
    for (Long64_t i = 0; i < N; ++i) {
        ch.GetEntry(i);
        if (!ok) continue;
        double shg = 0.; for (int k = 0; k < 8; ++k) shg += hg[k];
        mHG->Fill(x, y, shg); mMC->Fill(x, y, mcp); mPB->Fill(x, y, spb);
        if (std::fabs(y - seedY) < band) { hHGx.Fill(x, shg); hMCx.Fill(x, mcp); }
        if (std::fabs(x - seedX) < band) { hHGy.Fill(y, shg); hMCy.Fill(y, mcp); }
        bx += x; by += y; bxx += x*x; byy += y*y; bn += 1;
    }

    double rxL, rxR, ryL, ryR, mxL, mxR, myL, myR;
    const double radX = edgeCentre(&hHGx, minN, rxL, rxR);
    const double radY = edgeCentre(&hHGy, minN, ryL, ryR);
    const double mcX  = edgeCentre(&hMCx, minN, mxL, mxR);
    const double mcY  = edgeCentre(&hMCy, minN, myL, myR);
    const double beamX = bx/bn, beamY = by/bn;
    const double beamSx = std::sqrt(std::max(0., bxx/bn - beamX*beamX));
    const double beamSy = std::sqrt(std::max(0., byy/bn - beamY*beamY));
    const double dMC_x = mcX - radX,  dMC_y = mcY - radY;
    const double dB_x  = beamX - radX, dB_y = beamY - radY;

    std::printf("\n  Transverse alignment (WC track frame):\n");
    std::printf("    RADiCAL (HG edges) : (%.2f, %.2f) mm  footprint %.1f x %.1f mm\n",
                radX, radY, rxR-rxL, ryR-ryL);
    std::printf("    MCP    (FWHM mid)  : (%.2f, %.2f) mm  FWHM %.1f x %.1f mm\n",
                mcX, mcY, mxR-mxL, myR-myL);
    std::printf("    Beam   (centroid)  : (%.2f, %.2f) mm  sigma %.1f x %.1f mm\n",
                beamX, beamY, beamSx, beamSy);
    std::printf("    offset MCP-RADiCAL : (%+.2f, %+.2f) mm  -> %.2f mm\n",
                dMC_x, dMC_y, std::hypot(dMC_x, dMC_y));
    std::printf("    offset beam-RADiCAL: (%+.2f, %+.2f) mm  -> %.2f mm\n\n",
                dB_x, dB_y, std::hypot(dB_x, dB_y));

    // ── Figure: 2x2 (three response maps + alignment schematic) ────────────────
    TString sumDir = "output/Summary"; gSystem->mkdir(sumDir, kTRUE);
    TCanvas c("c_align", "", 1500, 1300);
    c.Divide(2, 2, 0.005, 0.005);
    auto drawMap = [&](int pad, TProfile2D* m, double cx, double cy, bool mark){
        c.cd(pad); gPad->SetRightMargin(0.15); gPad->SetLeftMargin(0.13); gPad->SetTicks(1,1);
        m->SetMinimum(0.); m->Draw("COLZ");
        if (mark) { TMarker* mk = new TMarker(cx, cy, 29); mk->SetMarkerColor(kRRed);
                    mk->SetMarkerSize(2.2); mk->Draw(); }
    };
    drawMap(1, mHG, radX, radY, true);
    drawMap(2, mMC, mcX,  mcY,  true);
    drawMap(3, mPB, 0, 0, false);

    // schematic
    c.cd(4); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.12); gPad->SetTicks(1,1);
    TH2F* fr = new TH2F("aFrame", "Detector alignment (WC frame);x track (mm);y track (mm)",
                        10, -8, 24, 10, -10, 18);
    fr->SetDirectory(nullptr); fr->Draw();
    TBox* rb = new TBox(rxL, ryL, rxR, ryR);                 // RADiCAL footprint
    rb->SetLineColor(kRBlue); rb->SetLineWidth(3); rb->SetFillStyle(0); rb->Draw("l");
    TEllipse* me = new TEllipse(mcX, mcY, 0.5*(mxR-mxL), 0.5*(myR-myL));  // MCP FWHM
    me->SetLineColor(kRGreen); me->SetLineWidth(3); me->SetFillStyle(0); me->Draw();
    TEllipse* be = new TEllipse(beamX, beamY, beamSx, beamSy);            // beam 1-sigma
    be->SetLineColor(kROrange); be->SetLineWidth(3); be->SetLineStyle(2); be->SetFillStyle(0); be->Draw();
    auto cmark = [&](double cx, double cy, int col){ TMarker* m = new TMarker(cx, cy, 29);
                  m->SetMarkerColor(col); m->SetMarkerSize(2.0); m->Draw(); };
    cmark(radX, radY, kRBlue); cmark(mcX, mcY, kRGreen); cmark(beamX, beamY, kROrange);
    TLegend* lg = new TLegend(0.58, 0.74, 0.97, 0.93); lg->SetBorderSize(0); lg->SetTextSize(0.034);
    lg->AddEntry(rb, "RADiCAL (edges)", "l");
    lg->AddEntry(me, "MCP (FWHM)", "l");
    lg->AddEntry(be, "beam (1#sigma)", "l");
    lg->Draw();
    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.032);
    tx.DrawLatex(0.16, 0.27, Form("MCP#minusRADiCAL: %.2f mm", std::hypot(dMC_x, dMC_y)));
    tx.DrawLatex(0.16, 0.21, Form("beam#minusRADiCAL: %.2f mm", std::hypot(dB_x, dB_y)));

    TString png = sumDir + "/alignment.png";
    c.Print(png);
    c.Print(sumDir + "/alignment.pdf");
    std::printf("[alignment] wrote %s (+ .pdf)\n", png.Data());

    // ── Persist ────────────────────────────────────────────────────────────────
    TFile fo(sumDir + "/alignment.root", "RECREATE");
    TParameter<double>("rad_center_x", radX).Write();  TParameter<double>("rad_center_y", radY).Write();
    TParameter<double>("mcp_center_x", mcX ).Write();  TParameter<double>("mcp_center_y", mcY ).Write();
    TParameter<double>("beam_center_x", beamX).Write(); TParameter<double>("beam_center_y", beamY).Write();
    TParameter<double>("off_mcp_rad",  std::hypot(dMC_x, dMC_y)).Write();
    TParameter<double>("off_beam_rad", std::hypot(dB_x,  dB_y )).Write();
    mHG->Write(); mMC->Write(); mPB->Write();
    fo.Close();
    std::printf("[alignment] wrote %s/alignment.root\n", sumDir.Data());
}

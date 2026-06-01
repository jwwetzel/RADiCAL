// ============================================================================
// compareEnergies.C — cross-energy collated plots for RADiCAL analysis
// ============================================================================
//
// Pure reader: opens processed ntuples (Analysis/Output/<E>/ntuple.root) and
// summary ROOT files in Analysis/Output/Summary/.  Never reads raw DRS4 data.
//
// Produces four multi-page PDFs in Analysis/Output/Summary/:
//
//   cross_energy_beam_quality.pdf    — Group 1: Beam & selection purity
//   cross_energy_containment.pdf     — Group 2: Shower containment
//   cross_energy_channels.pdf        — Group 3: Channel performance
//   cross_energy_timing_methods.pdf  — Group 4: Timing methods comparison
//
// Each PDF has 4 pages.  The layout is: 2×3 grid for per-energy panels,
// single overlaid plots for cross-energy comparisons.
//
// Usage (from repository root):
//   root -l -b -q 'Analysis/compareEnergies.C+'
// ============================================================================

#include "ChannelConfig.h"   // kRuns, kNRuns, kCap, kNCap, kCalo_x0/y0
#include "PlotUtils.h"       // StylePad, FitGaussCore, ScanRunCenters

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TArrow.h"
#include "TEllipse.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TAxis.h"
#include "TF1.h"
#include "TPad.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

// ---------------------------------------------------------------------------
// Colour and marker scheme for the 6 beam energies (25 → 150 GeV).
// Colors from RADiCALStyle.h (kREnergyCols, set by ApplyRADiCALStyle).
// ---------------------------------------------------------------------------
// kREnergyCols[6] and kRChannelCols[8] are provided by RADiCALStyle.h.
static const int kEMark[6] = { 20, 21, 22, 23, 24, 25 };

static const char* kSumDir = "Analysis/Output/Summary/";

// ---------------------------------------------------------------------------
// OpenNtuple — open the processed ntuple for energy index r.
// Returns nullptr (with message) if the file or tree is missing.
// Caller owns the returned TFile* and must Close()+delete it.
// ---------------------------------------------------------------------------
static TFile* OpenNtuple(int r, TTree*& tree)
{
    TString path = Form("Analysis/Output/%s/ntuple.root",
                        kRuns[r].label.Data());
    TFile* f = TFile::Open(path);
    if (!f || f->IsZombie()) {
        std::cerr << "compareEnergies: cannot open " << path << "\n";
        tree = nullptr;
        return nullptr;
    }
    tree = static_cast<TTree*>(f->Get("rad"));
    if (!tree) {
        std::cerr << "compareEnergies: tree 'rad' missing in " << path << "\n";
        f->Close();
        delete f;
        tree = nullptr;
        return nullptr;
    }
    return f;
}

// ---------------------------------------------------------------------------
// NormToUnity — scale h so its integral == 1.0 (in-place).
// Enables shape comparison between energies with different event counts.
// ---------------------------------------------------------------------------
static void NormToUnity(TH1F* h)
{
    double integral = h->Integral();
    if (integral > 0.) h->Scale(1.0 / integral);
}

// ---------------------------------------------------------------------------
// ApplyLineStyle — set color, width, marker style for a graph or histo.
// Called once per energy index r to enforce the standard scheme.
// ---------------------------------------------------------------------------
static void ApplyLineStyle(TH1F* h, int r)
{
    h->SetLineColor(kREnergyCols[r]);
    h->SetLineWidth(2);
    h->SetMarkerColor(kREnergyCols[r]);
    h->SetMarkerStyle(kEMark[r]);
}

static void ApplyGraphStyle(TGraph* g, int r)
{
    g->SetLineColor(kREnergyCols[r]);
    g->SetLineWidth(2);
    g->SetMarkerColor(kREnergyCols[r]);
    g->SetMarkerStyle(kEMark[r]);
    g->SetMarkerSize(1.0);
}

// ---------------------------------------------------------------------------
// DrawEnergyLegend — add one entry per energy to a legend.
// ---------------------------------------------------------------------------
static void DrawEnergyLegend(TObject* obj[6], bool valid[6],
                              double x1, double y1, double x2, double y2,
                              const char* opt = "lp")
{
    TLegend* leg = new TLegend(x1, y1, x2, y2);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.040);
    for (int r = 0; r < kNRuns; ++r)
        if (valid[r])
            leg->AddEntry(obj[r], Form("%.0f GeV", kRuns[r].energy_GeV), opt);
    leg->Draw();
    // Owned by current pad; do not delete here.
}

// ---------------------------------------------------------------------------
// PageTitle — draw a centred title at the top of the current canvas.
// ---------------------------------------------------------------------------
static void PageTitle(const char* title)
{
    TLatex t;
    t.SetNDC();
    t.SetTextSize(0.030);
    t.SetTextAlign(23);   // centre-top: text hangs BELOW the baseline, not above the canvas edge
    t.DrawLatex(0.50, 0.992, title);
}

// ===========================================================================
// Group 1 — Beam & Selection Purity
//
//   Page 1: Hit maps (2×3 grid) — is the beam centred consistently at all E?
//   Page 2: MCP1 amplitude (overlaid, normalised) — saturation trend vs E
//   Page 3: ΣLG energy distributions (overlaid, normalised) — linearity check
//   Page 4: Fiducial efficiency vs radius — where the tail is relative to cuts
// ===========================================================================
static void PlotGroup1_BeamQuality()
{
    std::cout << "compareEnergies: Group 1 — Beam & Selection Purity\n";

    // Histogram ranges
    const double kHitXlo = -2.,  kHitXhi = 16.;  // WC x range [mm]
    const double kHitYlo = -4.,  kHitYhi = 14.;  // WC y range [mm]
    const double kMCPMax = 1100.;                 // MCP amplitude range [mV]
    // sum_lg common axis: 55 mV/GeV × 150 GeV covers ≈ 3σ above peak at all energies
    const double kLGMax  = 55. * 150.;

    // Per-energy beam radius scan (WC-OK + MCP-quality events)
    struct BeamPt { float r_beam; };
    std::vector<BeamPt> scanPts[6];

    TH2F* hHit[6]    = {};
    TH1F* hMCPAmp[6] = {};
    TH1F* hSumLG[6]  = {};
    double xCen[6] = {}, yCen[6] = {};

    // ── Fill loop ────────────────────────────────────────────────────────────
    for (int r = 0; r < kNRuns; ++r) {
        TTree* t = nullptr;
        TFile* fin = OpenNtuple(r, t);
        if (!fin) continue;

        // Derive beam centroid for this run (used for r_beam calculation)
        double xc, yc, tcfd, trms;
        ScanRunCenters(t, xc, yc, tcfd, trms);
        xCen[r] = xc;
        yCen[r] = yc;

        TString tag = kRuns[r].label;
        // 1 mm/bin matches WC delay-line wire pitch (kWC_resBin from SelectionCuts.h)
        int nHitX = static_cast<int>(std::round((kHitXhi - kHitXlo) / kWC_resBin));
        int nHitY = static_cast<int>(std::round((kHitYhi - kHitYlo) / kWC_resBin));
        hHit[r] = new TH2F(Form("hHit_%s",    tag.Data()), "",
                            nHitX, kHitXlo, kHitXhi, nHitY, kHitYlo, kHitYhi);
        hHit[r]->SetDirectory(nullptr);
        hMCPAmp[r] = new TH1F(Form("hMCPAmp_%s", tag.Data()), "",
                               110, 0., kMCPMax);
        hMCPAmp[r]->SetDirectory(nullptr);
        hSumLG[r]  = new TH1F(Form("hSumLG_%s",  tag.Data()), "",
                               150, 0., kLGMax);
        hSumLG[r]->SetDirectory(nullptr);

        Float_t x_trk, y_trk, mcp_peak, sum_lg;
        Bool_t  wc_ok;
        t->SetBranchAddress("x_trk",    &x_trk);
        t->SetBranchAddress("y_trk",    &y_trk);
        t->SetBranchAddress("wc_ok",    &wc_ok);
        t->SetBranchAddress("mcp_peak", &mcp_peak);
        t->SetBranchAddress("sum_lg",   &sum_lg);

        Long64_t nEv = t->GetEntries();
        for (Long64_t ev = 0; ev < nEv; ++ev) {
            t->GetEntry(ev);
            if (!wc_ok) continue;

            hHit[r]->Fill(x_trk, y_trk);
            hMCPAmp[r]->Fill(mcp_peak);

            float dx     = x_trk - static_cast<float>(xc);
            float dy     = y_trk - static_cast<float>(yc);
            float r_beam = std::sqrt(dx*dx + dy*dy);

            // For the efficiency scan: keep events with clean MCP reference
            bool mcp_ok = (mcp_peak > kMCP1_minPeak && mcp_peak < kMCP1_maxPeak);
            if (mcp_ok) {
                BeamPt bp; bp.r_beam = r_beam;
                scanPts[r].push_back(bp);

                // Energy distribution: timing-fiducial events only
                if (r_beam < static_cast<float>(kFiducial_r_timing))
                    hSumLG[r]->Fill(sum_lg);
            }
        }
        t->ResetBranchAddresses();
        fin->Close();
        delete fin;
        std::cout << "  " << kRuns[r].label
                  << ": " << hHit[r]->GetEntries() << " WC-OK  "
                  << scanPts[r].size() << " MCP-ok\n";
    }

    TString outPDF = Form("%scross_energy_beam_quality.pdf", kSumDir);
    TCanvas c("c_bq", "", 1200, 800);

    // ── Page 1: Hit maps (2×3 grid) ─────────────────────────────────────────
    c.Divide(3, 2, 0.003, 0.035);  // 3.5% top/bottom margin reserves space for PageTitle
    for (int r = 0; r < kNRuns; ++r) {
        c.cd(r + 1);
        gPad->SetLeftMargin(0.11); gPad->SetBottomMargin(0.11);
        gPad->SetRightMargin(0.14); gPad->SetTopMargin(0.10);
        gPad->SetTickx(1); gPad->SetTicky(1);
        gPad->SetLogz();   // log z reveals the beam spot at low-stat energies
        if (!hHit[r]) continue;

        hHit[r]->GetXaxis()->SetTitle("x_{trk} (mm)");
        hHit[r]->GetYaxis()->SetTitle("y_{trk} (mm)");
        hHit[r]->GetXaxis()->SetTitleSize(0.07); hHit[r]->GetXaxis()->SetLabelSize(0.07);
        hHit[r]->GetYaxis()->SetTitleSize(0.07); hHit[r]->GetYaxis()->SetLabelSize(0.07);
        hHit[r]->GetYaxis()->SetTitleOffset(0.85);
        hHit[r]->SetMinimum(0.5);   // suppress zero-count bins (rendered as white)
        hHit[r]->Draw("COLZ");

        // Timing fiducial circle (orange)
        TEllipse* e3 = new TEllipse(xCen[r], yCen[r],
                                    kFiducial_r_timing, kFiducial_r_timing);
        e3->SetLineColor(kOrange+1); e3->SetLineWidth(2);
        e3->SetFillStyle(0); e3->Draw("SAME");

        // Energy fiducial circle (dashed red, tighter)
        TEllipse* e2 = new TEllipse(xCen[r], yCen[r],
                                    kFiducial_r_energy, kFiducial_r_energy);
        e2->SetLineColor(kRed+1); e2->SetLineStyle(2); e2->SetLineWidth(2);
        e2->SetFillStyle(0); e2->Draw("SAME");

        // Labels placed at y=0.80/0.71 — clear of the top axis tick at y≈0.90
        TLatex lab; lab.SetNDC(); lab.SetTextSize(0.075); lab.SetTextFont(72);
        lab.DrawLatex(0.16, 0.80, Form("%.0f GeV", kRuns[r].energy_GeV));
        lab.SetTextFont(42); lab.SetTextSize(0.062);
        lab.DrawLatex(0.16, 0.71,
                      Form("N_{WC} = %.0fk", hHit[r]->GetEntries()/1000.));
    }
    c.cd(0);
    PageTitle("Beam hit maps (WC-OK) | orange: timing fiducial r=3 mm | red dashed: energy fiducial r=2 mm");
    c.Print(outPDF + "(");

    // ── Page 2: MCP1 amplitude distributions (overlaid, normalised) ─────────
    c.Clear(); c.cd(); StylePad();
    double yMaxMCP = 0.;
    for (int r = 0; r < kNRuns; ++r) {
        if (!hMCPAmp[r]) continue;
        NormToUnity(hMCPAmp[r]);
        yMaxMCP = std::max(yMaxMCP, hMCPAmp[r]->GetMaximum());
    }
    bool first = true;
    TObject* mObj[6] = {};
    bool     mVal[6] = {};
    for (int r = 0; r < kNRuns; ++r) {
        if (!hMCPAmp[r]) continue;
        ApplyLineStyle(hMCPAmp[r], r);
        hMCPAmp[r]->GetXaxis()->SetTitle("MCP1 amplitude (mV)");
        hMCPAmp[r]->GetYaxis()->SetTitle("Events (normalised)");
        hMCPAmp[r]->GetYaxis()->SetRangeUser(0., 1.18 * yMaxMCP);
        if (first) { hMCPAmp[r]->Draw("HIST"); first = false; }
        else          hMCPAmp[r]->Draw("HIST SAME");
        mObj[r] = hMCPAmp[r]; mVal[r] = true;
    }
    TLine lLo(kMCP1_minPeak, 0., kMCP1_minPeak, 1.05*yMaxMCP);
    lLo.SetLineColor(kGray+2); lLo.SetLineStyle(2); lLo.SetLineWidth(2); lLo.Draw();
    TLine lHi(kMCP1_maxPeak, 0., kMCP1_maxPeak, 1.05*yMaxMCP);
    lHi.SetLineColor(kOrange+1); lHi.SetLineStyle(2); lHi.SetLineWidth(2); lHi.Draw();
    {
        TLatex ann; ann.SetNDC(); ann.SetTextSize(0.040);
        ann.SetTextColor(kGray+2);
        ann.DrawLatex(0.18, 0.76, Form("lower cut: %.0f mV", (double)kMCP1_minPeak));
        ann.SetTextColor(kOrange+1);
        ann.DrawLatex(0.18, 0.69, Form("sat. cut:  %.0f mV", (double)kMCP1_maxPeak));
    }
    DrawEnergyLegend(mObj, mVal, 0.63, 0.55, 0.93, 0.88);
    PageTitle("MCP1 amplitude distributions (WC-OK events, normalised to unit area)");
    c.Print(outPDF);

    // ── Page 3: ΣLG energy distributions (overlaid, normalised) ─────────────
    c.Clear(); c.cd(); StylePad();
    double yMaxLG = 0.;
    for (int r = 0; r < kNRuns; ++r) {
        if (!hSumLG[r]) continue;
        NormToUnity(hSumLG[r]);
        yMaxLG = std::max(yMaxLG, hSumLG[r]->GetMaximum());
    }
    first = true;
    TObject* lgObj[6] = {};
    bool     lgVal[6] = {};
    for (int r = 0; r < kNRuns; ++r) {
        if (!hSumLG[r]) continue;
        ApplyLineStyle(hSumLG[r], r);
        hSumLG[r]->GetXaxis()->SetTitle("#Sigma LG (mV)");
        hSumLG[r]->GetYaxis()->SetTitle("Events (normalised)");
        hSumLG[r]->GetYaxis()->SetRangeUser(0., 1.18 * yMaxLG);
        if (first) { hSumLG[r]->Draw("HIST"); first = false; }
        else          hSumLG[r]->Draw("HIST SAME");
        lgObj[r] = hSumLG[r]; lgVal[r] = true;
    }
    DrawEnergyLegend(lgObj, lgVal, 0.63, 0.55, 0.93, 0.88);
    PageTitle("#Sigma LG signal (WC-OK + MCP-quality + r < 3 mm)");
    c.Print(outPDF);

    // ── Page 4: Fiducial efficiency vs radius cut ────────────────────────────
    // Shows what fraction of good beam events survive each fiducial tightness.
    // Denominator: WC-OK + MCP-quality events (no fiducial cut applied yet).
    c.Clear(); c.cd(); StylePad();

    const int    kNrPts  = 20;
    const double kRstep  = 0.30;   // mm per step
    const double kRstart = 0.30;   // mm, first scan point

    TGraph*    gEff[6] = {};
    TObject* effObj[6] = {};
    bool     effVal[6] = {};
    for (int r = 0; r < kNRuns; ++r) {
        if (scanPts[r].empty()) continue;
        long nTot = static_cast<long>(scanPts[r].size());
        gEff[r] = new TGraph(kNrPts);
        for (int ip = 0; ip < kNrPts; ++ip) {
            double rThresh = kRstart + ip * kRstep;
            long nPass = 0;
            for (size_t ev = 0; ev < scanPts[r].size(); ++ev)
                if (scanPts[r][ev].r_beam < static_cast<float>(rThresh)) ++nPass;
            gEff[r]->SetPoint(ip, rThresh, static_cast<double>(nPass) / nTot);
        }
        ApplyGraphStyle(gEff[r], r);
        effObj[r] = gEff[r]; effVal[r] = true;
    }

    TMultiGraph* mgEff = new TMultiGraph();
    for (int r = 0; r < kNRuns; ++r)
        if (gEff[r]) mgEff->Add(gEff[r], "LP");
    mgEff->Draw("A");
    mgEff->GetXaxis()->SetTitle("Radius cut r_{cut} (mm)");
    mgEff->GetYaxis()->SetTitle("Fraction of MCP-quality events passing r < r_{cut}");
    mgEff->GetYaxis()->SetRangeUser(0., 1.08);
    mgEff->GetXaxis()->SetLimits(0., kRstart + (kNrPts-1)*kRstep + 0.3);

    TLine lTim(kFiducial_r_timing, 0., kFiducial_r_timing, 1.04);
    lTim.SetLineColor(kOrange+1); lTim.SetLineStyle(2); lTim.SetLineWidth(2);
    lTim.Draw("SAME");
    TLine lEn(kFiducial_r_energy, 0., kFiducial_r_energy, 1.04);
    lEn.SetLineColor(kRed+1); lEn.SetLineStyle(2); lEn.SetLineWidth(2);
    lEn.Draw("SAME");
    {
        TLatex ann; ann.SetNDC(); ann.SetTextSize(0.038);
        ann.SetTextColor(kRed+1);
        ann.DrawLatex(0.22, 0.35, Form("energy fid: %.1f mm", kFiducial_r_energy));
        ann.SetTextColor(kOrange+1);
        ann.DrawLatex(0.22, 0.28, Form("timing fid: %.1f mm", kFiducial_r_timing));
    }
    DrawEnergyLegend(effObj, effVal, 0.60, 0.25, 0.92, 0.58);
    PageTitle("Fiducial efficiency vs radius cut (WC-OK + MCP-quality events)");
    c.Print(outPDF + ")");

    // Cleanup.
    // TMultiGraph destructor calls SetOwner(kTRUE) on its internal TList before
    // deleting it, so deleting mgEff also deletes all contained gEff[r] graphs.
    // Delete mgEff first; do NOT delete gEff[r] separately (double-free).
    delete mgEff;   // owns and deletes gEff[0..5]
    for (int r = 0; r < kNRuns; ++r) {
        delete hHit[r]; delete hMCPAmp[r]; delete hSumLG[r];
    }
    std::cout << "compareEnergies: Group 1 -> " << outPDF << "\n";
}

// ===========================================================================
// Group 2 — Shower Containment
//
//   Page 1: sum_pb / sum_lg ratio (overlaid) — hadronic fraction vs energy
//   Page 2: PbGlass scatter 6-panel — two-population structure at each energy
//   Page 3: Contained fraction vs beam radius — edge-shower effect vs energy
//   Page 4: Contained fraction vs energy — headline containment number
// ===========================================================================
static void PlotGroup2_Containment()
{
    std::cout << "compareEnergies: Group 2 — Shower Containment\n";

    // r-bin scan for pages 3 and 4
    const int    kNrBins = 10;
    const double kRbinW  = 0.5;   // mm per bin

    TH1F* hRatio[6]   = {};       // sum_pb / sum_lg for in-fiducial events
    TH2F* hScatter[6] = {};       // sum_pb vs sum_lg scatter

    // Contained/total counts per r-bin (for page 3)
    long nInBin[6][kNrBins]   = {};  // events with good LG signal in this r-bin
    long nContBin[6][kNrBins] = {};  // above + passing containment cut

    // Timing-fiducial totals for page 4
    long nFidTot[6]  = {};
    long nFidCont[6] = {};

    double xCen[6] = {}, yCen[6] = {};

    for (int r = 0; r < kNRuns; ++r) {
        TTree* t = nullptr;
        TFile* fin = OpenNtuple(r, t);
        if (!fin) continue;

        double xc, yc, tcfd, trms;
        ScanRunCenters(t, xc, yc, tcfd, trms);
        xCen[r] = xc; yCen[r] = yc;

        // Energy-derived scatter axis: 30 mV/GeV × E_GeV, plus headroom
        double xmax_lg = 55. * kRuns[r].energy_GeV;
        double ymax_pb = 0.50 * xmax_lg;

        TString tag = kRuns[r].label;
        hRatio[r] = new TH1F(Form("hRatio_%s",   tag.Data()), "",
                              100, 0., 1.2);
        hRatio[r]->SetDirectory(nullptr);
        hScatter[r] = new TH2F(Form("hScatter_%s", tag.Data()), "",
                                80, 0., xmax_lg, 60, 0., ymax_pb);
        hScatter[r]->SetDirectory(nullptr);

        Float_t x_trk, y_trk, mcp_peak, sum_lg, sum_pb;
        Bool_t  wc_ok;
        t->SetBranchAddress("x_trk",    &x_trk);
        t->SetBranchAddress("y_trk",    &y_trk);
        t->SetBranchAddress("wc_ok",    &wc_ok);
        t->SetBranchAddress("mcp_peak", &mcp_peak);
        t->SetBranchAddress("sum_lg",   &sum_lg);
        t->SetBranchAddress("sum_pb",   &sum_pb);

        Long64_t nEv = t->GetEntries();
        for (Long64_t ev = 0; ev < nEv; ++ev) {
            t->GetEntry(ev);
            if (!wc_ok) continue;
            if (mcp_peak < kMCP_minPeak_E) continue;   // loose MCP quality

            float dx     = x_trk - static_cast<float>(xc);
            float dy     = y_trk - static_cast<float>(yc);
            float r_beam = std::sqrt(dx*dx + dy*dy);

            // Scatter: all events near fiducial (for population visibility)
            if (r_beam < static_cast<float>(kFiducial_r_timing))
                hScatter[r]->Fill(sum_lg, sum_pb);

            // Ratio + r-bin analysis: need a genuine LG signal
            bool has_signal = (sum_lg > kSumLG_centroid);
            if (!has_signal) continue;

            // Containment ratio plot (timing fiducial events)
            if (r_beam < static_cast<float>(kFiducial_r_timing))
                hRatio[r]->Fill(sum_pb / sum_lg);

            // r-bin scan: which radial shell are we in?
            int ib = static_cast<int>(r_beam / kRbinW);
            if (ib >= 0 && ib < kNrBins) {
                ++nInBin[r][ib];
                if (sum_pb < kPb_maxRatio * sum_lg) ++nContBin[r][ib];
            }

            // Timing-fiducial total for page 4
            if (r_beam < static_cast<float>(kFiducial_r_timing)) {
                ++nFidTot[r];
                if (sum_pb < kPb_maxRatio * sum_lg) ++nFidCont[r];
            }
        }
        t->ResetBranchAddresses();
        fin->Close();
        delete fin;
        std::cout << "  " << kRuns[r].label << ": containment "
                  << Form("%.1f%%", nFidTot[r] > 0
                          ? 100.*nFidCont[r]/nFidTot[r] : 0.) << "\n";
    }

    TString outPDF = Form("%scross_energy_containment.pdf", kSumDir);
    TCanvas c("c_cont", "", 1200, 800);

    // ── Page 1: Containment ratio distributions (overlaid, normalised) ───────
    c.Clear(); c.cd(); StylePad();
    double yMaxR = 0.;
    for (int r = 0; r < kNRuns; ++r) {
        if (!hRatio[r]) continue;
        NormToUnity(hRatio[r]);
        yMaxR = std::max(yMaxR, hRatio[r]->GetMaximum());
    }
    bool first = true;
    TObject* rObj[6] = {};
    bool     rVal[6] = {};
    for (int r = 0; r < kNRuns; ++r) {
        if (!hRatio[r]) continue;
        ApplyLineStyle(hRatio[r], r);
        hRatio[r]->GetXaxis()->SetTitle("#Sigma PbGlass / #Sigma LG");
        hRatio[r]->GetYaxis()->SetTitle("Events (normalised)");
        // Zoom to the populated range: the ratio peaks at ~0 and the cut is at
        // 0.30; the original 0-1.2 axis left ~75% empty.  (Normalised tail
        // beyond 0.5 is negligible.)
        hRatio[r]->GetXaxis()->SetRangeUser(0., 0.5);
        hRatio[r]->GetYaxis()->SetRangeUser(0., 1.20 * yMaxR);
        if (first) { hRatio[r]->Draw("HIST"); first = false; }
        else          hRatio[r]->Draw("HIST SAME");
        rObj[r] = hRatio[r]; rVal[r] = true;
    }
    TLine lCut(kPb_maxRatio, 0., kPb_maxRatio, 1.10*yMaxR);
    lCut.SetLineColor(kRed+1); lCut.SetLineStyle(2); lCut.SetLineWidth(2);
    lCut.Draw("SAME");
    {
        TLatex ann; ann.SetNDC(); ann.SetTextSize(0.038); ann.SetTextColor(kRed+1);
        ann.DrawLatex(0.60, 0.84,
                      Form("containment cut: %.2f", (double)kPb_maxRatio));
    }
    DrawEnergyLegend(rObj, rVal, 0.63, 0.55, 0.93, 0.86);
    PageTitle("#Sigma_{PbGlass} / #Sigma_{LG} ratio (in timing fiducial, normalised)");
    c.Print(outPDF + "(");

    // ── Page 2: PbGlass scatter 6-panel (2×3 grid) ──────────────────────────
    c.Clear();
    c.Divide(3, 2, 0.003, 0.035);  // 3.5% top/bottom margin reserves space for PageTitle
    gStyle->SetPalette(kCherry); TColor::InvertPalette();
    for (int r = 0; r < kNRuns; ++r) {
        c.cd(r + 1);
        gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12);
        gPad->SetRightMargin(0.14); gPad->SetTopMargin(0.10);
        gPad->SetTickx(1); gPad->SetTicky(1);
        gPad->SetLogz();
        if (!hScatter[r]) continue;

        hScatter[r]->GetXaxis()->SetTitle("#Sigma LG (mV)");
        hScatter[r]->GetYaxis()->SetTitle("#Sigma PbGlass (mV)");
        hScatter[r]->GetXaxis()->SetTitleSize(0.07); hScatter[r]->GetXaxis()->SetLabelSize(0.06);
        hScatter[r]->GetYaxis()->SetTitleSize(0.07); hScatter[r]->GetYaxis()->SetLabelSize(0.06);
        hScatter[r]->GetYaxis()->SetTitleOffset(0.90);
        hScatter[r]->Draw("COLZ");

        // Containment cut line: sum_pb = kPb_maxRatio * sum_lg
        double xhi = hScatter[r]->GetXaxis()->GetXmax();
        TLine* lc = new TLine(0., 0., xhi, kPb_maxRatio * xhi);
        lc->SetLineColor(kRed+1); lc->SetLineStyle(2); lc->SetLineWidth(2);
        lc->Draw("SAME");

        TLatex lab; lab.SetNDC(); lab.SetTextSize(0.08); lab.SetTextFont(72);
        lab.DrawLatex(0.17, 0.85, Form("%.0f GeV", kRuns[r].energy_GeV));
    }
    c.cd(0);
    PageTitle("#Sigma_{PbGlass} vs #Sigma_{LG} scatter (timing fiducial) | dashed: 30% containment cut");
    c.Print(outPDF);

    // ── Page 3: Contained fraction vs beam radius (6 energy curves) ─────────
    // Reveals whether edge showers have more hadronic punch-through.
    c.Clear(); c.cd(); StylePad();
    TGraph* gContR[6] = {};
    TObject* crObj[6] = {};
    bool     crVal[6] = {};
    for (int r = 0; r < kNRuns; ++r) {
        gContR[r] = new TGraph();
        for (int ib = 0; ib < kNrBins; ++ib) {
            if (nInBin[r][ib] < 20) continue;  // too few events in this bin
            double rCenter = (ib + 0.5) * kRbinW;
            double frac    = static_cast<double>(nContBin[r][ib]) / nInBin[r][ib];
            gContR[r]->SetPoint(gContR[r]->GetN(), rCenter, frac);
        }
        ApplyGraphStyle(gContR[r], r);
        crObj[r] = gContR[r]; crVal[r] = (gContR[r]->GetN() > 0);
    }
    TMultiGraph* mgCR = new TMultiGraph();
    for (int r = 0; r < kNRuns; ++r)
        if (crVal[r]) mgCR->Add(gContR[r], "LP");
    mgCR->Draw("A");
    mgCR->GetXaxis()->SetTitle("Beam radius r (mm)");
    mgCR->GetYaxis()->SetTitle(Form("Fraction with #Sigma_{PbGlass} < %.0f%% #times #Sigma_{LG}",
                                     100.*kPb_maxRatio));
    mgCR->GetYaxis()->SetRangeUser(0.4, 1.05);
    TLine lFid(kFiducial_r_timing, 0.40, kFiducial_r_timing, 1.02);
    lFid.SetLineColor(kOrange+1); lFid.SetLineStyle(2); lFid.SetLineWidth(2);
    lFid.Draw("SAME");
    DrawEnergyLegend(crObj, crVal, 0.15, 0.20, 0.50, 0.53);
    PageTitle("Shower containment fraction vs beam radius (events with #Sigma_{LG} > 300 mV)");
    c.Print(outPDF);

    // ── Page 4: Contained fraction vs beam energy (headline number) ──────────
    c.Clear(); c.cd(); StylePad();
    TGraphErrors* gCont = new TGraphErrors();
    for (int r = 0; r < kNRuns; ++r) {
        if (nFidTot[r] < 10) continue;
        double frac = static_cast<double>(nFidCont[r]) / nFidTot[r];
        // Binomial error on fraction
        double err  = std::sqrt(frac * (1. - frac) / nFidTot[r]);
        int np = gCont->GetN();
        gCont->SetPoint(np, kRuns[r].energy_GeV, 100. * frac);
        gCont->SetPointError(np, 0., 100. * err);
    }
    gCont->SetLineColor(kBlue+1); gCont->SetLineWidth(2);
    gCont->SetMarkerColor(kBlue+1); gCont->SetMarkerStyle(21); gCont->SetMarkerSize(1.4);
    gCont->Draw("APL");
    gCont->GetXaxis()->SetTitle("Beam energy (GeV)");
    gCont->GetYaxis()->SetTitle(Form("Contained fraction (%%): #Sigma_{Pb} < %.0f%% #times #Sigma_{LG}",
                                     100.*kPb_maxRatio));
    gCont->GetYaxis()->SetRangeUser(60., 105.);
    gCont->GetXaxis()->SetLimits(15., 165.);
    {
        TLatex ann; ann.SetNDC(); ann.SetTextSize(0.038); ann.SetTextColor(kGray+1);
        ann.DrawLatex(0.18, 0.20,
                      "timing fiducial (r < 3 mm) + #Sigma_{LG} > 300 mV");
    }
    PageTitle("Shower containment efficiency vs beam energy");
    c.Print(outPDF + ")");

    // Persist the headline contained-fraction graph for the data-driven report.
    {
        TFile fout(Form("%scross_energy.root", kSumDir), "RECREATE");
        gCont->Write("gContainedFrac");
        fout.Close();
        std::cout << "compareEnergies: wrote cross_energy.root (gContainedFrac)\n";
    }

    // Cleanup.  mgCR destructor owns and deletes all gContR[r] graphs.
    delete mgCR;    // owns and deletes gContR[0..5]
    delete gCont;
    for (int r = 0; r < kNRuns; ++r) {
        delete hRatio[r]; delete hScatter[r];
    }
    std::cout << "compareEnergies: Group 2 -> " << outPDF << "\n";
}

// ===========================================================================
// Group 3 — Channel Performance
//
//   Page 1: σ_t per channel vs energy — which capillary drives the resolution?
//   Page 2: Mean HG peak per channel (6-panel) — signal level vs energy
//   Page 3: Channel active fraction vs energy — stability of each capillary
//   Page 4: σ_t heat map (8 channels × 6 energies) — at-a-glance overview
// ===========================================================================
static void PlotGroup3_Channels()
{
    std::cout << "compareEnergies: Group 3 — Channel Performance\n";

    // Channel names matching the gTRes_<name> keys in summary.root
    static const char* kChName[8] = {
        "NW-D","NE-D","SE-D","SW-D","NW-U","NE-U","SE-U","SW-U"
    };

    // ── Pages 1 and 4: read per-channel timing from summary.root ─────────────
    TString sumPath = Form("%ssummary.root", kSumDir);
    TFile* fSum = TFile::Open(sumPath);
    if (!fSum || fSum->IsZombie()) {
        std::cerr << "compareEnergies: cannot open " << sumPath << "\n";
        return;
    }

    TGraphErrors* gTRes[8] = {};
    for (int ic = 0; ic < kNCap; ++ic) {
        gTRes[ic] = static_cast<TGraphErrors*>(
            fSum->Get(Form("gTRes_%s", kChName[ic])));
        if (!gTRes[ic])
            std::cerr << "  WARNING: gTRes_" << kChName[ic] << " not found\n";
    }
    TGraphErrors* gTResAvg = static_cast<TGraphErrors*>(
        fSum->Get("gTimingResolution"));

    // Build heat-map matrix: heatSig[ich][iE] = sigma_t [ps]
    double heatSig[8][6] = {};
    for (int ic = 0; ic < kNCap; ++ic) {
        if (!gTRes[ic]) continue;
        for (int ip = 0; ip < gTRes[ic]->GetN(); ++ip) {
            double E = gTRes[ic]->GetX()[ip];
            for (int ie = 0; ie < kNRuns; ++ie) {
                if (std::fabs(kRuns[ie].energy_GeV - E) < 1.) {
                    heatSig[ic][ie] = gTRes[ic]->GetY()[ip];
                    break;
                }
            }
        }
    }

    // ── Pages 2 and 3: read ntuples for HG peak data ─────────────────────────
    // nActive[r][ic]  = events with hg_peak[ic] > kHG_minPeak (in timing fiducial)
    // sumHGPeak[r][ic] = sum of hg_peak values for active events
    // nFidEvt[r]       = total timing-fiducial MCP-quality events
    long   nActive[6][8]   = {};
    double sumHGPeak[6][8] = {};
    long   nFidEvt[6]      = {};
    double xCen[6] = {}, yCen[6] = {};

    for (int r = 0; r < kNRuns; ++r) {
        TTree* t = nullptr;
        TFile* fin = OpenNtuple(r, t);
        if (!fin) continue;

        double xc, yc, tcfd, trms;
        ScanRunCenters(t, xc, yc, tcfd, trms);
        xCen[r] = xc; yCen[r] = yc;

        Float_t x_trk, y_trk, mcp_peak, hg_peak[8];
        Bool_t  wc_ok;
        t->SetBranchAddress("x_trk",    &x_trk);
        t->SetBranchAddress("y_trk",    &y_trk);
        t->SetBranchAddress("wc_ok",    &wc_ok);
        t->SetBranchAddress("mcp_peak", &mcp_peak);
        t->SetBranchAddress("hg_peak",   hg_peak);

        Long64_t nEv = t->GetEntries();
        for (Long64_t ev = 0; ev < nEv; ++ev) {
            t->GetEntry(ev);
            if (!wc_ok) continue;
            if (mcp_peak < kMCP1_minPeak || mcp_peak > kMCP1_maxPeak) continue;

            float dx     = x_trk - static_cast<float>(xc);
            float dy     = y_trk - static_cast<float>(yc);
            float r_beam = std::sqrt(dx*dx + dy*dy);
            if (r_beam >= static_cast<float>(kFiducial_r_timing)) continue;

            ++nFidEvt[r];
            for (int ic = 0; ic < kNCap; ++ic) {
                if (hg_peak[ic] > kHG_minPeak) {
                    ++nActive[r][ic];
                    sumHGPeak[r][ic] += hg_peak[ic];
                }
            }
        }
        t->ResetBranchAddresses();
        fin->Close();
        delete fin;
    }

    TString outPDF = Form("%scross_energy_channels.pdf", kSumDir);
    TCanvas c("c_ch", "", 1200, 800);

    // ── Page 1: σ_t per channel vs energy (8 lines) ──────────────────────────
    // Use sidebar (right margin = 0.27) — 9 legend entries need the space.
    c.Clear(); c.cd(); StylePad(false, true);
    double yMaxSig = 0.;
    for (int ic = 0; ic < kNCap; ++ic) {
        if (!gTRes[ic]) continue;
        for (int ip = 0; ip < gTRes[ic]->GetN(); ++ip)
            yMaxSig = std::max(yMaxSig, gTRes[ic]->GetY()[ip]);
    }
    if (yMaxSig < 10.) yMaxSig = 400.;
    TH1F* frame1 = static_cast<TH1F*>(
        c.DrawFrame(15., 50., 165., 1.30 * yMaxSig,
                    ";Beam energy (GeV);#sigma_{t} (ps)"));
    frame1->GetXaxis()->SetTitleSize(0.050);
    frame1->GetYaxis()->SetTitleSize(0.050);
    frame1->GetYaxis()->SetTitleOffset(1.30);

    TLegend* pleg1 = MakeLegend(9);
    for (int ic = 0; ic < kNCap; ++ic) {
        if (!gTRes[ic]) continue;
        gTRes[ic]->SetLineColor(kRChannelCols[ic]);
        gTRes[ic]->SetMarkerColor(kRChannelCols[ic]);
        gTRes[ic]->SetMarkerStyle(20 + ic);
        gTRes[ic]->SetMarkerSize(1.1);
        gTRes[ic]->SetLineWidth(2);
        gTRes[ic]->Draw("PL SAME");
        pleg1->AddEntry(gTRes[ic], kChName[ic], "lp");
    }
    if (gTResAvg) {
        gTResAvg->SetLineColor(kBlack); gTResAvg->SetLineStyle(2);
        gTResAvg->SetLineWidth(2); gTResAvg->SetMarkerStyle(29);
        gTResAvg->SetMarkerColor(kBlack);
        gTResAvg->Draw("PL SAME");
        pleg1->AddEntry(gTResAvg, "8-ch average", "lp");
    }
    pleg1->Draw();
    PageTitle("Per-channel CFD-5% timing resolution vs beam energy");
    c.Print(outPDF + "(");

    // ── Page 2: Mean HG peak per channel, 6-panel grid ───────────────────────
    // hBar objects must remain alive until after c.Print() — deleting a drawn
    // object before Print() removes it from the pad's primitives list, leaving
    // a blank panel.  Collect all six and delete after printing.
    TH1F* hBar[6] = {};
    c.Clear();
    c.Divide(3, 2, 0.003, 0.035);  // 3.5% top/bottom margin reserves space for PageTitle
    for (int r = 0; r < kNRuns; ++r) {
        c.cd(r + 1);
        gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.20);
        gPad->SetRightMargin(0.05); gPad->SetTopMargin(0.10);
        gPad->SetTickx(1); gPad->SetTicky(1);
        if (nFidEvt[r] == 0) continue;

        hBar[r] = new TH1F(Form("hBar_%s", kRuns[r].label.Data()), "",
                            kNCap, 0, kNCap);
        hBar[r]->SetDirectory(nullptr);
        for (int ic = 0; ic < kNCap; ++ic) {
            hBar[r]->GetXaxis()->SetBinLabel(ic + 1, kChName[ic]);
            double meanHG = (nActive[r][ic] > 0)
                            ? sumHGPeak[r][ic] / nActive[r][ic] : 0.;
            hBar[r]->SetBinContent(ic + 1, meanHG);
        }
        hBar[r]->SetFillColor(kREnergyCols[r]);
        hBar[r]->SetLineColor(kREnergyCols[r]);
        hBar[r]->GetYaxis()->SetTitle("Mean HG peak (mV)");
        hBar[r]->GetXaxis()->SetLabelSize(0.085);
        hBar[r]->GetYaxis()->SetTitleSize(0.075);
        hBar[r]->GetYaxis()->SetTitleOffset(0.95);
        // Common y-axis from 0 across all 6 panels.  Per-panel auto-ranging
        // started the bars near their own minimum (e.g. 790-820 mV), turning a
        // ~4% channel spread into dramatic-looking bars and hiding the real
        // energy evolution.  A shared 0-900 mV axis shows both honestly.
        hBar[r]->GetYaxis()->SetRangeUser(0., 900.);
        hBar[r]->Draw("BAR");

        TLatex lab; lab.SetNDC(); lab.SetTextSize(0.08); lab.SetTextFont(72);
        lab.DrawLatex(0.60, 0.84, Form("%.0f GeV", kRuns[r].energy_GeV));
    }
    c.cd(0);
    PageTitle("Mean HG peak amplitude per capillary (timing fiducial events)");
    c.Print(outPDF);
    for (int r = 0; r < kNRuns; ++r) delete hBar[r];  // safe after Print

    // ── Page 3: Channel active fraction vs energy (8 lines) ──────────────────
    // Active = hg_peak > kHG_minPeak. Measures how often each capillary
    // produces a usable signal — drops near edges or in noisy channels.
    c.Clear(); c.cd(); StylePad();

    TGraph* gAct[8] = {};
    TLegend leg3(0.62, 0.15, 0.93, 0.53);
    leg3.SetBorderSize(0); leg3.SetFillStyle(0); leg3.SetTextSize(0.038);
    double yMinAct = 1.0;
    for (int ic = 0; ic < kNCap; ++ic) {
        gAct[ic] = new TGraph();
        for (int r = 0; r < kNRuns; ++r) {
            if (nFidEvt[r] == 0) continue;
            double frac = static_cast<double>(nActive[r][ic]) / nFidEvt[r];
            gAct[ic]->SetPoint(gAct[ic]->GetN(), kRuns[r].energy_GeV, 100. * frac);
            yMinAct = std::min(yMinAct, frac);
        }
        gAct[ic]->SetLineColor(kRChannelCols[ic]);
        gAct[ic]->SetMarkerColor(kRChannelCols[ic]);
        gAct[ic]->SetMarkerStyle(20 + ic);
        gAct[ic]->SetLineWidth(2);
    }
    TH1F* frame3 = static_cast<TH1F*>(
        c.DrawFrame(15., std::max(0., (yMinAct - 0.05) * 100.),
                    165., 108.,
                    ";Beam energy (GeV);Active fraction (%)  [hg_peak > 20 mV]"));
    frame3->GetXaxis()->SetTitleSize(0.050);
    frame3->GetYaxis()->SetTitleSize(0.050);
    frame3->GetYaxis()->SetTitleOffset(1.30);
    for (int ic = 0; ic < kNCap; ++ic) {
        if (gAct[ic]->GetN() > 0) {
            gAct[ic]->Draw("PL SAME");
            leg3.AddEntry(gAct[ic], kChName[ic], "lp");
        }
    }
    leg3.Draw();
    PageTitle("Per-channel hit efficiency vs beam energy (timing fiducial, MCP-quality)");
    c.Print(outPDF);

    // ── Page 4: Heat map — σ_t [ps] for all channels × energies ─────────────
    // Colour encodes σ_t: read as a table, each row is a capillary,
    // each column is an energy.  Bright = worse timing.
    c.Clear(); c.cd();
    gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.11);
    gPad->SetRightMargin(0.18); gPad->SetTopMargin(0.10);
    gPad->SetTickx(0); gPad->SetTicky(0);
    gStyle->SetPalette(kCherry); TColor::InvertPalette();

    TH2F* hHeat = new TH2F("hHeat", "", kNRuns, 0, kNRuns, kNCap, 0, kNCap);
    hHeat->SetDirectory(nullptr);
    for (int ie = 1; ie <= kNRuns; ++ie)
        hHeat->GetXaxis()->SetBinLabel(ie, Form("%.0fG", kRuns[ie-1].energy_GeV));
    for (int ic = 1; ic <= kNCap; ++ic)
        hHeat->GetYaxis()->SetBinLabel(ic, kChName[ic-1]);
    for (int ic = 0; ic < kNCap; ++ic)
        for (int ie = 0; ie < kNRuns; ++ie)
            if (heatSig[ic][ie] > 0.)
                hHeat->SetBinContent(ie + 1, ic + 1, heatSig[ic][ie]);

    hHeat->GetXaxis()->SetTitle("Beam energy");
    hHeat->GetZaxis()->SetTitle("#sigma_{t} (ps)");
    hHeat->GetXaxis()->SetLabelSize(0.055);
    hHeat->GetYaxis()->SetLabelSize(0.055);
    hHeat->GetZaxis()->SetTitleOffset(1.5);
    // SetPaintTextFormat is GLOBAL and applied at PAINT (Print) time, not at
    // Draw time — so the format must still be "%.0f" when c.Print() runs.
    // Restoring "%g" before Print made every cell paint as the literal "%g".
    gStyle->SetPaintTextFormat("%.0f");   // integer ps values in cells
    hHeat->Draw("COLZ TEXT");  // TEXT overlays the numeric value in each cell

    PageTitle("#sigma_{t} (ps) per capillary per energy  [CFD-5%, timing fiducial]");
    c.Print(outPDF + ")");
    gStyle->SetPaintTextFormat("%g");     // restore default AFTER printing

    // Cleanup
    for (int ic = 0; ic < kNCap; ++ic) delete gAct[ic];
    delete hHeat;
    fSum->Close();
    delete fSum;
    std::cout << "compareEnergies: Group 3 -> " << outPDF << "\n";
}

// ===========================================================================
// Group 4 — Timing Methods Comparison
//
//   Page 1: Combination method σ_t vs energy — which combination wins?
//   Page 2: CFD fraction scan (10/20/30/50%) — is 20% optimal?
//   Page 3: Walk correction benefit (M1 vs M5/M6/M7) — is correction worth it?
//   Page 4: σ_t vs 1/√E with stochastic fits — the headline physics plot
// ===========================================================================
static void PlotGroup4_TimingMethods()
{
    std::cout << "compareEnergies: Group 4 — Timing Methods Comparison\n";

    // Open summary ROOT files
    TString tSumPath     = Form("%stiming_summary.root",      kSumDir);
    TString tMethPath    = Form("%stiming_methods.root",      kSumDir);
    TString tBinsPath    = Form("%stiming_energy_bins.root",  kSumDir);

    TFile* fTSum  = TFile::Open(tSumPath);
    TFile* fTMeth = TFile::Open(tMethPath);
    TFile* fTBins = TFile::Open(tBinsPath);

    if (!fTSum  || fTSum->IsZombie())  {
        std::cerr << "compareEnergies: missing " << tSumPath  << "\n"; return; }
    if (!fTMeth || fTMeth->IsZombie()) {
        std::cerr << "compareEnergies: missing " << tMethPath << "\n";
        fTSum->Close(); delete fTSum; return; }
    if (!fTBins || fTBins->IsZombie()) {
        std::cerr << "compareEnergies: missing " << tBinsPath << "\n";
        fTSum->Close(); delete fTSum; fTMeth->Close(); delete fTMeth; return; }

    // Combination methods from timing_summary.root
    // M0=best single, M1=4-corner, M2=mean-all, M3=A²-wgt-all, M4=A²-wgt-4corner
    static const char* kCombLabel[5] = {
        "Best single channel",
        "4-corner (U+D)/2",
        "Mean all 8",
        "A^{2}-wgt all 8",
        "A^{2}-wgt 4 corners"
    };
    static const int kCombCol[5] = {
        kCyan+2, kGreen+2, kBlue+1, kRed+1, kOrange+1
    };
    TGraphErrors* gComb[5] = {};
    for (int m = 0; m < 5; ++m) {
        gComb[m] = static_cast<TGraphErrors*>(
            fTSum->Get(Form("gTiming_M%d", m)));
    }

    // CFD fraction methods from timing_methods.root:
    // m0=CFD10%, m1=CFD20%, m2=CFD30%, m3=CFD50%
    static const char* kCFDLabel[4] = {
        "CFD 10%", "CFD 20%", "CFD 30%", "CFD 50%"
    };
    static const int kCFDCol[4] = { kCyan+2, kBlack, kGreen+2, kOrange+1 };
    TGraph* gCFD[4] = {};
    for (int m = 0; m < 4; ++m)
        gCFD[m] = static_cast<TGraph*>(fTMeth->Get(Form("gBest_m%d", m)));

    // Walk correction methods: m1=CFD20% (baseline), m5=TOT, m6=1/A, m7=HG/LG
    static const char* kWalkLabel[4] = {
        "CFD 20% (baseline)",
        "LED + TOT walk corr.",
        "CFD + 1/A walk corr.",
        "CFD + HG/LG ratio corr."
    };
    static const int kWalkCol[4] = { kBlack, kMagenta+1, kBlue+1, kViolet+2 };
    int kWalkIdx[4] = { 1, 5, 6, 7 };
    TGraph* gWalk[4] = {};
    for (int i = 0; i < 4; ++i)
        gWalk[i] = static_cast<TGraph*>(fTMeth->Get(Form("gBest_m%d", kWalkIdx[i])));

    // Energy-binned best estimator from timing_energy_bins.root
    TGraphErrors* gTEB[3] = {};
    gTEB[0] = static_cast<TGraphErrors*>(fTBins->Get("gBestSigma_teb_m0"));
    gTEB[1] = static_cast<TGraphErrors*>(fTBins->Get("gBestSigma_teb_m1"));
    gTEB[2] = static_cast<TGraphErrors*>(fTBins->Get("gBestSigma_teb_m2"));
    TGraphErrors* gPaper = static_cast<TGraphErrors*>(fTBins->Get("gPaper_teb"));

    static const char* kTEBLabel[3] = {
        "(DW#minusUP)/2  CFD-20%",
        "(DW+UP)/2   CFD-20%",
        "(DW#minusUP)/2  M7-corr."
    };
    static const int kTEBCol[3] = { kBlue+1, kRed+1, kViolet+2 };

    TString outPDF = Form("%scross_energy_timing_methods.pdf", kSumDir);
    TCanvas c("c_tm", "", 960, 720);

    // ── Page 1: Combination methods σ_t vs energy ────────────────────────────
    c.Clear(); c.cd(); StylePad();
    double yMaxC = 0.;
    for (int m = 0; m < 5; ++m)
        if (gComb[m])
            for (int p = 0; p < gComb[m]->GetN(); ++p)
                yMaxC = std::max(yMaxC, gComb[m]->GetY()[p]);
    if (yMaxC < 10.) yMaxC = 300.;

    TH1F* fC = static_cast<TH1F*>(
        c.DrawFrame(15., 50., 165., 1.35 * yMaxC,
                    ";Beam energy (GeV);#sigma_{t} (ps)"));
    fC->GetXaxis()->SetTitleSize(0.050);
    fC->GetYaxis()->SetTitleSize(0.050);
    fC->GetYaxis()->SetTitleOffset(1.30);

    TLegend legC(0.65, 0.62, 0.93, 0.88);
    legC.SetBorderSize(0); legC.SetFillStyle(0); legC.SetTextSize(0.040);
    for (int m = 0; m < 5; ++m) {
        if (!gComb[m]) continue;
        gComb[m]->SetLineColor(kCombCol[m]);
        gComb[m]->SetMarkerColor(kCombCol[m]);
        gComb[m]->SetMarkerStyle(20 + m);
        gComb[m]->SetMarkerSize(1.2);
        gComb[m]->SetLineWidth(2);
        gComb[m]->Draw("PL SAME");
        legC.AddEntry(gComb[m], kCombLabel[m], "lp");
    }
    legC.Draw();
    PageTitle("Channel-combination methods: #sigma_{t} vs beam energy");
    c.Print(outPDF + "(");

    // ── Page 2: CFD fraction scan (10/20/30/50%) ─────────────────────────────
    c.Clear(); c.cd(); StylePad();
    double yMaxCFD = 0.;
    for (int m = 0; m < 4; ++m)
        if (gCFD[m])
            for (int p = 0; p < gCFD[m]->GetN(); ++p)
                yMaxCFD = std::max(yMaxCFD, gCFD[m]->GetY()[p]);
    if (yMaxCFD < 10.) yMaxCFD = 300.;

    TH1F* fCFD = static_cast<TH1F*>(
        c.DrawFrame(15., 50., 165., 1.35 * yMaxCFD,
                    ";Beam energy (GeV);Best-channel #sigma_{t} (ps)"));
    fCFD->GetXaxis()->SetTitleSize(0.050);
    fCFD->GetYaxis()->SetTitleSize(0.050);
    fCFD->GetYaxis()->SetTitleOffset(1.30);

    TLegend legCFD(0.65, 0.64, 0.93, 0.88);
    legCFD.SetBorderSize(0); legCFD.SetFillStyle(0); legCFD.SetTextSize(0.042);
    for (int m = 0; m < 4; ++m) {
        if (!gCFD[m]) continue;
        gCFD[m]->SetLineColor(kCFDCol[m]);
        gCFD[m]->SetMarkerColor(kCFDCol[m]);
        gCFD[m]->SetMarkerStyle(20 + m);
        gCFD[m]->SetMarkerSize(1.2);
        gCFD[m]->SetLineWidth(m == 1 ? 3 : 1);  // highlight CFD-20% (m=1)
        gCFD[m]->Draw("PL SAME");
        legCFD.AddEntry(gCFD[m], kCFDLabel[m], "lp");
    }
    legCFD.Draw();
    PageTitle("CFD fraction comparison: best-channel #sigma_{t} vs beam energy");
    c.Print(outPDF);

    // ── Page 3: Walk correction benefit ──────────────────────────────────────
    c.Clear(); c.cd(); StylePad();
    double yMaxW = 0.;
    for (int i = 0; i < 4; ++i)
        if (gWalk[i])
            for (int p = 0; p < gWalk[i]->GetN(); ++p)
                yMaxW = std::max(yMaxW, gWalk[i]->GetY()[p]);
    if (yMaxW < 10.) yMaxW = 300.;

    TH1F* fW = static_cast<TH1F*>(
        c.DrawFrame(15., 50., 165., 1.35 * yMaxW,
                    ";Beam energy (GeV);Best-channel #sigma_{t} (ps)"));
    fW->GetXaxis()->SetTitleSize(0.050);
    fW->GetYaxis()->SetTitleSize(0.050);
    fW->GetYaxis()->SetTitleOffset(1.30);

    TLegend legW(0.65, 0.64, 0.93, 0.88);
    legW.SetBorderSize(0); legW.SetFillStyle(0); legW.SetTextSize(0.040);
    for (int i = 0; i < 4; ++i) {
        if (!gWalk[i]) continue;
        gWalk[i]->SetLineColor(kWalkCol[i]);
        gWalk[i]->SetMarkerColor(kWalkCol[i]);
        gWalk[i]->SetMarkerStyle(20 + i);
        gWalk[i]->SetMarkerSize(1.2);
        gWalk[i]->SetLineWidth(i == 0 ? 1 : 2);
        gWalk[i]->SetLineStyle(i == 0 ? 1 : 2);  // baseline solid, corrected dashed
        gWalk[i]->Draw("PL SAME");
        legW.AddEntry(gWalk[i], kWalkLabel[i], "lp");
    }
    legW.Draw();
    PageTitle("Walk correction benefit: best-channel #sigma_{t} vs beam energy");
    c.Print(outPDF);

    // ── Page 4: σ_t vs 1/√E — the headline physics plot ─────────────────────
    // Converts to 1/√E axis so stochastic scaling appears as a straight line.
    // Fit σ_t = √(a²/E + b²) to extract stochastic (a) and constant (b) terms.
    c.Clear(); c.cd(); StylePad();

    // Build 1/√E graphs from the energy-binned best-estimator data
    TGraph* gTEB_sqrtE[3] = {};
    TGraph* gPaper_sqrtE  = new TGraph();
    for (int m = 0; m < 3; ++m) {
        gTEB_sqrtE[m] = new TGraph();
        if (!gTEB[m]) continue;
        for (int p = 0; p < gTEB[m]->GetN(); ++p) {
            double E = gTEB[m]->GetX()[p];
            if (E > 0.)
                gTEB_sqrtE[m]->SetPoint(gTEB_sqrtE[m]->GetN(),
                                        1.0 / std::sqrt(E),
                                        gTEB[m]->GetY()[p]);
        }
    }
    if (gPaper) {
        for (int p = 0; p < gPaper->GetN(); ++p) {
            double E = gPaper->GetX()[p];
            if (E > 0.)
                gPaper_sqrtE->SetPoint(gPaper_sqrtE->GetN(),
                                       1.0 / std::sqrt(E),
                                       gPaper->GetY()[p]);
        }
    }

    double yMaxB = 0.;
    for (int m = 0; m < 3; ++m)
        if (gTEB[m])
            for (int p = 0; p < gTEB[m]->GetN(); ++p)
                yMaxB = std::max(yMaxB, gTEB[m]->GetY()[p]);
    for (int p = 0; p < gPaper_sqrtE->GetN(); ++p)
        yMaxB = std::max(yMaxB, gPaper_sqrtE->GetY()[p]);
    if (yMaxB < 10.) yMaxB = 100.;

    // x range: 1/√150 ≈ 0.082 to 1/√25 = 0.200
    TH1F* fB = static_cast<TH1F*>(
        c.DrawFrame(0.075, 0., 0.215, 1.40 * yMaxB,
                    ";1/#sqrt{E_{beam}}  (GeV^{-1/2});#sigma_{t} (ps)"));
    fB->GetXaxis()->SetTitleSize(0.050);
    fB->GetYaxis()->SetTitleSize(0.050);
    fB->GetYaxis()->SetTitleOffset(1.30);

    // Stochastic fit on 1/√E axis.
    // x = 1/√E  →  x² = 1/E  →  σ_t = √(a²/E + b²) = √(a²·x² + b²)
    // TF1 formula: sqrt([0]*[0]*x*x + [1]*[1])   ([0]=a, [1]=b)
    // Seed: a≈200 ps√GeV, b≈30 ps  (from rough fit to the data)
    TF1* fFit[3] = {};
    double fitA[3] = {}, fitB[3] = {};
    for (int m = 0; m < 3; ++m) {
        fFit[m] = new TF1(Form("fFit_ce_m%d", m),
                          "sqrt([0]*[0]*x*x + [1]*[1])", 0.075, 0.215);
        fFit[m]->SetParameters(200., 30.);
        fFit[m]->SetLineColor(kTEBCol[m]);
        fFit[m]->SetLineStyle(2);
        fFit[m]->SetLineWidth(2);
        if (gTEB_sqrtE[m] && gTEB_sqrtE[m]->GetN() >= 3) {
            gTEB_sqrtE[m]->Fit(fFit[m], "RQ");
            fitA[m] = std::fabs(fFit[m]->GetParameter(0));
            fitB[m] = std::fabs(fFit[m]->GetParameter(1));
            fFit[m]->DrawCopy("SAME");
        }
    }

    // Paper reference fit (same axis convention: x = 1/√E)
    TF1* fPaper = new TF1("fPaper_ce",
                          "sqrt(256.*256.*x*x + 17.5*17.5)", 0.075, 0.215);
    fPaper->SetLineColor(kGray+2); fPaper->SetLineStyle(3); fPaper->SetLineWidth(2);
    fPaper->DrawCopy("SAME");

    // Data points on top of fit lines
    for (int m = 0; m < 3; ++m) {
        if (!gTEB_sqrtE[m] || gTEB_sqrtE[m]->GetN() == 0) continue;
        gTEB_sqrtE[m]->SetLineColor(kTEBCol[m]);
        gTEB_sqrtE[m]->SetMarkerColor(kTEBCol[m]);
        gTEB_sqrtE[m]->SetMarkerStyle(20 + m);
        gTEB_sqrtE[m]->SetMarkerSize(1.3);
        gTEB_sqrtE[m]->SetLineWidth(2);
        gTEB_sqrtE[m]->Draw("PL SAME");
    }
    if (gPaper_sqrtE->GetN() > 0) {
        gPaper_sqrtE->SetLineColor(kGray+2);
        gPaper_sqrtE->SetMarkerColor(kGray+2);
        gPaper_sqrtE->SetMarkerStyle(25);
        gPaper_sqrtE->SetMarkerSize(1.2);
        gPaper_sqrtE->Draw("P SAME");
    }

    // CMS Phase II Barrel Timing Layer requirement: 30 ps at MIP equivalent
    // At test beam energies, the BTL works with MIPs (~200 MeV deposit),
    // so this is a constant requirement line, not energy-dependent.
    // For comparison, we plot it as a horizontal dashed line.
    TLine* lBTL = new TLine(0.075, 30., 0.215, 30.);
    lBTL->SetLineColor(kGray+1); lBTL->SetLineStyle(2); lBTL->SetLineWidth(3);
    lBTL->Draw("SAME");
    TLatex btlLab; btlLab.SetNDC(); btlLab.SetTextSize(0.034);
    btlLab.SetTextColor(kGray+1);
    btlLab.DrawLatex(0.20, 0.28, "CMS BTL Phase II requirement: 30 ps");

    // Gap annotation between our best result and published 27 ps
    if (gTEB[0] && gTEB[0]->GetN() > 0) {
        // Find our 150 GeV sigma_t: last point in gTEB_sqrtE[0] (highest x = lowest E),
        // but 150 GeV has x = 1/sqrt(150) ~ 0.0816 -- scan for it
        double ourVal = -1.;
        for (int p = 0; p < gTEB_sqrtE[0]->GetN(); ++p) {
            double xp = gTEB_sqrtE[0]->GetX()[p];
            if (std::fabs(xp - 1.0/std::sqrt(150.)) < 0.005) {
                ourVal = gTEB_sqrtE[0]->GetY()[p];
                break;
            }
        }
        // Fallback: use the point with smallest x (highest energy = 150 GeV)
        if (ourVal < 0. && gTEB_sqrtE[0]->GetN() > 0) {
            int minIdx = 0;
            for (int p = 1; p < gTEB_sqrtE[0]->GetN(); ++p)
                if (gTEB_sqrtE[0]->GetX()[p] < gTEB_sqrtE[0]->GetX()[minIdx])
                    minIdx = p;
            ourVal = gTEB_sqrtE[0]->GetY()[minIdx];
        }
        if (ourVal > 27.) {
            TArrow* arr = new TArrow(0.0816, ourVal, 0.0816, 27., 0.015, "|>");
            arr->SetLineColor(kBlue+1); arr->SetFillColor(kBlue+1);
            arr->SetLineWidth(2);
            arr->Draw();
            TLatex gapLab; gapLab.SetTextSize(0.030); gapLab.SetTextColor(kBlue+1);
            gapLab.DrawLatex(0.0826, 0.5*(ourVal + 27.),
                Form("gap: %.0f ps", ourVal - 27.));
        }
    }

    // Legend in the EMPTY upper-left region (the data band rises left->right,
    // so top-left is clear).  Short entries to avoid right-edge truncation.
    TLegend legB(0.16, 0.70, 0.55, 0.90);
    legB.SetBorderSize(0); legB.SetFillStyle(0); legB.SetTextSize(0.032);
    for (int m = 0; m < 3; ++m) {
        if (!gTEB_sqrtE[m] || gTEB_sqrtE[m]->GetN() == 0) continue;
        legB.AddEntry(gTEB_sqrtE[m], kTEBLabel[m], "lp");
    }
    legB.AddEntry(gPaper_sqrtE, "arXiv:2401.01747 (LYSO/SiPM)", "p");
    legB.AddEntry(lBTL, "CMS BTL req. 30 ps", "l");
    legB.Draw();

    // Fit-parameter annotations: compact block just below the legend (left side,
    // clear of the data band and the legend).  Colour-matched per method.
    {
        TLatex fann; fann.SetNDC(); fann.SetTextSize(0.029);
        for (int m = 0; m < 3; ++m) {
            if (fitA[m] <= 0. && fitB[m] <= 0.) continue;
            fann.SetTextColor(kTEBCol[m]);
            fann.DrawLatex(0.17, 0.66 - m * 0.050,
                Form("a = %.0f ps#sqrt{GeV},  b = %.1f ps", fitA[m], fitB[m]));
        }
        fann.SetTextColor(kGray+2);
        fann.DrawLatex(0.17, 0.66 - 3 * 0.050,
            "arXiv: a = 256 ps#sqrt{GeV}, b = 17.5 ps");
    }

    // Beam-energy labels near the bottom axis.  The data band sits at >=37 ps,
    // so labels at y ~ 0.07*frameTop are well clear of both the data and the
    // (now top-left) legend that the top placement collided with.
    if (gTEB[0] && gTEB[0]->GetN() > 0) {
        const double yLab = 0.06 * (1.40 * yMaxB);
        TLatex elab; elab.SetTextSize(0.026); elab.SetTextAlign(21);
        elab.SetTextColor(kGray+1);
        for (int r = 0; r < kNRuns; ++r) {
            double E = kRuns[r].energy_GeV;
            double x = 1.0 / std::sqrt(E);
            elab.DrawLatex(x, yLab, Form("%.0f", E));
        }
        elab.SetTextAlign(11);
        elab.DrawLatex(0.203, yLab, "GeV");
    }
    PageTitle("Energy-binned #sigma_{t} vs 1/#sqrt{E_{beam}}  (stochastic scaling)");
    c.Print(outPDF + ")");

    // Cleanup
    for (int m = 0; m < 3; ++m) { delete gTEB_sqrtE[m]; delete fFit[m]; }
    delete gPaper_sqrtE; delete fPaper;
    fTSum->Close();  delete fTSum;
    fTMeth->Close(); delete fTMeth;
    fTBins->Close(); delete fTBins;
    std::cout << "compareEnergies: Group 4 -> " << outPDF << "\n";
}

// ===========================================================================
// Entry point
// ===========================================================================
void compareEnergies()
{
    TH1::AddDirectory(kFALSE);   // prevent TFile from owning histograms

    ApplyRADiCALStyle();   // RADiCALStyle.h via PlotUtils.h

    gSystem->mkdir(kSumDir, kTRUE);

    PlotGroup1_BeamQuality();
    PlotGroup2_Containment();
    PlotGroup3_Channels();
    PlotGroup4_TimingMethods();

    std::cout << "\ncompareEnergies: all done.  Output in "
              << kSumDir << "\n";
}

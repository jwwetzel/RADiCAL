// ============================================================================
// positionCorrection.C — charge-sharing position reconstruction + position-
//                        binned timing correction  (Gap-closing item #1)
// ============================================================================
//
// MOTIVATION (Ledovskoy, Dec 2021 + Jan-2021 simulation)
// ------------------------------------------------------
// The combination timing has a position-dependent walk: where the shower lands
// in the cell changes the relative light arrival and biases the measured time.
// Ledovskoy reconstructed the impact position from the calorimeter's own charge
// sharing (asymmetry of the corner amplitudes) and removed the walk by working
// at fixed position — 154 ps -> 113 ps.  His simulation predicted this
// "position correction" is worth ~25 ps.  We MEASURE the non-uniformity
// (uniformityScan) but never CORRECT for it.  This macro does.
//
// METHOD (no assumptions — everything is validated)
// -------------------------------------------------
//   1. Build charge-sharing asymmetries from the 4 corner LG (energy) amplitudes:
//        asym_x = (East - West)/Sum   East = NE+SE,  West = NW+SW
//        asym_y = (North - South)/Sum  North = NW+NE, South = SW+SE
//   2. VALIDATE that asym tracks the wire-chamber position (Pearson rho).
//      If charge sharing does not track position for our geometry, the
//      technique does not apply and we say so.
//   3. Train a 2D position correction (mean A^2-combo timing per position bin)
//      on EVEN events; apply to ODD events; measure the core sigma_t before and
//      after OUT-OF-SAMPLE.  Done both with the calorimeter asymmetry and with
//      the (independent) wire-chamber position, for cross-check.
//
// Output:
//   Analysis/Output/Summary/position_correction.pdf  (4 pages)
//   Analysis/Output/Summary/position_correction.root (sigma vs energy graphs)
//
// Usage:  root -l -b -q 'Analysis/positionCorrection.C+'
// ============================================================================

#include "ChannelConfig.h"
#include "PlotUtils.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile2D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TSystem.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

static const float kNoTime = -1e9f;

// ---------------------------------------------------------------------------
// A^2-weighted combination of the 8 HG CFD-20% channels (MCP-referenced).
// Returns kNoTime if fewer than 2 valid channels.  (Same estimator as M3.)
// ---------------------------------------------------------------------------
static float A2Combo(const Float_t hg_peak[8], const Float_t hg_cfd[8])
{
    double sw = 0., swt = 0.; int n = 0;
    for (int i = 0; i < kNCap; ++i) {
        if (hg_peak[i] < kHG_minPeak) continue;
        if (hg_cfd[i]  < -1e5f)       continue;
        const double w = static_cast<double>(hg_peak[i]) * hg_peak[i];
        sw += w; swt += w * hg_cfd[i]; ++n;
    }
    return (n >= 2 && sw > 0.) ? static_cast<float>(swt / sw) : kNoTime;
}

// ---------------------------------------------------------------------------
// Corr2D — reproducible 2D-binned offset correction (the position analogue of
// drs4::StopCellCorrection).  Train on one event set, apply to another.
// ---------------------------------------------------------------------------
class Corr2D {
public:
    Corr2D(int nb, double lo, double hi)
        : nb_(nb), lo_(lo), hi_(hi), ready_(false), gmean_(0.),
          sum_(nb*nb, 0.), cnt_(nb*nb, 0), off_(nb*nb, 0.) {}

    void Accumulate(double ax, double ay, double v) {
        const int b = Bin(ax, ay);
        if (b >= 0) { sum_[b] += v; cnt_[b] += 1; }
    }
    void Finalize(long minCount = 20) {
        double gs = 0.; long gc = 0;
        for (size_t b = 0; b < sum_.size(); ++b) { gs += sum_[b]; gc += cnt_[b]; }
        gmean_ = (gc > 0) ? gs / gc : 0.;
        for (size_t b = 0; b < sum_.size(); ++b)
            off_[b] = (cnt_[b] >= minCount) ? sum_[b]/cnt_[b] - gmean_ : 0.;
        ready_ = true;
    }
    double Offset(double ax, double ay) const {
        if (!ready_) return 0.;
        const int b = Bin(ax, ay);
        return (b < 0) ? 0. : off_[b];
    }
    double OffsetRMS() const {
        double s = 0., s2 = 0.; int n = 0;
        for (size_t b = 0; b < off_.size(); ++b)
            if (cnt_[b] > 0) { s += off_[b]; s2 += off_[b]*off_[b]; ++n; }
        if (n < 2) return 0.;
        const double m = s/n; const double v = s2/n - m*m;
        return v > 0. ? std::sqrt(v) : 0.;
    }
private:
    int Bin(double ax, double ay) const {
        if (ax < lo_ || ax >= hi_ || ay < lo_ || ay >= hi_) return -1;
        int ix = static_cast<int>((ax - lo_) / (hi_ - lo_) * nb_);
        int iy = static_cast<int>((ay - lo_) / (hi_ - lo_) * nb_);
        if (ix < 0 || ix >= nb_ || iy < 0 || iy >= nb_) return -1;
        return ix * nb_ + iy;
    }
    int nb_; double lo_, hi_; bool ready_; double gmean_;
    std::vector<double> sum_; std::vector<long> cnt_; std::vector<double> off_;
};

// ---------------------------------------------------------------------------
// Core Gaussian sigma [ps] of a value list, fit within +-win ns of the mean.
// ---------------------------------------------------------------------------
static double CoreSigmaPs(const std::vector<float>& v, double win)
{
    if (static_cast<int>(v.size()) < 100) return -1.;
    double m = 0.; for (float x : v) m += x; m /= v.size();
    TH1F h("_pc_core", "", 200, m - win, m + win);
    h.SetDirectory(nullptr);
    for (float x : v) h.Fill(x);
    double mu, muE, s, sE;
    FitGaussCore(&h, 2.0, mu, muE, s, sE);
    return (s > 0.) ? s * 1000. : -1.;
}

// ---------------------------------------------------------------------------
// One per-event record (kept so we can train/apply the split-half correction).
// ---------------------------------------------------------------------------
struct PCEvent { float ax, ay, xt, yt, tcombo, tcorner; };

// ---------------------------------------------------------------------------
// 4-corner (DW-UP)/2 mean — the MCP-FREE headline-type estimator.
// For each corner k, (t_HG,D - t_HG,U)/2 cancels the common MCP reference.
// Returns kNoTime if no corner pair is valid.
// ---------------------------------------------------------------------------
static float CornerDiff(const Float_t hg_peak[8], const Float_t hg_cfd[8])
{
    double sum = 0.; int n = 0;
    for (int k = 0; k < 4; ++k) {           // NW,NE,SE,SW ; U = k+4
        if (hg_peak[k] < kHG_minPeak || hg_peak[k+4] < kHG_minPeak) continue;
        if (hg_cfd[k]  < -1e5f || hg_cfd[k+4] < -1e5f)              continue;
        sum += 0.5 * (hg_cfd[k] - hg_cfd[k+4]); ++n;
    }
    return (n > 0) ? static_cast<float>(sum / n) : kNoTime;
}

// ===========================================================================
// Main
// ===========================================================================
void positionCorrection()
{
    TH1::AddDirectory(kFALSE);
    ApplyRADiCALStyle();
    gSystem->mkdir("Analysis/Output/Summary", kTRUE);

    // Per-energy out-of-sample results
    std::vector<double> vE, vSigIncl, vSigAsym, vSigWC;

    // Detailed objects for the highest-energy page
    TProfile2D* pTimeAsym = nullptr;   // combo timing vs (asym_x, asym_y)
    TH2F*       hAsymVsXtrk = nullptr;  // validation: asym_x vs x_trk
    TH2F*       hAsymVsYtrk = nullptr;
    double rhoX = 0., rhoY = 0.;
    std::vector<float> detInclOdd, detAsymOdd;  // 150 GeV before/after (odd set)

    for (int iRun = 0; iRun < kNRuns; ++iRun) {
        const RunCfg& rc = kRuns[iRun];
        TString path = TString("Analysis/Output/") + rc.label + "/ntuple.root";
        TFile* f = TFile::Open(path);
        if (!f || f->IsZombie()) { std::cout << "  skip " << rc.label << "\n"; continue; }
        TTree* t = (TTree*)f->Get("rad");
        if (!t || t->GetEntries() == 0) { f->Close(); continue; }

        double xc, yc, tcfd, trms;
        ScanRunCenters(t, xc, yc, tcfd, trms);

        Bool_t  wc_ok;
        Float_t x_trk, y_trk, mcp_peak, sum_lg, sum_pb;
        Float_t hg_peak[8], hg_cfd[8], lg_peak[8];
        t->SetBranchAddress("wc_ok",   &wc_ok);
        t->SetBranchAddress("x_trk",   &x_trk);
        t->SetBranchAddress("y_trk",   &y_trk);
        t->SetBranchAddress("mcp_peak",&mcp_peak);
        t->SetBranchAddress("sum_lg",  &sum_lg);
        t->SetBranchAddress("sum_pb",  &sum_pb);
        t->SetBranchAddress("hg_peak",  hg_peak);
        t->SetBranchAddress("hg_cfd",   hg_cfd);
        t->SetBranchAddress("lg_peak",  lg_peak);

        std::vector<PCEvent> ev;
        ev.reserve(20000);

        const bool isTop = (iRun == kNRuns - 1);
        if (isTop) {
            hAsymVsXtrk = new TH2F("hAsymX",";x_{trk} (mm);asym_{x}",
                                   60, xc-6, xc+6, 60, -0.6, 0.6);
            hAsymVsYtrk = new TH2F("hAsymY",";y_{trk} (mm);asym_{y}",
                                   60, yc-6, yc+6, 60, -0.6, 0.6);
            pTimeAsym = new TProfile2D("pTimeAsym",
                ";asym_{x};asym_{y};#LT t_{combo} #GT (ns)",
                12,-0.6,0.6, 12,-0.6,0.6);
            hAsymVsXtrk->SetDirectory(nullptr);
            hAsymVsYtrk->SetDirectory(nullptr);
            pTimeAsym->SetDirectory(nullptr);
        }

        Long64_t nEv = t->GetEntries();
        for (Long64_t i = 0; i < nEv; ++i) {
            t->GetEntry(i);
            if (!wc_ok || mcp_peak < kMCP1_minPeak || mcp_peak > kMCP1_maxPeak) continue;
            float dx = x_trk - (float)xc, dy = y_trk - (float)yc;
            if (std::sqrt(dx*dx + dy*dy) >= (float)kFiducial_r_timing) continue;
            if (sum_lg > kSumLG_centroid && sum_pb >= kPb_maxRatio * sum_lg) continue;

            // Charge-sharing position from the 4 corner LG amplitudes
            // kCap: 0 NW-D,1 NE-D,2 SE-D,3 SW-D,4 NW-U,5 NE-U,6 SE-U,7 SW-U
            double A_NW = lg_peak[0] + lg_peak[4];
            double A_NE = lg_peak[1] + lg_peak[5];
            double A_SE = lg_peak[2] + lg_peak[6];
            double A_SW = lg_peak[3] + lg_peak[7];
            double sum  = A_NW + A_NE + A_SE + A_SW;
            if (sum <= 0.) continue;
            float ax = (float)(((A_NE + A_SE) - (A_NW + A_SW)) / sum);  // East-West
            float ay = (float)(((A_NW + A_NE) - (A_SW + A_SE)) / sum);  // North-South

            float tc = A2Combo(hg_peak, hg_cfd);
            if (tc <= kNoTime + 1e5f) continue;
            // keep the physics peak only
            if (std::fabs(tc - tcfd) > 5.f) continue;

            float tk = CornerDiff(hg_peak, hg_cfd);   // MCP-free corner estimator
            if (tk <= kNoTime + 1e5f) tk = kNoTime;    // may be invalid for some events

            ev.push_back({ax, ay, x_trk, y_trk, tc, tk});

            if (isTop) {
                hAsymVsXtrk->Fill(x_trk, ax);
                hAsymVsYtrk->Fill(y_trk, ay);
                pTimeAsym->Fill(ax, ay, tc);
            }
        }
        t->ResetBranchAddresses();
        f->Close();

        if (ev.size() < 400) { std::cout << "  " << rc.label << ": too few events\n"; continue; }

        // Train on EVEN events (both asym-based and WC-based corrections).
        // WC corrector works in centred coords (x_trk - xc, y_trk - yc), so its
        // window is symmetric about 0 (+-6 mm, 12 bins -> 1 mm/bin; fiducial is
        // +-3 mm so it is comfortably covered).
        Corr2D corrAsym(12, -0.6, 0.6);
        Corr2D corrWC  (12, -6.0, 6.0);
        Corr2D corrKnr (12, -0.6, 0.6);   // asym correction for the corner estimator
        for (size_t i = 0; i < ev.size(); i += 2) {
            corrAsym.Accumulate(ev[i].ax, ev[i].ay, ev[i].tcombo);
            corrWC  .Accumulate(ev[i].xt - xc, ev[i].yt - yc, ev[i].tcombo);
            if (ev[i].tcorner > kNoTime + 1e5f)
                corrKnr.Accumulate(ev[i].ax, ev[i].ay, ev[i].tcorner);
        }
        corrAsym.Finalize(20);
        corrWC.Finalize(20);
        corrKnr.Finalize(20);

        // Evaluate on ODD events
        std::vector<float> incl, cAsym, cWC, kIncl, kAsym;
        for (size_t i = 1; i < ev.size(); i += 2) {
            incl .push_back(ev[i].tcombo);
            cAsym.push_back(ev[i].tcombo - (float)corrAsym.Offset(ev[i].ax, ev[i].ay));
            cWC  .push_back(ev[i].tcombo - (float)corrWC.Offset(ev[i].xt - xc, ev[i].yt - yc));
            if (ev[i].tcorner > kNoTime + 1e5f) {
                kIncl.push_back(ev[i].tcorner);
                kAsym.push_back(ev[i].tcorner - (float)corrKnr.Offset(ev[i].ax, ev[i].ay));
            }
        }
        double sIncl = CoreSigmaPs(incl, 2.0);
        double sAsym = CoreSigmaPs(cAsym, 2.0);
        double sWC   = CoreSigmaPs(cWC,   2.0);
        // MCP-free corner estimator: inclusive vs asym-corrected (headline test)
        double sKincl = CoreSigmaPs(kIncl, 2.0);
        double sKasym = CoreSigmaPs(kAsym, 2.0);
        std::cout << Form("           (DW-UP)/2 MCP-free: incl %.1f -> asym %.1f ps\n",
                          sKincl, sKasym);

        if (isTop) {
            // recompute correlation factors
            rhoX = hAsymVsXtrk->GetCorrelationFactor();
            rhoY = hAsymVsYtrk->GetCorrelationFactor();
            detInclOdd = incl;
            detAsymOdd = cAsym;
        }

        vE.push_back(rc.energy_GeV);
        vSigIncl.push_back(sIncl);
        vSigAsym.push_back(sAsym);
        vSigWC.push_back(sWC);

        std::cout << Form("  %-7s  asym pattern RMS %.1f ps | "
                          "OOS sigma: incl %.1f -> asym %.1f -> WC %.1f ps\n",
                          rc.label.Data(), corrAsym.OffsetRMS()*1000.,
                          sIncl, sAsym, sWC);
    }

    if (vE.empty()) { std::cout << "positionCorrection: no data.\n"; return; }

    // =========================================================================
    // Output PDF
    // =========================================================================
    TString outPDF  = "Analysis/Output/Summary/position_correction.pdf";
    TString outROOT = "Analysis/Output/Summary/position_correction.root";
    TCanvas c("c_pc", "", 1100, 760);

    // -- Page 1: validation -- asymmetry vs wire-chamber position --------------
    c.Clear(); c.Divide(2, 1, 0.02, 0.02);
    if (hAsymVsXtrk) {
        c.cd(1); StylePad(true);
        hAsymVsXtrk->Draw("COLZ");
        DrawPadTitle("asym_{x} vs x_{trk}", 0.060);
        TLatex a; a.SetNDC(); a.SetTextSize(0.045); a.SetTextColor(kRRed);
        a.DrawLatex(0.18, 0.85, Form("#rho = %.2f", rhoX));
    }
    if (hAsymVsYtrk) {
        c.cd(2); StylePad(true);
        hAsymVsYtrk->Draw("COLZ");
        DrawPadTitle("asym_{y} vs y_{trk}", 0.060);
        TLatex a; a.SetNDC(); a.SetTextSize(0.045); a.SetTextColor(kRRed);
        a.DrawLatex(0.18, 0.85, Form("#rho = %.2f", rhoY));
    }
    c.cd(0);
    DrawPageTitle("Validation: calorimeter charge-sharing asymmetry vs wire-chamber position");
    c.Print(outPDF + "(");

    // -- Page 2: combo timing across the asymmetry plane (the walk) ------------
    c.Clear(); c.cd(); StylePad(true);
    if (pTimeAsym) {
        pTimeAsym->Draw("COLZ");
        DrawPadTitle("#LT A^{2}-combo timing #GT vs charge-sharing position (150 GeV)", 0.050);
    }
    DrawPageTitle("Position-dependent timing walk across the cell");
    c.Print(outPDF);

    // -- Page 3: before/after timing (150 GeV, odd/held-out set) ---------------
    c.Clear(); c.cd(); StylePad();
    if (!detInclOdd.empty()) {
        double m = 0.; for (float x : detInclOdd) m += x; m /= detInclOdd.size();
        TH1F* hI = new TH1F("hInclOdd",";A^{2}-combo  t_{HG}#minust_{MCP} (ns);events",120,m-1.2,m+1.2);
        TH1F* hA = new TH1F("hAsymOdd","",120,m-1.2,m+1.2);
        hI->SetDirectory(nullptr); hA->SetDirectory(nullptr);
        for (float x : detInclOdd) hI->Fill(x);
        for (float x : detAsymOdd) hA->Fill(x);
        hI->SetLineColor(kGray+2); hI->SetLineWidth(2);
        hA->SetLineColor(kRGreen);  hA->SetLineWidth(2);
        double yM = std::max(hI->GetMaximum(), hA->GetMaximum());
        hI->GetYaxis()->SetRangeUser(0., 1.25*yM);
        hI->Draw("HIST"); hA->Draw("HIST SAME");
        double sI = CoreSigmaPs(detInclOdd, 2.0), sA = CoreSigmaPs(detAsymOdd, 2.0);
        TLegend leg(0.60,0.72,0.93,0.86); leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.038);
        leg.AddEntry(hI, Form("inclusive  #sigma=%.0f ps", sI), "l");
        leg.AddEntry(hA, Form("pos-corrected  #sigma=%.0f ps", sA), "l");
        leg.Draw();
        TLatex a; a.SetNDC(); a.SetTextSize(0.034); a.SetTextColor(kGray+2);
        a.DrawLatex(0.16, 0.30, "150 GeV, held-out (odd) events; asym correction trained on even events");
    }
    DrawPageTitle("A^{2}-combo timing: inclusive vs position-corrected (out-of-sample)");
    c.Print(outPDF);

    // -- Page 4: sigma_t vs energy, inclusive vs position-corrected ------------
    c.Clear(); c.cd(); StylePad(false, true);
    double yMax = 0.;
    for (double s : vSigIncl) yMax = std::max(yMax, s);
    TH1F* fr = (TH1F*)c.DrawFrame(15., 0., 165., 1.30*yMax,
                                  ";Beam energy (GeV);A^{2}-combo #sigma_{t} (ps)");
    fr->GetXaxis()->SetTitleSize(0.048); fr->GetYaxis()->SetTitleSize(0.048);
    fr->GetYaxis()->SetTitleOffset(1.30);
    auto mkG = [&](std::vector<double>& y, int col, int mk){
        TGraph* g = new TGraph();
        for (size_t i = 0; i < vE.size(); ++i) if (y[i] > 0.) g->SetPoint(g->GetN(), vE[i], y[i]);
        g->SetMarkerStyle(mk); g->SetMarkerSize(1.4); g->SetMarkerColor(col);
        g->SetLineColor(col); g->SetLineWidth(2); return g;
    };
    TGraph* gI = mkG(vSigIncl, kGray+2, 20);
    TGraph* gA = mkG(vSigAsym, kRGreen,  21);
    TGraph* gW = mkG(vSigWC,   kRBlue,   22);
    gI->Draw("PL SAME"); gA->Draw("PL SAME"); gW->Draw("PL SAME");
    TLegend* Lc = MakeLegend(3);
    Lc->AddEntry(gI, "inclusive", "lp");
    Lc->AddEntry(gA, "charge-sharing", "lp");
    Lc->AddEntry(gW, "wire-chamber", "lp");
    Lc->Draw();
    DrawPageTitle("A^{2}-combo #sigma_{t} vs energy: inclusive vs position-corrected (out-of-sample)");
    c.Print(outPDF + ")");

    // -- ROOT output -----------------------------------------------------------
    TFile* fo = new TFile(outROOT, "RECREATE");
    gI->Write("gSigInclusive");
    gA->Write("gSigAsymCorrected");
    gW->Write("gSigWCCorrected");
    if (pTimeAsym)   pTimeAsym->Write();
    if (hAsymVsXtrk) hAsymVsXtrk->Write();
    if (hAsymVsYtrk) hAsymVsYtrk->Write();
    fo->Close(); delete fo;

    std::cout << "positionCorrection: done -> " << outPDF << "\n";
}

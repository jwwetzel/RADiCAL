// ============================================================================
// averageWaveforms.C — average DRS4 waveforms for RADiCAL CERN May 2023
// ============================================================================
//
// Reads raw "pulse" TChain data for each energy and computes per-channel
// average waveforms (mean ± RMS band) using TProfile, aligned on the
// CFD-20% crossing time of each channel.
//
// Channels: 8 HG capillary channels (kCap[0..7]) + MCP1 reference.
//
// Per-energy output (3×3 canvas):
//   output/<label>/average_waveforms.pdf
//
// Summary output (4×2 canvas, all energies overlaid):
//   output/Summary/average_waveforms_summary.pdf
//
// Usage (from repository root):
//   root -l -b -q 'Analysis/averageWaveforms.C+'
// ============================================================================

#include "ChannelConfig.h"   // kRuns, kNRuns, kCap, kNCap, kMCP1, kMCP1_t
#include "WaveformUtils.h"   // ExtractPulse, Pulse, kNoTime
#include "PlotUtils.h"       // StylePad

#include "TChain.h"
#include "TFile.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLine.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TAxis.h"
#include "TColor.h"
#include "TROOT.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TObjArray.h"
#include "TObjString.h"

#include <cmath>
#include <iostream>

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

// Number of channels: 8 HG capillaries + 1 MCP1
static const int kNWF = 9;

// TProfile binning: 300 bins, −6 to +18 ns around CFD crossing (HG: fast pulse)
static const int    kNBins   = 300;
static const double kTMin    = -6.0;
static const double kTMax    = 18.0;

// LG pulses (full-length fibre → energy) are much slower/longer than HG, so
// they need a wider window to show the full rise, peak and fall.
static const int    kNBins_LG = 320;
static const double kTMin_LG  = -10.0;
static const double kTMax_LG  =  1000.0;

// Maximum events per energy (raw data is large; keep runtime reasonable)
static const Long64_t kNMax  = 5000LL;

// Energy colors from RADiCALStyle.h (kREnergyCols, set by ApplyRADiCALStyle)

static const char* kSumDir = "output/Summary/";

// ---------------------------------------------------------------------------
// Helper: channel display name (index 0-7 → capillary, 8 → MCP1)
// ---------------------------------------------------------------------------
static const char* ChanName(int i)
{
    if (i < 8) return kCap[i].name;
    return "MCP1";
}

// ---------------------------------------------------------------------------
// Helper: book a TProfile for waveform averaging
// ---------------------------------------------------------------------------
static TProfile* BookWF(const char* label, const char* name,
                        int nbins = kNBins, double tmin = kTMin, double tmax = kTMax)
{
    TProfile* p = new TProfile(
        Form("wf_%s_%s", label, name),
        Form(";t #minus t_{CFD} (ns);Amplitude (mV)"),
        nbins, tmin, tmax);
    p->SetErrorOption("S");   // store RMS (std-dev) in error bars
    p->SetDirectory(nullptr); // keep in memory; do not attach to a TFile
    return p;
}

// ---------------------------------------------------------------------------
// Helper: draw one waveform sub-pad
//   pad     — the TPad to draw into (already cd()'d by caller)
//   hWF     — filled TProfile
//   title   — e.g. "150 GeV  NW-D"
//   iEnergy — energy index for colour (−1 → draw in black)
// ---------------------------------------------------------------------------
static void DrawWFPad(TProfile* hWF, const char* title,
                      int iEnergy = -1)
{
    StylePad();

    // Choose colour
    int col = (iEnergy >= 0 && iEnergy < kNRuns) ? kREnergyCols[iEnergy] : kBlack;

    // RMS band (E3 = filled area between mean ± error = RMS)
    TProfile* hBand = (TProfile*)hWF->Clone(Form("%s_band", hWF->GetName()));
    hBand->SetFillColorAlpha(col, 0.25);
    hBand->SetFillStyle(1001);
    hBand->SetLineColor(0);
    hBand->SetMarkerSize(0);
    hBand->Draw("E3");

    // Mean line on top
    TProfile* hLine = (TProfile*)hWF->Clone(Form("%s_line", hWF->GetName()));
    hLine->SetLineColor(col);
    hLine->SetLineWidth(2);
    hLine->SetMarkerSize(0);
    hLine->Draw("HIST SAME");

    // CFD crossing reference line at t_rel = 0
    double ymin = hWF->GetMinimum();
    double ymax = hWF->GetMaximum();
    // Guard against empty profiles
    if (ymin >= ymax) { ymin = 0.; ymax = 1.; }
    // Add a little headroom above max
    double ylo = (ymin < 0.) ? ymin * 1.1 : ymin * 0.9;
    double yhi = ymax * 1.2;
    hBand->GetYaxis()->SetRangeUser(ylo, yhi);

    TLine* lCFD = new TLine(0., ylo, 0., yhi);
    lCFD->SetLineStyle(2);
    lCFD->SetLineColor(kRed);
    lCFD->SetLineWidth(1);
    lCFD->Draw("SAME");

    // Annotate mean peak amplitude
    double meanPeak = 0.;
    int    nBins    = hWF->GetNbinsX();
    double binMax   = -1e9;
    for (int b = 1; b <= nBins; ++b) {
        double v = hWF->GetBinContent(b);
        if (v > binMax) { binMax = v; }
    }
    meanPeak = binMax;

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.050);
    lat.SetTextColor(kGray+2);
    if (meanPeak > 0.)
        lat.DrawLatex(0.62, 0.72, Form("%.0f mV", meanPeak));

    DrawPadTitle(title);   // from RADiCALStyle.h — margin-aware centering

    gPad->Modified();
    gPad->Update();
}

// ---------------------------------------------------------------------------
// ProcessEnergy — fill TProfile[kNWF] from raw "pulse" tree for one energy
//
// Returns true if at least one file was found and events were processed.
// hWF[0..7] → HG capillaries, hWF[8] → MCP1
// ---------------------------------------------------------------------------
static bool ProcessEnergy(const RunCfg& run, TProfile* hWF[], TProfile* hLG[])
{
    // Build TChain: split run.inFiles on ';'
    TChain* chain = new TChain("pulse");
    TString patterns(run.inFiles);
    TObjArray* parts = patterns.Tokenize(";");
    int nAdded = 0;
    for (int i = 0; i < parts->GetEntries(); ++i) {
        TString pat = ((TObjString*)parts->At(i))->GetString().Strip();
        // Check that the path is not empty before adding
        if (pat.IsNull()) continue;
        int n = chain->Add(pat);
        nAdded += n;
    }
    delete parts;

    if (nAdded == 0) {
        std::cerr << "[averageWaveforms] No files matched for "
                  << run.label << ": " << run.inFiles << "\n";
        delete chain;
        return false;
    }

    Long64_t nTotal = chain->GetEntries();
    Long64_t nProc  = std::min(nTotal, kNMax);

    std::cout << "[averageWaveforms] " << run.label
              << ": processing " << nProc
              << " / " << nTotal << " events\n";

    TTreeReader reader(chain);
    TTreeReaderArray<float> time_v(reader, "timevalue");
    TTreeReaderArray<float> amp_v (reader, "amplitude");

    Long64_t nProcessed = 0;

    while (reader.Next() && nProcessed < kNMax) {
        const float* A = &amp_v[0];
        const float* T = &time_v[0];

        // ── MCP1 quality gate ────────────────────────────────────────────
        Pulse mcp = ExtractPulse(T + kMCP1_t, A + kMCP1, 0.20f, kMCP1_minPeak);
        if (!mcp.valid) { ++nProcessed; continue; }
        if (mcp.peak < kMCP1_minPeak || mcp.peak > kMCP1_maxPeak)
            { ++nProcessed; continue; }

        // ── Wire chamber: require all four planes for fiducial check ─────
        // (Loose check — skip WC fiducial here; we want shape not position bias)
        // We do require WC reconstruction to suppress beam halo/noise
        Pulse wcR = ExtractPulse(T + kWC_t, A + kWC_R, 0.50f, kWC_minPeak);
        Pulse wcL = ExtractPulse(T + kWC_t, A + kWC_L, 0.50f, kWC_minPeak);
        Pulse wcD = ExtractPulse(T + kWC_t, A + kWC_D, 0.50f, kWC_minPeak);
        Pulse wcU = ExtractPulse(T + kWC_t, A + kWC_U, 0.50f, kWC_minPeak);
        bool wc_ok = wcR.valid && wcL.valid && wcD.valid && wcU.valid;
        if (!wc_ok) { ++nProcessed; continue; }

        // ── Fill MCP1 waveform (index 8) ─────────────────────────────────
        {
            Pulse mcp_cfd = ExtractPulse(T + kMCP1_t, A + kMCP1, 0.20f, kMCP1_minPeak);
            if (mcp_cfd.valid && mcp_cfd.crossingTime != kNoTime) {
                float t_cfd = mcp_cfd.crossingTime;
                float ped   = mcp_cfd.pedestal;
                const float* Aa = A + kMCP1;
                const float* Ta = T + kMCP1_t;
                for (int s = 0; s < 1024; ++s) {
                    float t_rel = Ta[s] - t_cfd;
                    if (t_rel < (float)kTMin || t_rel > (float)kTMax) continue;
                    float signal = ped - Aa[s];
                    hWF[8]->Fill(t_rel, signal);
                }
            }
        }

        // ── Fill HG capillary waveforms (indices 0-7) ────────────────────
        for (int i = 0; i < 8; ++i) {
            const float* Aa = A + kCap[i].hg;
            const float* Ta = T + kCap[i].hg_t;

            Pulse p = ExtractPulse(Ta, Aa, 0.20f, kHG_minPeak);
            if (!p.valid) continue;
            if (p.crossingTime == kNoTime) continue;

            float t_cfd = p.crossingTime;
            float ped   = p.pedestal;

            for (int s = 0; s < 1024; ++s) {
                float t_rel = Ta[s] - t_cfd;
                if (t_rel < (float)kTMin || t_rel > (float)kTMax) continue;
                float signal = ped - Aa[s];
                hWF[i]->Fill(t_rel, signal);
            }
        }

        // ── Fill LG capillary waveforms (indices 0-7) ────────────────────
        //    Persisted for the Layer-1 hero glance; not drawn in this macro's
        //    own PDFs.  Aligned on each LG channel's own CFD-20% crossing.
        for (int i = 0; i < 8; ++i) {
            const float* Aa = A + kCap[i].lg;
            const float* Ta = T + kCap[i].lg_t;

            Pulse p = ExtractPulse(Ta, Aa, 0.20f, kLG_minPeak);
            if (!p.valid) continue;
            if (p.crossingTime == kNoTime) continue;

            float t_cfd = p.crossingTime;
            float ped   = p.pedestal;

            for (int s = 0; s < 1024; ++s) {
                float t_rel = Ta[s] - t_cfd;
                if (t_rel < (float)kTMin_LG || t_rel > (float)kTMax_LG) continue;
                hLG[i]->Fill(t_rel, ped - Aa[s]);
            }
        }

        ++nProcessed;
    }

    std::cout << "[averageWaveforms] " << run.label
              << ": filled " << nProcessed << " events\n";

    delete chain;
    return (nProcessed > 0);
}

// ---------------------------------------------------------------------------
// DrawEnergyPage — 3×3 canvas for one energy, saved to PDF
// ---------------------------------------------------------------------------
static void DrawEnergyPage(TProfile* hWF[], const RunCfg& run, int iRun,
                            const char* pdfPath)
{
    TCanvas c("cWF", "", 2400, 2400);
    c.Divide(3, 3, 0.003, 0.035);

    // Channels 0-7: HG capillaries; channel 8: MCP1
    for (int i = 0; i < kNWF; ++i) {
        c.cd(i + 1);
        TString title = Form("%.0f GeV  %s",
                             run.energy_GeV, ChanName(i));
        DrawWFPad(hWF[i], title.Data(), iRun);
    }

    // Page title on the top pad (pad 0 == the full canvas frame)
    c.cd(0);
    DrawPageTitle(Form("Average HG Waveforms "
                       "#scale[0.85]{#font[52]{(mean #pm RMS)}}  "
                       "#bf{%.0f GeV}", run.energy_GeV));

    c.Print(pdfPath);
}

// ---------------------------------------------------------------------------
// DrawSummaryPage — 4×2 canvas, all energies overlaid, one panel per channel
// ---------------------------------------------------------------------------
static void DrawSummaryPage(TProfile* hAllWF[][kNWF], bool valid[],
                             const char* pdfPath)
{
    // 8 HG capillary channels only in the summary (skip MCP1 = index 8)
    TCanvas c("cSum", "", 3200, 1600);
    c.Divide(4, 2, 0.003, 0.035);

    for (int i = 0; i < 8; ++i) {
        c.cd(i + 1);
        StylePad();

        // Find global y range for this channel across energies
        double yhi_all = 0.;
        for (int r = 0; r < kNRuns; ++r) {
            if (!valid[r]) continue;
            int nb = hAllWF[r][i]->GetNbinsX();
            for (int b = 1; b <= nb; ++b) {
                double v = hAllWF[r][i]->GetBinContent(b);
                if (v > yhi_all) yhi_all = v;
            }
        }
        if (yhi_all <= 0.) yhi_all = 1.;

        bool first = true;
        for (int r = 0; r < kNRuns; ++r) {
            if (!valid[r]) continue;
            TProfile* h = hAllWF[r][i];

            // RMS band
            TProfile* hBand = (TProfile*)h->Clone(
                Form("sum_band_%s_%d", kCap[i].name, r));
            hBand->SetFillColorAlpha(kREnergyCols[r], 0.18);
            hBand->SetFillStyle(1001);
            hBand->SetLineColor(0);
            hBand->SetMarkerSize(0);
            hBand->GetYaxis()->SetRangeUser(-yhi_all * 0.1, yhi_all * 1.25);

            // Mean line
            TProfile* hLine = (TProfile*)h->Clone(
                Form("sum_line_%s_%d", kCap[i].name, r));
            hLine->SetLineColor(kREnergyCols[r]);
            hLine->SetLineWidth(2);
            hLine->SetMarkerSize(0);

            if (first) {
                hBand->Draw("E3");
                hLine->Draw("HIST SAME");
                first = false;
            } else {
                hBand->Draw("E3 SAME");
                hLine->Draw("HIST SAME");
            }
        }

        // CFD reference line
        TLine* lCFD = new TLine(0., -yhi_all * 0.1, 0., yhi_all * 1.25);
        lCFD->SetLineStyle(2);
        lCFD->SetLineColor(kRed);
        lCFD->SetLineWidth(1);
        lCFD->Draw("SAME");

        DrawPadTitle(kCap[i].name);   // from RADiCALStyle.h

        gPad->Modified();
        gPad->Update();
    }

    // Legend in the 8th pad (pad 8, lower-right corner)
    // We need to ensure pad 8 has at least the legend even if empty
    // Build legend separately on top of last pad
    c.cd(8);

    TLegend* leg = new TLegend(0.14, 0.18, 0.92, 0.82);
    leg->SetBorderSize(1);
    leg->SetTextSize(0.060);
    leg->SetFillStyle(0);
    for (int r = 0; r < kNRuns; ++r) {
        if (!valid[r]) continue;
        // Dummy line for legend entry
        TLine* dLine = new TLine(0, 0, 1, 1);
        dLine->SetLineColor(kREnergyCols[r]);
        dLine->SetLineWidth(2);
        leg->AddEntry(dLine,
            Form("%.0f GeV", kRuns[r].energy_GeV), "l");
    }
    leg->Draw();

    // Page title
    c.cd(0);
    DrawPageTitle("Average HG Waveforms  "
                  "#scale[0.85]{#font[52]{(mean #pm RMS)}}  "
                  "#bf{all energies}  "
                  "#scale[0.85]{150 GeV in red}");

    c.Print(pdfPath);
}

// ---------------------------------------------------------------------------
// Main entry point
// ---------------------------------------------------------------------------
void averageWaveforms()
{
    // ── Global style ─────────────────────────────────────────────────────────
    ApplyRADiCALStyle();   // RADiCALStyle.h via PlotUtils.h

    // ── Ensure output directories exist ─────────────────────────────────────
    for (int r = 0; r < kNRuns; ++r) {
        TString dir = Form("output/%s", kRuns[r].label.Data());
        gSystem->mkdir(dir, kTRUE);
    }
    gSystem->mkdir(kSumDir, kTRUE);

    // ── Storage for all-energy summary: hAllWF[energy][channel] ─────────────
    // Allocated on heap because 6×9 TProfiles are each ~300 bins
    TProfile* hAllWF[6][kNWF];
    TProfile* hAllLG[6][8];          // 8 LG capillaries (persisted for the hero glance)
    bool      validRun[6] = {};

    for (int r = 0; r < kNRuns; ++r) {
        for (int i = 0; i < kNWF; ++i) hAllWF[r][i] = nullptr;
        for (int i = 0; i < 8;    ++i) hAllLG[r][i] = nullptr;
    }

    // ── Per-energy loop ───────────────────────────────────────────────────────
    for (int r = 0; r < kNRuns; ++r) {
        const RunCfg& run = kRuns[r];

        std::cout << "\n[averageWaveforms] ====== " << run.label << " ======\n";

        // Book profiles for this energy
        TProfile* hWF[kNWF];
        for (int i = 0; i < 8; ++i)
            hWF[i] = BookWF(run.label.Data(), kCap[i].name);
        hWF[8] = BookWF(run.label.Data(), "MCP1");
        TProfile* hLG[8];
        for (int i = 0; i < 8; ++i)
            hLG[i] = BookWF(run.label.Data(), Form("%s_LG", kCap[i].name),
                            kNBins_LG, kTMin_LG, kTMax_LG);

        // Fill from raw data
        bool ok = ProcessEnergy(run, hWF, hLG);
        if (!ok) {
            for (int i = 0; i < kNWF; ++i) delete hWF[i];
            for (int i = 0; i < 8;    ++i) delete hLG[i];
            continue;
        }
        validRun[r] = true;

        // Save clones for the summary / hero (HG+MCP1 drawn here; LG persisted only)
        for (int i = 0; i < kNWF; ++i) {
            hAllWF[r][i] = (TProfile*)hWF[i]->Clone(
                Form("sum_%s_%s", run.label.Data(), ChanName(i)));
            hAllWF[r][i]->SetDirectory(nullptr);
        }
        for (int i = 0; i < 8; ++i) {
            hAllLG[r][i] = (TProfile*)hLG[i]->Clone(
                Form("sum_%s_%s_LG", run.label.Data(), kCap[i].name));
            hAllLG[r][i]->SetDirectory(nullptr);
        }

        // Write per-energy PDF (HG + MCP1 only — appendix layout unchanged)
        TString pdfPath = Form("output/%s/average_waveforms.pdf",
                               run.label.Data());
        DrawEnergyPage(hWF, run, r, pdfPath.Data());
        std::cout << "[averageWaveforms] Written: " << pdfPath << "\n";

        // Clean up per-energy profiles
        for (int i = 0; i < kNWF; ++i) delete hWF[i];
        for (int i = 0; i < 8;    ++i) delete hLG[i];
    }

    // ── Summary page ──────────────────────────────────────────────────────────
    bool anyValid = false;
    for (int r = 0; r < kNRuns; ++r) anyValid |= validRun[r];

    if (anyValid) {
        TString sumPath = Form("%saverage_waveforms_summary.pdf", kSumDir);
        DrawSummaryPage(hAllWF, validRun, sumPath.Data());
        std::cout << "[averageWaveforms] Written: " << sumPath << "\n";
    } else {
        std::cerr << "[averageWaveforms] No valid energies — summary skipped.\n";
    }

    // ── Persist average-waveform profiles for the Layer-1 hero plot ─────────
    //   Keys: sum_<label>_<chan> (e.g. sum_150GeV_NW-D).  layer1Summary.C reads
    //   these to build the clean "average pulse shape per capillary" figure.
    if (anyValid) {
        TFile fout(Form("%saverage_waveforms.root", kSumDir), "RECREATE");
        for (int r = 0; r < kNRuns; ++r) {
            for (int i = 0; i < kNWF; ++i) if (hAllWF[r][i]) hAllWF[r][i]->Write();
            for (int i = 0; i < 8;    ++i) if (hAllLG[r][i]) hAllLG[r][i]->Write();
        }
        fout.Close();
        std::cout << "[averageWaveforms] Wrote " << kSumDir
                  << "average_waveforms.root\n";
    }

    // ── Cleanup ───────────────────────────────────────────────────────────────
    for (int r = 0; r < kNRuns; ++r) {
        for (int i = 0; i < kNWF; ++i) if (hAllWF[r][i]) delete hAllWF[r][i];
        for (int i = 0; i < 8;    ++i) if (hAllLG[r][i]) delete hAllLG[r][i];
    }

    std::cout << "\n[averageWaveforms] Done.\n";
}

// ===========================================================================
// transverseMaps.C  —  per-channel transverse (x,y) amplitude maps, BEFORE and
// AFTER the 1x1-cm trigger selection.
//
// For every capillary (and the 1x1 trigger itself) this fills the mean pulse
// amplitude versus beam-track impact point:
//     map[ch]->Fill(x_trk, y_trk, peak[ch])
// exactly as the raw scan does, with two event selections:
//     BEFORE : good track only (all four WC planes valid)
//     AFTER  : good track AND the 1x1-cm trigger fired (ch15 peak > kCh15_trig)
//
// The 1x1 trigger (DRS0 grp1 ch6, the "1x1 cm" detector) tags a clean, single
// beam particle through the trigger aperture; requiring it removes the no-hit
// pedestal population (~20%) and the MCP BNC/SMA connector-shower events that
// show up as fixed rings in the ch15 map at ~(7,5) and ~(20,3) mm.
//
// Binning: x at 0.25 mm (x_trk is continuous), y at 0.5 mm (suppresses the
// peak-sample comb in y) — "finest but no finer".  Palette: inverted kRust
// (perceptually uniform, colour-blind safe).
//
// Reads raw waveforms (Data/ via kRuns, like drs4TimeBase.C).  No reprocessing.
//
// Run:  ROOT_INCLUDE_PATH=Analysis root -l -b -q 'Analysis/transverseMaps.C+'
//       (optional arg: energy index 0..5, default -1 = all energies)
// ===========================================================================
#include "ChannelConfig.h"
#include "WaveformUtils.h"
#include "PlotUtils.h"      // StylePad, RADiCAL style

#include "TChain.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TProfile2D.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TColor.h"

#include <vector>

namespace {
const int    kCh15    = chanOff(0,1,6);   // "Ch 15 - DRS 0" = 1x1 cm trigger
const int    kCh15_t  = kT_D0G1;
const float  kCh15_trig = 60.0f;          // mV: 1x1 trigger threshold (valley
                                          // between no-hit pedestal ~40 mV and
                                          // the fired-signal continuum)
// Map axes — symmetric 40 mm window CENTRED ON THE MODULE (kCalo_x0/y0, the
// calorimeter face centre in WC coords), so the module sits at the frame centre.
// Absolute WC coordinates are kept (connector positions stay (7,5)/(20,3) mm).
const double kHalfWin = 20.;                                       // mm
const int kNX = 160; const double kX0 = kCalo_x0 - kHalfWin, kX1 = kCalo_x0 + kHalfWin;  // 0.25 mm/bin
const int kNY =  80; const double kY0 = kCalo_y0 - kHalfWin, kY1 = kCalo_y0 + kHalfWin;  // 0.50 mm/bin

void drawGrid(TCanvas& c, TProfile2D* m[9], const char* tag, const char* outPDF)
{
    // 3x3 grid: 8 capillaries + the 1x1 trigger (cell 9).
    c.Clear();
    c.Divide(3, 3, 0.004, 0.004);
    for (int i = 0; i < 9; ++i) {
        c.cd(i + 1);
        gPad->SetRightMargin(0.15);
        gPad->SetLeftMargin(0.13);
        gPad->SetTicks(1, 1);
        if (!m[i]) continue;
        m[i]->SetMinimum(0.);
        m[i]->Draw("COLZ");
    }
    c.cd(0);
    TLatex t; t.SetNDC(); t.SetTextSize(0.022); t.SetTextAlign(22);
    t.DrawLatex(0.5, 0.985, tag);
    c.Print(outPDF);
}
} // namespace

void transverseMaps(int energyIndex = -1)
{
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kRust); TColor::InvertPalette();  // inverted rust: signal dark on white
    TString sumDir = "output/Summary";
    gSystem->mkdir(sumDir, kTRUE);

    const int e0 = (energyIndex >= 0) ? energyIndex : 0;
    const int e1 = (energyIndex >= 0) ? energyIndex + 1 : kNRuns;

    for (int e = e0; e < e1; ++e) {
        TChain ch("pulse");
        TObjArray* parts = kRuns[e].inFiles.Tokenize(";");
        for (int i = 0; i < parts->GetEntries(); ++i)
            ch.Add(((TObjString*)parts->At(i))->GetString());
        delete parts;
        if (ch.GetEntries() == 0) {
            std::printf("[transverseMaps] %s: no raw files — skipping\n", kRuns[e].label.Data());
            continue;
        }

        TTreeReader r(&ch);
        TTreeReaderArray<float> tv(r, "timevalue"), av(r, "amplitude");

        // 9 maps (8 caps + 1x1 trigger), before and after the trigger cut.
        TProfile2D *bef[9], *aft[9];
        const char* nm[9];
        for (int i = 0; i < 8; ++i) nm[i] = kCap[i].name;
        nm[8] = "1x1 trigger";
        for (int i = 0; i < 9; ++i) {
            bef[i] = new TProfile2D(Form("bef_%s_%d", kRuns[e].label.Data(), i),
                Form("%s;x Track (mm);y Track (mm)", nm[i]), kNX, kX0, kX1, kNY, kY0, kY1);
            aft[i] = new TProfile2D(Form("aft_%s_%d", kRuns[e].label.Data(), i),
                Form("%s;x Track (mm);y Track (mm)", nm[i]), kNX, kX0, kX1, kNY, kY0, kY1);
            bef[i]->SetDirectory(nullptr); aft[i]->SetDirectory(nullptr);
            bef[i]->GetZaxis()->SetTitle("mean amplitude (mV)");
            aft[i]->GetZaxis()->SetTitle("mean amplitude (mV)");
        }

        long ngood = 0, ntrig = 0;
        while (r.Next()) {
            const float* T = &tv[0];
            const float* A = &av[0];
            Pulse R = ExtractPulse(T + kWC_t, A + kWC_R, 0.50f, kWC_minPeak);
            Pulse L = ExtractPulse(T + kWC_t, A + kWC_L, 0.50f, kWC_minPeak);
            Pulse D = ExtractPulse(T + kWC_t, A + kWC_D, 0.50f, kWC_minPeak);
            Pulse U = ExtractPulse(T + kWC_t, A + kWC_U, 0.50f, kWC_minPeak);
            if (!(R.valid && L.valid && D.valid && U.valid)) continue;   // good track
            ++ngood;
            float x = kWC_Scale * (R.peakTime - L.peakTime);
            float y = kWC_Scale * (D.peakTime - U.peakTime);

            float pk[9];
            for (int i = 0; i < 8; ++i)
                pk[i] = ExtractPulse(T + kCap[i].hg_t, A + kCap[i].hg, 0.20f, 5.0f).peak;
            pk[8] = ExtractPulse(T + kCh15_t, A + kCh15, 0.20f, 5.0f).peak;

            bool trig = (pk[8] > kCh15_trig);
            if (trig) ++ntrig;
            for (int i = 0; i < 9; ++i) {
                bef[i]->Fill(x, y, pk[i]);
                if (trig) aft[i]->Fill(x, y, pk[i]);
            }
        }
        std::printf("[transverseMaps] %s: good tracks=%ld, 1x1-trigger=%ld (%.1f%%)\n",
                    kRuns[e].label.Data(), ngood, ntrig,
                    ngood ? 100.*ntrig/ngood : 0.);

        TString pdf = sumDir + Form("/transverse_maps_%s.pdf", kRuns[e].label.Data());
        TCanvas c("c_tmap", "", 1800, 1700);
        c.Print(pdf + "[");                                   // open multipage
        drawGrid(c, bef, Form("%s  --  BEFORE 1x1 trigger (good track only)", kRuns[e].label.Data()), pdf);
        drawGrid(c, aft, Form("%s  --  AFTER 1x1 trigger (ch15 > %.0f mV)", kRuns[e].label.Data(), kCh15_trig), pdf);
        c.Print(pdf + "]");                                   // close multipage
        std::printf("[transverseMaps] wrote %s\n", pdf.Data());

        for (int i = 0; i < 9; ++i) { delete bef[i]; delete aft[i]; }
    }
}

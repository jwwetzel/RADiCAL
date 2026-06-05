// shadowMap.C — THE money plot. Instead of just counting where particles went,
// we colour each position by the MEAN pulse amplitude seen there (a TProfile2D).
// That turns the beam into an X-ray: dense material in front of a detector
// absorbs/scatters the beam, leaving a darker "shadow"; the calorimeter's edges
// and the two connector rings on the 1x1 scintillator jump straight out.
//
// We plot two channels side by side:
//   left  = the 1x1 cm scintillator (slot 15) -> two connector rings + the calo shadow
//   right = an NE-D timing capillary (slot 2) -> the calorimeter lit up in the beam
//
// Run:  root -l 'student/shadowMap.C'   (compile for speed: '...shadowMap.C+')
// Saves shadow_map.png.
#include "radlab.h"
#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TProfile2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TColor.h"
#include <cstdio>

void shadowMap(const char* fname = "datasets/2023/raw/RUN1258_150_GeV.root")
{
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);

    TFile* f = TFile::Open(fname);
    if (!f || f->IsZombie()){ printf("could not open %s\n", fname); return; }
    TTree* pulse = (TTree*)f->Get("pulse");

    TTreeReader r(pulse);
    TTreeReaderArray<float> amp(r, "amplitude");
    TTreeReaderArray<float> tim(r, "timevalue");

    // 40 mm window CENTRED on the calorimeter face, so the module sits in the middle.
    const double W = 20.0;
    TProfile2D* m1x1 = new TProfile2D("m1x1",
        "1x1 scintillator: mean amplitude;x (mm);y (mm)",
        160, CALO_X0-W, CALO_X0+W, 80, CALO_Y0-W, CALO_Y0+W);
    TProfile2D* mcap = new TProfile2D("mcap",
        "NE-D capillary: mean amplitude;x (mm);y (mm)",
        160, CALO_X0-W, CALO_X0+W, 80, CALO_Y0-W, CALO_Y0+W);

    long ngood = 0;
    while (r.Next()){
        const float* A = &amp[0];
        const float* T = &tim[0];

        const float* aR = A + ampOffset(WC_R); const float* aL = A + ampOffset(WC_L);
        const float* aD = A + ampOffset(WC_D); const float* aU = A + ampOffset(WC_U);
        const float* tw = T + timeOffset(WC_R);
        if (findPeak(aR) < 20 || findPeak(aL) < 20 ||
            findPeak(aD) < 20 || findPeak(aU) < 20) continue;        // good track only

        float x = WC_SCALE * (findPeakTime(aR, tw) - findPeakTime(aL, tw));
        float y = WC_SCALE * (findPeakTime(aD, tw) - findPeakTime(aU, tw));

        m1x1->Fill(x, y, findPeak(A + ampOffset(SLOT_1x1)));         // mean 1x1 amplitude here
        mcap->Fill(x, y, findPeak(A + ampOffset(HG_SLOT[1])));       // mean NE-D amplitude here
        ++ngood;
    }
    printf("filled maps from %ld good tracks\n", ngood);

    TCanvas* c = new TCanvas("c_sh", "", 1400, 620);
    c->Divide(2, 1, 0.005, 0.005);
    c->cd(1); gPad->SetRightMargin(0.15); m1x1->SetMinimum(0); m1x1->Draw("COLZ");
    c->cd(2); gPad->SetRightMargin(0.15); mcap->SetMinimum(0); mcap->Draw("COLZ");
    c->SaveAs("shadow_map.png");
}

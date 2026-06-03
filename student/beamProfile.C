// beamProfile.C — reconstruct the beam impact point from the wire chamber and
// histogram WHERE the particles went (a 2D beam spot).
//
// The wire chamber has 4 planes (Right, Left, Down, Up). A particle's position
// is encoded in the DIFFERENCE of arrival times between opposite planes:
//     x = scale * (t_Right - t_Left),   y = scale * (t_Down - t_Up)
//
// Run:  root -l 'student/beamProfile.C'   (compile for speed: '...beamProfile.C+')
// Saves beam_profile.png.
#include "radlab.h"
#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include <cstdio>

void beamProfile(const char* fname = "Data/RUN1258_150_GeV.root")
{
    gStyle->SetOptStat(0);
    TFile* f = TFile::Open(fname);
    if (!f || f->IsZombie()){ printf("could not open %s\n", fname); return; }
    TTree* pulse = (TTree*)f->Get("pulse");

    TTreeReader r(pulse);
    TTreeReaderArray<float> amp(r, "amplitude");
    TTreeReaderArray<float> tim(r, "timevalue");

    TH2F* h = new TH2F("h", "beam profile (wire chamber);x (mm);y (mm)",
                       120, -25, 25, 120, -25, 25);

    long ngood = 0;
    while (r.Next()){
        const float* A = &amp[0];
        const float* T = &tim[0];

        const float* aR = A + ampOffset(WC_R); const float* aL = A + ampOffset(WC_L);
        const float* aD = A + ampOffset(WC_D); const float* aU = A + ampOffset(WC_U);
        const float* tw = T + timeOffset(WC_R);   // all four planes share DRS1-group1 time axis

        // require a real signal (> 20 mV) in all four planes = a good track
        if (findPeak(aR) < 20 || findPeak(aL) < 20 ||
            findPeak(aD) < 20 || findPeak(aU) < 20) continue;

        float x = WC_SCALE * (findPeakTime(aR, tw) - findPeakTime(aL, tw));
        float y = WC_SCALE * (findPeakTime(aD, tw) - findPeakTime(aU, tw));
        h->Fill(x, y);
        ++ngood;
    }
    printf("good tracks: %ld\n", ngood);

    TCanvas* c = new TCanvas("c_bp", "", 700, 650);
    c->SetRightMargin(0.14);
    h->Draw("COLZ");
    c->SaveAs("beam_profile.png");
}

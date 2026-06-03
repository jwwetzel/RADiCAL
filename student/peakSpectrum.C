// peakSpectrum.C — your first event loop. Read every event, measure the peak
// height of one capillary, and histogram it.
//
// Run:  root -l 'student/peakSpectrum.C(2)'        (slot 2 = NE-D capillary)
//   For speed on a 2 GB file you can compile it:  root -l 'student/peakSpectrum.C+(2)'
// Saves peak_spectrum.png.
#include "radlab.h"
#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TH1F.h"
#include "TCanvas.h"
#include <cstdio>

void peakSpectrum(int slot = 2, const char* fname = "Data/RUN1258_150_GeV.root")
{
    TFile* f = TFile::Open(fname);
    if (!f || f->IsZombie()){ printf("could not open %s\n", fname); return; }
    TTree* pulse = (TTree*)f->Get("pulse");

    TTreeReader r(pulse);
    TTreeReaderArray<float> amp(r, "amplitude");

    TH1F* h = new TH1F("h", "capillary peak height;peak amplitude (mV);events", 100, 0, 1000);

    long n = 0;
    while (r.Next()){
        const float* a = &amp[0] + ampOffset(slot);
        h->Fill(findPeak(a));
        ++n;
    }
    printf("looped over %ld events\n", n);

    TCanvas* c = new TCanvas("c_pk", "", 800, 550);
    h->SetLineColor(kAzure+1); h->SetLineWidth(2); h->Draw();
    c->SaveAs("peak_spectrum.png");
}

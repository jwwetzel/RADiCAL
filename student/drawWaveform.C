// drawWaveform.C — draw ONE raw waveform so you can see what the detector records.
//
// Run (interpreted):  root -l 'student/drawWaveform.C(2, 0)'
//   first arg = slot (2 = the NE-D timing capillary), second arg = event number.
// Saves waveform.png and leaves it on screen.
#include "radlab.h"
#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TGraph.h"
#include "TCanvas.h"
#include <cstdio>

void drawWaveform(int slot = 2, long event = 0,
                  const char* fname = "Data/RUN1258_150_GeV.root")
{
    TFile* f = TFile::Open(fname);
    if (!f || f->IsZombie()){ printf("could not open %s\n", fname); return; }
    TTree* pulse = (TTree*)f->Get("pulse");

    TTreeReader r(pulse);
    TTreeReaderArray<float> amp(r, "amplitude");
    TTreeReaderArray<float> tim(r, "timevalue");

    // walk forward to the event we want
    long i = 0; bool found = false;
    while (r.Next()){ if (i == event){ found = true; break; } ++i; }
    if (!found){ printf("event %ld not found\n", event); return; }

    const float* a = &amp[0] + ampOffset(slot);    // this slot's 1024 samples
    const float* t = &tim[0] + timeOffset(slot);   // its time axis

    TGraph* g = new TGraph(1024);
    for (int s = 0; s < 1024; ++s) g->SetPoint(s, t[s], a[s]);

    g->SetTitle(Form("raw waveform  -  slot %d, event %ld;time (ns);amplitude (mV)", slot, event));
    g->SetLineColor(kAzure+1); g->SetLineWidth(2);
    TCanvas* c = new TCanvas("c_wf", "", 900, 500);
    g->Draw("AL");
    c->SaveAs("waveform.png");
    printf("pedestal = %.1f mV,  peak height = %.1f mV\n", findPedestal(a), findPeak(a));
}

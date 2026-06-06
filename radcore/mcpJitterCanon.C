// ============================================================================
// mcpJitterCanon.C — what is the "100 ps MCP jitter", really?  (canonical data)
// ----------------------------------------------------------------------------
// MCP1 and MCP2 are ONE MCP signal passively split into the two DT5742 groups
// (corr=1.000). So in  D = mcp_time - mcp2_time  the MCP's OWN jitter (Hamamatsu
// ~10 ps) and the beam arrival T both CANCEL (common-mode). What is left is the
// INTER-GROUP DRS4 timing jitter:
//     D = [d0(stopcell0) - d1(stopcell1)] + [n0 - n1] + cable_offset
//   d_k = group-k domino-wave / stop-cell digitisation timing
//   n_k = CFD/baseline noise on the MCP copy in group k
// Tests, vs energy:
//   (1) sigma(D)/sqrt2  = per-group reference jitter           (expect ~70 ps, FLAT)
//   (2) is it FLAT vs MCP amplitude?  -> clock/stop-cell, not slew
//   (3) OOS additive stop-cell correction d_k(stopcell_k):  does it COLLAPSE?
//       If 70 ps -> ~15 ps, the floor is the DRS4 domino phase (correctable,
//       and same-group-cancelling), NOT the MCP.  The residual bounds the MCP
//       copy's CFD/readout noise; the true 10 ps MCP jitter is invisible here
//       because it is common-mode in the split.
//   ROOT_INCLUDE_PATH=radcore:Analysis root -l -b -q radcore/mcpJitterCanon.C+
// ============================================================================
#include "RadTiming.h"       // rad::tebSigma  (ns vector -> ps core sigma)
#include "DataPaths.h"
#include "SelectionCuts.h"
#include "PlotUtils.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>

struct Mev { float d; int sc0, sc1; float amp; char par; };

static double medianOf(std::vector<float>& v){
    if(v.empty()) return 0; size_t k=v.size()/2; std::nth_element(v.begin(),v.begin()+k,v.end());
    return v[k];
}
static double sigPs(std::vector<float>& v){ return rad::tebSigma(v); }  // /sqrt2 applied by caller

void mcpJitterCanon(){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    const double Es[6]={25,50,75,100,125,150};
    double gRaw[6], gCor[6], gAmpFlat[6];
    const int NB=32;   // stop-cell bins (1024/32 = 32 cells/bin)

    for(int e=0;e<6;++e){
        TFile* fp=TFile::Open(radReduced("DSB1",Es[e])); TTree* t=(TTree*)fp->Get("rad");
        Bool_t wc; Float_t mp,mp2,mt,mt2; Int_t sc[4];
        t->SetBranchAddress("wc_ok",&wc);
        t->SetBranchAddress("mcp_peak",&mp);   t->SetBranchAddress("mcp2_peak",&mp2);
        t->SetBranchAddress("mcp_time",&mt);   t->SetBranchAddress("mcp2_time",&mt2);
        t->SetBranchAddress("stopcell",sc);
        long N=t->GetEntries(); std::vector<Mev> ev; ev.reserve(N); long ie=0;
        for(long i=0;i<N;++i){ t->GetEntry(i);
            if(!wc) continue;
            if(mp<80||mp>450||mp2<80||mp2>450) continue;          // valid, unsaturated MCP copies
            if(mt<-1e4||mt2<-1e4) continue;
            Mev m; m.d=mt-mt2; m.sc0=sc[0]; m.sc1=sc[1]; m.amp=0.5f*(mp+mp2); m.par=(char)(ie++&1);
            if(std::fabs(m.d)>5.0f) continue;                     // ns; reject pathological
            ev.push_back(m);
        }
        // (1) raw per-group jitter
        std::vector<float> draw; for(auto&m:ev) draw.push_back(m.d);
        gRaw[e]= sigPs(draw)/std::sqrt(2.0);

        // (3) OOS additive stop-cell correction: train on par==0, apply to par==1
        std::vector<float> off0(NB,0), off1(NB,0);
        { std::vector<std::vector<float>> b0(NB);
          for(auto&m:ev) if(m.par==0){ int b=m.sc0*NB/1024; if(b<0)b=0; if(b>=NB)b=NB-1; b0[b].push_back(m.d); }
          for(int b=0;b<NB;++b) off0[b]= b0[b].size()>30 ? medianOf(b0[b]) : 0.0; }
        // residual after off0, then train off1
        { std::vector<std::vector<float>> b1(NB);
          for(auto&m:ev) if(m.par==0){ int b0i=m.sc0*NB/1024; if(b0i<0)b0i=0; if(b0i>=NB)b0i=NB-1;
                float r=m.d-off0[b0i]; int b=m.sc1*NB/1024; if(b<0)b=0; if(b>=NB)b=NB-1; b1[b].push_back(r); }
          for(int b=0;b<NB;++b) off1[b]= b1[b].size()>30 ? medianOf(b1[b]) : 0.0; }
        std::vector<float> dcor;
        for(auto&m:ev) if(m.par==1){ int a=m.sc0*NB/1024,b=m.sc1*NB/1024;
            if(a<0)a=0;if(a>=NB)a=NB-1;if(b<0)b=0;if(b>=NB)b=NB-1;
            dcor.push_back(m.d-off0[a]-off1[b]); }
        gCor[e]= sigPs(dcor)/std::sqrt(2.0);

        // (2) amplitude dependence: per-group jitter in the bright half vs dim half
        std::vector<Mev> s=ev; std::sort(s.begin(),s.end(),[](const Mev&A,const Mev&B){return A.amp<B.amp;});
        std::vector<float> dHi; for(size_t i=s.size()/2;i<s.size();++i) dHi.push_back(s[i].d);
        gAmpFlat[e]= sigPs(dHi)/std::sqrt(2.0);   // bright-half per-group jitter (compare to gRaw -> flat?)

        printf("E=%3.0f  per-group sigma: raw=%.1f ps  bright-half=%.1f ps  ->  stop-cell corrected=%.1f ps   (N=%zu)\n",
               Es[e], gRaw[e], gAmpFlat[e], gCor[e], ev.size());
        fp->Close();
    }

    TCanvas* c=new TCanvas("c_mcp","",900,680);
    TH1F* fr=c->DrawFrame(0,0,160,90);
    fr->SetTitle("DT5742 inter-group timing jitter (the \"MCP jitter\" is really the DRS4 domino phase);beam energy (GeV);#sigma(MCP1#minusMCP2)/#sqrt{2}  (ps)");
    TGraph* gr=new TGraph(6,Es,gRaw); gr->SetMarkerStyle(24); gr->SetMarkerColor(kGray+2); gr->SetLineColor(kGray+2); gr->SetMarkerSize(1.6); gr->SetLineWidth(2); gr->Draw("PL SAME");
    TGraph* gc=new TGraph(6,Es,gCor); gc->SetMarkerStyle(20); gc->SetMarkerColor(kRed+1);  gc->SetLineColor(kRed+1);  gc->SetMarkerSize(1.7); gc->SetLineWidth(2); gc->Draw("PL SAME");
    TLegend* lg=new TLegend(0.16,0.74,0.74,0.88); lg->SetBorderSize(0);
    lg->AddEntry(gr,"raw per-group jitter (uncorrected)","pl");
    lg->AddEntry(gc,"after OOS stop-cell correction","pl"); lg->Draw();
    TLatex tx; tx.SetNDC(); tx.SetTextFont(42); tx.SetTextSize(0.030); tx.SetTextColor(kGray+3);
    tx.DrawLatex(0.16,0.20,"MCP1/MCP2 = one split signal: the MCP's own ~10 ps jitter is common-mode and cancels here.");
    tx.DrawLatex(0.16,0.16,"What remains is the per-group DRS4 stop-cell / domino-wave phase -- correctable, and");
    tx.DrawLatex(0.16,0.12,"cancelled exactly in same-group estimators (the headline (DW#minusUP)/2 corners).");
    gSystem->mkdir("radcore/figs",kTRUE); c->Print("radcore/figs/mcp_jitter_canon.png");
    printf("wrote radcore/figs/mcp_jitter_canon.png\n");
}

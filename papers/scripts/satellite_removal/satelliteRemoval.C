// ============================================================================
// satelliteRemoval.C — the satellite/tail-removal demonstration for §5.3.
// Protocol: AUDIT.md (this directory). Identical-event sample (method-gain
// design): intersection validity+veto across {srCFD, cfd05, LED}, brightest-1000
// once. Observable: per-event (DW-UP)/2 (the production-width quantity).
// Metrics: core Gaussian sigma, production tebSigma, robust IQR sigma, and the
// COMMON-window tail fraction f_tail(|x-med| > 2.5*sigma_core(srCFD)).
// Panels: A distributions @150 (log-y, unit area); B metrics @150;
//         C clipped (nsat>=7) vs less-clipped (nsat<=2) tail fractions @75 GeV.
//   source setup.sh
//   root -l -b -q 'papers/scripts/satellite_removal/satelliteRemoval.C+'
// Output: papers/figures/satellite_removal/satellite_removal.{png,pdf} + caption.
// ============================================================================
#include "RadView.h"
#include "RadTiming.h"
#include "DataPaths.h"
#include "PlotUtils.h"
#include "SelectionCuts.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TBox.h"
#include "TStyle.h"
#include "TSystem.h"
#include <vector>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <cstdio>
using namespace rad;

static const int NS=3;
static const int  SRC[NS]={RadView::kLGCFD,RadView::kCFD05,RadView::kLED};
static const char* SNM[NS]={"srCFD","cfd05","LED"};

struct Ev { float dwup[NS]; float slg; int nsat; };

static std::vector<Ev> gather(double E, BuildConfig& cfg){
    std::vector<Ev> out;
    TString p=radReduced("DSB1",E); if(gSystem->AccessPathName(p.Data())) return out;
    TFile* f=TFile::Open(p.Data()); if(!f||f->IsZombie()){if(f)f->Close();return out;}
    TTree* t=(TTree*)f->Get("rad"); if(!t){f->Close();return out;}
    RadView v; v.attach(t,&cfg); double xc,yc; v.beamCenter(xc,yc);
    double R=TimingFiducialR(E), r2=R*R;
    Bool_t sat[8]; t->SetBranchAddress("hg_saturated",sat);
    Long64_t N=v.entries(); out.reserve(N/8);
    for(Long64_t i=0;i<N;++i){ v.get(i);
        if(!v.wc_ok()||v.mcp1_peak()<kMCP1_minPeak||v.mcp1_peak()>kMCP1_maxPeak) continue;
        double dx=v.x_trk()-xc,dy=v.y_trk()-yc; if(dx*dx+dy*dy>=r2) continue;
        Ev e; bool allok=true;
        for(int s=0;s<NS&&allok;++s){
            float tt[8];
            for(int c=0;c<8;++c){ if(v.hg_peak(c)<kHG_minPeak){allok=false;break;}
                float tc=v.timeOf(c,SRC[s]); if(tc<=-1e5f){allok=false;break;} tt[c]=tc; }
            if(!allok) break;
            float a[8]; for(int c=0;c<8;++c)a[c]=tt[c];
            std::nth_element(a,a+4,a+8); float med=a[4];
            for(int c=0;c<8;++c) if(std::fabs(tt[c]-med)>=kTimingChanConsistency_ns){allok=false;break;}
            if(!allok) break;
            double ds=0,us=0; for(int c=0;c<4;++c)ds+=tt[c]; for(int c=4;c<8;++c)us+=tt[c];
            e.dwup[s]=0.5f*(float)(ds/4.0-us/4.0);
        }
        if(!allok) continue;
        e.slg=(float)v.sum_lg(); e.nsat=0; for(int c=0;c<8;++c) if(sat[c])++e.nsat;
        out.push_back(e);
    }
    f->Close(); return out;
}

struct M { double core=0,prod=0,rob=0,ftail=0,med=0; long n=0,ntail=0; };
static double medianOf(std::vector<float>& v){ std::vector<float> a=v; std::nth_element(a.begin(),a.begin()+a.size()/2,a.end()); return a[a.size()/2]; }
static M metrics(std::vector<float>& v,double W){
    M m; if(v.size()<200) return m; m.n=v.size(); m.med=medianOf(v);
    std::vector<float> t=v; m.prod=tebSigma(t);
    std::vector<float> a=v; std::sort(a.begin(),a.end());
    m.rob=(a[(size_t)(0.75*a.size())]-a[(size_t)(0.25*a.size())])/1.349*1000.0;
    double lo=m.med-0.4,hi=m.med+0.4;
    TH1F h("hsm","",160,lo,hi); h.SetDirectory(nullptr); for(float x:v)h.Fill(x);
    double mu,muE,s,sE; FitGaussCore(&h,2.0,mu,muE,s,sE); m.core=s*1000.0;
    for(float x:v) if(std::fabs(x-m.med)>W) ++m.ntail;
    m.ftail=100.0*m.ntail/m.n;
    return m;
}

void satelliteRemoval(){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0); gStyle->SetGridStyle(3); gStyle->SetGridColor(kGray);
    gSystem->mkdir("papers/figures/satellite_removal",kTRUE);
    BuildConfig cfg=BuildConfig::Load(radConfig("DSB1").Data());
    printf("\n===== SATELLITE / TAIL REMOVAL DEMONSTRATION (identical events) =====\n");

    // ---------- 150 GeV: panels A + B ----------
    std::vector<Ev> ev=gather(150,cfg);
    std::nth_element(ev.begin(),ev.begin()+1000,ev.end(),[](const Ev&a,const Ev&b){return a.slg>b.slg;});
    ev.resize(1000);
    std::vector<float> vt[NS]; for(int s=0;s<NS;++s) for(Ev&e:ev) vt[s].push_back(e.dwup[s]);
    // common window from srCFD core
    std::vector<float> t0=vt[0]; double med0=medianOf(t0);
    TH1F h0("h0","",160,med0-0.4,med0+0.4); h0.SetDirectory(nullptr); for(float x:vt[0])h0.Fill(x);
    double mu,muE,s0,sE; FitGaussCore(&h0,2.0,mu,muE,s0,sE);
    const double W=2.5*s0;   // ns
    M m150[NS]; for(int s=0;s<NS;++s) m150[s]=metrics(vt[s],W);
    printf("  150 GeV, N=1000 identical; common tail window |x-med| > 2.5*sigma_core(srCFD) = %.1f ps\n",W*1000);
    printf("  %-6s core(ps) prod(ps) robust(ps) f_tail(%%) n_tail\n","");
    for(int s=0;s<NS;++s) printf("  %-6s %7.1f %8.1f %9.1f %8.2f %6ld\n",
        SNM[s],m150[s].core,m150[s].prod,m150[s].rob,m150[s].ftail,m150[s].ntail);
    printf("  Gaussian expectation for f_tail: ~1.24%%\n");

    // ---------- 75 GeV: panel C (clip split; pre-registered fallback) ----------
    // FULL-FIDUCIAL sample at 75 GeV: in the brightest-1000 the less-clipped subset is
    // empty (N=2 -- bright showers clip by construction), so the split uses all fiducial
    // events (documented deviation under the audit's too-small-subset rule).
    std::vector<Ev> e75=gather(75,cfg);
    std::vector<float> c75[NS],u75[NS];
    long nC=0,nU=0;
    for(Ev&e:e75){ if(e.nsat>=7){++nC; for(int s=0;s<NS;++s)c75[s].push_back(e.dwup[s]);}
                   else if(e.nsat<=2){++nU; for(int s=0;s<NS;++s)u75[s].push_back(e.dwup[s]);} }
    // common windows per subset from srCFD core of that subset
    auto coreOf=[&](std::vector<float>& v)->double{ if(v.size()<200) return -1;
        double md=medianOf(v); TH1F h("hh","",160,md-0.4,md+0.4); h.SetDirectory(nullptr);
        for(float x:v)h.Fill(x); double mu2,mE,s2,se2; FitGaussCore(&h,2.0,mu2,mE,s2,se2); return s2; };
    double Wc=2.5*coreOf(c75[0]), Wu=2.5*coreOf(u75[0]);
    M mc[NS],mu75[NS];
    for(int s=0;s<NS;++s){ mc[s]=metrics(c75[s],Wc); mu75[s]=metrics(u75[s],Wu); }
    printf("\n  75 GeV split (FULL FIDUCIAL; brightest-1000 left N=2 unclipped): clipped nsat>=7 N=%ld, less-clipped nsat<=2 N=%ld\n",nC,nU);
    printf("  f_tail clipped:      srCFD %.2f%%  cfd05 %.2f%%  LED %.2f%%\n",mc[0].ftail,mc[1].ftail,mc[2].ftail);
    printf("  f_tail less-clipped: srCFD %.2f%%  cfd05 %.2f%%  LED %.2f%%\n",mu75[0].ftail,mu75[1].ftail,mu75[2].ftail);

    // ---------- figure ----------
    TCanvas* c=new TCanvas("sr","",1640,620); c->Divide(3,1,0.005,0.005);
    // A: distributions @150, log-y
    c->cd(1); gPad->SetLeftMargin(0.14);gPad->SetRightMargin(0.03);gPad->SetTopMargin(0.085);gPad->SetBottomMargin(0.13);gPad->SetLogy();
    { double lo=-0.35,hi=0.35;
      double medC=medianOf(vt[1]);
      TH1F* hA=new TH1F("hA",";(DW#minusUP)/2 #minus median  (ns);normalized events",120,lo,hi); hA->SetDirectory(nullptr);
      TH1F* hB=new TH1F("hB","",120,lo,hi); hB->SetDirectory(nullptr);
      for(float x:vt[1])hA->Fill(x-medC); for(float x:vt[0])hB->Fill(x-med0);
      if(hA->Integral()>0)hA->Scale(1.0/hA->Integral()); if(hB->Integral()>0)hB->Scale(1.0/hB->Integral());
      hA->SetLineColor(kBlack); hA->SetLineWidth(2);
      hB->SetLineColor(kAzure+2); hB->SetFillColorAlpha(kAzure+2,0.25); hB->SetLineWidth(2);
      double ym=1.6*std::max(hA->GetMaximum(),hB->GetMaximum()); hA->SetMaximum(ym); hA->SetMinimum(2e-5);
      hA->GetYaxis()->SetTitleSize(0.05);hA->GetYaxis()->SetTitleOffset(1.35);hA->GetXaxis()->SetTitleSize(0.05);
      hA->Draw("HIST"); hB->Draw("HIST SAME");
      TLine* l1=new TLine(-W,2e-5,-W,ym*0.25); TLine* l2=new TLine(W,2e-5,W,ym*0.25);
      for(TLine* l:{l1,l2}){ l->SetLineColor(kRed+1);l->SetLineStyle(2);l->SetLineWidth(2);l->Draw(); }
      TLegend* lg=new TLegend(0.16,0.74,0.62,0.90); lg->SetBorderSize(0);lg->SetFillStyle(0);lg->SetTextSize(0.038);
      lg->AddEntry(hA,"cfd05 (clipped-peak CFD)","l"); lg->AddEntry(hB,"srCFD (recovered edge)","f"); lg->Draw();
      TLatex tx; tx.SetNDC(); tx.SetTextSize(0.032); tx.SetTextColor(kRed+2);
      tx.DrawLatex(0.16,0.69,Form("tail window: |x#minusmed| > %.0f ps",W*1000));
      DrawPadTitle("A  identical events, 150 GeV (log)"); }
    // B: metrics @150
    c->cd(2); gPad->SetLeftMargin(0.14);gPad->SetRightMargin(0.03);gPad->SetTopMargin(0.085);gPad->SetBottomMargin(0.13);
    { TLatex t4; t4.SetNDC(); t4.SetTextSize(0.042); double yy=0.87;
      auto L=[&](const char* s){ t4.DrawLatex(0.08,yy,s); yy-=0.072; };
      L("metrics @150 GeV (same 1000 events):");
      L(Form("core #sigma:    srCFD %.1f   cfd05 %.1f ps",m150[0].core,m150[1].core));
      L(Form("production:  srCFD %.1f   cfd05 %.1f ps",m150[0].prod,m150[1].prod));
      L(Form("robust IQR: srCFD %.1f   cfd05 %.1f ps",m150[0].rob,m150[1].rob));
      L(Form("f_{tail}:  srCFD %.2f%%   cfd05 %.2f%%",m150[0].ftail,m150[1].ftail));
      L("   (Gaussian expectation 1.24%)");
      L(Form("tail events: %ld #rightarrow %ld of 1000",m150[1].ntail,m150[0].ntail));
      L("#Rightarrow the gain lives in the non-Gaussian");
      L("    tail, not the Gaussian core");
      DrawPadTitle("B  width and tail metrics"); }
    // C: clip split @75
    c->cd(3); gPad->SetLeftMargin(0.14);gPad->SetRightMargin(0.03);gPad->SetTopMargin(0.085);gPad->SetBottomMargin(0.13);gPad->SetGridy();
    { TH1F* fr=gPad->DrawFrame(-0.5,0,1.5,1.35*std::max({mc[0].ftail,mc[1].ftail,mu75[0].ftail,mu75[1].ftail})+1);
      fr->SetTitle(";event subset;f_{tail} (%)");
      fr->GetXaxis()->SetBinLabel(fr->GetXaxis()->FindBin(0.0),"clipped (nsat#geq7)");
      fr->GetXaxis()->SetBinLabel(fr->GetXaxis()->FindBin(1.0),"less-clipped (nsat#leq2)");
      fr->GetXaxis()->SetLabelSize(0.055); fr->GetXaxis()->SetNdivisions(2);
      fr->GetYaxis()->SetTitleSize(0.05);
      double xs[2]={0,1};
      double yc[2]={mc[1].ftail,mu75[1].ftail}, ys[2]={mc[0].ftail,mu75[0].ftail};
      TGraph* g1=new TGraph(2,xs,yc); g1->SetMarkerStyle(24);g1->SetMarkerColor(kBlack);g1->SetMarkerSize(2.0);g1->Draw("P SAME");
      TGraph* g2=new TGraph(2,xs,ys); g2->SetMarkerStyle(20);g2->SetMarkerColor(kAzure+2);g2->SetMarkerSize(2.0);g2->Draw("P SAME");
      TLegend* lg=new TLegend(0.45,0.74,0.95,0.90); lg->SetBorderSize(0);lg->SetFillStyle(0);lg->SetTextSize(0.04);
      lg->AddEntry(g1,"cfd05","p"); lg->AddEntry(g2,"srCFD","p"); lg->Draw();
      TLatex tx; tx.SetNDC(); tx.SetTextSize(0.034); tx.SetTextColor(kGray+3);
      tx.DrawLatex(0.18,0.26,Form("75 GeV, full fiducial: N=%ld / %ld",nC,nU));
      tx.DrawLatex(0.18,0.20,"less-clipped (dim) events: srCFD tail RISES");
      tx.DrawLatex(0.18,0.15,"(noisy LG anchor) #rightarrow srCFD is for bright/clipped");
      DrawPadTitle("C  tail excess concentrates in clipped events"); }
    c->cd(0); DrawSuperTitle("Satellite/tail removal: srCFD vs cfd05 on identical events #minus the method gain is the suppression of clipping-induced tails",0.019f);
    c->Print("papers/figures/satellite_removal/satellite_removal.png");
    c->Print("papers/figures/satellite_removal/satellite_removal.pdf");
    printf("\n  wrote papers/figures/satellite_removal/satellite_removal.{png,pdf}\n");

    // sidecar caption
    std::ofstream cap("papers/figures/satellite_removal/satellite_removal_CAPTION.txt");
    cap << "Approved caption (satellite/tail-removal demonstration):\n\n";
    char b[900];
    snprintf(b,sizeof(b),
      "\"Tail suppression by the saturation-recovered CFD on identical events (DSB1, brightest-1000,\n"
      "production selection). (A) The per-event (DW-UP)/2 distribution -- the quantity whose width is\n"
      "the time resolution -- at 150 GeV under cfd05 (line) and srCFD (filled), unit-normalized,\n"
      "logarithmic scale; the common tail window |x-median| > %.0f ps (2.5 sigma of the srCFD core) is\n"
      "marked. (B) Width and tail metrics on the same events: the Gaussian-core widths differ by\n"
      "%.1f ps while the tail fraction drops from %.2f%% to %.2f%% (Gaussian expectation 1.24%%) --\n"
      "the improvement from srCFD is concentrated in the non-Gaussian tail component, supporting its\n"
      "interpretation as recovery from estimator bias on clipped pulses rather than a change in\n"
      "detector response. (C) At 75 GeV (full-fiducial sample; in the brightest-1000 the\n"
      "less-clipped subset is empty, N=2, since bright showers clip by construction): the cfd05\n"
      "tail excess concentrates in the clipped subset and largely vanishes in the less-clipped\n"
      "subset, while srCFD shows the OPPOSITE pattern -- on dim, unclipped pulses its low-gain\n"
      "amplitude anchor is noisy and its tail fraction rises. Each estimator is biased in the\n"
      "regime it was not designed for, which is the basis of the per-regime adopted-source rule\n"
      "(srCFD for bright/clipped, fixed-threshold for dim).\"\n\nFORBIDDEN in any version: 'improves intrinsic detector timing',\n"
      "'saturation fully corrected', 'all method dependence eliminated', universal CFD claims.\n",
      W*1000,m150[0].prod,m150[1].prod,m150[1].ftail,m150[0].ftail);
    cap << b; cap.close();
    printf("  wrote papers/figures/satellite_removal/satellite_removal_CAPTION.txt\n");
}

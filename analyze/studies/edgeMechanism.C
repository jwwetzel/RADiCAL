// ===========================================================================
// edgeMechanism.C  —  the *cause* behind the CFD-fraction timing trend.
//
//   Two square panels:
//   layer4_edge_shape.png  — mean normalized rising edge of a Down (SE-D) vs an
//       Up (SE-U) capillary at 100 GeV, aligned at the 20% crossing, with the
//       5% and 20% levels marked.  KEY POINT: the 20% level sits on a STEEP part
//       of the edge (slope annotated) and the Down edge is, if anything, steeper
//       there than the Up edge -- so the Down-capillary timing shoulder is NOT a
//       mean-slope / electronic-noise effect (that would make 20% *better*).
//
//   layer4_edge_jitter.png — the actual mechanism: the pulse-to-pulse leading-
//       edge SHAPE jitter sigma(t_f - t_5%) (MCP cancels; purely intra-pulse)
//       vs CFD fraction, for the Down vs Up capillary groups.  It rises the
//       higher up the edge one times, and is systematically larger for the Down
//       capillaries -- so timing low on the edge (CFD-5%) is the most reproducible
//       point.  Energy-independent (50-150 GeV) => an intrinsic pulse-shape
//       property, not shower physics.
//
// Run: ROOT_INCLUDE_PATH=Analysis root -l -b -q 'Analysis/edgeMechanism.C+'
// ===========================================================================
#include <cstdio>
#include <vector>
#include <algorithm>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TLine.h"
#include "TMarker.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TSystem.h"

#include "RADiCALStyle.h"

namespace {

const char* kSumDir = "output/Summary/";

double medianOf(std::vector<float> v){
    if (v.empty()) return 0.;
    std::nth_element(v.begin(), v.begin()+v.size()/2, v.end());
    return v[v.size()/2];
}
double winRMS(const std::vector<float>& v, double win, long& n){
    double med = medianOf(v);
    double s=0,s2=0; n=0;
    for (float x : v) if (std::fabs(x-med)<=win){ s+=x; s2+=x*x; ++n; }
    if (n<2) return -1.;
    double m=s/n, var=s2/n-m*m; return var>0.?std::sqrt(var):0.;
}

// crossing time where the normalized profile first reaches frac (linear interp)
double crossTime(TProfile* p, double pk, double frac){
    double tgt = frac*pk; int nb=p->GetNbinsX();
    for (int b=2;b<=nb;++b){ double v0=p->GetBinContent(b-1), v1=p->GetBinContent(b);
        if (v0<tgt && v1>=tgt){ double t0=p->GetBinCenter(b-1), t1=p->GetBinCenter(b);
            return t0 + (t1-t0)*(tgt-v0)/(v1-v0); } }
    return 0.;
}
double normSlopeAt(TProfile* p, double pk, double frac){
    double tc = crossTime(p, pk, frac);
    int b = p->FindBin(tc);
    if (b<2 || b>=p->GetNbinsX()) return 0.;
    double dt = p->GetBinCenter(b+1)-p->GetBinCenter(b-1);
    return (p->GetBinContent(b+1)-p->GetBinContent(b-1))/dt/pk;   // (1/ns)
}

// ── Panel 1: mean rising edge, Down vs Up ───────────────────────────────────
void EdgeShape()
{
    TFile f(Form("%saverage_waveforms.root", kSumDir));
    if (f.IsZombie()){ printf("[edge] no average_waveforms.root — skip shape\n"); return; }
    TProfile* pD = (TProfile*)f.Get("sum_100GeV_SE-D");
    TProfile* pU = (TProfile*)f.Get("sum_100GeV_SE-U");
    if (!pD || !pU){ printf("[edge] profiles missing — skip shape\n"); return; }

    double pkD=0, pkU=0; int nb=pD->GetNbinsX();
    for (int b=1;b<=nb;++b){ pkD=std::max(pkD,pD->GetBinContent(b)); pkU=std::max(pkU,pU->GetBinContent(b)); }
    double t20D = crossTime(pD,pkD,0.20), t20U = crossTime(pU,pkU,0.20);
    double slD = normSlopeAt(pD,pkD,0.20), slU = normSlopeAt(pU,pkU,0.20);

    // build normalized, 20%-aligned edge graphs
    TGraph* gD = new TGraph(); TGraph* gU = new TGraph();
    int iD=0, iU=0;
    for (int b=1;b<=nb;++b){
        double t=pD->GetBinCenter(b)-t20D, y=pD->GetBinContent(b)/pkD;
        if (t>-1.2 && t<1.8) gD->SetPoint(iD++, t, y);
    }
    for (int b=1;b<=nb;++b){
        double t=pU->GetBinCenter(b)-t20U, y=pU->GetBinContent(b)/pkU;
        if (t>-1.2 && t<1.8) gU->SetPoint(iU++, t, y);
    }

    TCanvas* c = NewSquareCanvas("layer4_edge_shape", 660);
    c->cd();
    gD->SetLineColor(kRRed);   gD->SetLineWidth(3);
    gU->SetLineColor(kRData);  gU->SetLineWidth(3);
    gD->Draw("AL");
    gD->GetXaxis()->SetTitle("t #minus t_{20%}  (ns)");
    gD->GetYaxis()->SetTitle("a(t) / a_{max}");
    gD->GetXaxis()->SetLimits(-1.2, 1.8);
    gD->GetYaxis()->SetRangeUser(0., 1.08);
    gD->GetXaxis()->SetTitleSize(0.046);
    gD->GetYaxis()->SetTitleSize(0.046);
    gU->Draw("L SAME");

    // 5% and 20% level lines
    for (double fr : {0.05, 0.20}){
        TLine* l = new TLine(-1.2, fr, 1.8, fr);
        l->SetLineColor(kGray+1); l->SetLineStyle(3); l->SetLineWidth(1); l->Draw("SAME");
    }
    { TLatex t; t.SetTextColor(kGray+2); t.SetTextSize(0.030);
      t.DrawLatex(-1.12, 0.075, "5%"); t.DrawLatex(-1.12, 0.225, "20%"); }

    // markers at the 20% crossings (both at t=0 by construction)
    TMarker* mD = new TMarker(0., 0.20, 20); mD->SetMarkerColor(kRRed);  mD->SetMarkerSize(1.4); mD->Draw();
    TMarker* mU = new TMarker(0., 0.20, 20); mU->SetMarkerColor(kRData); mU->SetMarkerSize(1.4); mU->Draw();

    TLegend* leg = MakeCornerLegend(2, "br", 0.038);
    leg->AddEntry(gD, "SE-D (Down)", "l");
    leg->AddEntry(gU, "SE-U (Up)",   "l");
    leg->Draw();

    { TLatex t; t.SetNDC(); t.SetTextSize(0.030); t.SetTextColor(kGray+3);
      t.DrawLatex(0.16, 0.78, Form("norm. slope @20%%:"));
      t.SetTextColor(kRRed);  t.DrawLatex(0.16, 0.735, Form("Down %.2f /ns", slD));
      t.SetTextColor(kRData); t.DrawLatex(0.16, 0.690, Form("Up   %.2f /ns", slU));
      t.SetTextColor(kGray+3);t.DrawLatex(0.16, 0.63, "20% is on a STEEP part");
      t.DrawLatex(0.16, 0.59, "#Rightarrow not a slope/noise effect"); }

    DrawPadTitle("Mean rising edge (100 GeV)", 0.052f);
    c->Print(Form("%slayer4_edge_shape.png", kSumDir));
    printf("[edge] wrote layer4_edge_shape.png  (slope@20%% D=%.2f U=%.2f /ns)\n", slD, slU);
}

// ── Panel 2: edge-shape jitter sigma(t_f - t_5%) vs fraction ────────────────
void EdgeJitter()
{
    const char* br[5]  = {"hg_cfd05","hg_cfd10","hg_cfd","hg_cfd30","hg_cfd50"};
    const double frac[5]= {5,10,20,30,50};
    const int down[4]={0,1,2,3}, up[4]={4,5,6,7};

    // representative energy 100 GeV (verified energy-independent 50-150)
    TFile f("output/100GeV/ntuple.root");
    if (f.IsZombie()){ printf("[edge] no 100GeV ntuple — skip jitter\n"); return; }
    TTree* t=(TTree*)f.Get("rad");
    Float_t cfd[5][8], hg_peak[8]; Bool_t in_fid;
    // distinct arrays; note cfd[0] (hg_cfd05) is the reference
    for (int j=0;j<5;++j) t->SetBranchAddress(br[j], cfd[j]);
    t->SetBranchAddress("hg_peak", hg_peak);
    t->SetBranchAddress("in_fiducial", &in_fid);

    TGraphErrors* g[2] = { new TGraphErrors(), new TGraphErrors() };
    for (int grp=0; grp<2; ++grp){
        const int* ch = grp? up : down;
        int pt=0;
        // anchor at (5%, 0) by definition
        g[grp]->SetPoint(pt, 5., 0.); g[grp]->SetPointError(pt, 0., 0.); ++pt;
        for (int j=1;j<5;++j){           // f = 10,20,30,50
            std::vector<float> d;
            for (Long64_t e=0;e<t->GetEntries();++e){ t->GetEntry(e); if(!in_fid) continue;
                for (int lc=0; lc<4; ++lc){ int c=ch[lc]; if (hg_peak[c]<=20.f) continue;
                    float a=cfd[j][c], b=cfd[0][c];
                    if (a>-60.f&&a<60.f&&b>-60.f&&b<60.f) d.push_back(a-b); } }
            long n; double s = winRMS(d, 1.0, n);
            if (s>=0.){ g[grp]->SetPoint(pt, frac[j], s*1000.);
                        g[grp]->SetPointError(pt, 0., (s*1000.)/std::sqrt(2.0*n)); ++pt; }
        }
    }
    f.Close();

    TCanvas* c = NewSquareCanvas("layer4_edge_jitter", 660);
    c->cd();
    g[0]->SetMarkerStyle(20); g[0]->SetMarkerSize(1.3); g[0]->SetMarkerColor(kRRed);
    g[0]->SetLineColor(kRRed);  g[0]->SetLineWidth(3);
    g[1]->SetMarkerStyle(21); g[1]->SetMarkerSize(1.3); g[1]->SetMarkerColor(kRData);
    g[1]->SetLineColor(kRData); g[1]->SetLineWidth(3);

    g[0]->Draw("APL");
    g[0]->GetXaxis()->SetTitle("CFD fraction (%)");
    g[0]->GetYaxis()->SetTitle("edge-shape jitter  #sigma(t_{f}#minus t_{5%})  (ps)");
    g[0]->GetXaxis()->SetLimits(0., 54.);
    g[0]->GetYaxis()->SetRangeUser(0., 360.);
    g[0]->GetXaxis()->SetTitleSize(0.044);
    g[0]->GetYaxis()->SetTitleSize(0.042);
    g[1]->Draw("PL SAME");

    // mark 20% (where CFD-20% times)
    TLine* l20 = new TLine(20., 0., 20., 335.);
    l20->SetLineColor(kGray+2); l20->SetLineStyle(2); l20->SetLineWidth(2); l20->Draw("SAME");
    { TLatex a; a.SetTextColor(kGray+3); a.SetTextSize(0.032); a.SetTextAngle(90);
      a.DrawLatex(22.0, 120., "CFD-20%"); }

    TLegend* leg = MakeCornerLegend(2, "tl", 0.040);
    leg->AddEntry(g[0], "Down capillaries", "lp");
    leg->AddEntry(g[1], "Up capillaries",   "lp");
    leg->Draw();

    { TLatex t; t.SetNDC(); t.SetTextSize(0.029); t.SetTextColor(kGray+3);
      t.DrawLatex(0.46, 0.30, "Edge shape jitters more");
      t.DrawLatex(0.46, 0.255, "the higher you time, and");
      t.DrawLatex(0.46, 0.21, "more for the Down edge.");
      t.DrawLatex(0.46, 0.155, "The 5% edge = most reproducible");
      t.DrawLatex(0.46, 0.11,  "(srCFD's operating point)."); }

    DrawPadTitle("Leading-edge shape jitter (100 GeV)", 0.050f);
    c->Print(Form("%slayer4_edge_jitter.png", kSumDir));
    printf("[edge] wrote layer4_edge_jitter.png\n");
}

} // namespace

void edgeMechanism()
{
    ApplyRADiCALStyle();
    gSystem->mkdir(kSumDir, kTRUE);
    EdgeShape();
    EdgeJitter();
    printf("\n[edgeMechanism] done — heroes in %s\n", kSumDir);
}

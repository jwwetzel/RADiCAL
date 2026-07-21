// ===========================================================================
// elbowFractionTrend.C  —  publication-quality version of the "elbow" finding.
//
//   The per-channel timing resolution of the Down vs Up capillary GROUPS as a
//   function of CFD fraction, at all six beam energies, with error bars.
//
//   This is the quantitative statement behind the per-channel page-3 diagnostic:
//   on the DOWN capillaries sigma_t rises steeply with CFD fraction (the leading-
//   edge SHAPE jitters more high on the edge — not a mean-slope effect; see
//   edgeMechanism.C), while the UP capillaries depend far less on fraction.
//   The 5% edge is the optimum among plain CFD fractions; the production
//   headline estimator (srCFD) times exactly this low edge, with LG-guided
//   recovery of clipped pulses (see hg_lgcfd in the reduction).
//
//   Estimator (stated on the figure caption in the report):
//     For each (energy, group, fraction) the four capillaries of the group are
//     each shifted to their own median Delta-t (MCP-referenced; fiducial events,
//     hg_peak > 20 mV, |Delta-t| < 60 ns) and POOLED.  sigma_t is the RMS within
//     +/-1.0 ns of the pooled median (truncated to exclude failed-crossing
//     outliers); the statistical error is sigma/sqrt(2N).
//
// Output: Summary/layer4_cfd_fraction_down.png and _up.png (square heroes).
// Run: ROOT_INCLUDE_PATH=Analysis root -l -b -q 'Analysis/elbowFractionTrend.C+'
// ===========================================================================
#include <cstdio>
#include <vector>
#include <algorithm>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TSystem.h"

#include "RADiCALStyle.h"

namespace {

const char* kSumDir = "output/Summary/";

const char* kEneLabel[6] = {"25GeV","50GeV","75GeV","100GeV","125GeV","150GeV"};
const double kEneGeV[6]   = {25,50,75,100,125,150};
const char* kCapName[8]   = {"NW-D","NE-D","SE-D","SW-D","NW-U","NE-U","SE-U","SW-U"};

// CFD fractions and their ntuple branches (ascending fraction)
const int    kNF = 6;
const double kFrac[kNF]      = {3, 5, 10, 20, 30, 50};
const char*  kFracBr[kNF]    = {"hg_cfd03","hg_cfd05","hg_cfd10","hg_cfd","hg_cfd30","hg_cfd50"};

double median(std::vector<float> v){
    if (v.empty()) return 0.;
    std::nth_element(v.begin(), v.begin()+v.size()/2, v.end());
    return v[v.size()/2];
}
// windowed RMS within +/- win of center c; returns sigma and fills nUsed
double winRMS(const std::vector<float>& v, double c, double win, long& nUsed){
    double s=0,s2=0; long k=0;
    for (float x : v) if (std::fabs(x-c) <= win) { s+=x; s2+=x*x; ++k; }
    nUsed = k;
    if (k < 2) return -1.;
    double m=s/k, var=s2/k - m*m;
    return var > 0. ? std::sqrt(var) : 0.;
}

// Build one square panel for a capillary group (chans = 4 channel indices).
void DrawGroup(const char* pngName, const char* title,
               const int chans[4], double yMax)
{
    // sig[ie][jf], err[ie][jf]  in ps
    double sig[6][kNF], err[6][kNF];
    for (int ie=0; ie<6; ++ie) for (int jf=0; jf<kNF; ++jf){ sig[ie][jf]=-1; err[ie][jf]=0; }

    for (int ie=0; ie<6; ++ie){
        TString fn = TString(Form("output/%s/ntuple.root", kEneLabel[ie]));
        TFile f(fn);
        if (f.IsZombie()) { printf("[trend] missing %s — skip\n", fn.Data()); continue; }
        TTree* t = (TTree*)f.Get("rad");
        if (!t) continue;

        Float_t hg_peak[8]; Bool_t in_fid;
        Float_t cfd[kNF][8];
        t->SetBranchAddress("hg_peak", hg_peak);
        t->SetBranchAddress("in_fiducial", &in_fid);
        for (int jf=0; jf<kNF; ++jf) t->SetBranchAddress(kFracBr[jf], cfd[jf]);

        // raw[localChan][frac] -> values
        std::vector<float> raw[4][kNF];
        Long64_t N = t->GetEntries();
        for (Long64_t e=0; e<N; ++e){
            t->GetEntry(e);
            if (!in_fid) continue;
            for (int lc=0; lc<4; ++lc){
                int ch = chans[lc];
                if (hg_peak[ch] <= 20.f) continue;
                for (int jf=0; jf<kNF; ++jf){
                    float v = cfd[jf][ch];
                    if (v > -60.f && v < 60.f) raw[lc][jf].push_back(v);
                }
            }
        }
        // median-subtract each channel, pool the 4 channels per fraction
        for (int jf=0; jf<kNF; ++jf){
            std::vector<float> pooled;
            for (int lc=0; lc<4; ++lc){
                double med = median(raw[lc][jf]);
                for (float v : raw[lc][jf]) pooled.push_back(v - med);
            }
            long n=0;
            double med0 = median(pooled);
            double s = winRMS(pooled, med0, 1.0, n);
            if (s >= 0. && n > 10){ sig[ie][jf] = s*1000.; err[ie][jf] = (s*1000.)/std::sqrt(2.0*n); }
        }
        f.Close();
    }

    // ---- draw ----
    TCanvas* c = NewSquareCanvas(pngName, 660);
    c->cd();

    TGraphErrors* g[6] = {};
    double drawnMax = 0.;
    for (int ie=0; ie<6; ++ie){
        g[ie] = new TGraphErrors();
        g[ie]->SetName(Form("g_%s_%s", title, kEneLabel[ie]));
        int p=0;
        double dodge = (ie - 2.5) * 0.45;   // small x-dodge so error bars don't overlap
        for (int jf=0; jf<kNF; ++jf){
            if (sig[ie][jf] < 0.) continue;
            g[ie]->SetPoint(p, kFrac[jf] + dodge, sig[ie][jf]);
            g[ie]->SetPointError(p, 0., err[ie][jf]);
            if (sig[ie][jf] + err[ie][jf] > drawnMax) drawnMax = sig[ie][jf] + err[ie][jf];
            ++p;
        }
        g[ie]->SetMarkerStyle(20);
        g[ie]->SetMarkerSize(1.15);
        g[ie]->SetMarkerColor(kREnergyCols[ie]);
        g[ie]->SetLineColor(kREnergyCols[ie]);
        g[ie]->SetLineWidth(2);
    }
    double hi = (yMax > 0.) ? yMax : drawnMax * 1.12;

    // frame via the first graph
    bool first = true;
    for (int ie=0; ie<6; ++ie){
        if (g[ie]->GetN() == 0) continue;
        if (first){
            g[ie]->Draw("APL");
            g[ie]->GetXaxis()->SetTitle("CFD fraction (%)");
            g[ie]->GetYaxis()->SetTitle("group #sigma_{t}  (ps)");
            g[ie]->GetXaxis()->SetLimits(0., 54.);
            g[ie]->GetYaxis()->SetRangeUser(0., hi);
            g[ie]->GetXaxis()->SetTitleSize(0.046);
            g[ie]->GetYaxis()->SetTitleSize(0.046);
            first = false;
        } else {
            g[ie]->Draw("PL SAME");
        }
    }

    // adopted-5% marker
    TLine* adopt = new TLine(5., 0., 5., hi*0.93);
    adopt->SetLineColor(kGray+2); adopt->SetLineStyle(2); adopt->SetLineWidth(2);
    adopt->Draw("SAME");
    { TLatex a; a.SetTextColor(kGray+3); a.SetTextSize(0.034); a.SetTextAngle(90);
      a.DrawLatex(7.2, hi*0.40, "5% edge #minus srCFD's operating point"); }

    // legend (energies)
    TLegend* leg = MakeCornerLegend(6, "tl", 0.034);
    for (int ie=0; ie<6; ++ie)
        if (g[ie]->GetN()) leg->AddEntry(g[ie], Form("%.0f GeV", kEneGeV[ie]), "lp");
    leg->Draw();

    DrawPadTitle(title, 0.052f);
    c->Print(Form("%s%s.png", kSumDir, pngName));
    printf("[trend] wrote %s.png\n", pngName);
}

} // namespace

void elbowFractionTrend()
{
    ApplyRADiCALStyle();
    gSystem->mkdir(kSumDir, kTRUE);

    const int down[4] = {0,1,2,3};   // NW-D NE-D SE-D SW-D
    const int up[4]   = {4,5,6,7};   // NW-U NE-U SE-U SW-U

    // Shared y-axis so the Down rise and Up flatness are directly comparable.
    const double yShared = 520.;
    DrawGroup("layer4_cfd_fraction_down", "Down capillaries", down, yShared);
    DrawGroup("layer4_cfd_fraction_up",   "Up capillaries",   up,   yShared);

    printf("\n[elbowFractionTrend] done — heroes in %s\n", kSumDir);
}

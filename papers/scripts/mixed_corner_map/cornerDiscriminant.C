// ============================================================================
// cornerDiscriminant.C — GATE 1 (data-only half): brightness-INDEPENDENT pulse-
// shape classification of the MIXED corners against labeled references.
// Protocol pre-registered in AUDIT.md (this directory). Discriminants:
//   D1 = lg_charge/lg_peak  (LG width; never clips; PRIMARY)
//   D2 = hg_cfd50 - hg_cfd05 (5%->50% rise time; UNCLIPPED events only)
//   D3 = hg_charge/hg_peak  (HG width; UNCLIPPED only)
// References: pure DSB1 build (8 labeled DSB1 ends) and pure LUAG build (8
// labeled LuAG ends), same reduction batch. Classification: per MIXED end,
// compare the median discriminant to the two reference medians (distance in
// units of reference robust width); both ends of a capillary must agree.
//   source setup.sh
//   root -l -b -q 'papers/scripts/mixed_corner_map/cornerDiscriminant.C+'
// Output: papers/figures/mixed_corner_map/corner_discriminant.png + stdout table.
// ============================================================================
#include "RadView.h"
#include "DataPaths.h"
#include "PlotUtils.h"
#include "SelectionCuts.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>
using namespace rad;

static const char* ENM[8]={"NW-D","NE-D","SE-D","SW-D","NW-U","NE-U","SE-U","SW-U"};
// inferred map to TEST (1 = DSB1-like expected, 0 = LuAG-like expected)
static const int EXPECT_DSB1[8]={0,1,0,1, 0,1,0,1};   // NE,SW DSB1 ; NW,SE LuAG

struct Med { double med=0, w=0; long n=0; };   // median + robust width (1.4826*MAD)
static Med medw(std::vector<float>& v){
    Med m; if(v.size()<500) return m;
    std::sort(v.begin(),v.end()); m.n=v.size(); m.med=v[v.size()/2];
    std::vector<float> ad(v.size()); for(size_t i=0;i<v.size();++i) ad[i]=std::fabs(v[i]-m.med);
    std::nth_element(ad.begin(),ad.begin()+ad.size()/2,ad.end()); m.w=1.4826*ad[ad.size()/2];
    if(m.w<1e-6) m.w=1e-6; return m;
}

// collect D1 (all events), D2+D3 (unclipped) per end for one build over its energies
struct EndData { std::vector<float> d1,d2,d3; };
static void collect(const char* build, EndData ed[8]){
    BuildConfig cfg=BuildConfig::Load(radConfig(build).Data());
    const double Es[6]={25,50,75,100,125,150};
    for(double e:Es){ TString p=radReduced(build,e); if(gSystem->AccessPathName(p.Data()))continue;
        TFile* f=TFile::Open(p.Data()); if(!f||f->IsZombie()){if(f)f->Close();continue;}
        TTree* t=(TTree*)f->Get("rad"); if(!t){f->Close();continue;}
        // direct branch access (shape ratios need raw arrays incl. non-timing E-type ends)
        Float_t lgq[8],lgp[8],hgq[8],hgp[8],c05[8],c50[8],mcp; Bool_t sat[8]; Bool_t wc;
        t->SetBranchStatus("*",0);
        for(const char* b:{"lg_charge","lg_peak","hg_charge","hg_peak","hg_cfd05","hg_cfd50","hg_saturated","mcp1_peak","wc_ok"}) t->SetBranchStatus(b,1);
        t->SetBranchAddress("lg_charge",lgq); t->SetBranchAddress("lg_peak",lgp);
        t->SetBranchAddress("hg_charge",hgq); t->SetBranchAddress("hg_peak",hgp);
        t->SetBranchAddress("hg_cfd05",c05);  t->SetBranchAddress("hg_cfd50",c50);
        t->SetBranchAddress("hg_saturated",sat); t->SetBranchAddress("mcp1_peak",&mcp); t->SetBranchAddress("wc_ok",&wc);
        Long64_t N=t->GetEntries(); Long64_t cap=600000; if(N>cap)N=cap;   // plenty for medians
        for(Long64_t i=0;i<N;++i){ t->GetEntry(i);
            if(!wc||mcp<kMCP1_minPeak||mcp>kMCP1_maxPeak) continue;
            for(int c=0;c<8;++c){
                if(lgp[c]>10 && lgq[c]>0) ed[c].d1.push_back(lgq[c]/lgp[c]);
                if(!sat[c] && hgp[c]>60 && hgp[c]<700){            // clean unclipped HG pulse
                    if(c05[c]>-1e5f&&c50[c]>-1e5f) ed[c].d2.push_back(c50[c]-c05[c]);
                    if(hgq[c]>0) ed[c].d3.push_back(hgq[c]/hgp[c]); } }
        }
        f->Close();
    }
}

void cornerDiscriminant(){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    gSystem->mkdir("papers/figures/mixed_corner_map",kTRUE);
    EndData dsb[8], lua[8], mix[8];
    printf("\n===== GATE 1 (data half): MIXED corner pulse-shape classification =====\n");
    printf("collecting references + MIXED ...\n");
    collect("DSB1",dsb); collect("LUAG",lua); collect("MIXED",mix);

    // pooled reference templates (all 8 ends of each pure build)
    auto pool=[&](EndData* ed,int which)->std::vector<float>{ std::vector<float> o;
        for(int c=0;c<8;++c){ auto&v=(which==1?ed[c].d1:which==2?ed[c].d2:ed[c].d3); o.insert(o.end(),v.begin(),v.end()); } return o; };
    std::vector<float> rd1=pool(dsb,1), rl1=pool(lua,1), rd2=pool(dsb,2), rl2=pool(lua,2), rd3=pool(dsb,3), rl3=pool(lua,3);
    Med Rd1=medw(rd1), Rl1=medw(rl1), Rd2=medw(rd2), Rl2=medw(rl2), Rd3=medw(rd3), Rl3=medw(rl3);
    printf("\n  reference templates (median +- robust width):\n");
    printf("    D1 lgQ/lgP : DSB1 %.2f+-%.2f   LuAG %.2f+-%.2f\n",Rd1.med,Rd1.w,Rl1.med,Rl1.w);
    printf("    D2 cfd50-05: DSB1 %.3f+-%.3f  LuAG %.3f+-%.3f (ns)\n",Rd2.med,Rd2.w,Rl2.med,Rl2.w);
    printf("    D3 hgQ/hgP : DSB1 %.2f+-%.2f   LuAG %.2f+-%.2f\n",Rd3.med,Rd3.w,Rl3.med,Rl3.w);

    // classify each MIXED end: signed score s = (|d-LuAG| - |d-DSB1|) / (w_D+w_L)/2  (>0 => DSB1-like)
    printf("\n  per-end classification (score>0 => DSB1-like; expected map in [..]):\n");
    printf("  %-6s %10s %10s %10s   verdict  vs expected\n","end","D1 score","D2 score","D3 score");
    int agree=0, total=0; double s1all[8];
    for(int c=0;c<8;++c){
        Med m1=medw(mix[c].d1), m2=medw(mix[c].d2), m3=medw(mix[c].d3);
        auto score=[&](Med& m,Med& RD,Med& RL)->double{ if(m.n<500) return 0.0/0.0;
            return (std::fabs(m.med-RL.med)-std::fabs(m.med-RD.med))/(0.5*(RD.w+RL.w)); };
        double s1=score(m1,Rd1,Rl1), s2=score(m2,Rd2,Rl2), s3=score(m3,Rd3,Rl3); s1all[c]=s1;
        int votes=0,nv=0;
        for(double s:{s1,s2,s3}) if(s==s){ ++nv; if(s>0)++votes; }
        bool isDSB1 = (nv>0)&&(votes*2>nv);
        bool exp_=EXPECT_DSB1[c]; bool ok=(isDSB1==exp_);
        ++total; if(ok)++agree;
        printf("  %-6s %10.1f %10.1f %10.1f   %-9s [%s]  %s\n",ENM[c],s1,s2,s3,
               isDSB1?"DSB1":"LuAG", exp_?"DSB1":"LuAG", ok?"OK":"** CONTRADICTS **");
    }
    printf("\n  capillary-pair consistency: ");
    bool pairOK=true; for(int k=0;k<4;++k){ bool a=s1all[k]>0,b=s1all[k+4]>0; if(a!=b){pairOK=false;printf("%s/%s DISAGREE  ",ENM[k],ENM[k+4]);} }
    printf(pairOK?"all 4 capillaries: D and U ends agree\n":"\n");
    printf("  VERDICT: %d/%d ends match the inferred map (NE+SW=DSB1, NW+SE=LuAG)\n",agree,total);

    // figure: D1 reference distributions + the 8 MIXED end medians
    TCanvas* c=new TCanvas("cd","",1000,700); c->SetLeftMargin(0.12);c->SetRightMargin(0.04);c->SetTopMargin(0.08);c->SetBottomMargin(0.13);
    double xmax=std::max(Rl1.med+5*Rl1.w, Rd1.med+5*Rd1.w)*1.15;
    TH1F* hD=new TH1F("hD",";lg_charge / lg_peak  (effective LG pulse width);normalized",120,0,xmax);
    TH1F* hL=new TH1F("hL","",120,0,xmax); hD->SetDirectory(nullptr); hL->SetDirectory(nullptr);
    for(float x:rd1)hD->Fill(x); for(float x:rl1)hL->Fill(x);
    if(hD->Integral()>0)hD->Scale(1.0/hD->Integral()); if(hL->Integral()>0)hL->Scale(1.0/hL->Integral());
    hD->SetLineColor(kRed+1); hD->SetFillColorAlpha(kRed+1,0.25); hD->SetLineWidth(2);
    hL->SetLineColor(kGreen+3); hL->SetFillColorAlpha(kGreen+3,0.25); hL->SetLineWidth(2);
    double ym=1.35*std::max(hD->GetMaximum(),hL->GetMaximum()); hD->SetMaximum(ym);
    hD->Draw("HIST"); hL->Draw("HIST SAME");
    TLegend* lg=new TLegend(0.55,0.70,0.95,0.90); lg->SetBorderSize(0);lg->SetFillStyle(0);lg->SetTextSize(0.034);
    lg->AddEntry(hD,"reference: pure DSB1 build (8 ends)","f");
    lg->AddEntry(hL,"reference: pure LuAG build (8 ends)","f");
    lg->Draw();
    TLatex tx; tx.SetTextSize(0.030);
    for(int k=0;k<8;++k){ Med m=medw(mix[k].d1); if(m.n<500) continue;
        bool isD=s1all[k]>0;
        TLine* l=new TLine(m.med,0,m.med,0.55*ym); l->SetLineColor(isD?kRed+2:kGreen+3); l->SetLineWidth(2); l->SetLineStyle(2); l->Draw();
        tx.SetTextColor(isD?kRed+2:kGreen+3); tx.SetTextAngle(90);
        tx.DrawLatex(m.med+xmax*0.004,0.57*ym,ENM[k]); }
    tx.SetTextAngle(0); tx.SetNDC(); tx.SetTextColor(kGray+3); tx.SetTextSize(0.030);
    tx.DrawLatex(0.14,0.86,"dashed = MIXED end medians");
    DrawSuperTitle("MIXED corner identification by pulse SHAPE (brightness-independent): LG width D1 vs labeled references",0.020f);
    c->Print("papers/figures/mixed_corner_map/corner_discriminant.png");
    printf("  wrote papers/figures/mixed_corner_map/corner_discriminant.png\n");
}

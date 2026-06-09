// ============================================================================
// methodCompare.C — motivation for the choice of timing estimator, per build.
// ----------------------------------------------------------------------------
// On IDENTICAL events (same brightest-1000 selection, same (DW-UP)/2 estimator,
// same a/sqrt(E)(+)b photostat fit -- only the per-channel TIME SOURCE differs),
// compares the three candidate timing methods:
//   * cfd05   — CFD at 5% of the CLIPPED measured peak  (the PUBLISHED method)
//   * led     — leading edge at a FIXED absolute threshold (robust on dim/slow pulses)
//   * hg_lgcfd— CFD on the LG-recovered TRUE peak, on the steep edge below the clip
// The method that wins is set by the LIGHT REGIME: a bright, heavily-clipped build
// (DSB1) must recover its edge (lgcfd); a dim build (LuAG) cannot use the published
// 5% foot (it sits in noise) and wants led. So the estimator choice is itself a
// fingerprint of the light-yield thesis. cfd05 is noisy / non-photostatistic on the
// dim builds (its a/sqrt(E)+b fit is poor) so its fit line is suppressed there.
//
//   source setup.sh
//   root -l -b -q -e '.L analyze/studies/methodCompare.C+' -e 'methodSurvey()'  // 3-panel + table
//   root -l -b -q 'analyze/studies/methodCompare.C+("DSB1")'                    // one build
// Output: figures/<year>/narrative/method_compare_<BUILD>.png, method_survey.png,
//         papers/timing/tab_methods.tex
// ============================================================================
#include "RadView.h"
#include "RadTiming.h"
#include "DataPaths.h"
#include "FigPaths.h"
#include "PlotUtils.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include <vector>
#include <string>
#include <cmath>
#include <cstdio>
#include <fstream>
using namespace rad;

static const int    SRCS[3] = { RadView::kCFD05, RadView::kLED, RadView::kLGCFD };
static const char*  SLAB[3] = { "cfd05", "led", "hg_lgcfd" };
static int          adoptedSrc(const char* b){ std::string s=b; return (s=="LUAG"||s=="TENERGY")?RadView::kLED:RadView::kLGCFD; }
static const char*  titleFor(const char* b){ std::string s=b;
    if(s=="DSB1")  return "DSB1 #minus high light, heavily clipped";
    if(s=="LUAG")  return "LuAG #minus low light, little clipping";
    if(s=="MIXED") return "MIXED #minus bright DSB1 + dim LuAG corners";
    if(s=="TENERGY") return "TENERGY #minus high light, segmented";
    return b; }

struct Ladder { std::vector<double> E,S,Se,ze; double a=0,b=0,be=0,chi=0; double s150=-1; };

static Ladder ladder(const char* build, BuildConfig& cfg, int src){
    Ladder L; const double Es[]={25,50,75,100,125,150};
    for(double e:Es){ TString p=radReduced(build,e); if(gSystem->AccessPathName(p.Data())) continue;
        TFile* fp=TFile::Open(p.Data()); if(!fp||fp->IsZombie()){ if(fp)fp->Close(); continue; }
        TTree* t=(TTree*)fp->Get("rad"); if(!t){fp->Close();continue;}
        RadView v; v.attach(t,&cfg);
        TimingResult r=timingBrightestK(v,e,src,1000);
        if(r.sigma_ps>0){ L.E.push_back(e); L.S.push_back(r.sigma_ps); L.Se.push_back(r.sigma_ps/std::sqrt(2000.0)); L.ze.push_back(0);
                          if(e>149) L.s150=r.sigma_ps; }
        fp->Close(); }
    if(L.E.size()>=3){ TGraphErrors g(L.E.size(),&L.E[0],&L.S[0],&L.ze[0],&L.Se[0]);
        TF1 f("f","sqrt([0]*[0]/x+[1]*[1])",20,160); f.SetParameters(300,18); g.Fit(&f,"QN");
        L.a=std::fabs(f.GetParameter(0)); L.b=std::fabs(f.GetParameter(1));
        L.chi=f.GetChisquare()/std::max(1,f.GetNDF()); L.be=f.GetParError(1); if(L.chi>1) L.be*=std::sqrt(L.chi); }
    return L;
}

static void computeLadders(const char* build, Ladder L[3]){
    BuildConfig cfg=BuildConfig::Load(radConfig(build).Data());
    for(int i=0;i<3;++i) L[i]=ladder(build,cfg,SRCS[i]);
}

// draw one build's 3-method comparison on the current pad, with shared [lo,hi].
// lo>=hi => auto per-build range. Off-scale points (>hi) are annotated, not clipped silently.
// A photostat fit line is drawn ONLY when the fit is good (chi2/ndf < 8); cfd05 on the dim
// builds is non-photostatistic so it shows as points (no misleading curve).
static void drawPanel(const char* build, Ladder L[3], double lo, double hi){
    int adopt=adoptedSrc(build);
    int col[3]={kRed+1,kGreen+3,kAzure+2}; int mk[3]={24,25,26};
    if(lo>=hi){ double ymn=1e9,ymx=0;
        for(int i=0;i<3;++i) for(double s:L[i].S){ ymn=std::min(ymn,s); if(s<=80) ymx=std::max(ymx,s); }
        if(ymn>ymx){ymn=20;ymx=60;} double r=ymx-ymn; if(r<5)r=5; lo=std::max(0.0,ymn-0.10*r); hi=ymx+0.12*r; }
    TH1F* fr=gPad->DrawFrame(0,lo,165,hi);
    fr->SetTitle(Form("%s;beam energy E (GeV);brightest-1000  (DW#minusUP)/2  #sigma_{t} (ps)",titleFor(build)));
    fr->GetYaxis()->SetTitleSize(0.045); fr->GetYaxis()->SetTitleOffset(1.2); fr->GetXaxis()->SetTitleSize(0.045);
    TLegend* lg=new TLegend(0.54,0.66,0.95,0.89); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.037);
    for(int i=0;i<3;++i){ if(L[i].E.empty())continue; bool A=(SRCS[i]==adopt);
        TGraphErrors* g=new TGraphErrors(L[i].E.size(),&L[i].E[0],&L[i].S[0],&L[i].ze[0],&L[i].Se[0]);
        g->SetMarkerStyle(A?(20+i):mk[i]); g->SetMarkerColor(col[i]); g->SetLineColor(col[i]); g->SetMarkerSize(A?1.7:1.2);
        if(L[i].a>0 && L[i].chi<8.0){ TF1* f=new TF1(Form("f%s%d",build,i),"sqrt([0]*[0]/x+[1]*[1])",20,160);
            f->SetParameters(L[i].a,L[i].b); f->SetLineColor(col[i]); f->SetLineWidth(A?3:2); f->SetLineStyle(A?1:2); f->Draw("SAME"); }
        g->Draw("P SAME");
        const char* leglab=SLAB[i];
        lg->AddEntry(g,Form("%s%s",leglab,A?"  (adopted)":""),"lp");
        for(size_t j=0;j<L[i].S.size();++j) if(L[i].S[j]>hi){ TLatex tx; tx.SetTextSize(0.032); tx.SetTextColor(col[i]);
            tx.DrawLatex(L[i].E[j]-7, hi-0.07*(hi-lo), Form("#uparrow%.0f",L[i].S[j])); } }
    lg->Draw();
}

void methodCompare(const char* build="DSB1"){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    Ladder L[3]; computeLadders(build,L);
    int adopt=adoptedSrc(build);
    printf("\n=== %s: cfd05 (published) vs led vs hg_lgcfd, brightest-1000, IDENTICAL events ===\n",build);
    printf("   E      cfd05     led    lgcfd\n");
    const double Es[]={25,50,75,100,125,150};
    for(double e:Es){ bool any=false; for(int i=0;i<3;++i) for(double x:L[i].E) if(std::fabs(x-e)<1) any=true; if(!any)continue;
        printf("  %3.0f  ",e); for(int i=0;i<3;++i){ double s=-1; for(size_t j=0;j<L[i].E.size();++j) if(std::fabs(L[i].E[j]-e)<1) s=L[i].S[j];
            printf("%7.1f",s);} printf("\n"); }
    for(int i=0;i<3;++i) printf("  %-9s a=%4.0f b=%4.1f+-%.1f chi2/ndf=%.1f  sigma150=%.1f%s\n",
        SLAB[i],L[i].a,L[i].b,L[i].be,L[i].chi,L[i].s150, SRCS[i]==adopt?"  <- ADOPTED":"");
    TCanvas* c=new TCanvas("mc","",900,650); c->SetLeftMargin(0.12); c->SetRightMargin(0.04);
    c->SetTopMargin(0.07); c->SetBottomMargin(0.13); c->SetGridy();
    drawPanel(build,L,1,-1);   // auto per-build range
    c->Print(radFigP(Form("figures/narrative/method_compare_%s.png",build)));
    printf("  wrote figures/narrative/method_compare_%s.png\n",build);
}

void methodSurvey(){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    const char* builds[3]={"DSB1","LUAG","MIXED"};
    Ladder LL[3][3];
    for(int b=0;b<3;++b) computeLadders(builds[b], LL[b]);
    // shared y-range: ignore >80 ps outliers (the worst cfd05 points) for the upper bound; annotate those.
    double gmin=1e9, gmax=0;
    for(int b=0;b<3;++b) for(int i=0;i<3;++i) for(double s:LL[b][i].S){ gmin=std::min(gmin,s); if(s<=80) gmax=std::max(gmax,s); }
    double lo=std::max(0.0,gmin-5.0), hi=gmax+0.10*(gmax-gmin)+4.0;

    TCanvas* c=new TCanvas("ms","",1500,560); c->Divide(3,1,0.004,0.004);
    for(int b=0;b<3;++b){ c->cd(b+1); gPad->SetLeftMargin(0.13); gPad->SetRightMargin(0.03);
        gPad->SetTopMargin(0.09); gPad->SetBottomMargin(0.14); gPad->SetGridy();
        drawPanel(builds[b], LL[b], lo, hi); }
    c->cd(0); TLatex t; t.SetNDC(); t.SetTextFont(62); t.SetTextSize(0.029);
    t.DrawLatex(0.045,0.96,"Timing-method survey: cfd05 (published) vs led vs hg_lgcfd on identical brightest-1000 events #minus the estimator that wins tracks the light regime");
    c->Print(radFigP("figures/narrative/method_survey.png"));
    printf("  wrote figures/narrative/method_survey.png\n");

    for(int b=0;b<3;++b){ TCanvas* cc=new TCanvas(Form("mc%d",b),"",900,650);
        cc->SetLeftMargin(0.12); cc->SetRightMargin(0.04); cc->SetTopMargin(0.07); cc->SetBottomMargin(0.13); cc->SetGridy();
        drawPanel(builds[b], LL[b], 1, -1);
        cc->Print(radFigP(Form("figures/narrative/method_compare_%s.png",builds[b]))); }

    // LaTeX summary table
    std::ofstream o("papers/timing/tab_methods.tex");
    o<<"% auto-generated by analyze/studies/methodCompare.C (methodSurvey)\n";
    o<<"\\begin{table}[t]\n\\centering\n";
    o<<"\\caption{Brightest-1000 $(\\mathrm{DW}-\\mathrm{UP})/2$ time resolution under the three candidate\n";
    o<<"timing estimators, on identical events per build: cfd05 (the published 5\\%-of-clipped-peak method),\n";
    o<<"led (fixed-threshold leading edge), and hg\\_lgcfd (CFD on the LG-recovered edge). $\\sigma_t(150)$ in ps,\n";
    o<<"and the fitted floor $b$ of $\\sigma_t=a/\\sqrt{E}\\oplus b$. The adopted estimator per build is in bold.}\n";
    o<<"\\label{tab:methods}\n\\small\n\\begin{tabular}{lccc}\n\\toprule\n";
    o<<"Build & cfd05 & led & hg\\_lgcfd \\\\\n\\midrule\n";
    o<<"\\multicolumn{4}{l}{\\emph{$\\sigma_t(150)$ [ps]}}\\\\\n";
    for(int b=0;b<3;++b){ int ad=adoptedSrc(builds[b]); o<<builds[b];
        for(int i=0;i<3;++i){ double s=LL[b][i].s150; bool A=(SRCS[i]==ad);
            if(s>0) o<<Form(" & %s%.1f%s",A?"\\textbf{":"",s,A?"}":""); else o<<" & --"; } o<<" \\\\\n"; }
    o<<"\\midrule\n\\multicolumn{4}{l}{\\emph{floor }$b$\\emph{ [ps]}}\\\\\n";
    for(int b=0;b<3;++b){ int ad=adoptedSrc(builds[b]); o<<builds[b];
        for(int i=0;i<3;++i){ double bb=LL[b][i].b; bool A=(SRCS[i]==ad);
            if(LL[b][i].a>0) o<<Form(" & %s%.1f%s",A?"\\textbf{":"",bb,A?"}":""); else o<<" & --"; } o<<" \\\\\n"; }
    o<<"\\bottomrule\n\\end{tabular}\n\\end{table}\n"; o.close();
    printf("  wrote papers/timing/tab_methods.tex\n");

    printf("\n========== METHOD SURVEY: sigma_t(150) per build per method ==========\n");
    printf("build    cfd05     led    lgcfd\n");
    for(int b=0;b<3;++b){ printf("%-6s ",builds[b]); for(int i=0;i<3;++i) printf("%7.1f",LL[b][i].s150); printf("\n"); }
}

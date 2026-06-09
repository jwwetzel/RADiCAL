// ============================================================================
// methodCompare.C — motivation for the choice of timing estimator, per build.
// ----------------------------------------------------------------------------
// On IDENTICAL events (same brightest-1000 selection, same (DW-UP)/2 estimator,
// same a/sqrt(E)(+)b fit -- only the per-channel TIME SOURCE differs), compares
// the three candidate timing methods:
//   * cfd05   — CFD at 5% of the CLIPPED measured peak  (the PUBLISHED method)
//   * led     — leading edge at a FIXED absolute threshold (robust on dim/slow pulses)
//   * hg_lgcfd— CFD on the LG-recovered TRUE peak, on the steep edge below the clip
// The method that wins is set by the LIGHT REGIME: a bright, heavily-clipped build
// (DSB1/MIXED) must recover its edge (lgcfd); a dim build (LuAG) cannot use the
// published 5% foot (it sits in noise) and wants led. So the estimator choice is
// itself a fingerprint of the light-yield thesis.
//
//   source setup.sh
//   root -l -b -q 'analyze/studies/methodCompare.C+'                 // -> survey (3 builds)
//   root -l -b -q 'analyze/studies/methodCompare.C+("DSB1")'         // -> one build
// Output: figures/<year>/narrative/method_compare_<BUILD>.png (per build),
//         figures/<year>/narrative/method_survey.png (3-panel), papers/timing/tab_methods.tex
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

// the candidate methods, in display order
static const int    SRCS[3] = { RadView::kCFD05, RadView::kLED, RadView::kLGCFD };
static const char*  SLAB[3] = { "cfd05", "led", "hg_lgcfd" };
static int          adoptedSrc(const char* b){ std::string s=b; return (s=="LUAG"||s=="TENERGY")?RadView::kLED:RadView::kLGCFD; }

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

// draw the 3-method comparison for one build on the current pad; fills L[3].
static void drawBuild(const char* build, Ladder L[3]){
    BuildConfig cfg=BuildConfig::Load(radConfig(build).Data());
    int adopt=adoptedSrc(build);
    int col[3]={kRed+1,kGreen+3,kAzure+2}; int mk[3]={24,25,26};
    double ymax=55;
    for(int i=0;i<3;++i){ L[i]=ladder(build,cfg,SRCS[i]);
        for(double s:L[i].S) ymax=std::max(ymax,s); }
    ymax=std::min(ymax*1.12, 110.0);
    TH1F* fr=gPad->DrawFrame(0,0,165,ymax);
    fr->SetTitle(Form("%s;beam energy E (GeV);brightest-1000  (DW#minusUP)/2  #sigma_{t} (ps)",build));
    fr->GetYaxis()->SetTitleSize(0.045); fr->GetYaxis()->SetTitleOffset(1.2); fr->GetXaxis()->SetTitleSize(0.045);
    TLegend* lg=new TLegend(0.40,0.70,0.96,0.92); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.032);
    for(int i=0;i<3;++i){ if(L[i].E.empty())continue; bool A=(SRCS[i]==adopt);
        TGraphErrors* g=new TGraphErrors(L[i].E.size(),&L[i].E[0],&L[i].S[0],&L[i].ze[0],&L[i].Se[0]);
        g->SetMarkerStyle(A?(20+i):mk[i]); g->SetMarkerColor(col[i]); g->SetLineColor(col[i]); g->SetMarkerSize(A?1.7:1.2);
        if(L[i].a>0){ TF1* f=new TF1(Form("f%s%d",build,i),"sqrt([0]*[0]/x+[1]*[1])",20,160);
            f->SetParameters(L[i].a,L[i].b); f->SetLineColor(col[i]); f->SetLineWidth(A?3:1); f->SetLineStyle(A?1:2); f->Draw("SAME"); }
        g->Draw("P SAME");
        lg->AddEntry(g,Form("%s%s:  a=%.0f, b=%.0f#pm%.0f ps",SLAB[i],A?" (adopted)":"",L[i].a,L[i].b,L[i].be),"lp");
    }
    lg->Draw();
}

void methodCompare(const char* build="DSB1"){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    TCanvas* c=new TCanvas("mc","",900,650); c->SetLeftMargin(0.12); c->SetRightMargin(0.04);
    c->SetTopMargin(0.07); c->SetBottomMargin(0.13); c->SetGridy();
    Ladder L[3]; drawBuild(build,L);
    int adopt=adoptedSrc(build);
    printf("\n=== %s: cfd05 (published) vs led vs hg_lgcfd, brightest-1000, IDENTICAL events ===\n",build);
    printf("   E      cfd05     led    lgcfd\n");
    for(size_t k=0;k<6;++k){ double e=L[0].E.size()>k?L[0].E[k]:(L[2].E.size()>k?L[2].E[k]:0); if(e==0)continue;
        printf("  %3.0f  ",e);
        for(int i=0;i<3;++i){ double s=-1; for(size_t j=0;j<L[i].E.size();++j) if(std::fabs(L[i].E[j]-e)<1) s=L[i].S[j];
            printf("%7.1f",s); } printf("\n"); }
    for(int i=0;i<3;++i) printf("  %-9s a=%4.0f b=%4.1f+-%.1f chi2/ndf=%.1f  sigma150=%.1f%s\n",
        SLAB[i],L[i].a,L[i].b,L[i].be,L[i].chi,L[i].s150, SRCS[i]==adopt?"  <- ADOPTED":"");
    c->Print(radFigP(Form("figures/narrative/method_compare_%s.png",build)));
    printf("  wrote figures/narrative/method_compare_%s.png\n",build);
}

// the full survey: 3 builds side by side + a LaTeX summary table.
void methodSurvey(){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    const char* builds[3]={"DSB1","LUAG","MIXED"};
    Ladder LL[3][3];   // [build][method]
    TCanvas* c=new TCanvas("ms","",1500,540); c->Divide(3,1,0.004,0.004);
    for(int b=0;b<3;++b){ c->cd(b+1); gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.03);
        gPad->SetTopMargin(0.08); gPad->SetBottomMargin(0.14); gPad->SetGridy();
        drawBuild(builds[b], LL[b]);
        // also write the per-build single figure
    }
    c->cd(0); TLatex t; t.SetNDC(); t.SetTextFont(62); t.SetTextSize(0.030);
    t.DrawLatex(0.06,0.955,"Timing-method survey: cfd05 (published) vs led vs hg_lgcfd, identical brightest-1000 events");
    c->Print(radFigP("figures/narrative/method_survey.png"));
    printf("  wrote figures/narrative/method_survey.png\n");

    // per-build single figures (for the paper, one per build)
    for(int b=0;b<3;++b) methodCompare(builds[b]);

    // LaTeX summary table: sigma_t(150) + floor b per build per method
    std::ofstream o("papers/timing/tab_methods.tex");
    o<<"% auto-generated by analyze/studies/methodCompare.C (methodSurvey)\n";
    o<<"\\begin{table}[t]\n\\centering\n";
    o<<"\\caption{Brightest-1000 $(\\mathrm{DW}-\\mathrm{UP})/2$ time resolution under the three candidate\n";
    o<<"timing estimators, on identical events per build: cfd05 (the published 5\\%-of-clipped-peak method),\n";
    o<<"led (fixed-threshold leading edge), and hg\\_lgcfd (CFD on the LG-recovered edge). $\\sigma_t(150)$ in ps,\n";
    o<<"and the fitted floor $b$ of $\\sigma_t=a/\\sqrt{E}\\oplus b$. The adopted estimator per build is in bold.}\n";
    o<<"\\label{tab:methods}\n\\small\n\\begin{tabular}{lccc}\n\\toprule\n";
    o<<"Build & cfd05 & led & hg\\_lgcfd \\\\\n";
    o<<"\\midrule\n";
    o<<"\\multicolumn{4}{l}{\\emph{$\\sigma_t(150)$ [ps]}}\\\\\n";
    for(int b=0;b<3;++b){ int ad=adoptedSrc(builds[b]); o<<builds[b];
        for(int i=0;i<3;++i){ double s=LL[b][i].s150; bool A=(SRCS[i]==ad);
            if(s>0) o<<Form(" & %s%.1f%s",A?"\\textbf{":"",s,A?"}":""); else o<<" & --"; } o<<" \\\\\n"; }
    o<<"\\midrule\n";
    o<<"\\multicolumn{4}{l}{\\emph{floor }$b$\\emph{ [ps]}}\\\\\n";
    for(int b=0;b<3;++b){ int ad=adoptedSrc(builds[b]); o<<builds[b];
        for(int i=0;i<3;++i){ double bb=LL[b][i].b; bool A=(SRCS[i]==ad);
            if(LL[b][i].a>0) o<<Form(" & %s%.1f%s",A?"\\textbf{":"",bb,A?"}":""); else o<<" & --"; } o<<" \\\\\n"; }
    o<<"\\bottomrule\n\\end{tabular}\n\\end{table}\n";
    o.close();
    printf("  wrote papers/timing/tab_methods.tex\n");

    // console summary
    printf("\n========== METHOD SURVEY: sigma_t(150) per build per method ==========\n");
    printf("build    cfd05     led    lgcfd   adopted\n");
    for(int b=0;b<3;++b){ printf("%-6s ",builds[b]);
        for(int i=0;i<3;++i) printf("%7.1f",LL[b][i].s150);
        printf("   %s\n", SLAB[adoptedSrc(builds[b])==RadView::kLED?1:2]); }
}

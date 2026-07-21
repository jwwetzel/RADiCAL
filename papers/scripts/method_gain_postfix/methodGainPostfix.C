// ============================================================================
// methodGainPostfix.C — §5.3 recompute: srCFD vs cfd05 vs LED on IDENTICAL events
// with the full POST-FIX production chain. Protocol: AUDIT.md (this directory).
//   * intersection event set: all 8 ends valid + per-source median veto passed,
//     for ALL three sources; brightest-1000 by sum_lg selected ONCE per energy.
//   * width = production tebSigma (post-fix) per source on the same 1000 events.
//   * paired Poisson bootstrap (2000) of the fixed-window truncated-RMS difference.
//   * floors: a/sqrtE (+) b per source over the identical-event sigma(E).
//   * 150 GeV saturation split: nsat>=7 vs nsat<=5 (hg_saturated).
//   source setup.sh
//   root -l -b -q 'papers/scripts/method_gain_postfix/methodGainPostfix.C+'
// Output: papers/figures/method_gain_postfix/method_gain_postfix.png,
//         papers/tables/method_gain_postfix_2026-06-09.md,
//         papers/timing/tab_methods.tex (regenerated post-fix).
// ============================================================================
#include "RadView.h"
#include "RadTiming.h"
#include "DataPaths.h"
#include "PlotUtils.h"
#include "SelectionCuts.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TRandom3.h"
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
struct W { double lo,hi; };
// KERNEL CLONE: the iterative 2.5-sigma truncated window of rad::tebSigma
// (lib/physics/RadTiming.h); used here to define a COMMON window for the
// identical-event width conventions (production widths come from tebSigma itself).
static W cw(std::vector<float>& v){ W w{0,0}; if(v.size()<200)return w;
    double mu=0; for(float x:v)mu+=x; mu/=v.size();
    double sd=0; for(float x:v)sd+=(x-mu)*(x-mu); sd=std::sqrt(sd/v.size()); if(sd<1e-4)sd=0.1;
    for(int it=0;it<6;++it){ double s=0,ss=0;long n=0; for(float x:v)if(std::fabs(x-mu)<2.5*sd){s+=x;ss+=x*x;++n;}
        if(n<100)break; mu=s/n; double q=ss/n-mu*mu; if(q>0)sd=std::sqrt(q); }
    w.lo=mu-2.5*sd; w.hi=mu+2.5*sd; return w; }
static double wrms(std::vector<float>& v,W w,std::vector<int>* wt=nullptr){
    double s=0,ss=0,n=0;
    for(size_t i=0;i<v.size();++i){ float x=v[i]; if(x<w.lo||x>w.hi)continue;
        double k=wt?(*wt)[i]:1.0; if(k<=0)continue; s+=k*x; ss+=k*x*x; n+=k; }
    if(n<50)return -1; double mu=s/n,q=ss/n-mu*mu; return q>0?std::sqrt(q):-1; }

// gather the per-energy intersection set (per-source dwup with the production veto)
static std::vector<Ev> gather(double E, BuildConfig& cfg){
    std::vector<Ev> out;
    TString p=radReduced("DSB1",E); if(gSystem->AccessPathName(p.Data())) return out;
    TFile* f=TFile::Open(p.Data()); if(!f||f->IsZombie()){if(f)f->Close();return out;}
    TTree* t=(TTree*)f->Get("rad"); if(!t){f->Close();return out;}
    RadView v; v.attach(t,&cfg); double xc,yc; v.beamCenter(xc,yc);
    double R=TimingFiducialR(E), r2=R*R;
    // saturation flags via direct branch (RadView has no accessor)
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

void methodGainPostfix(){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0); gStyle->SetGridStyle(3); gStyle->SetGridColor(kGray);
    gSystem->mkdir("papers/figures/method_gain_postfix",kTRUE);
    gSystem->mkdir("papers/tables",kTRUE);
    BuildConfig cfg=BuildConfig::Load(radConfig("DSB1").Data());
    const double Es[6]={25,50,75,100,125,150};
    TRandom3 rng(20260610);
    double S[NS][6], Se[NS][6]; bool ok[6]={false};
    double d150_cfd=-1,d150_lo=-1,d150_hi=-1, d150_led=-1,d150_led_lo=0,d150_led_hi=0;
    double d150_tail=-1, d150_led_tail=-1;   // un-resampled fixed-window (tail-sensitive) centrals
    double clip150[NS]={-1,-1,-1}, unclip150[NS]={-1,-1,-1}; long nClip=0,nUnclip=0;
    printf("\n===== METHOD GAIN, POST-FIX (DSB1, identical events, production chain) =====\n");
    for(int ie=0;ie<6;++ie){
        std::vector<Ev> ev=gather(Es[ie],cfg);
        if(ev.size()<2000){ printf("  %3.0f GeV: intersection too small (%zu)\n",Es[ie],ev.size()); continue; }
        // brightest-1000 ONCE (same events for all three sources)
        std::nth_element(ev.begin(),ev.begin()+1000,ev.end(),[](const Ev&a,const Ev&b){return a.slg>b.slg;});
        ev.resize(1000);
        std::vector<float> vt[NS]; for(int s=0;s<NS;++s){ vt[s].reserve(1000); for(Ev&e:ev) vt[s].push_back(e.dwup[s]); }
        for(int s=0;s<NS;++s){ std::vector<float> tmp=vt[s]; S[s][ie]=tebSigma(tmp); Se[s][ie]=S[s][ie]/std::sqrt(2000.0); }
        ok[ie]=true;
        printf("  %3.0f GeV (N=1000 identical): srCFD %.1f  cfd05 %.1f  LED %.1f ps  | Delta(cfd05-srCFD)=%.1f\n",
               Es[ie],S[0][ie],S[1][ie],S[2][ie],S[1][ie]-S[0][ie]);
        if(Es[ie]>149){
            // paired Poisson bootstrap of the fixed-window RMS differences
            W w0=cw(vt[0]), w1=cw(vt[1]), w2=cw(vt[2]);
            // un-resampled tail-sensitive centrals (the quantity the bootstrap CI brackets;
            // convention: fixed-window truncated RMS, NOT the production Gaussian-core width)
            { double t0=wrms(vt[0],w0), t1=wrms(vt[1],w1), t2=wrms(vt[2],w2);
              if(t0>0&&t1>0) d150_tail=(t1-t0)*1000.0;
              if(t0>0&&t2>0) d150_led_tail=(t2-t0)*1000.0;
              printf("    tail-sensitive centrals @150: Delta(cfd05-srCFD)=%.1f ps, Delta(LED-srCFD)=%.1f ps\n",
                     d150_tail,d150_led_tail); }
            std::vector<double> bo1, bo2; bo1.reserve(2000); bo2.reserve(2000);
            std::vector<int> wt(1000);
            for(int b=0;b<2000;++b){ for(auto&x:wt)x=rng.Poisson(1.0);
                double s0=wrms(vt[0],w0,&wt), s1=wrms(vt[1],w1,&wt), s2=wrms(vt[2],w2,&wt);
                if(s0>0&&s1>0) bo1.push_back((s1-s0)*1000.0);
                if(s0>0&&s2>0) bo2.push_back((s2-s0)*1000.0); }
            std::sort(bo1.begin(),bo1.end()); std::sort(bo2.begin(),bo2.end());
            auto q=[&](std::vector<double>&v,double f){ return v[(size_t)(f*(v.size()-1))]; };
            d150_cfd=S[1][ie]-S[0][ie]; d150_lo=q(bo1,0.16); d150_hi=q(bo1,0.84);
            d150_led=S[2][ie]-S[0][ie]; d150_led_lo=q(bo2,0.16); d150_led_hi=q(bo2,0.84);
            printf("    paired bootstrap @150: Delta(cfd05-srCFD) 68%% [%.1f, %.1f] ps; Delta(LED-srCFD) 68%% [%.1f, %.1f] ps\n",
                   d150_lo,d150_hi,d150_led_lo,d150_led_hi);
            // saturation split
            std::vector<float> c[NS], u[NS];
            for(Ev&e:ev){ for(int s=0;s<NS;++s){ if(e.nsat>=7)c[s].push_back(e.dwup[s]); else if(e.nsat<=5)u[s].push_back(e.dwup[s]); } }
            nClip=c[0].size(); nUnclip=u[0].size();
            for(int s=0;s<NS;++s){ std::vector<float> t1=c[s],t2=u[s];
                clip150[s]= (t1.size()>=300)? tebSigma(t1):-1; unclip150[s]=(t2.size()>=300)? tebSigma(t2):-1; }
            printf("    saturation split @150: nsat>=7 (N=%ld): srCFD %.1f cfd05 %.1f | nsat<=5 (N=%ld): srCFD %.1f cfd05 %.1f\n",
                   nClip,clip150[0],clip150[1],nUnclip,unclip150[0],unclip150[1]);
        }
    }
    // floors per source
    double A[NS],Ae[NS],B[NS],Be[NS],Chi[NS]; int Ndf[NS];
    for(int s=0;s<NS;++s){ std::vector<double> E,Y,Yer,z;
        for(int ie=0;ie<6;++ie) if(ok[ie]&&S[s][ie]>0){E.push_back(Es[ie]);Y.push_back(S[s][ie]);Yer.push_back(Se[s][ie]);z.push_back(0);}
        TGraphErrors g(E.size(),&E[0],&Y[0],&z[0],&Yer[0]);
        TF1 f("f","sqrt([0]*[0]/x+[1]*[1])",E.front(),E.back()); f.SetParameters(250,18); g.Fit(&f,"QN");
        A[s]=std::fabs(f.GetParameter(0)); Ae[s]=f.GetParError(0);
        B[s]=std::fabs(f.GetParameter(1)); Be[s]=f.GetParError(1);
        Chi[s]=f.GetChisquare(); Ndf[s]=f.GetNDF();
        if(Ndf[s]>0&&Chi[s]/Ndf[s]>1){double k=std::sqrt(Chi[s]/Ndf[s]);Ae[s]*=k;Be[s]*=k;}
        printf("  floor %-6s: a=%.0f+-%.0f  b=%.1f+-%.1f  chi2/ndf=%.1f/%d\n",SNM[s],A[s],Ae[s],B[s],Be[s],Chi[s],Ndf[s]);
    }
    double gainPct=100.0*(S[1][5]-S[0][5])/S[1][5];
    printf("\n  HEADLINE: srCFD improves identical-event sigma(150) by %.1f ps (%.0f%%) vs cfd05; floor %.1f->%.1f ps\n",
           S[1][5]-S[0][5],gainPct,B[1],B[0]);

    // ---- markdown table ----
    { std::ofstream md("papers/tables/method_gain_postfix_2026-06-09.md");
      md << "# Method gain, post-fix (DSB1, identical events) — 2026-06-09\n\n";
      md << "Generated by `papers/scripts/method_gain_postfix/methodGainPostfix.C` (audit: AUDIT.md).\n";
      md << "Identical-event design: intersection of per-source validity+veto; brightest-1000 selected once.\n";
      md << "Width = production post-fix tebSigma. SUPERSEDES the pre-fix 3.1 ps / 22.1->19.5 / tab_methods values.\n\n";
      md << "| E (GeV) | srCFD (ps) | cfd05 (ps) | LED (ps) | cfd05−srCFD (ps) |\n|---|---|---|---|---|\n";
      char b[320];
      for(int ie=0;ie<6;++ie){ if(!ok[ie])continue;
          snprintf(b,sizeof(b),"| %.0f | %.1f | %.1f | %.1f | %+.1f |\n",Es[ie],S[0][ie],S[1][ie],S[2][ie],S[1][ie]-S[0][ie]); md<<b; }
      md << "\nNOTE — two width conventions (see docs/STATS_CONVENTIONS.md): table widths = production\n";
      md << "tebSigma (Gaussian-core, the convention of every paper resolution); the bootstrap CI below\n";
      md << "brackets the TAIL-SENSITIVE (fixed-window truncated-RMS) difference, whose central value is\n";
      md << "listed beside it. The two are different observables and are quoted separately in Sec. 5.3.\n";
      snprintf(b,sizeof(b),"\n- 150 GeV core-width gain (cfd05−srCFD, tebSigma): **%.1f ps (%.0f%%)**\n",
               S[1][5]-S[0][5],gainPct); md<<b;
      snprintf(b,sizeof(b),"- 150 GeV tail-sensitive gain (cfd05−srCFD, fixed-window RMS): %.1f ps, paired-bootstrap 68%% CI [%.1f, %.1f] ps\n",
               d150_tail,d150_lo,d150_hi); md<<b;
      snprintf(b,sizeof(b),"- 150 GeV LED−srCFD: core %+.1f ps; tail-sensitive %+.1f ps, 68%% CI [%.1f, %.1f] ps\n",
               d150_led,d150_led_tail,d150_led_lo,d150_led_hi); md<<b;
      { char u0[24],u1[24];
        if(unclip150[0]<0) snprintf(u0,sizeof(u0),"n/a"); else snprintf(u0,sizeof(u0),"%.1f",unclip150[0]);
        if(unclip150[1]<0) snprintf(u1,sizeof(u1),"n/a"); else snprintf(u1,sizeof(u1),"%.1f",unclip150[1]);
        snprintf(b,sizeof(b),"- Saturation split @150: fully clipped (nsat≥7, N=%ld): srCFD %.1f vs cfd05 %.1f ps; less clipped (nsat≤5, N=%ld): %s vs %s (subset < 300: width not quoted; the less-clipped comparison lives in the 75 GeV full-fiducial split)\n",
               nClip,clip150[0],clip150[1],nUnclip,u0,u1); md<<b; }
      md << "- Floors (identical-event fits):\n";
      for(int s=0;s<NS;++s){ snprintf(b,sizeof(b),"  - %s: a=%.0f±%.0f ps·√GeV, b=%.1f±%.1f ps (χ²/ndf %.1f/%d)\n",
               SNM[s],A[s],Ae[s],B[s],Be[s],Chi[s],Ndf[s]); md<<b; }
      md << "- Closure: at 25 GeV (4% clipped) srCFD vs cfd05 difference is the table's first row.\n";
      md.close(); printf("  wrote papers/tables/method_gain_postfix_2026-06-09.md\n"); }

    // ---- regenerate tab_methods.tex (post-fix) ----
    { std::ofstream tx("papers/timing/tab_methods.tex");
      tx << "% auto-generated by papers/scripts/method_gain_postfix/methodGainPostfix.C (POST-FIX, 2026-06-09)\n";
      tx << "% identical-event comparison, DSB1, production estimator chain. Not \\input by the manuscript\n";
      tx << "% by default; numbers match Sec. 5.3 (method gain).\n";
      tx << "\\begin{table}[t]\n\\centering\n";
      tx << "\\caption{Identical-event timing comparison of the three per-channel estimators (DSB1,\n";
      tx << "brightest-1000, production selection). $\\sigma_t(150)$ and the fitted floor $b$ per source.}\n";
      tx << "\\label{tab:methods}\n\\small\n\\begin{tabular}{lccc}\n\\toprule\n";
      tx << " & srCFD & cfd05 & LED \\\\\n\\midrule\n";
      char b2[256];
      snprintf(b2,sizeof(b2),"$\\sigma_t(150)$ (ps) & $%.1f$ & $%.1f$ & $%.1f$ \\\\\n",S[0][5],S[1][5],S[2][5]); tx<<b2;
      snprintf(b2,sizeof(b2),"floor $b$ (ps) & $%.1f\\pm%.1f$ & $%.1f\\pm%.1f$ & $%.1f\\pm%.1f$ \\\\\n",B[0],Be[0],B[1],Be[1],B[2],Be[2]); tx<<b2;
      snprintf(b2,sizeof(b2),"$a$ (ps$\\sqrt{\\mathrm{GeV}}$) & $%.0f\\pm%.0f$ & $%.0f\\pm%.0f$ & $%.0f\\pm%.0f$ \\\\\n",A[0],Ae[0],A[1],Ae[1],A[2],Ae[2]); tx<<b2;
      tx << "\\bottomrule\n\\end{tabular}\n\\end{table}\n";
      tx.close(); printf("  regenerated papers/timing/tab_methods.tex (post-fix)\n"); }

    // ---- figure ----
    TCanvas* c=new TCanvas("mg","",900,860);
    TPad* p1=new TPad("p1","",0,0.34,1,1); TPad* p2=new TPad("p2","",0,0,1,0.34);
    p1->SetBottomMargin(0.02); p1->SetLeftMargin(0.13); p1->SetRightMargin(0.04); p1->SetTopMargin(0.05); p1->SetGridy(); p1->SetLogx();
    p2->SetTopMargin(0.04); p2->SetBottomMargin(0.32); p2->SetLeftMargin(0.13); p2->SetRightMargin(0.04); p2->SetGridy(); p2->SetLogx();
    p1->Draw(); p2->Draw();
    p1->cd();
    { TH1F* fr=gPad->DrawFrame(20,18,200,60);
      fr->GetYaxis()->SetTitle("identical-event #sigma_{t} (ps)"); fr->GetYaxis()->SetTitleSize(0.055); fr->GetYaxis()->SetTitleOffset(1.1);
      fr->GetXaxis()->SetLabelSize(0);
      int cols[NS]={kAzure+2,kBlack,kRed+1}; int mks[NS]={20,24,25};
      TLegend* lg=new TLegend(0.52,0.62,0.92,0.90); lg->SetBorderSize(0);lg->SetFillStyle(0);lg->SetTextSize(0.045);
      for(int s=0;s<NS;++s){ std::vector<double> E,Y,Yer,z;
          for(int ie=0;ie<6;++ie) if(ok[ie]){E.push_back(Es[ie]);Y.push_back(S[s][ie]);Yer.push_back(Se[s][ie]);z.push_back(0);}
          TGraphErrors* g=new TGraphErrors(E.size(),&E[0],&Y[0],&z[0],&Yer[0]);
          g->SetMarkerStyle(mks[s]);g->SetMarkerColor(cols[s]);g->SetLineColor(cols[s]);g->SetMarkerSize(1.4);
          TF1* f=new TF1(Form("fmg%d",s),"sqrt([0]*[0]/x+[1]*[1])",25,150); f->SetParameters(A[s],B[s]);
          f->SetLineColor(cols[s]); f->SetLineWidth(2); f->SetLineStyle(s==0?1:2); f->Draw("SAME");
          g->Draw("P SAME");
          lg->AddEntry(g,Form("%s: b=%.1f#pm%.1f ps",SNM[s],B[s],Be[s]),"lp"); }
      lg->Draw();
      TLatex tx; tx.SetNDC(); tx.SetTextSize(0.040); tx.SetTextColor(kGray+3);
      tx.DrawLatex(0.16,0.16,"DSB1, brightest-1000, SAME events for all estimators"); }
    p2->cd();
    { TH1F* fr=gPad->DrawFrame(20,-3,200,8);
      fr->GetYaxis()->SetTitle("#Delta#sigma vs srCFD (ps)"); fr->GetYaxis()->SetTitleSize(0.085); fr->GetYaxis()->SetTitleOffset(0.70);
      fr->GetYaxis()->SetLabelSize(0.075); fr->GetYaxis()->SetNdivisions(505);
      fr->GetXaxis()->SetTitle("beam energy E (GeV)"); fr->GetXaxis()->SetTitleSize(0.11); fr->GetXaxis()->SetLabelSize(0.075);
      fr->GetXaxis()->SetMoreLogLabels(); fr->GetXaxis()->SetNoExponent();
      TLine* l0=new TLine(20,0,200,0); l0->SetLineColor(kGray+2); l0->SetLineStyle(2); l0->Draw();
      std::vector<double> E,D1,D2,z;
      for(int ie=0;ie<6;++ie) if(ok[ie]){E.push_back(Es[ie]);D1.push_back(S[1][ie]-S[0][ie]);D2.push_back(S[2][ie]-S[0][ie]);z.push_back(0);}
      TGraph* g1=new TGraph(E.size(),&E[0],&D1[0]); g1->SetMarkerStyle(24);g1->SetMarkerColor(kBlack);g1->SetLineColor(kBlack);g1->SetMarkerSize(1.3);g1->Draw("PL SAME");
      TGraph* g2=new TGraph(E.size(),&E[0],&D2[0]); g2->SetMarkerStyle(25);g2->SetMarkerColor(kRed+1);g2->SetLineColor(kRed+1);g2->SetMarkerSize(1.3);g2->Draw("PL SAME");
      TLatex tx; tx.SetNDC(); tx.SetTextSize(0.075);
      // core-width gain only; the tail-sensitive gain + its bootstrap CI live in Sec. 5.3 (pairing
      // the core central value with the tail-convention CI in one stamp read as self-contradictory)
      tx.DrawLatex(0.33,0.88,Form("@150: cfd05#minussrCFD = %+.1f ps (core width)",S[1][5]-S[0][5])); }
    // paper convention (format pass 2026-06-09): no internal super-title; the LaTeX caption carries it
    c->cd(0);
    c->Print("papers/figures/method_gain_postfix/method_gain_postfix.png");
    printf("  wrote papers/figures/method_gain_postfix/method_gain_postfix.png\n");
}

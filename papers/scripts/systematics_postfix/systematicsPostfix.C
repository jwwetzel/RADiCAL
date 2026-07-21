// ============================================================================
// systematicsPostfix.C — regenerate papers/timing/tab_systematics.tex with the
// FULL POST-FIX production chain. Protocol: AUDIT.md (this directory).
//   * per build @150 GeV: one pass storing per-event end times/peaks (adopted
//     source); variations applied post-hoc on identical stored events:
//     K=500/2000, r_fid=2.5/3.5, MCP<700, HG>=30 mV, veto W=1.5/3.0 ns.
//   * width = production post-fix tebSigma; nominal = production cuts, K=1000.
//   * DSB1 floor block: fresh production timingBrightestK 6-energy scan ->
//     b +- stat (sqrt(chi2/ndf) inflated) and the fit-range (drop 25 GeV) shift.
//   * methodological caveats go to table NOTES, not fake-precision rows.
//   source setup.sh
//   root -l -b -q 'papers/scripts/systematics_postfix/systematicsPostfix.C+'
// Output: papers/timing/tab_systematics.tex (regenerated),
//         papers/tables/systematics_postfix_2026-06-09.md, log (this dir).
// ============================================================================
#include "RadView.h"
#include "RadTiming.h"
#include "DataPaths.h"
#include "PlotUtils.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TF1.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TBox.h"
#include "TLine.h"
#include "TLatex.h"
#include <vector>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <cstdio>
using namespace rad;

struct Ev { float t[8], pk[8], r2, slg, mcp; bool istim[8]; };

static std::vector<Ev> gather150(const char* build,int src,BuildConfig& cfg){
    std::vector<Ev> out;
    TString p=radReduced(build,150); if(gSystem->AccessPathName(p.Data())) return out;
    TFile* f=TFile::Open(p.Data()); if(!f||f->IsZombie()){if(f)f->Close();return out;}
    TTree* t=(TTree*)f->Get("rad"); if(!t){f->Close();return out;}
    RadView v; v.attach(t,&cfg); double xc,yc; v.beamCenter(xc,yc);
    Long64_t N=v.entries(); out.reserve(N/4);
    for(Long64_t i=0;i<N;++i){ v.get(i);
        if(!v.wc_ok()) continue;
        Ev e; e.mcp=v.mcp1_peak();
        if(e.mcp<kMCP1_minPeak||e.mcp>kMCP1_maxPeak) continue;   // widest window; tighter applied post-hoc
        double dx=v.x_trk()-xc,dy=v.y_trk()-yc; e.r2=(float)(dx*dx+dy*dy);
        if(e.r2>=3.5*3.5) continue;                               // widest radius; tighter post-hoc
        e.slg=(float)v.sum_lg(); bool any=false;
        for(int c=0;c<8;++c){ e.istim[c]=v.is_timing(c); e.pk[c]=v.hg_peak(c);
            float tc=e.istim[c]? v.timeOf(c,src) : -1e9f; e.t[c]=tc; if(e.istim[c]&&tc>-1e5f)any=true; }
        if(any) out.push_back(e);
    }
    f->Close(); return out;
}

// production-equivalent sigma under a variant (post-hoc cuts on the stored events)
// KERNEL CLONE: rad::eventDWUP + brightest-K + rad::tebSigma re-run on STORED events with
// cuts as arguments (the whole point: post-hoc variants on identical events). Nominal
// args (3.0,750,20,2.0,1000) reproduce the production headline 25.7 — the equivalence proof.
static double variantSigma(std::vector<Ev>& EV,double rmax,double mcpHi,float hgThr,float W,int K){
    std::vector<std::pair<float,float>> sd; sd.reserve(EV.size());
    for(Ev& e:EV){
        if(e.r2>=rmax*rmax||e.mcp>mcpHi) continue;
        float tt[8]; bool okc[8]; int m=0; float a[8];
        for(int c=0;c<8;++c){ okc[c]=e.istim[c]&&e.pk[c]>=hgThr&&e.t[c]>-1e5f; if(okc[c]){tt[c]=e.t[c];a[m++]=e.t[c];} }
        if(m<2) continue;
        std::nth_element(a,a+m/2,a+m); float med=a[m/2];
        double ds=0,us=0;int dn=0,un=0;
        for(int c=0;c<8;++c){ if(!okc[c]||std::fabs(tt[c]-med)>=W) continue;
            if(c<4){ds+=tt[c];++dn;} else {us+=tt[c];++un;} }
        if(dn<1||un<1) continue;
        sd.push_back({e.slg,0.5f*(float)(ds/dn-us/un)});
    }
    if((int)sd.size()<K) return -1;
    std::nth_element(sd.begin(),sd.begin()+K,sd.end(),[](const std::pair<float,float>&a,const std::pair<float,float>&b){return a.first>b.first;});
    std::vector<float> vt; vt.reserve(K); for(int i=0;i<K;++i)vt.push_back(sd[i].second);
    return tebSigma(vt);
}

void systematicsPostfix(){
    gSystem->mkdir("papers/tables",kTRUE);
    struct B { const char* build; int src; } BS[4]={
        {"DSB1",RadView::kLGCFD},{"TENERGY",RadView::kLED},{"MIXED",RadView::kLGCFD},{"LUAG",RadView::kLED}};
    const char* VN[8]={"$K{=}500$","$K{=}2000$","$r_{\\mathrm{fid}}{=}2.5$~mm","$r_{\\mathrm{fid}}{=}3.5$~mm",
                       "MCP$\\,{<}700$~mV","HG$\\,{\\ge}30$~mV","veto $1.5$~ns","veto $3.0$~ns"};
    double nom[4], shift[4][8], tot[4];
    printf("\n===== POST-FIX systematics @150 GeV (production chain; post-hoc variants) =====\n");
    for(int ib=0;ib<4;++ib){
        BuildConfig cfg=BuildConfig::Load(radConfig(BS[ib].build).Data());
        std::vector<Ev> EV=gather150(BS[ib].build,BS[ib].src,cfg);
        nom[ib]=variantSigma(EV,3.0,750,20,2.0f,1000);
        double v[8]={
            variantSigma(EV,3.0,750,20,2.0f,500),  variantSigma(EV,3.0,750,20,2.0f,2000),
            variantSigma(EV,2.5,750,20,2.0f,1000), variantSigma(EV,3.5,750,20,2.0f,1000),
            variantSigma(EV,3.0,700,20,2.0f,1000), variantSigma(EV,3.0,750,30,2.0f,1000),
            variantSigma(EV,3.0,750,20,1.5f,1000), variantSigma(EV,3.0,750,20,3.0f,1000)};
        double rms=0; int n=0;
        for(int k=0;k<8;++k){
            if(v[k]<=0||nom[ib]<=0)   // insurance (never trips on committed data): a failed
                printf("  WARNING: %s variant %d FAILED (sigma<0) — recorded as 0 shift; do NOT quote without investigating\n",BS[ib].build,k);
            shift[ib][k]=(v[k]>0&&nom[ib]>0)?v[k]-nom[ib]:0; rms+=shift[ib][k]*shift[ib][k]; ++n; }
        tot[ib]=std::sqrt(rms/n);
        printf("  %-8s nominal %.1f ps | shifts:",BS[ib].build,nom[ib]);
        for(int k=0;k<8;++k) printf(" %+0.1f",shift[ib][k]);
        printf(" | total RMS %.1f\n",tot[ib]);
    }
    // DSB1 floor block (fresh production 6-E scan)
    double bAll=0,beAll=0,bDrop=0;
    { BuildConfig cfg=BuildConfig::Load(radConfig("DSB1").Data());
      const double Es[6]={25,50,75,100,125,150};
      std::vector<double> E,S,Se,z;
      for(double e:Es){ TString p=radReduced("DSB1",e); if(gSystem->AccessPathName(p.Data()))continue;
          TFile* f=TFile::Open(p.Data()); TTree* t=(TTree*)f->Get("rad"); RadView v; v.attach(t,&cfg);
          TimingResult r=timingBrightestK(v,e,RadView::kLGCFD,1000);
          if(r.sigma_ps>0){E.push_back(e);S.push_back(r.sigma_ps);Se.push_back(r.sigma_ps/std::sqrt(2000.0));z.push_back(0);}
          f->Close(); }
      auto fitb=[&](size_t lo)->std::pair<double,double>{
          std::vector<double> E2(E.begin()+lo,E.end()),S2(S.begin()+lo,S.end()),Se2(Se.begin()+lo,Se.end()),z2(E2.size(),0);
          TGraphErrors g(E2.size(),&E2[0],&S2[0],&z2[0],&Se2[0]);
          TF1 f("f","sqrt([0]*[0]/x+[1]*[1])",E2.front(),E2.back()); f.SetParameters(250,18); g.Fit(&f,"QN");
          double b=std::fabs(f.GetParameter(1)),be=f.GetParError(1);
          if(f.GetNDF()>0&&f.GetChisquare()/f.GetNDF()>1) be*=std::sqrt(f.GetChisquare()/f.GetNDF());
          return {b,be}; };
      auto r1=fitb(0), r2=fitb(1);
      bAll=r1.first; beAll=r1.second; bDrop=r2.first;
      printf("  DSB1 floor: b(25-150)=%.1f+-%.1f, b(50-150)=%.1f  (fit-range shift %+0.1f)\n",bAll,beAll,bDrop,bDrop-bAll);
    }

    // ---- regenerate tab_systematics.tex ----
    std::ofstream tx("papers/timing/tab_systematics.tex");
    tx << "% auto-generated by papers/scripts/systematics_postfix/systematicsPostfix.C (POST-FIX, 2026-06-09)\n";
    tx << "% production chain: post-fix tebSigma + in-event veto + continuous fiducial; variants post-hoc\n";
    tx << "% on identical stored events. Supersedes the pre-fix paperSystematics.C table.\n";
    tx << "\\begin{table}[t]\n\\centering\n";
    tx << "\\caption{Selection-systematic budget for the brightest-1000 $(\\mathrm{DW}-\\mathrm{UP})/2$\n";
    tx << "time resolution at $150$~GeV, recomputed with the production estimator chain. Each row is the\n";
    tx << "shift of $\\sigma_t(150)$ under one selection variation on identical events; the total is the\n";
    tx << "RMS of the shifts. The floor block gives the DSB1 fit floor and its fit-range stability.\n";
    tx << "The timing-method choice (srCFD vs.\\ cfd05/LED) is treated in Sec.~\\ref{sec:method-gain}\n";
    tx << "with its own paired uncertainty and is not folded in here. Notes: the MIXED module-wide\n";
    tx << "estimator is ill-posed (Sec.~\\ref{sec:mixed}; its same-shower ratio carries $\\pm0.05$\n";
    tx << "energy-period scatter and $<0.01$ run-jackknife spread); TENERGY times with three\n";
    tx << "capillaries ($\\sqrt{4/3}$ statistical penalty); LuAG and TENERGY floors are extrapolations\n";
    tx << "(Sec.~\\ref{sec:thesis}).}\n";
    tx << "\\label{tab:syst}\n\\footnotesize\n\\setlength{\\tabcolsep}{4pt}\n\\begin{tabular}{@{}lcccc@{}}\n\\toprule\n";
    tx << "Variation & DSB1 & TENERGY & MIXED & LUAG \\\\\n\\midrule\n";
    tx << "\\multicolumn{5}{l}{\\emph{shift in }$\\sigma_t(150)$\\emph{ [ps]}}\\\\\n";
    char b[256];
    for(int k=0;k<8;++k){ snprintf(b,sizeof(b),"%s & $%+0.1f$ & $%+0.1f$ & $%+0.1f$ & $%+0.1f$ \\\\\n",
        VN[k],shift[0][k],shift[1][k],shift[2][k],shift[3][k]); tx<<b; }
    tx << "\\midrule\n";
    snprintf(b,sizeof(b),"nominal $\\sigma_t(150)$ & $%.1f$ & $%.1f$ & $%.1f$ & $%.1f$ \\\\\n",nom[0],nom[1],nom[2],nom[3]); tx<<b;
    snprintf(b,sizeof(b),"total syst. & $\\pm%.1f$ & $\\pm%.1f$ & $\\pm%.1f$ & $\\pm%.1f$ \\\\\n",tot[0],tot[1],tot[2],tot[3]); tx<<b;
    tx << "\\midrule\n";
    snprintf(b,sizeof(b),"DSB1 floor $b$ (stat) & \\multicolumn{4}{l}{$%.1f\\pm%.1f$~ps} \\\\\n",bAll,beAll); tx<<b;
    snprintf(b,sizeof(b),"DSB1 $b$, fit $50$--$150$ only & \\multicolumn{4}{l}{$%+0.1f$~ps shift} \\\\\n",bDrop-bAll); tx<<b;
    tx << "\\bottomrule\n\\end{tabular}\n\\end{table}\n";
    tx.close();
    printf("  regenerated papers/timing/tab_systematics.tex\n");

    // md mirror
    std::ofstream md("papers/tables/systematics_postfix_2026-06-09.md");
    md << "# Post-fix selection systematics @150 GeV (2026-06-09)\n\nGenerated by systematicsPostfix.C; mirrors papers/timing/tab_systematics.tex.\n\n";
    md << "| variation | DSB1 | TENERGY | MIXED | LUAG |\n|---|---|---|---|---|\n";
    const char* MN[8]={"K=500","K=2000","r=2.5","r=3.5","MCP<700","HG>=30","veto 1.5ns","veto 3.0ns"};
    for(int k=0;k<8;++k){ snprintf(b,sizeof(b),"| %s | %+0.1f | %+0.1f | %+0.1f | %+0.1f |\n",MN[k],shift[0][k],shift[1][k],shift[2][k],shift[3][k]); md<<b; }
    snprintf(b,sizeof(b),"| **nominal** | %.1f | %.1f | %.1f | %.1f |\n| **total RMS** | %.1f | %.1f | %.1f | %.1f |\n",
             nom[0],nom[1],nom[2],nom[3],tot[0],tot[1],tot[2],tot[3]); md<<b;
    snprintf(b,sizeof(b),"\nDSB1 floor b = %.1f ± %.1f ps (25–150); fit-range (50–150) shift %+0.1f ps.\n",bAll,beAll,bDrop-bAll); md<<b;
    md.close();
    printf("  wrote papers/tables/systematics_postfix_2026-06-09.md\n");

    // ---- stability figure: POST-FIX replacement for the pre-fix paperSystematics.C A.10 ----
    // Paper convention: no internal super-title (the LaTeX caption carries it). No baked floor
    // numbers — Tables 1-2 are authoritative; the only annotated number is tot[] computed above,
    // so the figure agrees with tab_systematics.tex by construction.
    { const char* XL[9]={"nominal","K=500","K=2000","r_{fid}=2.5","r_{fid}=3.5","MCP<700","HG#geq30","veto 1.5","veto 3.0"};
      const char* PT[4]={"DSB1 (srCFD)","TENERGY (LED)","MIXED (module-wide ref)","LUAG (LED)"};
      int cols[4]={kRData,kOrange+8,kGray+2,kRGreen+1};
      TCanvas* c=new TCanvas("syst","",1400,1000); c->Divide(2,2,0.004,0.004);
      for(int ib=0;ib<4;++ib){ c->cd(ib+1);
          gPad->SetLeftMargin(0.13); gPad->SetRightMargin(0.04); gPad->SetTopMargin(0.09); gPad->SetBottomMargin(0.22); gPad->SetGridy();
          double mx=0; for(int k=0;k<8;++k) mx=std::max(mx,std::fabs(shift[ib][k]));
          double pad=std::max(2.2*tot[ib],1.25*mx)+0.6;
          TH1F* fr=gPad->DrawFrame(-0.7,nom[ib]-pad,8.7,nom[ib]+pad);
          fr->SetTitle(";;brightest-1000 #sigma_{t}(150)  (ps)");
          fr->GetYaxis()->SetTitleSize(0.055); fr->GetYaxis()->SetTitleOffset(1.15); fr->GetYaxis()->SetLabelSize(0.05);
          fr->GetYaxis()->SetNdivisions(506);
          for(int k=0;k<9;++k) fr->GetXaxis()->SetBinLabel(fr->GetXaxis()->FindBin((double)k),XL[k]);
          fr->GetXaxis()->LabelsOption("v"); fr->GetXaxis()->SetLabelSize(0.055);
          TBox* bx=new TBox(-0.7,nom[ib]-tot[ib],8.7,nom[ib]+tot[ib]);
          bx->SetFillColorAlpha(cols[ib],0.18); bx->Draw();
          TLine* ln=new TLine(-0.7,nom[ib],8.7,nom[ib]); ln->SetLineColor(cols[ib]); ln->SetLineWidth(2); ln->Draw();
          double xs[8],ys[8];
          for(int k=0;k<8;++k){ xs[k]=k+1; ys[k]=nom[ib]+shift[ib][k]; }
          TGraph* g=new TGraph(8,xs,ys); g->SetMarkerStyle(20); g->SetMarkerColor(cols[ib]); g->SetMarkerSize(1.5); g->Draw("P SAME");
          TGraph* gn=new TGraph(1); gn->SetPoint(0,0,nom[ib]); gn->SetMarkerStyle(29); gn->SetMarkerSize(2.6); gn->SetMarkerColor(kBlack); gn->Draw("P SAME");
          TLatex tl; tl.SetNDC(); tl.SetTextSize(0.055); tl.SetTextColor(cols[ib]);
          tl.DrawLatex(0.58,0.84,Form("total syst #pm%.1f ps",tot[ib]));
          DrawPadTitle(PT[ib]);
      }
      gSystem->mkdir("papers/figures/systematics_postfix",kTRUE);
      c->Print("papers/figures/systematics_postfix/systematics_stability.png");
      c->Print("papers/timing/figs/systematics.png");
      printf("  wrote POST-FIX stability figure -> papers/timing/figs/systematics.png (replaces pre-fix A.10)\n");
    }
}

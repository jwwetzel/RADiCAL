// ============================================================================
// paperSystematics.C — Paper 1 systematic-uncertainty budget, all four builds.
// ----------------------------------------------------------------------------
// The headline brightest-1000 (DW-UP)/2 sigma_t = a/sqrt(E) (+) b is recomputed
// under deliberate variations of every selection knob, and the spread is quoted
// as the systematic on sigma_t(150), on a, and on the floor b. The variations:
//   * brightest-N   K = 500 / 1000(nom) / 2000      (selection tightness)
//   * fiducial r    x0.83 / x1.0(nom) / x1.17       (2.5 / 3.0 / 3.5 mm)
//   * timing source robust(nom) <-> cfd05           (the published method)
//   * MCP window    upper 750(nom) -> 700           (reference saturation)
//   * HG threshold  20(nom) -> 30 mV                (per-channel noise floor)
//   * fit range     25-150(nom) -> 50-150           (low-E lever arm)
//
// Method fidelity: sigmaKnob() at NOMINAL knobs reproduces rad::timingBrightestK
// bit-for-bit (same cuts, same nth_element brightest-K, same tebSigma core fit) --
// validated and printed per build. One event pass per (build,energy) caches every
// per-event quantity the knobs need, so all variations are applied in-memory.
//
// The reduced files are the COMPLETE statistics (the reducer processed all raw
// events), so this runs at full stats locally; identical on Argon.
//
//   source setup.sh
//   root -l -b -q 'analyze/studies/paperSystematics.C+'
// Output: figures/narrative/systematics.png, papers/timing/tab_systematics.tex
// ============================================================================
#include "RadView.h"
#include "RadTiming.h"        // rad::tebSigma, rad::timingBrightestK (validation)
#include "DataPaths.h"        // radReduced
#include "BuildConfig.h"
#include "SelectionCuts.h"    // TimingFiducialR, kMCP1_*, kHG_minPeak
#include "PlotUtils.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TBox.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <fstream>
using namespace rad;

static int  robustSrc(const char* b){ std::string s=b; if(s=="LUAG"||s=="TENERGY") return RadView::kLED; return RadView::kLGCFD; }
static const char* srcName(int s){ return s==RadView::kLED?"led":(s==RadView::kLGCFD?"lgcfd":"cfd05"); }

// per-event cache: everything any knob needs (timing-role baked into hg[c]=-1)
struct Ev { float slg, x, y, mcp; float hg[8]; float tR[8]; float tC[8]; };

// one pass over a (build,energy) file: cache events + the sum_lg-weighted center
static void loadFile(const char* path, BuildConfig& cfg, int srcR,
                     std::vector<Ev>& out, double& xc, double& yc){
    out.clear(); xc=yc=0;
    TFile* fp=TFile::Open(path); if(!fp||fp->IsZombie()){ if(fp)fp->Close(); return; }
    TTree* t=(TTree*)fp->Get("rad"); if(!t){ fp->Close(); return; }
    RadView v; v.attach(t,&cfg);
    Long64_t N=v.entries(); out.reserve(N);
    double wx=0,wy=0,w=0;
    for(Long64_t i=0;i<N;++i){ v.get(i);
        if(!v.wc_ok()) continue;
        double slg=v.sum_lg();
        if(v.mcp1_peak()>=50.0f && slg>300.0){ wx+=v.x_trk()*slg; wy+=v.y_trk()*slg; w+=slg; }
        Ev e; e.slg=(float)slg; e.x=v.x_trk(); e.y=v.y_trk(); e.mcp=v.mcp1_peak();
        for(int c=0;c<8;++c){ bool tim=v.is_timing(c);
            e.hg[c]= tim? v.hg_peak(c) : -1.f;
            e.tR[c]= tim? v.timeOf(c,srcR)         : kNoTime;
            e.tC[c]= tim? v.timeOf(c,RadView::kCFD05) : kNoTime; }
        out.push_back(e);
    }
    xc=(w>0)?wx/w:0; yc=(w>0)?wy/w:0;
    fp->Close();
}

// in-memory mirror of timingBrightestK with every cut exposed.
//  useC: false=robust source (tR), true=cfd05 (tC).
static double sigmaKnob(const std::vector<Ev>& ev, double xc, double yc,
                        double rFid, int K, bool useC,
                        float mcpHi, float hgMin){
    double r2=rFid*rFid;
    std::vector<std::pair<float,float>> sd; sd.reserve(ev.size()/4+16);
    for(const Ev& e : ev){
        if(e.mcp<kMCP1_minPeak || e.mcp>mcpHi) continue;
        double dx=e.x-xc, dy=e.y-yc; if(dx*dx+dy*dy>=r2) continue;
        double ds=0,us=0; int dn=0,un=0;
        for(int c=0;c<4;++c) if(e.hg[c]>=hgMin){ float tc=useC?e.tC[c]:e.tR[c]; if(tc>-1e5f){ds+=tc;++dn;} }
        for(int c=4;c<8;++c) if(e.hg[c]>=hgMin){ float tc=useC?e.tC[c]:e.tR[c]; if(tc>-1e5f){us+=tc;++un;} }
        if(dn<1||un<1) continue;
        sd.push_back({ e.slg, 0.5f*(float)(ds/dn-us/un) });
    }
    if((int)sd.size()<K) return -1;
    std::nth_element(sd.begin(), sd.begin()+K, sd.end(),
                     [](const std::pair<float,float>&a,const std::pair<float,float>&b){ return a.first>b.first; });
    std::vector<float> vt; vt.reserve(K); for(int i=0;i<K;++i) vt.push_back(sd[i].second);
    return tebSigma(vt);
}

// fit sigma(E) -> a, b (PDG-scaled b error); optionally restrict to E>=Emin
static void fitAB(const std::vector<double>& E, const std::vector<double>& S,
                  const std::vector<double>& Se, double Emin,
                  double& a, double& b, double& be){
    std::vector<double> e,s,se,ze;
    for(size_t i=0;i<E.size();++i) if(E[i]>=Emin-1){ e.push_back(E[i]); s.push_back(S[i]); se.push_back(Se[i]); ze.push_back(0); }
    a=b=be=0; if(e.size()<3) return;
    TGraphErrors g(e.size(),&e[0],&s[0],&ze[0],&se[0]);
    TF1 f("fab","sqrt([0]*[0]/x+[1]*[1])",20,160); f.SetParameters(220,18);
    g.Fit(&f,"QN");
    a=std::fabs(f.GetParameter(0)); b=std::fabs(f.GetParameter(1));
    double cn=f.GetChisquare()/std::max(1,f.GetNDF()); be=f.GetParError(1); if(cn>1) be*=std::sqrt(cn);
}

void paperSystematics(){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    const char* builds[]={"DSB1","TENERGY","MIXED","LUAG"};
    int col[4]={kAzure+2,kViolet+1,kOrange+8,kGreen+3};
    const double Es[]={25,50,75,100,125,150};

    // variation menu (applied to the per-energy sigma curve). idx 0 = NOMINAL.
    struct Var { const char* tag; double rScale; int K; bool useC; float mcpHi; float hgMin; double Emin; };
    // NOTE: the timing-method choice (cfd05 vs the recovered edge) is the paper's
    // thesis, quantified separately in Sec. method-gain -- NOT a selection
    // systematic, so it is deliberately excluded here. These probe robustness of
    // the result to reasonable analysis choices, all on the adopted source.
    std::vector<Var> V = {
        {"nominal",        1.00, 1000, false, 750.f, 20.f,  25},
        {"K=500",          1.00,  500, false, 750.f, 20.f,  25},
        {"K=2000",         1.00, 2000, false, 750.f, 20.f,  25},
        {"r_{fid}=2.5",    0.833,1000, false, 750.f, 20.f,  25},
        {"r_{fid}=3.5",    1.167,1000, false, 750.f, 20.f,  25},
        {"MCP<700",        1.00, 1000, false, 700.f, 20.f,  25},
        {"HG>30 mV",       1.00, 1000, false, 750.f, 30.f,  25},
        {"fit 50-150",     1.00, 1000, false, 750.f, 20.f,  50},
    };
    const int NV=(int)V.size();

    // results[build][var] -> {sigma150, a, b, be}
    double R150[4][32]={}, Ra[4][32]={}, Rb[4][32]={}, Rbe[4][32]={};

    printf("\n========== Paper 1 systematic budget (brightest-1000 (DW-UP)/2) ==========\n");
    for(int bi=0;bi<4;++bi){ const char* B=builds[bi]; int srcR=robustSrc(B);
        BuildConfig cfg=BuildConfig::Load(radConfig(B).Data()); if(!cfg.valid()){printf("%s: no config\n",B);continue;}
        printf("\n--- %s (nominal src=%s) ---\n",B,srcName(srcR));

        // sigma(E) for every variation, energy by energy (cache once per energy)
        std::vector<double> curveE[32], curveS[32], curveSe[32];
        bool validated=false;
        for(double e : Es){ TString path=radReduced(B,e);
            if(gSystem->AccessPathName(path.Data())) continue;
            std::vector<Ev> ev; double xc,yc; loadFile(path.Data(),cfg,srcR,ev,xc,yc);
            if(ev.empty()) continue;
            double rNom=TimingFiducialR(e);
            // one-time fidelity check vs the locked rad::timingBrightestK
            if(!validated){ TFile* fp=TFile::Open(path.Data()); TTree* t=(TTree*)fp->Get("rad"); RadView vv; vv.attach(t,&cfg);
                TimingResult tr=timingBrightestK(vv,e,srcR,1000);
                double mine=sigmaKnob(ev,xc,yc,rNom,1000,false,750.f,20.f);
                printf("   [fidelity @ %.0f GeV] mirror=%.2f ps  locked=%.2f ps  (delta %.2f)\n",
                       e,mine,tr.sigma_ps,mine-tr.sigma_ps);
                fp->Close(); validated=true; }
            for(int vi=0;vi<NV;++vi){ const Var& vv=V[vi];
                double s=sigmaKnob(ev,xc,yc,rNom*vv.rScale,vv.K,vv.useC,vv.mcpHi,vv.hgMin);
                if(s>0){ curveE[vi].push_back(e); curveS[vi].push_back(s); curveSe[vi].push_back(s/std::sqrt(2.0*vv.K)); }
            }
        }
        // fit each variation; record sigma@150 + (a,b)
        for(int vi=0;vi<NV;++vi){ double a,b,be; fitAB(curveE[vi],curveS[vi],curveSe[vi],V[vi].Emin,a,b,be);
            double s150=-1; for(size_t i=0;i<curveE[vi].size();++i) if(curveE[vi][i]>149) s150=curveS[vi][i];
            R150[bi][vi]=s150; Ra[bi][vi]=a; Rb[bi][vi]=b; Rbe[bi][vi]=be;
            printf("   %-12s  sigma150=%5.1f  a=%5.0f  b=%4.1f +- %.1f\n",V[vi].tag,s150,a,b,be);
        }
    }

    // ---- systematic = RMS of (variation - nominal) over the non-nominal knobs ----
    auto budget=[&](double R[4][32], int bi, double& nom, double& sys){
        nom=R[bi][0]; double q=0; int n=0;
        for(int vi=1;vi<NV;++vi){ if(R[bi][vi]<=0) continue; double d=R[bi][vi]-nom; q+=d*d; ++n; }
        sys=n? std::sqrt(q/n):0; };

    // ---------- LaTeX table ----------
    {std::ofstream o("papers/timing/tab_systematics.tex");
     o<<"% auto-generated by analyze/studies/paperSystematics.C\n";
     o<<"\\begin{table}[t]\n\\centering\n";
     o<<"\\caption{Systematic budget for the brightest-1000 $(\\mathrm{DW}-\\mathrm{UP})/2$ time\n";
     o<<"resolution. Each row recomputes $\\sigma_t(150)$ and the fit floor $b$ under one\n";
     o<<"selection variation; the entry is the shift from the nominal value (ps). The total\n";
     o<<"systematic is the RMS of the shifts; the statistical floor error is from the fit. The\n";
     o<<"timing-method choice (cfd05 vs.\\ the recovered edge) is the subject of\n";
     o<<"Sec.~\\ref{sec:method-gain} and is not folded in here.}\n";
     o<<"\\label{tab:syst}\n\\small\n\\begin{tabular}{lcccc}\n\\toprule\n";
     o<<"Variation & DSB1 & TENERGY & MIXED & LUAG \\\\\n";
     o<<"\\midrule\n";
     o<<"\\multicolumn{5}{l}{\\emph{shift in }$\\sigma_t(150)$\\emph{ [ps]}}\\\\\n";
     // clean LaTeX labels (table order == V order)
     const char* TL[8]={ "nominal",
        "$K{=}500$", "$K{=}2000$",
        "$r_{\\mathrm{fid}}{=}2.5$~mm", "$r_{\\mathrm{fid}}{=}3.5$~mm",
        "MCP${}<700$~mV", "HG${}>30$~mV", "fit $50$--$150$ only" };
     for(int vi=1;vi<NV;++vi){ o<<TL[vi];
        for(int bi=0;bi<4;++bi){ double d=(R150[bi][vi]>0&&R150[bi][0]>0)?R150[bi][vi]-R150[bi][0]:0;
            o<<Form(" & $%+.1f$",d); } o<<" \\\\\n"; }
     o<<"\\midrule\n";
     o<<"nominal $\\sigma_t(150)$";
     for(int bi=0;bi<4;++bi) o<<Form(" & $%.1f$",R150[bi][0]); o<<" \\\\\n";
     o<<"total syst.";
     for(int bi=0;bi<4;++bi){ double nom,sys; budget(R150,bi,nom,sys); o<<Form(" & $\\pm%.1f$",sys);} o<<" \\\\\n";
     o<<"\\midrule\n";
     o<<"floor $b$ (nom)";
     for(int bi=0;bi<4;++bi) o<<Form(" & $%.1f$",Rb[bi][0]); o<<" \\\\\n";
     o<<"$b$ syst.";
     for(int bi=0;bi<4;++bi){ double nom,sys; budget(Rb,bi,nom,sys); o<<Form(" & $\\pm%.1f$",sys);} o<<" \\\\\n";
     o<<"$b$ stat.";
     for(int bi=0;bi<4;++bi) o<<Form(" & $\\pm%.1f$",Rbe[bi][0]); o<<" \\\\\n";
     o<<"\\bottomrule\n\\end{tabular}\n\\end{table}\n";
     o.close(); printf("\n  wrote papers/timing/tab_systematics.tex\n"); }

    // ---------- figure: 2x2, per-build sigma150 stability under the variations ----------
    TCanvas* c=new TCanvas("syst","",1100,860); c->Divide(2,2,0.005,0.005);
    for(int bi=0;bi<4;++bi){ c->cd(bi+1); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.20); gPad->SetTopMargin(0.09); gPad->SetRightMargin(0.04);
        double nom,sys; budget(R150,bi,nom,sys);
        double bnom,bsys; budget(Rb,bi,bnom,bsys);
        double ymin=nom-std::max(4.0,3*sys), ymax=nom+std::max(4.0,3*sys);
        TH1F* fr=gPad->DrawFrame(-0.5,ymin,NV-0.5,ymax);
        fr->SetTitle(Form("%s;;brightest-1000  #sigma_{t}(150)  (ps)",builds[bi]));
        fr->GetYaxis()->SetTitleSize(0.05); fr->GetYaxis()->SetTitleOffset(1.25); fr->GetYaxis()->SetLabelSize(0.045);
        fr->GetXaxis()->SetNdivisions(0);
        // systematic band + nominal line
        TBox* bx=new TBox(-0.5,nom-sys,NV-0.5,nom+sys); bx->SetFillColorAlpha(col[bi],0.18); bx->SetLineWidth(0); bx->Draw();
        TLine* ln=new TLine(-0.5,nom,NV-0.5,nom); ln->SetLineColor(col[bi]); ln->SetLineWidth(2); ln->SetLineStyle(1); ln->Draw();
        // variation points
        for(int vi=0;vi<NV;++vi){ if(R150[bi][vi]<=0) continue;
            TGraph* g=new TGraph(1); g->SetPoint(0,vi,R150[bi][vi]);
            g->SetMarkerStyle(vi==0?29:20); g->SetMarkerSize(vi==0?2.4:1.4);
            g->SetMarkerColor(vi==0?kBlack:col[bi]); g->Draw("P SAME");
            TLatex tx; tx.SetTextAngle(90); tx.SetTextAlign(32); tx.SetTextSize(0.040); tx.SetTextColor(kGray+3);
            tx.DrawLatex(vi, ymin+0.04*(ymax-ymin), V[vi].tag); }
        TLatex t; t.SetNDC(); t.SetTextSize(0.052); t.SetTextColor(col[bi]);
        t.DrawLatex(0.55,0.84,Form("syst #pm%.1f ps",sys));
        t.SetTextColor(kBlack); t.SetTextSize(0.044);
        t.DrawLatex(0.55,0.78,Form("b=%.1f#pm%.1f(st)#pm%.1f(sy)",Rb[bi][0],Rbe[bi][0],bsys));
    }
    c->cd(0); TLatex tt; tt.SetNDC(); tt.SetTextFont(62); tt.SetTextSize(0.026);
    tt.DrawLatex(0.06,0.965,"Systematic stability of #sigma_{t}(150): nominal #pm band, each selection variation as a point");
    gSystem->mkdir("figures/narrative",kTRUE);
    c->Print("figures/narrative/systematics.png");
    printf("  wrote figures/narrative/systematics.png\n");
}

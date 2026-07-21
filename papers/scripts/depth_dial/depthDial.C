// ============================================================================
// depthDial.C — GATE 3: mean (t_DW - t_UP)/2 vs ln(E): the longitudinal depth
// observable. Pre-registered method in AUDIT.md (read it first):
//   primary = DSB1 (full 25-150 GeV scan), source = lgcfd (amplitude-independent
//   CFD; led/cfd05 as walk-sensitive cross-checks), FIXED r=2.5 mm fiducial at
//   all E, ALL timing ends required valid (fixed channel composition), full
//   fiducial (no brightness cut), robust truncated mean, Delta relative to the
//   lowest energy, fit p1*ln(E/Eref) in ps per e-fold.
// Prediction: p1 = -X0/v_g = -5.4 mm / 205 mm/ns = -26 ps per e-fold (PASS band
// -(26+-10), negative sign required; high-E compression allowed - filament edge).
//   source setup.sh
//   root -l -b -q 'papers/scripts/depth_dial/depthDial.C+'
// Output: papers/figures/depth_dial/depth_dial.png (+ _diag.png) + stdout table.
// ============================================================================
#include "RadView.h"
#include "DataPaths.h"
#include "PlotUtils.h"       // ApplyRADiCALStyle, DrawSuperTitle
#include "SelectionCuts.h"   // kMCP1_*, kHG_minPeak, kTimingChanConsistency_ns
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
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>
using namespace rad;

static const double kX0_mm   = 5.4;    // effective radiation length of the stack
static const double kVg_mmns = 205.0;  // quartz-core group velocity @495 nm (+-5)
static const double kPredSlope_ps = -kX0_mm / kVg_mmns * 1000.0;   // -26.3 ps / e-fold

// robust (2.5-sigma truncated, iterated) mean of t_diff; returns {mean, errOnMean, rms, n}
struct RM { double mu=0, err=0, rms=0; long n=0; };
// KERNEL CLONE: the iterative 2.5-sigma truncated core of rad::tebSigma
// (lib/physics/RadTiming.h), returning the MEAN (this gate reads a mean shift,
// not a width) — hence no Gaussian-core fit stage and no debias factor.
static RM robustMean(std::vector<float>& v){
    RM r; if(v.size()<100) return r;
    double mu=0; for(float x:v)mu+=x; mu/=v.size();
    double sd=0; for(float x:v)sd+=(x-mu)*(x-mu); sd=std::sqrt(sd/v.size()); if(sd<1e-4)sd=0.1;
    for(int it=0;it<6;++it){ double s=0,ss=0; long n=0;
        for(float x:v) if(std::fabs(x-mu)<2.5*sd){ s+=x; ss+=x*x; ++n; }
        if(n<100) break; mu=s/n; double w2=ss/n-mu*mu; if(w2>0) sd=std::sqrt(w2); r.n=n; }
    r.mu=mu; r.rms=sd; r.err= r.n>0 ? sd/std::sqrt((double)r.n) : 0;
    return r;
}

// gather per-event t_diff for (build, E, src): pre-registered selection (AUDIT.md)
static std::vector<float> gather(const char* build, double E, int src, BuildConfig& cfg){
    std::vector<float> out;
    TString p=radReduced(build,E); if(gSystem->AccessPathName(p.Data())) return out;
    TFile* f=TFile::Open(p.Data()); if(!f||f->IsZombie()){ if(f)f->Close(); return out; }
    TTree* t=(TTree*)f->Get("rad"); if(!t){ f->Close(); return out; }
    RadView v; v.attach(t,&cfg); if(!v.hasSrc(src)){ f->Close(); return out; }
    double xc,yc; v.beamCenter(xc,yc); const double r2=2.5*2.5;   // FIXED r=2.5 at ALL E
    // the timing-end set for this build (fixed composition requirement)
    std::vector<int> tend; for(int c=0;c<8;++c) if(v.is_timing(c)) tend.push_back(c);
    Long64_t N=v.entries(); out.reserve(N/4);
    for(Long64_t i=0;i<N;++i){ v.get(i);
        if(!v.wc_ok()||v.mcp1_peak()<kMCP1_minPeak||v.mcp1_peak()>kMCP1_maxPeak) continue;
        double dx=v.x_trk()-xc,dy=v.y_trk()-yc; if(dx*dx+dy*dy>=r2) continue;
        // all timing ends must be valid (peak + crossing) — fixed channel mix vs E
        float tt[8]; bool ok=true;
        for(int c:tend){ if(v.hg_peak(c)<kHG_minPeak){ok=false;break;}
            float tc=v.timeOf(c,src); if(tc<=-1e5f){ok=false;break;} tt[c]=tc; }
        if(!ok) continue;
        // in-event consistency veto (drop EVENT if any end is broken — keeps composition fixed)
        float med; { float a[8]; int m=0; for(int c:tend)a[m++]=tt[c];
            std::nth_element(a,a+m/2,a+m); med=a[m/2]; }
        for(int c:tend) if(std::fabs(tt[c]-med)>=kTimingChanConsistency_ns){ ok=false; break; }
        if(!ok) continue;
        double ds=0,us=0; int dn=0,un=0;
        for(int c:tend){ if(c<4){ds+=tt[c];++dn;} else {us+=tt[c];++un;} }
        if(dn<1||un<1) continue;
        out.push_back(0.5f*(float)(ds/dn-us/un));
    }
    f->Close(); return out;
}

struct Series { std::vector<double> E,M,Me,R; double slope=0, slopeErr=0, chi2=0; int ndf=0; double ref=0; };
static Series buildSeries(const char* build,int src){
    Series S; BuildConfig cfg=BuildConfig::Load(radConfig(build).Data());
    const double Es[6]={25,50,75,100,125,150};
    for(double e:Es){ std::vector<float> v=gather(build,e,src,cfg); RM r=robustMean(v);
        if(r.n<1000) continue;
        S.E.push_back(e); S.M.push_back(r.mu*1000.0); S.Me.push_back(r.err*1000.0); S.R.push_back(r.rms*1000.0);
        printf("  %-7s %-6s %3.0f GeV: N=%6ld  <t_diff>=%9.2f ps  err=%5.2f  RMS=%5.1f\n",
               build,RadView::srcName(src),e,r.n,r.mu*1000,r.err*1000,r.rms*1000);
    }
    if(S.E.size()>=3){ S.ref=S.M.front();
        std::vector<double> lx(S.E.size()); for(size_t i=0;i<S.E.size();++i) lx[i]=std::log(S.E[i]/S.E.front());
        TGraphErrors g(S.E.size(),&lx[0],&S.M[0],nullptr,&S.Me[0]);
        TF1 f("f","[0]+[1]*x"); f.SetParameters(S.M.front(),-25);
        g.Fit(&f,"QN");
        S.slope=f.GetParameter(1); S.slopeErr=f.GetParError(1);
        S.chi2=f.GetChisquare(); S.ndf=f.GetNDF();
        // inflate error for systematic scatter if chi2/ndf>1 (sub-ps stat errors)
        if(S.ndf>0 && S.chi2/S.ndf>1) S.slopeErr*=std::sqrt(S.chi2/S.ndf);
    }
    return S;
}

void depthDial(){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    gStyle->SetGridStyle(3); gStyle->SetGridColor(kGray);
    gSystem->mkdir("papers/figures/depth_dial",kTRUE);
    printf("\n===== GATE 3 DEPTH DIAL: mean (DW-UP)/2 vs ln(E) =====\n");
    printf("prediction: slope = -X0/v_g = %.1f ps per e-fold (PASS band -(26+-10))\n\n",kPredSlope_ps);

    // primary: DSB1 lgcfd; method cross-checks: DSB1 led + cfd05; build cross-checks (adopted src)
    Series dP = buildSeries("DSB1",  RadView::kLGCFD);
    Series dL = buildSeries("DSB1",  RadView::kLED);
    Series dC = buildSeries("DSB1",  RadView::kCFD05);
    Series tE = buildSeries("TENERGY",RadView::kLED);
    Series lG = buildSeries("LUAG",  RadView::kLED);

    printf("\n  ---- slopes (ps per e-fold of E) ----\n");
    auto rep=[&](const char* n,Series& S){ if(S.E.size()>=3)
        printf("  %-16s %+7.1f +- %4.1f   chi2/ndf=%.1f/%d\n",n,S.slope,S.slopeErr,S.chi2,S.ndf); };
    rep("DSB1 lgcfd",dP); rep("DSB1 led",dL); rep("DSB1 cfd05",dC); rep("TENERGY led",tE); rep("LUAG led",lG);

    // ---------- publication figure: DSB1 lgcfd, Delta vs E (log x) ----------
    TCanvas* c=new TCanvas("dd","",900,680);
    c->SetLeftMargin(0.13); c->SetRightMargin(0.05); c->SetTopMargin(0.04); c->SetBottomMargin(0.13);
    c->SetLogx(); c->SetGridy();
    TH1F* fr=c->DrawFrame(20,-70,200,25);
    fr->SetTitle(";beam energy E (GeV);#Delta#LT(t_{DW}#minust_{UP})/2#GT  (ps)");
    fr->GetYaxis()->SetTitleSize(0.045); fr->GetYaxis()->SetTitleOffset(1.3);
    fr->GetXaxis()->SetTitleSize(0.045); fr->GetXaxis()->SetMoreLogLabels(); fr->GetXaxis()->SetNoExponent();
    fr->GetXaxis()->SetLabelSize(0.033);   // keep the 90 and 100 log labels from colliding
    // prediction band: slope -26.3 +- (v_g 205+-5 & X0 uncertainty ~ +-2 ps/efold band)
    { double e0=dP.E.empty()?25:dP.E.front();
      TGraph* up=new TGraph(); TGraph* dn=new TGraph(); int k=0;
      for(double e=e0;e<=160;e*=1.05){ double l=std::log(e/e0);
          up->SetPoint(k,e,(kPredSlope_ps+2.5)*l); dn->SetPoint(k,e,(kPredSlope_ps-2.5)*l); ++k; }
      up->SetLineColor(kGray+1); up->SetLineStyle(2); up->SetLineWidth(2);
      dn->SetLineColor(kGray+1); dn->SetLineStyle(2); dn->SetLineWidth(2);
      up->Draw("L SAME"); dn->Draw("L SAME"); }
    if(dP.E.size()>=3){
        std::vector<double> dM(dP.E.size()),ze(dP.E.size(),0);
        for(size_t i=0;i<dP.E.size();++i) dM[i]=dP.M[i]-dP.ref;
        TGraphErrors* g=new TGraphErrors(dP.E.size(),&dP.E[0],&dM[0],&ze[0],&dP.Me[0]);
        g->SetMarkerStyle(20); g->SetMarkerColor(kAzure+2); g->SetLineColor(kAzure+2); g->SetMarkerSize(1.7);
        TF1* fit=new TF1("fitp",Form("[0]*log(x/%.1f)",dP.E.front()),dP.E.front(),155);
        fit->SetParameter(0,dP.slope); fit->SetLineColor(kAzure+2); fit->SetLineWidth(3);
        fit->Draw("SAME"); g->Draw("P SAME");
        TLegend* lg=new TLegend(0.30,0.70,0.84,0.92); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.034);
        lg->AddEntry(g,"DSB1, srCFD (lgcfd), full fiducial r=2.5 mm","lp");
        lg->AddEntry(fit,Form("fit: %+.1f #pm %.1f ps per e-fold",dP.slope,dP.slopeErr),"l");
        lg->AddEntry((TObject*)nullptr,Form("prediction #minusX_{0}/v_{g} = %.1f ps per e-fold (dashed)",kPredSlope_ps),"");
        lg->Draw();
    }
    // paper convention (format pass 2026-06-09): no internal super-title; the LaTeX caption carries it
    c->Print("papers/figures/depth_dial/depth_dial.png");

    // ---------- diagnostic figure: method + build cross-checks + RMS ----------
    TCanvas* d=new TCanvas("ddg","",1500,560); d->Divide(2,1,0.005,0.005);
    d->cd(1); gPad->SetLeftMargin(0.13);gPad->SetRightMargin(0.03);gPad->SetTopMargin(0.08);gPad->SetBottomMargin(0.14);gPad->SetLogx();gPad->SetGridy();
    TH1F* f1=gPad->DrawFrame(20,-80,200,30);
    f1->SetTitle("method & build cross-checks;beam energy E (GeV);#Delta#LTt_{diff}#GT (ps)");
    f1->GetXaxis()->SetMoreLogLabels(); f1->GetXaxis()->SetNoExponent();
    TLegend* l1=new TLegend(0.15,0.16,0.62,0.42); l1->SetBorderSize(0);l1->SetFillStyle(0);l1->SetTextSize(0.033);
    auto drawS=[&](Series& S,int col,int mk,const char* lab){ if(S.E.size()<3) return;
        std::vector<double> dM(S.E.size()),ze(S.E.size(),0);
        for(size_t i=0;i<S.E.size();++i) dM[i]=S.M[i]-S.ref;
        TGraphErrors* g=new TGraphErrors(S.E.size(),&S.E[0],&dM[0],&ze[0],&S.Me[0]);
        g->SetMarkerStyle(mk);g->SetMarkerColor(col);g->SetLineColor(col);g->SetMarkerSize(1.3);g->Draw("PL SAME");
        l1->AddEntry(g,Form("%s: %+.1f#pm%.1f",lab,S.slope,S.slopeErr),"lp"); };
    drawS(dP,kAzure+2,20,"DSB1 srCFD");
    drawS(dC,kBlack,24,"DSB1 cfd05");
    drawS(dL,kRed+1,25,"DSB1 led (walk-prone)");
    drawS(tE,kGreen+3,21,"TENERGY led");
    drawS(lG,kOrange+8,26,"LuAG led (#tau~70 ns)");
    l1->Draw();
    d->cd(2); gPad->SetLeftMargin(0.13);gPad->SetRightMargin(0.03);gPad->SetTopMargin(0.08);gPad->SetBottomMargin(0.14);gPad->SetLogx();gPad->SetGridy();
    TH1F* f2=gPad->DrawFrame(20,0,200,120);
    f2->SetTitle("event-by-event spread (the depth-fluctuation floor);beam energy E (GeV);RMS of t_{diff} (ps)");
    f2->GetXaxis()->SetMoreLogLabels(); f2->GetXaxis()->SetNoExponent();
    if(dP.E.size()>=3){ std::vector<double> ze(dP.E.size(),0);
        TGraphErrors* g=new TGraphErrors(dP.E.size(),&dP.E[0],&dP.R[0],&ze[0],&ze[0]);
        g->SetMarkerStyle(20);g->SetMarkerColor(kAzure+2);g->SetLineColor(kAzure+2);g->SetMarkerSize(1.4);g->Draw("PL SAME");
        TLatex tx; tx.SetNDC(); tx.SetTextSize(0.035); tx.SetTextColor(kGray+3);
        tx.DrawLatex(0.40,0.80,"full-fiducial RMS: photostat #oplus depth");
        tx.DrawLatex(0.40,0.74,"#sigma_{z}=v_{g}#upoint#sigma_{t}: 1 X_{0} #leftrightarrow 26 ps"); }
    d->cd(0); DrawSuperTitle("Depth dial diagnostics: method/build consistency (left), t_{diff} spread (right)",0.022f);
    d->Print("papers/figures/depth_dial/depth_dial_diag.png");

    // verdict
    bool pass = dP.E.size()>=4 && dP.slope<0 && std::fabs(dP.slope-kPredSlope_ps)<=10.0+2*dP.slopeErr;
    printf("\n  GATE 3 VERDICT: slope %+.1f +- %.1f ps/e-fold vs predicted %.1f  ->  %s\n",
           dP.slope,dP.slopeErr,kPredSlope_ps, pass?"PASS":"REVIEW (outside pre-registered band)");
    printf("  wrote papers/figures/depth_dial/depth_dial.png + depth_dial_diag.png\n");
}

// ============================================================================
// positionReconcile.C — GATE 2: what is the position number, really?
// Pre-registered protocol in AUDIT.md (this directory). Reproduces the
// showerLocalization.C estimator on the CURRENT schema with:
//   * train/test split (fit linear estimator on even events, evaluate on odd),
//   * unbinned held-out residual distributions (Gaussian core, RMS, robust sigma),
//   * binned closure (mean x_reco vs x_WC, slope, bin-mean residuals + errors),
//   * comparator decomposition with the explicit imaginary-unfolding check,
//   * fiducial variants (|x|,|y|<6 vs r<2.5) and energies (150 primary, 50, 25).
//   source setup.sh
//   root -l -b -q 'papers/scripts/position_reconciliation/positionReconcile.C+'
// Output: papers/figures/position_reconciliation/*.png + stdout log.
// ============================================================================
#include "RadView.h"
#include "DataPaths.h"
#include "PlotUtils.h"
#include "SelectionCuts.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TLinearFitter.h"
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

static const double kWC_sumUB = 3.6;   // mm — the (inapplicable) sum-side upper bound, for the check

struct Stats { double rms=0, gaus=0, gausErr=0, rob=0, mean=0; long n=0; };
static Stats stats(std::vector<float>& v,const char* hn){
    Stats s; if(v.size()<200) return s; s.n=v.size();
    double mu=0; for(float x:v)mu+=x; mu/=v.size(); s.mean=mu;
    double m2=0; for(float x:v)m2+=(x-mu)*(x-mu); s.rms=std::sqrt(m2/v.size());
    std::vector<float> a=v; std::sort(a.begin(),a.end());
    double q25=a[(size_t)(0.25*a.size())], q75=a[(size_t)(0.75*a.size())];
    s.rob=(q75-q25)/1.349;
    TH1F h(hn,"",160,mu-5*s.rms,mu+5*s.rms); h.SetDirectory(nullptr);
    for(float x:v)h.Fill(x);
    double m,mE,sg,sE; FitGaussCore(&h,2.0,m,mE,sg,sE); s.gaus=sg; s.gausErr=sE;
    return s;
}

struct Ev { float f[3]; float xw,yw; float r2; };

// gatherer (direct branch access; module centre from the LG-weighted beam centroid)
static std::vector<Ev> gather2(const char* build,double E){
    std::vector<Ev> out;
    TString p=radReduced(build,E); if(gSystem->AccessPathName(p.Data())) return out;
    TFile* f=TFile::Open(p.Data()); if(!f||f->IsZombie()){if(f)f->Close();return out;}
    TTree* t=(TTree*)f->Get("rad"); if(!t){f->Close();return out;}
    Bool_t wc; Float_t x,y,mp,lg[8];
    t->SetBranchStatus("*",0);
    for(const char* b:{"wc_ok","x_trk","y_trk","mcp1_peak","lg_peak"}) t->SetBranchStatus(b,1);
    t->SetBranchAddress("wc_ok",&wc); t->SetBranchAddress("x_trk",&x); t->SetBranchAddress("y_trk",&y);
    t->SetBranchAddress("mcp1_peak",&mp); t->SetBranchAddress("lg_peak",lg);
    // pass 1: LG-weighted beam centroid (selection-matched, replaces hardcoded 6.6/4.7)
    Long64_t N=t->GetEntries(); double sx=0,sy=0; long nc=0;
    for(Long64_t i=0;i<N;++i){ t->GetEntry(i);
        if(!wc||mp<kMCP1_minPeak||mp>kMCP1_maxPeak) continue;
        double s=0; for(int k=0;k<8;++k)s+=lg[k]; if(s<500) continue;
        sx+=x; sy+=y; ++nc; }
    if(nc<5000){ f->Close(); return out; }
    double cx0=sx/nc, cy0=sy/nc;
    // pass 2: events
    out.reserve(nc);
    for(Long64_t i=0;i<N;++i){ t->GetEntry(i);
        if(!wc||mp<kMCP1_minPeak||mp>kMCP1_maxPeak) continue;
        double A[4]={(double)lg[0]+lg[4],(double)lg[1]+lg[5],(double)lg[2]+lg[6],(double)lg[3]+lg[7]};
        double s=A[0]+A[1]+A[2]+A[3]; if(s<500) continue;
        double xw=x-cx0, yw=y-cy0; if(std::fabs(xw)>6||std::fabs(yw)>6) continue;
        Ev e; e.f[0]=A[0]/s; e.f[1]=A[1]/s; e.f[2]=A[2]/s; e.xw=xw; e.yw=yw; e.r2=xw*xw+yw*yw;
        out.push_back(e);
    }
    f->Close(); return out;
}

// fit linear estimator (const + 3 fractions) on TRAIN, return coefficients
static void fitLin(std::vector<Ev>& tr,bool useY,double c[4]){
    TLinearFitter lf(3,"hyp3");
    for(Ev& e:tr){ double ff[3]={e.f[0],e.f[1],e.f[2]}; lf.AddPoint(ff, useY? e.yw : e.xw, 1.0); }
    lf.Eval(); for(int k=0;k<4;++k)c[k]=lf.GetParameter(k);
}
static inline double predict(const double c[4],const Ev& e){ return c[0]+c[1]*e.f[0]+c[2]*e.f[1]+c[3]*e.f[2]; }

void positionReconcile(){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0); gStyle->SetGridStyle(3); gStyle->SetGridColor(kGray);
    gSystem->mkdir("papers/figures/position_reconciliation",kTRUE);
    printf("\n===== GATE 2: POSITION ARITHMETIC RECONCILIATION =====\n");
    printf("estimator: x_reco = c0 + c.f (3 corner light fractions), trained on EVEN events,\n");
    printf("evaluated on ODD events (held-out). All residuals are unbinned, event-level.\n");

    struct Cfg { const char* build; double E; const char* tag; } CFGS[3]={
        {"DSB1",150,"DSB1 150 GeV"},{"DSB1",50,"DSB1 50 GeV"},{"DSB1",25,"DSB1 25 GeV"}};
    std::vector<float> rx150, ry150;  double slope150=0; std::vector<Ev> ev150;
    for(auto& C:CFGS){
        std::vector<Ev> ev=gather2(C.build,C.E);
        if(ev.size()<20000){ printf("  %s: insufficient events (%zu)\n",C.tag,ev.size()); continue; }
        std::vector<Ev> tr,te; tr.reserve(ev.size()/2); te.reserve(ev.size()/2);
        for(size_t i=0;i<ev.size();++i) ((i&1)?te:tr).push_back(ev[i]);
        double cx[4],cy[4]; fitLin(tr,false,cx); fitLin(tr,true,cy);
        // residuals: train (for the overfit check) and held-out test
        std::vector<float> rxT,ryT,rxE,ryE,rxTight,ryTight;
        for(Ev& e:tr){ rxT.push_back(predict(cx,e)-e.xw); ryT.push_back(predict(cy,e)-e.yw); }
        for(Ev& e:te){ float dx=predict(cx,e)-e.xw, dy=predict(cy,e)-e.yw;
            rxE.push_back(dx); ryE.push_back(dy);
            if(e.r2<2.5*2.5){ rxTight.push_back(dx); ryTight.push_back(dy); } }
        Stats sxT=stats(rxT,"sxT"), sxE=stats(rxE,"sxE"), syE=stats(ryE,"syE");
        Stats sxTt=stats(rxTight,"sxTt"), syTt=stats(ryTight,"syTt");
        // beam-spot (no-information) baseline on the test sample
        std::vector<float> bx; for(Ev& e:te) bx.push_back(e.xw);
        Stats sb=stats(bx,"sb");
        printf("\n  --- %s (N=%zu; train %zu / test %zu) ---\n",C.tag,ev.size(),tr.size(),te.size());
        printf("  beam-spot spread (no-info baseline):           sigma_x = %.2f mm (RMS)\n",sb.rms);
        printf("  TRAIN residual x:  RMS=%.2f  (the historical in-sample quantity)\n",sxT.rms);
        printf("  HELD-OUT residual x: RMS=%.2f  Gaus-core=%.2f+-%.2f  robust=%.2f  (N=%ld)\n",
               sxE.rms,sxE.gaus,sxE.gausErr,sxE.rob,sxE.n);
        printf("  HELD-OUT residual y: RMS=%.2f  Gaus-core=%.2f+-%.2f  robust=%.2f\n",syE.rms,syE.gaus,syE.gausErr,syE.rob);
        printf("  tight fiducial r<2.5: x RMS=%.2f  y RMS=%.2f\n",sxTt.rms,syTt.rms);
        printf("  overfit check (train vs test RMS): %.3f vs %.3f mm  (delta=%.1f um)\n",
               sxT.rms,sxE.rms,1000.0*(sxE.rms-sxT.rms));
        // comparator decomposition (x, held-out RMS)
        double meas=sxE.rms;
        printf("  comparator decomposition (x): measured held-out RMS = %.2f mm\n",meas);
        double q=meas*meas-kWC_sumUB*kWC_sumUB;
        if(q<0) printf("    vs 3.6 mm sum-side UB: %.2f^2 - %.2f^2 < 0  -> UNFOLDING IMAGINARY/INVALID:\n"
                       "    the 3.6 mm sum-side bound CANNOT be the event-level comparator term\n"
                       "    (it is t0-inflated; t0 cancels in the position difference). The residual\n"
                       "    itself bounds the applicable term: sigma_WC(diff) <= %.2f mm.\n",meas,kWC_sumUB,meas);
        printf("    joint bound: sigma_capillary <= %.2f mm AND sigma_WC(diff) <= %.2f mm;\n"
               "    intrinsic unfolding NOT performed (no measured applicable comparator term).\n",meas,meas);
        if(C.E>100){ rx150=rxE; ry150=ryE; ev150=te;
            // closure slope on held-out: profile of x_reco vs x_WC
            TProfile pf("pf","",6,-6,6); for(Ev& e:te) pf.Fill(e.xw,predict(cx,e));
            TF1 lin("lin","[0]+[1]*x",-6,6); pf.Fit(&lin,"QN"); slope150=lin.GetParameter(1);
        }
    }

    // ---------- figures (150 GeV, held-out) ----------
    if(rx150.size()>10000){
        std::vector<Ev>& te=ev150;
        // recompute coefficients for drawing (train half again — same split)
        // (figure uses stored residuals; closure profile rebuilt from te + stored predictions via residual+truth)
        TCanvas* c=new TCanvas("pr","",1500,1050); c->Divide(2,2,0.006,0.010);
        // pad 1: unbinned residual x
        c->cd(1); gPad->SetLeftMargin(0.13);gPad->SetBottomMargin(0.12);gPad->SetTopMargin(0.08);gPad->SetGridy();
        Stats sx=stats(rx150,"fsx");
        TH1F* hx=new TH1F("hxr",";x_{reco} #minus x_{WC}  (mm);events",160,-8,8); hx->SetDirectory(nullptr);
        for(float r:rx150)hx->Fill(r);
        hx->SetLineColor(kAzure+2); hx->SetFillColorAlpha(kAzure+2,0.30); hx->Draw("HIST");
        TF1* g=new TF1("gg","gaus",-8,8); g->SetParameters(hx->GetMaximum(),sx.mean,sx.gaus);
        g->SetLineColor(kRed+1); g->Draw("SAME");
        TLatex tx; tx.SetNDC(); tx.SetTextSize(0.038);
        tx.DrawLatex(0.16,0.86,Form("HELD-OUT: RMS=%.2f, core #sigma=%.2f, robust=%.2f mm",sx.rms,sx.gaus,sx.rob));
        tx.DrawLatex(0.16,0.80,Form("N=%ld",sx.n));
        DrawPadTitle("unbinned event-level residual, x (DSB1 150 GeV)");
        // pad 2: y
        c->cd(2); gPad->SetLeftMargin(0.13);gPad->SetBottomMargin(0.12);gPad->SetTopMargin(0.08);gPad->SetGridy();
        Stats sy=stats(ry150,"fsy");
        TH1F* hy=new TH1F("hyr",";y_{reco} #minus y_{WC}  (mm);events",160,-8,8); hy->SetDirectory(nullptr);
        for(float r:ry150)hy->Fill(r);
        hy->SetLineColor(kGreen+3); hy->SetFillColorAlpha(kGreen+3,0.30); hy->Draw("HIST");
        TF1* g2=new TF1("gg2","gaus",-8,8); g2->SetParameters(hy->GetMaximum(),sy.mean,sy.gaus);
        g2->SetLineColor(kRed+1); g2->Draw("SAME");
        tx.DrawLatex(0.16,0.86,Form("HELD-OUT: RMS=%.2f, core #sigma=%.2f, robust=%.2f mm",sy.rms,sy.gaus,sy.rob));
        DrawPadTitle("unbinned event-level residual, y (DSB1 150 GeV)");
        // pad 3: binned closure (rebuild from te using residual+truth)
        c->cd(3); gPad->SetLeftMargin(0.13);gPad->SetBottomMargin(0.12);gPad->SetTopMargin(0.08);gPad->SetGridy();
        TProfile* pf=new TProfile("pfc",";x_{WC} (mm);#LT x_{reco} #GT (mm)",12,-6,6); pf->SetDirectory(nullptr);
        for(size_t i=0;i<te.size();++i) pf->Fill(te[i].xw, te[i].xw + rx150[i]);
        pf->SetMarkerStyle(20); pf->SetMarkerColor(kAzure+2); pf->SetLineColor(kAzure+2); pf->SetMarkerSize(1.3);
        pf->SetMinimum(-6); pf->SetMaximum(6); pf->Draw();
        TLine* dg=new TLine(-6,-6,6,6); dg->SetLineStyle(2); dg->SetLineColor(kGray+2); dg->Draw();
        tx.DrawLatex(0.16,0.86,Form("closure slope (held-out) = %.3f (ideal 1)",slope150));
        tx.DrawLatex(0.16,0.80,"bin-mean errors #approx RMS/#sqrt{N_{bin}} #ll mm (NOT a resolution)");
        DrawPadTitle("binned closure: #LTx_{reco}#GT vs x_{WC}");
        // pad 4: decomposition summary text
        c->cd(4); gPad->SetTopMargin(0.08);
        TLatex t4; t4.SetNDC(); t4.SetTextSize(0.040); double yy=0.86;
        auto L=[&](const char* s){ t4.DrawLatex(0.06,yy,s); yy-=0.075; };
        L("comparator decomposition (x, held-out RMS):");
        L(Form("measured residual width: %.2f mm",sx.rms));
        L(Form("#sigma_{WC}(sum-side UB) = %.1f mm: INAPPLICABLE (t_{0}-inflated;",kWC_sumUB));
        L("   t_{0} cancels in the position difference)");
        L(Form("unfolding vs 3.6: %.1f^{2}#minus%.1f^{2} < 0 #Rightarrow imaginary, REFUSED",sx.rms,kWC_sumUB));
        L(Form("joint bound: #sigma_{cap} #leq %.2f mm AND #sigma_{WC}(diff) #leq %.2f mm",sx.rms,sx.rms));
        L("no intrinsic resolution is quoted from this data");
        DrawPadTitle("what may be claimed");
        c->cd(0); DrawSuperTitle("GATE 2 reconciliation: held-out event-level residuals, closure, and the comparator decomposition (DSB1 150 GeV)",0.019f);
        c->Print("papers/figures/position_reconciliation/position_reconcile.png");
        printf("\n  wrote papers/figures/position_reconciliation/position_reconcile.png\n");
    }
}

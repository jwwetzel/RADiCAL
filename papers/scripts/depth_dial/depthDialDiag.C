// ============================================================================
// depthDialDiag.C — GATE 3 hardening: stability diagnostics for the depth dial
// (companion to depthDial.C; see AUDIT.md + DEPTH_DIAL_REVIEW.md).
// Single pass per DSB1 energy file collects per-event:
//   t_diff for FIVE sources (lgcfd, cfd05, led, cfd20, cfd30) each with its own
//   all-ends-valid + median-veto flag (replicates depthDial.C per-source logic),
//   per-capillary t_diff (lgcfd, 4 capillaries), r^2 to beam centroid, sum_lg, run.
// Variants computed post-hoc from the same event store:
//   V1 fiducial: r<2.5 (baseline) vs r<1.5 (tight)
//   V2 subset fits: all 6 E / drop 150 / drop {125,150}   (baseline events)
//   V3 per-capillary slopes (NW/NE/SE/SW independently)
//   V4 brightness quartiles (sum_lg quartiles per energy)
//   V5 per-run means at every E (table; spread at 125/150 re-steepening check)
//   V6 method table: slopes for all five sources
//   source setup.sh
//   root -l -b -q 'papers/scripts/depth_dial/depthDialDiag.C+'
// Output: papers/figures/depth_dial/diagnostics/depth_dial_stability.png + stdout log.
// ============================================================================
#include "RadView.h"
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
#include "TStyle.h"
#include "TSystem.h"
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include <cstdio>
using namespace rad;

static const int    NSRC=5;
static const int    SRC [NSRC]={RadView::kLGCFD,RadView::kCFD05,RadView::kLED,RadView::kCFD20,RadView::kCFD30};
static const char*  SNM [NSRC]={"lgcfd","cfd05","led","cfd20","cfd30"};
static const double kPred=-26.3;

struct Ev { float r2, slg, td[NSRC], cap[4]; bool ok[NSRC]; int run; };
struct RM { double mu=0,err=0,rms=0; long n=0; };

static RM rmean(std::vector<float>& v){
    RM r; if(v.size()<100) return r;
    double mu=0; for(float x:v)mu+=x; mu/=v.size();
    double sd=0; for(float x:v)sd+=(x-mu)*(x-mu); sd=std::sqrt(sd/v.size()); if(sd<1e-4)sd=0.1;
    for(int it=0;it<6;++it){ double s=0,ss=0; long n=0;
        for(float x:v) if(std::fabs(x-mu)<2.5*sd){s+=x;ss+=x*x;++n;}
        if(n<100)break; mu=s/n; double w2=ss/n-mu*mu; if(w2>0)sd=std::sqrt(w2); r.n=n; }
    r.mu=mu; r.rms=sd; r.err=r.n>0?sd/std::sqrt((double)r.n):0; return r;
}
// weighted ln(E) fit -> slope (ps/e-fold), err inflated by sqrt(chi2/ndf)
struct Fit { double s=0,se=0,chi=0; int ndf=0; };
static Fit lfit(std::vector<double>& E,std::vector<double>& M,std::vector<double>& Me){
    Fit F; if(E.size()<3) return F;
    std::vector<double> lx(E.size()); for(size_t i=0;i<E.size();++i) lx[i]=std::log(E[i]/E[0]);
    TGraphErrors g(E.size(),&lx[0],&M[0],nullptr,&Me[0]);
    TF1 f("lf","[0]+[1]*x"); f.SetParameters(M[0],-25); g.Fit(&f,"QN");
    F.s=f.GetParameter(1); F.se=f.GetParError(1); F.chi=f.GetChisquare(); F.ndf=f.GetNDF();
    if(F.ndf>0&&F.chi/F.ndf>1) F.se*=std::sqrt(F.chi/F.ndf);
    return F;
}

static std::vector<Ev> gatherAll(double E, BuildConfig& cfg){
    std::vector<Ev> out;
    TString p=radReduced("DSB1",E); if(gSystem->AccessPathName(p.Data())) return out;
    TFile* f=TFile::Open(p.Data()); if(!f||f->IsZombie()){if(f)f->Close();return out;}
    TTree* t=(TTree*)f->Get("rad"); if(!t){f->Close();return out;}
    RadView v; v.attach(t,&cfg); double xc,yc; v.beamCenter(xc,yc);
    std::vector<int> tend; for(int c=0;c<8;++c) if(v.is_timing(c)) tend.push_back(c);
    Long64_t N=v.entries(); out.reserve(N/4);
    for(Long64_t i=0;i<N;++i){ v.get(i);
        if(!v.wc_ok()||v.mcp1_peak()<kMCP1_minPeak||v.mcp1_peak()>kMCP1_maxPeak) continue;
        double dx=v.x_trk()-xc,dy=v.y_trk()-yc; double r2=dx*dx+dy*dy; if(r2>=2.5*2.5) continue;
        Ev e; e.r2=(float)r2; e.slg=(float)v.sum_lg(); e.run=v.run();
        bool peaksOK=true; for(int c:tend) if(v.hg_peak(c)<kHG_minPeak){peaksOK=false;break;}
        if(!peaksOK) continue;
        for(int s=0;s<NSRC;++s){ e.ok[s]=false; float tt[8]; bool good=true;
            for(int c:tend){ float tc=v.timeOf(c,SRC[s]); if(tc<=-1e5f){good=false;break;} tt[c]=tc; }
            if(!good) continue;
            float a[8]; int m=0; for(int c:tend)a[m++]=tt[c];
            std::nth_element(a,a+m/2,a+m); float med=a[m/2];
            for(int c:tend) if(std::fabs(tt[c]-med)>=kTimingChanConsistency_ns){good=false;break;}
            if(!good) continue;
            double ds=0,us=0;int dn=0,un=0;
            for(int c:tend){ if(c<4){ds+=tt[c];++dn;} else {us+=tt[c];++un;} }
            if(dn<1||un<1) continue;
            e.td[s]=0.5f*(float)(ds/dn-us/un); e.ok[s]=true;
            if(s==0) for(int k=0;k<4;++k) e.cap[k]=0.5f*(tt[k]-tt[k+4]);   // per-capillary, lgcfd
        }
        if(e.ok[0]||e.ok[1]||e.ok[2]) out.push_back(e);
    }
    f->Close(); return out;
}

void depthDialDiag(){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0); gStyle->SetGridStyle(3); gStyle->SetGridColor(kGray);
    gSystem->mkdir("papers/figures/depth_dial/diagnostics",kTRUE);
    BuildConfig cfg=BuildConfig::Load(radConfig("DSB1").Data());
    const double Es[6]={25,50,75,100,125,150};
    std::vector<std::vector<Ev>> ev(6);
    printf("\n===== DEPTH-DIAL STABILITY DIAGNOSTICS (DSB1) =====\n");
    for(int i=0;i<6;++i){ ev[i]=gatherAll(Es[i],cfg); printf("  %3.0f GeV: %zu baseline events (r<2.5, peaks ok)\n",Es[i],ev[i].size()); }

    // ---------- V1 + V2: radius variants & subset fits (lgcfd) ----------
    auto series=[&](double rmax,int qlo,int qhi)->std::vector<RM>{   // brightness quartile range [qlo,qhi] of 4; 0-3=all
        std::vector<RM> out(6);
        for(int i=0;i<6;++i){ std::vector<Ev>& V=ev[i]; if(V.empty())continue;
            std::vector<float> sl; for(Ev&e:V) if(e.ok[0]&&e.r2<rmax*rmax) sl.push_back(e.slg);
            if(sl.size()<400) continue;
            std::sort(sl.begin(),sl.end());
            float lo=sl[(size_t)(qlo*sl.size()/4)], hi=(qhi>=3)?1e12f:sl[(size_t)((qhi+1)*sl.size()/4)];
            std::vector<float> td; for(Ev&e:V) if(e.ok[0]&&e.r2<rmax*rmax&&e.slg>=lo&&e.slg<hi) td.push_back(e.td[0]);
            out[i]=rmean(td); }
        return out;
    };
    auto toFit=[&](std::vector<RM>& R,int nE)->Fit{
        std::vector<double> E,M,Me;
        for(int i=0;i<nE;++i) if(R[i].n>500){E.push_back(Es[i]);M.push_back(R[i].mu*1000);Me.push_back(R[i].err*1000);}
        return lfit(E,M,Me);
    };
    std::vector<RM> base=series(2.5,0,3), tight=series(1.5,0,3);
    Fit fAll=toFit(base,6), fD150=toFit(base,5), fD125=toFit(base,4), fTight=toFit(tight,6);
    printf("\n  V1/V2 fiducial + subsets (lgcfd, ps per e-fold; pred %.1f):\n",kPred);
    printf("    r<2.5 all 6 E    : %+7.1f +- %4.1f  (chi2/ndf %.0f/%d)\n",fAll.s,fAll.se,fAll.chi,fAll.ndf);
    printf("    r<2.5 drop 150   : %+7.1f +- %4.1f  (chi2/ndf %.0f/%d)\n",fD150.s,fD150.se,fD150.chi,fD150.ndf);
    printf("    r<2.5 drop125+150: %+7.1f +- %4.1f  (chi2/ndf %.0f/%d)\n",fD125.s,fD125.se,fD125.chi,fD125.ndf);
    printf("    r<1.5 all 6 E    : %+7.1f +- %4.1f  (chi2/ndf %.0f/%d)\n",fTight.s,fTight.se,fTight.chi,fTight.ndf);

    // ---------- V3: per-capillary ----------
    const char* CNM[4]={"NW","NE","SE","SW"}; Fit fCap[4]; std::vector<RM> capR[4];
    for(int k=0;k<4;++k){ capR[k].resize(6);
        for(int i=0;i<6;++i){ std::vector<float> td; for(Ev&e:ev[i]) if(e.ok[0]) td.push_back(e.cap[k]); capR[k][i]=rmean(td); }
        fCap[k]=toFit(capR[k],6); }
    printf("\n  V3 per-capillary slopes (lgcfd):\n");
    for(int k=0;k<4;++k) printf("    %s: %+7.1f +- %4.1f\n",CNM[k],fCap[k].s,fCap[k].se);

    // ---------- V4: brightness quartiles ----------
    Fit fQ[4];
    for(int q=0;q<4;++q){ std::vector<RM> R=series(2.5,q,q); fQ[q]=toFit(R,6); }
    printf("\n  V4 brightness-quartile slopes (lgcfd, Q1=dimmest):\n");
    for(int q=0;q<4;++q) printf("    Q%d: %+7.1f +- %4.1f\n",q+1,fQ[q].s,fQ[q].se);

    // ---------- V5: per-run means ----------
    printf("\n  V5 per-run robust means (lgcfd, runs with N>700):\n");
    std::vector<double> runX125,runY125,runX150,runY150;
    for(int i=0;i<6;++i){ std::map<int,std::vector<float>> byrun;
        for(Ev&e:ev[i]) if(e.ok[0]) byrun[e.run].push_back(e.td[0]);
        double mn=1e9,mx=-1e9; int nr=0;
        for(auto& kv:byrun){ if(kv.second.size()<700) continue; RM r=rmean(kv.second);
            if(r.n<500) continue;
            ++nr; mn=std::min(mn,r.mu*1000); mx=std::max(mx,r.mu*1000);
            if(i==4){runX125.push_back(kv.first);runY125.push_back(r.mu*1000);}
            if(i==5){runX150.push_back(kv.first);runY150.push_back(r.mu*1000);} }
        if(nr>1) printf("    %3.0f GeV: %2d runs with N>700, run-to-run spread (max-min) = %5.1f ps\n",Es[i],nr,mx-mn);
    }

    // ---------- V6: method table ----------
    printf("\n  V6 method slopes (all 6 E, own-validity samples):\n");
    Fit fS[NSRC];
    for(int s=0;s<NSRC;++s){ std::vector<RM> R(6);
        for(int i=0;i<6;++i){ std::vector<float> td; for(Ev&e:ev[i]) if(e.ok[s]) td.push_back(e.td[s]); R[i]=rmean(td); }
        fS[s]=toFit(R,6); printf("    %-6s: %+7.1f +- %4.1f  (chi2/ndf %.0f/%d)\n",SNM[s],fS[s].s,fS[s].se,fS[s].chi,fS[s].ndf); }

    // ---------- figure: 2x2 diagnostics ----------
    TCanvas* c=new TCanvas("ddd","",1500,1100); c->Divide(2,2,0.006,0.010);
    auto frame=[&](int pad,const char* title,double ylo,double yhi){ c->cd(pad);
        gPad->SetLeftMargin(0.13);gPad->SetRightMargin(0.04);gPad->SetTopMargin(0.09);gPad->SetBottomMargin(0.13);
        gPad->SetLogx();gPad->SetGridy();
        TH1F* fr=gPad->DrawFrame(20,ylo,200,yhi);
        fr->SetTitle(Form("%s;beam energy E (GeV);#Delta#LTt_{diff}#GT (ps)",title));
        fr->GetXaxis()->SetMoreLogLabels(); fr->GetXaxis()->SetNoExponent(); };
    auto drawRM=[&](std::vector<RM>& R,int col,int mk,TLegend* lg,const char* lab,Fit& F){
        std::vector<double> E,M,Me; for(int i=0;i<6;++i) if(R[i].n>500){E.push_back(Es[i]);M.push_back(R[i].mu*1000);Me.push_back(R[i].err*1000);}
        if(E.size()<3) return; double ref=M[0];
        for(auto&m:M)m-=ref; std::vector<double> ze(E.size(),0);
        TGraphErrors* g=new TGraphErrors(E.size(),&E[0],&M[0],&ze[0],&Me[0]);
        g->SetMarkerStyle(mk);g->SetMarkerColor(col);g->SetLineColor(col);g->SetMarkerSize(1.3);g->Draw("PL SAME");
        if(lg) lg->AddEntry(g,Form("%s: %+.1f#pm%.1f",lab,F.s,F.se),"lp"); };
    // pad 1: radius + prediction
    frame(1,"V1 fiducial radius",-85,25);
    { TLegend* l=new TLegend(0.15,0.16,0.62,0.34); l->SetBorderSize(0);l->SetFillStyle(0);l->SetTextSize(0.034);
      drawRM(base,kAzure+2,20,l,"r<2.5 (baseline)",fAll); drawRM(tight,kRed+1,24,l,"r<1.5 (tight)",fTight);
      TGraph* pr=new TGraph(); int k=0; for(double e=25;e<=160;e*=1.06) pr->SetPoint(k++,e,kPred*std::log(e/25));
      pr->SetLineColor(kGray+1);pr->SetLineStyle(2);pr->SetLineWidth(2);pr->Draw("L SAME");
      l->AddEntry(pr,"prediction -26.3","l"); l->Draw(); }
    // pad 2: per-capillary
    frame(2,"V3 per-capillary (lgcfd)",-85,25);
    { TLegend* l=new TLegend(0.15,0.16,0.62,0.40); l->SetBorderSize(0);l->SetFillStyle(0);l->SetTextSize(0.034);
      int cc[4]={kAzure+2,kGreen+3,kOrange+8,kMagenta+1}; int mm[4]={20,21,22,23};
      for(int k=0;k<4;++k) drawRM(capR[k],cc[k],mm[k],l,CNM[k],fCap[k]); l->Draw(); }
    // pad 3: brightness quartiles
    frame(3,"V4 brightness quartiles (lgcfd)",-85,25);
    { TLegend* l=new TLegend(0.15,0.16,0.62,0.40); l->SetBorderSize(0);l->SetFillStyle(0);l->SetTextSize(0.034);
      int cc[4]={kBlue+1,kCyan+2,kOrange+7,kRed+1};
      for(int q=0;q<4;++q){ std::vector<RM> R=series(2.5,q,q); drawRM(R,cc[q],20+q,l,Form("Q%d",q+1),fQ[q]); } l->Draw(); }
    // pad 4: per-run at 125/150
    c->cd(4); gPad->SetLeftMargin(0.13);gPad->SetRightMargin(0.04);gPad->SetTopMargin(0.09);gPad->SetBottomMargin(0.13);gPad->SetGridy();
    { double xmn=1e9,xmx=-1e9,ymn=1e9,ymx=-1e9;
      for(size_t i=0;i<runX125.size();++i){xmn=std::min(xmn,runX125[i]);xmx=std::max(xmx,runX125[i]);ymn=std::min(ymn,runY125[i]);ymx=std::max(ymx,runY125[i]);}
      for(size_t i=0;i<runX150.size();++i){xmn=std::min(xmn,runX150[i]);xmx=std::max(xmx,runX150[i]);ymn=std::min(ymn,runY150[i]);ymx=std::max(ymx,runY150[i]);}
      if(xmn>xmx){xmn=0;xmx=1;ymn=-300;ymx=-200;}
      TH1F* fr=gPad->DrawFrame(xmn-3,ymn-8,xmx+3,ymx+12);
      fr->SetTitle("V5 per-run means, 125 & 150 GeV;run number;#LTt_{diff}#GT (ps)");
      TGraph* g1=new TGraph(runX125.size(),&runX125[0],&runY125[0]); g1->SetMarkerStyle(20);g1->SetMarkerColor(kOrange+8);g1->SetMarkerSize(1.4);g1->Draw("P SAME");
      TGraph* g2=new TGraph(runX150.size(),&runX150[0],&runY150[0]); g2->SetMarkerStyle(21);g2->SetMarkerColor(kRed+1);g2->SetMarkerSize(1.4);g2->Draw("P SAME");
      TLegend* l=new TLegend(0.65,0.74,0.95,0.90); l->SetBorderSize(0);l->SetFillStyle(0);l->SetTextSize(0.036);
      l->AddEntry(g1,"125 GeV","p"); l->AddEntry(g2,"150 GeV","p"); l->Draw(); }
    c->cd(0); DrawSuperTitle("Depth-dial stability diagnostics (DSB1, srCFD unless noted): fiducial, per-capillary, brightness, run-period",0.020f);
    c->Print("papers/figures/depth_dial/diagnostics/depth_dial_stability.png");
    printf("\n  wrote papers/figures/depth_dial/diagnostics/depth_dial_stability.png\n");
}

// ============================================================================
// mixedKillshotBootstrap.C — GATE 6: harden the MIXED in-event DSB1-vs-LuAG
// kill-shot. Protocol pre-registered in AUDIT.md (this directory).
//   M1 direct event-paired width ratio (Down layer: DSB1 slots 1,3 vs LuAG 0,2)
//      on a COMMON all-7-ends-valid event set; fixed-window truncated RMS.
//   M2 pairwise-solve reproduction with GATE-1 labels (original quantity).
//   M3 Poisson event bootstrap (1000 replicas): median, 68%, 95% CI.
//   M4 run jackknife (leave-one-run-out).
//   M5 method dependence: srCFD (primary), cfd05 (original), LED.
//   M6 wrong-map diagnostics: swapped + scrambled, with the amplitude control.
// Corner map (GATE 1, pulse-shape confirmed): DSB1 = NE,SW (slots 1,3,5);
// LuAG = NW,SE (slots 0,2,4,6). SW-U excluded (cross-group).
//   source setup.sh
//   root -l -b -q 'papers/scripts/mixed_killshot_bootstrap/mixedKillshotBootstrap.C+'
// ============================================================================
#include "RadView.h"
#include "DataPaths.h"
#include "PlotUtils.h"
#include "SelectionCuts.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TRandom3.h"
#include <vector>
#include <map>
#include <algorithm>
#include <random>
#include <cmath>
#include <cstdio>
using namespace rad;

static const int NSRC=3;
static const char* SNM[NSRC]={"srCFD","cfd05","led"};
static const char* HBR[NSRC]={"hg_lgcfd","hg_cfd05","hg_led"};

struct Ev { float t[7]; float amp[7]; int run; };
struct Win { double lo,hi; };

// full-sample 2.5-sigma core window (frozen -> enables Poisson bootstrap)
// KERNEL CLONE: the iterative 2.5-sigma truncated window of rad::tebSigma
// (lib/physics/RadTiming.h); window FROZEN per dataset before bootstrap resampling.
static Win coreWindow(std::vector<float>& v){
    Win w{0,0}; if(v.size()<300) return w;
    double mu=0; for(float x:v)mu+=x; mu/=v.size();
    double sd=0; for(float x:v)sd+=(x-mu)*(x-mu); sd=std::sqrt(sd/v.size()); if(sd<1e-4)sd=0.1;
    for(int it=0;it<6;++it){ double s=0,ss=0; long n=0;
        for(float x:v) if(std::fabs(x-mu)<2.5*sd){s+=x;ss+=x*x;++n;}
        if(n<100)break; mu=s/n; double v2=ss/n-mu*mu; if(v2>0)sd=std::sqrt(v2); }
    w.lo=mu-2.5*sd; w.hi=mu+2.5*sd; return w;
}
// truncated RMS within a FIXED window, with optional per-event weights
static double winRMS(std::vector<float>& v, Win w, std::vector<int>* wt=nullptr){
    double s=0,ss=0,n=0;
    for(size_t i=0;i<v.size();++i){ float x=v[i]; if(x<w.lo||x>w.hi) continue;
        double k=wt?(*wt)[i]:1.0; if(k<=0)continue; s+=k*x; ss+=k*x*x; n+=k; }
    if(n<50) return -1; double mu=s/n, var=ss/n-mu*mu; return var>0?std::sqrt(var):-1;
}

static std::vector<Ev> gather(double E,int src){
    std::vector<Ev> out;
    TString p=radReduced("MIXED",E); if(gSystem->AccessPathName(p.Data())) return out;
    TFile* f=TFile::Open(p.Data()); if(!f||f->IsZombie()){if(f)f->Close();return out;}
    TTree* t=(TTree*)f->Get("rad"); if(!t){f->Close();return out;}
    Bool_t wc; Float_t x,y,mp,hp[8],ht[8]; Int_t run;
    t->SetBranchStatus("*",0);
    for(const char* b:{"wc_ok","x_trk","y_trk","mcp1_peak","hg_peak","run"}) t->SetBranchStatus(b,1);
    t->SetBranchStatus(HBR[src],1);
    t->SetBranchAddress("wc_ok",&wc); t->SetBranchAddress("x_trk",&x); t->SetBranchAddress("y_trk",&y);
    t->SetBranchAddress("mcp1_peak",&mp); t->SetBranchAddress("hg_peak",hp);
    t->SetBranchAddress(HBR[src],ht); t->SetBranchAddress("run",&run);
    Long64_t N=t->GetEntries(); double sx=0,sy=0; long nc=0;
    for(Long64_t i=0;i<N&&nc<40000;++i){ t->GetEntry(i); if(wc&&x>-100&&x<100){sx+=x;sy+=y;++nc;} }
    if(nc<1000){f->Close();return out;}
    double xc=sx/nc, yc=sy/nc;
    out.reserve(N/4);
    for(Long64_t i=0;i<N;++i){ t->GetEntry(i);
        if(!wc||mp<kMCP1_minPeak||mp>kMCP1_maxPeak) continue;
        double dx=x-xc,dy=y-yc; if(dx*dx+dy*dy>=9.0) continue;       // r<3 (matches original)
        bool ok=true; Ev e; e.run=run;
        for(int s=0;s<7;++s){ if(hp[s]<20||ht[s]<=-1e5f){ok=false;break;} e.t[s]=ht[s]; e.amp[s]=hp[s]; }
        if(ok) out.push_back(e);                                      // COMMON all-7-valid event set
    }
    f->Close(); return out;
}

// direct paired ratio on an event list: R = winRMS(t1-t3)/winRMS(t0-t2) with frozen windows
struct Direct { std::vector<float> dD,dL; Win wD,wL; double sD=-1,sL=-1,R=-1; };
static Direct direct(std::vector<Ev>& ev,int a1,int a3,int b0,int b2){
    Direct d; d.dD.reserve(ev.size()); d.dL.reserve(ev.size());
    for(Ev& e:ev){ d.dD.push_back(e.t[a1]-e.t[a3]); d.dL.push_back(e.t[b0]-e.t[b2]); }
    d.wD=coreWindow(d.dD); d.wL=coreWindow(d.dL);
    d.sD=winRMS(d.dD,d.wD); d.sL=winRMS(d.dL,d.wL);
    if(d.sD>0&&d.sL>0) d.R=d.sD/d.sL;
    return d;
}

// pairwise solve (original quantity) with GATE-1 labels; returns group ratio D/L
static double solveRatio(std::vector<Ev>& ev){
    static const bool isD[7]={false,true,false,true,false,true,false};
    auto pairSig=[&](int i,int j)->double{ std::vector<float> v; v.reserve(ev.size());
        for(Ev& e:ev) v.push_back(e.t[i]-e.t[j]); Win w=coreWindow(v); return winRMS(v,w); };
    double xd[4]; { double sij[4][4]={{0}};
        for(int a=0;a<4;++a)for(int b=a+1;b<4;++b){ double s=pairSig(a,b); sij[a][b]=sij[b][a]=s>0?s*s:0; }
        double R[4],S=0; for(int a=0;a<4;++a){R[a]=0;for(int b=0;b<4;++b)if(b!=a)R[a]+=sij[a][b];S+=R[a];} S/=6.0;
        for(int a=0;a<4;++a) xd[a]=(R[a]-S)/2.0; }
    double xu[3]; { double sij[3][3]={{0}};
        for(int a=0;a<3;++a)for(int b=a+1;b<3;++b){ double s=pairSig(4+a,4+b); sij[a][b]=sij[b][a]=s>0?s*s:0; }
        double R[3],S=0; for(int a=0;a<3;++a){R[a]=0;for(int b=0;b<3;++b)if(b!=a)R[a]+=sij[a][b];S+=R[a];} S/=4.0;
        for(int a=0;a<3;++a) xu[a]=R[a]-S; }
    double qD=0,qL=0; int nD=0,nL=0;
    for(int a=0;a<4;++a){ if(xd[a]<=0)continue; if(isD[a]){qD+=xd[a];++nD;} else {qL+=xd[a];++nL;} }
    for(int a=0;a<3;++a){ if(xu[a]<=0)continue; if(isD[4+a]){qD+=xu[a];++nD;} else {qL+=xu[a];++nL;} }
    if(nD<1||nL<1||qD<=0||qL<=0) return -1;
    return std::sqrt((qD/nD)/(qL/nL));
}

void mixedKillshotBootstrap(){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0); gStyle->SetGridStyle(3); gStyle->SetGridColor(kGray);
    gSystem->mkdir("papers/figures/mixed_killshot_bootstrap",kTRUE);
    const double Es[5]={50,75,100,125,150};
    TRandom3 rng(20260609);
    printf("\n===== GATE 6: MIXED KILL-SHOT BOOTSTRAP =====\n");
    printf("map (GATE 1): DSB1 = NE,SW (slots 1,3,5); LuAG = NW,SE (slots 0,2,4,6); SW-U excluded.\n");
    printf("M1 direct paired ratio convention: R = sigma(NE-D minus SW-D)/sigma(NW-D minus SE-D) = DSB1/LuAG.\n");

    // ---------- primary (srCFD) per energy: M1 + M2 + M3 + amplitudes ----------
    const int NB=1000;
    std::vector<double> Rdir(5,-1), Rsol(5,-1); std::vector<long> Nev(5,0);
    std::vector<std::vector<double>> boots(5);
    std::vector<double> sDv(5,-1), sLv(5,-1);
    std::vector<std::vector<Ev>> EV(5);
    double ampD=0, ampL=0; long ampN=0;
    for(int ie=0;ie<5;++ie){
        EV[ie]=gather(Es[ie],0); std::vector<Ev>& ev=EV[ie]; Nev[ie]=ev.size();
        if(ev.size()<5000){ printf("  %3.0f GeV: too few common events (%zu)\n",Es[ie],ev.size()); continue; }
        Direct d=direct(ev,1,3,0,2);
        Rdir[ie]=d.R; sDv[ie]=d.sD/std::sqrt(2.0)*1000.0; sLv[ie]=d.sL/std::sqrt(2.0)*1000.0;
        Rsol[ie]=solveRatio(ev);
        for(Ev& e:ev){ ampD+=0.5*(e.amp[1]+e.amp[3]); ampL+=0.5*(e.amp[0]+e.amp[2]); ++ampN; }
        // Poisson bootstrap of the direct ratio (frozen windows d.wD/d.wL)
        boots[ie].reserve(NB);
        std::vector<int> wt(ev.size());
        for(int b=0;b<NB;++b){ for(size_t i=0;i<wt.size();++i) wt[i]=rng.Poisson(1.0);
            double sD=winRMS(d.dD,d.wD,&wt), sL=winRMS(d.dL,d.wL,&wt);
            if(sD>0&&sL>0) boots[ie].push_back(sD/sL); }
        std::sort(boots[ie].begin(),boots[ie].end());
        auto q=[&](double f)->double{ return boots[ie][(size_t)(f*(boots[ie].size()-1))]; };
        printf("  %3.0f GeV: N=%6ld  sigmaCap(DSB1)=%5.0f ps  sigmaCap(LuAG)=%5.0f ps  "
               "R(direct)=%.3f [68%%: %.3f-%.3f] [95%%: %.3f-%.3f]  R(solve)=%.3f\n",
               Es[ie],Nev[ie],sDv[ie],sLv[ie],Rdir[ie],q(0.16),q(0.84),q(0.025),q(0.975),Rsol[ie]);
    }
    // 5-energy mean ratio bootstrap (replica-wise mean across energies)
    std::vector<double> meanBoot; long nbmin=1e9;
    for(int ie=0;ie<5;++ie) if(!boots[ie].empty()) nbmin=std::min(nbmin,(long)boots[ie].size());
    // re-shuffle each energy's replicas to decorrelate sorting, then average index-wise
    for(int ie=0;ie<5;++ie) std::shuffle(boots[ie].begin(),boots[ie].end(),std::mt19937(7+ie));
    for(long b=0;b<nbmin;++b){ double m=0;int n=0; for(int ie=0;ie<5;++ie) if(!boots[ie].empty()){m+=boots[ie][b];++n;}
        if(n) meanBoot.push_back(m/n); }
    std::sort(meanBoot.begin(),meanBoot.end());
    auto Q=[&](double f)->double{ return meanBoot.empty()?-1:meanBoot[(size_t)(f*(meanBoot.size()-1))]; };
    double Rmean=0; int nR=0; for(double r:Rdir) if(r>0){Rmean+=r;++nR;} Rmean/=std::max(1,nR);
    double RmeanSol=0; int nS=0; for(double r:Rsol) if(r>0){RmeanSol+=r;++nS;} RmeanSol/=std::max(1,nS);
    printf("\n  5-E mean direct ratio = %.3f   bootstrap median %.3f  [68%%: %.3f-%.3f]  [95%%: %.3f-%.3f]\n",
           Rmean,Q(0.5),Q(0.16),Q(0.84),Q(0.025),Q(0.975));
    printf("  5-E mean solve ratio (original quantity, GATE-1 labels) = %.3f\n",RmeanSol);
    printf("  amplitude control (correct map): <amp>_DSB1/<amp>_LuAG = %.2f\n",(ampD/ampN)/(ampL/ampN));

    // ---------- M4 run jackknife (srCFD, direct ratio) ----------
    printf("\n  M4 run jackknife (leave-one-run-out, runs with N>=700):\n");
    double jkMin=1e9,jkMax=-1e9,jkSum=0; int jkN=0; int worstRun=-1; double worstDev=0;
    std::vector<double> jkX,jkY;
    for(int ie=0;ie<5;++ie){ if(EV[ie].size()<5000) continue;
        std::map<int,long> nrun; for(Ev& e:EV[ie]) nrun[e.run]++;
        for(auto& kv:nrun){ if(kv.second<700) continue;
            std::vector<Ev> sub; sub.reserve(EV[ie].size()-kv.second);
            for(Ev& e:EV[ie]) if(e.run!=kv.first) sub.push_back(e);
            Direct d=direct(sub,1,3,0,2); if(d.R<0) continue;
            double m=0;int n=0; for(int je=0;je<5;++je){ double r=(je==ie)?d.R:Rdir[je]; if(r>0){m+=r;++n;} }
            m/=n; jkX.push_back(kv.first); jkY.push_back(m);
            jkSum+=m; ++jkN; jkMin=std::min(jkMin,m); jkMax=std::max(jkMax,m);
            if(std::fabs(m-Rmean)>worstDev){worstDev=std::fabs(m-Rmean);worstRun=kv.first;}
        } }
    printf("    %d leave-one-out variants: mean %.3f, spread (max-min) %.3f, worst run %d (|dev|=%.3f)\n",
           jkN,jkN?jkSum/jkN:-1,jkMax-jkMin,worstRun,worstDev);
    printf("    single-run dependence: %s\n",(jkMax-jkMin)<0.02?"NO (spread < 0.02)":"check worst run");

    // ---------- M5 method dependence (direct ratio per source) ----------
    printf("\n  M5 method dependence (5-E mean direct ratio):\n");
    double Rsrc[NSRC]={Rmean,0,0};
    for(int s=1;s<NSRC;++s){ double m=0;int n=0;
        for(int ie=0;ie<5;++ie){ std::vector<Ev> ev=gather(Es[ie],s); if(ev.size()<5000)continue;
            Direct d=direct(ev,1,3,0,2); if(d.R>0){m+=d.R;++n;} }
        Rsrc[s]=n?m/n:-1; }
    printf("    srCFD (primary): %.3f   cfd05 (original): %.3f   led: %.3f\n",Rsrc[0],Rsrc[1],Rsrc[2]);

    // ---------- M6 wrong-map diagnostics (srCFD, Down layer) ----------
    printf("\n  M6 wrong-map diagnostics (150 GeV):\n");
    { std::vector<Ev>& ev=EV[4];
      Direct dC=direct(ev,1,3,0,2);            // correct
      Direct dX=direct(ev,0,2,1,3);            // swapped: ratio must invert
      Direct dS=direct(ev,1,0,3,2);            // scrambled: cross-material "pairs"
      double aD=0,aL=0,aS1=0,aS2=0;
      for(Ev& e:ev){ aD+=0.5*(e.amp[1]+e.amp[3]); aL+=0.5*(e.amp[0]+e.amp[2]);
                     aS1+=0.5*(e.amp[1]+e.amp[0]); aS2+=0.5*(e.amp[3]+e.amp[2]); }
      printf("    correct map : R=%.3f   amplitude control=%.2f (expect >~2)\n",dC.R,aD/aL);
      printf("    swapped map : R=%.3f   (expect ~1/correct=%.3f)\n",dX.R,1.0/dC.R);
      printf("    scrambled   : R=%.3f   amplitude control=%.2f (expect ~1)\n",dS.R,aS1/aS2);
      printf("    NOTE: scrambled-pair widths include the material-difference variance and lose\n");
      printf("    the same-capillary correlations -> not directly comparable to the correct-map sigma.\n");
    }

    // ---------- figures ----------
    TCanvas* c=new TCanvas("ks","",1500,1050); c->Divide(2,2,0.006,0.010);
    // pad 1: bootstrap distribution of the 5-E mean ratio
    c->cd(1); gPad->SetLeftMargin(0.13);gPad->SetBottomMargin(0.12);gPad->SetTopMargin(0.08);gPad->SetGridy();
    { TH1F* hb=new TH1F("hb",";5-energy mean ratio  #sigma_{DSB1}/#sigma_{LuAG};bootstrap replicas",60,
          std::max(0.85,Q(0.001)-0.02), Q(0.999)+0.02); hb->SetDirectory(nullptr);
      for(double r:meanBoot)hb->Fill(r);
      hb->SetLineColor(kAzure+2);hb->SetFillColorAlpha(kAzure+2,0.30); hb->Draw("HIST");
      TLine* l1=new TLine(1.0,0,1.0,hb->GetMaximum()*1.02); l1->SetLineColor(kRed+1);l1->SetLineStyle(2);l1->SetLineWidth(2);l1->Draw();
      TLatex tx; tx.SetNDC(); tx.SetTextSize(0.037);
      tx.DrawLatex(0.16,0.86,Form("median %.3f",Q(0.5)));
      tx.DrawLatex(0.16,0.80,Form("68%%: [%.3f, %.3f]",Q(0.16),Q(0.84)));
      tx.DrawLatex(0.16,0.74,Form("95%%: [%.3f, %.3f]",Q(0.025),Q(0.975)));
      DrawPadTitle("event-bootstrap of the paired ratio (srCFD)"); }
    // pad 2: jackknife
    c->cd(2); gPad->SetLeftMargin(0.13);gPad->SetBottomMargin(0.12);gPad->SetTopMargin(0.08);gPad->SetGridy();
    { double xmn=1e9,xmx=-1e9; for(double x:jkX){xmn=std::min(xmn,x);xmx=std::max(xmx,x);}
      if(xmn>xmx){xmn=0;xmx=1;}
      TH1F* fr=gPad->DrawFrame(xmn-5,Rmean-0.05,xmx+5,Rmean+0.05);
      fr->SetTitle(";excluded run;5-E mean ratio (leave-one-out)");
      TGraph* g=new TGraph(jkX.size(),jkX.data(),jkY.data());
      g->SetMarkerStyle(20);g->SetMarkerColor(kGreen+3);g->SetMarkerSize(1.1);g->Draw("P SAME");
      TLine* l=new TLine(xmn-5,Rmean,xmx+5,Rmean); l->SetLineColor(kGray+2);l->SetLineStyle(2);l->Draw();
      TLine* l1=new TLine(xmn-5,1.0,xmx+5,1.0); l1->SetLineColor(kRed+1);l1->SetLineStyle(3);l1->Draw();
      DrawPadTitle("run jackknife"); }
    // pad 3: widths vs E
    c->cd(3); gPad->SetLeftMargin(0.13);gPad->SetBottomMargin(0.12);gPad->SetTopMargin(0.08);gPad->SetGridy();
    { TH1F* fr=gPad->DrawFrame(40,250,160,450);
      fr->SetTitle(";beam energy E (GeV);per-capillary #sigma_{t} (ps)");
      std::vector<double> ze(5,0), eD(5),eL(5);
      for(int i=0;i<5;++i){ double n=std::sqrt(2.0*std::max(1L,Nev[i])); eD[i]=sDv[i]/n*3; eL[i]=sLv[i]/n*3; }
      TGraphErrors* gD=new TGraphErrors(5,Es,sDv.data(),ze.data(),eD.data());
      gD->SetMarkerStyle(20);gD->SetMarkerColor(kRData);gD->SetLineColor(kRData);gD->SetMarkerSize(1.5);gD->SetLineWidth(2);gD->Draw("PL SAME");
      TGraphErrors* gL=new TGraphErrors(5,Es,sLv.data(),ze.data(),eL.data());
      gL->SetMarkerStyle(21);gL->SetMarkerColor(kGreen+3);gL->SetLineColor(kGreen+3);gL->SetMarkerSize(1.5);gL->SetLineWidth(2);gL->Draw("PL SAME");
      TLegend* lg=new TLegend(0.45,0.74,0.95,0.90); lg->SetBorderSize(0);lg->SetFillStyle(0);lg->SetTextSize(0.036);
      lg->AddEntry(gD,"DSB1 pair (NE#minusSW, Down)","lp"); lg->AddEntry(gL,"LuAG pair (NW#minusSE, Down)","lp"); lg->Draw();
      DrawPadTitle("same-shower per-capillary widths (srCFD)"); }
    // pad 4: method + map summary
    c->cd(4); gPad->SetTopMargin(0.08);
    { TLatex t4; t4.SetNDC(); t4.SetTextSize(0.040); double yy=0.88;
      auto L=[&](const char* s){ t4.DrawLatex(0.06,yy,s); yy-=0.073; };
      L("method dependence (5-E mean direct ratio):");
      L(Form("   srCFD %.3f    cfd05 %.3f    led %.3f",Rsrc[0],Rsrc[1],Rsrc[2]));
      L(Form("solve-method reproduction (GATE-1 labels): %.3f",RmeanSol));
      L("wrong-map diagnostics (150 GeV):");
      L("   swapped map inverts the ratio (1/R)  #checkmark");
      L("   scrambled map collapses the #times2#minus3 amplitude");
      L("   control to #approx1  #checkmark");
      DrawPadTitle("method + map controls"); }
    c->cd(0); DrawSuperTitle("GATE 6: MIXED in-event kill-shot hardening #minus paired bootstrap, jackknife, methods, map controls",0.019f);
    c->Print("papers/figures/mixed_killshot_bootstrap/killshot_bootstrap.png");
    printf("\n  wrote papers/figures/mixed_killshot_bootstrap/killshot_bootstrap.png\n");
}

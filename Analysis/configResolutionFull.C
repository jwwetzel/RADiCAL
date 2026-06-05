// configResolutionFull.C — follow-on paper gap #1+#2:
//   per-config sigma_t(E) (OOS best-bin, run-folded) AND sigma_E/E(E), all six
//   energies, all four builds (DSB1 / LuAG / MIXED / TENERGY), with a/sqrtE(+b/E)+c
//   fits and a comparison figure.  DSB1 reads the processRun ntuples; the other
//   three read the config-agnostic reduced ntuples.
//   root -l 'Analysis/configResolutionFull.C+'
#include "ChannelConfig.h"
#include "PlotUtils.h"
#include "TFile.h"
#include "DataPaths.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include <cstdio>
#include <vector>
#include <algorithm>
#include <cmath>

struct Ev { int run; float slg; float t; };

static double coreS(std::vector<float> v){
    if(v.size()<150) return -1;
    std::sort(v.begin(),v.end());
    double mu=v[v.size()/2], s=0.3;
    for(int it=0;it<5;++it){ double a=0,a2=0; long n=0;
        for(float x:v) if(std::fabs(x-mu)<2.5*s){a+=x;a2+=x*x;++n;}
        if(n<40) break; mu=a/n; double var=a2/n-mu*mu; s=var>0?std::sqrt(var):s; }
    return s;
}
static double coreSrel(std::vector<float>& v){ // sigma/mean of an energy distribution
    if(v.size()<150) return -1;
    std::vector<float> t=v; std::sort(t.begin(),t.end());
    double md=t[t.size()/2],a=0,a2=0; long n=0;
    for(float x:t) if(std::fabs(x-md)<0.4*md){a+=x;a2+=x*x;++n;}
    if(n<50) return -1; double mu=a/n, var=a2/n-mu*mu; return var>0? std::sqrt(var)/mu : -1;
}
// Gaussian-core sigma [ps] of a set of corner times (matches the official pipeline's FitGaussCore).
static double gcoreSig(std::vector<float>& v){
    if(v.size()<150) return -1;
    std::vector<float> s=v; std::sort(s.begin(),s.end()); double md=s[s.size()/2];
    TH1F h("_gc","",200,md-1.0,md+1.0); for(float x:v) h.Fill(x);
    double mu,muE,sg,sgErr; FitGaussCore(&h,2.0,mu,muE,sg,sgErr);
    return sg>0? sg*1000.0 : -1;   // ns -> ps
}
// 5-fold run-folded out-of-sample best-bin sigma_t [ps]; select min-sigma E_meas bin
// on training runs, MEASURE on held-out runs, pool across folds, fit once. Matches
// timingEnergyBins.C (Method A, #G5).
static double oosBest(std::vector<Ev>& ev, double& insamp, double& errOut, double& effOut){
    insamp=-1; errOut=0; effOut=0; if(ev.size()<2500) return -1;
    std::vector<float> sl; for(auto&e:ev) sl.push_back(e.slg); std::sort(sl.begin(),sl.end());
    double md=sl[sl.size()/2],a=0,a2=0; long n=0;
    for(float v:sl) if(std::fabs(v-md)<0.5*md){a+=v;a2+=v*v;++n;}
    double mE=a/n, sE=std::sqrt(a2/n-mE*mE), lo=mE-2*sE, bw=4*sE/9.0;
    std::vector<int> ur; for(auto&e:ev) ur.push_back(e.run); std::sort(ur.begin(),ur.end());
    ur.erase(std::unique(ur.begin(),ur.end()),ur.end());
    const int kF=5; auto fold=[&](int r){ int idx=(int)(std::lower_bound(ur.begin(),ur.end(),r)-ur.begin()); return idx%kF; };
    const int kTrainMin=(500*(kF-1))/kF;
    { double best=1e9; for(int b=0;b<9;++b){ double blo=lo+b*bw,bhi=blo+bw; std::vector<float> vt;
        for(auto&e:ev) if(e.slg>=blo&&e.slg<bhi) vt.push_back(e.t);
        if((long)vt.size()<500) continue; double s=gcoreSig(vt); if(s>0&&s<best) best=s; } insamp=best; }
    std::vector<float> pool;
    for(int f=0;f<kF;++f){ int bb=-1; double btr=1e9;
        for(int b=0;b<9;++b){ double blo=lo+b*bw,bhi=blo+bw; std::vector<float> vt;
            for(auto&e:ev) if(fold(e.run)!=f && e.slg>=blo && e.slg<bhi) vt.push_back(e.t);
            if((long)vt.size()<kTrainMin) continue; double s=gcoreSig(vt); if(s>0&&s<btr){btr=s;bb=b;} }
        if(bb<0) continue; double blo=lo+bb*bw,bhi=blo+bw;
        for(auto&e:ev) if(fold(e.run)==f && e.slg>=blo && e.slg<bhi) pool.push_back(e.t); }
    if(pool.size()<150) return -1;
    effOut = 100.0*pool.size()/ev.size();
    double s = gcoreSig(pool); errOut = (s>0)? s/std::sqrt(2.0*pool.size()) : 0;
    return s;
}

// ---- loaders: fill ev with fiducial, MCP-quality events; set sErel = sigma_E/E ----
static void loadReduced(const char* dir, double E, std::vector<Ev>& ev, double& sErel){
    int dw[4],up[4]; bool upM2[4]; int en[8];
    for(int i=0;i<4;++i) dw[i]=kCap[i].hg/1024;
    for(int i=4;i<8;++i){ up[i-4]=kCap[i].hg/1024; upM2[i-4]=kCap[i].use_mcp2; }
    for(int i=0;i<8;++i) en[i]=kCap[i].lg/1024;
    TFile* fp=TFile::Open(radReduced(dir,E)); if(!fp||fp->IsZombie()) return;
    TTree* t=(TTree*)fp->Get("rad"); if(!t){ fp->Close(); return; }
    Int_t run; Bool_t wc; Float_t x,y,m1t,m2t,m1p,sp[36],sc[36];
    t->SetBranchAddress("run",&run); t->SetBranchAddress("wc_ok",&wc);
    t->SetBranchAddress("x_trk",&x); t->SetBranchAddress("y_trk",&y);
    t->SetBranchAddress("mcp1_time",&m1t); t->SetBranchAddress("mcp2_time",&m2t);
    t->SetBranchAddress("mcp1_peak",&m1p); t->SetBranchAddress("s_peak",sp); t->SetBranchAddress("s_cfd05",sc);
    long N=t->GetEntries(); double xs=0,ys=0; long nw=0;
    for(long i=0;i<N&&nw<50000;++i){ t->GetEntry(i); if(wc&&x>-100&&x<100){xs+=x;ys+=y;++nw;} }
    double xc=nw?xs/nw:0, yc=nw?ys/nw:0; std::vector<float> elist;
    for(long i=0;i<N;++i){ t->GetEntry(i);
        if(!wc) continue; double dx=x-xc,dy=y-yc; if(dx*dx+dy*dy>=9.0) continue;
        if(m1p<200||m1p>750||m1t<-1e5) continue;
        double dsum=0,usum=0; int dn=0,un=0;
        for(int k=0;k<4;++k){ int s=dw[k]; if(sp[s]>20&&sc[s]>-1e5){ dsum+=sc[s]-m1t; ++dn; } }
        for(int k=0;k<4;++k){ int s=up[k]; double ref=upM2[k]?m2t:m1t; if(ref<-1e5)continue;
            if(sp[s]>20&&sc[s]>-1e5){ usum+=sc[s]-ref; ++un; } }
        double slg=0; for(int k=0;k<8;++k) slg+=sp[en[k]];
        if(dn>=1&&un>=1){ Ev e; e.run=run; e.slg=(float)slg; e.t=0.5f*(float)(dsum/dn-usum/un); ev.push_back(e); }
        elist.push_back((float)slg);
    }
    sErel = coreSrel(elist); fp->Close();
}
static void loadDSB1(double E, std::vector<Ev>& ev, double& sErel){
    TFile* fp=TFile::Open(radReduced("DSB1",E)); if(!fp||fp->IsZombie()) return;
    TTree* t=(TTree*)fp->Get("rad"); if(!t){ fp->Close(); return; }
    Int_t run; Bool_t wc; Float_t x,y,slg,mp,cfd[8];
    t->SetBranchAddress("run",&run); t->SetBranchAddress("wc_ok",&wc);
    t->SetBranchAddress("x_trk",&x); t->SetBranchAddress("y_trk",&y);
    t->SetBranchAddress("sum_lg",&slg); t->SetBranchAddress("mcp_peak",&mp); t->SetBranchAddress("hg_cfd05",cfd);
    long N=t->GetEntries(); double xs=0,ys=0; long nw=0;
    for(long i=0;i<N&&nw<50000;++i){ t->GetEntry(i); if(wc&&x>-100&&x<100){xs+=x;ys+=y;++nw;} }
    double xc=nw?xs/nw:0, yc=nw?ys/nw:0; std::vector<float> elist;
    for(long i=0;i<N;++i){ t->GetEntry(i);
        if(!wc) continue; double dx=x-xc,dy=y-yc; if(dx*dx+dy*dy>=9.0) continue;
        if(mp<200||mp>750) continue;
        double dsum=0,usum=0; int dn=0,un=0;
        for(int k=0;k<4;++k) if(cfd[k]>-1e5){ dsum+=cfd[k]; ++dn; }
        for(int k=4;k<8;++k) if(cfd[k]>-1e5){ usum+=cfd[k]; ++un; }
        if(dn>=1&&un>=1){ Ev e; e.run=run; e.slg=slg; e.t=0.5f*(float)(dsum/dn-usum/un); ev.push_back(e); }
        elist.push_back(slg);
    }
    sErel = coreSrel(elist); fp->Close();
}

void configResolutionFull(){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    const char* lab[4]={"DSB1","LuAG","MIXED","TENERGY"};
    const char* dir[4]={"","LUAG","MIXED","TENERGY"};
    int col[4]={kRRed, kAzure+1, 30, kOrange+7};
    double Es[6]={25,50,75,100,125,150};
    double sig[4][6], sigErr[4][6], sigE[4][6], eff[4][6];
    for(int c=0;c<4;++c) for(int e=0;e<6;++e){ sig[c][e]=-1; sigE[c][e]=-1; eff[c][e]=0; sigErr[c][e]=0; }

    for(int c=0;c<4;++c){
        printf("\n=== %s : sigma_t(E) OOS best-bin, sigma_E/E ===\n", lab[c]);
        printf("%5s %9s %12s %10s %9s %8s\n","E","N","sigt_OOS[ps]","(insamp)","eff[%]","sE/E[%]");
        for(int e=0;e<6;++e){
            std::vector<Ev> ev; double sErel=-1;
            if(c==0) loadDSB1(Es[e],ev,sErel); else loadReduced(dir[c],Es[e],ev,sErel);
            if(ev.size()<2000){ printf("%5.0f %9zu  (insufficient)\n",Es[e],ev.size()); continue; }
            double in,err,ef; double oos=oosBest(ev,in,err,ef);
            sig[c][e]=oos; sigErr[c][e]=err; sigE[c][e]=sErel*100.0; eff[c][e]=ef;
            printf("%5.0f %9zu %12.1f %10.1f %9.1f %8.1f\n",Es[e],ev.size(),oos,in,ef,sErel*100.0);
        }
    }

    // ---- sigma_t(E) figure with a/sqrtE (+) b fits ----
    TCanvas* c1=new TCanvas("c_st","",900,640); c1->SetGridx(); c1->SetGridy();
    TLegend* lg=new TLegend(0.55,0.62,0.88,0.88);
    TH1F* fr=c1->DrawFrame(15,15,165,80); fr->SetTitle("Per-build timing resolution;Beam energy (GeV);#sigma_{t} (OOS best-bin) [ps]");
    printf("\n=== sigma_t(E) = a/sqrtE (+) b  fits ===\n");
    for(int c=0;c<4;++c){
        std::vector<double> ex,ey,exe,eye;
        for(int e=0;e<6;++e) if(sig[c][e]>0){ ex.push_back(Es[e]); ey.push_back(sig[c][e]); exe.push_back(0); eye.push_back(sigErr[c][e]); }
        if(ex.size()<3) continue;
        TGraphErrors* g=new TGraphErrors(ex.size(),&ex[0],&ey[0],&exe[0],&eye[0]);
        g->SetMarkerStyle(20+c); g->SetMarkerColor(col[c]); g->SetLineColor(col[c]); g->SetMarkerSize(1.3);
        TF1* f=new TF1(Form("f%d",c),"sqrt([0]*[0]/x+[1]*[1])",20,160); f->SetParameters(200,25); f->SetLineColor(col[c]); f->SetLineWidth(2);
        g->Fit(f,"RQN"); g->Draw("P SAME"); f->Draw("SAME");
        printf("  %-8s a=%6.1f ps*sqrtGeV   b=%5.1f ps   (sigt@150=%.1f)\n",lab[c],f->GetParameter(0),fabs(f->GetParameter(1)),sig[c][5]);
        lg->AddEntry(g,Form("%s  (b=%.0f ps)",lab[c],fabs(f->GetParameter(1))),"p");
    }
    lg->Draw();
    TLatex tl; tl.SetNDC(); tl.SetTextSize(0.03); tl.DrawLatex(0.13,0.16,"OOS best-bin (DW-UP)/2, CFD-5%, r<3mm; fit #sigma_{t}=a/#sqrt{E}#oplus b");
    c1->Print("Analysis/capillary_figs/config_sigmat_vs_E.png");

    // ---- sigma_E/E(E) figure with 52%/sqrtE+...  fits ----
    TCanvas* c2=new TCanvas("c_se","",900,640); c2->SetGridx(); c2->SetGridy();
    TLegend* lg2=new TLegend(0.55,0.62,0.88,0.88);
    TH1F* fr2=c2->DrawFrame(15,8,165,24); fr2->SetTitle("Per-build shower-max energy resolution;Beam energy (GeV);#sigma_{E}/E [%]");
    printf("\n=== sigma_E/E(E) = a/sqrtE (+) b/E (+) c  fits ===\n");
    for(int c=0;c<4;++c){
        std::vector<double> ex,ey,exe,eye;
        for(int e=0;e<6;++e) if(sigE[c][e]>0){ ex.push_back(Es[e]); ey.push_back(sigE[c][e]); exe.push_back(0); eye.push_back(0.3); }
        if(ex.size()<3) continue;
        TGraphErrors* g=new TGraphErrors(ex.size(),&ex[0],&ey[0],&exe[0],&eye[0]);
        g->SetMarkerStyle(20+c); g->SetMarkerColor(col[c]); g->SetLineColor(col[c]); g->SetMarkerSize(1.3);
        TF1* f=new TF1(Form("fe%d",c),"sqrt([0]*[0]/x+[1]*[1]/(x*x)+[2]*[2])",20,160); f->SetParameters(50,30,9); f->SetLineColor(col[c]); f->SetLineWidth(2);
        g->Fit(f,"RQN"); g->Draw("P SAME"); f->Draw("SAME");
        printf("  %-8s a=%5.1f%%/sqrtGeV  b=%5.1f%%/GeV  c=%4.1f%%\n",lab[c],fabs(f->GetParameter(0)),fabs(f->GetParameter(1)),fabs(f->GetParameter(2)));
        lg2->AddEntry(g,lab[c],"p");
    }
    lg2->Draw();
    c2->Print("Analysis/capillary_figs/config_sigmaE_vs_E.png");
    printf("\nwrote Analysis/capillary_figs/config_sigmat_vs_E.png, config_sigmaE_vs_E.png\n");
}

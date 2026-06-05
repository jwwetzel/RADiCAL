// whyDiff.C — isolate why sigmaT gives 29.0 ps vs the published 27.4 ps (DSB1/150).
// Identical event selection (matches SelectionCuts.h: MCP1 in [200,750], r<3mm,
// (DW-UP)/2 over 4 corners, CFD-5%); vary ONLY the sigma-estimation method.
//   ROOT_INCLUDE_PATH=radcore:Analysis root -l -b -q radcore/whyDiff.C+
#include "RadView.h"
#include "PlotUtils.h"     // FitGaussCore
#include "DataPaths.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>

struct Ev { int run; float slg; float t; };

// FitGaussCore sigma [ps] with a configurable half-window [ns] and bin count
static double gfit(std::vector<float>& v, double halfWin, int nb) {
    if (v.size() < 150) return -1;
    std::vector<float> s = v; std::sort(s.begin(), s.end()); double md = s[s.size()/2];
    TH1F h("_g","",nb, md-halfWin, md+halfWin); for (float x : v) h.Fill(x);
    double mu,muE,sg,sgE; FitGaussCore(&h,2.0,mu,muE,sg,sgE); return sg>0 ? sg*1000.0 : -1;
}
// truncated-RMS core [ps] (the coreS estimator used by dsb1OOSandCorr.C)
static double trms(std::vector<float>& v) {
    if (v.size() < 150) return -1;
    std::sort(v.begin(), v.end()); double mu=v[v.size()/2], s=0.3;
    for (int it=0;it<6;++it){ double a=0,a2=0; long n=0; for(float x:v) if(std::fabs(x-mu)<2.5*s){a+=x;a2+=x*x;++n;} if(n<40)break; mu=a/n; double var=a2/n-mu*mu; s=var>0?std::sqrt(var):s; }
    return s*1000.0;
}

// best single energy-bin, IN-SAMPLE, with a given sigma estimator
template<class F>
static double inSampleBest(std::vector<Ev>& ev, F sig, int nb=9, int minN=500, double* bestE=nullptr) {
    std::vector<float> sl; for(auto&e:ev) sl.push_back(e.slg); std::sort(sl.begin(),sl.end());
    double md=sl[sl.size()/2],a=0,a2=0; long n=0; for(float v:sl) if(std::fabs(v-md)<0.5*md){a+=v;a2+=v*v;++n;}
    double mE=a/n,sE=std::sqrt(a2/n-mE*mE),lo=mE-2*sE,bw=4*sE/nb; double best=1e9;
    for(int b=0;b<nb;++b){ double blo=lo+b*bw,bhi=blo+bw; std::vector<float> vt;
        for(auto&e:ev) if(e.slg>=blo&&e.slg<bhi) vt.push_back(e.t);
        if((long)vt.size()<minN) continue; double s=sig(vt); if(s>0&&s<best){ best=s; if(bestE)*bestE=0.5*(blo+bhi); } }
    return best;
}
// OOS best-bin; foldByRun chooses run-level (true) vs event-index (false) folds;
// pool=true pools held-out events then one fit, else averages per-fold sigma.
template<class F>
static double oos(std::vector<Ev>& ev, F sig, int kF, bool foldByRun, bool pool) {
    std::vector<float> sl; for(auto&e:ev) sl.push_back(e.slg); std::sort(sl.begin(),sl.end());
    double md=sl[sl.size()/2],a=0,a2=0; long n=0; for(float v:sl) if(std::fabs(v-md)<0.5*md){a+=v;a2+=v*v;++n;}
    double mE=a/n,sE=std::sqrt(a2/n-mE*mE),lo=mE-2*sE,bw=4*sE/9.0;
    std::vector<int> ur; for(auto&e:ev) ur.push_back(e.run); std::sort(ur.begin(),ur.end()); ur.erase(std::unique(ur.begin(),ur.end()),ur.end());
    auto foldOf=[&](const Ev&e,long idx)->int{
        if(!foldByRun) return (int)(idx % kF);
        int ri=(int)(std::lower_bound(ur.begin(),ur.end(),e.run)-ur.begin()); return ri%kF; };
    std::vector<float> poolv; double acc=0; int nf=0;
    for(int f=0; f<kF; ++f){
        int bb=-1; double btr=1e9;
        for(int b=0;b<9;++b){ double blo=lo+b*bw,bhi=blo+bw; std::vector<float> vt; long idx=0;
            for(auto&e:ev){ if(foldOf(e,idx)!=f && e.slg>=blo && e.slg<bhi) vt.push_back(e.t); ++idx; }
            if((long)vt.size()<400) continue; double s=sig(vt); if(s>0&&s<btr){btr=s;bb=b;} }
        if(bb<0) continue; double blo=lo+bb*bw,bhi=blo+bw; std::vector<float> ho; long idx=0;
        for(auto&e:ev){ if(foldOf(e,idx)==f && e.slg>=blo && e.slg<bhi) ho.push_back(e.t); ++idx; }
        if(pool) poolv.insert(poolv.end(),ho.begin(),ho.end());
        else { double s=sig(ho); if(s>0){acc+=s;++nf;} }
    }
    if(pool) return poolv.size()<150?-1:sig(poolv);
    return nf?acc/nf:-1;
}

void whyDiff() {
    rad::BuildConfig cfg = rad::BuildConfig::Load("datasets/2023/configs/DSB1.json");
    TFile* fp = TFile::Open(radReduced("DSB1",150)); TTree* t=(TTree*)fp->Get("rad");
    rad::RadView v; v.attach(t,&cfg);
    long N=v.entries(); double xs=0,ys=0; long nw=0;
    for(long i=0;i<N&&nw<50000;++i){ v.get(i); if(v.wc_ok()&&v.x_trk()>-100&&v.x_trk()<100){xs+=v.x_trk();ys+=v.y_trk();++nw;} }
    double xc=xs/nw, yc=ys/nw;
    std::vector<Ev> ev;
    for(long i=0;i<N;++i){ v.get(i);
        if(!v.wc_ok()||v.mcp1_peak()<200||v.mcp1_peak()>750) continue;
        double dx=v.x_trk()-xc,dy=v.y_trk()-yc; if(dx*dx+dy*dy>=9.0) continue;
        double ds=0,us=0; int dn=0,un=0;
        for(int c=0;c<4;++c){ float tc=v.cfd05(c); if(tc>-1e5){ds+=tc;++dn;} }
        for(int c=4;c<8;++c){ float tc=v.cfd05(c); if(tc>-1e5){us+=tc;++un;} }
        if(dn<1||un<1) continue; Ev e; e.run=v.run(); e.slg=v.sum_lg(); e.t=0.5f*(float)(ds/dn-us/un); ev.push_back(e);
    }
    printf("\nDSB1/150 identical selection: %zu events  (published headline = 27.4 ps; in-sample 25.6)\n\n", ev.size());

    // refined event set: per-channel hg_peak>=20 quality cut (published line 453)
    std::vector<Ev> evQ;
    for(long i=0;i<N;++i){ v.get(i);
        if(!v.wc_ok()||v.mcp1_peak()<200||v.mcp1_peak()>750) continue;
        double dx=v.x_trk()-xc,dy=v.y_trk()-yc; if(dx*dx+dy*dy>=9.0) continue;
        double ds=0,us=0; int dn=0,un=0;
        for(int c=0;c<4;++c){ float tc=v.cfd05(c); if(tc>-1e5 && v.hg_peak(c)>=20.f){ds+=tc;++dn;} }
        for(int c=4;c<8;++c){ float tc=v.cfd05(c); if(tc>-1e5 && v.hg_peak(c)>=20.f){us+=tc;++un;} }
        if(dn<1||un<1) continue; Ev e; e.run=v.run(); e.slg=v.sum_lg(); e.t=0.5f*(float)(ds/dn-us/un); evQ.push_back(e);
    }
    // per-run alignment: subtract each run's median t (the ScanRunCenters CFD offset)
    auto alignPerRun=[](std::vector<Ev> in){ std::map<int,std::vector<float>> by;
        for(auto&e:in) by[e.run].push_back(e.t);
        std::map<int,float> med; for(auto&kv:by){ auto&w=kv.second; std::sort(w.begin(),w.end()); med[kv.first]=w[w.size()/2]; }
        for(auto&e:in) e.t-=med[e.run]; return in; };

    auto gF = [](std::vector<float>& x){ return gfit(x, 0.4, 200); };   // FitGaussCore +-0.4ns
    std::vector<Ev> evA  = alignPerRun(ev);
    std::vector<Ev> evQA = alignPerRun(evQ);

    printf("IN-SAMPLE best energy-bin (FitGaussCore), building up the published refinements:\n");
    printf("   sigmaT (loose, cfd05 valid only)                 = %.1f ps\n", inSampleBest(ev , gF));
    printf("   + per-channel hg_peak>=20 quality cut            = %.1f ps   (%zu evt)\n", inSampleBest(evQ, gF), evQ.size());
    printf("   + per-run CFD-offset alignment (ScanRunCenters)  = %.1f ps\n", inSampleBest(evA , gF));
    printf("   + BOTH quality cut AND per-run alignment         = %.1f ps   <== published in-sample 25.6\n", inSampleBest(evQA, gF));
    printf("\nOOS best-bin (5-fold by RUN, pool) of the fully-refined set = %.1f ps   <== published headline 27.4\n",
           oos(evQA, gF, 5, true, true));

    // per-RUN beam centroid (ScanRunCenters) instead of one global centroid
    std::map<int,double> rx,ry; std::map<int,long> rn;
    for(long i=0;i<N;++i){ v.get(i); if(v.wc_ok()&&v.x_trk()>-100&&v.x_trk()<100){ rx[v.run()]+=v.x_trk(); ry[v.run()]+=v.y_trk(); rn[v.run()]++; } }
    std::map<int,double> cx,cy; for(auto&kv:rn){ cx[kv.first]=rx[kv.first]/kv.second; cy[kv.first]=ry[kv.first]/kv.second; }
    std::vector<Ev> evR;
    for(long i=0;i<N;++i){ v.get(i);
        if(!v.wc_ok()||v.mcp1_peak()<200||v.mcp1_peak()>750) continue;
        double dx=v.x_trk()-cx[v.run()],dy=v.y_trk()-cy[v.run()]; if(dx*dx+dy*dy>=9.0) continue;
        double ds=0,us=0; int dn=0,un=0;
        for(int c=0;c<4;++c){ float tc=v.cfd05(c); if(tc>-1e5 && v.hg_peak(c)>=20.f){ds+=tc;++dn;} }
        for(int c=4;c<8;++c){ float tc=v.cfd05(c); if(tc>-1e5 && v.hg_peak(c)>=20.f){us+=tc;++un;} }
        if(dn<1||un<1) continue; Ev e; e.run=v.run(); e.slg=v.sum_lg(); e.t=0.5f*(float)(ds/dn-us/un); evR.push_back(e);
    }
    // per-bin sigma_t profile across the 9 published-style energy bins (muE +- 2 sigE)
    auto gW = [](std::vector<float>& x){ return gfit(x, 1.0, 200); };   // == published timing fit window
    std::vector<float> sl; for(auto&e:evQ) sl.push_back(e.slg); std::sort(sl.begin(),sl.end());
    double md=sl[sl.size()/2],a=0,a2=0; long nn=0; for(float x:sl) if(std::fabs(x-md)<0.5*md){a+=x;a2+=x*x;++nn;}
    double mE=a/nn,sE=std::sqrt(a2/nn-mE*mE),lo=mE-2*sE,bw=4*sE/9.0;
    printf("\nPER-BIN sigma_t profile (9 bins muE+-2sigE, FitGaussCore +-1ns) [quality set]:\n");
    double best=1e9; std::vector<float> upperHalf;
    for(int b=0;b<9;++b){ double blo=lo+b*bw,bhi=blo+bw; std::vector<float> vt;
        for(auto&e:evQ) if(e.slg>=blo&&e.slg<bhi) vt.push_back(e.t);
        double s = vt.size()>=400 ? gW(vt) : -1; if(s>0&&s<best)best=s;
        if(b>=5) upperHalf.insert(upperHalf.end(),vt.begin(),vt.end());   // published "best region" = bins 5..8
        printf("   bin %d  sum_lg~%5.0f mV  N=%6zu  sigma_t = %5.1f ps%s\n",
               b, 0.5*(blo+bhi), vt.size(), s, b>=5?"  <-- published pools these":""); }
    printf("\n   sigmaT single-min-bin = %.1f ps   (best = top energy bin)\n", best);

    // TIMING-FIT HISTOGRAM RESOLUTION scan on the brightest events (bins 7-8)
    std::vector<float> bright; for(auto&e:evQ) if(e.slg >= lo+7*bw) bright.push_back(e.t);
    printf("\nFIT-RESOLUTION scan on the brightest slice (bins 7-8, N=%zu, FitGaussCore):\n", bright.size());
    printf("   +-1.0ns/200b (10ps bins, sigmaT) = %.1f ps\n", gfit(bright,1.0,200));
    printf("   +-0.4ns/200b ( 4ps bins)         = %.1f ps\n", gfit(bright,0.4,200));
    printf("   +-0.2ns/200b ( 2ps bins)         = %.1f ps\n", gfit(bright,0.2,200));
    printf("   +-0.4ns/400b ( 2ps bins)         = %.1f ps\n", gfit(bright,0.4,400));
    printf("   trunc-RMS (binning-free)         = %.1f ps\n", trms(bright));

    // EXACT published bin definition: FitGaussCore on the sum_lg histogram -> muE,sigE
    { std::vector<float> sv; for(auto&e:evQ) sv.push_back(e.slg);
      double smin=*std::min_element(sv.begin(),sv.end()), smax=*std::max_element(sv.begin(),sv.end());
      TH1F hS("hS","",150,smin,smax); for(float x:sv) hS.Fill(x);
      double muE,muEe,sigE,sigEe; FitGaussCore(&hS,2.0,muE,muEe,sigE,sigEe);
      double blo2=muE-2*sigE, bw2=4*sigE/9.0, best2=1e9; int bb=-1;
      for(int b=0;b<9;++b){ double a1=blo2+b*bw2,a2v=a1+bw2; std::vector<float> vt;
          for(auto&e:evQ) if(e.slg>=a1&&e.slg<a2v) vt.push_back(e.t);
          if(vt.size()<400) continue; double s=gfit(vt,1.0,200); if(s>0&&s<best2){best2=s;bb=b;} }
      printf("\nPUBLISHED-EXACT bin def (FitGaussCore sum_lg: muE=%.0f sigE=%.0f mV; 9 bins muE+-2sigE):\n", muE, sigE);
      printf("   best-bin sigma_t = %.1f ps  (bin %d)\n", best2, bb);
    }

    // THE FIX: per-run LG-light-weighted beam centroid (ScanRunCenters), not a global position mean
    std::map<int,double> wx,wy,ww;
    for(long i=0;i<N;++i){ v.get(i); if(v.wc_ok()&&v.x_trk()>-100&&v.x_trk()<100){ double w=v.sum_lg(); wx[v.run()]+=w*v.x_trk(); wy[v.run()]+=w*v.y_trk(); ww[v.run()]+=w; } }
    std::map<int,double> lcx,lcy; for(auto&kv:ww){ lcx[kv.first]=wx[kv.first]/kv.second; lcy[kv.first]=wy[kv.first]/kv.second; }
    std::vector<Ev> evL;
    for(long i=0;i<N;++i){ v.get(i);
        if(!v.wc_ok()||v.mcp1_peak()<200||v.mcp1_peak()>750) continue;
        double dx=v.x_trk()-lcx[v.run()],dy=v.y_trk()-lcy[v.run()]; if(dx*dx+dy*dy>=9.0) continue;
        double ds=0,us=0; int dn=0,un=0;
        for(int c=0;c<4;++c){ float tc=v.cfd05(c); if(tc>-1e5 && v.hg_peak(c)>=20.f){ds+=tc;++dn;} }
        for(int c=4;c<8;++c){ float tc=v.cfd05(c); if(tc>-1e5 && v.hg_peak(c)>=20.f){us+=tc;++un;} }
        if(dn<1||un<1) continue; Ev e; e.run=v.run(); e.slg=v.sum_lg(); e.t=0.5f*(float)(ds/dn-us/un); evL.push_back(e);
    }
    printf("\n=== THE FIX: per-run LG-weighted centroid (ScanRunCenters) ===\n");
    printf("   fiducial events: %zu   (was %zu with crude global mean; published collects 155596)\n", evL.size(), evQ.size());
    printf("   in-sample best-bin sigma_t = %.1f ps   <== published headline 27.4 ps\n", inSampleBest(evL, gW, 9, 500));
    fp->Close();
}

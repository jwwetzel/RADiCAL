// ============================================================================
// sigmaProbe.C — exhaustive sigma_t(E) diagnostic to isolate the CAUSE of any
// non-monotonicity (premise: a correct sigma_t-vs-E curve must be monotonically
// decreasing to a floor; any wiggle is a FIT or SELECTION artifact, not physics).
//
// For one (build, cornerSet, src, fidMode), at each energy it gathers the SAME
// pool of valid (DW-UP)/2 timing values, then crosses two axes on identical data:
//   SELECTION axis (which events): K500 / K1000 / K2000 brightest-by-sum_lg,
//      PCT2 = brightest fixed 2% (matched percentile across energies), and the
//      mean/median selected light + N (so a stats/percentile artifact is visible).
//   FIT axis (how sigma is read): gaus = FitGaussCore (current), trunc = 2.5sigma
//      truncated RMS (current fallback), iqr = IQR/1.349, rms3 = RMS within 3sigma.
//   Plus per-point skew & excess-kurtosis of the selected values (non-Gaussianity).
// If the wiggle moves with the FIT estimator -> fit problem. If it moves with the
// SELECTION scheme -> selection problem. If it survives every cell -> deeper.
//
//   cornerSet: 0=all4  1=DSB1 diag(NE,SW)  2=LuAG diag(NW,SE)
//   src (RadView enum): 1=cfd05  6=led  7=lgcfd
//   fidMode: 0=default(TimingFiducialR + per-file beamCenter)  1=fixed r=2.5
//   source setup.sh
//   root -l -b -q -e '.L analyze/studies/sigmaProbe.C+' -e 'sigmaProbe("MIXED",0,7,0)'
// Pure stdout (machine-parseable PROBE/ROW lines); no figure, no file writes.
// ============================================================================
#include "RadView.h"
#include "DataPaths.h"
#include "PlotUtils.h"       // FitGaussCore
#include "SelectionCuts.h"   // TimingFiducialR, kMCP1_*, kHG_minPeak
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TSystem.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>
using namespace rad;

static void cornersFor(int cs, std::vector<int>& dn, std::vector<int>& up){
    if(cs==1){ dn={1,3}; up={5,7}; }       // DSB1 diagonal NE,SW
    else if(cs==2){ dn={0,2}; up={4,6}; }  // LuAG diagonal NW,SE
    else { dn={0,1,2,3}; up={4,5,6,7}; }   // all 4 corners
}
static const char* csName(int cs){ return cs==1?"DSB1diag":cs==2?"LuAGdiag":"all4"; }
static const char* srcName(int s){ return s==7?"lgcfd":s==6?"led":s==1?"cfd05":"src?"; }

// ---- sigma estimators on a value vector (ns in, ps out) ----
static double estTrunc(const std::vector<float>& v){            // tebSigma's robust core
    if(v.size()<20) return -1; double mu=0; for(float x:v)mu+=x; mu/=v.size();
    double ms=0; for(float x:v)ms+=(x-mu)*(x-mu); ms=std::sqrt(ms/v.size()); if(ms<0.008)ms=0.1;
    double rc=mu,rw=ms; for(int it=0;it<5;++it){double s=0,ss=0;long n=0;for(float x:v)if(std::fabs(x-rc)<2.5*rw){s+=x;ss+=x*x;++n;}if(n<20)break;rc=s/n;double w2=ss/n-rc*rc;if(w2>0)rw=std::sqrt(w2);}
    return rw/0.9546*1000.0;
}
static double estGaus(const std::vector<float>& v){             // FitGaussCore (current headline)
    if(v.size()<50) return -1; double mu1=0; for(float x:v)mu1+=x; mu1/=v.size();
    double ms1=0; for(float x:v)ms1+=(x-mu1)*(x-mu1); ms1=std::sqrt(ms1/v.size()); if(ms1<0.008)ms1=0.1;
    double mu2=0;int n2=0;for(float x:v)if(std::fabs(x-mu1)<5*ms1){mu2+=x;++n2;} double ms2=ms1;
    if(n2>0){mu2/=n2;ms2=0;for(float x:v)if(std::fabs(x-mu1)<5*ms1)ms2+=(x-mu2)*(x-mu2);ms2=std::sqrt(ms2/n2);if(ms2<0.008)ms2=0.1;}else mu2=mu1;
    TH1F h("hpb","",120,mu2-4*ms2,mu2+4*ms2); h.SetDirectory(nullptr); for(float x:v)h.Fill(x);
    double mu,muE,s,sE; FitGaussCore(&h,2.0,mu,muE,s,sE); return s>0? s*1000.0 : -1;
}
static double estIQR(std::vector<float> v){
    if(v.size()<20) return -1; std::sort(v.begin(),v.end());
    double q1=v[(size_t)(0.25*v.size())], q3=v[(size_t)(0.75*v.size())]; return (q3-q1)/1.349*1000.0;
}
static double estRMS3(const std::vector<float>& v){
    if(v.size()<20) return -1; double mu=0; for(float x:v)mu+=x; mu/=v.size();
    double ms=0; for(float x:v)ms+=(x-mu)*(x-mu); ms=std::sqrt(ms/v.size()); if(ms<1e-4)ms=0.1;
    double s=0,ss=0;long n=0;for(float x:v)if(std::fabs(x-mu)<3*ms){s+=x;ss+=x*x;++n;} if(n<20)return -1;
    double m=s/n,w2=ss/n-m*m; return w2>0? std::sqrt(w2)*1000.0 : -1;
}
static void moments(const std::vector<float>& v,double& sk,double& ku){
    sk=ku=0; if(v.size()<20)return; double mu=0;for(float x:v)mu+=x;mu/=v.size();
    double m2=0,m3=0,m4=0;for(float x:v){double d=x-mu;m2+=d*d;m3+=d*d*d;m4+=d*d*d*d;}
    m2/=v.size();m3/=v.size();m4/=v.size(); if(m2<=0)return; sk=m3/std::pow(m2,1.5); ku=m4/(m2*m2)-3.0;
}

void sigmaProbe(const char* build, int cornerSet=0, int src=7, int fidMode=0, int cleanMode=0, double W=2.0, int rankMode=0){
    std::vector<int> dn,up; cornersFor(cornerSet,dn,up);
    BuildConfig cfg=BuildConfig::Load(radConfig(build).Data());
    const double Es[6]={25,50,75,100,125,150};
    printf("PROBE build=%s corners=%s src=%s fid=%s clean=%s%s rank=%s\n",build,csName(cornerSet),srcName(src),fidMode?"fixedR2.5":"default",
           cleanMode?"inEventConsistency":"none", cleanMode?Form("(W=%.1fns)",W):"", rankMode?"HGtimingSum":"sum_lg");
    printf("# E   Nfid  selK1000light medlight | est{gaus trunc iqr rms3} per sel{K500 K1000 K2000 PCT2}  skew kurt(@K1000)\n");
    for(double E:Es){ TString p=radReduced(build,E); if(gSystem->AccessPathName(p.Data())){ printf("ROW E=%.0f N=0 (missing)\n",E); continue; }
        TFile* f=TFile::Open(p.Data()); if(!f||f->IsZombie()){ if(f)f->Close(); printf("ROW E=%.0f N=0 (zombie)\n",E); continue; }
        TTree* t=(TTree*)f->Get("rad"); if(!t){ f->Close(); printf("ROW E=%.0f N=0 (notree)\n",E); continue; }
        RadView v; v.attach(t,&cfg); double xc,yc; v.beamCenter(xc,yc);
        double R=fidMode?2.5:TimingFiducialR(E); double r2=R*R;
        std::vector<std::pair<float,float>> sd; Long64_t N=v.entries();
        for(Long64_t i=0;i<N;++i){ v.get(i);
            if(!v.wc_ok()||v.mcp1_peak()<kMCP1_minPeak||v.mcp1_peak()>kMCP1_maxPeak) continue;
            double dx=v.x_trk()-xc,dy=v.y_trk()-yc; if(dx*dx+dy*dy>=r2) continue;
            // optional in-event consistency cleaning: drop any channel whose crossing time
            // disagrees with the event's median channel time by > W ns (a broken/wrong-feature
            // crossing or a near-noise channel). One rule kills both pathologies.
            float medt=0;
            if(cleanMode){ std::vector<float> at;
                for(int c:dn) if(v.is_timing(c)&&v.hg_peak(c)>=kHG_minPeak){float tc=v.timeOf(c,src);if(tc>-1e5f)at.push_back(tc);}
                for(int c:up) if(v.is_timing(c)&&v.hg_peak(c)>=kHG_minPeak){float tc=v.timeOf(c,src);if(tc>-1e5f)at.push_back(tc);}
                if(at.empty()) continue; std::nth_element(at.begin(),at.begin()+at.size()/2,at.end()); medt=at[at.size()/2]; }
            double ds=0,us=0,hgs=0;int nd=0,nu=0;
            for(int c:dn) if(v.is_timing(c)&&v.hg_peak(c)>=kHG_minPeak){float tc=v.timeOf(c,src);if(tc>-1e5f&&(!cleanMode||std::fabs(tc-medt)<W)){ds+=tc;hgs+=v.hg_peak(c);++nd;}}
            for(int c:up) if(v.is_timing(c)&&v.hg_peak(c)>=kHG_minPeak){float tc=v.timeOf(c,src);if(tc>-1e5f&&(!cleanMode||std::fabs(tc-medt)<W)){us+=tc;hgs+=v.hg_peak(c);++nu;}}
            if(nd<1||nu<1) continue;
            // rank by HG timing-pulse sum (clean-pulse selection) or by LG shower sum
            sd.push_back({ rankMode? (float)hgs : (float)v.sum_lg(), 0.5f*(float)(ds/nd-us/nu)});
        }
        long Nfid=sd.size();
        // sort by light descending once
        std::sort(sd.begin(),sd.end(),[](const std::pair<float,float>&a,const std::pair<float,float>&b){return a.first>b.first;});
        auto topK=[&](int K)->std::vector<float>{ std::vector<float> o; int kk=std::min((long)K,Nfid); for(int i=0;i<kk;++i)o.push_back(sd[i].second); return o; };
        int Kpct=std::max(200,(int)(0.02*Nfid));
        std::vector<float> vK500=topK(500), vK1000=topK(1000), vK2000=topK(2000), vPCT=topK(Kpct);
        double l1000=0,med=0; { int kk=std::min(1000L,Nfid); for(int i=0;i<kk;++i)l1000+=sd[i].first; if(kk)l1000/=kk; if(Nfid>0)med=sd[Nfid/2].first; }
        double sk=0,ku=0; moments(vK1000,sk,ku);
        auto S=[&](const std::vector<float>& vv)->std::string{ char b[96];
            snprintf(b,sizeof(b),"%.1f/%.1f/%.1f/%.1f",estGaus(vv),estTrunc(vv),estIQR(vv),estRMS3(vv)); return std::string(b); };
        std::string s5=S(vK500), s1=S(vK1000), s2=S(vK2000), sp=S(vPCT);
        printf("ROW E=%.0f N=%ld L1000=%.0f Lmed=%.0f | K500[%s] K1000[%s] K2000[%s] PCT2(K=%d)[%s] | skew=%.2f kurt=%.2f\n",
               E,Nfid,l1000,med, s5.c_str(), s1.c_str(), s2.c_str(), Kpct, sp.c_str(), sk,ku);
        f->Close();
    }
    printf("ENDPROBE\n");
}

// walkCorrTest.C — measure the CFD-20% amplitude-walk curve on the Down
//   capillaries and test whether a data-driven walk correction collapses the
//   shoulder.  Honest train/test split (odd events train the correction, even
//   events are evaluated).  Compares raw CFD-20%, walk-corrected CFD-20%,
//   CFD-5%, and LED.  50 GeV.
//
// Run: ROOT_INCLUDE_PATH=Analysis root -l -b -q 'Analysis/walkCorrTest.C+'
#include <cstdio>
#include <vector>
#include <algorithm>
#include <cmath>
#include "TFile.h"
#include "TTree.h"

namespace {
inline bool ok(float t){ return t>-60.f && t<60.f; }
double median(std::vector<float> v){ if(v.empty())return 0; std::sort(v.begin(),v.end());
    size_t n=v.size(); return (n&1)?v[n/2]:0.5*(v[n/2-1]+v[n/2]); }
double trms(const std::vector<float>&v,double c,double w){ double s=0,s2=0;int k=0;
    for(float x:v) if(std::fabs(x-c)<=w){s+=x;s2+=x*x;++k;} if(k<2)return 0;
    double m=s/k; double var=s2/k-m*m; return var>0?std::sqrt(var):0; }
}

void walkCorrTest(){
    TFile f("output/50GeV/ntuple.root");
    TTree* t=(TTree*)f.Get("rad");
    Float_t hg_peak[8],hg_cfd[8],hg_cfd05[8],hg_led[8]; Bool_t in_fid;
    t->SetBranchAddress("hg_peak",hg_peak);
    t->SetBranchAddress("hg_cfd",hg_cfd);
    t->SetBranchAddress("hg_cfd05",hg_cfd05);
    t->SetBranchAddress("hg_led",hg_led);
    t->SetBranchAddress("in_fiducial",&in_fid);

    const char* nm[8]={"NW-D","NE-D","SE-D","SW-D","NW-U","NE-U","SE-U","SW-U"};
    const int chans[4]={0,1,2,6};
    const Long64_t N=t->GetEntries();

    // walk variable: u = 1/sqrt(A)  (leading-edge walk ~ noise/slope ~ 1/(dV/dt) ~ 1/A
    //  near threshold; 1/sqrt(A) is a robust monotone proxy that also linearises
    //  the rounded-edge case).  We fit t20 = a + b*u by OLS on the TRAIN set,
    //  then also do a data-driven profile (bin-and-subtract) as a model-free check.
    for(int ci=0;ci<4;++ci){
        const int ch=chans[ci];
        double gC05raw=0, gC05cor=0;
        std::vector<float> A_tr,t_tr, A5_tr,t5_tr;             // train (odd)
        std::vector<float> A_te,t20_te,t05_te,tled_te, A5_te; // test (even)
        for(Long64_t e=0;e<N;++e){ t->GetEntry(e);
            if(!in_fid||!ok(hg_cfd[ch])||hg_peak[ch]<=20) continue;
            float A=hg_peak[ch];
            if(e&1){ A_tr.push_back(A); t_tr.push_back(hg_cfd[ch]);
                     if(ok(hg_cfd05[ch])){ A5_tr.push_back(A); t5_tr.push_back(hg_cfd05[ch]); } }
            else   { A_te.push_back(A); t20_te.push_back(hg_cfd[ch]);
                     t05_te.push_back(ok(hg_cfd05[ch])?hg_cfd05[ch]:999);
                     tled_te.push_back(ok(hg_led[ch])?hg_led[ch]:999);
                     A5_te.push_back(A); }
        }
        // OLS walk fit for CFD-5%: t5 = a5 + b5*u
        { double su=0,st=0,suu=0,sut=0; int n5=A5_tr.size();
          for(int i=0;i<n5;++i){ double u=1.0/std::sqrt(A5_tr[i]); su+=u; st+=t5_tr[i]; suu+=u*u; sut+=u*t5_tr[i]; }
          double b5=(n5*sut-su*st)/(n5*suu-su*su); double a5=(st-b5*su)/n5;
          // evaluate CFD-5% raw vs corrected on test
          std::vector<float> c05r, c05c;
          for(size_t i=0;i<A5_te.size();++i){ if(t05_te[i]<900){ c05r.push_back(t05_te[i]);
              double u=1.0/std::sqrt(A5_te[i]); c05c.push_back(t05_te[i]-(a5+b5*u)); } }
          gC05raw=trms(c05r,median(c05r),1.5)*1000; gC05cor=trms(c05c,median(c05c),1.5)*1000; }
        // OLS fit t = a + b*u, u=1/sqrt(A), on train
        double su=0,st=0,suu=0,sut=0; int n=A_tr.size();
        for(int i=0;i<n;++i){ double u=1.0/std::sqrt(A_tr[i]); su+=u; st+=t_tr[i]; suu+=u*u; sut+=u*t_tr[i]; }
        double b=(n*sut-su*st)/(n*suu-su*su); double a=(st-b*su)/n;

        // data-driven profile on train: 12 quantile-ish amplitude bins, mean t per bin
        const int NB=12; std::vector<float> As=A_tr; std::sort(As.begin(),As.end());
        double q[NB+1]; for(int i=0;i<=NB;++i) q[i]=As[std::min((size_t)(As.size()*i/NB),As.size()-1)];
        double pm[NB]={0}; int pc[NB]={0}; double pa[NB]={0};
        for(int i=0;i<n;++i){ for(int b2=0;b2<NB;++b2) if(A_tr[i]>=q[b2]&&A_tr[i]<=q[b2+1]){ pm[b2]+=t_tr[i]; pa[b2]+=A_tr[i]; ++pc[b2]; break; } }
        for(int b2=0;b2<NB;++b2){ if(pc[b2]){pm[b2]/=pc[b2]; pa[b2]/=pc[b2];} }
        double pglob=median(t_tr);
        auto profCorr=[&](double A)->double{ // interpolate mean-t(A), return offset to subtract
            if(A<=pa[0])return pm[0]-pglob; if(A>=pa[NB-1])return pm[NB-1]-pglob;
            for(int b2=0;b2<NB-1;++b2) if(A>=pa[b2]&&A<pa[b2+1]){
                double f=(A-pa[b2])/(pa[b2+1]-pa[b2]); return (pm[b2]+f*(pm[b2+1]-pm[b2]))-pglob; }
            return 0; };

        // evaluate on test
        std::vector<float> raw=t20_te, corrU, corrP, c05, led;
        for(size_t i=0;i<A_te.size();++i){
            double u=1.0/std::sqrt(A_te[i]);
            corrU.push_back(t20_te[i]-(a+b*u)+ (a+b*1.0/std::sqrt(500.0)) ); // re-center at A=500
            corrP.push_back(t20_te[i]-profCorr(A_te[i]));
            if(t05_te[i]<900) c05.push_back(t05_te[i]);
            if(tled_te[i]<900) led.push_back(tled_te[i]);
        }
        printf("\n==== %s (ch%d) 50 GeV  [test N=%zu] ====\n",nm[ch],ch,A_te.size());
        printf("  walk fit: t20 = %.3f + %.3f/sqrt(A)   (b in ns*sqrt(mV))\n",a,b);
        printf("  CFD-20%% raw            coreSig = %6.1f ps\n", trms(raw,  median(raw),1.5)*1000);
        printf("  CFD-20%% + 1/sqrtA fit  coreSig = %6.1f ps\n", trms(corrU,median(corrU),1.5)*1000);
        printf("  CFD-20%% + profile corr coreSig = %6.1f ps\n", trms(corrP,median(corrP),1.5)*1000);
        printf("  CFD-5%% (no corr)       coreSig = %6.1f ps\n", trms(c05, median(c05),1.5)*1000);
        printf("  CFD-5%% + 1/sqrtA fit   coreSig = %6.1f ps\n", gC05cor);
        printf("  LED-20mV (no corr)     coreSig = %6.1f ps\n", trms(led, median(led),1.5)*1000);
    }
    printf("\n[walkCorrTest] done.\n");
}

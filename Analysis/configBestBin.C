// configBestBin.C — best-contained-bin timing (the DSB1 headline method) on a
// reduced config dataset. Bins fiducial events by SumLG into 9 equal-width bins
// over <E>+/-2sigma, computes the corner (DW-UP)/2 CFD-5% core sigma per bin
// (N>=500), reports the best (min) bin + its efficiency.
// NOTE: in-sample best-bin — needs run-folded OOS validation before quoting as a
// headline (see the DSB1 overfitting lesson). First look only.
//   root -l 'Analysis/configBestBin.C+("reduced/LUAG","LuAG")'
#include "TFile.h"
#include "TTree.h"
#include "ChannelConfig.h"
#include <cstdio>
#include <vector>
#include <algorithm>
#include <cmath>

static double coreSig(std::vector<float> v){
    if(v.size()<200) return -1;
    std::sort(v.begin(),v.end()); double mu=v[v.size()/2], s=0.2;
    for(int it=0;it<5;++it){ double a=0,a2=0; long n=0;
        for(float x:v){ if(std::fabs(x-mu)<2.5*s){a+=x;a2+=x*x;++n;} }
        if(n<50)break; mu=a/n; double var=a2/n-mu*mu; s=var>0?std::sqrt(var):s; }
    return s;
}

void configBestBin(const char* dir, const char* label="")
{
    int dw[4],up[4]; bool upM2[4];
    for(int i=0;i<4;++i) dw[i]=kCap[i].hg/1024;
    for(int i=4;i<8;++i){ up[i-4]=kCap[i].hg/1024; upM2[i-4]=kCap[i].use_mcp2; }
    int en[8]; for(int i=0;i<8;++i) en[i]=kCap[i].lg/1024;

    double Es[6]={25,50,75,100,125,150};
    printf("\n=== %s  (best-contained-bin, in-sample) ===\n", label[0]?label:dir);
    printf("%5s %12s %10s %8s %8s\n","E","sigma_t[ps]","bin<SumLG>","eff%","Nbin");
    for(int e=0;e<6;++e){
        TFile* fp=TFile::Open(Form("%s/%.0fGeV.root",dir,Es[e])); if(!fp||fp->IsZombie())continue;
        TTree* t=(TTree*)fp->Get("rad"); if(!t){fp->Close();continue;}
        Bool_t wc; Float_t x,y,m1t,m2t,m1p,sp[36],sc[36];
        t->SetBranchAddress("wc_ok",&wc);t->SetBranchAddress("x_trk",&x);t->SetBranchAddress("y_trk",&y);
        t->SetBranchAddress("mcp1_time",&m1t);t->SetBranchAddress("mcp2_time",&m2t);t->SetBranchAddress("mcp1_peak",&m1p);
        t->SetBranchAddress("s_peak",sp);t->SetBranchAddress("s_cfd05",sc);
        long N=t->GetEntries();
        double xs=0,ys=0; long nw=0;
        for(long i=0;i<N&&nw<50000;++i){t->GetEntry(i); if(wc&&x>-100&&x<100){xs+=x;ys+=y;++nw;}}
        double xc=nw?xs/nw:0, yc=nw?ys/nw:0;
        std::vector<float> slg, tc;
        for(long i=0;i<N;++i){ t->GetEntry(i);
            if(!wc)continue; double dx=x-xc,dy=y-yc; if(dx*dx+dy*dy>=9.0)continue;
            if(m1p<200||m1p>750||m1t<-1e5)continue;
            double ds=0,us=0; int dn=0,un=0;
            for(int k=0;k<4;++k){int s=dw[k]; if(sp[s]>20&&sc[s]>-1e5){ds+=sc[s]-m1t;++dn;}}
            for(int k=0;k<4;++k){int s=up[k];double ref=upM2[k]?m2t:m1t; if(ref<-1e5)continue; if(sp[s]>20&&sc[s]>-1e5){us+=sc[s]-ref;++un;}}
            if(dn<1||un<1)continue;
            double E=0; for(int k=0;k<8;++k)E+=sp[en[k]];
            slg.push_back((float)E); tc.push_back(0.5f*(float)(ds/dn-us/un));
        }
        long nt=slg.size(); if(nt<2000){printf("%5.0f  (low stats)\n",Es[e]);fp->Close();continue;}
        // mean,sigma of SumLG
        std::vector<float> s2=slg; std::sort(s2.begin(),s2.end()); double md=s2[s2.size()/2];
        double a=0,a2=0; long n=0; for(float v:slg){if(std::fabs(v-md)<0.5*md){a+=v;a2+=v*v;++n;}} double mE=a/n,sE=std::sqrt(a2/n-mE*mE);
        double lo=mE-2*sE, hi=mE+2*sE, bw=(hi-lo)/9.0;
        double best=1e9,bc=0; long bn=0;
        for(int b=0;b<9;++b){ double blo=lo+b*bw,bhi=blo+bw; std::vector<float> vt;
            for(long i=0;i<nt;++i) if(slg[i]>=blo&&slg[i]<bhi) vt.push_back(tc[i]);
            if((long)vt.size()<500)continue; double s=coreSig(vt)*1000.0;
            if(s>0&&s<best){best=s; bc=blo+0.5*bw; bn=vt.size();} }
        printf("%5.0f %12.1f %10.0f %8.1f %8ld\n",Es[e],best,bc,100.0*bn/nt,bn);
        fp->Close();
    }
}

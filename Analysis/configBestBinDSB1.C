// configBestBinDSB1.C — best-contained-bin (DW-UP)/2 on the DSB1 processRun
// ntuples, to CHECK that this simple pipeline reproduces the ~27 ps headline
// when the same best-bin selection is applied (vs 58 ps all-fiducial).
//   root -l 'Analysis/configBestBinDSB1.C+'
#include "TFile.h"
#include "TTree.h"
#include <cstdio>
#include <vector>
#include <algorithm>
#include <cmath>
static double coreSig(std::vector<float> v){
    if(v.size()<200)return -1; std::sort(v.begin(),v.end()); double mu=v[v.size()/2],s=0.2;
    for(int it=0;it<5;++it){double a=0,a2=0;long n=0;for(float x:v)if(std::fabs(x-mu)<2.5*s){a+=x;a2+=x*x;++n;}
        if(n<50)break;mu=a/n;double var=a2/n-mu*mu;s=var>0?std::sqrt(var):s;} return s;}
void configBestBinDSB1(){
    double Es[6]={25,50,75,100,125,150};
    printf("\n=== DSB1 (processRun) : all-fiducial vs best-bin ===\n");
    printf("%5s %12s %12s %8s\n","E","allfid[ps]","bestbin[ps]","eff%");
    for(int e=0;e<6;++e){
        TFile* fp=TFile::Open(Form("Analysis/Output/%.0fGeV/ntuple.root",Es[e])); if(!fp||fp->IsZombie())continue;
        TTree* t=(TTree*)fp->Get("rad"); if(!t){fp->Close();continue;}
        Bool_t wc;Float_t x,y,slg,cfd[8];
        t->SetBranchAddress("wc_ok",&wc);t->SetBranchAddress("x_trk",&x);t->SetBranchAddress("y_trk",&y);
        t->SetBranchAddress("sum_lg",&slg);t->SetBranchAddress("hg_cfd05",cfd);
        long N=t->GetEntries(); double xs=0,ys=0;long nw=0;
        for(long i=0;i<N&&nw<50000;++i){t->GetEntry(i);if(wc&&x>-100&&x<100){xs+=x;ys+=y;++nw;}}
        double xc=nw?xs/nw:0,yc=nw?ys/nw:0;
        std::vector<float> sl,tc;
        for(long i=0;i<N;++i){t->GetEntry(i); if(!wc)continue; double dx=x-xc,dy=y-yc; if(dx*dx+dy*dy>=9.0)continue;
            double ds=0,us=0;int dn=0,un=0; for(int k=0;k<4;++k)if(cfd[k]>-1e5){ds+=cfd[k];++dn;}
            for(int k=4;k<8;++k)if(cfd[k]>-1e5){us+=cfd[k];++un;} if(dn<1||un<1)continue;
            sl.push_back(slg); tc.push_back(0.5f*(float)(ds/dn-us/un));}
        long nt=sl.size(); if(nt<2000){fp->Close();continue;}
        double allf=coreSig(tc)*1000.0;
        std::vector<float> s2=sl;std::sort(s2.begin(),s2.end());double md=s2[s2.size()/2],a=0,a2=0;long n=0;
        for(float v:sl)if(std::fabs(v-md)<0.5*md){a+=v;a2+=v*v;++n;} double mE=a/n,sE=std::sqrt(a2/n-mE*mE);
        double lo=mE-2*sE,hi=mE+2*sE,bw=(hi-lo)/9.0,best=1e9;long bn=0;
        for(int b=0;b<9;++b){double blo=lo+b*bw,bhi=blo+bw;std::vector<float> vt;
            for(long i=0;i<nt;++i)if(sl[i]>=blo&&sl[i]<bhi)vt.push_back(tc[i]);
            if((long)vt.size()<500)continue;double s=coreSig(vt)*1000.0;if(s>0&&s<best){best=s;bn=vt.size();}}
        printf("%5.0f %12.1f %12.1f %8.1f\n",Es[e],allf,best,100.0*bn/nt);
        fp->Close();
    }
}

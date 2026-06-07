// configResolutionDSB1.C — same simple method as configResolution.C, but reading
// the original DSB1 processRun ntuples (output/<E>GeV/ntuple.root),
// where hg_cfd05[8] is ALREADY MCP-referenced.  Gives the DSB1 baseline in the
// identical method for a fair cross-config comparison.
//   root -l 'Analysis/configResolutionDSB1.C+'
#include "TFile.h"
#include "TTree.h"
#include <cstdio>
#include <vector>
#include <algorithm>
#include <cmath>

static double coreSigma(std::vector<float> v){
    if(v.size()<200) return -1;
    std::sort(v.begin(),v.end()); double mu=v[v.size()/2], s=0.2;
    for(int it=0;it<5;++it){ double a=0,a2=0; long n=0;
        for(float x:v){ if(std::fabs(x-mu)<2.5*s){a+=x;a2+=x*x;++n;} }
        if(n<50)break; mu=a/n; double var=a2/n-mu*mu; s=var>0?std::sqrt(var):s; }
    return s;
}

void configResolutionDSB1()
{
    double Es[6]={25,50,75,100,125,150};
    printf("\n=== DSB1 (processRun ntuples, same simple method) ===\n");
    printf("%5s %9s %10s %10s %8s\n","E","N_fid","sigma_t[ps]","sE/E[%]","<SumLG>");
    for(int e=0;e<6;++e){
        TFile* fp=TFile::Open(Form("output/%.0fGeV/ntuple.root",Es[e]));
        if(!fp||fp->IsZombie())continue; TTree* t=(TTree*)fp->Get("rad"); if(!t){fp->Close();continue;}
        Bool_t wc; Float_t x,y,slg,cfd[8];
        t->SetBranchAddress("wc_ok",&wc);t->SetBranchAddress("x_trk",&x);t->SetBranchAddress("y_trk",&y);
        t->SetBranchAddress("sum_lg",&slg);t->SetBranchAddress("hg_cfd05",cfd);
        long N=t->GetEntries();
        double xs=0,ys=0; long nw=0;
        for(long i=0;i<N&&nw<50000;++i){t->GetEntry(i); if(wc&&x>-100&&x<100){xs+=x;ys+=y;++nw;}}
        double xc=nw?xs/nw:0, yc=nw?ys/nw:0;
        std::vector<float> vt, ve; long nf=0;
        for(long i=0;i<N;++i){ t->GetEntry(i);
            if(!wc)continue; double dx=x-xc,dy=y-yc; if(dx*dx+dy*dy>=9.0)continue; ++nf;
            double ds=0,us=0; int dn=0,un=0;
            for(int k=0;k<4;++k) if(cfd[k]>-1e5){ds+=cfd[k];++dn;}
            for(int k=4;k<8;++k) if(cfd[k]>-1e5){us+=cfd[k];++un;}
            if(dn>=1&&un>=1) vt.push_back(0.5f*(float)(ds/dn-us/un));
            ve.push_back(slg);
        }
        double st=coreSigma(vt)*1000.0, sE=0,mE=0;
        { std::vector<float> tmp=ve; std::sort(tmp.begin(),tmp.end()); double md=tmp[tmp.size()/2],a=0,a2=0;long n=0;
          for(float v:tmp) if(std::fabs(v-md)<0.5*md){a+=v;a2+=v*v;++n;} if(n>50){mE=a/n;sE=std::sqrt(a2/n-mE*mE);} }
        printf("%5.0f %9ld %10.1f %10.1f %8.0f\n",Es[e],nf,st,100.0*sE/mE,mE);
        fp->Close();
    }
}

// dsb1OOSandCorr.C — (a) DSB1 HG best-bin OOS in the SAME method as the reduced
// configs (processRun ntuples; run-folded), for apples-to-apples; and (b) a
// viability check for the LG->HG de-saturation idea: per-event correlation of
// hg_peak vs lg_peak in the UNSATURATED regime (low E), per channel.
//   root -l 'Analysis/dsb1OOSandCorr.C+'
#include "TFile.h"
#include "DataPaths.h"
#include "TTree.h"
#include <cstdio>
#include <vector>
#include <algorithm>
#include <cmath>
struct Ev{int run;float slg;float t;};
static double coreS(std::vector<float>&v){ if(v.size()<200)return -1;
    std::sort(v.begin(),v.end());double mu=v[v.size()/2],s=0.3;
    for(int it=0;it<5;++it){double a=0,a2=0;long n=0;for(float x:v)if(std::fabs(x-mu)<2.5*s){a+=x;a2+=x*x;++n;}
        if(n<40)break;mu=a/n;double var=a2/n-mu*mu;s=var>0?std::sqrt(var):s;}return s;}
static double oosBest(std::vector<Ev>&ev,double&insamp){ insamp=-1; if(ev.size()<3000)return -1;
    std::vector<float> sl;for(auto&e:ev)sl.push_back(e.slg);std::sort(sl.begin(),sl.end());
    double md=sl[sl.size()/2],a=0,a2=0;long n=0;for(float v:sl)if(std::fabs(v-md)<0.5*md){a+=v;a2+=v*v;++n;}
    double mE=a/n,sE=std::sqrt(a2/n-mE*mE),lo=mE-2*sE,bw=4*sE/9.0;
    {double best=1e9;for(int b=0;b<9;++b){double blo=lo+b*bw,bhi=blo+bw;std::vector<float>vt;
        for(auto&e:ev)if(e.slg>=blo&&e.slg<bhi)vt.push_back(e.t);if((long)vt.size()<500)continue;
        double s=coreS(vt)*1000;if(s>0&&s<best)best=s;}insamp=best;}
    double acc=0;int nf=0;for(int f=0;f<2;++f){int bb=-1;double btr=1e9;
        for(int b=0;b<9;++b){double blo=lo+b*bw,bhi=blo+bw;std::vector<float>vt;
            for(auto&e:ev)if((e.run%2)==f&&e.slg>=blo&&e.slg<bhi)vt.push_back(e.t);if((long)vt.size()<400)continue;
            double s=coreS(vt)*1000;if(s>0&&s<btr){btr=s;bb=b;}}
        if(bb<0)continue;double blo=lo+bb*bw,bhi=blo+bw;std::vector<float>vt;
        for(auto&e:ev)if((e.run%2)!=f&&e.slg>=blo&&e.slg<bhi)vt.push_back(e.t);
        if((long)vt.size()<200)continue;double s=coreS(vt)*1000;if(s>0){acc+=s;++nf;}}
    return nf?acc/nf:-1;}

void dsb1OOSandCorr(){
    double Es[6]={25,50,75,100,125,150};
    printf("\n=== (a) DSB1 HG best-bin (DW-UP)/2 sigma_t [ps]: OOS (in-sample) ===\n");
    for(int e=0;e<6;++e){
        TFile* fp=TFile::Open(radReduced("DSB1",Es[e]));if(!fp||fp->IsZombie())continue;
        TTree* t=(TTree*)fp->Get("rad");if(!t){fp->Close();continue;}
        Int_t run;Bool_t wc;Float_t x,y,slg,cfd[8];
        t->SetBranchAddress("run",&run);t->SetBranchAddress("wc_ok",&wc);t->SetBranchAddress("x_trk",&x);t->SetBranchAddress("y_trk",&y);
        t->SetBranchAddress("sum_lg",&slg);t->SetBranchAddress("hg_cfd05",cfd);
        long N=t->GetEntries();double xs=0,ys=0;long nw=0;
        for(long i=0;i<N&&nw<40000;++i){t->GetEntry(i);if(wc&&x>-100&&x<100){xs+=x;ys+=y;++nw;}}
        double xc=nw?xs/nw:0,yc=nw?ys/nw:0;std::vector<Ev>hg;
        for(long i=0;i<N;++i){t->GetEntry(i);if(!wc)continue;double dx=x-xc,dy=y-yc;if(dx*dx+dy*dy>=9.0)continue;
            double ds=0,us=0;int dn=0,un=0;for(int k=0;k<4;++k)if(cfd[k]>-1e5){ds+=cfd[k];++dn;}
            for(int k=4;k<8;++k)if(cfd[k]>-1e5){us+=cfd[k];++un;}if(dn<1||un<1)continue;
            Ev ev;ev.run=run;ev.slg=slg;ev.t=0.5f*(float)(ds/dn-us/un);hg.push_back(ev);}
        double in,oo=oosBest(hg,in);printf("  E=%3.0f : %6.1f (%5.1f)\n",Es[e],oo,in);fp->Close();}

    // (b) per-event hg_peak vs lg_peak correlation in the UNSATURATED regime (25 & 50 GeV)
    printf("\n=== (b) LG->HG de-saturation viability: per-event hg_peak vs lg_peak ===\n");
    printf("    (unsaturated regime; correlation r and relative scatter of HG about the LG-fit)\n");
    for(double E : {25.,50.}){
        TFile* fp=TFile::Open(radReduced("DSB1",E));if(!fp||fp->IsZombie())continue;
        TTree* t=(TTree*)fp->Get("rad");if(!t){fp->Close();continue;}
        Bool_t wc;Float_t hg[8],lg[8];
        t->SetBranchAddress("wc_ok",&wc);t->SetBranchAddress("hg_peak",hg);t->SetBranchAddress("lg_peak",lg);
        // channel 2 (NE-D, a bright timing cap); accumulate for linear fit hg = a + b*lg
        double sx=0,sy=0,sxx=0,sxy=0,syy=0;long n=0;long N=t->GetEntries();
        for(long i=0;i<N;++i){t->GetEntry(i);if(!wc)continue;double X=lg[2],Y=hg[2];
            if(X<20||Y<20||Y>780)continue; sx+=X;sy+=Y;sxx+=X*X;sxy+=X*Y;syy+=Y*Y;++n;}
        if(n>500){double mx=sx/n,my=sy/n,cxy=sxy/n-mx*my,vx=sxx/n-mx*mx,vy=syy/n-my*my;
            double b=cxy/vx,r=cxy/std::sqrt(vx*vy);
            double resid2=vy-b*cxy; double rms=std::sqrt(resid2>0?resid2:0);
            printf("  E=%3.0f  NE-D: r=%.3f  slope=%.2f  HG scatter about LG-fit = %.1f mV (%.1f%% of <HG>=%.0f)\n",
                   E,r,b,rms,100*rms/my,my);}
        fp->Close();}
}

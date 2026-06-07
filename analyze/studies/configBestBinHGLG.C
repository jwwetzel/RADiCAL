// configBestBinHGLG.C — best-contained-bin (DW-UP)/2 timing, run-folded OOS, for
// the HG (fast/saturated, timing) vs LG (slow/clean, energy) capillary chains of
// a reduced config.  Apples-to-apples with the DSB1 headline method (best ΣLG
// bin), plus 2-fold run-parity OOS so the number is defensible (no overfitting).
//   root -l 'Analysis/configBestBinHGLG.C+("reduced/LUAG","LuAG")'
#include "ChannelConfig.h"
#include "TFile.h"
#include "TTree.h"
#include <cstdio>
#include <vector>
#include <algorithm>
#include <cmath>
struct Ev{int run;float slg;float t;};
static double coreS(std::vector<float>&v){ if(v.size()<200)return -1;
    std::sort(v.begin(),v.end()); double mu=v[v.size()/2],s=0.3;
    for(int it=0;it<5;++it){double a=0,a2=0;long n=0;for(float x:v)if(std::fabs(x-mu)<2.5*s){a+=x;a2+=x*x;++n;}
        if(n<40)break;mu=a/n;double var=a2/n-mu*mu;s=var>0?std::sqrt(var):s;} return s;}

static double oosBest(std::vector<Ev>&ev,double&insamp){
    insamp=-1; if(ev.size()<4000)return -1;
    std::vector<float> sl; for(auto&e:ev)sl.push_back(e.slg);
    std::sort(sl.begin(),sl.end()); double md=sl[sl.size()/2],a=0,a2=0;long n=0;
    for(float v:sl)if(std::fabs(v-md)<0.5*md){a+=v;a2+=v*v;++n;} double mE=a/n,sE=std::sqrt(a2/n-mE*mE);
    double lo=mE-2*sE,bw=(4*sE)/9.0;
    // in-sample best bin
    {double best=1e9;for(int b=0;b<9;++b){double blo=lo+b*bw,bhi=blo+bw;std::vector<float>vt;
        for(auto&e:ev)if(e.slg>=blo&&e.slg<bhi)vt.push_back(e.t);
        if((long)vt.size()<500)continue;double s=coreS(vt)*1000;if(s>0&&s<best)best=s;} insamp=best;}
    // 2-fold run-parity OOS
    double acc=0;int nf=0;
    for(int f=0;f<2;++f){ int bestb=-1;double btr=1e9;
        for(int b=0;b<9;++b){double blo=lo+b*bw,bhi=blo+bw;std::vector<float>vt;
            for(auto&e:ev)if((e.run%2)==f&&e.slg>=blo&&e.slg<bhi)vt.push_back(e.t);
            if((long)vt.size()<400)continue;double s=coreS(vt)*1000;if(s>0&&s<btr){btr=s;bestb=b;}}
        if(bestb<0)continue; double blo=lo+bestb*bw,bhi=blo+bw;std::vector<float>vt;
        for(auto&e:ev)if((e.run%2)!=f&&e.slg>=blo&&e.slg<bhi)vt.push_back(e.t);
        if((long)vt.size()<200)continue; double s=coreS(vt)*1000; if(s>0){acc+=s;++nf;}}
    return nf?acc/nf:-1;
}

void configBestBinHGLG(const char* dir,const char* label=""){
    int hgD[4],hgU[4];bool m2[4];int lgD[4],lgU[4],lgAll[8];
    for(int i=0;i<4;++i)hgD[i]=kCap[i].hg/1024;
    for(int i=4;i<8;++i){hgU[i-4]=kCap[i].hg/1024;m2[i-4]=kCap[i].use_mcp2;}
    for(int i=0;i<4;++i)lgD[i]=kCap[i].lg/1024;
    for(int i=4;i<8;++i)lgU[i-4]=kCap[i].lg/1024;
    for(int i=0;i<8;++i)lgAll[i]=kCap[i].lg/1024;
    double Es[6]={25,50,75,100,125,150};
    printf("\n=== %s : best-bin (DW-UP)/2 sigma_t [ps], OOS (in-sample) ===\n",label[0]?label:dir);
    printf("%5s %18s %18s\n","E","HG (timing chain)","LG (energy chain)");
    for(int e=0;e<6;++e){
        TFile* fp=TFile::Open(Form("%s/%.0fGeV.root",dir,Es[e]));if(!fp||fp->IsZombie())continue;
        TTree* t=(TTree*)fp->Get("rad");if(!t){fp->Close();continue;}
        Int_t run;Bool_t wc;Float_t x,y,m1t,m2t,m1p,sp[36],sc[36];
        t->SetBranchAddress("run",&run);t->SetBranchAddress("wc_ok",&wc);t->SetBranchAddress("x_trk",&x);t->SetBranchAddress("y_trk",&y);
        t->SetBranchAddress("mcp1_time",&m1t);t->SetBranchAddress("mcp2_time",&m2t);t->SetBranchAddress("mcp1_peak",&m1p);
        t->SetBranchAddress("s_peak",sp);t->SetBranchAddress("s_cfd05",sc);
        long N=t->GetEntries();double xs=0,ys=0;long nw=0;
        for(long i=0;i<N&&nw<40000;++i){t->GetEntry(i);if(wc&&x>-100&&x<100){xs+=x;ys+=y;++nw;}}
        double xc=nw?xs/nw:0,yc=nw?ys/nw:0;
        std::vector<Ev> hg,lg;
        for(long i=0;i<N;++i){t->GetEntry(i);if(!wc)continue;double dx=x-xc,dy=y-yc;if(dx*dx+dy*dy>=9.0)continue;
            if(m1p<200||m1p>750||m1t<-1e5)continue;
            double slg=0;for(int k=0;k<8;++k)slg+=sp[lgAll[k]];
            // HG (DW-UP)/2, MCP-referenced
            {double ds=0,us=0;int dn=0,un=0;
             for(int k=0;k<4;++k){int s=hgD[k];if(sp[s]>20&&sc[s]>-1e5){ds+=sc[s]-m1t;++dn;}}
             for(int k=0;k<4;++k){int s=hgU[k];double r=m2[k]?m2t:m1t;if(r<-1e5)continue;if(sp[s]>20&&sc[s]>-1e5){us+=sc[s]-r;++un;}}
             if(dn>=1&&un>=1){Ev ev;ev.run=run;ev.slg=(float)slg;ev.t=0.5f*(float)(ds/dn-us/un);hg.push_back(ev);} }
            // LG (DW-UP)/2, same DRS1 G0 group -> no MCP needed
            {double ds=0,us=0;int dn=0,un=0;
             for(int k=0;k<4;++k){int s=lgD[k];if(sp[s]>15&&sc[s]>-1e5){ds+=sc[s];++dn;}}
             for(int k=0;k<4;++k){int s=lgU[k];if(sp[s]>15&&sc[s]>-1e5){us+=sc[s];++un;}}
             if(dn>=1&&un>=1){Ev ev;ev.run=run;ev.slg=(float)slg;ev.t=0.5f*(float)(ds/dn-us/un);lg.push_back(ev);} }
        }
        double hi,li,ho=oosBest(hg,hi),lo2=oosBest(lg,li);
        printf("%5.0f   %7.1f (%5.1f)      %7.1f (%6.1f)\n",Es[e],ho,hi,lo2,li);
        fp->Close();
    }
}

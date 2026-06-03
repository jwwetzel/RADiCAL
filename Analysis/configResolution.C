// configResolution.C — first-look timing + energy resolution for a reduced
// config dataset (reduceRaw.C output), reusing the DSB1 corner/energy geometry.
//   timing : (DW-UP)/2 corner estimator, CFD-5%, MCP-referenced (cancels in diff)
//   energy : Sigma LG-peak sum  ->  sigma_E/E
// Assumes the config shares the DSB1 wiring (verified by discoverReduced.C).
//   root -l 'Analysis/configResolution.C+("reduced/LUAG")'
#include "TFile.h"
#include "TTree.h"
#include "ChannelConfig.h"
#include <cstdio>
#include <vector>
#include <algorithm>
#include <cmath>

static double coreSigma(std::vector<float>& v, double win){ // iterative truncated RMS [same units]
    if(v.size()<200) return -1;
    std::sort(v.begin(),v.end());
    double mu=v[v.size()/2], s=win;
    for(int it=0;it<5;++it){ double a=0,a2=0; long n=0;
        for(float x:v){ double d=x-mu; if(std::fabs(d)<2.5*s){a+=x;a2+=x*x;++n;} }
        if(n<50) break; mu=a/n; double var=a2/n-mu*mu; s=var>0?std::sqrt(var):s; }
    return s;
}

void configResolution(const char* dir, const char* label="")
{
    // slot indices from the DSB1 map
    int dw[4], up[4]; bool upM2[4]; int en[8];
    for(int i=0;i<4;++i) dw[i]=kCap[i].hg/1024;
    for(int i=4;i<8;++i){ up[i-4]=kCap[i].hg/1024; upM2[i-4]=kCap[i].use_mcp2; }
    for(int i=0;i<8;++i) en[i]=kCap[i].lg/1024;

    double Es[6]={25,50,75,100,125,150};
    printf("\n=== %s ===\n", label[0]?label:dir);
    printf("%5s %9s %10s %10s %8s\n","E","N_fid","sigma_t[ps]","sE/E[%]","<SumLG>");
    for(int e=0;e<6;++e){
        TString f=Form("%s/%.0fGeV.root",dir,Es[e]);
        TFile* fp=TFile::Open(f); if(!fp||fp->IsZombie()){ continue; }
        TTree* t=(TTree*)fp->Get("rad"); if(!t){ fp->Close(); continue; }
        Bool_t wc; Float_t x,y,m1t,m2t,m1p; Float_t sp[36],sc[36];
        t->SetBranchAddress("wc_ok",&wc); t->SetBranchAddress("x_trk",&x); t->SetBranchAddress("y_trk",&y);
        t->SetBranchAddress("mcp1_time",&m1t); t->SetBranchAddress("mcp2_time",&m2t);
        t->SetBranchAddress("mcp1_peak",&m1p);
        t->SetBranchAddress("s_peak",sp); t->SetBranchAddress("s_cfd05",sc);
        long N=t->GetEntries();
        // beam centroid from wc_ok
        double xs=0,ys=0; long nw=0;
        for(long i=0;i<N && nw<50000;++i){ t->GetEntry(i); if(wc&&x>-100&&x<100){xs+=x;ys+=y;++nw;} }
        double xc=nw?xs/nw:0, yc=nw?ys/nw:0;
        std::vector<float> vt, ve; long nfid=0;
        for(long i=0;i<N;++i){ t->GetEntry(i);
            if(!wc) continue; double dx=x-xc,dy=y-yc; if(dx*dx+dy*dy>=9.0) continue;  // r<3mm
            if(m1p<200||m1p>750) continue; if(m1t<-1e5) continue;
            ++nfid;
            // timing
            double dsum=0,usum=0; int dn=0,un=0;
            for(int k=0;k<4;++k){ int s=dw[k]; if(sp[s]>20&&sc[s]>-1e5){ dsum+=sc[s]-m1t; ++dn; } }
            for(int k=0;k<4;++k){ int s=up[k]; double ref=upM2[k]?m2t:m1t; if(ref<-1e5)continue;
                                  if(sp[s]>20&&sc[s]>-1e5){ usum+=sc[s]-ref; ++un; } }
            if(dn>=1&&un>=1) vt.push_back(0.5f*(float)(dsum/dn-usum/un));
            // energy
            double slg=0; for(int k=0;k<8;++k) slg+=sp[en[k]];
            ve.push_back((float)slg);
        }
        double st=coreSigma(vt,0.2)*1000.0;        // ns->ps
        // energy: mean + core sigma
        std::vector<float> ve2=ve; double sE=coreSigma(ve2,0.0);
        double meanE=0; for(float v:ve) meanE+=v; meanE/= (ve.empty()?1:ve.size());
        // recompute sE about mean properly
        { std::vector<float> tmp=ve; std::sort(tmp.begin(),tmp.end()); double md=tmp[tmp.size()/2],a=0,a2=0; long n=0;
          for(float v:tmp){ if(std::fabs(v-md)<0.4*md){a+=v;a2+=v*v;++n;} } if(n>50){ double mu=a/n; sE=std::sqrt(a2/n-mu*mu); meanE=mu; } }
        printf("%5.0f %9ld %10.1f %10.1f %8.0f\n", Es[e], nfid, st, 100.0*sE/meanE, meanE);
        fp->Close();
    }
}

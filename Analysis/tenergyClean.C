// tenergyClean.C — TENERGY contains an E-type ("energy") capillary at the NW
//   corner (logbook "3xDSB1, 1xEnergy" / "E,52,54,57"; data: NW ~4-5x dimmer).
//   Our (DW-UP)/2 estimator wrongly averages that energy cap in as a timing
//   channel. This recomputes TENERGY's OOS best-bin sigma_t (a) with all 4 corners
//   (contaminated) and (b) EXCLUDING the NW corner (clean 3xDSB1 timing), and
//   shows DSB1 for reference.
//   root -l 'Analysis/tenergyClean.C+'
#include "ChannelConfig.h"
#include "PlotUtils.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include <cstdio>
#include <vector>
#include <algorithm>
#include <cmath>

struct Ev { int run; float slg; float t; };
static double gcoreSig(std::vector<float>& v){
    if(v.size()<150) return -1;
    std::vector<float> s=v; std::sort(s.begin(),s.end()); double md=s[s.size()/2];
    TH1F h("_g","",240,md-1.5,md+1.5); for(float x:v) h.Fill(x);
    double mu,muE,sg,sgE; FitGaussCore(&h,2.0,mu,muE,sg,sgE); return sg>0? sg*1000.0:-1;
}
static double oosBest(std::vector<Ev>& ev){
    if(ev.size()<2500) return -1;
    std::vector<float> sl; for(auto&e:ev) sl.push_back(e.slg); std::sort(sl.begin(),sl.end());
    double md=sl[sl.size()/2],a=0,a2=0; long n=0;
    for(float v:sl) if(std::fabs(v-md)<0.5*md){a+=v;a2+=v*v;++n;}
    double mE=a/n, sE=std::sqrt(a2/n-mE*mE), lo=mE-2*sE, bw=4*sE/9.0;
    std::vector<int> ur; for(auto&e:ev) ur.push_back(e.run); std::sort(ur.begin(),ur.end());
    ur.erase(std::unique(ur.begin(),ur.end()),ur.end());
    const int kF=5; auto fold=[&](int r){ int idx=(int)(std::lower_bound(ur.begin(),ur.end(),r)-ur.begin()); return idx%kF; };
    const int kTrainMin=(500*(kF-1))/kF; std::vector<float> pool;
    for(int f=0;f<kF;++f){ int bb=-1; double btr=1e9;
        for(int b=0;b<9;++b){ double blo=lo+b*bw,bhi=blo+bw; std::vector<float> vt;
            for(auto&e:ev) if(fold(e.run)!=f && e.slg>=blo && e.slg<bhi) vt.push_back(e.t);
            if((long)vt.size()<kTrainMin) continue; double s=gcoreSig(vt); if(s>0&&s<btr){btr=s;bb=b;} }
        if(bb<0) continue; double blo=lo+bb*bw,bhi=blo+bw;
        for(auto&e:ev) if(fold(e.run)==f && e.slg>=blo && e.slg<bhi) pool.push_back(e.t); }
    return pool.size()<150? -1 : gcoreSig(pool);
}
// excl = corner to drop (-1 none; 0 = NW). Drops that corner's DOWN and UP readout.
static void loadCfg(const char* dir, bool isDSB1, double E, int excl, std::vector<Ev>& ev){
    TString fn = isDSB1 ? Form("Analysis/Output/%.0fGeV/ntuple.root",E) : Form("%s/%.0fGeV.root",dir,E);
    TFile* fp=TFile::Open(fn); if(!fp||fp->IsZombie()) return;
    TTree* t=(TTree*)fp->Get("rad"); if(!t){ fp->Close(); return; }
    Int_t run; Bool_t wc; Float_t x,y,m1t,m2t,m1p,sp[36],sc[36],mp,cfd[8],slgD;
    t->SetBranchAddress("wc_ok",&wc); t->SetBranchAddress("x_trk",&x); t->SetBranchAddress("y_trk",&y);
    if(isDSB1){ t->SetBranchAddress("mcp_peak",&mp); t->SetBranchAddress("hg_cfd05",cfd); t->SetBranchAddress("sum_lg",&slgD); }
    else { t->SetBranchAddress("mcp1_time",&m1t); t->SetBranchAddress("mcp2_time",&m2t);
           t->SetBranchAddress("mcp1_peak",&m1p); t->SetBranchAddress("s_peak",sp); t->SetBranchAddress("s_cfd05",sc); }
    int dw[4],up[4]; bool upM2[4]; int en[8];
    for(int i=0;i<4;++i) dw[i]=kCap[i].hg/1024;
    for(int i=4;i<8;++i){ up[i-4]=kCap[i].hg/1024; upM2[i-4]=kCap[i].use_mcp2; }
    for(int i=0;i<8;++i) en[i]=kCap[i].lg/1024;
    long N=t->GetEntries(); double xs=0,ys=0; long nw=0;
    for(long i=0;i<N&&nw<50000;++i){ t->GetEntry(i); if(wc&&x>-100&&x<100){xs+=x;ys+=y;++nw;} }
    double xc=nw?xs/nw:0, yc=nw?ys/nw:0;
    for(long i=0;i<N;++i){ t->GetEntry(i);
        if(!wc) continue; double dx=x-xc,dy=y-yc; if(dx*dx+dy*dy>=9.0) continue;
        if(isDSB1){ if(mp<200||mp>750) continue; } else { if(m1p<200||m1p>750||m1t<-1e5) continue; }
        double dsum=0,usum=0; int dn=0,un=0, slg=0; double slgv=0;
        for(int k=0;k<4;++k){ if(k==excl) continue;
            if(isDSB1){ if(cfd[k]>-1e5){ dsum+=cfd[k]; ++dn; } }
            else { int s=dw[k]; if(sp[s]>20&&sc[s]>-1e5){ dsum+=sc[s]-m1t; ++dn; } } }
        for(int k=0;k<4;++k){ if(k==excl) continue;
            if(isDSB1){ if(cfd[k+4]>-1e5){ usum+=cfd[k+4]; ++un; } }
            else { int s=up[k]; double ref=upM2[k]?m2t:m1t; if(ref<-1e5)continue;
                   if(sp[s]>20&&sc[s]>-1e5){ usum+=sc[s]-ref; ++un; } } }
        if(isDSB1) slgv=slgD; else { for(int k=0;k<8;++k) slgv+=sp[en[k]]; }
        if(dn>=1&&un>=1){ Ev e; e.run=run; e.slg=(float)slgv; e.t=0.5f*(float)(dsum/dn-usum/un); ev.push_back(e); }
    }
    fp->Close();
}
void tenergyClean(){
    double Es[6]={25,50,75,100,125,150};
    printf("\n=== TENERGY timing: all-4-corner (contaminated by NW E-type) vs NW-excluded (clean 3xDSB1) ===\n");
    printf("%5s %14s %16s %12s\n","E","TEN_all4[ps]","TEN_exclNW[ps]","DSB1[ps]");
    for(int e=0;e<6;++e){
        std::vector<Ev> a,b,d;
        loadCfg("reduced/TENERGY",false,Es[e],-1,a);
        loadCfg("reduced/TENERGY",false,Es[e], 0,b);   // drop NW corner
        loadCfg("",true,Es[e],-1,d);
        double sa=oosBest(a), sb=oosBest(b), sd=oosBest(d);
        printf("%5.0f %14.1f %16.1f %12.1f\n",Es[e],sa,sb,sd);
    }
    printf("(NW corner = the E-type 'energy' capillary; excluding it gives a clean 3xDSB1-timing estimate)\n");
}

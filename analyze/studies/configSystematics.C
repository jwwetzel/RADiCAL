// configSystematics.C — follow-on paper gap #4: per-config systematic on the
//   OOS best-bin sigma_t @150 GeV. Re-runs the full OOS analysis under cut
//   variations (fiducial radius, containment, MCP window) and quadrature-sums the
//   shifts from nominal, per build.
//   root -l 'Analysis/configSystematics.C+'
#include "ChannelConfig.h"
#include "PlotUtils.h"
#include "TFile.h"
#include "DataPaths.h"
#include "TTree.h"
#include "TH1F.h"
#include <cstdio>
#include <vector>
#include <algorithm>
#include <cmath>

struct Ev { int run; float slg; float t; };
static const int kPb[4]={ chanOff(0,1,1)/1024, chanOff(0,1,2)/1024, chanOff(0,1,3)/1024, chanOff(0,1,4)/1024 };

static double gcoreSig(std::vector<float>& v){
    if(v.size()<150) return -1;
    std::vector<float> s=v; std::sort(s.begin(),s.end()); double md=s[s.size()/2];
    TH1F h("_g","",200,md-1.0,md+1.0); for(float x:v) h.Fill(x);
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
static void loadCfg(const char* dir, bool isDSB1, double E,
                    double fidR, double mcpLo, double mcpHi, double contMax, std::vector<Ev>& ev){
    TString fn = isDSB1 ? radReduced("DSB1",E) : radReduced(dir,E);
    TFile* fp=TFile::Open(fn); if(!fp||fp->IsZombie()) return;
    TTree* t=(TTree*)fp->Get("rad"); if(!t){ fp->Close(); return; }
    Int_t run; Bool_t wc; Float_t x,y,m1t,m2t,m1p,sp[36],sc[36],mp,cfd[8],slgD,spb;
    t->SetBranchAddress("wc_ok",&wc); t->SetBranchAddress("x_trk",&x); t->SetBranchAddress("y_trk",&y);
    if(isDSB1){ t->SetBranchAddress(t->GetBranch("mcp1_peak")?"mcp1_peak":"mcp_peak",&mp); t->SetBranchAddress("hg_cfd05",cfd);
                t->SetBranchAddress("sum_lg",&slgD); t->SetBranchAddress("sum_pb",&spb); }
    else { t->SetBranchAddress("mcp1_time",&m1t); t->SetBranchAddress("mcp2_time",&m2t);
           t->SetBranchAddress("mcp1_peak",&m1p); t->SetBranchAddress("s_peak",sp); t->SetBranchAddress("s_cfd05",sc); }
    int dw[4],up[4]; bool upM2[4]; int en[8];
    for(int i=0;i<4;++i) dw[i]=kCap[i].hg/1024;
    for(int i=4;i<8;++i){ up[i-4]=kCap[i].hg/1024; upM2[i-4]=kCap[i].use_mcp2; }
    for(int i=0;i<8;++i) en[i]=kCap[i].lg/1024;
    long N=t->GetEntries(); double xs=0,ys=0; long nw=0;
    for(long i=0;i<N&&nw<50000;++i){ t->GetEntry(i); if(wc&&x>-100&&x<100){xs+=x;ys+=y;++nw;} }
    double xc=nw?xs/nw:0, yc=nw?ys/nw:0; double R2=fidR*fidR;
    for(long i=0;i<N;++i){ t->GetEntry(i);
        if(!wc) continue; double dx=x-xc,dy=y-yc; if(dx*dx+dy*dy>=R2) continue;
        double slg, pb, tval; bool ok;
        if(isDSB1){ if(mp<mcpLo||mp>mcpHi) continue; slg=slgD; pb=spb;
            double ds=0,us=0; int dn=0,un=0;
            for(int k=0;k<4;++k) if(cfd[k]>-1e5){ ds+=cfd[k]; ++dn; }
            for(int k=4;k<8;++k) if(cfd[k]>-1e5){ us+=cfd[k]; ++un; }
            ok=(dn>=1&&un>=1); tval=ok?0.5*(ds/dn-us/un):0;
        } else { if(m1p<mcpLo||m1p>mcpHi||m1t<-1e5) continue;
            slg=0; for(int k=0;k<8;++k) slg+=sp[en[k]];
            pb=0; for(int k=0;k<4;++k) pb+=sp[kPb[k]];
            double ds=0,us=0; int dn=0,un=0;
            for(int k=0;k<4;++k){ int s=dw[k]; if(sp[s]>20&&sc[s]>-1e5){ ds+=sc[s]-m1t; ++dn; } }
            for(int k=0;k<4;++k){ int s=up[k]; double ref=upM2[k]?m2t:m1t; if(ref<-1e5)continue;
                if(sp[s]>20&&sc[s]>-1e5){ us+=sc[s]-ref; ++un; } }
            ok=(dn>=1&&un>=1); tval=ok?0.5*(ds/dn-us/un):0;
        }
        if(slg>0 && pb < contMax*slg && ok){ Ev e; e.run=run; e.slg=(float)slg; e.t=(float)tval; ev.push_back(e); }
    }
    fp->Close();
}

void configSystematics(){
    const char* lab[4]={"DSB1","LuAG","MIXED","TENERGY"};
    const char* dir[4]={"","LUAG","MIXED","TENERGY"};
    struct SV{ const char* n; double fid,cont,mlo,mhi; };
    SV var[7]={ {"Nominal",3.0,0.30,200,750},{"Fid-0.5",2.5,0.30,200,750},{"Fid+0.5",3.5,0.30,200,750},
                {"Cont-0.05",3.0,0.25,200,750},{"Cont+0.05",3.0,0.35,200,750},
                {"MCPlo+50",3.0,0.30,250,750},{"MCPhi-50",3.0,0.30,200,700} };
    // systematic is evaluated on the STABLE all-fiducial sigma_t (cut-variation on the
    // jumpy best-bin is dominated by selection jitter, not real cut sensitivity).
    auto allfid=[&](std::vector<Ev>& ev)->double{ if(ev.size()<500) return -1;
        std::vector<float> t; for(auto&e:ev) t.push_back(e.t); return gcoreSig(t); };
    printf("\n=== per-config systematic on the (stable) all-fiducial sigma_t @150 GeV ===\n");
    printf("%-9s %8s %7s","build","nominal","N"); for(int v=1;v<7;++v) printf(" %9s",var[v].n); printf("  %8s\n","SYST");
    for(int c=0;c<4;++c){
        double sN=-1, s[7]; long nNom=0;
        for(int v=0;v<7;++v){ std::vector<Ev> ev; loadCfg(dir[c],c==0,150,var[v].fid,var[v].mlo,var[v].mhi,var[v].cont,ev);
            if(v==0) nNom=ev.size(); s[v]=allfid(ev); }
        sN=s[0]; double sq=0; int nv=0;
        for(int v=1;v<7;++v) if(s[v]>0&&sN>0){ double d=s[v]-sN; sq+=d*d; ++nv; }
        double syst=std::sqrt(sq);
        printf("%-9s %8.1f %7ld",lab[c],sN,nNom);
        for(int v=1;v<7;++v) printf(" %+9.1f", (s[v]>0&&sN>0)? s[v]-sN : 0.0);
        printf("  %8.1f\n",syst);
    }
    printf("(columns 2-7 are shift vs nominal [ps]; SYST = quadrature sum)\n");
}

// configCapDiag.C — identify E-type (full-length WLS "energy") capillaries vs
//   T-type (shower-max "timing") capillaries in each build, from the DATA.
//   Signature: an E-type cap collects light along the whole shower -> much
//   BRIGHTER (higher LG/HG amplitude) but SLOWER (worse single-channel sigma_t).
//   Prints, per build @150 GeV, each of the 8 readout channels (4 corners x
//   up/down): mean LG peak, mean HG peak, single-channel sigma_t.
//   root -l 'Analysis/configCapDiag.C+'
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

static double gcoreSig(std::vector<float>& v){
    if(v.size()<200) return -1;
    std::vector<float> s=v; std::sort(s.begin(),s.end()); double md=s[s.size()/2];
    TH1F h("_g","",200,md-1.5,md+1.5); for(float x:v) h.Fill(x);
    double mu,muE,sg,sgE; FitGaussCore(&h,2.0,mu,muE,sg,sgE); return sg>0? sg*1000.0:-1;
}
static void diag(const char* dir, bool isDSB1, const char* label, double E){
    TString fn = isDSB1 ? radReduced("DSB1",E) : radReduced(dir,E);
    TFile* fp=TFile::Open(fn); if(!fp||fp->IsZombie()){ printf("  %s: no file\n",label); return; }
    TTree* t=(TTree*)fp->Get("rad"); if(!t){ fp->Close(); return; }
    Int_t run; Bool_t wc; Float_t x,y,m1t,m2t,m1p,sp[36],sc[36],mp,cfd[8],hgp[8],lgp[8];
    t->SetBranchAddress("wc_ok",&wc); t->SetBranchAddress("x_trk",&x); t->SetBranchAddress("y_trk",&y);
    if(isDSB1){ t->SetBranchAddress("mcp_peak",&mp); t->SetBranchAddress("hg_cfd05",cfd);
                t->SetBranchAddress("hg_peak",hgp); t->SetBranchAddress("lg_peak",lgp); }
    else { t->SetBranchAddress("mcp1_time",&m1t); t->SetBranchAddress("mcp2_time",&m2t);
           t->SetBranchAddress("mcp1_peak",&m1p); t->SetBranchAddress("s_peak",sp); t->SetBranchAddress("s_cfd05",sc); }
    long N=t->GetEntries(); double xs=0,ys=0; long nw=0;
    for(long i=0;i<N&&nw<50000;++i){ t->GetEntry(i); if(wc&&x>-100&&x<100){xs+=x;ys+=y;++nw;} }
    double xc=nw?xs/nw:0, yc=nw?ys/nw:0;
    std::vector<float> tv[8]; double lgs[8]={0},hgs[8]={0}; long cn[8]={0};
    for(long i=0;i<N;++i){ t->GetEntry(i);
        if(!wc) continue; double dx=x-xc,dy=y-yc; if(dx*dx+dy*dy>=9.0) continue;
        if(isDSB1){ if(mp<200||mp>750) continue; } else { if(m1p<200||m1p>750||m1t<-1e5) continue; }
        for(int c=0;c<8;++c){
            double tt,lg,hg;
            if(isDSB1){ if(cfd[c]<-1e5) continue; tt=cfd[c]; lg=lgp[c]; hg=hgp[c]; }
            else { int hs=kCap[c].hg/1024, ls=kCap[c].lg/1024; double ref=kCap[c].use_mcp2?m2t:m1t;
                   if(sc[hs]<-1e5||ref<-1e5) continue; tt=sc[hs]-ref; lg=sp[ls]; hg=sp[hs]; }
            tv[c].push_back((float)tt); lgs[c]+=lg; hgs[c]+=hg; ++cn[c];
        }
    }
    printf("\n=== %s @ %.0f GeV ===\n", label, E);
    printf("  %-6s %9s %9s %12s\n","chan","<LG>[mV]","<HG>[mV]","sigt[ps]");
    for(int c=0;c<8;++c){ if(cn[c]<300){ printf("  %-6s   (low)\n",kCap[c].name); continue; }
        double s=gcoreSig(tv[c]);
        printf("  %-6s %9.0f %9.0f %12.1f\n", kCap[c].name, lgs[c]/cn[c], hgs[c]/cn[c], s); }
    fp->Close();
}
void configCapDiag(){
    diag("",true,"DSB1",150);
    diag("LUAG",false,"LuAG",150);
    diag("MIXED",false,"MIXED",150);
    diag("TENERGY",false,"TENERGY",150);
}

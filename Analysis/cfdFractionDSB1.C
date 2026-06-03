// cfdFractionDSB1.C — (DW-UP)/2 sigma_t vs CFD fraction on the SATURATED DSB1
// channels.  Tests where on the rising edge the timing is best: low fractions
// sit on the sharp early edge (well below the clip), high fractions climb toward
// the saturated peak.  If sigma_t rises with fraction, the saturated upper part
// of the pulse times WORSE -> you rely on the sharp low edge, and the "blown up"
// (clipped) top is not usable for timing.
//   root -l 'Analysis/cfdFractionDSB1.C+'
#include "TFile.h"
#include "TTree.h"
#include <cstdio>
#include <vector>
#include <algorithm>
#include <cmath>
static double coreSig(std::vector<float> v){
    if(v.size()<300)return -1; std::sort(v.begin(),v.end()); double mu=v[v.size()/2],s=0.2;
    for(int it=0;it<5;++it){double a=0,a2=0;long n=0;for(float x:v)if(std::fabs(x-mu)<2.5*s){a+=x;a2+=x*x;++n;}
        if(n<50)break;mu=a/n;double var=a2/n-mu*mu;s=var>0?std::sqrt(var):s;} return s;}
void cfdFractionDSB1(){
    const char* br[5]={"hg_cfd05","hg_cfd10","hg_cfd","hg_cfd30","hg_cfd50"};
    double frac[5]={5,10,20,30,50};
    double Es[3]={50,100,150};
    printf("\n=== DSB1 (DW-UP)/2 sigma_t [ps] vs CFD fraction (all-fiducial) ===\n");
    printf("%6s %8s %8s %8s %8s %8s\n","E\\frac","5%","10%","20%","30%","50%");
    for(int e=0;e<3;++e){
        TFile* fp=TFile::Open(Form("Analysis/Output/%.0fGeV/ntuple.root",Es[e]));if(!fp||fp->IsZombie())continue;
        TTree* t=(TTree*)fp->Get("rad");if(!t){fp->Close();continue;}
        Bool_t wc;Float_t x,y,cf[5][8];
        t->SetBranchAddress("wc_ok",&wc);t->SetBranchAddress("x_trk",&x);t->SetBranchAddress("y_trk",&y);
        for(int k=0;k<5;++k) t->SetBranchAddress(br[k],cf[k]);
        long N=t->GetEntries();double xs=0,ys=0;long nw=0;
        for(long i=0;i<N&&nw<40000;++i){t->GetEntry(i);if(wc&&x>-100&&x<100){xs+=x;ys+=y;++nw;}}
        double xc=nw?xs/nw:0,yc=nw?ys/nw:0;
        std::vector<float> vt[5];
        for(long i=0;i<N;++i){t->GetEntry(i);if(!wc)continue;double dx=x-xc,dy=y-yc;if(dx*dx+dy*dy>=9.0)continue;
            for(int k=0;k<5;++k){double ds=0,us=0;int dn=0,un=0;
                for(int c=0;c<4;++c)if(cf[k][c]>-1e5){ds+=cf[k][c];++dn;}
                for(int c=4;c<8;++c)if(cf[k][c]>-1e5){us+=cf[k][c];++un;}
                if(dn>=1&&un>=1)vt[k].push_back(0.5f*(float)(ds/dn-us/un));}}
        printf("%6.0f",Es[e]); for(int k=0;k<5;++k)printf("%8.1f",coreSig(vt[k])*1000.0); printf("\n");
        fp->Close();
    }
    printf("(fraction is %% of the [clipped] peak; 5%% = sharp low edge, 50%% climbs toward saturation)\n");
}

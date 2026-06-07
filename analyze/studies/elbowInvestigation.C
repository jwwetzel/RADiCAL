// ===========================================================================
// elbowInvestigation.C  —  diagnose the bimodal "elbow" in the per-channel
//   CFD-20% timing distributions on the Down capillaries (50 GeV NW-D/NE-D/SE-D).
//
//   (1) Is the elbow caused by poor selection?
//   (2) Can tighter cuts narrow the peaks?
//   (3) What is the origin?  (amplitude walk vs DRS4 satellite vs position)
//
// Sentinel-safe: a CFD/LED value is only used if finite and |t| < 60 ns.
// Robust statistics: median + truncated RMS within +/- 1.5 ns of the median,
//   so a handful of failed crossings cannot blow up the width.
//
// Run:  ROOT_INCLUDE_PATH=Analysis root -l -b -q 'Analysis/elbowInvestigation.C+'
// ===========================================================================
#include <cstdio>
#include <vector>
#include <algorithm>
#include <cmath>
#include "TFile.h"
#include "TTree.h"

namespace {

inline bool ok(float t){ return t > -60.f && t < 60.f; }   // sentinel-safe window

double median(std::vector<float> v){
    if(v.empty()) return 0;
    std::sort(v.begin(),v.end());
    size_t n=v.size(); return (n&1)? v[n/2] : 0.5*(v[n/2-1]+v[n/2]);
}
// truncated RMS within +/- win of a center
double trms(const std::vector<float>& v, double c, double win, double* nused=nullptr){
    double s=0,s2=0; int k=0;
    for(float x:v){ if(std::fabs(x-c)<=win){ s+=x; s2+=x*x; ++k; } }
    if(nused)*nused=k;
    if(k<2) return 0;
    double m=s/k; double var=s2/k-m*m; return var>0? std::sqrt(var):0;
}

} // namespace

void elbowInvestigation()
{
    const char* fn = "output/50GeV/ntuple.root";
    TFile f(fn);
    if (f.IsZombie()){ printf("cannot open %s\n", fn); return; }
    TTree* t = (TTree*)f.Get("rad");
    if (!t){ printf("no tree\n"); return; }

    Float_t hg_peak[8], hg_cfd[8], hg_cfd05[8], hg_cfd03[8], hg_cfd10[8],
            hg_cfd30[8], hg_cfd50[8], hg_led[8], hg_ped_rms[8];
    Bool_t  hg_sat[8], hg_spike[8], in_fid;
    Int_t   stopcell[4];
    Float_t x_trk, y_trk, mcp_peak;
    t->SetBranchAddress("hg_peak", hg_peak);
    t->SetBranchAddress("hg_cfd",  hg_cfd);
    t->SetBranchAddress("hg_cfd05",hg_cfd05);
    t->SetBranchAddress("hg_cfd03",hg_cfd03);
    t->SetBranchAddress("hg_cfd10",hg_cfd10);
    t->SetBranchAddress("hg_cfd30",hg_cfd30);
    t->SetBranchAddress("hg_cfd50",hg_cfd50);
    t->SetBranchAddress("hg_led",  hg_led);
    t->SetBranchAddress("hg_ped_rms",hg_ped_rms);
    t->SetBranchAddress("hg_saturated",hg_sat);
    t->SetBranchAddress("hg_spike",hg_spike);
    t->SetBranchAddress("in_fiducial",&in_fid);
    t->SetBranchAddress("stopcell", stopcell);
    t->SetBranchAddress("x_trk",&x_trk);
    t->SetBranchAddress("y_trk",&y_trk);
    t->SetBranchAddress("mcp_peak",&mcp_peak);

    const char* nm[8]={"NW-D","NE-D","SE-D","SW-D","NW-U","NE-U","SE-U","SW-U"};
    const Long64_t N=t->GetEntries();
    const int chans[4]={0,1,2,6};   // 3 Down (in the plot) + SE-U control

    for(int ci=0;ci<4;++ci){
        const int ch=chans[ci];
        printf("\n================  channel %d  (%s)  50 GeV  ================\n",ch,nm[ch]);

        // gather per-channel arrays under the standard selection (fiducial, valid CFD-20)
        std::vector<float> v20,v05,v03,v10,v30,v50,vled,vamp,vx,vy,vmcp;
        for(Long64_t e=0;e<N;++e){
            t->GetEntry(e);
            if(!in_fid) continue;
            if(!ok(hg_cfd[ch])) continue;
            if(hg_peak[ch]<=20) continue;
            v20.push_back(hg_cfd[ch]);
            v05.push_back(ok(hg_cfd05[ch])? hg_cfd05[ch] : 999);
            v03.push_back(ok(hg_cfd03[ch])? hg_cfd03[ch] : 999);
            v10.push_back(ok(hg_cfd10[ch])? hg_cfd10[ch] : 999);
            v30.push_back(ok(hg_cfd30[ch])? hg_cfd30[ch] : 999);
            v50.push_back(ok(hg_cfd50[ch])? hg_cfd50[ch] : 999);
            vled.push_back(ok(hg_led[ch])? hg_led[ch] : 999);
            vamp.push_back(hg_peak[ch]);
            vx.push_back(x_trk); vy.push_back(y_trk); vmcp.push_back(mcp_peak);
        }
        double med = median(v20);

        // ---- (3a) method comparison: core sigma (truncated +/-1.5 ns) + valid% ----
        printf("\n (3a) per-method core resolution (truncated RMS within +/-1.5 ns of median)\n");
        printf("   %-9s  coreSig(ps)  valid%%\n","method");
        auto report=[&](const char* l, std::vector<float>& v){
            double n=0; double s=trms(v,median(v),1.5,&n);
            int valid=0; for(float x:v) if(x<900) ++valid;
            printf("   %-9s  %8.1f     %5.1f\n", l, s*1000.0, 100.0*valid/v.size());
        };
        report("CFD-3%",v03); report("CFD-5%",v05); report("CFD-10%",v10);
        report("CFD-20%",v20); report("CFD-30%",v30); report("CFD-50%",v50);
        report("LED-20mV",vled);

        // ---- (1)/(2)/(3b): core sigma & shoulder fraction vs amplitude bin ----
        // shoulder = |t20 - med| in (0.4, 3.0] ns  (real shoulder, excludes garbage tail)
        const int NB=6; double edg[NB+1]={20,40,60,100,200,400,5000};
        printf("\n (1)/(2)/(3b) CFD-20%% vs amplitude bin\n");
        printf("   amp(mV)      N    coreSig(ps)  shoulder%%   meanX   meanY\n");
        for(int b=0;b<NB;++b){
            std::vector<float> sub; double mx=0,my=0; int nx=0; int nsh=0;
            for(size_t i=0;i<v20.size();++i){
                if(vamp[i]>edg[b] && vamp[i]<=edg[b+1]){
                    sub.push_back(v20[i]);
                    mx+=vx[i]; my+=vy[i]; ++nx;
                    if(std::fabs(v20[i]-med)>0.4 && std::fabs(v20[i]-med)<=3.0) ++nsh;
                }
            }
            if(sub.empty()){ printf("   %4.0f-%-4.0f    %5d   --\n",edg[b],edg[b+1],0); continue; }
            double sig=trms(sub,median(sub),1.5);
            printf("   %4.0f-%-4.0f   %5zu   %8.1f    %6.1f   %6.2f  %6.2f\n",
                   edg[b],edg[b+1],sub.size(),sig*1000.0,100.0*nsh/sub.size(),mx/nx,my/nx);
        }

        // ---- (3c) shoulder events: are they low-amplitude? off-center? -------
        std::vector<size_t> mainI, shI;
        for(size_t i=0;i<v20.size();++i){
            double d=std::fabs(v20[i]-med);
            if(d<=0.4) mainI.push_back(i);
            else if(d<=3.0) shI.push_back(i);   // shoulder, not garbage tail
        }
        auto mean=[&](std::vector<size_t>&id,std::vector<float>&a){double m=0;for(size_t i:id)m+=a[i];return id.empty()?0:m/id.size();};
        auto radius=[&](std::vector<size_t>&id){double m=0;for(size_t i:id)m+=std::sqrt(vx[i]*vx[i]+vy[i]*vy[i]);return id.empty()?0:m/id.size();};
        printf("\n (3c) main(|dt|<0.4ns) vs shoulder(0.4-3ns) populations\n");
        printf("   pop       N     meanAmp(mV)  meanMCP(mV)  mean|r_trk|(mm)\n");
        printf("   main   %5zu     %8.1f     %8.1f      %8.2f\n",
               mainI.size(),mean(mainI,vamp),mean(mainI,vmcp),radius(mainI));
        printf("   shoulder %3zu     %8.1f     %8.1f      %8.2f\n",
               shI.size(),mean(shI,vamp),mean(shI,vmcp),radius(shI));
        printf("   shoulder fraction = %.1f%% of events\n",100.0*shI.size()/v20.size());
    }
    printf("\n[elbowInvestigation] done.\n");
}

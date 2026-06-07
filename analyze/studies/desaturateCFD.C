// desaturateCFD.C — DEFINITIVE test of the LG-referenced de-saturation idea.
// For saturated HG pulses, standard CFD-5% takes 5% of the CLIPPED peak (a fixed
// voltage -> residual time-walk). Here we instead reconstruct the TRUE HG peak
// from the (linear, unsaturated) LG via the per-channel dual-gain slope b_i, and
// place the CFD threshold at 5% of that TRUE peak on the recorded rising edge.
// Compare (DW-UP)/2 sigma_t (best-bin) head-to-head: standard vs de-saturated,
// in identical events, on the raw DSB1 waveforms (150 GeV = most saturated).
//   root -l 'Analysis/desaturateCFD.C+'
#include "WaveformUtils.h"
#include "ChannelConfig.h"
#include "TFile.h"
#include "TChain.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include <cstdio>
#include <vector>
#include <algorithm>
#include <cmath>

static inline int ampOff(int s){return s*1024;}
static inline int timOff(int s){int drs=s/18,grp=(s/9)%2;return (drs*2+grp)*1024;}

// rising-edge crossing at an ABSOLUTE threshold V (mV), same conventions as ExtractPulse
static float crossAbs(const float* time,const float* amp,float V){
    double ped=0;for(int i=3;i<53;++i)ped+=amp[i];ped/=50.0;
    float vmin=amp[3];int imin=3;for(int i=4;i<1024;++i)if(amp[i]<vmin){vmin=amp[i];imin=i;}
    for(int i=3;i<imin;++i){float si=ped-amp[i],sip=ped-amp[i+1];
        if(si<V&&sip>=V){float sl=(sip-si)/(time[i+1]-time[i]);return sl>0?time[i]+(V-si)/sl:kNoTime;}}
    return kNoTime;
}
static double coreS(std::vector<float> v){if(v.size()<200)return -1;std::sort(v.begin(),v.end());
    double mu=v[v.size()/2],s=0.3;for(int it=0;it<5;++it){double a=0,a2=0;long n=0;
    for(float x:v)if(std::fabs(x-mu)<2.5*s){a+=x;a2+=x*x;++n;}if(n<40)break;mu=a/n;double var=a2/n-mu*mu;s=var>0?std::sqrt(var):s;}return s;}

void desaturateCFD(){
    int HG[8]={1,2,3,0,5,4,6,9}, LG[8]={19,20,21,18,23,22,24,25};
    bool useM2[8]={0,0,0,0,0,0,0,1}; int M1=7,M2=16;
    int WCr=28,WCl=29,WCd=30,WCu=32; double WSc=7.0/36.0;

    // ---- Pass 1: calibrate dual-gain slope b_i on UNSATURATED real-signal events ----
    double sHL[8]={0},sLL[8]={0};
    for(const char* f : {"Data/RUN1211_25_GeV.root","Data/RUN1148_50_GeV.root","Data/RUN1112_75_GeV.root"}){
        TFile* fp=TFile::Open(f);if(!fp||fp->IsZombie())continue;TTree* t=(TTree*)fp->Get("pulse");if(!t){fp->Close();continue;}
        TTreeReader r(t);TTreeReaderArray<float> A(r,"amplitude"),T(r,"timevalue");
        while(r.Next()){const float* a=&A[0];const float* tt=&T[0];
            for(int i=0;i<8;++i){Pulse hg=ExtractPulse(tt+timOff(HG[i]),a+ampOff(HG[i]),0.05f,5.f);
                Pulse lg=ExtractPulse(tt+timOff(LG[i]),a+ampOff(LG[i]),0.20f,5.f);
                if(hg.peak>60&&hg.peak<700&&lg.peak>25){sHL[i]+=hg.peak*lg.peak;sLL[i]+=lg.peak*lg.peak;}}}
        fp->Close();}
    double b[8];printf("\n=== dual-gain slope b_i = HG_true/LG (unsaturated fit) ===\n");
    for(int i=0;i<8;++i){b[i]=sLL[i]>0?sHL[i]/sLL[i]:0;printf("  cap%d  b=%.2f\n",i,b[i]);}

    // ---- Pass 2: test on 150 GeV (most saturated) ----
    TChain ch("pulse");
    for(const char* f:{"Data/RUN1258_150_GeV.root","Data/RUN1259_150_GeV.root","Data/RUN1260_150_GeV.root","Data/RUN1261_150_GeV.root"}) ch.Add(f);
    // beam center
    double xs=0,ys=0;long nw=0;
    {TTreeReader r(&ch);TTreeReaderArray<float> A(r,"amplitude"),T(r,"timevalue");long c=0;
     while(r.Next()&&c<40000){++c;const float* a=&A[0];const float* tt=&T[0];
        Pulse R=ExtractPulse(tt+timOff(WCr),a+ampOff(WCr),0.5f,8.f),L=ExtractPulse(tt+timOff(WCl),a+ampOff(WCl),0.5f,8.f);
        Pulse D=ExtractPulse(tt+timOff(WCd),a+ampOff(WCd),0.5f,8.f),U=ExtractPulse(tt+timOff(WCu),a+ampOff(WCu),0.5f,8.f);
        if(R.valid&&L.valid&&D.valid&&U.valid){xs+=WSc*(R.peakTime-L.peakTime);ys+=WSc*(D.peakTime-U.peakTime);++nw;}}}
    double xc=nw?xs/nw:0,yc=nw?ys/nw:0;
    printf("\n150 GeV beam center (%.1f,%.1f) from %ld events\n",xc,yc,nw);

    // store corner times at CFD 5% and 20%, standard (frac of CLIPPED peak) and
    // de-saturated (frac of LG-reconstructed TRUE peak)
    std::vector<float> slg, t5s,t5d, t20s,t20d;
    {TTreeReader r(&ch);TTreeReaderArray<float> A(r,"amplitude"),T(r,"timevalue");
     while(r.Next()){const float* a=&A[0];const float* tt=&T[0];
        Pulse R=ExtractPulse(tt+timOff(WCr),a+ampOff(WCr),0.5f,8.f),L=ExtractPulse(tt+timOff(WCl),a+ampOff(WCl),0.5f,8.f);
        Pulse D=ExtractPulse(tt+timOff(WCd),a+ampOff(WCd),0.5f,8.f),U=ExtractPulse(tt+timOff(WCu),a+ampOff(WCu),0.5f,8.f);
        if(!(R.valid&&L.valid&&D.valid&&U.valid))continue;
        double X=WSc*(R.peakTime-L.peakTime),Y=WSc*(D.peakTime-U.peakTime);
        if((X-xc)*(X-xc)+(Y-yc)*(Y-yc)>=9.0)continue;
        Pulse m1=ExtractPulse(tt+timOff(M1),a+ampOff(M1),0.20f,30.f),m2=ExtractPulse(tt+timOff(M2),a+ampOff(M2),0.20f,30.f);
        if(!m1.valid||m1.peak>750)continue;
        double dw5s=0,up5s=0,dw5d=0,up5d=0,dw20s=0,up20s=0,dw20d=0,up20d=0;int dn=0,un=0;bool ok=true;double SLG=0;
        for(int i=0;i<8;++i){
            const float* ht=tt+timOff(HG[i]); const float* ha=a+ampOff(HG[i]);
            Pulse hg=ExtractPulse(ht,ha,0.05f,5.f);
            Pulse lg=ExtractPulse(tt+timOff(LG[i]),a+ampOff(LG[i]),0.20f,5.f);
            if(!hg.valid||!lg.valid){ok=false;break;}
            float mref=useM2[i]?m2.crossingTime:m1.crossingTime; if(mref<-1e5){ok=false;break;}
            float Vtrue=b[i]*lg.peak;                          // reconstructed true HG peak
            float c5s =crossAbs(ht,ha,0.05f*hg.peak), c5d =crossAbs(ht,ha,0.05f*Vtrue);
            float c20s=crossAbs(ht,ha,0.20f*hg.peak), c20d=crossAbs(ht,ha,0.20f*Vtrue);
            if(c5s<-1e5||c5d<-1e5||c20s<-1e5||c20d<-1e5){ok=false;break;}
            if(i<4){dw5s+=c5s-mref;dw5d+=c5d-mref;dw20s+=c20s-mref;dw20d+=c20d-mref;++dn;}
            else   {up5s+=c5s-mref;up5d+=c5d-mref;up20s+=c20s-mref;up20d+=c20d-mref;++un;}
            SLG+=lg.peak;
        }
        if(!ok||dn<4||un<4)continue;
        slg.push_back((float)SLG);
        t5s.push_back(0.5f*(float)(dw5s/dn-up5s/un));   t5d.push_back(0.5f*(float)(dw5d/dn-up5d/un));
        t20s.push_back(0.5f*(float)(dw20s/dn-up20s/un));t20d.push_back(0.5f*(float)(dw20d/dn-up20d/un));
    }}
    long nt=slg.size();printf("150 GeV fiducial events with all 8 caps valid: %ld\n",nt);
    if(nt<2000){printf("low stats\n");return;}
    // best-bin helper over the same SumLG bins
    std::vector<float> s2=slg;std::sort(s2.begin(),s2.end());double md=s2[s2.size()/2],a=0,a2=0;long n=0;
    for(float v:slg)if(std::fabs(v-md)<0.5*md){a+=v;a2+=v*v;++n;}double mE=a/n,sE=std::sqrt(a2/n-mE*mE),lo=mE-2*sE,bw=4*sE/9;
    auto bestbin=[&](std::vector<float>&tv){double best=1e9;for(int bn=0;bn<9;++bn){double blo=lo+bn*bw,bhi=blo+bw;std::vector<float> v;
        for(long i=0;i<nt;++i)if(slg[i]>=blo&&slg[i]<bhi)v.push_back(tv[i]);if((long)v.size()<500)continue;double s=coreS(v)*1000;if(s>0&&s<best)best=s;}return best;};
    printf("\n=== 150 GeV (DW-UP)/2 sigma_t [ps]: standard (frac of clip) vs de-saturated (frac of LG-true) ===\n");
    printf("  %-14s %10s %10s\n","","standard","de-sat");
    printf("  CFD-5%%  best-bin %10.1f %10.1f   (all-fid %.1f / %.1f)\n",bestbin(t5s),bestbin(t5d),coreS(t5s)*1000,coreS(t5d)*1000);
    printf("  CFD-20%% best-bin %10.1f %10.1f   (all-fid %.1f / %.1f)\n",bestbin(t20s),bestbin(t20d),coreS(t20s)*1000,coreS(t20d)*1000);
}

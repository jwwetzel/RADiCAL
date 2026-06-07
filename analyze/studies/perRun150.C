// perRun150.C — per-run MCP/SiPM leading-edge timing for the four 150 GeV runs.
// Are they identical (as expected) or does the absolute MCP time differ run-to-run?
//   ROOT_INCLUDE_PATH=lib/waveform:lib/io:lib/physics:lib/viz root -l -b -q analyze/studies/perRun150.C+
#include "ChannelConfig.h"
#include "WaveformUtils.h"
#include "DataPaths.h"
#include "TChain.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>

static void stat(std::vector<float> v, double&mean,double&rms,double&md){
    if(v.empty()){mean=rms=md=0;return;} double s=0,ss=0; for(float x:v){s+=x;ss+=x*x;}
    mean=s/v.size(); rms=std::sqrt(std::max(0.0,ss/v.size()-mean*mean));
    std::sort(v.begin(),v.end()); md=v[v.size()/2];
}
void perRun150(){
    const char* runs[4]={"RUN1258_150_GeV.root","RUN1259_150_GeV.root","RUN1260_150_GeV.root","RUN1261_150_GeV.root"};
    for(int r=0;r<4;++r){
        TChain ch("pulse"); ch.Add(radRaw(runs[r])); TTreeReader rd(&ch);
        TTreeReaderArray<float> amp(rd,"amplitude"), tim(rd,"timevalue");
        std::vector<float> mt, ht, diff; long bad=0,tot=0;
        while(rd.Next()){ const float*A=&amp[0]; const float*T=&tim[0]; ++tot;
            PulseMulti mc=ExtractPulseMulti(T+kMCP1_t,A+kMCP1,20.f,30.f,2000.f);
            if(mc.peak<kMCP1_minPeak||mc.peak>kMCP1_maxPeak){++bad;continue;}
            PulseMulti h=ExtractPulseMulti(T+kCap[0].hg_t,A+kCap[0].hg,20.f,5.f,820.f);
            if(mc.ledTime>-1e4) mt.push_back(mc.ledTime);
            if(h.ledTime>-1e4&&h.peak>=20){ ht.push_back(h.ledTime); if(mc.ledTime>-1e4) diff.push_back(h.ledTime-mc.ledTime); }
        }
        double mm,mr,mmd,hm,hr,hmd,dm,dr,dmd;
        stat(mt,mm,mr,mmd); stat(ht,hm,hr,hmd); stat(diff,dm,dr,dmd);
        printf("%s  N=%ld badMCP=%ld\n", runs[r], tot, bad);
        printf("     mcp_led  mean=%8.2f rms=%8.2f med=%8.2f ns\n", mm,mr,mmd);
        printf("     hg_led   mean=%8.2f rms=%8.2f med=%8.2f ns\n", hm,hr,hmd);
        printf("     hg-mcp   mean=%8.2f rms=%8.2f med=%8.2f ns\n", dm,dr,dmd);
    }
}

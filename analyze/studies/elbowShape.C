// elbowShape.C — dump the signed shoulder shape for SE-D (idx2) at 50 GeV.
//   Tests "+0.2 ns DRS4 satellite" (would be a spike at med+0.2) vs amplitude
//   walk (broad one-sided skew, amplitude-correlated).
#include <cstdio>
#include <vector>
#include <algorithm>
#include <cmath>
#include "TFile.h"
#include "TTree.h"

void elbowShape(){
    TFile f("output/50GeV/ntuple.root");
    TTree* t=(TTree*)f.Get("rad");
    Float_t hg_peak[8],hg_cfd[8],hg_cfd05[8]; Bool_t in_fid;
    t->SetBranchAddress("hg_peak",hg_peak);
    t->SetBranchAddress("hg_cfd",hg_cfd);
    t->SetBranchAddress("hg_cfd05",hg_cfd05);
    t->SetBranchAddress("in_fiducial",&in_fid);
    const int ch=2;
    std::vector<float> dt, amp;
    std::vector<float> tmp;
    for(Long64_t e=0;e<t->GetEntries();++e){ t->GetEntry(e);
        if(!in_fid||hg_cfd[ch]<-60||hg_cfd[ch]>60||hg_peak[ch]<=20) continue;
        tmp.push_back(hg_cfd[ch]); }
    std::sort(tmp.begin(),tmp.end()); double med=tmp[tmp.size()/2];
    for(Long64_t e=0;e<t->GetEntries();++e){ t->GetEntry(e);
        if(!in_fid||hg_cfd[ch]<-60||hg_cfd[ch]>60||hg_peak[ch]<=20) continue;
        dt.push_back(hg_cfd[ch]-med); amp.push_back(hg_peak[ch]); }

    printf("SE-D 50 GeV  median(CFD20)=%.3f ns  N=%zu\n",med,dt.size());
    printf("\n signed offset histogram (t20 - median), bin=0.1 ns:\n");
    printf("  offset(ns)  count   meanAmp(mV)\n");
    for(double lo=-1.0; lo<2.0; lo+=0.1){
        int c=0; double a=0;
        for(size_t i=0;i<dt.size();++i) if(dt[i]>=lo&&dt[i]<lo+0.1){++c; a+=amp[i];}
        char bar[81]; int nb=std::min(70,c/3); for(int k=0;k<nb;++k)bar[k]='#'; bar[nb]=0;
        printf("  %+5.1f      %4d   %7.1f  %s\n", lo, c, c?a/c:0.0, bar);
    }
    // counts in key windows
    auto inwin=[&](double a,double b){int c=0;for(float x:dt)if(x>=a&&x<b)++c;return c;};
    printf("\n core |dt|<0.2 : %d   |0.2-0.4| : %d   late 0.4-1.5 : %d   early <-0.4 : %d\n",
        inwin(-0.2,0.2), inwin(0.2,0.4)+inwin(-0.4,-0.2), inwin(0.4,1.5), inwin(-1.5,-0.4));
    printf(" satellite test: bump localized at +0.2 ns? compare 0.1-0.3 (%d) vs spread 0.4-1.5 (%d)\n",
        inwin(0.1,0.3), inwin(0.4,1.5));
}

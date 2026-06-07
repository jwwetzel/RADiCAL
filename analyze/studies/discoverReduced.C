// discoverReduced.C — derive a config's slot->role map from a REDUCED ntuple
// (the config-agnostic output of reduceRaw.C).  Classifies all 36 DRS slots by
// occupancy, mean amplitude, and pulse width (charge/peak), so the live timing
// vs energy capillaries of any capillary config (LuAG, mixed, ...) can be read
// off without prior knowledge.
//   root -l 'Analysis/discoverReduced.C+("reduced/LUAG/150GeV.root")'
#include "TFile.h"
#include "TTree.h"
#include <cstdio>
#include <vector>

void discoverReduced(const char* file, long maxEv=300000)
{
    TFile f(file);
    if(f.IsZombie()){ printf("cannot open %s\n",file); return; }
    TTree* t=(TTree*)f.Get("rad");
    const int N=36;
    Float_t speak[36], scharge[36];
    t->SetBranchAddress("s_peak",  speak);
    t->SetBranchAddress("s_charge",scharge);

    std::vector<long>   nv(N,0);
    std::vector<double> sp(N,0), sw(N,0);   // sum peak, sum width(=charge/peak)
    long ne=std::min((long)t->GetEntries(),maxEv);
    for(long i=0;i<ne;++i){ t->GetEntry(i);
        for(int s=0;s<N;++s){ if(speak[s]>15.f){ ++nv[s]; sp[s]+=speak[s];
            if(speak[s]>0) sw[s]+=scharge[s]/speak[s]; } } }

    printf("\n%s   (%ld events)\n",file,ne);
    printf("%-9s %5s %8s %8s  %-16s\n","slot(D,G,c)","occ%","peak","width","role");
    printf("%s\n",std::string(56,'-').c_str());
    for(int s=0;s<N;++s){
        int drs=s/18, grp=(s/9)%2, ch=s%9;
        double occ=100.*nv[s]/ne, pk=nv[s]?sp[s]/nv[s]:0, wd=nv[s]?sw[s]/nv[s]:0;
        const char* role;
        if(ch==8)                          role="TR/sync (DAQ)";
        else if(occ<5)                     role="dead/unused";
        else if(drs==1&&grp==1)            role="WC / scint";
        else if(drs==1&&grp==0)            role="LG energy (wide)";
        else if(drs==0&&pk<260&&wd<8)      role="MCP / fast ref";
        else if(drs==0&&pk>=300)           role="HG timing";
        else if(drs==0&&grp==1)            role="PbGlass / aux";
        else                               role="(weak/ambiguous)";
        printf("D%dG%dc%-3d %5.0f %8.0f %8.1f  %-16s\n",drs,grp,ch,occ,pk,wd,role);
    }
}

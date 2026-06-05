// etypeEnergy.C — energy resolution from the ENERGY (E-type) capillary.
//   The 2023 run has NO 4xE-type configuration (max 1 E-type per build), so a true
//   "energy-capillary-only" resolution is not available. This shows the best the data
//   allows: sigma_E/E vs beam energy from
//     (a) the single E-type cap (TENERGY NW, full-length WLS, the energy element),
//     (b) a single T-type cap (avg of the 3 DSB1 timing caps, localized at shower max),
//     (c) the full 8-low-gain sum (the standard shower-max energy estimate),
//   so one can compare the per-capillary energy performance of the two cap types.
//   root -l 'Analysis/etypeEnergy.C+'
#include "ChannelConfig.h"
#include "PlotUtils.h"
#include "DataPaths.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include <cstdio>
#include <vector>
#include <algorithm>
#include <cmath>

// sigma/mean (%) of a distribution, median-anchored truncated core
static double sErel(std::vector<float> v){
    if(v.size()<150) return -1; std::sort(v.begin(),v.end());
    double md=v[v.size()/2], a=0,a2=0; long n=0;
    for(float x:v) if(std::fabs(x-md)<0.45*md){a+=x;a2+=x*x;++n;}
    if(n<60) return -1; double mu=a/n, var=a2/n-mu*mu; return var>0? 100.0*std::sqrt(var)/mu : -1;
}
void etypeEnergy(){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    int lg[8]; for(int i=0;i<8;++i) lg[i]=kCap[i].lg/1024;
    double Es[5]={50,75,100,125,150};
    std::vector<double> ex, eE, eT, eS;        // E-type(1), T-type(1 avg), 8-sum
    printf("\n=== TENERGY energy resolution sigma_E/E [%%] ===\n");
    printf("%5s %14s %16s %12s\n","E","E-type(1 cap)","T-type(1 cap avg)","8-LG sum");
    for(int e=0;e<5;++e){
        TFile* fp=TFile::Open(radReduced("TENERGY",Es[e])); if(!fp||fp->IsZombie()) continue;
        TTree* t=(TTree*)fp->Get("rad"); if(!t){fp->Close();continue;}
        Bool_t wc; Float_t x,y,m1p,sp[36];
        t->SetBranchAddress("wc_ok",&wc); t->SetBranchAddress("x_trk",&x); t->SetBranchAddress("y_trk",&y);
        t->SetBranchAddress("mcp1_peak",&m1p); t->SetBranchAddress("s_peak",sp);
        long N=t->GetEntries(); double xs=0,ys=0; long nw=0;
        for(long i=0;i<N&&nw<40000;++i){ t->GetEntry(i); if(wc&&x>-100&&x<100){xs+=x;ys+=y;++nw;} }
        double xc=nw?xs/nw:0, yc=nw?ys/nw:0;
        std::vector<float> ve, vt1,vt2,vt3, vs;
        for(long i=0;i<N;++i){ t->GetEntry(i);
            if(!wc||m1p<200||m1p>750) continue; double dx=x-xc,dy=y-yc; if(dx*dx+dy*dy>=9.0) continue;
            ve.push_back(sp[lg[0]]+sp[lg[4]]);                       // E-type NW (down+up)
            vt1.push_back(sp[lg[1]]+sp[lg[5]]);                      // T-type NE
            vt2.push_back(sp[lg[2]]+sp[lg[6]]);                      // T-type SE
            vt3.push_back(sp[lg[3]]+sp[lg[7]]);                      // T-type SW
            double s=0; for(int k=0;k<8;++k) s+=sp[lg[k]]; vs.push_back((float)s);  // 8-sum
        }
        double rE=sErel(ve), rT=(sErel(vt1)+sErel(vt2)+sErel(vt3))/3.0, rS=sErel(vs);
        printf("%5.0f %14.1f %16.1f %12.1f\n",Es[e],rE,rT,rS);
        ex.push_back(Es[e]); eE.push_back(rE); eT.push_back(rT); eS.push_back(rS);
        fp->Close();
    }
    if(ex.size()<3){ printf("low stats\n"); return; }
    TCanvas* c=new TCanvas("c_ee","",900,640); c->SetGridx(); c->SetGridy();
    TH1F* fr=c->DrawFrame(40,10,160,19);   // zoom so the gentle a/sqrtE falloff is visible
    fr->SetTitle("Per-capillary energy resolution: E-type vs T-type (TENERGY);Beam energy (GeV);#sigma_{E}/E  [%]");
    printf("\n=== fits  sigma_E/E = a/sqrtE (+) c  ===\n");
    auto mk=[&](std::vector<double>&y,int col,int ms,const char* nm)->TGraph*{
        TGraph* g=new TGraph(ex.size(),&ex[0],&y[0]); g->SetMarkerStyle(ms); g->SetMarkerColor(col); g->SetLineColor(col); g->SetMarkerSize(1.6);
        TF1* f=new TF1(Form("f_%s",nm),"sqrt([0]*[0]/x+[1]*[1])",45,158); f->SetParameters(40,13); f->SetLineColor(col); f->SetLineStyle(2); f->SetLineWidth(2);
        g->Fit(f,"RQN"); g->Draw("P SAME"); f->Draw("SAME");
        printf("  %-7s  a = %4.0f %%#sqrt{GeV}   c = %4.1f %% (constant)\n",nm,fabs(f->GetParameter(0)),fabs(f->GetParameter(1)));
        return g; };
    TGraph* gE=mk(eE,kAzure+1,21,"E-type"); TGraph* gT=mk(eT,kRRed,20,"T-type"); TGraph* gS=mk(eS,kGreen+3,33,"8sum");
    TLegend* lg2=new TLegend(0.13,0.135,0.55,0.33); lg2->SetTextSize(0.027); lg2->SetFillColorAlpha(0,0);
    lg2->AddEntry(gE,"E-type cap (full-length, energy)","pl");
    lg2->AddEntry(gT,"T-type cap (shower-max, single)","pl");
    lg2->AddEntry(gS,"all 8 low-gain summed","pl");
    lg2->Draw();
    TLatex tl; tl.SetNDC(); tl.SetTextSize(0.027);
    tl.DrawLatex(0.56,0.305,"Nearly flat = #bf{constant-term dominated}:");
    tl.DrawLatex(0.56,0.255,"the a/#sqrt{E} stochastic term is small over");
    tl.DrawLatex(0.56,0.205,"50-150 GeV; the large sampling / position");
    tl.DrawLatex(0.56,0.155,"constant term c sets #sigma_{E}/E (dashed fit).");
    c->Print("Analysis/capillary_figs/etype_energy_resolution.png");
    printf("wrote Analysis/capillary_figs/etype_energy_resolution.png\n");
}

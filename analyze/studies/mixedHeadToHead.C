// mixedHeadToHead.C — the clean DSB1-vs-LuAG crystal comparison, IN-EVENT.
// ----------------------------------------------------------------------------
// The MIXED module ("2xDSB1, 2xLuAG") has both crystals reading the SAME showers.
// Timing slots 0-6 all live in DRS0 Group 0 (one time axis), so the pairwise
// difference of any two of their CFD-5% times cancels EVERYTHING common to the
// event -- arrival time, the group/run timebase offset, the MCP reference --
// leaving exactly sqrt(s_i^2 + s_j^2), the quadrature of the two capillaries'
// intrinsic jitters.  Restricting to SAME-LAYER pairs (Down: 0,1,2,3 ; Up:
// 4,5,6) also removes the longitudinal shower-time spread.  Solving the pairwise
// system per layer gives each capillary's intrinsic sigma_t with no proxy, no
// saturation confound, no run-offset -- the definitive crystal comparison.
//
// Crystal id is read from the data (bright caps = DSB1, dim = LuAG).
//   root -l 'Analysis/mixedHeadToHead.C+'
#include "PlotUtils.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TLegend.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>

static double coreSig(std::vector<float> v){
    if(v.size()<300)return -1; std::sort(v.begin(),v.end()); double mu=v[v.size()/2],s=0.2;
    for(int it=0;it<5;++it){double a=0,a2=0;long n=0;for(float x:v)if(std::fabs(x-mu)<2.5*s){a+=x;a2+=x*x;++n;}
        if(n<50)break;mu=a/n;double var=a2/n-mu*mu;s=var>0?std::sqrt(var):s;} return s;}

void mixedHeadToHead(){
    ApplyRADiCALStyle();
    double Es[5]={50,75,100,125,150};
    int down[4]={0,1,2,3}, up[3]={4,5,6};   // DRS0 G0 slot indices, by layer
    std::vector<double> ED, sD, sL;          // energy, DSB1 cap sigma, LuAG cap sigma
    printf("\n=== MIXED in-event DSB1-vs-LuAG (per-capillary intrinsic sigma_t) ===\n");
    for(int e=0;e<5;++e){
        TFile* fp=TFile::Open(Form("reduced/MIXED/%.0fGeV.root",Es[e])); if(!fp||fp->IsZombie())continue;
        TTree* t=(TTree*)fp->Get("rad"); if(!t){fp->Close();continue;}
        Bool_t wc; Float_t x,y,sp[36],sc[36];
        t->SetBranchAddress("wc_ok",&wc);t->SetBranchAddress("x_trk",&x);t->SetBranchAddress("y_trk",&y);
        t->SetBranchAddress("s_peak",sp);t->SetBranchAddress("s_cfd05",sc);
        long N=t->GetEntries();
        double xs=0,ys=0;long nw=0; for(long i=0;i<N&&nw<40000;++i){t->GetEntry(i);if(wc&&x>-100&&x<100){xs+=x;ys+=y;++nw;}}
        double xc=nw?xs/nw:0,yc=nw?ys/nw:0;
        // collect pairwise CFD-5% differences (same-layer) + per-slot amplitude
        std::vector<float> dpair[4][4]; std::vector<float> upair[7][7];   // sparse use
        double amp[7]={0}; long namp[7]={0};
        for(long i=0;i<N;++i){t->GetEntry(i); if(!wc)continue; double dx=x-xc,dy=y-yc; if(dx*dx+dy*dy>=9.0)continue;
            for(int a=0;a<7;++a){int s=(a<4)?down[a]:up[a-4]; if(sp[s]>20){amp[a]+=sp[s];++namp[a];}}
            for(int a=0;a<4;++a)for(int b=a+1;b<4;++b){int si=down[a],sj=down[b];
                if(sp[si]>20&&sp[sj]>20&&sc[si]>-1e5&&sc[sj]>-1e5) dpair[a][b].push_back(sc[si]-sc[sj]);}
            for(int a=0;a<3;++a)for(int b=a+1;b<3;++b){int si=up[a],sj=up[b];
                if(sp[si]>20&&sp[sj]>20&&sc[si]>-1e5&&sc[sj]>-1e5) upair[a][b].push_back(sc[si]-sc[sj]);}
        }
        for(int a=0;a<7;++a) amp[a]= namp[a]? amp[a]/namp[a] : 0;
        double mx=0; for(int a=0;a<7;++a) mx=std::max(mx,amp[a]); double thr=0.65*mx;
        // classify slots 0..6 (index a): DSB1 if amp>thr
        bool isD[7]; for(int a=0;a<7;++a) isD[a]=amp[a]>thr;
        // solve per-layer x_i = s_i^2 (ns^2)
        double x_d[4], x_u[3];
        { double sij[4][4]={{0}}; for(int a=0;a<4;++a)for(int b=a+1;b<4;++b){double s=coreSig(dpair[a][b]); double v=s*s; sij[a][b]=sij[b][a]=v;}
          double R[4],S=0; for(int a=0;a<4;++a){R[a]=0;for(int b=0;b<4;++b)if(b!=a)R[a]+=sij[a][b];S+=R[a];} S/=6.0;
          for(int a=0;a<4;++a) x_d[a]=(R[a]-S)/2.0; }
        { double sij[3][3]={{0}}; for(int a=0;a<3;++a)for(int b=a+1;b<3;++b){double s=coreSig(upair[a][b]); double v=s*s; sij[a][b]=sij[b][a]=v;}
          double R[3],S=0; for(int a=0;a<3;++a){R[a]=0;for(int b=0;b<3;++b)if(b!=a)R[a]+=sij[a][b];S+=R[a];} S/=4.0;
          for(int a=0;a<3;++a) x_u[a]=R[a]-S; }
        // group DSB1 / LuAG (sigma in ps), averaging x_i (variance) then sqrt
        double sumD=0,sumL=0; int nD=0,nL=0;
        for(int a=0;a<4;++a){ double xi=x_d[a]; if(xi<=0)continue; if(isD[a]){sumD+=xi;++nD;} else {sumL+=xi;++nL;} }
        for(int a=0;a<3;++a){ double xi=x_u[a]; if(xi<=0)continue; if(isD[a+4]){sumD+=xi;++nD;} else {sumL+=xi;++nL;} }
        double sigD=nD?std::sqrt(sumD/nD)*1000.0:-1, sigL=nL?std::sqrt(sumL/nL)*1000.0:-1;
        printf(" E=%3.0f GeV  DSB1 caps: %.0f ps (n=%d)   LuAG caps: %.0f ps (n=%d)   thr=%.0f mV\n",
               Es[e],sigD,nD,sigL,nL,thr);
        if(sigD>0&&sigL>0){ED.push_back(Es[e]);sD.push_back(sigD);sL.push_back(sigL);}
        fp->Close();
    }
    // figure
    TCanvas* c=new TCanvas("c_h2h","",900,680);c->SetLeftMargin(0.13);c->SetBottomMargin(0.13);c->SetTopMargin(0.10);
    TH1F* fr=gPad->DrawFrame(0,100,165,320);
    fr->GetXaxis()->SetTitle("beam energy (GeV)");fr->GetYaxis()->SetTitle("per-capillary intrinsic #sigma_{t} (ps)");
    fr->GetXaxis()->SetTitleSize(0.045);fr->GetYaxis()->SetTitleSize(0.045);
    TGraph* gD=new TGraph(ED.size(),ED.data(),sD.data()); gD->SetLineColor(kRData);gD->SetMarkerColor(kRData);gD->SetMarkerStyle(20);gD->SetMarkerSize(1.6);gD->SetLineWidth(3);
    TGraph* gL=new TGraph(ED.size(),ED.data(),sL.data()); gL->SetLineColor(kRGreen+1);gL->SetMarkerColor(kRGreen+1);gL->SetMarkerStyle(21);gL->SetMarkerSize(1.6);gL->SetLineWidth(3);
    gD->Draw("PL");gL->Draw("PL");
    TLegend* L=new TLegend(0.50,0.74,0.93,0.88);L->SetBorderSize(0);L->SetFillStyle(0);L->SetTextSize(0.036);
    L->AddEntry(gD,"DSB1 capillaries","pl");L->AddEntry(gL,"LuAG capillaries","pl");L->Draw();
    {TLatex t;t.SetNDC();t.SetTextSize(0.029);t.SetTextColor(kGray+3);
     t.DrawLatex(0.16,0.30,"Same module, same showers, same MCP, same DRS group.");
     t.DrawLatex(0.16,0.255,"Pairwise-difference method -> offset/MCP/shower-time cancel.");
     t.DrawLatex(0.16,0.21,"The ONLY difference is the crystal.");}
    DrawPageTitle("MIXED in-event head-to-head: DSB1 vs LuAG capillary timing");
    c->Print("/tmp/mixed_h2h.png");
    printf("\nwrote /tmp/mixed_h2h.png\n");
}

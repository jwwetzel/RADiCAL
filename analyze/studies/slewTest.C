// slewTest.C — sigma_t vs light yield (SumLG), all energies pooled, per config.
// Tests "it's light yield, not crystal speed": if timing follows noise/slew with
// slew proportional to light, DSB1 and LuAG fall on the SAME sigma_t-vs-SumLG
// curve, DSB1 simply extending to higher light (it saturates the HG because it
// makes more light).  SumLG is the LOW-gain energy sum -> unsaturated, ~ light.
#include "PlotUtils.h"
#include "ChannelConfig.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TLegend.h"
#include <vector>
#include <algorithm>
#include <cmath>

static double coreSig(std::vector<float> v){
    if(v.size()<300)return -1; std::sort(v.begin(),v.end()); double mu=v[v.size()/2],s=0.2;
    for(int it=0;it<5;++it){double a=0,a2=0;long n=0;for(float x:v)if(std::fabs(x-mu)<2.5*s){a+=x;a2+=x*x;++n;}
        if(n<50)break;mu=a/n;double var=a2/n-mu*mu;s=var>0?std::sqrt(var):s;} return s;}

static void center(TTree*t,Float_t&x,Float_t&y,Bool_t&wc,double&xc,double&yc){
    long N=t->GetEntries();double xs=0,ys=0;long nw=0;
    for(long i=0;i<N&&nw<40000;++i){t->GetEntry(i);if(wc&&x>-100&&x<100){xs+=x;ys+=y;++nw;}}
    xc=nw?xs/nw:0;yc=nw?ys/nw:0;}

// reduced format
static void collectRed(const char* dir,std::vector<float>&SL,std::vector<float>&TC){
    int dw[4],up[4];bool m2[4];int en[8];
    for(int i=0;i<4;++i)dw[i]=kCap[i].hg/1024;
    for(int i=4;i<8;++i){up[i-4]=kCap[i].hg/1024;m2[i-4]=kCap[i].use_mcp2;}
    for(int i=0;i<8;++i)en[i]=kCap[i].lg/1024;
    for(double E:{50.,75.,100.,125.,150.}){
        TFile*fp=TFile::Open(Form("%s/%.0fGeV.root",dir,E));if(!fp||fp->IsZombie())continue;
        TTree*t=(TTree*)fp->Get("rad");if(!t){fp->Close();continue;}
        Bool_t wc;Float_t x,y,m1t,m2t,m1p,sp[36],sc[36];
        t->SetBranchAddress("wc_ok",&wc);t->SetBranchAddress("x_trk",&x);t->SetBranchAddress("y_trk",&y);
        t->SetBranchAddress("mcp1_time",&m1t);t->SetBranchAddress("mcp2_time",&m2t);t->SetBranchAddress("mcp1_peak",&m1p);
        t->SetBranchAddress("s_peak",sp);t->SetBranchAddress("s_cfd05",sc);
        double xc,yc;center(t,x,y,wc,xc,yc);long N=t->GetEntries();
        for(long i=0;i<N;++i){t->GetEntry(i);if(!wc)continue;double dx=x-xc,dy=y-yc;if(dx*dx+dy*dy>=9.0)continue;
            if(m1p<200||m1p>750||m1t<-1e5)continue;
            double ds=0,us=0;int dn=0,un=0;
            for(int k=0;k<4;++k){int s=dw[k];if(sp[s]>20&&sc[s]>-1e5){ds+=sc[s]-m1t;++dn;}}
            for(int k=0;k<4;++k){int s=up[k];double r=m2[k]?m2t:m1t;if(r<-1e5)continue;if(sp[s]>20&&sc[s]>-1e5){us+=sc[s]-r;++un;}}
            if(dn<1||un<1)continue;
            // x = mean HG timing-pulse amplitude (mV) over the 8 timing slots = slew proxy
            double hg=0;int hn=0; for(int k=0;k<4;++k){if(sp[dw[k]]>20){hg+=sp[dw[k]];++hn;}}
            for(int k=0;k<4;++k){if(sp[up[k]]>20){hg+=sp[up[k]];++hn;}}
            if(hn<1)continue;
            SL.push_back((float)(hg/hn));TC.push_back(0.5f*(float)(ds/dn-us/un));}
        fp->Close();}
}
// DSB1 processRun
static void collectDSB1(std::vector<float>&SL,std::vector<float>&TC){
    for(double E:{25.,50.,75.,100.,125.,150.}){
        TFile*fp=TFile::Open(Form("output/%.0fGeV/ntuple.root",E));if(!fp||fp->IsZombie())continue;
        TTree*t=(TTree*)fp->Get("rad");if(!t){fp->Close();continue;}
        Bool_t wc;Float_t x,y,hgp[8],cfd[8];
        t->SetBranchAddress("wc_ok",&wc);t->SetBranchAddress("x_trk",&x);t->SetBranchAddress("y_trk",&y);
        t->SetBranchAddress("hg_peak",hgp);t->SetBranchAddress("hg_cfd05",cfd);
        double xc,yc;center(t,x,y,wc,xc,yc);long N=t->GetEntries();
        for(long i=0;i<N;++i){t->GetEntry(i);if(!wc)continue;double dx=x-xc,dy=y-yc;if(dx*dx+dy*dy>=9.0)continue;
            double ds=0,us=0;int dn=0,un=0;for(int k=0;k<4;++k)if(cfd[k]>-1e5){ds+=cfd[k];++dn;}
            for(int k=4;k<8;++k)if(cfd[k]>-1e5){us+=cfd[k];++un;}if(dn<1||un<1)continue;
            double hg=0;int hn=0;for(int k=0;k<8;++k){if(hgp[k]>20){hg+=hgp[k];++hn;}} if(hn<1)continue;
            SL.push_back((float)(hg/hn));TC.push_back(0.5f*(float)(ds/dn-us/un));}
        fp->Close();}
}
static TGraph* curve(std::vector<float>&SL,std::vector<float>&TC,int col,int mst){
    double bw=40; std::vector<double> bx,by;
    for(double lo=100;lo<820;lo+=bw){ std::vector<float> tv; double sx=0;long n=0;
        for(size_t i=0;i<SL.size();++i)if(SL[i]>=lo&&SL[i]<lo+bw){tv.push_back(TC[i]);sx+=SL[i];++n;}
        if((long)tv.size()<500)continue; double s=coreSig(tv)*1000.0; if(s<=0)continue;
        bx.push_back(sx/n); by.push_back(s);}
    TGraph*g=new TGraph(bx.size(),bx.data(),by.data());
    g->SetLineColor(col);g->SetMarkerColor(col);g->SetMarkerStyle(mst);g->SetMarkerSize(1.5);g->SetLineWidth(3);return g;
}
void slewTest(){
    ApplyRADiCALStyle();
    std::vector<float> dSL,dTC,lSL,lTC,mSL,mTC,tSL,tTC;
    collectDSB1(dSL,dTC); collectRed("reduced/LUAG",lSL,lTC);
    collectRed("reduced/MIXED",mSL,mTC); collectRed("reduced/TENERGY",tSL,tTC);
    TCanvas*c=new TCanvas("c_slew","",900,680); c->SetLeftMargin(0.13);c->SetBottomMargin(0.13);c->SetTopMargin(0.09);
    TH1F*fr=gPad->DrawFrame(0,20,820,90);
    fr->GetXaxis()->SetTitle("mean HG amplitude (mV)  #propto slew");
    fr->GetYaxis()->SetTitle("#sigma_{t} (DW#minusUP)/2 (ps)");
    fr->GetXaxis()->SetTitleSize(0.045);fr->GetYaxis()->SetTitleSize(0.045);
    auto gD=curve(dSL,dTC,kRData,20), gL=curve(lSL,lTC,kRGreen+1,21),
         gM=curve(mSL,mTC,kROrange,22), gT=curve(tSL,tTC,kRRed,23);
    gM->Draw("PL");gT->Draw("PL");gD->Draw("PL");gL->Draw("PL");
    TLegend*L=new TLegend(0.50,0.66,0.92,0.88);L->SetBorderSize(0);L->SetFillStyle(0);L->SetTextSize(0.034);
    L->AddEntry(gD,"DSB1","pl");L->AddEntry(gL,"LuAG","pl");
    L->AddEntry(gM,"2xDSB1+2xLuAG","pl");L->AddEntry(gT,"3xDSB1+1xEnergy","pl");L->Draw();
    { TLatex t;t.SetNDC();t.SetTextSize(0.029);t.SetTextColor(kGray+3);
      t.DrawLatex(0.40,0.55,"At the SAME HG amplitude (same slew),");
      t.DrawLatex(0.40,0.51,"do the crystals agree (light-limited)");
      t.DrawLatex(0.40,0.47,"or does one sit lower (faster pulse)?");
      t.SetTextColor(kRRed+1);t.DrawLatex(0.40,0.41,"DSB1 HG saturates near ~780 mV (right edge).");}
    DrawPageTitle("#sigma_{t} vs HG amplitude: is timing light-limited or crystal-limited?");
    c->Print("/tmp/slew_test.png");
    printf("wrote /tmp/slew_test.png\n");
}

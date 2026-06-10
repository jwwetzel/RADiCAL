// timingAllMethods.C — brightest-1000 (DW-UP)/2 sigma_t(E) for ALL available timing
// methods (cfd03/05/10/20/30/50, led, lgcfd) on one plot per build, via the PRODUCTION
// (fixed) rad::timingBrightestK. a/sqrt(E)+b photostat fits; floor b in the legend.
// Shows the full method landscape and which method is adopted per build.
//   source setup.sh; root -l -b -q 'analyze/studies/timingAllMethods.C+("DSB1")'
// Output: figures/<year>/narrative/timing_allmethods_<BUILD>.png
#include "RadView.h"
#include "RadTiming.h"
#include "DataPaths.h"
#include "FigPaths.h"
#include "PlotUtils.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include <vector>
#include <cstdio>
#include <cmath>
using namespace rad;

// distinct hues + distinct marker shapes (survive in greyscale / colourblind)
static const int MCOL[8]={kGray+1,kBlue+1,kAzure+7,kCyan+2,kGreen+2,kOrange+7,kViolet+1,kRed+1};
static const int MMK [8]={24,20,21,22,23,33,34,29};

void timingAllMethods(const char* build="DSB1"){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    BuildConfig cfg=BuildConfig::Load(radConfig(build).Data());
    const double Es[6]={25,50,75,100,125,150};
    bool avail[8]={};
    std::vector<double> E[8],S[8],Se[8],z[8]; double bfl[8],ber[8],chi[8]; bool mono[8]; bool anyNM=false;
    double minE=1e9,maxE=0;
    for(int e=0;e<6;++e){ int En=(int)Es[e]; TString p=radReduced(build,En); if(gSystem->AccessPathName(p.Data()))continue;
        TFile* f=TFile::Open(p.Data()); if(!f||f->IsZombie()){if(f)f->Close();continue;} TTree* t=(TTree*)f->Get("rad"); if(!t){f->Close();continue;}
        RadView v; v.attach(t,&cfg);
        for(int s=0;s<8;++s){ if(!v.hasSrc(s))continue; avail[s]=true;
            TimingResult r=timingBrightestK(v,En,s,1000);
            if(r.sigma_ps>0){ E[s].push_back(En);S[s].push_back(r.sigma_ps);Se[s].push_back(r.sigma_ps/std::sqrt(2000.0));z[s].push_back(0);
                minE=std::min(minE,(double)En); maxE=std::max(maxE,(double)En);} }
        f->Close();
    }
    if(minE>maxE){minE=25;maxE=150;}
    TCanvas* c=new TCanvas(Form("tam_%s",build),"",1000,760); c->SetLeftMargin(0.12);c->SetRightMargin(0.04);c->SetTopMargin(0.07);c->SetBottomMargin(0.13);c->SetGridy();
    double ymax=0,ymin=1e9; for(int s=0;s<8;++s)for(double y:S[s]){ymax=std::max(ymax,y);ymin=std::min(ymin,y);} if(ymin>ymax){ymin=20;ymax=70;}
    double rng=ymax-ymin; if(rng<5)rng=5;
    TH1F* fr=c->DrawFrame(minE-12,std::max(0.0,ymin-0.10*rng),maxE+12,ymax+0.34*rng);   // tight to data; top headroom
    fr->SetTitle(Form("%s  #minus  brightest-1000 (DW#minusUP)/2  #sigma_{t}(E), all methods;beam energy E (GeV);#sigma_{t} (ps)",build));
    fr->GetYaxis()->SetTitleSize(0.045);fr->GetYaxis()->SetTitleOffset(1.25);fr->GetXaxis()->SetTitleSize(0.045);
    TLegend* lg=new TLegend(0.40,0.605,0.955,0.92); lg->SetBorderSize(0);lg->SetFillStyle(0);lg->SetTextSize(0.030); lg->SetNColumns(2); lg->SetMargin(0.16); lg->SetColumnSeparation(0.04);
    for(int s=0;s<8;++s){ if(!avail[s]||E[s].size()<3)continue;
        TGraphErrors g(E[s].size(),&E[s][0],&S[s][0],&z[s][0],&Se[s][0]);
        TF1 f("f","sqrt([0]*[0]/x+[1]*[1])",minE,maxE); f.SetParameters(300,18); g.Fit(&f,"QN");
        bfl[s]=std::fabs(f.GetParameter(1)); chi[s]=f.GetChisquare()/std::max(1,f.GetNDF()); ber[s]=f.GetParError(1); if(chi[s]>1)ber[s]*=std::sqrt(chi[s]);
        mono[s]=true; for(size_t i=1;i<S[s].size();++i) if(S[s][i]>S[s][i-1]+0.1) mono[s]=false; if(!mono[s])anyNM=true;
        // small horizontal dodge so the 8 co-located markers fan out
        std::vector<double> Ed(E[s].size()); for(size_t i=0;i<E[s].size();++i) Ed[i]=E[s][i]+(s-3.5)*1.1;
        TGraphErrors* gg=new TGraphErrors(Ed.size(),&Ed[0],&S[s][0],&z[s][0],&Se[s][0]);
        gg->SetMarkerStyle(MMK[s]);gg->SetMarkerColor(MCOL[s]);gg->SetLineColor(MCOL[s]);gg->SetMarkerSize(1.35);
        if(chi[s]<10){TF1* ff=new TF1(Form("ff%d%s",s,build),"sqrt([0]*[0]/x+[1]*[1])",minE,maxE);ff->SetParameters(f.GetParameter(0),bfl[s]);ff->SetLineColor(MCOL[s]);ff->SetLineWidth(2);ff->Draw("SAME");}
        gg->Draw("P SAME");
        lg->AddEntry(gg,Form("%s  b=%.0f%s",RadView::srcName(s),bfl[s],mono[s]?"":" *"),"lp");
    }
    lg->Draw();
    if(anyNM){ TLatex nt; nt.SetNDC();nt.SetTextSize(0.025);nt.SetTextColor(kGray+2); nt.DrawLatex(0.40,0.575,"b = floor (ps);  * = non-monotonic, floor unreliable"); }
    c->Print(radFigP(Form("figures/narrative/timing_allmethods_%s.png",build)));
    printf("wrote figures/narrative/timing_allmethods_%s.png\n",build);
}

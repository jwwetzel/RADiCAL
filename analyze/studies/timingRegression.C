// timingRegression.C — validate the sigma-monotonicity fix using the PRODUCTION
// RadTiming.h functions (patched tebSigma robust window + eventDWUP veto + continuous
// TimingFiducialR). Checks: (1) every adopted-method sigma_t(E) is monotonic;
// (2) the published DSB1 cfd05 best-bin headline (27.4 ps @150) + a/sqrt(E)+b floor
// (~17.5 ps) are preserved.
//   source setup.sh; root -l -b -q 'analyze/studies/timingRegression.C+'
#include "RadView.h"
#include "RadTiming.h"
#include "DataPaths.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TSystem.h"
#include <vector>
#include <cstdio>
#include <cmath>
using namespace rad;

static void series(const char* build,int src,bool bestbin,const char* tag){
    BuildConfig cfg=BuildConfig::Load(radConfig(build).Data());
    const double Es[6]={25,50,75,100,125,150};
    std::vector<double> E,S,Se,ze;
    printf("  %-22s :",tag);
    double prev=-1; bool mono=true;
    for(double e:Es){ TString p=radReduced(build,e); if(gSystem->AccessPathName(p.Data())){printf("   --  ");continue;}
        TFile* f=TFile::Open(p.Data()); TTree* t=(TTree*)f->Get("rad"); if(!t){f->Close();printf("   --  ");continue;}
        RadView v; v.attach(t,&cfg);
        TimingResult r = bestbin ? timingBestBin(v,e,src) : timingBrightestK(v,e,src,1000);
        double s=r.sigma_ps; f->Close();
        if(s>0){ printf(" %5.1f",s); if(prev>0 && s>prev+0.05) mono=false; prev=s;
                 E.push_back(e);S.push_back(s);Se.push_back(s/std::sqrt(2000.0));ze.push_back(0);}
        else printf("   --  ");
    }
    double b=-1,berr=-1,chi=-1;
    if(E.size()>=3){ TGraphErrors g(E.size(),&E[0],&S[0],&ze[0],&Se[0]);
        TF1 fit("fit","sqrt([0]*[0]/x+[1]*[1])",20,160); fit.SetParameters(300,18); g.Fit(&fit,"QN");
        b=std::fabs(fit.GetParameter(1)); chi=fit.GetChisquare()/std::max(1,fit.GetNDF()); berr=fit.GetParError(1); if(chi>1)berr*=std::sqrt(chi); }
    printf("   | floor b=%.1f+-%.1f ps  chi2/ndf=%.1f  %s\n",b,berr,chi, mono?"MONOTONIC":"*** NON-MONOTONIC ***");
}

void timingRegression(){
    printf("\n===== PRODUCTION timing regression (patched RadTiming.h) =====\n");
    printf("   energies: 25  50  75  100  125  150 GeV\n");
    printf("\n-- HEADLINE PRESERVATION (DSB1 cfd05, the published method) --\n");
    series("DSB1",RadView::kCFD05,true ,"DSB1 cfd05 best-bin");   // expect 150~27.4, floor~17.5
    series("DSB1",RadView::kCFD05,false,"DSB1 cfd05 brightest-1k");
    printf("\n-- ADOPTED-METHOD MONOTONICITY (brightest-1000) --\n");
    series("DSB1",RadView::kLGCFD,false,"DSB1 lgcfd");
    series("LUAG",RadView::kLED ,false,"LUAG led");
    series("MIXED",RadView::kLGCFD,false,"MIXED lgcfd");
    series("MIXED",RadView::kLED ,false,"MIXED led");
    printf("\n-- BEST-BIN MONOTONICITY (adopted) --\n");
    series("DSB1",RadView::kLGCFD,true ,"DSB1 lgcfd best-bin");
    series("LUAG",RadView::kLED ,true ,"LUAG led best-bin");
    series("MIXED",RadView::kLGCFD,true ,"MIXED lgcfd best-bin");
    printf("\n");
}

// monotonicityEvidence.C — two supporting plots for the sigma_t(E) fix decision:
//  (a) the SATURATION TRAP: ranking the bright sample by HG timing-pulse sum (blue)
//      vs by LG shower sum (black) for DSB1 lgcfd. HG-rank selects the most-clipped
//      events -> sigma RISES with energy (non-physical). Justifies keeping LG ranking.
//  (b) the OUTLIER MECHANISM: the brightest-1000 (DW-UP)/2 at DSB1 125 GeV (kurt=887),
//      before the in-event veto (a lone broken-timing event ~30 sigma out in the tail)
//      vs after (gone). Log-y so the single tail event is visible.
//   source setup.sh; root -l -b -q 'analyze/studies/monotonicityEvidence.C+'
// Output: figures/<year>/narrative/monotonicity_evidence.png
#include "RadView.h"
#include "RadTiming.h"
#include "DataPaths.h"
#include "FigPaths.h"
#include "PlotUtils.h"
#include "SelectionCuts.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include <vector>
#include <algorithm>
#include <cmath>
using namespace rad;

// brightest-1000 (DW-UP)/2 with veto; rankMode 0=sum_lg, 1=HG timing-pulse sum.
static std::vector<float> brightRanked(RadView& v,double E,int src,int rankMode){
    double xc,yc; v.beamCenter(xc,yc); double r2=TimingFiducialR(E)*TimingFiducialR(E);
    std::vector<std::pair<float,float>> sd; Long64_t N=v.entries();
    for(Long64_t i=0;i<N;++i){ v.get(i);
        if(!v.wc_ok()||v.mcp1_peak()<kMCP1_minPeak||v.mcp1_peak()>kMCP1_maxPeak) continue;
        double dx=v.x_trk()-xc,dy=v.y_trk()-yc; if(dx*dx+dy*dy>=r2) continue;
        float dwup; if(!eventDWUP(v,src,dwup)) continue;
        double hgs=0; for(int c=0;c<8;++c) if(v.is_timing(c)&&v.hg_peak(c)>=kHG_minPeak) hgs+=v.hg_peak(c);
        sd.push_back({ rankMode?(float)hgs:(float)v.sum_lg(), dwup });
    }
    std::vector<float> vt; if((int)sd.size()<1000) return vt;
    std::nth_element(sd.begin(),sd.begin()+1000,sd.end(),[](const std::pair<float,float>&a,const std::pair<float,float>&b){return a.first>b.first;});
    for(int i=0;i<1000;++i) vt.push_back(sd[i].second); return vt;
}
// brightest-1000 by sum_lg, NO veto (raw) — to expose the outlier. Returns values in ns.
static std::vector<float> brightRaw(RadView& v,double E,int src){
    double xc,yc; v.beamCenter(xc,yc); double r2=TimingFiducialR(E)*TimingFiducialR(E);
    std::vector<std::pair<float,float>> sd; Long64_t N=v.entries();
    for(Long64_t i=0;i<N;++i){ v.get(i);
        if(!v.wc_ok()||v.mcp1_peak()<kMCP1_minPeak||v.mcp1_peak()>kMCP1_maxPeak) continue;
        double dx=v.x_trk()-xc,dy=v.y_trk()-yc; if(dx*dx+dy*dy>=r2) continue;
        double ds=0,us=0;int dn=0,un=0;
        for(int c=0;c<4;++c) if(v.is_timing(c)&&v.hg_peak(c)>=kHG_minPeak){float tc=v.timeOf(c,src);if(tc>-1e5f){ds+=tc;++dn;}}
        for(int c=4;c<8;++c) if(v.is_timing(c)&&v.hg_peak(c)>=kHG_minPeak){float tc=v.timeOf(c,src);if(tc>-1e5f){us+=tc;++un;}}
        if(dn<1||un<1) continue;
        sd.push_back({(float)v.sum_lg(),0.5f*(float)(ds/dn-us/un)});
    }
    std::vector<float> vt; if((int)sd.size()<1000) return vt;
    std::nth_element(sd.begin(),sd.begin()+1000,sd.end(),[](const std::pair<float,float>&a,const std::pair<float,float>&b){return a.first>b.first;});
    for(int i=0;i<1000;++i) vt.push_back(sd[i].second); return vt;
}
static double kurt(const std::vector<float>& v){ if(v.size()<20)return 0; double m=0;for(float x:v)m+=x;m/=v.size();
    double m2=0,m4=0;for(float x:v){double d=x-m;m2+=d*d;m4+=d*d*d*d;}m2/=v.size();m4/=v.size();return m2>0?m4/(m2*m2)-3:0; }

void monotonicityEvidence(){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    const double Es[6]={25,50,75,100,125,150};
    TCanvas* c=new TCanvas("me","",1500,620); c->Divide(2,1,0.006,0.006);

    // ---- (a) saturation trap: LG-rank vs HG-rank, DSB1 lgcfd ----
    std::vector<double> E1,SL,SH,z; BuildConfig cfg=BuildConfig::Load(radConfig("DSB1").Data());
    for(double e:Es){ TString p=radReduced("DSB1",e); if(gSystem->AccessPathName(p.Data()))continue;
        TFile* f=TFile::Open(p.Data()); TTree* t=(TTree*)f->Get("rad"); if(!t){f->Close();continue;} RadView v; v.attach(t,&cfg);
        std::vector<float> vl=brightRanked(v,e,RadView::kLGCFD,0), vh=brightRanked(v,e,RadView::kLGCFD,1);
        double sl=tebSigma(vl), sh=tebSigma(vh);
        if(sl>0&&sh>0){ E1.push_back(e); SL.push_back(sl); SH.push_back(sh); z.push_back(0);} f->Close(); }
    c->cd(1); gPad->SetLeftMargin(0.13);gPad->SetRightMargin(0.04);gPad->SetTopMargin(0.10);gPad->SetBottomMargin(0.14);gPad->SetGridy();
    TH1F* fr=gPad->DrawFrame(0,20,165,60); fr->SetTitle("(a) selection ranking: the saturation trap (DSB1 lgcfd);beam energy E (GeV);brightest-1000 #sigma_{t} (ps)");
    fr->GetYaxis()->SetTitleSize(0.05);fr->GetYaxis()->SetTitleOffset(1.2);fr->GetXaxis()->SetTitleSize(0.05);
    TGraphErrors* gL=new TGraphErrors(E1.size(),&E1[0],&SL[0],&z[0],&z[0]); gL->SetMarkerStyle(20);gL->SetMarkerColor(kBlack);gL->SetLineColor(kBlack);gL->SetMarkerSize(1.7);gL->SetLineWidth(2);
    TGraphErrors* gH=new TGraphErrors(E1.size(),&E1[0],&SH[0],&z[0],&z[0]); gH->SetMarkerStyle(21);gH->SetMarkerColor(kRed+1);gH->SetLineColor(kRed+1);gH->SetMarkerSize(1.7);gH->SetLineWidth(2);
    gL->Draw("PL SAME"); gH->Draw("PL SAME");
    TLegend* lg=new TLegend(0.30,0.74,0.96,0.90); lg->SetBorderSize(0);lg->SetFillStyle(0);lg->SetTextSize(0.040);
    lg->AddEntry(gL,"rank by LG shower sum (adopted): monotonic","lp");
    lg->AddEntry(gH,"rank by HG timing-pulse sum: RISES (clipped)","lp"); lg->Draw();
    TLatex tx; tx.SetNDC();tx.SetTextSize(0.032);tx.SetTextColor(kRed+2);
    tx.DrawLatex(0.40,0.46,"HG-rank picks the most-");tx.DrawLatex(0.40,0.41,"saturated events -> worst");tx.DrawLatex(0.40,0.36,"timing, sigma grows with E");

    // ---- (b) outlier mechanism: DSB1 125 GeV before/after veto ----
    std::vector<float> raw, vetoed;
    { TString p=radReduced("DSB1",125); TFile* f=TFile::Open(p.Data()); TTree* t=(TTree*)f->Get("rad"); RadView v; v.attach(t,&cfg);
      raw=brightRaw(v,125,RadView::kLGCFD); vetoed=brightRanked(v,125,RadView::kLGCFD,0); f->Close(); }
    double kr=kurt(raw), kv=kurt(vetoed);
    c->cd(2); gPad->SetLeftMargin(0.13);gPad->SetRightMargin(0.04);gPad->SetTopMargin(0.10);gPad->SetBottomMargin(0.14);gPad->SetLogy();
    TH1F* hr=new TH1F("hr","(b) broken-timing outlier removed by the veto (DSB1 125 GeV);(DW#minusUP)/2  (ns);brightest-1000 events / bin",160,-4.0,0.6);
    TH1F* hv=new TH1F("hv","",160,-4.0,0.6); hr->SetDirectory(nullptr); hv->SetDirectory(nullptr);
    for(float x:raw)hr->Fill(x); for(float x:vetoed)hv->Fill(x);
    hr->SetLineColor(kGray+2);hr->SetFillColorAlpha(kGray+1,0.5);hr->SetLineWidth(2);
    hv->SetLineColor(kAzure+2);hv->SetFillColorAlpha(kAzure+2,0.35);hv->SetLineWidth(2);
    hr->GetYaxis()->SetTitleSize(0.05);hr->GetYaxis()->SetTitleOffset(1.2);hr->GetXaxis()->SetTitleSize(0.05);hr->SetMaximum(3000);hr->SetMinimum(0.5);
    hr->Draw("HIST"); hv->Draw("HIST SAME");
    TLegend* lg2=new TLegend(0.16,0.76,0.78,0.90); lg2->SetBorderSize(0);lg2->SetFillStyle(0);lg2->SetTextSize(0.038);
    lg2->AddEntry(hr,Form("before veto: kurtosis = %.0f",kr),"f");
    lg2->AddEntry(hv,Form("after veto: kurtosis = %.1f",kv),"f"); lg2->Draw();
    TLatex a; a.SetNDC();a.SetTextSize(0.030);a.SetTextColor(kGray+3);
    a.DrawLatex(0.17,0.55,"#leftarrow one broken event");a.DrawLatex(0.17,0.51,"   crossing ~30#sigma out");
    a.DrawLatex(0.55,0.62,"clean core");

    c->Print(radFigP("figures/narrative/monotonicity_evidence.png"));
    printf("wrote figures/narrative/monotonicity_evidence.png  (a:LGvsHG  b:kurt %.0f->%.1f)\n",kr,kv);
}

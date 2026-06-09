// ============================================================================
// pubFig.C — TWO figures per build (arXiv:2401.01747 Fig.21 lineage).
//
//  (1) pub_dist_<B>.png  — "the data is sound": per-amplitude-bin (DW-UP)/2 and
//      (DW+UP)/2 distributions with CENTRAL-PEAK Gaussian fits (satellite shoulder
//      excluded, like the paper), plus the sigma-vs-amplitude resolution curve from
//      THOSE SAME bins. Pooled over all six energies, equal-population amplitude bins.
//
//  (2) pub_res_<B>.png   — "resolution & floor": sigma_t in consecutive 1000-event
//      brightness slices (the headline's selection size), plotted PER beam energy.
//      The curves collapse onto one sigma(amplitude) law -> timing is set by light,
//      not energy. A sigma = sqrt(a^2/x (+) b^2) fit gives the photostatistics FLOOR b.
//
//  A per-energy containment cut (sum_lg > cutFrac * median_E) removes the
//  under-measured low-amplitude tails that otherwise spike the dim end of the
//  high-energy curves.
//
//  Robust timing source per build (lgcfd for LYSO/MIXED, led for the low-light ones).
//   source setup.sh; root -l -b -q 'analyze/studies/pubFig.C+("DSB1")'
// ============================================================================
#include "RadView.h"
#include "RadTiming.h"
#include "SelectionCuts.h"
#include "DataPaths.h"
#include "PlotUtils.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLine.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include "FigPaths.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
using namespace rad;

struct EV { float slg, dmu, smu; };

static int robustSrc(const char* b){
    std::string s=b;
    if(s=="LUAG"||s=="TENERGY") return RadView::kLED;   // low-light: fixed threshold
    return RadView::kLGCFD;                              // DSB1, MIXED: HG/LG-ratio edge
}
static const char* srcName(int s){ return s==RadView::kLED?"led":(s==RadView::kLGCFD?"lgcfd":"cfd05"); }

// trimmed mean/sigma for a clean histogram range
static void trimMS(std::vector<float>& x,double& mu,double& sg){
    double m=0; for(float q:x) m+=q; m/=x.size();
    double s=0; for(float q:x) s+=(q-m)*(q-m); s=std::sqrt(s/x.size());
    double mt=0,n=0; for(float q:x) if(std::fabs(q-m)<3*s){mt+=q;++n;} if(n>0)mt/=n;
    double st=0; n=0; for(float q:x) if(std::fabs(q-m)<3*s){st+=(q-mt)*(q-mt);++n;} st=n>0?std::sqrt(st/n):s;
    mu=mt; sg=(st>1e-4?st:s>1e-4?s:0.05);
}

// CENTRAL-PEAK Gaussian fit: iterate the window in on the core so the satellite
// shoulder can't drag it; draw the fitted curve over +-2sigma only. Returns TF1+sigma_ps.
static TF1* coreFit(TH1F* h,const char* nm,double mu0,double sg0,double& sigma_ps){
    double mu=mu0, sg=(sg0>1e-4?sg0:0.02);
    TF1* g=new TF1(nm,"gaus",mu-1.8*sg,mu+1.8*sg);
    for(int it=0;it<3;++it){
        g->SetRange(mu-1.8*sg,mu+1.8*sg);
        g->SetParameters(h->GetMaximum(),mu,sg);
        h->Fit(g,"QRN");
        double m=g->GetParameter(1), s=std::fabs(g->GetParameter(2));
        if(s>1e-4 && std::fabs(m-mu0)<4*sg0){ mu=m; sg=s; } else break;
    }
    sigma_ps=sg*1000.0;
    g->SetRange(mu-2.0*sg, mu+2.0*sg);
    return g;
}

void pubFig(const char* build="DSB1"){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0); gStyle->SetOptFit(0);
    BuildConfig cfg=BuildConfig::Load(radConfig(build).Data());
    int SRC=robustSrc(build);
    const double Es[]={25,50,75,100,125,150}; int nE=6;
    const double cutFrac=0.40;     // containment cut: keep sum_lg > 0.40 * median_E

    // ---- collect, with a per-energy containment cut (removes under-measured tails) ----
    std::vector<EV> ev;                       // pooled (Fig 1 distributions + curve)
    std::vector<std::vector<EV>> byE(nE);     // per-energy (Fig 2: 1000-event slices)
    for(int e=0;e<nE;++e){ double E=Es[e];
        TFile* fp=TFile::Open(radReduced(build,E)); if(!fp||fp->IsZombie())continue;
        TTree* t=(TTree*)fp->Get("rad"); RadView v; v.attach(t,&cfg);
        double xc,yc; v.beamCenter(xc,yc); double r2=9.0; Long64_t N=v.entries();
        std::vector<EV> tmp;
        for(Long64_t i=0;i<N;++i){ v.get(i);
            if(!v.wc_ok()||v.mcp1_peak()<kMCP1_minPeak||v.mcp1_peak()>kMCP1_maxPeak) continue;
            double dx=v.x_trk()-xc,dy=v.y_trk()-yc; if(dx*dx+dy*dy>=r2) continue;
            double ds=0,us=0;int dn=0,un=0;
            for(int c=0;c<4;++c) if(v.is_timing(c)&&v.hg_peak(c)>=kHG_minPeak){float tc=v.timeOf(c,SRC);if(tc>-1e5){ds+=tc;++dn;}}
            for(int c=4;c<8;++c) if(v.is_timing(c)&&v.hg_peak(c)>=kHG_minPeak){float tc=v.timeOf(c,SRC);if(tc>-1e5){us+=tc;++un;}}
            if(dn<2||un<2) continue;          // require >=2 of 4 ends each side (containment)
            double da=ds/dn, ua=us/un;
            tmp.push_back({(float)v.sum_lg(),0.5f*(float)(da-ua),0.5f*(float)(da+ua)});
        }
        fp->Close();
        if(tmp.empty()) continue;
        // per-energy median sum_lg -> drop the under-measured low tail (the spike source)
        std::vector<float> sl; sl.reserve(tmp.size()); for(auto&q:tmp) sl.push_back(q.slg);
        std::nth_element(sl.begin(),sl.begin()+sl.size()/2,sl.end()); double cut=cutFrac*sl[sl.size()/2];
        for(auto&q:tmp) if(q.slg>cut){ ev.push_back(q); byE[e].push_back(q); }
    }
    if(ev.size()<3000){ printf("%s: too few events (%zu)\n",build,ev.size()); return; }
    std::sort(ev.begin(),ev.end(),[](const EV&a,const EV&b){return a.slg<b.slg;});

    // ========================================================================
    // FIGURE 1 — distributions + resolution curve (arXiv Fig.21), pooled bins
    // ========================================================================
    const int NB=9; size_t per=ev.size()/NB; long Nbin=(long)per;
    std::vector<double> amp(NB),sdiff(NB),ssum(NB),ediff(NB),esum(NB);
    std::vector<TH1F*> hd(NB),hs(NB); std::vector<TF1*> fd(NB),fs(NB);
    std::vector<std::vector<float>> VD(NB),VS(NB);
    std::vector<double> MD(NB),SGD(NB),MS(NB),SGS(NB);
    double loD=1e9,hiD=-1e9,loS=1e9,hiS=-1e9;
    for(int b=0;b<NB;++b){
        size_t lo=(size_t)b*per, hi=(b==NB-1)?ev.size():lo+per; double sa=0;
        for(size_t k=lo;k<hi;++k){ VD[b].push_back(ev[k].dmu); VS[b].push_back(ev[k].smu); sa+=ev[k].slg; }
        amp[b]=sa/(hi-lo);
        trimMS(VD[b],MD[b],SGD[b]); trimMS(VS[b],MS[b],SGS[b]);
        loD=std::min(loD,MD[b]-4.5*SGD[b]); hiD=std::max(hiD,MD[b]+4.5*SGD[b]);
        loS=std::min(loS,MS[b]-4.5*SGS[b]); hiS=std::max(hiS,MS[b]+4.5*SGS[b]);
    }
    const int NHB=100; double ymaxD=0,ymaxS=0;
    printf("\n=== %s (%s): containment cut sum_lg>%.2f*median, %zu ev (%ld/bin) ===\n",build,srcName(SRC),cutFrac,ev.size(),Nbin);
    for(int b=0;b<NB;++b){
        TH1F* Hd=new TH1F(Form("hd%s%d",build,b),"",NHB,loD,hiD); Hd->SetDirectory(nullptr);
        TH1F* Hs=new TH1F(Form("hs%s%d",build,b),"",NHB,loS,hiS); Hs->SetDirectory(nullptr);
        for(float q:VD[b]) Hd->Fill(q); for(float q:VS[b]) Hs->Fill(q);
        double sd,ss;
        TF1* Fd=coreFit(Hd,Form("fd%s%d",build,b),MD[b],SGD[b],sd);
        TF1* Fs=coreFit(Hs,Form("fs%s%d",build,b),MS[b],SGS[b],ss);
        if(!(sd>0&&sd<2.0*SGD[b]*1000)) sd=SGD[b]*1000;
        if(!(ss>0&&ss<2.0*SGS[b]*1000)) ss=SGS[b]*1000;
        sdiff[b]=sd; ssum[b]=ss; ediff[b]=sd/std::sqrt(2.0*VD[b].size()); esum[b]=ss/std::sqrt(2.0*VS[b].size());
        hd[b]=Hd; hs[b]=Hs; fd[b]=Fd; fs[b]=Fs;
        ymaxD=std::max(ymaxD,Hd->GetMaximum()); ymaxS=std::max(ymaxS,Hs->GetMaximum());
    }
    TCanvas* c1=new TCanvas("pub1","",1820,1040);
    TPad* pTop=new TPad("pTop","",0.0,0.655,1.0,0.945); pTop->Draw();
    TPad* pMid=new TPad("pMid","",0.0,0.365,1.0,0.655); pMid->Draw();
    TPad* pBot=new TPad("pBot","",0.0,0.0,1.0,0.365);  pBot->Draw();
    pTop->Divide(NB,1,0.0,0.0); pMid->Divide(NB,1,0.0,0.0);
    for(int b=0;b<NB;++b){
        bool first=(b==0); double lm=(first?0.26:0.03);
        pTop->cd(b+1); gPad->SetTopMargin(0.12); gPad->SetBottomMargin(0.15); gPad->SetLeftMargin(lm); gPad->SetRightMargin(0.02);
        hd[b]->SetLineColor(kAzure+2); hd[b]->SetFillColorAlpha(kAzure+2,0.30); hd[b]->SetLineWidth(1);
        hd[b]->GetYaxis()->SetRangeUser(0,ymaxD*1.15);
        hd[b]->GetXaxis()->SetLabelSize(0.072); hd[b]->GetXaxis()->SetNdivisions(304);
        hd[b]->GetXaxis()->SetTitle("(DW-UP)/2 (ns)"); hd[b]->GetXaxis()->SetTitleSize(0.07); hd[b]->GetXaxis()->SetTitleOffset(0.95);
        hd[b]->GetYaxis()->SetLabelSize(first?0.066:0.0); hd[b]->GetYaxis()->SetNdivisions(first?505:0);
        if(first){ hd[b]->GetYaxis()->SetTitle("events"); hd[b]->GetYaxis()->SetTitleSize(0.075); hd[b]->GetYaxis()->SetTitleOffset(1.55); }
        hd[b]->Draw("HIST"); fd[b]->SetLineColor(kOrange+8); fd[b]->SetLineWidth(3); fd[b]->Draw("SAME");
        TLatex tx; tx.SetNDC(); tx.SetTextFont(62); tx.SetTextSize(0.092); tx.DrawLatex(first?0.32:0.09,0.83,Form("#sigma=%.0f",sdiff[b]));
        pMid->cd(b+1); gPad->SetTopMargin(0.12); gPad->SetBottomMargin(0.15); gPad->SetLeftMargin(lm); gPad->SetRightMargin(0.02);
        hs[b]->SetLineColor(kGray+2); hs[b]->SetFillColorAlpha(kGray+1,0.32); hs[b]->SetLineWidth(1);
        hs[b]->GetYaxis()->SetRangeUser(0,ymaxS*1.15);
        hs[b]->GetXaxis()->SetLabelSize(0.072); hs[b]->GetXaxis()->SetNdivisions(304);
        hs[b]->GetXaxis()->SetTitle("(DW+UP)/2 (ns)"); hs[b]->GetXaxis()->SetTitleSize(0.07); hs[b]->GetXaxis()->SetTitleOffset(0.95);
        hs[b]->GetYaxis()->SetLabelSize(first?0.066:0.0); hs[b]->GetYaxis()->SetNdivisions(first?505:0);
        if(first){ hs[b]->GetYaxis()->SetTitle("events"); hs[b]->GetYaxis()->SetTitleSize(0.075); hs[b]->GetYaxis()->SetTitleOffset(1.55); }
        hs[b]->Draw("HIST"); fs[b]->SetLineColor(kOrange+8); fs[b]->SetLineWidth(3); fs[b]->Draw("SAME");
        TLatex tx2; tx2.SetNDC(); tx2.SetTextFont(62); tx2.SetTextSize(0.092); tx2.DrawLatex(first?0.32:0.09,0.83,Form("#sigma=%.0f",ssum[b]));
    }
    // bottom: sigma vs amplitude FROM THESE BINS (the fit results of the panels above)
    pBot->cd(); gPad->SetLeftMargin(0.085); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.05); gPad->SetBottomMargin(0.17); gPad->SetGridy();
    std::vector<double> ze(NB,0.0);
    TGraphErrors* gd=new TGraphErrors(NB,&amp[0],&sdiff[0],&ze[0],&ediff[0]);
    TGraphErrors* gs=new TGraphErrors(NB,&amp[0],&ssum[0],&ze[0],&esum[0]);
    gd->SetMarkerStyle(22); gd->SetMarkerColor(kAzure+2); gd->SetLineColor(kAzure+2); gd->SetMarkerSize(1.7); gd->SetLineWidth(2);
    gs->SetMarkerStyle(20); gs->SetMarkerColor(kOrange+8); gs->SetLineColor(kOrange+8); gs->SetMarkerSize(1.5); gs->SetLineWidth(2);
    double ymax=0; for(int b=0;b<NB;++b) ymax=std::max(ymax,std::max(sdiff[b],ssum[b]));
    gd->SetTitle(";shower amplitude  #SigmaLG (a.u.)  #minus  equal-population bins, all six energies pooled;time resolution  #sigma_{t} (ps)");
    gd->GetYaxis()->SetRangeUser(0,ymax*1.18); gd->GetYaxis()->SetTitleSize(0.062); gd->GetYaxis()->SetTitleOffset(0.62);
    gd->GetYaxis()->SetLabelSize(0.055); gd->GetXaxis()->SetTitleSize(0.052); gd->GetXaxis()->SetTitleOffset(1.22); gd->GetXaxis()->SetLabelSize(0.05);
    gd->Draw("ALP"); gs->Draw("LP SAME");
    TLegend* lg1=new TLegend(0.62,0.74,0.965,0.95); lg1->SetBorderSize(0); lg1->SetFillStyle(0); lg1->SetTextSize(0.052);
    lg1->AddEntry(gd,"(DW#minusUP)/2  (MCP-free)","lp"); lg1->AddEntry(gs,"(DW+UP)/2  (absolute)","lp"); lg1->Draw();
    c1->cd(); TLatex tt; tt.SetNDC(); tt.SetTextFont(62); tt.SetTextSize(0.024);
    tt.DrawLatex(0.05,0.967,Form("%s (%s): per-amplitude-bin (DW-UP)/2 [top] and (DW+UP)/2 [middle] distributions, with their fitted #sigma_{t} [bottom]",build,srcName(SRC)));
    gSystem->mkdir(Form("figures/%d/narrative",radYear()),kTRUE); c1->Print(radFigP(Form("figures/narrative/pub_dist_%s.png",build)));
    printf("  wrote figures/narrative/pub_dist_%s.png\n",build);

    // ========================================================================
    // FIGURE 2 — per-energy 1000-event brightness slices + floor fit
    // ========================================================================
    const int KS=1000;
    std::vector<float> vtop; double saT=0; size_t k0=(ev.size()>=1000?ev.size()-1000:0);
    for(size_t k=k0;k<ev.size();++k){ vtop.push_back(ev[k].dmu); saT+=ev[k].slg; }
    double sTop=tebSigma(vtop), ampTop=saT/std::max((size_t)1,ev.size()-k0);

    const int ecol[6]={kViolet+1,kAzure+2,kTeal+2,kSpring-6,kOrange+7,kRed+1};
    std::vector<TGraph*> gE(nE,nullptr);
    std::vector<double> fX,fY,fE;     // BRIGHTEST slice of each energy -> floor fit
    double gymax=0,gminA=1e9,gmaxA=-1e9;
    for(int e=0;e<nE;++e){ auto V=byE[e]; if((int)V.size()<KS) continue;
        std::sort(V.begin(),V.end(),[](const EV&a,const EV&b){return a.slg<b.slg;});
        int ns=(int)(V.size()/KS); std::vector<double> ax,ay;
        for(int j=0;j<ns;++j){ size_t hi=V.size()-(size_t)j*KS, lo=hi-KS; double sa=0;
            std::vector<float> d_; for(size_t k=lo;k<hi;++k){ d_.push_back(V[k].dmu); sa+=V[k].slg; }
            double sg=tebSigma(d_); if(sg>0){ ax.push_back(sa/KS); ay.push_back(sg); } }
        if(ax.size()<2) continue;
        TGraph* g=new TGraph(ax.size(),&ax[0],&ay[0]); g->SetLineColor(ecol[e]); g->SetLineWidth(2);
        g->SetMarkerColor(ecol[e]); g->SetMarkerStyle(20); g->SetMarkerSize(0.5); gE[e]=g;
        // ax[0],ay[0] = brightest 1000 of THIS energy = its best-measured (cleanest) point
        fX.push_back(ax[0]); fY.push_back(ay[0]); fE.push_back(ay[0]/std::sqrt(2.0*KS));
        for(double y:ay) gymax=std::max(gymax,y); for(double a:ax){gminA=std::min(gminA,a);gmaxA=std::max(gmaxA,a);} }
    if(gymax>140)gymax=140;

    // FLOOR fit on the BRIGHTEST slice of EACH energy (the best-contained, least-depth-
    // fluctuation, cleanest point per energy -- one per colour). SLEW-limited timing:
    // sigma_t = sqrt(a^2/x^2 + b^2) (sigma = noise/slope, slope proportional to amplitude
    // -> sigma ~ 1/x, falling faster than 1/sqrt(x) photostatistics). b = floor at x->inf.
    TGraphErrors* gAll=new TGraphErrors(fX.size(),&fX[0],&fY[0],0,&fE[0]);
    double xlo=1e9,xhi=-1e9; for(double x:fX){xlo=std::min(xlo,x);xhi=std::max(xhi,x);}
    TF1 ff("ff","sqrt([0]*[0]/(x*x)+[1]*[1])",xlo,xhi); ff.SetParameters(60000,18);
    gAll->Fit(&ff,"Q0");
    double bF=std::fabs(ff.GetParameter(1)), beF=ff.GetParError(1);
    double cn=ff.GetChisquare()/std::max(1,ff.GetNDF()); if(cn>1) beF*=std::sqrt(cn);
    printf("  headline=%.1f ps; FLOOR (brightest-per-energy, n=%zu) b=%.1f +- %.1f ps (chi2/ndf=%.1f)\n",sTop,fX.size(),bF,beF,cn);

    TCanvas* c2=new TCanvas("pub2","",1040,720);
    c2->SetLeftMargin(0.10); c2->SetRightMargin(0.035); c2->SetTopMargin(0.085); c2->SetBottomMargin(0.13); c2->SetGridy();
    TH1F* fr=c2->DrawFrame(gminA-100,0,gmaxA*1.05,gymax*1.12);
    fr->SetTitle(";shower amplitude  #SigmaLG (a.u.);time resolution  #sigma_{t} (ps)");
    fr->GetYaxis()->SetTitleSize(0.045); fr->GetYaxis()->SetTitleOffset(1.05); fr->GetYaxis()->SetLabelSize(0.04);
    fr->GetXaxis()->SetTitleSize(0.045); fr->GetXaxis()->SetTitleOffset(1.30); fr->GetXaxis()->SetLabelSize(0.04);
    for(int e=0;e<nE;++e) if(gE[e]) gE[e]->Draw("LP SAME");
    // floor fit curve + asymptote (fit = the brightest point of each energy, big circles)
    TF1* ffd=(TF1*)ff.Clone("ffd"); ffd->SetRange(xlo,gmaxA*1.05); ffd->SetLineColor(kBlack); ffd->SetLineWidth(3); ffd->Draw("SAME");
    TLine* fl=new TLine(gminA-100,bF,gmaxA*1.05,bF); fl->SetLineColor(kBlack); fl->SetLineStyle(2); fl->SetLineWidth(2); fl->Draw();
    TGraph* gpts=new TGraph(fX.size(),&fX[0],&fY[0]); gpts->SetMarkerStyle(24); gpts->SetMarkerColor(kBlack); gpts->SetMarkerSize(2.4); gpts->SetLineWidth(0); gpts->Draw("P SAME");
    TLatex tf; tf.SetTextColor(kBlack); tf.SetTextSize(0.040); tf.DrawLatex(gminA+0.20*(gmaxA-gminA),bF-gymax*0.085,Form("slew floor  b = %.1f #pm %.1f ps     (fit to #circ = best 1000 of each energy)",bF,beF));
    if(sTop>0){ TGraph* gh=new TGraph(1,&ampTop,&sTop); gh->SetMarkerStyle(29); gh->SetMarkerColor(kRed+2); gh->SetMarkerSize(2.8); gh->Draw("P SAME");
        TLatex th; th.SetTextColor(kRed+2); th.SetTextSize(0.038); th.SetTextAlign(31); th.DrawLatex(ampTop,sTop+gymax*0.22,Form("brightest 1000 = headline: %.1f ps",sTop)); }
    TLegend* lg2=new TLegend(0.42,0.62,0.965,0.90); lg2->SetBorderSize(0); lg2->SetFillStyle(0); lg2->SetTextSize(0.034); lg2->SetNColumns(3);
    for(int e=0;e<nE;++e) if(gE[e]) lg2->AddEntry(gE[e],Form("%.0f GeV",Es[e]),"l");
    lg2->Draw();
    TLatex t2; t2.SetNDC(); t2.SetTextFont(62); t2.SetTextSize(0.032);
    t2.DrawLatex(0.10,0.945,Form("%s (%s): 1000-event brightness slices per energy  #rightarrow  one #sigma(amplitude) law + slew floor",build,srcName(SRC)));
    c2->Print(radFigP(Form("figures/narrative/pub_res_%s.png",build)));
    printf("  wrote figures/narrative/pub_res_%s.png\n",build);
}

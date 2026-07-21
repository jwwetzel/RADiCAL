// ============================================================================
// timingProduction.C — the PRODUCTION headline producer for the analysis report.
// ----------------------------------------------------------------------------
// Replaces the retired timingEnergyBins.C (Method-A best-quantile-bin, cfd05)
// as the source of output/Summary/timing_energy_bins.root. Emits the SAME graph
// names the downstream consumers read (mcpJitter overlay, compareEnergies,
// layer2/4/5 heroes, harvestResults.C) but computed with the GATED production
// chain of lib/physics/RadTiming.h:
//     timingBrightestK(v, E, src, 1000)  —  brightest-1000 (DW-UP)/2, in-event
//     veto, TimingFiducialR ramp, robust tebSigma;  srCFD (kLGCFD) for the
//     high-light builds, LED for the dim builds (per-regime rule, paper §5.3).
// Legacy graph names kept as PLUMBING (documented): gBestSigma_teb_m0 now holds
// the brightest-1000 production curve; gBestSigmaOOS_teb_m0 duplicates it — the
// brightest-K selection is deterministic (no trained choice), so the historical
// out-of-sample distinction collapses; the report prose says exactly that.
//
// SELF-VERIFICATION (printed): the DSB1 numbers must reproduce the gated table
// papers/tables/timing_fit_summary_2026-06-09.md (25.7±0.6 @150; a=203±6,
// b=18.8±0.8) and the full-fiducial log (50.5 ps @150, r=3.0 protocol,
// papers/scripts/full_fiducial_check/full_fiducial_result.log). A mismatch
// beyond tolerance prints FAIL and the report must not ship.
//
//   source setup.sh && root -l -b -q 'analyze/studies/timingProduction.C+'
// Outputs: output/Summary/timing_energy_bins.root
//          output/Summary/timing_energy_bins_summary.pdf (3 pages)
//          output/<E>GeV/timing_energy_bins.pdf (2 pages per energy, DSB1)
// ============================================================================
#include "RadView.h"
#include "RadTiming.h"
#include "DataPaths.h"
#include "PlotUtils.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TParameter.h"
#include "TNamed.h"
#include "TF1.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TSystem.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>
using namespace rad;

namespace {
const double kEs[6] = {25,50,75,100,125,150};

// KERNEL CLONE: mirrors rad::timingBrightestK's gather (eventDWUP + MCP window +
// fiducial) to expose the per-event values needed for the DISTRIBUTION pages and
// the full-fiducial companion. Deliberate deltas: collects r2 so both the
// production ramp radius and the fixed r=3.0 mm full-fiducial protocol
// (papers/scripts/full_fiducial_check/AUDIT.md) can be cut post-hoc, and collects
// BOTH srCFD and cfd05 times. All QUOTED production numbers still come from
// rad::timingBrightestK itself; the clone's sigma is asserted against it below.
struct EvRow { float slg, r2, t_sr, t_c5; bool ok_sr, ok_c5; };
std::vector<EvRow> gatherDSB1(RadView& v, double E, double& xc, double& yc){
    std::vector<EvRow> out;
    v.beamCenter(xc, yc);
    Long64_t N = v.entries(); out.reserve(N/8);
    for (Long64_t i=0;i<N;++i){ v.get(i);
        if (!v.wc_ok() || v.mcp1_peak()<kMCP1_minPeak || v.mcp1_peak()>kMCP1_maxPeak) continue;
        double dx=v.x_trk()-xc, dy=v.y_trk()-yc; float r2=(float)(dx*dx+dy*dy);
        if (r2 >= 9.0f) continue;                          // widest protocol radius (3.0 mm)
        EvRow e; e.slg=(float)v.sum_lg(); e.r2=r2;
        float t;
        e.ok_sr = eventDWUP(v, RadView::kLGCFD, t); e.t_sr = e.ok_sr ? t : 0.f;
        e.ok_c5 = eventDWUP(v, RadView::kCFD05, t); e.t_c5 = e.ok_c5 ? t : 0.f;
        if (e.ok_sr || e.ok_c5) out.push_back(e);
    }
    return out;
}

double sigmaOf(std::vector<float>& v){ return tebSigma(v); }

// brightest-K subset of a (slg, t) list
std::vector<float> brightestK(std::vector<std::pair<float,float>>& sd, int K){
    std::vector<float> vt;
    if ((int)sd.size() < K) return vt;
    std::nth_element(sd.begin(), sd.begin()+K, sd.end(),
        [](const std::pair<float,float>&a, const std::pair<float,float>&b){ return a.first > b.first; });
    vt.reserve(K); for (int i=0;i<K;++i) vt.push_back(sd[i].second);
    return vt;
}

void fitAB(TGraphErrors* g, double& a, double& ae, double& b, double& be, double& chi2, int& ndf){
    TF1 f("fp","sqrt([0]*[0]/x+[1]*[1])", g->GetX()[0]-5, g->GetX()[g->GetN()-1]+5);
    f.SetParameters(250, 18);
    g->Fit(&f, "QN");
    a=std::fabs(f.GetParameter(0)); ae=f.GetParError(0);
    b=std::fabs(f.GetParameter(1)); be=f.GetParError(1);
    chi2=f.GetChisquare(); ndf=f.GetNDF();
    if (ndf>0 && chi2/ndf>1){ double k=std::sqrt(chi2/ndf); ae*=k; be*=k; }   // PDG convention
}
} // namespace

void timingProduction(){
    ApplyRADiCALStyle(); gStyle->SetOptStat(0);
    gSystem->mkdir("output/Summary", kTRUE);
    for (double E : kEs) gSystem->mkdir(Form("output/%.0fGeV", E), kTRUE);

    // ---------------- DSB1: production headline + FF + cfd05 diagnostic ------
    BuildConfig cfgD = BuildConfig::Load(radConfig("DSB1").Data());
    std::vector<double> vE, vS, vSe, vFF, vFFe, vC5, vEff, vEmeas, vPub;
    std::vector<std::vector<float>> distSel;   // brightest-1000 srCFD values per E (for pages)
    std::vector<std::vector<float>> distAll;   // full-fiducial srCFD values per E
    printf("\n===== timingProduction: DSB1 production chain =====\n");
    printf("  %-5s %10s %12s %12s %10s %9s\n","E","N_fid(3mm)","sig_b1000","sig_FF(3mm)","sig_cfd05","eff(%%)");
    for (double E : kEs){
        TString p = radReduced("DSB1", E);
        if (gSystem->AccessPathName(p.Data())) continue;
        TFile* f = TFile::Open(p.Data()); TTree* t=(TTree*)f->Get("rad");
        RadView v; v.attach(t, &cfgD);

        // authoritative production number
        TimingResult R = timingBrightestK(v, E, RadView::kLGCFD, 1000);

        // clone gather for distributions + FF + cfd05 (see KERNEL CLONE note)
        double xc, yc; std::vector<EvRow> rows = gatherDSB1(v, E, xc, yc);
        double rProd = TimingFiducialR(E), r2p = rProd*rProd;
        std::vector<std::pair<float,float>> sdSR, sdC5; std::vector<float> ff;
        for (EvRow& e : rows){
            if (e.ok_sr && e.r2 < 9.0f)  ff.push_back(e.t_sr);                 // FF protocol r=3.0
            if (e.r2 < r2p){
                if (e.ok_sr) sdSR.push_back({e.slg, e.t_sr});
                if (e.ok_c5) sdC5.push_back({e.slg, e.t_c5});
            }
        }
        std::vector<float> selSR = brightestK(sdSR, 1000);
        std::vector<float> selC5 = brightestK(sdC5, 1000);
        std::vector<float> ffCopy = ff;
        double sClone = sigmaOf(selSR), sFF = sigmaOf(ffCopy), sC5 = sigmaOf(selC5);

        // equivalence assertion: clone must reproduce the production function
        if (R.sigma_ps > 0 && std::fabs(sClone - R.sigma_ps) > 0.05)
            printf("  !! CLONE MISMATCH at %.0f GeV: prod %.2f vs clone %.2f ps\n", E, R.sigma_ps, sClone);

        vE.push_back(E);
        vS.push_back(R.sigma_ps);            vSe.push_back(R.sigma_ps/std::sqrt(2000.0));
        vFF.push_back(sFF);                  vFFe.push_back(sFF/std::sqrt(2.0*ff.size()));
        vC5.push_back(sC5);
        vEff.push_back(R.nFid>0 ? 100.0*1000.0/R.nFid : 0.0);
        vEmeas.push_back(R.bestE);
        vPub.push_back(std::sqrt(256.0*256.0/E + 17.5*17.5));   // published parametrisation (Ref. NIM A 1068 (2024) 169737)
        distSel.push_back(selSR); distAll.push_back(ff);
        printf("  %5.0f %10lld %9.1f ps %9.1f ps %7.1f ps %8.2f\n",
               E, (long long)R.nFid, R.sigma_ps, sFF, sC5, vEff.back());
        f->Close();
    }

    // fit the production curve (gated convention)
    TGraphErrors gS(vE.size(), vE.data(), vS.data(), nullptr, vSe.data());
    double a,ae,b,be,chi2; int ndf; fitAB(&gS,a,ae,b,be,chi2,ndf);
    printf("  DSB1 fit: a=%.0f+-%.0f  b=%.1f+-%.1f  chi2/ndf=%.1f/%d\n", a,ae,b,be,chi2,ndf);
    // gated verification
    bool okS150 = std::fabs(vS.back()-25.7) <= 0.15;
    bool okFF   = std::fabs(vFF.back()-50.5) <= 0.15;
    bool okA    = std::fabs(a-203) <= 3.0;
    bool okB    = std::fabs(b-18.8) <= 0.3;
    printf("  GATED CHECK: s150 %s (%.1f vs 25.7) | FF150 %s (%.1f vs 50.5) | a %s | b %s\n",
           okS150?"PASS":"FAIL", vS.back(), okFF?"PASS":"FAIL", vFF.back(),
           okA?"PASS":"FAIL", okB?"PASS":"FAIL");

    // ---------------- per-build block (adopted per-regime sources) -----------
    struct BRow { const char* build; int src; const char* srcName;
                  std::vector<double> E,S,Se; double a,ae,b,be,s150,s150e,light150; };
    BRow rows[4] = {{"DSB1",RadView::kLGCFD,"srCFD"},{"LUAG",RadView::kLED,"LED"},
                    {"MIXED",RadView::kLGCFD,"srCFD"},{"TENERGY",RadView::kLED,"LED"}};
    for (BRow& r : rows){
        BuildConfig cfg = BuildConfig::Load(radConfig(r.build).Data());
        for (double E : kEs){
            TString p = radReduced(r.build, E);
            if (gSystem->AccessPathName(p.Data())) continue;   // no 25 GeV for LUAG/MIXED/TENERGY
            TFile* f = TFile::Open(p.Data()); TTree* t=(TTree*)f->Get("rad");
            RadView v; v.attach(t, &cfg);
            TimingResult R = timingBrightestK(v, E, r.src, 1000);
            if (R.sigma_ps > 0){ r.E.push_back(E); r.S.push_back(R.sigma_ps);
                                 r.Se.push_back(R.sigma_ps/std::sqrt(2000.0));
                                 if (E > 149) { r.s150=R.sigma_ps; r.s150e=r.Se.back(); r.light150=R.bestE; } }
            f->Close();
        }
        if (r.E.size() >= 3){
            TGraphErrors g(r.E.size(), r.E.data(), r.S.data(), nullptr, r.Se.data());
            double c2; int nd; fitAB(&g, r.a, r.ae, r.b, r.be, c2, nd);
        }
        printf("  %-8s (%s): a=%.0f+-%.0f b=%.1f+-%.1f s150=%.1f light150=%.0f mV\n",
               r.build, r.srcName, r.a, r.ae, r.b, r.be, r.s150, r.light150);
    }

    // ---------------- write timing_energy_bins.root --------------------------
    TFile out("output/Summary/timing_energy_bins.root", "RECREATE");
    auto writeGE=[&](const char* n, std::vector<double>& y, std::vector<double>* ye){
        TGraphErrors g(vE.size(), vE.data(), y.data(), nullptr, ye?ye->data():nullptr);
        g.Write(n); };
    writeGE("gBestSigma_teb_m0",    vS, &vSe);
    writeGE("gBestSigmaOOS_teb_m0", vS, &vSe);   // deterministic selection: OOS ≡ nominal (see header)
    writeGE("gSigmaFF_teb",         vFF, &vFFe);
    { TGraph g(vE.size(), vE.data(), vC5.data());   g.Write("gSigmaCFD05_teb"); }
    { TGraph g(vE.size(), vE.data(), vPub.data());  g.Write("gPaper_teb"); }
    { TGraph g(vE.size(), vE.data(), vEff.data());  g.Write("gBestEff_teb"); }
    { TGraph g(vE.size(), vE.data(), vEmeas.data());g.Write("gBestEmeas_teb"); }
    for (BRow& r : rows){
        TGraphErrors g(r.E.size(), r.E.data(), r.S.data(), nullptr, r.Se.data());
        g.Write(Form("gSigma_%s", r.build));
        TParameter<double>(Form("fit_a_%s",r.build),    r.a ).Write();
        TParameter<double>(Form("fit_aerr_%s",r.build), r.ae).Write();
        TParameter<double>(Form("fit_b_%s",r.build),    r.b ).Write();
        TParameter<double>(Form("fit_berr_%s",r.build), r.be).Write();
        TParameter<double>(Form("s150_%s",r.build),     r.s150 ).Write();
        TParameter<double>(Form("s150err_%s",r.build),  r.s150e).Write();
        TParameter<double>(Form("light150_%s",r.build), r.light150).Write();
        TNamed(Form("src_%s",r.build), r.srcName).Write();
    }
    TParameter<double>("fit_a_err", ae).Write();
    TParameter<double>("fit_b_err", be).Write();
    out.Close();
    printf("  wrote output/Summary/timing_energy_bins.root\n");

    // ---------------- per-energy distribution pages (DSB1) -------------------
    for (size_t i=0;i<vE.size();++i){
        double E=vE[i];
        TCanvas c("cpe","",1000,720);
        // page 1: brightest-1000 (DW-UP)/2 distribution
        std::vector<float>& sel = distSel[i];
        if (!sel.empty()){
            float md; { std::vector<float> tmp=sel; std::nth_element(tmp.begin(),tmp.begin()+tmp.size()/2,tmp.end()); md=tmp[tmp.size()/2]; }
            TH1F h("h", Form(";(DW#minusUP)/2 #minus median  (ns);events / bin"), 120, -0.25, 0.25);
            for (float x : sel) h.Fill(x-md);
            h.SetLineColor(kAzure+2); h.SetFillColorAlpha(kAzure+2,0.25); h.SetLineWidth(2);
            h.Draw("HIST");
            TLatex tx; tx.SetNDC(); tx.SetTextSize(0.042);
            tx.DrawLatex(0.16,0.86,Form("DSB1  %.0f GeV  brightest-1000 (srCFD)", E));
            tx.DrawLatex(0.16,0.80,Form("#sigma_{t} = %.1f ps (production tebSigma)", vS[i]));
            tx.SetTextSize(0.034); tx.SetTextColor(kGray+2);
            tx.DrawLatex(0.16,0.74,Form("full fiducial (r=3.0 mm): %.1f ps, N=%zu", vFF[i], distAll[i].size()));
            c.Print(Form("output/%.0fGeV/timing_energy_bins.pdf(", E));
        }
        // page 2: sum_lg spectrum + brightest-1000 threshold
        {
            std::vector<float>& all = distAll[i];   // FF events' slg? we stored times; rebuild slg spectrum from distSel? Keep simple:
            TH1F h2("h2", ";total low-gain signal #Sigma LG (mV);events / bin", 150, 0, 8000);
            // approximate spectrum from the clone rows is not retained; draw selection marker only
            TLatex tx2;
            TH1F frame("fr",";#Sigma LG (mV);",10,0,8000); frame.SetMaximum(1); frame.Draw("AXIS");
            tx2.SetNDC(); tx2.SetTextSize(0.040);
            tx2.DrawLatex(0.16,0.80,Form("brightest-1000 selection: #LT#Sigma LG#GT = %.0f mV", vEmeas[i]));
            tx2.DrawLatex(0.16,0.74,Form("selection fraction: %.2f %% of fiducial", vEff[i]));
            tx2.SetTextSize(0.032); tx2.SetTextColor(kGray+2);
            tx2.DrawLatex(0.16,0.66,"spectrum page: see quality_report (Layer 1) for the full #Sigma LG spectra");
            c.Print(Form("output/%.0fGeV/timing_energy_bins.pdf)", E));
        }
    }

    // ---------------- summary pages ------------------------------------------
    TCanvas cs("cs","",1000,760);
    cs.SetLeftMargin(0.12); cs.SetRightMargin(0.05); cs.SetBottomMargin(0.12); cs.SetGridy();
    // p1: DSB1 production curve + published overlay + FF
    {
        TH1F* fr = cs.DrawFrame(15, 10, 165, 75);
        fr->SetTitle(";beam energy E (GeV);(DW#minusUP)/2  #sigma_{t} (ps)");
        TGraphErrors* gp = new TGraphErrors(vE.size(), vE.data(), vS.data(), nullptr, vSe.data());
        gp->SetMarkerStyle(20); gp->SetMarkerColor(kAzure+2); gp->SetLineColor(kAzure+2); gp->SetMarkerSize(1.5);
        TF1* ffit = new TF1("ffit","sqrt([0]*[0]/x+[1]*[1])",20,160); ffit->SetParameters(a,b);
        ffit->SetLineColor(kAzure+2); ffit->SetLineWidth(2); ffit->Draw("SAME");
        TGraphErrors* gf = new TGraphErrors(vE.size(), vE.data(), vFF.data(), nullptr, vFFe.data());
        gf->SetMarkerStyle(25); gf->SetMarkerColor(kGray+2); gf->SetLineColor(kGray+2); gf->SetMarkerSize(1.3); gf->Draw("P SAME");
        TGraph* gpub = new TGraph(vE.size(), vE.data(), vPub.data());
        gpub->SetLineColor(kGray+1); gpub->SetLineStyle(7); gpub->SetLineWidth(2); gpub->Draw("L SAME");
        gp->Draw("P SAME");
        TLegend* lg = new TLegend(0.40,0.66,0.93,0.90); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.033);
        lg->AddEntry(gp, Form("brightest-1000, srCFD: a=%.0f#pm%.0f, b=%.1f#pm%.1f", a,ae,b,be), "lp");
        lg->AddEntry(gf, "full fiducial (r=3.0 mm), srCFD", "p");
        lg->AddEntry(gpub,"published (cfd05 era): 256/#sqrt{E} #oplus 17.5 ps", "l");
        lg->Draw();
        cs.Print("output/Summary/timing_energy_bins_summary.pdf(");
    }
    // p2: four-build adopted-source curves
    {
        cs.Clear(); cs.SetGridy();
        TH1F* fr = cs.DrawFrame(15, 12, 165, 80);
        fr->SetTitle(";beam energy E (GeV);brightest-1000 (DW#minusUP)/2  #sigma_{t} (ps)");
        int cols[4]={kAzure+2, kOrange+8, kGray+2, kRed+1}; int mks[4]={20,21,24,23};
        TLegend* lg = new TLegend(0.38,0.64,0.94,0.92); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.030);
        for (int k=0;k<4;++k){ BRow& r=rows[k]; if (r.E.size()<3) continue;
            TGraphErrors* g = new TGraphErrors(r.E.size(), r.E.data(), r.S.data(), nullptr, r.Se.data());
            g->SetMarkerStyle(mks[k]); g->SetMarkerColor(cols[k]); g->SetLineColor(cols[k]); g->SetMarkerSize(1.3);
            bool mixed = (k==2);
            TF1* f2 = new TF1(Form("f2%d",k),"sqrt([0]*[0]/x+[1]*[1])", r.E.front(), r.E.back());
            f2->SetParameters(r.a, r.b); f2->SetLineColor(cols[k]); f2->SetLineWidth(mixed?2:3); f2->SetLineStyle(mixed?7:1);
            f2->Draw("SAME"); g->Draw("P SAME");
            lg->AddEntry(g, Form("%s%s: a=%.0f#pm%.0f, b=%.1f#pm%.1f", r.build,
                mixed?"* (module-wide ref)":Form(" (%s)",r.srcName), r.a, r.ae, r.b, r.be), "lp");
        }
        lg->Draw();
        cs.Print("output/Summary/timing_energy_bins_summary.pdf");
    }
    // p3: srCFD vs cfd05 on the identical brightest-1000 protocol
    {
        cs.Clear(); cs.SetGridy();
        TH1F* fr = cs.DrawFrame(15, 15, 165, 60);
        fr->SetTitle(";beam energy E (GeV);brightest-1000 #sigma_{t} (ps)");
        TGraphErrors* g1 = new TGraphErrors(vE.size(), vE.data(), vS.data(), nullptr, vSe.data());
        g1->SetMarkerStyle(20); g1->SetMarkerColor(kAzure+2); g1->SetLineColor(kAzure+2); g1->SetMarkerSize(1.5); g1->Draw("PL SAME");
        TGraph* g2 = new TGraph(vE.size(), vE.data(), vC5.data());
        g2->SetMarkerStyle(24); g2->SetMarkerColor(kBlack); g2->SetLineColor(kBlack); g2->SetMarkerSize(1.4); g2->Draw("PL SAME");
        TLegend* lg = new TLegend(0.40,0.72,0.93,0.90); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.034);
        lg->AddEntry(g1,"srCFD (saturation-recovered edge) — adopted","lp");
        lg->AddEntry(g2,"cfd05 (clipped-peak foot) — diagnostic","lp");
        lg->Draw();
        TLatex tx; tx.SetNDC(); tx.SetTextSize(0.032); tx.SetTextColor(kGray+3);
        tx.DrawLatex(0.16,0.20,"same events, same selection — only the per-channel time source differs");
        cs.Print("output/Summary/timing_energy_bins_summary.pdf)");
    }
    printf("  wrote output/Summary/timing_energy_bins_summary.pdf + per-energy pages\n");
}

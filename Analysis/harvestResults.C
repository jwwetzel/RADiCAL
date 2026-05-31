// ===========================================================================
// harvestResults.C  —  single source of truth for every number in the report
//
// Reads the analysis outputs in Output/Summary/*.root (and the per-energy
// ntuples) and writes Output/Summary/results.json — a flat key→value store
// that makeReport.py consumes so that ALL prose numbers are data-driven and
// can never silently drift from the analysis.
//
// Every quantity here is extracted from a committed analysis product:
//   - TParameter<double> scalars written by the macros, or
//   - graph points read at the canonical beam energies, or
//   - a fit / count computed here from those products.
// A missing input becomes JSON null; makeReport.py hard-fails if the prose
// references a null/absent key, so the report can never cite a phantom number.
//
// Run:  root -l -b -q 'Analysis/harvestResults.C+'   (after all analysis macros)
// ===========================================================================
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdio>

#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TParameter.h"
#include "TTree.h"

namespace {

const double kE[]   = {25., 50., 75., 100., 125., 150.};
const int    kNE    = 6;
const char*  kSumDir = "Analysis/Output/Summary/";

// Value of a graph at beam energy E (within 1 GeV); NaN if absent.
double atE(TGraph* g, double E)
{
    if (!g) return NAN;
    double x, y;
    for (int i = 0; i < g->GetN(); ++i) {
        g->GetPoint(i, x, y);
        if (std::fabs(x - E) < 1.0) return y;
    }
    return NAN;
}

double errAtE(TGraphErrors* g, double E)
{
    if (!g) return NAN;
    double x, y;
    for (int i = 0; i < g->GetN(); ++i) {
        g->GetPoint(i, x, y);
        if (std::fabs(x - E) < 1.0) return g->GetErrorY(i);
    }
    return NAN;
}

// Fill a per-energy array from a graph <gname> in Summary/<file>.
std::vector<double> arrFromGraph(const char* file, const char* gname,
                                 bool errors = false)
{
    std::vector<double> v(kNE, NAN);
    TFile* f = TFile::Open(Form("%s%s", kSumDir, file));
    if (!f || f->IsZombie()) { if (f) delete f; return v; }
    TGraph* g = dynamic_cast<TGraph*>(f->Get(gname));
    for (int i = 0; i < kNE; ++i)
        v[i] = errors ? errAtE(dynamic_cast<TGraphErrors*>(g), kE[i])
                      : atE(g, kE[i]);
    delete f;
    return v;
}

double scalarFromFile(const char* file, const char* name)
{
    TFile* f = TFile::Open(Form("%s%s", kSumDir, file));
    if (!f || f->IsZombie()) { if (f) delete f; return NAN; }
    TParameter<double>* p = dynamic_cast<TParameter<double>*>(f->Get(name));
    double v = p ? p->GetVal() : NAN;
    delete f;
    return v;
}

// Fit sigma_t = sqrt(a^2/E + b^2) to a graph; return a, b by reference.
void fitAB(const char* file, const char* gname, double& a, double& b)
{
    a = NAN; b = NAN;
    TFile* f = TFile::Open(Form("%s%s", kSumDir, file));
    if (!f || f->IsZombie()) { if (f) delete f; return; }
    TGraphErrors* g = dynamic_cast<TGraphErrors*>(f->Get(gname));
    if (g && g->GetN() >= 2) {
        TF1 fit("fitAB", "sqrt([0]*[0]/x + [1]*[1])", 20., 160.);
        fit.SetParameters(220., 30.);
        g->Fit(&fit, "QN");
        a = fit.GetParameter(0);
        b = fit.GetParameter(1);
    }
    delete f;
}

// High-energy timing FLOOR.  Two parametrisations of the sigma_t-vs-E curve:
//   paper  : sigma_t = sqrt(a^2/E + b^2)              -> floor = b   (matches arXiv:2401.01747)
//   timing : sigma_t = sqrt((a/E)^2 + (b/sqrt(E))^2 + c^2) -> floor = c
// The 1/E term is the physically correct scaling for the electronics-noise/slew
// contribution to TIMING resolution (slew rate dV/dt ∝ amplitude ∝ E), which the
// paper's pure 1/sqrt(E) form omits.  We report both so the extrapolated floor is
// not hostage to one parametrisation.  Returns floor value + its fit error.
void fitFloorPaper(const char* file, const char* gname, double& b, double& bErr)
{
    b = NAN; bErr = NAN;
    TFile* f = TFile::Open(Form("%s%s", kSumDir, file));
    if (!f || f->IsZombie()) { if (f) delete f; return; }
    TGraphErrors* g = dynamic_cast<TGraphErrors*>(f->Get(gname));
    if (g && g->GetN() >= 2) {
        TF1 fit("fitFP", "sqrt([0]*[0]/x + [1]*[1])", 20., 160.);
        fit.SetParameters(220., 25.);
        g->Fit(&fit, "QN");
        b    = std::fabs(fit.GetParameter(1));
        bErr = fit.GetParError(1);
    }
    delete f;
}

void fitFloorTiming(const char* file, const char* gname, double& c, double& cErr)
{
    c = NAN; cErr = NAN;
    TFile* f = TFile::Open(Form("%s%s", kSumDir, file));
    if (!f || f->IsZombie()) { if (f) delete f; return; }
    TGraphErrors* g = dynamic_cast<TGraphErrors*>(f->Get(gname));
    if (g && g->GetN() >= 3) {   // 3 free parameters
        TF1 fit("fitFT", "sqrt([0]*[0]/(x*x) + [1]*[1]/x + [2]*[2])", 20., 160.);
        fit.SetParameters(1500., 200., 20.);
        fit.SetParLimits(2, 0., 100.);
        g->Fit(&fit, "QN");
        c    = std::fabs(fit.GetParameter(2));
        cErr = fit.GetParError(2);
    }
    delete f;
}

double nEventsAt(double E)
{
    TFile* f = TFile::Open(Form("Analysis/Output/%dGeV/ntuple.root", (int)E));
    if (!f || f->IsZombie()) { if (f) delete f; return NAN; }
    TTree* t = dynamic_cast<TTree*>(f->Get("rad"));
    double n = t ? static_cast<double>(t->GetEntries()) : NAN;
    delete f;
    return n;
}

double meanIgnoringNaN(const std::vector<double>& v)
{
    double s = 0.; int n = 0;
    for (double x : v) if (!std::isnan(x)) { s += x; ++n; }
    return n ? s / n : NAN;
}

// ── JSON formatting ────────────────────────────────────────────────────────
std::string num(double v)
{
    if (std::isnan(v)) return "null";
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.6g", v);
    return buf;
}
std::string arr(const std::vector<double>& v)
{
    std::string s = "[";
    for (size_t i = 0; i < v.size(); ++i) {
        s += num(v[i]);
        if (i + 1 < v.size()) s += ", ";
    }
    return s + "]";
}

} // namespace

void harvestResults()
{
    // ── Physics extraction (headline) ───────────────────────────────────────
    std::vector<double> teb_sigma     = arrFromGraph("timing_energy_bins.root", "gBestSigma_teb_m0");
    std::vector<double> teb_sigma_err = arrFromGraph("timing_energy_bins.root", "gBestSigma_teb_m0", true);
    std::vector<double> paper_sigma   = arrFromGraph("timing_energy_bins.root", "gPaper_teb");
    // Best-bin disclosure: efficiency (% of fiducial) and E_meas (mV) of the single
    // best E_meas bin whose sigma is the advertised teb_sigma headline.
    std::vector<double> teb_eff       = arrFromGraph("timing_energy_bins.root", "gBestEff_teb");
    std::vector<double> teb_ebin_mV   = arrFromGraph("timing_energy_bins.root", "gBestEmeas_teb");
    double teb_a = NAN, teb_b = NAN;
    fitAB("timing_energy_bins.root", "gBestSigma_teb_m0", teb_a, teb_b);

    // #G5 bias-corrected (out-of-sample, run-folded selection CV) headline curve.
    std::vector<double> teb_sigma_oos     = arrFromGraph("timing_energy_bins.root", "gBestSigmaOOS_teb_m0");
    std::vector<double> teb_sigma_oos_err = arrFromGraph("timing_energy_bins.root", "gBestSigmaOOS_teb_m0", true);
    // High-energy floor, fit on the OOS curve, in BOTH parametrisations.
    double teb_floor_paper = NAN,  teb_floor_paper_err = NAN;
    double teb_floor_timing = NAN, teb_floor_timing_err = NAN;
    fitFloorPaper ("timing_energy_bins.root", "gBestSigmaOOS_teb_m0", teb_floor_paper,  teb_floor_paper_err);
    fitFloorTiming("timing_energy_bins.root", "gBestSigmaOOS_teb_m0", teb_floor_timing, teb_floor_timing_err);
    // Lowest MEASURED point (no extrapolation): the highest-energy OOS sigma.
    double teb_low_meas = NAN;
    for (double s : teb_sigma_oos) if (!std::isnan(s) && (std::isnan(teb_low_meas) || s < teb_low_meas)) teb_low_meas = s;

    // ── Channel-combination timing ──────────────────────────────────────────
    std::vector<double> combo_a2_8ch  = arrFromGraph("timing_summary.root", "gTiming_M3"); // A^2-wgt all-8 (CFD-5%)
    std::vector<double> scan_all8     = arrFromGraph("channel_combination_scan.root", "gAll8");
    std::vector<double> scan_best7    = arrFromGraph("channel_combination_scan.root", "gBest7");
    std::vector<double> scan_noSWU    = arrFromGraph("channel_combination_scan.root", "gNoSWU");
    std::vector<double> scan_bestN4   = arrFromGraph("channel_combination_scan.root", "gBestN4");

    // ── Reference / DRS4 ────────────────────────────────────────────────────
    std::vector<double> mcp_jitter    = arrFromGraph("mcp_jitter.root", "gMCP_jitter");
    double mcp_jitter_mean            = meanIgnoringNaN(mcp_jitter);
    double drs4_cellwidth_rms = scalarFromFile("drs4_timebase.root", "cellWidthRMS_D0G0_ps");
    double drs4_combo_before  = scalarFromFile("drs4_timebase.root", "sigma_combo_before_ps");
    double drs4_combo_after   = scalarFromFile("drs4_timebase.root", "sigma_combo_after_ps");

    // ── Energy resolution + systematics ─────────────────────────────────────
    std::vector<double> eres          = arrFromGraph("summary.root", "gEnergyResolution");
    std::vector<double> eres_detector = arrFromGraph("summary.root", "gEnergyResolution_detector");
    std::vector<double> syst_total    = arrFromGraph("systematic_uncertainties.root", "gSystTotal");

    // ── Beam characterisation (newly persisted by the gap macros) ───────────
    std::vector<double> punch_through = arrFromGraph("pbglass_investigation.root", "gPunchThroughFrac");
    std::vector<double> containment   = arrFromGraph("cross_energy.root", "gContainedFrac");

    // ── Dataset size ────────────────────────────────────────────────────────
    std::vector<double> n_events(kNE, NAN);
    for (int i = 0; i < kNE; ++i) n_events[i] = nEventsAt(kE[i]);

    // ── Write results.json ──────────────────────────────────────────────────
    std::ofstream o(Form("%sresults.json", kSumDir));
    o << "{\n";
    o << "  \"_note\": \"Auto-generated by harvestResults.C from Output/Summary/*.root and per-energy ntuples. Do not edit by hand.\",\n";
    o << "  \"_estimator\": \"headline sigma_t = (DW-UP)/2, CFD-5%, energy-binned\",\n";
    o << "  \"energies\":        " << arr(std::vector<double>(kE, kE + kNE)) << ",\n";
    o << "  \"teb_sigma\":       " << arr(teb_sigma)     << ",\n";
    o << "  \"teb_sigma_err\":   " << arr(teb_sigma_err) << ",\n";
    o << "  \"teb_fit_a\":       " << num(teb_a)         << ",\n";
    o << "  \"teb_fit_b\":       " << num(teb_b)         << ",\n";
    o << "  \"paper_sigma\":     " << arr(paper_sigma)   << ",\n";
    o << "  \"teb_eff\":         " << arr(teb_eff)       << ",\n";
    o << "  \"teb_ebin_mV\":     " << arr(teb_ebin_mV)   << ",\n";
    o << "  \"teb_sigma_oos\":     " << arr(teb_sigma_oos)     << ",\n";
    o << "  \"teb_sigma_oos_err\": " << arr(teb_sigma_oos_err) << ",\n";
    o << "  \"teb_floor_paper\":      " << num(teb_floor_paper)      << ",\n";
    o << "  \"teb_floor_paper_err\":  " << num(teb_floor_paper_err)  << ",\n";
    o << "  \"teb_floor_timing\":     " << num(teb_floor_timing)     << ",\n";
    o << "  \"teb_floor_timing_err\": " << num(teb_floor_timing_err) << ",\n";
    o << "  \"teb_low_meas\":         " << num(teb_low_meas)         << ",\n";
    o << "  \"combo_a2_8ch\":    " << arr(combo_a2_8ch)  << ",\n";
    o << "  \"scan_all8\":       " << arr(scan_all8)     << ",\n";
    o << "  \"scan_best7\":      " << arr(scan_best7)    << ",\n";
    o << "  \"scan_noSWU\":      " << arr(scan_noSWU)    << ",\n";
    o << "  \"scan_bestN4\":     " << arr(scan_bestN4)   << ",\n";
    o << "  \"mcp_jitter\":      " << arr(mcp_jitter)    << ",\n";
    o << "  \"mcp_jitter_mean\": " << num(mcp_jitter_mean) << ",\n";
    o << "  \"drs4_cellwidth_rms\": " << num(drs4_cellwidth_rms) << ",\n";
    o << "  \"drs4_combo_before\":  " << num(drs4_combo_before)  << ",\n";
    o << "  \"drs4_combo_after\":   " << num(drs4_combo_after)   << ",\n";
    o << "  \"eres\":            " << arr(eres)          << ",\n";
    o << "  \"eres_detector\":   " << arr(eres_detector) << ",\n";
    o << "  \"syst_total\":      " << arr(syst_total)    << ",\n";
    o << "  \"punch_through\":   " << arr(punch_through) << ",\n";
    o << "  \"containment\":     " << arr(containment)   << ",\n";
    o << "  \"n_events\":        " << arr(n_events)      << "\n";
    o << "}\n";
    o.close();

    std::printf("harvestResults: wrote %sresults.json\n", kSumDir);
    std::printf("  teb_sigma@150 = %.1f ps,  combo(M3)@150 = %.1f ps,  "
                "mcp_mean = %.1f ps\n",
                teb_sigma[5], combo_a2_8ch[5], mcp_jitter_mean);
    std::printf("  punch_through@150 = %.1f%%,  containment@150 = %.1f%%,  "
                "n_events@150 = %.0f\n",
                punch_through[5], containment[5], n_events[5]);
}

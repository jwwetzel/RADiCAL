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
// Run:  root -l -b -q 'analyze/studies/harvestResults.C+'   (after all analysis
// macros; the headline producer is analyze/studies/timingProduction.C — the
// PRODUCTION gated chain: brightest-1000 (DW-UP)/2, srCFD/LED per regime).
// MIGRATED 2026-07-21 from the retired best-bin/cfd05 chain.
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
const char*  kSumDir = "output/Summary/";

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

// Fit sigma_t = sqrt(a^2/E + b^2) to a graph; return a, b (+errors) by reference.
// Gated convention: parameter errors inflated by sqrt(chi2/ndf) when > 1 (PDG).
void fitAB(const char* file, const char* gname, double& a, double& b,
           double* aErr = nullptr, double* bErr = nullptr)
{
    a = NAN; b = NAN; if (aErr) *aErr = NAN; if (bErr) *bErr = NAN;
    TFile* f = TFile::Open(Form("%s%s", kSumDir, file));
    if (!f || f->IsZombie()) { if (f) delete f; return; }
    TGraphErrors* g = dynamic_cast<TGraphErrors*>(f->Get(gname));
    if (g && g->GetN() >= 2) {
        TF1 fit("fitAB", "sqrt([0]*[0]/x + [1]*[1])", 20., 160.);
        fit.SetParameters(220., 30.);
        g->Fit(&fit, "QN");
        a = std::fabs(fit.GetParameter(0));
        b = std::fabs(fit.GetParameter(1));
        double ae = fit.GetParError(0), be = fit.GetParError(1);
        if (fit.GetNDF() > 0 && fit.GetChisquare()/fit.GetNDF() > 1) {
            double k = std::sqrt(fit.GetChisquare()/fit.GetNDF()); ae *= k; be *= k; }
        if (aErr) *aErr = ae; if (bErr) *bErr = be;
    }
    delete f;
}

// Read a TParameter<double> / TNamed string from Summary/timing_energy_bins.root
double tebPar(const char* name)
{
    return scalarFromFile("timing_energy_bins.root", name);
}
std::string tebStr(const char* name)
{
    TFile* f = TFile::Open(Form("%stiming_energy_bins.root", kSumDir));
    if (!f || f->IsZombie()) { if (f) delete f; return ""; }
    TNamed* n = dynamic_cast<TNamed*>(f->Get(name));
    std::string s = n ? n->GetTitle() : "";
    delete f;
    return s;
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
    TFile* f = TFile::Open(Form("output/%dGeV/ntuple.root", (int)E));
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
    double teb_a = NAN, teb_b = NAN, teb_a_err = NAN, teb_b_err = NAN;
    fitAB("timing_energy_bins.root", "gBestSigma_teb_m0", teb_a, teb_b, &teb_a_err, &teb_b_err);

    // Full-fiducial companion (claims law: always quoted beside the brightest-1000
    // headline). Protocol: r = 3.0 mm at all energies, matching the gated
    // full_fiducial_check log (50.5 ps at 150 GeV).
    std::vector<double> sigma_ff      = arrFromGraph("timing_energy_bins.root", "gSigmaFF_teb");
    std::vector<double> sigma_ff_err  = arrFromGraph("timing_energy_bins.root", "gSigmaFF_teb", true);
    std::vector<double> sigma_cfd05   = arrFromGraph("timing_energy_bins.root", "gSigmaCFD05_teb");

    // Per-build production block (adopted per-regime sources), written by
    // timingProduction.C from live timingBrightestK passes — reproduces the gated
    // table papers/tables/timing_fit_summary_2026-06-09.md exactly (verified in
    // its stdout GATED CHECK).
    const char* kBuilds[4] = {"DSB1","LUAG","MIXED","TENERGY"};
    double bl_a[4], bl_ae[4], bl_b[4], bl_be[4], bl_s[4], bl_se[4], bl_l[4];
    std::string bl_src[4];
    for (int i = 0; i < 4; ++i) {
        bl_a[i]  = tebPar(Form("fit_a_%s",    kBuilds[i]));
        bl_ae[i] = tebPar(Form("fit_aerr_%s", kBuilds[i]));
        bl_b[i]  = tebPar(Form("fit_b_%s",    kBuilds[i]));
        bl_be[i] = tebPar(Form("fit_berr_%s", kBuilds[i]));
        bl_s[i]  = tebPar(Form("s150_%s",     kBuilds[i]));
        bl_se[i] = tebPar(Form("s150err_%s",  kBuilds[i]));
        bl_l[i]  = tebPar(Form("light150_%s", kBuilds[i]));
        bl_src[i]= tebStr(Form("src_%s",      kBuilds[i]));
    }

    // Gated companion results imported by value WITH provenance (committed gate
    // logs; not recomputed here): the MIXED same-shower ratio, the depth-dial
    // slope, and the per-build selection systematics.
    //   mixed ratio : papers/scripts/mixed_killshot_bootstrap (GATE 6 log)
    //   depth slope : papers/scripts/depth_dial (GATE 3 log)
    //   syst sel    : papers/tables/systematics_postfix_2026-06-09.md
    const double kMixedRatio = 1.04,  kMixedRatioErr = 0.05;
    const double kDepthSlope = -33.6, kDepthSlopeErr = 2.9;
    const double kSystSel[4] = {1.0, 1.9, 0.9, 1.1};   // DSB1, LUAG, MIXED, TENERGY

    // In the production chain the selection is deterministic (brightest-1000),
    // so the OOS graph is identical to the nominal one by construction — kept
    // under its legacy name for key compatibility.
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
    // Wire-chamber spatial resolution (data-driven self-calibration; reads raw
    // waveforms, so optional — null on a run without the raw files).
    double wc_res_x_mm   = scalarFromFile("wire_chamber_resolution.root", "wc_res_x_ub_um");
    double wc_res_y_mm   = scalarFromFile("wire_chamber_resolution.root", "wc_res_y_ub_um");
    double wc_res_cfd_mm = scalarFromFile("wire_chamber_resolution.root", "wc_res_x_ub_cfd_um");
    double wc_beam_sx    = scalarFromFile("wire_chamber_resolution.root", "wc_beam_sigma_x_mm");
    double wc_beam_sy    = scalarFromFile("wire_chamber_resolution.root", "wc_beam_sigma_y_mm");
    if (!std::isnan(wc_res_x_mm))   wc_res_x_mm   /= 1000.;   // um -> mm
    if (!std::isnan(wc_res_y_mm))   wc_res_y_mm   /= 1000.;
    if (!std::isnan(wc_res_cfd_mm)) wc_res_cfd_mm /= 1000.;

    // Shashlik module centre from the calorimeter edges (moduleCenter.C; optional).
    double mod_center_x = scalarFromFile("module_center.root", "module_center_x");
    double mod_center_y = scalarFromFile("module_center.root", "module_center_y");
    double mod_width_x  = scalarFromFile("module_center.root", "module_width_x");
    double mod_width_y  = scalarFromFile("module_center.root", "module_width_y");

    // Transverse alignment (alignmentAnalysis.C; optional).
    double mcp_center_x  = scalarFromFile("alignment.root", "mcp_center_x");
    double mcp_center_y  = scalarFromFile("alignment.root", "mcp_center_y");
    double beam_center_x = scalarFromFile("alignment.root", "beam_center_x");
    double beam_center_y = scalarFromFile("alignment.root", "beam_center_y");
    double off_mcp_rad   = scalarFromFile("alignment.root", "off_mcp_rad");
    double off_beam_rad  = scalarFromFile("alignment.root", "off_beam_rad");

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
    o << "  \"_note\": \"Auto-generated by harvestResults.C from output/Summary/*.root and the reduced ntuples. Do not edit by hand.\",\n";
    o << "  \"_estimator\": \"PRODUCTION chain: brightest-1000 (DW-UP)/2 with in-event consistency veto and robust core width (tebSigma); srCFD (saturation-recovered CFD, branch hg_lgcfd) for the high-light builds, LED for the dim builds (per-regime rule). Full-fiducial companion quoted alongside (r=3.0 mm protocol).\",\n";
    o << "  \"_provenance\": \"headline graphs: timingProduction.C (live timingBrightestK, verified vs papers/tables/timing_fit_summary_2026-06-09.md); mixed ratio / depth slope / selection systematics imported from the committed gate logs under papers/scripts/.\",\n";
    o << "  \"energies\":        " << arr(std::vector<double>(kE, kE + kNE)) << ",\n";
    o << "  \"teb_sigma\":       " << arr(teb_sigma)     << ",\n";
    o << "  \"teb_sigma_err\":   " << arr(teb_sigma_err) << ",\n";
    o << "  \"teb_fit_a\":       " << num(teb_a)         << ",\n";
    o << "  \"teb_fit_a_err\":   " << num(teb_a_err)     << ",\n";
    o << "  \"teb_fit_b\":       " << num(teb_b)         << ",\n";
    o << "  \"teb_fit_b_err\":   " << num(teb_b_err)     << ",\n";
    o << "  \"dsb1_sigma_ff\":       " << arr(sigma_ff)     << ",\n";
    o << "  \"dsb1_sigma_ff_err\":   " << arr(sigma_ff_err) << ",\n";
    o << "  \"dsb1_sigma_cfd05\":    " << arr(sigma_cfd05)  << ",\n";
    o << "  \"builds\": {\n";
    for (int i = 0; i < 4; ++i) {
        o << "    \"" << kBuilds[i] << "\": {"
          << "\"fit_a\": "     << num(bl_a[i])  << ", \"fit_a_err\": " << num(bl_ae[i])
          << ", \"fit_b\": "   << num(bl_b[i])  << ", \"fit_b_err\": " << num(bl_be[i])
          << ", \"s150\": "    << num(bl_s[i])  << ", \"s150_err\": "  << num(bl_se[i])
          << ", \"light150_mV\": " << num(bl_l[i])
          << ", \"source\": \"" << bl_src[i] << "\""
          << ", \"syst_sel\": " << num(kSystSel[i]) << "}" << (i < 3 ? "," : "") << "\n";
    }
    o << "  },\n";
    o << "  \"mixed_sameshower_ratio\":     " << num(kMixedRatio)    << ",\n";
    o << "  \"mixed_sameshower_ratio_err\": " << num(kMixedRatioErr) << ",\n";
    o << "  \"depth_slope_ps_per_efold\":   " << num(kDepthSlope)    << ",\n";
    o << "  \"depth_slope_err\":            " << num(kDepthSlopeErr) << ",\n";
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
    o << "  \"wc_res_x_mm\":     " << num(wc_res_x_mm)   << ",\n";
    o << "  \"wc_res_y_mm\":     " << num(wc_res_y_mm)   << ",\n";
    o << "  \"wc_res_cfd_mm\":   " << num(wc_res_cfd_mm) << ",\n";
    o << "  \"wc_beam_sx_mm\":   " << num(wc_beam_sx)    << ",\n";
    o << "  \"wc_beam_sy_mm\":   " << num(wc_beam_sy)    << ",\n";
    o << "  \"mod_center_x\":    " << num(mod_center_x) << ",\n";
    o << "  \"mod_center_y\":    " << num(mod_center_y) << ",\n";
    o << "  \"mod_width_x\":     " << num(mod_width_x)  << ",\n";
    o << "  \"mod_width_y\":     " << num(mod_width_y)  << ",\n";
    o << "  \"mcp_center_x\":    " << num(mcp_center_x)  << ",\n";
    o << "  \"mcp_center_y\":    " << num(mcp_center_y)  << ",\n";
    o << "  \"beam_center_x\":   " << num(beam_center_x) << ",\n";
    o << "  \"beam_center_y\":   " << num(beam_center_y) << ",\n";
    o << "  \"off_mcp_rad\":     " << num(off_mcp_rad)   << ",\n";
    o << "  \"off_beam_rad\":    " << num(off_beam_rad)  << ",\n";
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
    std::printf("  PRODUCTION: sigma@150 = %.1f ps (brightest-1000, srCFD), "
                "full-fiducial = %.1f ps, fit a=%.0f+-%.0f b=%.1f+-%.1f\n",
                teb_sigma.size() > 5 ? teb_sigma[5] : NAN,
                sigma_ff.size() > 5 ? sigma_ff[5] : NAN,
                teb_a, teb_a_err, teb_b, teb_b_err);
    std::printf("  teb_sigma@150 = %.1f ps,  combo(M3)@150 = %.1f ps,  "
                "mcp_mean = %.1f ps\n",
                teb_sigma[5], combo_a2_8ch[5], mcp_jitter_mean);
    std::printf("  punch_through@150 = %.1f%%,  containment@150 = %.1f%%,  "
                "n_events@150 = %.0f\n",
                punch_through[5], containment[5], n_events[5]);
}

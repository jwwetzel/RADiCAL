// ============================================================================
// qualityPlots.C — data quality, channel performance, and cut optimisation
// ============================================================================
//
// Reads Analysis/Output/<label>/ntuple.root and produces diagnostic plots to:
//   (a) confirm the detector and beam are healthy (data quality)
//   (b) characterise each SiPM / capillary channel (channel performance)
//   (c) show how the timing resolution and event efficiency depend on each
//       cut threshold (value-vs-cut optimisation)
//
// Run after processRun.C has been executed for all energies.
//
// ── Per-energy output (Analysis/Output/<label>/quality_report.pdf) ──────────
//
//   Page 1  Beam quality
//             2D WC hit map (both fiducial circles) · x/y beam profiles
//             MCP1 and MCP2 amplitude distributions with cut lines
//             Event statistics table
//
//   Page 2  HG channel amplitudes  (8 SiPM/capillary channels)
//             Peak amplitude distribution + hit efficiency table
//
//   Page 3  LG channel amplitudes + HG/LG amplitude correlation
//
//   Page 4  Time walk diagnostic
//             TProfile of t_CFD vs 1/A_HG (physics model: t = a + b/A is
//             linear in 1/A).  Slope b = M6 correction coefficient.
//
//   Page 5  TOT & LED diagnostics  (only if hg_led branch present)
//             TOT distribution + TProfile(t_LED vs TOT) showing M5 correction
//
//   Page 6  Cut optimisation — MCP amplitude threshold
//             Timing resolution and event efficiency vs MCP1 amplitude cut
//             Guides choice of kMCP1_minPeak (default 100 mV)
//
//   Page 7  Cut optimisation — fiducial radius
//             Timing resolution and event efficiency vs beam radius cut
//             Guides choice of kFiducial_r_timing (default 3.0 mm)
//
//   Page 8  Cut optimisation — HG amplitude threshold
//             Timing resolution and event efficiency vs minimum HG amplitude
//             Guides choice of kHG_minPeak (default 20 mV)
//
//   Page 9  Fiducial selection overview
//             WC hit map at each cut stage · radial beam profile · cut flow table
//             sum_LG before/after containment cut · map of rejected events
//
//   Page 10 PbGlass shower containment
//             sum_PbGlass vs sum_LG scatter (2D) with kPb_maxRatio cut line
//             PbGlass/RADiCAL ratio distribution · containment fraction vs radius
//             Guides choice of kPb_maxRatio (default 0.30)
//
// ── Summary output ────────────────────────────────────────────────────────────
//   Analysis/Output/Summary/quality_summary.pdf — efficiency vs cut, all energies
//
// Usage (from repository root):
//   root -l -b -q 'Analysis/qualityPlots.C+'
// ============================================================================

#include "ChannelConfig.h"
#include "PlotUtils.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TArc.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TAxis.h"
#include "TBox.h"

#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>

// ---------------------------------------------------------------------------
// Event data stored in memory for cut optimisation scans
// (only wc_ok + loose MCP events; fiducial cut applied during scan)
// ---------------------------------------------------------------------------
struct EvtData {
    float r_beam;     // distance from data-derived centroid [mm]
    float mcp_peak;   // MCP1 amplitude [mV]
    float mcp2_peak;  // MCP2 amplitude [mV]
    float hg_peak[8]; // HG amplitudes [mV]
    float hg_cfd[8];  // HG timing [ns] — CFD-5% (adopted headline fraction) when
                      // the hg_cfd05 branch is present, else CFD-20% fallback.
};

// ---------------------------------------------------------------------------
// Helper: fit Gaussian to a vector of timing values, return sigma [ns]
// Optionally restrict to events where r_beam < r_cut and amp[ch] > a_cut
// ---------------------------------------------------------------------------
static double TimingSigma(const std::vector<EvtData>& evts,
                           int ch, float r_cut, float a_cut)
{
    // Collect timing values passing both cuts
    std::vector<float> tv;
    tv.reserve(evts.size());
    for (size_t i = 0; i < evts.size(); ++i) {
        const EvtData& e = evts[i];
        if (e.r_beam   >= r_cut) continue;
        if (e.hg_peak[ch] < a_cut) continue;
        if (e.hg_cfd[ch]  < -1e5f) continue;
        tv.push_back(e.hg_cfd[ch]);
    }
    if ((int)tv.size() < 50) return -1.;

    // Auto-range histogram
    double mu = 0.;
    for (size_t i = 0; i < tv.size(); ++i) mu += tv[i];
    mu /= tv.size();
    double rms = 0.;
    for (size_t i = 0; i < tv.size(); ++i) rms += (tv[i]-mu)*(tv[i]-mu);
    rms = std::sqrt(rms / tv.size());
    if (rms < 0.010) rms = 0.200;

    TH1F htmp("_qp_ts", "", 200, mu - 4.*rms, mu + 4.*rms);
    for (size_t i = 0; i < tv.size(); ++i) htmp.Fill(tv[i]);
    double m2, mE, sig, sigE;
    FitGaussCore(&htmp, 2.0, m2, mE, sig, sigE);
    return sig;
}

// ---------------------------------------------------------------------------
// Helper: count events passing r_beam < r_cut and hg_peak[ch] >= a_cut
// ---------------------------------------------------------------------------
static int CountEvents(const std::vector<EvtData>& evts,
                        int ch, float r_cut, float a_cut)
{
    int n = 0;
    for (size_t i = 0; i < evts.size(); ++i) {
        if (evts[i].r_beam      >= r_cut) continue;
        if (evts[i].hg_peak[ch] <  a_cut) continue;
        if (evts[i].hg_cfd[ch]  < -1e5f)  continue;
        ++n;
    }
    return n;
}

// ---------------------------------------------------------------------------
// Helper: fit Gaussian to timing values with MCP amplitude threshold scan.
// Events must also pass r_beam < r_cut and hg_peak[ch] >= a_cut.
// ---------------------------------------------------------------------------
static double TimingSigmaMCP(const std::vector<EvtData>& evts,
                              int ch, float r_cut, float mcp_cut, float a_cut)
{
    std::vector<float> tv;
    tv.reserve(evts.size());
    for (size_t i = 0; i < evts.size(); ++i) {
        const EvtData& e = evts[i];
        if (e.r_beam      >= r_cut)  continue;
        if (e.mcp_peak    <  mcp_cut) continue;
        if (e.hg_peak[ch] <  a_cut)  continue;
        if (e.hg_cfd[ch]  < -1e5f)   continue;
        tv.push_back(e.hg_cfd[ch]);
    }
    if ((int)tv.size() < 50) return -1.;
    double mu = 0.;
    for (size_t i = 0; i < tv.size(); ++i) mu += tv[i];
    mu /= tv.size();
    double rms = 0.;
    for (size_t i = 0; i < tv.size(); ++i) rms += (tv[i]-mu)*(tv[i]-mu);
    rms = std::sqrt(rms / tv.size());
    if (rms < 0.010) rms = 0.200;
    TH1F htmp("_qp_mcps", "", 200, mu - 4.*rms, mu + 4.*rms);
    for (size_t i = 0; i < tv.size(); ++i) htmp.Fill(tv[i]);
    double m2, mE, sig, sigE;
    FitGaussCore(&htmp, 2.0, m2, mE, sig, sigE);
    return sig;
}

static int CountEventsMCP(const std::vector<EvtData>& evts,
                           int ch, float r_cut, float mcp_cut, float a_cut)
{
    int n = 0;
    for (size_t i = 0; i < evts.size(); ++i) {
        if (evts[i].r_beam      >= r_cut)  continue;
        if (evts[i].mcp_peak    <  mcp_cut) continue;
        if (evts[i].hg_peak[ch] <  a_cut)  continue;
        if (evts[i].hg_cfd[ch]  < -1e5f)   continue;
        ++n;
    }
    return n;
}

// ---------------------------------------------------------------------------
// Channel colors from RADiCALStyle.h (kRChannelCols, set by ApplyRADiCALStyle)

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------
void qualityPlots()
{
    // Prevent histograms from auto-registering to the current TFile directory.
    // Without this, all TH1F / TProfile objects booked while an input file is
    // the current directory are owned by that file and destroyed when
    // fin->Close() is called — before the drawing code runs.
    TH1::AddDirectory(kFALSE);

    ApplyRADiCALStyle();   // RADiCALStyle.h via PlotUtils.h

    gSystem->mkdir("Analysis/Output/Summary", kTRUE);

    // Per-energy efficiency at the nominal timing cuts — for summary page
    // vEff[ie][iCut][ch] where iCut: 0=r<3mm, 1=A>20mV
    std::vector<double> vEnergy;

    // ── Channel validation metrics: filled per-energy, drawn in summary PDF ──
    // Indexed [iRun][iCh].  iRun increments for each successfully processed energy.
    double sumHitEff [kNRuns][8] = {};  // hit efficiency [%]     (higher = better)
    double sumPedRMS [kNRuns][8] = {};  // mean HG pedestal noise [mV] (lower = better)
    double sumMeanHG [kNRuns][8] = {};  // mean HG peak amplitude [mV] (signal events)
    double sumTimeSig[kNRuns][8] = {};  // CFD-20% timing resolution [ps] (lower = better)
    double sumContain[kNRuns]    = {};  // PbGlass containment efficiency [%]
    int    nRunsDone = 0;

    // =========================================================================
    // Per-energy loop
    // =========================================================================
    for (int iRun = 0; iRun < kNRuns; ++iRun) {
        const RunCfg& rc = kRuns[iRun];
        TString ntFile   = TString("Analysis/Output/") + rc.label + "/ntuple.root";

        TFile* fin = TFile::Open(ntFile);
        if (!fin || fin->IsZombie()) {
            std::cout << "[qualityPlots] Skipping " << rc.label
                      << " — ntuple not found.\n";
            continue;
        }
        TTree* tree = (TTree*)fin->Get("rad");
        if (!tree || tree->GetEntries() == 0) { fin->Close(); continue; }

        bool hasNewBranches = (tree->GetBranch("hg_led") != nullptr);

        std::cout << "[qualityPlots] " << rc.label
                  << " — " << tree->GetEntries() << " events\n";

        // Derive per-run beam centroid
        double xc, yc, tOff, tRms;
        ScanRunCenters(tree, xc, yc, tOff, tRms);
        const float xcf = static_cast<float>(xc);
        const float ycf = static_cast<float>(yc);

        // Timing histogram window (centred on measured offset)
        double tWin = std::max(5.0 * tRms, 2.0);
        double tLo  = tOff - tWin;
        double tHi  = tOff + tWin;

        // ---------------------------------------------------------------------
        // Book histograms
        // ---------------------------------------------------------------------
        // Beam quality — bins matched to WC delay-line resolution (kWC_resBin = 1 mm)
        double hmHW   = 20.;
        // Occupancy maps: bin at 2 mm, an integer multiple (5x) of the ~0.4 mm WC
        // delay-line TDC quantum.  Binning at 1 mm (kWC_resBin) aliases that quantum
        // into a moire "comb" (the visible striping); 2 mm cancels it (bin-to-bin
        // comb roughness 10%->6%, measured on data) and still oversamples the
        // ~3.6 mm WC position resolution.
        const double kHitMapBin = 2.0;  // mm/bin for WC occupancy maps
        int    nBins2D = static_cast<int>(std::round(2.*hmHW / kHitMapBin)); // 20 bins, 2 mm
        TH2F* h2D = new TH2F(Form("h2D_%s", rc.label.Data()),
            ";x_{WC} (mm);y_{WC} (mm)",
            nBins2D, xc-hmHW, xc+hmHW, nBins2D, yc-hmHW, yc+hmHW);
        TH1F* hX  = new TH1F(Form("hX_%s",  rc.label.Data()),
            ";x_{WC} (mm);Events", nBins2D, xc-hmHW, xc+hmHW);
        TH1F* hY  = new TH1F(Form("hY_%s",  rc.label.Data()),
            ";y_{WC} (mm);Events", nBins2D, yc-hmHW, yc+hmHW);
        TH1F* hMCP  = new TH1F(Form("hMCP_%s",  rc.label.Data()),
            ";MCP1 amplitude (mV);Events", 150, 0, 1500);
        TH1F* hMCP2 = new TH1F(Form("hMCP2_%s", rc.label.Data()),
            ";MCP2 amplitude (mV);Events", 150, 0, 1500);

        // HG and LG amplitudes per channel (for fiducial events)
        TH1F* hHG[8], *hLG[8];
        for (int i = 0; i < 8; ++i) {
            hHG[i] = new TH1F(Form("hHG_%s_%d", rc.label.Data(), i),
                Form(";%s HG amplitude (mV);Events", kCap[i].name),
                150, 0, 500);
            hLG[i] = new TH1F(Form("hLG_%s_%d", rc.label.Data(), i),
                Form(";%s LG amplitude (mV);Events", kCap[i].name),
                150, 0, 3000);
        }

        // HG/LG amplitude correlation per channel (TProfile: mean HG vs LG bin)
        TProfile* pCorr[8];
        for (int i = 0; i < 8; ++i)
            pCorr[i] = new TProfile(Form("pCorr_%s_%d", rc.label.Data(), i),
                Form(";%s LG amplitude (mV);Mean HG amplitude (mV)", kCap[i].name),
                50, 0, 3000);

        // Time walk: mean t_cfd vs 1/hg_peak for TIMEABLE pulses.  First-order walk
        // gives t = a + b/A (linear on a 1/A axis).  Restricted to A > 100 mV
        // (1/A < 0.01): below that the CFD has no real leading edge and floors at
        // the search-window edge (~-90 ns) — mistiming, not walk.  See the fill cut.
        TProfile* pWalk[8];
        for (int i = 0; i < 8; ++i)
            pWalk[i] = new TProfile(Form("pWalk_%s_%d", rc.label.Data(), i),
                ";1/A_{HG} (1/mV);Mean t_{CFD} (ns)",
                40, 0.001, 0.011);

        // TOT distribution and LED vs TOT correlation (if new branches present)
        TH1F*    hTOT[8];
        TProfile* pLEDvsTOT[8];
        for (int i = 0; i < 8; ++i) {
            hTOT[i] = new TH1F(Form("hTOT_%s_%d", rc.label.Data(), i),
                Form(";%s TOT (ns);Events", kCap[i].name),
                100, 0, 30);
            pLEDvsTOT[i] = new TProfile(Form("pLedTot_%s_%d", rc.label.Data(), i),
                Form(";%s TOT (ns);Mean t_{LED} (ns)", kCap[i].name),
                50, 0, 30);
        }

        // Fiducial selection overview (Page 9)
        // hR: radial distance from centroid for all WC-ok events
        TH1F* hR = new TH1F(Form("hR_%s", rc.label.Data()),
            ";r_{beam} (mm);Events", 100, 0., 15.);
        // h2Dsel: hit map of events passing energy fiducial + containment
        TH2F* h2Dsel = new TH2F(Form("h2Dsel_%s", rc.label.Data()),
            ";x_{WC} (mm);y_{WC} (mm)",
            nBins2D, xc-hmHW, xc+hmHW, nBins2D, yc-hmHW, yc+hmHW);
        // h2DnoC: hit map of events inside energy fiducial but FAILING containment
        TH2F* h2DnoC = new TH2F(Form("h2DnoC_%s", rc.label.Data()),
            ";x_{WC} (mm);y_{WC} (mm)",
            nBins2D, xc-hmHW, xc+hmHW, nBins2D, yc-hmHW, yc+hmHW);
        // hSumLG_fid: sum_lg for events inside energy fiducial (before containment)
        TH1F* hSumLG_fid = new TH1F(Form("hSumLGfid_%s", rc.label.Data()),
            ";#SigmaLG_{RADiCAL} (mV);Events", 150, 0., 15000.);
        // hSumLG_cont: sum_lg for events also passing the containment cut
        TH1F* hSumLG_cont = new TH1F(Form("hSumLGcont_%s", rc.label.Data()),
            ";#SigmaLG_{RADiCAL} (mV);Events", 150, 0., 15000.);

        // PbGlass shower containment (Page 10)
        // Axis ranges are derived from beam energy: sum_lg peaks at ~30 mV/GeV,
        // so xmax = 55×E covers well beyond 3σ of the signal distribution at
        // every energy.  ymax = 0.50×xmax allows for extreme punch-through.
        const double xmax_lg = 55. * rc.energy_GeV;   // ~2× expected mean sum_lg
        const double ymax_pb = 0.50 * xmax_lg;        // PbGlass ≤ ~50% of LG range

        // hPbVsLG: 2D scatter of sum_pb vs sum_lg (shower containment diagnostic)
        TH2F* hPbVsLG = new TH2F(Form("hPbLG_%s", rc.label.Data()),
            ";#SigmaLG_{RADiCAL} (mV);#SigmaPbGlass (mV)",
            100, 0., xmax_lg, 100, 0., ymax_pb);
        // hPbRatio: distribution of sum_pb / sum_lg
        TH1F* hPbRatio = new TH1F(Form("hPbRat_%s", rc.label.Data()),
            ";#SigmaPbGlass / #SigmaLG;Events", 100, 0., 2.0);
        // pContVsR: mean containment fraction vs radial distance from centroid
        TProfile* pContVsR = new TProfile(Form("pContR_%s", rc.label.Data()),
            ";r_{beam} (mm);#SigmaPbGlass/(#SigmaPbGlass+#SigmaLG)",
            40, 0., 8.);

        // ALL-events PbGlass scatter (Page 10, Pad 1) — filled for every WC-ok
        // event with no fiducial or LG threshold.  Shows the full population
        // (beam halo at origin, full-energy showers in upper-right cluster).
        TH2F* hPbVsLG_all = new TH2F(Form("hPbLGall_%s", rc.label.Data()),
            ";#SigmaLG_{RADiCAL} (mV);#SigmaPbGlass (mV)",
            100, 0., xmax_lg, 100, 0., ymax_pb);
        TH1F* hSumLG_all  = new TH1F(Form("hSumLGall_%s", rc.label.Data()),
            ";#SigmaLG_{RADiCAL} (mV);Events", 150, 0., xmax_lg);

        // Pedestal RMS accumulator — averaged over WC-ok events.
        // Pedestal noise is a channel property (DRS4 capacitor cell variance),
        // not shower-dependent, so all WC-ok events contribute equally.
        double pedRMSsum[8]   = {};
        long   pedRMScount[8] = {};

        // ---------------------------------------------------------------------
        // Branch addresses
        // ---------------------------------------------------------------------
        Bool_t  wc_ok;
        Float_t x_trk, y_trk, mcp_peak, mcp2_peak, mcp_time, mcp2_time;
        Float_t hg_peak[8], hg_cfd[8], lg_peak[8];
        Float_t hg_led[8], hg_tot[8];

        tree->SetBranchAddress("wc_ok",    &wc_ok);
        tree->SetBranchAddress("x_trk",    &x_trk);
        tree->SetBranchAddress("y_trk",    &y_trk);
        tree->SetBranchAddress("mcp_peak", &mcp_peak);
        tree->SetBranchAddress("hg_peak",  hg_peak);
        // Per-channel timing uses CFD-5% (adopted headline fraction) when available:
        // it removes the broad CFD-20% shoulder on the Down capillaries (whose
        // leading-edge SHAPE jitters more high on the edge — not a mean-slope
        // effect).  Falls back to CFD-20% for pre-reprocess ntuples.  See
        // timingMethods.C page 3 and edgeMechanism.C.
        tree->SetBranchAddress(tree->GetBranch("hg_cfd05") ? "hg_cfd05" : "hg_cfd", hg_cfd);
        tree->SetBranchAddress("lg_peak",  lg_peak);

        // MCP2 and LED/TOT may be absent in older ntuples
        bool hasMCP2 = (tree->GetBranch("mcp2_peak") != nullptr);
        mcp2_peak = 0.f;
        mcp_time  = -1e6f;   // kNoTime sentinel (WaveformUtils.h not included here)
        mcp2_time = -1e6f;
        if (hasMCP2) {
            tree->SetBranchAddress("mcp2_peak", &mcp2_peak);
            // mcp_time and mcp2_time are present if mcp2_peak is
            if (tree->GetBranch("mcp_time"))  tree->SetBranchAddress("mcp_time",  &mcp_time);
            if (tree->GetBranch("mcp2_time")) tree->SetBranchAddress("mcp2_time", &mcp2_time);
        }

        if (hasNewBranches) {
            tree->SetBranchAddress("hg_led", hg_led);
            tree->SetBranchAddress("hg_tot", hg_tot);
        }

        // PbGlass — present in all ntuples produced by current processRun.C
        Float_t pb_peak_v[4], sum_pb_v = 0.f;
        bool hasPb = (tree->GetBranch("pb_peak") != nullptr);
        if (hasPb) {
            tree->SetBranchAddress("pb_peak", pb_peak_v);
            tree->SetBranchAddress("sum_pb",  &sum_pb_v);
        }

        // HG pedestal RMS per channel — added in last processRun.C update
        Float_t hg_ped_rms_v[8] = {};
        bool hasPedRMS = (tree->GetBranch("hg_ped_rms") != nullptr);
        if (hasPedRMS) tree->SetBranchAddress("hg_ped_rms", hg_ped_rms_v);

        // ---------------------------------------------------------------------
        // Fill histograms + load event cache for cut optimisation
        // ---------------------------------------------------------------------
        std::vector<EvtData> cache;
        cache.reserve(5000);

        long nWC = 0, nFid = 0, nTotal = 0, nContained = 0;
        long nMCP1sat = 0, nMCP2sat = 0;  // DRS4-saturated MCP events
        Long64_t nEntries = tree->GetEntries();

        for (Long64_t ev = 0; ev < nEntries; ++ev) {
            tree->GetEntry(ev);
            ++nTotal;

            // Beam quality histograms — all events
            if (wc_ok) {
                h2D->Fill(x_trk, y_trk);
                hX->Fill(x_trk);
                hY->Fill(y_trk);
                ++nWC;
                // Radial distance from data-derived centroid (Page 9)
                float dx_r = x_trk - xcf, dy_r = y_trk - ycf;
                hR->Fill(std::sqrt(dx_r*dx_r + dy_r*dy_r));

                // All-events PbGlass scatter — no fiducial or amplitude cuts
                if (hasPb) {
                    float slg_all = 0.f;
                    for (int i = 0; i < 8; ++i) slg_all += lg_peak[i];
                    hPbVsLG_all->Fill(slg_all, sum_pb_v);
                    hSumLG_all->Fill(slg_all);
                }
                // Pedestal RMS accumulation for summary metrics
                if (hasPedRMS) {
                    for (int i = 0; i < 8; ++i) {
                        pedRMSsum[i]  += hg_ped_rms_v[i];
                        ++pedRMScount[i];
                    }
                }
            }
            hMCP->Fill(mcp_peak);
            if (mcp_peak > kMCP1_maxPeak) ++nMCP1sat;
            if (hasMCP2) {
                hMCP2->Fill(mcp2_peak);
                if (mcp2_peak > kMCP2_maxPeak) ++nMCP2sat;
            }

            // Fiducial cut for channel diagnostics
            if (!wc_ok || mcp_peak < kMCP_minPeak_E) continue;
            float dx = x_trk - xcf, dy = y_trk - ycf;
            float r  = std::sqrt(dx*dx + dy*dy);

            // ── Containment and fiducial selection plots (Pages 9 & 10) ──────
            // These use the tighter energy fiducial (kFiducial_r_energy = 2 mm).
            // sum_lg is computed from the already-read lg_peak[] array.
            if (r < static_cast<float>(kFiducial_r_energy)) {
                float sum_lg_v = 0.f;
                for (int i = 0; i < 8; ++i) sum_lg_v += lg_peak[i];
                hSumLG_fid->Fill(sum_lg_v);

                if (hasPb && sum_lg_v > kSumLG_centroid) {
                    // Page 10: PbGlass vs RADiCAL scatter and ratio
                    hPbVsLG->Fill(sum_lg_v, sum_pb_v);
                    float ratio = sum_pb_v / sum_lg_v;
                    hPbRatio->Fill(ratio);
                    float frac  = (sum_pb_v + sum_lg_v) > 0.f
                                  ? sum_pb_v / (sum_pb_v + sum_lg_v) : 0.f;
                    pContVsR->Fill(r, frac);

                    // Page 9: hit maps of selected vs rejected events
                    bool contained = (sum_pb_v < kPb_maxRatio * sum_lg_v);
                    if (contained) {
                        h2Dsel->Fill(x_trk, y_trk);
                        hSumLG_cont->Fill(sum_lg_v);
                        ++nContained;
                    } else {
                        h2DnoC->Fill(x_trk, y_trk);
                    }
                }
            }
            // ─────────────────────────────────────────────────────────────────

            if (r >= static_cast<float>(kFiducial_r_timing)) continue;
            ++nFid;

            // Per-channel amplitude, correlation, walk, and TOT histograms
            for (int i = 0; i < 8; ++i) {
                hHG[i]->Fill(hg_peak[i]);
                hLG[i]->Fill(lg_peak[i]);
                if (lg_peak[i] > 0.) pCorr[i]->Fill(lg_peak[i], hg_peak[i]);
                // Walk diagnostic: fill t_CFD vs 1/A for TIMEABLE pulses only —
                // A > 100 mV AND a correctly-timed crossing (t > -50 ns).  Below
                // ~100 mV the CFD finds no real leading edge and floors near the
                // search-window edge (~-90 ns); that is mistiming, not walk, and it
                // would swamp the genuine (sub-ns) amplitude dependence.
                if (hg_peak[i] > 100.f && hg_cfd[i] > -50.f)
                    pWalk[i]->Fill(1.0f / hg_peak[i], hg_cfd[i]);
                if (hasNewBranches && hg_tot[i] > 0.f
                        && hg_peak[i] > kHG_minPeak && hg_led[i] > -1e5f) {
                    hTOT[i]->Fill(hg_tot[i]);
                    pLEDvsTOT[i]->Fill(hg_tot[i], hg_led[i]);
                }
            }

            // Cache for cut-optimisation scans.
            // Use the full TIMING MCP window (lower + upper) so that σ_t scans
            // reflect real timing-quality events only.  Saturated MCP events
            // have a biased CFD reference and must not enter the timing cache.
            if (mcp_peak >= kMCP1_minPeak && mcp_peak <= kMCP1_maxPeak) {
                EvtData ed;
                ed.r_beam    = r;
                ed.mcp_peak  = mcp_peak;
                ed.mcp2_peak = mcp2_peak;
                for (int i = 0; i < 8; ++i) {
                    ed.hg_peak[i] = hg_peak[i];
                    ed.hg_cfd[i]  = hg_cfd[i];
                }
                cache.push_back(ed);
            }
        }

        fin->Close();

        std::cout << "  Total: " << nTotal
                  << "  WC-ok: " << nWC
                  << "  Fiducial (" << kFiducial_r_timing << " mm): " << nFid
                  << Form("  (%.1f%%)", nFid > 0 ? 100.*nFid/nTotal : 0.)
                  << Form("  MCP1-sat: %ld (%.2f%%)", nMCP1sat, 100.*nMCP1sat/nTotal);
        if (hasPb) {
            // Report containment relative to energy fiducial (2mm) events
            long nEFid = h2Dsel->GetEntries() + h2DnoC->GetEntries();
            std::cout << Form("  Contained: %ld  (%.1f%% of E-fid)",
                nContained, nEFid > 0 ? 100.*nContained/nEFid : 0.);
        }
        std::cout << "\n";

        vEnergy.push_back(rc.energy_GeV);

        // Find the best single channel (lowest CFD-20% sigma at nominal cuts)
        int bestCh = 0;
        double bestSig = 9999.;
        for (int i = 0; i < 8; ++i) {
            double s = TimingSigma(cache, i,
                static_cast<float>(kFiducial_r_timing), kHG_minPeak);
            if (s > 0. && s < bestSig) { bestSig = s; bestCh = i; }
        }
        std::cout << "  Best channel for cut optimisation: "
                  << kCap[bestCh].name
                  << " sigma = " << bestSig*1000. << " ps\n";

        // ---------------------------------------------------------------------
        // Open per-energy PDF
        // ---------------------------------------------------------------------
        TString pdfPath = TString("Analysis/Output/") + rc.label
                          + "/quality_report.pdf";

        // =====================================================================
        // PAGE 1 — Beam quality
        // =====================================================================
        {
            TCanvas c("cQ1", "", 1600, 900);
            c.Divide(3, 2, 0.005, 0.035);

            // Pad 1: 2D hit map
            c.cd(1);
            gPad->SetRightMargin(0.16); StylePad(true);
            h2D->Draw("COLZ");
            // Fiducial circles at data-derived centroid
            TArc* arcE = new TArc(xc, yc, kFiducial_r_energy);
            arcE->SetLineColor(kRed); arcE->SetLineWidth(2);
            arcE->SetLineStyle(1);   arcE->SetFillStyle(0);
            arcE->Draw("SAME");
            TArc* arcT = new TArc(xc, yc, kFiducial_r_timing);
            arcT->SetLineColor(kRed); arcT->SetLineWidth(1);
            arcT->SetLineStyle(2);   arcT->SetFillStyle(0);
            arcT->Draw("SAME");
            { TLatex t; t.SetNDC(); t.SetTextSize(0.052); t.SetTextAlign(22);
              t.DrawLatex(0.44, 0.93,
                  Form("%.0f GeV  beam hit map", rc.energy_GeV)); }

            // Pad 2: x profile
            c.cd(2); StylePad();
            hX->SetLineColor(kBlue+1); hX->SetLineWidth(2); hX->Draw("HIST");
            { TLatex t; t.SetNDC(); t.SetTextSize(0.056); t.SetTextAlign(22);
              t.DrawLatex(0.54, 0.93, Form("%.0f GeV  x_{WC}", rc.energy_GeV)); }
            { TLine* ln = new TLine(xc, 0, xc, hX->GetMaximum());
              ln->SetLineColor(kRed); ln->SetLineStyle(2); ln->Draw("SAME"); }

            // Pad 3: y profile
            c.cd(3); StylePad();
            hY->SetLineColor(kBlue+1); hY->SetLineWidth(2); hY->Draw("HIST");
            { TLatex t; t.SetNDC(); t.SetTextSize(0.056); t.SetTextAlign(22);
              t.DrawLatex(0.54, 0.93, Form("%.0f GeV  y_{WC}", rc.energy_GeV)); }
            { TLine* ln = new TLine(yc, 0, yc, hY->GetMaximum());
              ln->SetLineColor(kRed); ln->SetLineStyle(2); ln->Draw("SAME"); }

            // Pad 4: MCP1 amplitude
            c.cd(4); StylePad();
            hMCP->SetLineColor(kBlue+1); hMCP->SetLineWidth(2); hMCP->Draw("HIST");
            // Lower cut line
            { TLine* ln = new TLine(kMCP1_minPeak, 0,
                                     kMCP1_minPeak, hMCP->GetMaximum());
              ln->SetLineColor(kRed); ln->SetLineStyle(2); ln->SetLineWidth(2);
              ln->Draw("SAME"); }
            // Upper (saturation) cut line
            { TLine* ln = new TLine(kMCP1_maxPeak, 0,
                                     kMCP1_maxPeak, hMCP->GetMaximum());
              ln->SetLineColor(kOrange+1); ln->SetLineStyle(2); ln->SetLineWidth(2);
              ln->Draw("SAME"); }
            { TLatex t; t.SetNDC(); t.SetTextSize(0.056); t.SetTextAlign(22);
              t.DrawLatex(0.54, 0.93, Form("%.0f GeV  MCP1", rc.energy_GeV)); }
            { TLatex t; t.SetNDC(); t.SetTextSize(0.044); t.SetTextColor(kRed);
              t.DrawLatex(0.55, 0.77,
                  Form("lo cut: %.0f mV", (double)kMCP1_minPeak)); }
            { TLatex t; t.SetNDC(); t.SetTextSize(0.044); t.SetTextColor(kOrange+1);
              t.DrawLatex(0.55, 0.68,
                  Form("sat cut: %.0f mV  (%.2f%%)",
                       (double)kMCP1_maxPeak, 100.*nMCP1sat/nTotal)); }

            // Pad 5: MCP2 amplitude (or note if absent)
            c.cd(5); StylePad();
            if (hasMCP2) {
                hMCP2->SetLineColor(kGreen+2); hMCP2->SetLineWidth(2);
                hMCP2->Draw("HIST");
                // Lower cut line
                TLine* ln2lo = new TLine(kMCP2_minPeak, 0,
                                          kMCP2_minPeak, hMCP2->GetMaximum());
                ln2lo->SetLineColor(kRed); ln2lo->SetLineStyle(2);
                ln2lo->SetLineWidth(2);    ln2lo->Draw("SAME");
                // Upper (saturation) cut line
                TLine* ln2hi = new TLine(kMCP2_maxPeak, 0,
                                          kMCP2_maxPeak, hMCP2->GetMaximum());
                ln2hi->SetLineColor(kOrange+1); ln2hi->SetLineStyle(2);
                ln2hi->SetLineWidth(2);         ln2hi->Draw("SAME");
                { TLatex t2; t2.SetNDC(); t2.SetTextSize(0.044); t2.SetTextColor(kRed);
                  t2.DrawLatex(0.55, 0.77,
                      Form("lo cut: %.0f mV", (double)kMCP2_minPeak)); }
                { TLatex t2; t2.SetNDC(); t2.SetTextSize(0.044); t2.SetTextColor(kOrange+1);
                  t2.DrawLatex(0.55, 0.68,
                      Form("sat cut: %.0f mV  (%.2f%%)",
                           (double)kMCP2_maxPeak, 100.*nMCP2sat/nTotal)); }
            }
            { TLatex t; t.SetNDC(); t.SetTextSize(0.056); t.SetTextAlign(22);
              t.DrawLatex(0.54, 0.93, Form("%.0f GeV  MCP2", rc.energy_GeV)); }

            // Pad 6: Run statistics text
            c.cd(6);
            gPad->SetLeftMargin(0.05); gPad->SetTopMargin(0.05);
            TLatex st; st.SetNDC(); st.SetTextSize(0.062);
            st.DrawLatex(0.05, 0.92, Form("%.0f GeV  run statistics", rc.energy_GeV));
            st.SetTextSize(0.052);
            st.DrawLatex(0.05, 0.82, Form("Total events:      %ld", nTotal));
            st.DrawLatex(0.05, 0.73, Form("WC-ok events:      %ld  (%.1f%%)",
                nWC, 100.*nWC/nTotal));
            st.DrawLatex(0.05, 0.64,
                Form("Fiducial (%.1fmm):   %ld  (%.1f%%)",
                     kFiducial_r_timing, nFid, 100.*nFid/nTotal));
            st.SetTextColor(kOrange+1);
            st.DrawLatex(0.05, 0.54,
                Form("MCP1 saturated:    %ld  (%.2f%%)",
                     nMCP1sat, 100.*nMCP1sat/nTotal));
            if (hasMCP2)
                st.DrawLatex(0.05, 0.45,
                    Form("MCP2 saturated:    %ld  (%.2f%%)",
                         nMCP2sat, 100.*nMCP2sat/nTotal));
            st.SetTextSize(0.040); st.SetTextColor(kGray+2);
            st.DrawLatex(0.05, 0.34, "Fiducial circles on hit map:");
            st.DrawLatex(0.05, 0.26, "  solid = energy cut (2.0 mm)");
            st.DrawLatex(0.05, 0.18, "  dashed = timing cut (3.0 mm)");
            st.DrawLatex(0.05, 0.10, "Red/orange lines = MCP lo/sat cuts");

            // Standalone "completeness overview" hero for the headline energy,
            // used as the single figure in the report's quality-report chapter.
            if (rc.energy_GeV > 149. && rc.energy_GeV < 151.)
                c.Print("Analysis/Output/Summary/quality_overview_150.png");
            c.Print(pdfPath + "(");
        }

        // =====================================================================
        // PAGE 2 — HG channel amplitudes + hit efficiency
        // =====================================================================
        {
            TCanvas c("cQ2", "", 2400, 1400);
            c.Divide(3, 3, 0.003, 0.003);

            // Per-channel amplitude distribution
            for (int i = 0; i < 8; ++i) {
                c.cd(i+1); StylePad();
                hHG[i]->SetLineColor(kRChannelCols[i]);
                hHG[i]->SetLineWidth(2); hHG[i]->Draw("HIST");
                // Mark the quality cut threshold
                TLine* cut = new TLine(kHG_minPeak, 0,
                                        kHG_minPeak, hHG[i]->GetMaximum());
                cut->SetLineColor(kRed); cut->SetLineStyle(2);
                cut->SetLineWidth(2);    cut->Draw("SAME");
                // Title and efficiency annotation
                { TLatex t; t.SetNDC(); t.SetTextSize(0.070);
                  t.SetTextAlign(22);
                  t.DrawLatex(0.54, 0.93,
                      Form("%.0f GeV  %s  HG", rc.energy_GeV, kCap[i].name)); }
                long nAbove = 0;
                for (int b = hHG[i]->FindBin(kHG_minPeak);
                         b <= hHG[i]->GetNbinsX(); ++b)
                    nAbove += static_cast<long>(hHG[i]->GetBinContent(b));
                double eff = nFid > 0 ? 100.*nAbove/nFid : 0.;
                { TLatex t; t.SetNDC(); t.SetTextSize(0.060);
                  t.SetTextColor(kRed);
                  t.DrawLatex(0.55, 0.76, Form("eff: %.1f%%", eff)); }
            }

            // Pad 9: efficiency summary table
            c.cd(9); gPad->SetLeftMargin(0.05); gPad->SetTopMargin(0.05);
            TLatex tb; tb.SetNDC(); tb.SetTextSize(0.060);
            tb.DrawLatex(0.05, 0.93, Form("%.0f GeV  HG hit efficiency", rc.energy_GeV));
            tb.SetTextSize(0.052);
            tb.DrawLatex(0.05, 0.83, Form("Cut: A_{HG} > %.0f mV", (double)kHG_minPeak));
            tb.SetTextSize(0.048);
            double yrow = 0.72;
            for (int i = 0; i < 8; ++i) {
                long nAbove = 0;
                for (int b = hHG[i]->FindBin(kHG_minPeak);
                         b <= hHG[i]->GetNbinsX(); ++b)
                    nAbove += static_cast<long>(hHG[i]->GetBinContent(b));
                double eff = nFid > 0 ? 100.*nAbove/nFid : 0.;
                tb.SetTextColor(kRChannelCols[i]);
                tb.DrawLatex(0.05, yrow,
                    Form("%-6s  %.1f%%", kCap[i].name, eff));
                yrow -= 0.085;
            }

            c.Print(pdfPath);
        }

        // =====================================================================
        // PAGE 3 — LG amplitudes + HG/LG correlation
        // =====================================================================
        {
            TCanvas c("cQ3", "", 2400, 1400);
            c.Divide(3, 3, 0.003, 0.003);

            for (int i = 0; i < 8; ++i) {
                c.cd(i+1); StylePad();
                hLG[i]->SetLineColor(kRChannelCols[i]);
                hLG[i]->SetLineWidth(2); hLG[i]->Draw("HIST");
                // Mark the LG quality cut
                TLine* cut = new TLine(kLG_minPeak, 0,
                                        kLG_minPeak, hLG[i]->GetMaximum());
                cut->SetLineColor(kRed); cut->SetLineStyle(2);
                cut->SetLineWidth(2);    cut->Draw("SAME");
                { TLatex t; t.SetNDC(); t.SetTextSize(0.070);
                  t.SetTextAlign(22);
                  t.DrawLatex(0.54, 0.93,
                      Form("%.0f GeV  %s  LG", rc.energy_GeV, kCap[i].name)); }
            }

            // Pad 9: HG/LG correlation for representative channel
            c.cd(9); StylePad();
            if (pCorr[bestCh]->GetEntries() > 20) {
                pCorr[bestCh]->SetMarkerColor(kRChannelCols[bestCh]);
                pCorr[bestCh]->SetMarkerStyle(20);
                pCorr[bestCh]->SetMarkerSize(0.6);
                pCorr[bestCh]->Draw("EP");
                { TLatex t; t.SetNDC(); t.SetTextSize(0.060);
                  t.SetTextAlign(22);
                  t.DrawLatex(0.54, 0.93,
                      Form("%.0f GeV  %s  HG vs LG",
                           rc.energy_GeV, kCap[bestCh].name)); }
            }

            c.Print(pdfPath);
        }

        // =====================================================================
        // PAGE 4 — Time walk diagnostic (t_CFD vs amplitude profiles)
        // =====================================================================
        {
            TCanvas c("cQ4", "", 2400, 1400);
            c.Divide(3, 3, 0.003, 0.003);

            for (int i = 0; i < 8; ++i) {
                c.cd(i+1); StylePad();
                if (pWalk[i]->GetEntries() < 10) continue;

                pWalk[i]->SetMarkerColor(kRChannelCols[i]);
                pWalk[i]->SetMarkerStyle(20);
                pWalk[i]->SetMarkerSize(0.8);
                // Zoom the y-axis around the good-timing band (all channels sit
                // near -26 ns) so the genuine, SUB-NS walk is visible.  The old
                // -20..-100 span was entirely the low-amplitude mistiming floor,
                // which is now excluded by the fill cut.
                pWalk[i]->GetYaxis()->SetRangeUser(-34., -18.);
                pWalk[i]->Draw("EP");

                // First-order walk model t = a + b/A is LINEAR on this 1/A axis.
                // Fit only the timeable region (A > 100 mV, 1/A < 0.01).  A flat
                // profile (b ~ 0) means negligible residual walk.
                TF1 fwalk(Form("fw_%d", i), "pol1", 0.001f, 0.010f);
                fwalk.SetLineColor(kRed+1);
                fwalk.SetLineWidth(1);
                fwalk.SetLineStyle(2);
                double bslope = 0.;
                if (pWalk[i]->GetEntries() > 20) {
                    pWalk[i]->Fit(&fwalk, "RQN");
                    bslope = fwalk.GetParameter(1);   // ns*mV
                }
                fwalk.DrawCopy("SAME");

                { TLatex t; t.SetNDC(); t.SetTextSize(0.060);
                  t.SetTextAlign(22);
                  t.DrawLatex(0.54, 0.93,
                      Form("%.0f GeV  %s  walk", rc.energy_GeV, kCap[i].name)); }
                { TLatex t; t.SetNDC(); t.SetTextSize(0.044);
                  t.SetTextColor(kGray+2);
                  t.DrawLatex(0.15, 0.80, "timeable pulses (A > 100 mV)");
                  t.DrawLatex(0.15, 0.74, Form("b = %.0f ns#upointmV  (pol1 fit)", bslope)); }
            }

            // Pad 9: interpretation note
            c.cd(9); gPad->SetLeftMargin(0.05); gPad->SetTopMargin(0.08);
            TLatex note; note.SetNDC(); note.SetTextSize(0.058);
            note.DrawLatex(0.05, 0.93, "Time-walk diagnostic");
            note.SetTextSize(0.043); note.SetTextColor(kGray+1);
            note.DrawLatex(0.05, 0.84, "Mean t_{CFD} vs 1/A_{HG} for TIMEABLE");
            note.DrawLatex(0.05, 0.78, "pulses (A > 100 mV).  Walk model t = a + b/A");
            note.DrawLatex(0.05, 0.72, "is linear on this axis; slope b is the");
            note.DrawLatex(0.05, 0.66, "amplitude->time coupling.");
            note.DrawLatex(0.05, 0.56, "Result: residual walk is SUB-NS over A = 100-");
            note.DrawLatex(0.05, 0.50, "1000 mV (slightly larger on the Down capillaries)");
            note.DrawLatex(0.05, 0.44, "— no explicit correction needed (headline = CFD-5%).");
            note.SetTextColor(kGray+2);
            note.DrawLatex(0.05, 0.32, "NB: below ~100 mV the CFD has no leading");
            note.DrawLatex(0.05, 0.26, "edge to find and floors near the window edge");
            note.DrawLatex(0.05, 0.20, "(~-90 ns) — mistiming, excluded here.");

            c.Print(pdfPath);
        }

        // =====================================================================
        // PAGE 5 — TOT & LED diagnostics
        // =====================================================================
        if (hasNewBranches)
        {
            TCanvas c("cQ5", "", 2400, 1400);
            c.Divide(3, 3, 0.003, 0.003);

            for (int i = 0; i < 8; ++i) {
                c.cd(i+1); StylePad();
                // TOT distribution
                hTOT[i]->SetLineColor(kRChannelCols[i]);
                hTOT[i]->SetLineWidth(2); hTOT[i]->Draw("HIST");
                { TLatex t; t.SetNDC(); t.SetTextSize(0.060); t.SetTextAlign(22);
                  t.DrawLatex(0.54, 0.93,
                      Form("%.0f GeV  %s  TOT", rc.energy_GeV, kCap[i].name)); }
                if (hTOT[i]->GetEntries() > 20) {
                    double mean = hTOT[i]->GetMean();
                    TLatex t; t.SetNDC(); t.SetTextSize(0.056); t.SetTextColor(kGray+2);
                    t.DrawLatex(0.55, 0.76, Form("<TOT> = %.1f ns", mean));
                }
            }

            // Pad 9: LED vs TOT profile for best channel
            c.cd(9); StylePad();
            if (pLEDvsTOT[bestCh]->GetEntries() > 20) {
                pLEDvsTOT[bestCh]->SetMarkerColor(kRChannelCols[bestCh]);
                pLEDvsTOT[bestCh]->SetMarkerStyle(20);
                pLEDvsTOT[bestCh]->SetMarkerSize(0.8);
                pLEDvsTOT[bestCh]->Draw("EP");
                // Fit: t_LED = a + b*TOT — the slope b is used by method M5
                TF1 flt(Form("flt_%d", bestCh), "pol1",
                        hTOT[bestCh]->GetXaxis()->GetXmin(),
                        hTOT[bestCh]->GetXaxis()->GetXmax());
                flt.SetLineColor(kRed+1); flt.SetLineWidth(1); flt.SetLineStyle(2);
                pLEDvsTOT[bestCh]->Fit(&flt, "RQN");
                flt.DrawCopy("SAME");
                { TLatex t; t.SetNDC(); t.SetTextSize(0.055); t.SetTextAlign(22);
                  t.DrawLatex(0.54, 0.93,
                      Form("%.0f GeV  %s  LED vs TOT",
                           rc.energy_GeV, kCap[bestCh].name)); }
                { TLatex t; t.SetNDC(); t.SetTextSize(0.048); t.SetTextColor(kGray+2);
                  t.DrawLatex(0.15, 0.76, "Slope used by M5 (TOT correction)"); }
            }

            c.Print(pdfPath);
        }

        // =====================================================================
        // PAGE 6 — Cut optimisation: MCP amplitude threshold
        // The MCP cut is the primary event-quality gate for timing — vary it
        // to find the trade-off between resolution improvement and efficiency loss.
        // =====================================================================
        {
            // Cache already contains all events with mcp_peak >= kMCP1_minPeak;
            // to scan below the nominal cut we need events down to ~0 mV.
            // Use the raw cache reference (0 mV cut) as the 100% efficiency baseline.
            // NOTE: the cache was filled with mcp_peak >= kMCP1_minPeak, so scans
            // BELOW that threshold cannot be evaluated.  The scan shows the gain
            // from tightening the cut above the nominal value.
            const float mcpCuts[] = {50.f, 75.f, 100.f, 125.f, 150.f,
                                      175.f, 200.f, 250.f, 300.f, 400.f};
            const int nMCP = 10;

            TGraph gSigVsMCP, gEffVsMCP;
            // Reference: events passing nominal timing cuts (fixed r and HG cuts)
            int nRefMCP = CountEvents(cache, bestCh,
                                       static_cast<float>(kFiducial_r_timing),
                                       kHG_minPeak);

            for (int im = 0; im < nMCP; ++im) {
                float M   = mcpCuts[im];
                double s  = TimingSigmaMCP(cache, bestCh,
                                            static_cast<float>(kFiducial_r_timing),
                                            M, kHG_minPeak);
                int    n  = CountEventsMCP(cache, bestCh,
                                            static_cast<float>(kFiducial_r_timing),
                                            M, kHG_minPeak);
                double eff = (nRefMCP > 0) ? 100.*n/nRefMCP : 0.;
                if (s > 0.) {
                    gSigVsMCP.SetPoint(gSigVsMCP.GetN(), M, s*1000.);
                    gEffVsMCP.SetPoint(gEffVsMCP.GetN(), M, eff);
                }
            }

            TCanvas c("cQ6mcp", "", 1400, 600);
            c.Divide(2, 1, 0.01, 0.01);

            // σ vs MCP threshold
            c.cd(1);
            gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14);
            gPad->SetRightMargin(0.05); gPad->SetTopMargin(0.12);
            gPad->SetTickx(1); gPad->SetTicky(1);

            if (gSigVsMCP.GetN() > 1) {
                double ylo = 9999., yhi = 0.;
                for (int p = 0; p < gSigVsMCP.GetN(); ++p) {
                    ylo = std::min(ylo, gSigVsMCP.GetY()[p]);
                    yhi = std::max(yhi, gSigVsMCP.GetY()[p]);
                }
                TH1F* fm1 = gPad->DrawFrame(30., ylo*0.85, 420., yhi*1.20,
                    ";MCP1 amplitude cut (mV);#sigma_{t} (ps)");
                fm1->GetXaxis()->SetTitleSize(0.052);
                fm1->GetYaxis()->SetTitleSize(0.052);
                gSigVsMCP.SetMarkerStyle(20); gSigVsMCP.SetMarkerSize(1.2);
                gSigVsMCP.SetLineWidth(2);
                gSigVsMCP.SetMarkerColor(kBlue+1); gSigVsMCP.SetLineColor(kBlue+1);
                gSigVsMCP.Draw("PL SAME");
                TLine* lM = new TLine(kMCP1_minPeak, ylo*0.85,
                                       kMCP1_minPeak, yhi*1.20);
                lM->SetLineColor(kRed); lM->SetLineStyle(2); lM->Draw("SAME");
            }
            { TLatex t; t.SetNDC(); t.SetTextSize(0.052); t.SetTextAlign(22);
              t.DrawLatex(0.54, 0.93,
                  Form("%.0f GeV  %s  MCP scan", rc.energy_GeV, kCap[bestCh].name)); }
            { TLatex t; t.SetNDC(); t.SetTextSize(0.044); t.SetTextColor(kRed);
              t.DrawLatex(0.55, 0.76,
                  Form("nominal cut: %.0f mV", (double)kMCP1_minPeak)); }

            // Efficiency vs MCP threshold
            c.cd(2);
            gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14);
            gPad->SetRightMargin(0.05); gPad->SetTopMargin(0.12);
            gPad->SetTickx(1); gPad->SetTicky(1);

            if (gEffVsMCP.GetN() > 1) {
                TH1F* fm2 = gPad->DrawFrame(30., 0., 420., 115.,
                    ";MCP1 amplitude cut (mV);Efficiency (%)");
                fm2->GetXaxis()->SetTitleSize(0.052);
                fm2->GetYaxis()->SetTitleSize(0.052);
                gEffVsMCP.SetMarkerStyle(20); gEffVsMCP.SetMarkerSize(1.2);
                gEffVsMCP.SetLineWidth(2);
                gEffVsMCP.SetMarkerColor(kGreen+2); gEffVsMCP.SetLineColor(kGreen+2);
                gEffVsMCP.Draw("PL SAME");
                TLine* lM2 = new TLine(kMCP1_minPeak, 0., kMCP1_minPeak, 115.);
                lM2->SetLineColor(kRed); lM2->SetLineStyle(2); lM2->Draw("SAME");
            }
            { TLatex t; t.SetNDC(); t.SetTextSize(0.052); t.SetTextAlign(22);
              t.DrawLatex(0.54, 0.93,
                  Form("%.0f GeV  MCP amplitude efficiency", rc.energy_GeV)); }

            c.Print(pdfPath);
        }

        // =====================================================================
        // PAGE 7 — Cut optimisation: fiducial radius
        // =====================================================================
        {
            // Scan fiducial radius for the best channel
            const float radii[] = {0.5f, 0.75f, 1.0f, 1.25f, 1.5f,
                                    1.75f, 2.0f, 2.5f, 3.0f, 3.5f, 4.0f, 5.0f};
            const int nR = 12;

            TGraph gSigVsR, gEffVsR;
            int nRef = CountEvents(cache, bestCh, 100.f, 0.f); // all within large radius

            for (int ir = 0; ir < nR; ++ir) {
                float R   = radii[ir];
                double s  = TimingSigma(cache, bestCh, R, kHG_minPeak);
                int    n  = CountEvents(cache, bestCh, R, kHG_minPeak);
                double eff = (nRef > 0) ? 100.*n/nRef : 0.;
                if (s > 0.) {
                    gSigVsR.SetPoint(gSigVsR.GetN(), R, s*1000.);  // [ps]
                    gEffVsR.SetPoint(gEffVsR.GetN(), R, eff);
                }
            }

            TCanvas c("cQ7", "", 1400, 600);
            c.Divide(2, 1, 0.01, 0.01);

            // σ vs radius
            c.cd(1);
            gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14);
            gPad->SetRightMargin(0.05); gPad->SetTopMargin(0.12);
            gPad->SetTickx(1); gPad->SetTicky(1);

            if (gSigVsR.GetN() > 1) {
                double ylo = 9999., yhi = 0.;
                for (int p = 0; p < gSigVsR.GetN(); ++p) {
                    ylo = std::min(ylo, gSigVsR.GetY()[p]);
                    yhi = std::max(yhi, gSigVsR.GetY()[p]);
                }
                TH1F* fr1 = gPad->DrawFrame(0., ylo*0.85, 5.5, yhi*1.20,
                    ";Fiducial radius (mm);#sigma_{t} (ps)");
                fr1->GetXaxis()->SetTitleSize(0.052);
                fr1->GetYaxis()->SetTitleSize(0.052);
                gSigVsR.SetMarkerStyle(20); gSigVsR.SetMarkerSize(1.2);
                gSigVsR.SetLineWidth(2);
                gSigVsR.SetMarkerColor(kBlue+1); gSigVsR.SetLineColor(kBlue+1);
                gSigVsR.Draw("PL SAME");
                // Mark nominal timing cut
                TLine* lR = new TLine(kFiducial_r_timing, ylo*0.85,
                                       kFiducial_r_timing, yhi*1.20);
                lR->SetLineColor(kRed); lR->SetLineStyle(2); lR->Draw("SAME");
            }
            { TLatex t; t.SetNDC(); t.SetTextSize(0.052); t.SetTextAlign(22);
              t.DrawLatex(0.54, 0.93,
                  Form("%.0f GeV  %s  fiducial scan", rc.energy_GeV, kCap[bestCh].name)); }
            { TLatex t; t.SetNDC(); t.SetTextSize(0.044); t.SetTextColor(kRed);
              t.DrawLatex(0.55, 0.76,
                  Form("dashed = %.1f mm cut", kFiducial_r_timing)); }

            // Efficiency vs radius
            c.cd(2);
            gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14);
            gPad->SetRightMargin(0.05); gPad->SetTopMargin(0.12);
            gPad->SetTickx(1); gPad->SetTicky(1);

            if (gEffVsR.GetN() > 1) {
                TH1F* fr2 = gPad->DrawFrame(0., 0., 5.5, 115.,
                    ";Fiducial radius (mm);Efficiency (%)");
                fr2->GetXaxis()->SetTitleSize(0.052);
                fr2->GetYaxis()->SetTitleSize(0.052);
                gEffVsR.SetMarkerStyle(20); gEffVsR.SetMarkerSize(1.2);
                gEffVsR.SetLineWidth(2);
                gEffVsR.SetMarkerColor(kGreen+2); gEffVsR.SetLineColor(kGreen+2);
                gEffVsR.Draw("PL SAME");
                TLine* lR2 = new TLine(kFiducial_r_timing, 0.,
                                        kFiducial_r_timing, 115.);
                lR2->SetLineColor(kRed); lR2->SetLineStyle(2); lR2->Draw("SAME");
            }
            { TLatex t; t.SetNDC(); t.SetTextSize(0.052); t.SetTextAlign(22);
              t.DrawLatex(0.54, 0.93,
                  Form("%.0f GeV  fiducial efficiency", rc.energy_GeV)); }

            c.Print(pdfPath);
        }

        // =====================================================================
        // PAGE 8 — Cut optimisation: HG amplitude threshold
        // =====================================================================
        {
            const float ampCuts[] = {0.f, 5.f, 10.f, 15.f, 20.f, 25.f,
                                      30.f, 40.f, 50.f, 75.f, 100.f, 150.f};
            const int nA = 12;

            TGraph gSigVsA, gEffVsA;
            // Reference: all fiducial events with no amplitude cut
            int nRef = CountEvents(cache, bestCh,
                                    static_cast<float>(kFiducial_r_timing), 0.f);

            for (int ia = 0; ia < nA; ++ia) {
                float A   = ampCuts[ia];
                double s  = TimingSigma(cache, bestCh,
                                         static_cast<float>(kFiducial_r_timing), A);
                int    n  = CountEvents(cache, bestCh,
                                         static_cast<float>(kFiducial_r_timing), A);
                double eff = (nRef > 0) ? 100.*n/nRef : 0.;
                if (s > 0.) {
                    gSigVsA.SetPoint(gSigVsA.GetN(), A, s*1000.);
                    gEffVsA.SetPoint(gEffVsA.GetN(), A, eff);
                }
            }

            TCanvas c("cQ8", "", 1400, 600);
            c.Divide(2, 1, 0.01, 0.01);

            // σ vs amplitude threshold
            c.cd(1);
            gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14);
            gPad->SetRightMargin(0.05); gPad->SetTopMargin(0.12);
            gPad->SetTickx(1); gPad->SetTicky(1);

            if (gSigVsA.GetN() > 1) {
                double ylo = 9999., yhi = 0.;
                for (int p = 0; p < gSigVsA.GetN(); ++p) {
                    ylo = std::min(ylo, gSigVsA.GetY()[p]);
                    yhi = std::max(yhi, gSigVsA.GetY()[p]);
                }
                TH1F* fa1 = gPad->DrawFrame(-5., ylo*0.85, 160., yhi*1.20,
                    ";Min HG amplitude (mV);#sigma_{t} (ps)");
                fa1->GetXaxis()->SetTitleSize(0.052);
                fa1->GetYaxis()->SetTitleSize(0.052);
                gSigVsA.SetMarkerStyle(20); gSigVsA.SetMarkerSize(1.2);
                gSigVsA.SetLineWidth(2);
                gSigVsA.SetMarkerColor(kBlue+1); gSigVsA.SetLineColor(kBlue+1);
                gSigVsA.Draw("PL SAME");
                TLine* lA = new TLine(kHG_minPeak, ylo*0.85,
                                       kHG_minPeak, yhi*1.20);
                lA->SetLineColor(kRed); lA->SetLineStyle(2); lA->Draw("SAME");
            }
            { TLatex t; t.SetNDC(); t.SetTextSize(0.052); t.SetTextAlign(22);
              t.DrawLatex(0.54, 0.93,
                  Form("%.0f GeV  %s  amplitude scan",
                       rc.energy_GeV, kCap[bestCh].name)); }
            { TLatex t; t.SetNDC(); t.SetTextSize(0.044); t.SetTextColor(kRed);
              t.DrawLatex(0.55, 0.76,
                  Form("dashed = %.0f mV cut", (double)kHG_minPeak)); }

            // Efficiency vs amplitude threshold
            c.cd(2);
            gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14);
            gPad->SetRightMargin(0.05); gPad->SetTopMargin(0.12);
            gPad->SetTickx(1); gPad->SetTicky(1);

            if (gEffVsA.GetN() > 1) {
                TH1F* fa2 = gPad->DrawFrame(-5., 0., 160., 115.,
                    ";Min HG amplitude (mV);Efficiency (%)");
                fa2->GetXaxis()->SetTitleSize(0.052);
                fa2->GetYaxis()->SetTitleSize(0.052);
                gEffVsA.SetMarkerStyle(20); gEffVsA.SetMarkerSize(1.2);
                gEffVsA.SetLineWidth(2);
                gEffVsA.SetMarkerColor(kGreen+2); gEffVsA.SetLineColor(kGreen+2);
                gEffVsA.Draw("PL SAME");
                TLine* lA2 = new TLine(kHG_minPeak, 0.,
                                        kHG_minPeak, 115.);
                lA2->SetLineColor(kRed); lA2->SetLineStyle(2); lA2->Draw("SAME");
            }
            { TLatex t; t.SetNDC(); t.SetTextSize(0.052); t.SetTextAlign(22);
              t.DrawLatex(0.54, 0.93,
                  Form("%.0f GeV  amplitude efficiency", rc.energy_GeV)); }

            c.Print(pdfPath);   // page 8 — more pages follow
        }

        // =====================================================================
        // PAGE 9 — Fiducial selection overview
        // Shows the complete event selection pipeline: WC hit map at each cut
        // stage, radial beam profile with fiducial lines, cut flow statistics,
        // and the impact of the PbGlass containment cut.
        // =====================================================================
        {
            TCanvas c("cQ9", "", 2400, 1400);
            c.Divide(3, 2, 0.005, 0.035);

            // Pad 1: WC hit map for all WC-ok events (same as Page 1 top-left)
            c.cd(1);
            gPad->SetRightMargin(0.16); StylePad(true);
            h2D->Draw("COLZ");
            TArc* a1E = new TArc(xc, yc, kFiducial_r_energy);
            a1E->SetLineColor(kRed); a1E->SetLineWidth(2);
            a1E->SetLineStyle(1);    a1E->SetFillStyle(0); a1E->Draw("SAME");
            TArc* a1T = new TArc(xc, yc, kFiducial_r_timing);
            a1T->SetLineColor(kRed); a1T->SetLineWidth(1);
            a1T->SetLineStyle(2);    a1T->SetFillStyle(0); a1T->Draw("SAME");
            { TLatex t; t.SetNDC(); t.SetTextSize(0.050); t.SetTextAlign(22);
              t.DrawLatex(0.44, 0.93,
                  Form("%.0f GeV  all WC-ok events", rc.energy_GeV)); }

            // Pad 2: 1D radial profile with fiducial cut lines
            c.cd(2); StylePad();
            hR->SetLineColor(kBlue+1); hR->SetLineWidth(2); hR->Draw("HIST");
            { TLine* lE = new TLine(kFiducial_r_energy, 0,
                                     kFiducial_r_energy, hR->GetMaximum());
              lE->SetLineColor(kRed); lE->SetLineStyle(1); lE->SetLineWidth(2);
              lE->Draw("SAME"); }
            { TLine* lT = new TLine(kFiducial_r_timing, 0,
                                     kFiducial_r_timing, hR->GetMaximum());
              lT->SetLineColor(kRed); lT->SetLineStyle(2); lT->SetLineWidth(2);
              lT->Draw("SAME"); }
            { TLatex t; t.SetNDC(); t.SetTextSize(0.050); t.SetTextAlign(22);
              t.DrawLatex(0.54, 0.93,
                  Form("%.0f GeV  radial beam profile", rc.energy_GeV)); }
            { TLatex t; t.SetNDC(); t.SetTextSize(0.044); t.SetTextColor(kRed);
              t.DrawLatex(0.55, 0.78,
                  Form("solid = %.1f mm (energy)", kFiducial_r_energy));
              t.DrawLatex(0.55, 0.70,
                  Form("dashed = %.1f mm (timing)", kFiducial_r_timing)); }

            // Pad 3: Cut flow table
            c.cd(3);
            gPad->SetLeftMargin(0.05); gPad->SetTopMargin(0.05);
            TLatex cf; cf.SetNDC(); cf.SetTextSize(0.064);
            cf.DrawLatex(0.05, 0.93, Form("%.0f GeV  cut flow", rc.energy_GeV));
            cf.SetTextSize(0.052);
            double yc3 = 0.82;
            auto cfRow = [&](const char* label, long n, long denom) {
                double pct = denom > 0 ? 100.*n/denom : 0.;
                cf.DrawLatex(0.05, yc3, Form("%-20s %7ld  (%.1f%%)", label, n, pct));
                yc3 -= 0.10;
            };
            cfRow("Total events",      nTotal, nTotal);
            cfRow("WC-ok",             nWC,    nTotal);
            cfRow("MCP (E cut)",
                  // count WC+MCP events directly from hR integral as proxy;
                  // use nFid-based estimate from timing fiducial
                  static_cast<long>(hR->GetEntries()), nWC);
            cfRow(Form("Fid T (%.1fmm)", kFiducial_r_timing), nFid, nTotal);
            long nEfid = static_cast<long>(h2Dsel->GetEntries() +
                                            h2DnoC->GetEntries());
            cfRow(Form("Fid E (%.1fmm)", kFiducial_r_energy), nEfid, nTotal);
            if (hasPb) {
                cfRow("Fid E + contained", nContained, nTotal);
                cf.SetTextSize(0.042); cf.SetTextColor(kGray+2);
                cf.DrawLatex(0.05, yc3 - 0.04,
                    Form("Containment cut: sum_pb < %.0f%% sum_lg",
                         100.*kPb_maxRatio));
            }

            // Pad 4: Hit map of events passing energy fiducial + containment
            c.cd(4);
            gPad->SetRightMargin(0.16); StylePad(true);
            if (h2Dsel->GetEntries() > 0) {
                h2Dsel->Draw("COLZ");
            } else {
                TLatex t; t.SetNDC(); t.SetTextSize(0.06); t.SetTextAlign(22);
                t.DrawLatex(0.5, 0.5, "No pb_peak branch");
            }
            TArc* a4E = new TArc(xc, yc, kFiducial_r_energy);
            a4E->SetLineColor(kRed); a4E->SetLineWidth(2);
            a4E->SetFillStyle(0); a4E->Draw("SAME");
            { TLatex t; t.SetNDC(); t.SetTextSize(0.050); t.SetTextAlign(22);
              t.DrawLatex(0.44, 0.93,
                  Form("%.0f GeV  fid+contained events", rc.energy_GeV)); }

            // Pad 5: sum_lg before vs after containment (overlay)
            c.cd(5); StylePad();
            hSumLG_fid->SetLineColor(kBlue+1); hSumLG_fid->SetLineWidth(2);
            hSumLG_fid->Draw("HIST");
            hSumLG_cont->SetLineColor(kGreen+2); hSumLG_cont->SetLineWidth(2);
            hSumLG_cont->SetLineStyle(2); hSumLG_cont->Draw("HIST SAME");
            { TLatex t; t.SetNDC(); t.SetTextSize(0.050); t.SetTextAlign(22);
              t.DrawLatex(0.54, 0.93,
                  Form("%.0f GeV  #SigmaLG (energy fiducial)", rc.energy_GeV)); }
            { TLegend* leg = new TLegend(0.45, 0.65, 0.88, 0.78);
              leg->SetBorderSize(0); leg->SetTextSize(0.044);
              leg->AddEntry(hSumLG_fid,  "All fid-E events", "l");
              leg->AddEntry(hSumLG_cont, "After containment", "l");
              leg->Draw(); }

            // Pad 6: Hit map of events failing containment (punch-through / edge)
            c.cd(6);
            gPad->SetRightMargin(0.16); StylePad(true);
            if (h2DnoC->GetEntries() > 0) {
                h2DnoC->Draw("COLZ");
            } else {
                TLatex t; t.SetNDC(); t.SetTextSize(0.06); t.SetTextAlign(22);
                t.DrawLatex(0.5, 0.5, hasPb ? "All events contained" : "No pb_peak branch");
            }
            TArc* a6E = new TArc(xc, yc, kFiducial_r_energy);
            a6E->SetLineColor(kRed); a6E->SetLineWidth(2);
            a6E->SetFillStyle(0); a6E->Draw("SAME");
            { TLatex t; t.SetNDC(); t.SetTextSize(0.050); t.SetTextAlign(22);
              t.DrawLatex(0.44, 0.93,
                  Form("%.0f GeV  containment-REJECTED", rc.energy_GeV)); }

            c.Print(pdfPath);
        }

        // =====================================================================
        // PAGE 10 — PbGlass shower containment  (3 × 2 layout)
        //
        // Pad 1 (ALL events):  sum_pb vs sum_lg for every WC-ok event.
        //   Shows the full population: beam halo/no-signal cluster near origin;
        //   full-energy showers in the upper-right; the cut line separates
        //   contained events (below) from punch-through (above).
        //
        // Pad 2 (energy fiducial):  same axes but only events inside the 2 mm
        //   energy fiducial where sum_lg > kSumLG_centroid.  Much cleaner
        //   population — directly shows the events used for energy analysis.
        //
        // Pads 3-6: ratio distribution, containment-fraction vs radius,
        //   sum_lg spectra, and annotation.
        // =====================================================================
        {
            TCanvas c("cQ10", "", 2400, 1200);
            c.Divide(3, 2, 0.008, 0.008);

            // Helper lambda: draw the standard cut line on the current pad
            auto drawCutLine = [&](TH2F* h) {
                if (!h || h->GetEntries() < 2) return;
                double xmx = h->GetXaxis()->GetXmax();
                TLine* lc = new TLine(0, 0, xmx, kPb_maxRatio * xmx);
                lc->SetLineColor(kRed); lc->SetLineStyle(2);
                lc->SetLineWidth(2);    lc->Draw("SAME");
            };

            // ── Pad 1: ALL WC-ok events ──────────────────────────────────────
            c.cd(1);
            gPad->SetRightMargin(0.16); gPad->SetLeftMargin(0.13);
            gPad->SetBottomMargin(0.13); StylePad(true);
            if (hasPb && hPbVsLG_all->GetEntries() > 10) {
                hPbVsLG_all->Draw("COLZ");
                drawCutLine(hPbVsLG_all);
            } else {
                TLatex t; t.SetNDC(); t.SetTextSize(0.06); t.SetTextAlign(22);
                t.DrawLatex(0.5, 0.5, hasPb ? "Too few events" : "No pb_peak branch");
            }
            { TLatex t; t.SetNDC(); t.SetTextSize(0.050); t.SetTextAlign(22);
              t.DrawLatex(0.50, 0.94,
                  Form("%.0f GeV  ALL WC-ok events", rc.energy_GeV)); }
            { TLatex t; t.SetNDC(); t.SetTextSize(0.042); t.SetTextColor(kRed);
              t.DrawLatex(0.14, 0.85,
                  Form("cut: sum_pb < %.0f%% #times sum_lg", 100.*kPb_maxRatio)); }

            // ── Pad 2: Energy fiducial events ────────────────────────────────
            c.cd(2);
            gPad->SetRightMargin(0.16); gPad->SetLeftMargin(0.13);
            gPad->SetBottomMargin(0.13); StylePad(true);
            if (hasPb && hPbVsLG->GetEntries() > 10) {
                hPbVsLG->Draw("COLZ");
                drawCutLine(hPbVsLG);
            } else {
                TLatex t; t.SetNDC(); t.SetTextSize(0.06); t.SetTextAlign(22);
                t.DrawLatex(0.5, 0.5, hasPb ? "Too few events" : "No pb_peak branch");
            }
            { TLatex t; t.SetNDC(); t.SetTextSize(0.050); t.SetTextAlign(22);
              t.DrawLatex(0.50, 0.94,
                  Form("%.0f GeV  energy fiducial (%.1f mm)", rc.energy_GeV,
                       kFiducial_r_energy)); }

            // ── Pad 3: PbGlass/RADiCAL ratio distribution ────────────────────
            c.cd(3); StylePad();
            if (hasPb && hPbRatio->GetEntries() > 10) {
                hPbRatio->SetLineColor(kBlue+1);
                hPbRatio->SetLineWidth(2); hPbRatio->Draw("HIST");
                TLine* lrr = new TLine(kPb_maxRatio, 0,
                                        kPb_maxRatio, hPbRatio->GetMaximum());
                lrr->SetLineColor(kRed); lrr->SetLineStyle(2);
                lrr->SetLineWidth(2);    lrr->Draw("SAME");
                { TLatex t; t.SetNDC(); t.SetTextSize(0.046); t.SetTextColor(kRed);
                  t.DrawLatex(0.52, 0.78,
                      Form("cut: ratio < %.2f", (double)kPb_maxRatio)); }
                int binCut    = hPbRatio->FindBin(kPb_maxRatio);
                double nPasd  = hPbRatio->Integral(1, binCut - 1);
                double nAllR  = hPbRatio->Integral();
                if (nAllR > 0) {
                    TLatex t; t.SetNDC(); t.SetTextSize(0.044); t.SetTextColor(kGray+2);
                    t.DrawLatex(0.52, 0.68, Form("%.1f%% pass", 100.*nPasd/nAllR));
                }
            } else {
                TLatex t; t.SetNDC(); t.SetTextSize(0.06); t.SetTextAlign(22);
                t.DrawLatex(0.5, 0.5, hasPb ? "Too few events" : "No pb_peak branch");
            }
            { TLatex t; t.SetNDC(); t.SetTextSize(0.050); t.SetTextAlign(22);
              t.DrawLatex(0.54, 0.94,
                  Form("%.0f GeV  containment ratio", rc.energy_GeV)); }

            // ── Pad 4: Containment fraction vs radial distance ───────────────
            c.cd(4); StylePad();
            if (hasPb && pContVsR->GetEntries() > 20) {
                pContVsR->SetMarkerColor(kBlue+1);
                pContVsR->SetMarkerStyle(20);
                pContVsR->SetMarkerSize(0.9);
                pContVsR->SetLineColor(kBlue+1);
                pContVsR->SetMinimum(0.);
                pContVsR->Draw("EP");
                { TLine* lE = new TLine(kFiducial_r_energy, 0,
                                         kFiducial_r_energy,
                                         pContVsR->GetMaximum() * 1.1);
                  lE->SetLineColor(kRed); lE->SetLineStyle(1); lE->Draw("SAME"); }
                { TLine* lT = new TLine(kFiducial_r_timing, 0,
                                         kFiducial_r_timing,
                                         pContVsR->GetMaximum() * 1.1);
                  lT->SetLineColor(kRed); lT->SetLineStyle(2); lT->Draw("SAME"); }
            } else {
                TLatex t; t.SetNDC(); t.SetTextSize(0.06); t.SetTextAlign(22);
                t.DrawLatex(0.5, 0.5, hasPb ? "Too few events" : "No pb_peak branch");
            }
            { TLatex t; t.SetNDC(); t.SetTextSize(0.050); t.SetTextAlign(22);
              t.DrawLatex(0.54, 0.94,
                  Form("%.0f GeV  frac vs radius", rc.energy_GeV)); }
            { TLatex t; t.SetNDC(); t.SetTextSize(0.040); t.SetTextColor(kGray+2);
              t.DrawLatex(0.14, 0.18, "solid = E fid  dashed = T fid"); }

            // ── Pad 5: sum_lg spectra — ALL vs energy-fiducial vs contained ──
            c.cd(5); StylePad();
            hSumLG_all->SetLineColor(kGray+1); hSumLG_all->SetLineWidth(1);
            // Rescale to same area for shape comparison
            double scaleAll  = 1.;
            double scaleFid  = 1.;
            double scaleCont = 1.;
            if (hSumLG_all->Integral()  > 0) scaleAll  = 1./hSumLG_all->Integral();
            if (hSumLG_fid->Integral()  > 0) scaleFid  = 1./hSumLG_fid->Integral();
            if (hSumLG_cont->Integral() > 0) scaleCont = 1./hSumLG_cont->Integral();
            TH1F* hA = (TH1F*)hSumLG_all->Clone("_tmpA");
            TH1F* hF = (TH1F*)hSumLG_fid->Clone("_tmpF");
            TH1F* hC = (TH1F*)hSumLG_cont->Clone("_tmpC");
            hA->Scale(scaleAll); hF->Scale(scaleFid); hC->Scale(scaleCont);
            hA->SetLineColor(kGray+1);  hA->SetLineWidth(1);
            hF->SetLineColor(kBlue+1);  hF->SetLineWidth(2);
            hC->SetLineColor(kGreen+2); hC->SetLineWidth(2); hC->SetLineStyle(2);
            double ymxS = std::max({hA->GetMaximum(), hF->GetMaximum(), hC->GetMaximum()});
            hA->SetMaximum(ymxS * 1.15);
            hA->Draw("HIST"); hF->Draw("HIST SAME"); hC->Draw("HIST SAME");
            { TLegend* leg = new TLegend(0.38, 0.68, 0.88, 0.84);
              leg->SetBorderSize(0); leg->SetTextSize(0.040);
              leg->AddEntry(hA, "All WC-ok",         "l");
              leg->AddEntry(hF, "Energy fiducial",   "l");
              leg->AddEntry(hC, "After containment", "l");
              leg->Draw(); }
            { TLatex t; t.SetNDC(); t.SetTextSize(0.048); t.SetTextAlign(22);
              t.DrawLatex(0.54, 0.94,
                  Form("%.0f GeV  #SigmaLG shape comparison", rc.energy_GeV)); }

            // ── Pad 6: Annotation ─────────────────────────────────────────────
            c.cd(6);
            gPad->SetLeftMargin(0.05); gPad->SetTopMargin(0.05);
            TLatex note; note.SetNDC(); note.SetTextSize(0.058);
            note.DrawLatex(0.05, 0.94, "PbGlass containment");
            note.SetTextSize(0.046); note.SetTextColor(kGray+1);
            note.DrawLatex(0.05, 0.85, "PbGlass is BEHIND the RADiCAL module.");
            note.DrawLatex(0.05, 0.77, "Pad 1: ALL WC-ok events (the full");
            note.DrawLatex(0.05, 0.70, "  beam population: halo + signal).");
            note.DrawLatex(0.05, 0.62, "Pad 2: Energy fiducial events:");
            note.DrawLatex(0.05, 0.55, "  central showers, clean signal region.");
            note.DrawLatex(0.05, 0.46, "Events BELOW cut line are contained");
            note.DrawLatex(0.05, 0.39, "in RADiCAL.  ABOVE = punch-through.");
            note.SetTextColor(kBlack);
            note.DrawLatex(0.05, 0.30,
                Form("Cut: sum_pb < %.0f%% #times sum_lg", 100.*kPb_maxRatio));
            note.DrawLatex(0.05, 0.22,
                Form("(kPb_maxRatio = %.2f)", (double)kPb_maxRatio));
            if (hasPb) {
                long nEfid2  = static_cast<long>(h2Dsel->GetEntries() +
                                                  h2DnoC->GetEntries());
                note.SetTextSize(0.050); note.SetTextColor(kBlue+1);
                note.DrawLatex(0.05, 0.11,
                    Form("E-fid pass: %ld/%ld  (%.1f%%)",
                         nContained, nEfid2,
                         nEfid2 > 0 ? 100.*nContained/nEfid2 : 0.));
            }

            c.Print(pdfPath + ")");   // close per-energy PDF (last page)
        }

        // ── Collect channel validation metrics for the summary PDF ────────────
        for (int i = 0; i < 8; ++i) {
            // Include the OVERFLOW bin (nbins+1): hHG is booked 0-500 mV, but at
            // high energy HG amplitudes exceed 500 mV and land in overflow.
            // Omitting it made the hit efficiency collapse to ~0 above ~40 GeV.
            long nAbove = (long)hHG[i]->Integral(hHG[i]->FindBin(kHG_minPeak),
                                                 hHG[i]->GetNbinsX() + 1);
            sumHitEff[nRunsDone][i] = nFid > 0 ? 100.*nAbove/nFid : 0.;

            hHG[i]->GetXaxis()->SetRangeUser(kHG_minPeak,
                                              hHG[i]->GetXaxis()->GetXmax());
            sumMeanHG[nRunsDone][i] = hHG[i]->GetMean();
            hHG[i]->GetXaxis()->SetRange(0, 0);

            sumPedRMS[nRunsDone][i] = pedRMScount[i] > 0
                ? pedRMSsum[i] / pedRMScount[i] : 0.;

            double s = TimingSigma(cache, i,
                static_cast<float>(kFiducial_r_timing), kHG_minPeak);
            sumTimeSig[nRunsDone][i] = (s > 0.) ? s * 1000. : 0.;
        }
        {
            long nEfidTot = static_cast<long>(h2Dsel->GetEntries() +
                                               h2DnoC->GetEntries());
            sumContain[nRunsDone] = nEfidTot > 0
                ? 100.*nContained/nEfidTot : 0.;
        }
        ++nRunsDone;

        std::cout << "[qualityPlots] " << rc.label
                  << ": wrote " << pdfPath << "\n";

    } // end per-energy loop

    // =========================================================================
    // SUMMARY PDF — Analysis/Output/Summary/quality_summary.pdf
    //
    // Page 1: Four colour-matrix heat maps (channel × energy):
    //   Hit efficiency [%] · Pedestal RMS [mV] ·
    //   Mean HG amplitude [mV] · Timing resolution [ps]
    //
    // Page 2: Trend plots — each metric vs beam energy, one line per channel.
    // =========================================================================
    if (nRunsDone > 0)
    {
        const TString sumPDF = "Analysis/Output/Summary/quality_summary.pdf";

        // Build axis labels
        const char* chNames[8];
        for (int i = 0; i < 8; ++i) chNames[i] = kCap[i].name;

        // ── Helper: fill a TH2F matrix from a 2D array ─────────────────────
        auto fillMatrix = [&](TH2F* h) {
            for (int ir = 0; ir < nRunsDone; ++ir)
                for (int ic = 0; ic < 8; ++ic)
                    h->GetXaxis()->SetBinLabel(ir + 1,
                        Form("%.0f", vEnergy[ir]));
            for (int ic = 0; ic < 8; ++ic)
                h->GetYaxis()->SetBinLabel(ic + 1, chNames[ic]);
            h->GetXaxis()->SetTitle("Beam energy (GeV)");
            h->GetXaxis()->SetTitleSize(0.055);
            h->GetYaxis()->SetTitleSize(0.000); // labels are enough
            h->GetXaxis()->SetLabelSize(0.065);
            h->GetYaxis()->SetLabelSize(0.065);
        };

        // ── Helper: paint each cell value by hand.  ROOT's "COLZ TEXT" painter
        //    renders the literal gStyle paint-format ("%.0f") instead of the
        //    value on this build, so we overlay the numbers ourselves — white on
        //    the dark high-value cells, black on the light low-value cells
        //    (inverted-kRust: high = dark).  zlo/zhi = the colour-scale range.
        auto overlayCells = [&](TH2F* h, const char* fmt, double zlo, double zhi) {
            const double zmid = 0.5 * (zlo + zhi);
            TLatex cell; cell.SetTextAlign(22); cell.SetTextFont(42); cell.SetTextSize(0.034);
            for (int ix = 1; ix <= h->GetNbinsX(); ++ix)
                for (int iy = 1; iy <= h->GetNbinsY(); ++iy) {
                    const double v = h->GetBinContent(ix, iy);
                    if (v <= 0.) continue;             // skip missing/failed cells
                    cell.SetTextColor(v > zmid ? kWhite : kBlack);
                    cell.DrawLatex(h->GetXaxis()->GetBinCenter(ix),
                                   h->GetYaxis()->GetBinCenter(iy), Form(fmt, v));
                }
        };

        // ── Book the four metric matrices ──────────────────────────────────
        TH2F* mEff = new TH2F("mEff", ";Beam energy (GeV);",
                               nRunsDone, 0., (double)nRunsDone,
                               8, 0., 8.);
        TH2F* mPed = new TH2F("mPed", ";Beam energy (GeV);",
                               nRunsDone, 0., (double)nRunsDone,
                               8, 0., 8.);
        TH2F* mHG  = new TH2F("mHG",  ";Beam energy (GeV);",
                               nRunsDone, 0., (double)nRunsDone,
                               8, 0., 8.);
        TH2F* mSig = new TH2F("mSig", ";Beam energy (GeV);",
                               nRunsDone, 0., (double)nRunsDone,
                               8, 0., 8.);

        fillMatrix(mEff); fillMatrix(mPed);
        fillMatrix(mHG);  fillMatrix(mSig);

        for (int ir = 0; ir < nRunsDone; ++ir) {
            for (int ic = 0; ic < 8; ++ic) {
                mEff->SetBinContent(ir + 1, ic + 1, sumHitEff [ir][ic]);
                mPed->SetBinContent(ir + 1, ic + 1, sumPedRMS [ir][ic]);
                mHG ->SetBinContent(ir + 1, ic + 1, sumMeanHG [ir][ic]);
                mSig->SetBinContent(ir + 1, ic + 1, sumTimeSig[ir][ic]);
            }
        }

        // ── Page 1: Heat map matrices ───────────────────────────────────────
        //  Cell values are painted by overlayCells() (ROOT's "COLZ TEXT" prints
        //  the literal "%.0f" on this build — see the helper above).
        {
            TCanvas cS1("cQsum1", "", 2400, 1600);
            cS1.Divide(2, 2, 0.01, 0.01);

            // Pad 1: Hit efficiency
            cS1.cd(1);
            gPad->SetRightMargin(0.14); gPad->SetLeftMargin(0.13);
            gPad->SetTopMargin(0.12);   gPad->SetBottomMargin(0.15);
            mEff->SetMinimum(0.); mEff->SetMaximum(100.);
            mEff->SetTitle(Form("Hit efficiency (%%)  cut: A_{HG} > %.0f mV",
                                 (double)kHG_minPeak));
            mEff->GetZaxis()->SetLabelSize(0.040);
            mEff->GetXaxis()->SetTitleOffset(1.05);
            mEff->Draw("COLZ");
            overlayCells(mEff, "%.0f", 0., 100.);
            { TLatex t; t.SetNDC(); t.SetTextSize(0.048); t.SetTextAlign(22);
              t.DrawLatex(0.50, 0.95, "Hit efficiency [%]  (higher = better)"); }

            // Pad 2: Pedestal RMS
            cS1.cd(2);
            gPad->SetRightMargin(0.14); gPad->SetLeftMargin(0.13);
            gPad->SetTopMargin(0.12);   gPad->SetBottomMargin(0.15);
            mPed->SetMinimum(0.); mPed->SetMaximum(6.);   // flag anything near the 5 mV floor
            mPed->SetTitle("HG pedestal RMS [mV]  (healthy < 3 mV)");
            mPed->GetZaxis()->SetLabelSize(0.040);
            mPed->GetXaxis()->SetTitleOffset(1.05);
            mPed->Draw("COLZ");
            overlayCells(mPed, "%.1f", 0., 6.);
            { TLatex t; t.SetNDC(); t.SetTextSize(0.048); t.SetTextAlign(22);
              t.DrawLatex(0.50, 0.95, "Pedestal RMS [mV]  (lower = better)"); }

            // Pad 3: Mean HG amplitude (auto colour range -> threshold on its data range)
            cS1.cd(3);
            gPad->SetRightMargin(0.14); gPad->SetLeftMargin(0.13);
            gPad->SetTopMargin(0.12);   gPad->SetBottomMargin(0.15);
            mHG->SetTitle(Form("Mean HG amplitude [mV] for A > %.0f mV",
                                (double)kHG_minPeak));
            mHG->GetZaxis()->SetLabelSize(0.040);
            mHG->GetXaxis()->SetTitleOffset(1.05);
            mHG->Draw("COLZ");
            overlayCells(mHG, "%.0f", mHG->GetMinimum(), mHG->GetMaximum());
            { TLatex t; t.SetNDC(); t.SetTextSize(0.048); t.SetTextAlign(22);
              t.DrawLatex(0.50, 0.95, "Mean HG amplitude [mV]"); }

            // Pad 4: Timing resolution
            cS1.cd(4);
            gPad->SetRightMargin(0.14); gPad->SetLeftMargin(0.13);
            gPad->SetTopMargin(0.12);   gPad->SetBottomMargin(0.15);
            mSig->SetMinimum(0.); mSig->SetMaximum(500.);
            mSig->SetTitle(Form("CFD-20%% timing resolution [ps]  "
                                "r < %.1f mm, A > %.0f mV",
                                kFiducial_r_timing, (double)kHG_minPeak));
            mSig->GetZaxis()->SetLabelSize(0.040);
            mSig->GetXaxis()->SetTitleOffset(1.05);
            mSig->Draw("COLZ");
            overlayCells(mSig, "%.0f", 0., 500.);
            { TLatex t; t.SetNDC(); t.SetTextSize(0.048); t.SetTextAlign(22);
              t.DrawLatex(0.50, 0.95, "Timing resolution [ps]  (lower = better)"); }

            cS1.Print(sumPDF + "(");
        }

        // ── Page 2: Trend plots (metric vs energy, one line per channel) ────
        {
            TCanvas cS2("cQsum2", "", 2400, 1600);
            cS2.Divide(2, 2, 0.01, 0.01);

            // Build per-channel TGraph arrays for each metric
            TGraph* gEff[8]; TGraph* gPedG[8];
            TGraph* gHGG[8]; TGraph* gSigG[8];
            for (int ic = 0; ic < 8; ++ic) {
                gEff [ic] = new TGraph(nRunsDone);
                gPedG[ic] = new TGraph(nRunsDone);
                gHGG [ic] = new TGraph(nRunsDone);
                gSigG[ic] = new TGraph(nRunsDone);
                for (int ir = 0; ir < nRunsDone; ++ir) {
                    gEff [ic]->SetPoint(ir, vEnergy[ir], sumHitEff [ir][ic]);
                    gPedG[ic]->SetPoint(ir, vEnergy[ir], sumPedRMS [ir][ic]);
                    gHGG [ic]->SetPoint(ir, vEnergy[ir], sumMeanHG [ir][ic]);
                    gSigG[ic]->SetPoint(ir, vEnergy[ir], sumTimeSig[ir][ic]);
                }
                for (auto* g : {gEff[ic], gPedG[ic], gHGG[ic], gSigG[ic]}) {
                    g->SetMarkerStyle(20 + ic % 8);
                    g->SetMarkerSize(1.1);
                    g->SetMarkerColor(kRChannelCols[ic]);
                    g->SetLineColor(kRChannelCols[ic]);
                    g->SetLineWidth(2);
                }
            }

            // Helper: draw a metric panel from 8 TGraph pointers
            auto drawTrend = [&](TGraph* gs[], const char* yTitle,
                                  double yMin, double yMax) {
                TMultiGraph* mg = new TMultiGraph();
                for (int ic = 0; ic < 8; ++ic) mg->Add(gs[ic], "PL");
                mg->Draw("A");
                mg->GetXaxis()->SetTitle("Beam energy (GeV)");
                mg->GetYaxis()->SetTitle(yTitle);
                mg->GetXaxis()->SetTitleSize(0.050);
                mg->GetYaxis()->SetTitleSize(0.050);
                mg->GetXaxis()->SetLabelSize(0.044);
                mg->GetYaxis()->SetLabelSize(0.044);
                if (yMin < yMax) { mg->SetMinimum(yMin); mg->SetMaximum(yMax); }
                // Legend
                TLegend* leg = new TLegend(0.70, 0.12, 0.92, 0.58);
                leg->SetBorderSize(0); leg->SetTextSize(0.036);
                for (int ic = 0; ic < 8; ++ic)
                    leg->AddEntry(gs[ic], kCap[ic].name, "lp");
                leg->Draw();
            };

            cS2.cd(1);
            gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
            gPad->SetRightMargin(0.05); gPad->SetTopMargin(0.10);
            gPad->SetTickx(1); gPad->SetTicky(1);
            drawTrend(gEff, "Hit efficiency (%)", 0., 105.);
            { TLatex t; t.SetNDC(); t.SetTextSize(0.050); t.SetTextAlign(22);
              t.DrawLatex(0.50, 0.94, "Hit efficiency vs energy"); }

            cS2.cd(2);
            gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
            gPad->SetRightMargin(0.05); gPad->SetTopMargin(0.10);
            gPad->SetTickx(1); gPad->SetTicky(1);
            drawTrend(gSigG, "#sigma_{t}  (ps)", 0., 0.);  // auto y-range
            { TLatex t; t.SetNDC(); t.SetTextSize(0.050); t.SetTextAlign(22);
              t.DrawLatex(0.50, 0.94, "Timing resolution (#sigma_{t}) vs energy"); }

            cS2.cd(3);
            gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
            gPad->SetRightMargin(0.05); gPad->SetTopMargin(0.10);
            gPad->SetTickx(1); gPad->SetTicky(1);
            drawTrend(gHGG, "Mean HG amplitude (mV)", 0., 0.);
            { TLatex t; t.SetNDC(); t.SetTextSize(0.050); t.SetTextAlign(22);
              t.DrawLatex(0.50, 0.94, "Mean HG amplitude vs energy"); }

            cS2.cd(4);
            gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
            gPad->SetRightMargin(0.05); gPad->SetTopMargin(0.10);
            gPad->SetTickx(1); gPad->SetTicky(1);
            drawTrend(gPedG, "Mean pedestal RMS (mV)", 0., 0.);
            // Draw a reference line at 3 mV (healthy threshold)
            { TLine* lH = new TLine(vEnergy.front(), 3., vEnergy.back(), 3.);
              lH->SetLineColor(kOrange+1); lH->SetLineStyle(2);
              lH->SetLineWidth(2); lH->Draw("SAME"); }
            { TLatex t; t.SetNDC(); t.SetTextSize(0.042); t.SetTextColor(kOrange+1);
              t.DrawLatex(0.14, 0.40, "3 mV healthy threshold"); }
            { TLatex t; t.SetNDC(); t.SetTextSize(0.050); t.SetTextAlign(22);
              t.DrawLatex(0.50, 0.94, "Pedestal RMS vs energy"); }

            cS2.Print(sumPDF + ")");
        }

        // Per-energy containment summary to stdout
        std::cout << "\n[qualityPlots] Summary: " << sumPDF << "\n";
        std::cout << "  Containment efficiency per energy:\n";
        for (int ir = 0; ir < nRunsDone; ++ir)
            std::cout << Form("    %.0f GeV:  %.1f%%\n",
                              vEnergy[ir], sumContain[ir]);
    }

    std::cout << "\n[qualityPlots] Done.\n"
              << "  Per-energy:  Analysis/Output/<label>/quality_report.pdf\n"
              << "  Summary:     Analysis/Output/Summary/quality_summary.pdf\n";
}

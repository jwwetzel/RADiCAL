#!/usr/bin/env bash
# =============================================================================
# runAll.sh — run the full RADiCAL analysis chain
# =============================================================================
#
# Run from the repository root:
#   bash Analysis/runAll.sh
#
# Step 1 runs all six energies in PARALLEL — each energy is an independent
# ROOT process, so all available CPU cores are used simultaneously.
# On an M4 Max (14 performance cores) with six ~2 GB input files, expect
# Step 1 to take roughly as long as the single slowest run (~5–10 min).
#
# Step 2 (resolution analysis) runs after all Step 1 jobs finish.
#
# Outputs:
#   Analysis/Output/25GeV/     ntuple.root + per-energy PDFs
#   Analysis/Output/50GeV/
#   Analysis/Output/75GeV/
#   Analysis/Output/100GeV/
#   Analysis/Output/125GeV/
#   Analysis/Output/150GeV/    (four runs chained: RUN1258–1261, ~120k events)
#   Analysis/Output/Summary/   resolution_summary.pdf  timing_per_capillary.pdf
#                               summary.root
# =============================================================================

set -euo pipefail

cd "$(dirname "$0")/.."

# ---------------------------------------------------------------------------
# Ensure Homebrew tools are reachable.  A non-interactive `bash runAll.sh`
# does NOT source the user's interactive rc files, so /opt/homebrew/bin
# (pdftoppm, gs, …) may be missing from PATH even when installed — this is
# what made Step 17 (makeReport) fail with "No PDF->PNG converter found".
# makeReport.py also probes these dirs directly, but exporting them here keeps
# every tool in the chain reachable regardless of how the script is launched.
# ---------------------------------------------------------------------------
for d in /opt/homebrew/bin /usr/local/bin; do
    [[ -d "$d" && ":$PATH:" != *":$d:"* ]] && export PATH="$d:$PATH"
done

# ---------------------------------------------------------------------------
# Ensure ROOT is in PATH (sources thisroot.sh if needed)
# ---------------------------------------------------------------------------
if ! command -v root &>/dev/null; then
    THISROOT_CANDIDATES=(
        "/opt/homebrew/opt/root/bin/thisroot.sh"
        "/usr/local/opt/root/bin/thisroot.sh"
        "/opt/homebrew/bin/thisroot.sh"
    )
    SOURCED=false
    for f in "${THISROOT_CANDIDATES[@]}"; do
        if [[ -f "$f" ]]; then
            # shellcheck disable=SC1090
            source "$f"; SOURCED=true; break
        fi
    done
    if ! $SOURCED; then
        echo "[runAll] ERROR: 'root' not found. Source thisroot.sh or run setup.sh first."
        exit 1
    fi
fi

ROOT_CMD="root -l -b -q"

echo "================================================================"
echo " RADiCAL analysis — $(date)"
echo " ROOT $(root --version 2>&1 | head -1)"
echo "================================================================"
echo ""
# SKIP_PROCESS=1 skips raw->ntuple processing and runs the analysis on
# pre-existing Analysis/Output/<E>/ntuple.root.  Used by the Argon HPC pipeline,
# where processRun is run per-run as a parallel array job and merged separately.
if [[ -n "$SKIP_PROCESS" ]]; then
    echo "Step 1: SKIPPED (SKIP_PROCESS set) — using existing Analysis/Output/<E>/ntuple.root"
    echo ""
else

echo "Step 1: Processing all energies in parallel"
echo "  Each energy is an independent process — all cores used at once."
echo ""

mkdir -p Analysis/Output

# ---------------------------------------------------------------------------
# Pre-compile processRun.C once before launching parallel jobs.
#
# Without this, all 6 background jobs start simultaneously and each tries
# to compile processRun_C.so via ACLiC.  The first job to reach the linker
# creates processRun_C_ACLiC_dict_rdict.pcm; every subsequent job sees the
# file already exists and the rootcling step fails, killing that energy.
#
# Pre-compiling here creates a fresh .so before any parallel work begins.
# The 6 parallel jobs then find the .so up-to-date and skip recompilation.
# ---------------------------------------------------------------------------
echo "  Pre-compiling processRun.C (avoids ACLiC race condition)..."
root -l -b -q -e '.L Analysis/processRun.C+' > Analysis/Output/precompile.log 2>&1
if [[ $? -ne 0 ]]; then
    # Some ROOT versions exit non-zero when .L finds no default function —
    # that is fine as long as the .so was actually created.
    if ls Analysis/processRun_C.so &>/dev/null || \
       ls Analysis/processRun_C*.so &>/dev/null 2>/dev/null; then
        echo "  Pre-compile: .so created (non-zero exit ignored — no default fn)"
    else
        echo "  Pre-compile FAILED — see Analysis/Output/precompile.log"
        exit 1
    fi
fi
echo "  Pre-compile done."
echo ""

# ---------------------------------------------------------------------------
# Launch all processRun jobs in the background (parallel)
# Each writes to its own Analysis/Output/<label>/ directory — no conflicts.
# ---------------------------------------------------------------------------
declare -a PIDS
declare -a LABELS

run_energy() {
    local files="$1" energy="$2" label="$3"
    echo "  [start] $label"
    $ROOT_CMD "Analysis/processRun.C+(\"$files\", $energy, \"$label\")" \
        > "Analysis/Output/${label}.log" 2>&1
    echo "  [done]  $label"
}

# 25 GeV
run_energy "Data/RUN1211_25_GeV.root"   25.  "25GeV"  &  PIDS+=($!); LABELS+=("25GeV")
# 50 GeV
run_energy "Data/RUN1148_50_GeV.root"   50.  "50GeV"  &  PIDS+=($!); LABELS+=("50GeV")
# 75 GeV
run_energy "Data/RUN1112_75_GeV.root"   75.  "75GeV"  &  PIDS+=($!); LABELS+=("75GeV")
# 100 GeV
run_energy "Data/RUN1075_100_GeV.root" 100. "100GeV"  &  PIDS+=($!); LABELS+=("100GeV")
# 125 GeV
run_energy "Data/RUN1034_125_GeV.root" 125. "125GeV"  &  PIDS+=($!); LABELS+=("125GeV")
# 150 GeV — four runs chained
run_energy "Data/RUN1258_150_GeV.root;Data/RUN1259_150_GeV.root;Data/RUN1260_150_GeV.root;Data/RUN1261_150_GeV.root" \
           150. "150GeV"  &  PIDS+=($!); LABELS+=("150GeV")

# Wait for all background jobs; report any failures
FAILED=0
for i in "${!PIDS[@]}"; do
    pid=${PIDS[$i]}; label=${LABELS[$i]}
    if wait "$pid"; then
        echo "  ✓ ${label}"
    else
        echo "  ✗ ${label} FAILED — see Analysis/Output/${label}.log"
        FAILED=$((FAILED+1))
    fi
done

if [[ $FAILED -gt 0 ]]; then
    echo ""
    echo "ERROR: $FAILED job(s) failed. Check the .log files above before proceeding."
    exit 1
fi

# Clean up per-energy log files on success
rm -f Analysis/Output/*.log

fi   # end SKIP_PROCESS guard

echo ""
echo "Step 2: Multi-energy resolution analysis"
echo ""

$ROOT_CMD 'Analysis/analyzeResolution.C+'

echo ""
echo "Step 3: Combined timing resolution (Up+Down corner averaging)"
echo "  Methods: best single · 4-corner (U+D)/2 · avg-8 · A²-wgt-8 · A²-wgt-corners"
echo ""

$ROOT_CMD 'Analysis/timingResolution.C+'

echo ""
echo "Step 4: Timing method comparison"
echo "  Methods: CFD 10/20/30/50% · LED · TOT-corrected LED · 1/A walk corr. · HG/LG ratio corr."
echo ""

$ROOT_CMD 'Analysis/timingMethods.C+'

echo ""
echo "Step 5: Data quality, channel performance, and cut optimisation plots"
echo "  Beam hit maps · MCP/HG/LG amplitude distributions · time walk · TOT"
echo "  Fiducial radius scan · HG amplitude threshold scan"
echo ""

$ROOT_CMD 'Analysis/qualityPlots.C+'

echo ""
echo "Step 6: Energy-binned timing resolution (paper arXiv:2401.01747 method)"
echo "  (DW-UP)/2 CFD-20% · (DW+UP)/2 CFD-20% · (DW-UP)/2 M7-corrected"
echo "  Best estimator from top-3 energy bins · fit σ = a/√E ⊕ b"
echo ""

$ROOT_CMD 'Analysis/timingEnergyBins.C+'

echo ""
echo "Step 7: Cross-energy collated plots (all 6 energies side-by-side)"
echo "  Beam quality · shower containment · channel performance · timing methods"
echo ""

$ROOT_CMD 'Analysis/compareEnergies.C+'

echo ""
echo "Step 8: PbGlass dual-band population analysis"
echo "  Classifies events into EM / hadronic / beam-halo / off-axis populations"
echo ""

$ROOT_CMD 'Analysis/investigatePbGlass.C+'

echo ""
echo "Step 9: MCP timing reference jitter decomposition"
echo "  σ(MCP1-MCP2)/√2 → σ_MCP_single; subtract from per-channel timing"
echo ""

$ROOT_CMD 'Analysis/mcpJitter.C+'

echo ""
echo "Step 10: Containment cut optimisation scan"
echo "  Scans sum_pb/sum_lg threshold 0.10→0.50; measures σ_t vs yield tradeoff"
echo ""

$ROOT_CMD 'Analysis/timingContainmentScan.C+'

echo ""
echo "Step 11: Channel combination brute-force scan (all 255 subsets)"
echo "  Finds the optimal N-channel combination; tests SW-Up (MCP2 ref) impact"
echo ""

$ROOT_CMD 'Analysis/channelCombinationScan.C+'

echo ""
echo "Step 11b: Charge-sharing position reconstruction + timing correction (Layer 5)"
echo "  Asymmetry position from corner LG amplitudes; position-binned timing"
echo "  correction validated split-half (out-of-sample)"
echo ""

$ROOT_CMD 'Analysis/positionCorrection.C+'

echo ""
echo "Step 12: DRS4 time-base verification (Layer 1)"
echo "  Confirms cell-width calibration state; validates stop-cell timing"
echo "  correction out-of-sample (reads raw data; uses highest-energy run)"
echo ""

$ROOT_CMD 'Analysis/drs4TimeBase.C+'

echo ""
echo "Step 12b: Wire-chamber spatial resolution (Layer 1)"
echo "  Delay-line self-calibration (sum-of-end-times); data-driven sigma_x/y"
echo "  for track binning (reads raw waveforms; uses highest-energy run)"
echo ""

$ROOT_CMD 'Analysis/wireChamberResolution.C+'

echo ""
echo "Step 12c: Per-channel transverse maps, before/after 1x1 trigger (Layer 1)"
echo "  Mean amplitude vs beam track (x,y) for each capillary + the 1x1 trigger,"
echo "  good-track vs 1x1-trigger selection (reads raw waveforms; all energies)"
echo ""

$ROOT_CMD 'Analysis/transverseMaps.C+'

echo ""
echo "Step 12d: Shashlik module centre from edges (Layer 1)"
echo "  X/Y edge profiles of mean sum_lg; centre = edge midpoint (beam-independent)"
echo ""

$ROOT_CMD 'Analysis/moduleCenter.C+'

echo ""
echo "Step 13: DRS4 hardware diagnostics (Layer 1)"
echo "  Cell-level noise floor, saturation fractions, spike rates"
echo ""

$ROOT_CMD 'Analysis/drs4Diagnostics.C+'

echo ""
echo "Step 14: Channel integrity report (Layer 1)"
echo "  Active fractions, HG/LG ratios, correlation matrix"
echo ""

$ROOT_CMD 'Analysis/channelIntegrity.C+'

echo ""
echo "Step 14b: Average waveform profiles (Layer 1)"
echo "  Mean +/- RMS HG pulse shapes per capillary (reads raw waveforms)"
echo ""

$ROOT_CMD 'Analysis/averageWaveforms.C+'

echo ""
echo "Step 14c: HG charge / amplitude profiles (Layer 1)"
echo "  Mean HG amplitude vs energy + position maps per capillary"
echo ""

$ROOT_CMD 'Analysis/chargeProfiles.C+'

echo ""
echo "Step 14d: Layer 1 hero figures (Ledovskoy-clean summary)"
echo "  layer1_pulse_shapes / vitals / linearity / timebase"
echo ""

$ROOT_CMD 'Analysis/layer1Summary.C+'

echo ""
echo "Step 14e: Layer 2 hero figures (Ledovskoy-clean summary)"
echo "  layer2_mcp_jitter / sub_mcp / per_channel"
echo ""

$ROOT_CMD 'Analysis/layer2Summary.C+'

echo ""
echo "Step 14f: Layer 3 hero figures (Ledovskoy-clean summary)"
echo "  layer3_beam_map / containment / quality"
echo ""

$ROOT_CMD 'Analysis/layer3Summary.C+'

echo ""
echo "Step 14g: Layer 4 hero figures (Ledovskoy-clean summary)"
echo "  layer4_ladder / walk"
echo ""

$ROOT_CMD 'Analysis/layer4Summary.C+'

echo ""
echo "Step 14h: CFD-fraction trend (Down vs Up groups, all energies)"
echo "  layer4_cfd_fraction_down / _up  — the timing-shoulder optimisation figure"
echo ""

$ROOT_CMD 'Analysis/elbowFractionTrend.C+'

echo ""
echo "Step 14i: CFD-fraction MECHANISM (mean edge + edge-shape jitter)"
echo "  layer4_edge_shape / layer4_edge_jitter  — the cause behind the trend"
echo ""

$ROOT_CMD 'Analysis/edgeMechanism.C+'

echo ""
echo "Step 15: Spatial uniformity scan (Layer 5)"
echo "  A2-weighted combo sigma_t vs beam position (x,y)"
echo ""

$ROOT_CMD 'Analysis/uniformityScan.C+'

echo ""
echo "Step 15b: Layer 5 hero figures (Ledovskoy-clean summary)"
echo "  layer5_timing / energy / uniformity  (needs uniformity_scan.root from Step 15)"
echo ""

$ROOT_CMD 'Analysis/layer5Summary.C+'

echo ""
echo "Step 16: Systematic uncertainty evaluation (Layer 6)"
echo "  Cut variations: fiducial, containment, MCP window, HG threshold"
echo ""

$ROOT_CMD 'Analysis/systematicUncertainties.C+'

echo ""
echo "Step 16b: Layer 6 hero figures (Ledovskoy-clean summary)"
echo "  layer6_budget / band"
echo ""

$ROOT_CMD 'Analysis/layer6Summary.C+'

echo ""
echo "Step 17: Harvest results -> Output/Summary/results.json"
echo "  Single source of truth for every number quoted in the report"
echo ""

$ROOT_CMD 'Analysis/harvestResults.C+'

echo ""
echo "Step 18: Generate interactive web report"
echo "  Converts all PDFs to PNG and assembles Analysis/Output/report.html"
echo ""

python3 Analysis/makeReport.py

echo ""
echo "================================================================"
echo " Done!  Key outputs:"
echo ""
echo " Per-energy (Analysis/Output/<energy>/):"
echo "   ntuple.root                 — analysis ntuple (tree 'rad')"
echo "   quality_report.pdf          — beam quality, cuts, MCP, PbGlass"
echo "   capillary_maps.pdf          — HG amplitude maps"
echo "   timing_distributions.pdf    — per-channel CFD-20% distributions"
echo "   combined_timing.pdf         — combination method distributions"
echo "   timing_methods.pdf          — all 8 timing methods"
echo "   timing_energy_bins.pdf      — energy-binned σ_t distributions"
echo "   energy_distribution.pdf     — ΣLG spectrum"
echo ""
echo " Cross-energy summary (Analysis/Output/Summary/):"
echo "   cross_energy_beam_quality.pdf      — beam maps + MCP + ΣLG trends"
echo "   cross_energy_containment.pdf       — PbGlass scatter + efficiency"
echo "   cross_energy_channels.pdf          — per-channel σ_t + active fraction"
echo "   cross_energy_timing_methods.pdf    — combination methods + σ=a/√E fit"
echo "   pbglass_investigation.pdf          — 4-population scatter analysis"
echo "   mcp_jitter.pdf                     — MCP1-MCP2 jitter, σ_MCP≈72 ps"
echo "   timing_containment_scan.pdf        — containment cut optimisation"
echo "   channel_combination_scan.pdf       — all-255-subset channel combo scan"
echo "   drs4_diagnostics.pdf               — cell-level noise, saturation, spike rates"
echo "   channel_integrity.pdf              — active fractions, HG/LG ratios, correlation matrix"
echo "   uniformity_scan.pdf                — A2-weighted sigma_t vs beam position (x,y)"
echo "   systematic_uncertainties.pdf       — cut-variation systematic uncertainty table"
echo "   resolution_summary.pdf / *.root    — all summary ROOT graphs"
echo ""
echo " Interactive report:"
echo "   Analysis/Output/report.html        — open in any browser"
echo "================================================================"

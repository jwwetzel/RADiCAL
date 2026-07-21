#!/usr/bin/env bash
# tools/run_gate.sh — safety wrapper for running a published gate macro.
# (REORG_PLAN B3: the gate macros' bytes stay untouched; the guard lives here.)
#
# ⚠ GATE MACROS REGENERATE THEIR COMMITTED OUTPUTS IN PLACE. Running one is a
# deliberate act of verification: afterwards `git diff` must be clean (byte-identical
# outputs) — a dirty diff is a FINDING, not a formality (ANALYSIS_GUIDE.md). The one
# known benign residual is thesis_postfix.pdf's embedded /CreationDate (see repro.sh,
# which automates the full check-and-restore protocol; prefer it for routine use).
#
# Usage: tools/run_gate.sh <gate> --yes     (from anywhere; runs from the repo root)
#        tools/run_gate.sh --list
# Refuses to run without data/2023/reduced/ (except apparatus_composite, tex-only)
# and without the explicit --yes acknowledgement of the clobber warning above.
# (bash-3.2-compatible on purpose: macOS ships bash 3.2, no associative arrays.)
set -u
cd "$(dirname "$0")/.." || exit 2

GATES="apparatus_composite depth_dial full_fiducial_check method_gain_postfix mixed_corner_map mixed_killshot_bootstrap position_reconciliation satellite_removal systematics_postfix timing_fit_summary"

gate_info() {  # $1 = gate; sets ENTRY, OUTPUTS, NEEDS_DATA
  NEEDS_DATA=yes
  case "$1" in
    timing_fit_summary)
      ENTRY="root -l -b -q 'papers/scripts/timing_fit_summary/timingFitSummary.C'"
      OUTPUTS="papers/tables/timing_fit_summary_2026-06-09.md, papers/figures/timing_fit_summary/timing_fit_summary.png, papers/timing/figs/thesis_postfix.{png,pdf}" ;;
    full_fiducial_check)
      ENTRY="root -l -b -q 'papers/scripts/full_fiducial_check/fullFiducialCheck.C+'"
      OUTPUTS="stdout only (committed log: papers/scripts/full_fiducial_check/full_fiducial_result.log — do NOT overwrite; diff your stdout against it)" ;;
    systematics_postfix)
      ENTRY="root -l -b -q 'papers/scripts/systematics_postfix/systematicsPostfix.C+'"
      OUTPUTS="papers/timing/tab_systematics.tex, papers/tables/systematics_postfix_2026-06-09.md, papers/figures/systematics_postfix/systematics_stability.png, papers/timing/figs/systematics.png" ;;
    method_gain_postfix)
      ENTRY="root -l -b -q 'papers/scripts/method_gain_postfix/methodGainPostfix.C+'"
      OUTPUTS="papers/figures/method_gain_postfix/method_gain_postfix.png, papers/tables/method_gain_postfix_2026-06-09.md, papers/timing/tab_methods.tex" ;;
    depth_dial)
      ENTRY="root -l -b -q 'papers/scripts/depth_dial/depthDial.C+'"
      OUTPUTS="papers/figures/depth_dial/depth_dial.png, papers/figures/depth_dial/depth_dial_diag.png" ;;
    mixed_killshot_bootstrap)
      ENTRY="root -l -b -q 'papers/scripts/mixed_killshot_bootstrap/mixedKillshotBootstrap.C+' && root -l -b -q 'papers/scripts/mixed_killshot_bootstrap/makeMixedKillshotFigure.C+'"
      OUTPUTS="papers/figures/mixed_killshot_bootstrap/{killshot_bootstrap.png, mixed_h2h_corrected.png/.pdf, makeMixedKillshotFigure_CAPTION.txt}" ;;
    mixed_corner_map)
      ENTRY="root -l -b -q 'papers/scripts/mixed_corner_map/cornerDiscriminant.C+'"
      OUTPUTS="papers/figures/mixed_corner_map/corner_discriminant.png" ;;
    position_reconciliation)
      ENTRY="root -l -b -q 'papers/scripts/position_reconciliation/positionReconcile.C+'"
      OUTPUTS="papers/figures/position_reconciliation/position_reconcile.png" ;;
    satellite_removal)
      ENTRY="root -l -b -q 'papers/scripts/satellite_removal/satelliteRemoval.C+'"
      OUTPUTS="papers/figures/satellite_removal/satellite_removal.{png,pdf} + caption sidecar" ;;
    apparatus_composite)
      ENTRY="cd papers/scripts/apparatus_composite && tectonic apparatus_composite.tex"
      OUTPUTS="papers/scripts/apparatus_composite/apparatus_composite.pdf (untracked build artifact; canonical copy lives in papers/timing/figs/)"
      NEEDS_DATA=no ;;
    *) return 1 ;;
  esac
}

if [[ "${1:-}" == "--list" || -z "${1:-}" ]]; then
  echo "Gates (see papers/scripts/INDEX.md for the full runbook):"
  for g in $GATES; do echo "  $g"; done
  echo "Usage: tools/run_gate.sh <gate> --yes"
  exit 0
fi

gate="$1"
gate_info "$gate" || { echo "ERROR: unknown gate '$gate' (tools/run_gate.sh --list)"; exit 2; }

if [[ "$NEEDS_DATA" == yes && ! -d data/2023/reduced/DSB1 ]]; then
  echo "REFUSING: data/2023/reduced/ is not present — the gate would either crash or"
  echo "regenerate committed evidence from nothing. Fetch the data first (data/2023/README.md,"
  echo "verify with: shasum -a 256 -c data/2023/MANIFEST.sha256)."
  exit 2
fi

if [[ "${2:-}" != "--yes" ]]; then
  echo "⚠ '$gate' REGENERATES COMMITTED OUTPUTS IN PLACE:"
  echo "   $OUTPUTS"
  echo "This is correct ONLY as a verification act: run on a clean tree, then 'git diff'"
  echo "must come back clean (dirty = finding). Re-run with:  tools/run_gate.sh $gate --yes"
  exit 2
fi

# shellcheck disable=SC1091
source setup.sh >/dev/null
echo "Running gate '$gate' from $(pwd)..."
eval "$ENTRY"
rc=$?
echo
echo "Gate exited $rc. Now check:  git status --porcelain --untracked-files=no"
echo "(expected: clean, or thesis_postfix.pdf timestamp-only — see tools/repro.sh)"
exit $rc

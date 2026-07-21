#!/usr/bin/env bash
# tools/repro.sh — one-command reproduction of the gated headline chain.
# (CODE_AUDIT_2026-07-21 deep item D1: collapses the manual reproduction protocol
#  to one command, with the known thesis_postfix.pdf /CreationDate exemption built in.)
#
# What it runs, from the repo root:
#   0. optional: shasum-verify data/2023/MANIFEST.sha256 (--check-data)
#   1. reduce/verify.C                     — data presence / staleness gate (stdout only)
#   2. papers/scripts/timing_fit_summary/timingFitSummary.C
#        regenerates its committed outputs IN PLACE (table md + 3 figures) — the
#        reproduction test is that git sees no change.
#   3. papers/scripts/full_fiducial_check/fullFiducialCheck.C
#        stdout-only; diffed against the committed full_fiducial_result.log.
#   4. verdict: `git status` must be clean EXCEPT papers/timing/figs/thesis_postfix.pdf,
#        which legitimately differs ONLY in its embedded PDF /CreationDate + /ModDate
#        (ROOT's PDF writer embeds the wall clock — environment drift, not value drift;
#        characterized byte-exactly in CODE_AUDIT_2026-07-21). If that is the only
#        difference, the file is restored from HEAD and the run PASSES.
#
# NOT covered here: the tectonic manuscript/apparatus builds (need tectonic; see
# REPRODUCE.md) and the remaining gates (see papers/scripts/INDEX.md for the full
# runbook). Requires data/2023/reduced/ (~13 GB, not in git) and ROOT.
#
# Usage: tools/repro.sh [--check-data]
# Exit:  0 = PASS (all outputs reproduced), 1 = FAIL, 2 = environment problem.
set -u

cd "$(dirname "$0")/.." || exit 2
command -v root >/dev/null 2>&1 || { echo "ERROR: ROOT not on PATH."; exit 2; }
[[ -d data/2023/reduced/DSB1 ]] || { echo "ERROR: data/2023/reduced/ missing — see data/2023/README.md."; exit 2; }

if ! git diff --quiet || ! git diff --cached --quiet; then
  echo "ERROR: working tree not clean. Commit or stash first — this script uses"
  echo "'git diff' as the reproduction verdict, so it needs a clean baseline."
  exit 2
fi

# shellcheck disable=SC1091
source setup.sh >/dev/null

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
overall=PASS
t0=$SECONDS

step() { echo; echo "== $* =="; }

if [[ "${1:-}" == "--check-data" ]]; then
  step "STEP 0: data manifest (shasum -a 256 -c)"
  if shasum -a 256 -c data/2023/MANIFEST.sha256 2>/dev/null | grep -v ': OK$' | grep -q .; then
    shasum -a 256 -c data/2023/MANIFEST.sha256 2>/dev/null | grep -v ': OK$'
    echo "STEP 0: FAIL (data drift)"; overall=FAIL
  else
    echo "STEP 0: PASS (29 files OK)"
  fi
fi

step "STEP 1: reduce/verify.C (data presence gate)"
if root -l -b -q 'reduce/verify.C+("data/2023/reduced")' > "$TMP/verify.log" 2>&1 \
   && ! grep -q 'STALE' "$TMP/verify.log"; then
  echo "STEP 1: PASS"
else
  tail -20 "$TMP/verify.log"; echo "STEP 1: FAIL"; overall=FAIL
fi

step "STEP 2: timing_fit_summary gate (headline 25.7 +- 0.6)"
if root -l -b -q 'papers/scripts/timing_fit_summary/timingFitSummary.C' > "$TMP/tfs.log" 2>&1; then
  echo "exit 0 ($((SECONDS-t0)) s elapsed total)"
else
  tail -20 "$TMP/tfs.log"; echo "STEP 2: FAIL (macro error)"; overall=FAIL
fi

step "STEP 3: full_fiducial_check gate (~50 ps companion)"
if root -l -b -q 'papers/scripts/full_fiducial_check/fullFiducialCheck.C+' > "$TMP/ffc_raw.log" 2>&1; then
  # Committed log = stdout capture; strip ROOT/ACLiC banner noise from both sides
  grep -v '^Processing\|^Info in\|^Warning in' "$TMP/ffc_raw.log" > "$TMP/ffc.log"
  grep -v '^Processing\|^Info in\|^Warning in' papers/scripts/full_fiducial_check/full_fiducial_result.log > "$TMP/ffc_committed.log"
  if diff -u "$TMP/ffc_committed.log" "$TMP/ffc.log" > "$TMP/ffc.diff"; then
    echo "STEP 3: PASS (stdout identical to committed log)"
  else
    head -30 "$TMP/ffc.diff"; echo "STEP 3: FAIL (numeric drift vs committed log)"; overall=FAIL
  fi
else
  tail -20 "$TMP/ffc_raw.log"; echo "STEP 3: FAIL (macro error)"; overall=FAIL
fi

step "STEP 4: git verdict (byte-identity of committed outputs)"
PDF="papers/timing/figs/thesis_postfix.pdf"
dirty="$(git status --porcelain --untracked-files=no)"
if [[ -z "$dirty" ]]; then
  echo "STEP 4: PASS (tree fully clean — even the PDF timestamp matched)"
elif [[ "$dirty" == " M $PDF" ]]; then
  # The one allowed residual. Verify it is ONLY /CreationDate + /ModDate metadata:
  git show "HEAD:$PDF" > "$TMP/head.pdf"
  scrub() { perl -pe 's{/(Creation|Mod)Date \(D:[^)]*\)}{/${1}Date (D:SCRUBBED)}g' "$1" > "$2"; }
  scrub "$TMP/head.pdf" "$TMP/head_scrubbed.pdf"
  scrub "$PDF"          "$TMP/tree_scrubbed.pdf"
  if cmp -s "$TMP/head_scrubbed.pdf" "$TMP/tree_scrubbed.pdf"; then
    git checkout -- "$PDF"
    echo "STEP 4: PASS ($PDF differed only in /CreationDate+/ModDate; restored from HEAD)"
  else
    echo "STEP 4: FAIL ($PDF differs beyond the timestamp exemption — investigate:"
    echo "         cmp -l $PDF <(git show HEAD:$PDF)   — tree left dirty on purpose)"
    overall=FAIL
  fi
else
  echo "$dirty"
  echo "STEP 4: FAIL (unexpected modified files above — a gated output drifted."
  echo "         Treat as a FINDING per ANALYSIS_GUIDE.md; tree left dirty on purpose)"
  overall=FAIL
fi

echo
echo "======================================"
echo "repro.sh: $overall  ($((SECONDS-t0)) s total)"
echo "======================================"
[[ "$overall" == PASS ]]

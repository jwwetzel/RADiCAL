# RADiCAL Analysis — Gap-Closing Plan

Living document. Tracks the work needed to bring the May-2023 test-beam analysis
to a publishable, industry-leading standard. Two sources feed this list:
 (A) the 6-layer CMS-reviewer assessment done earlier in the project, and
 (B) the comparison against A. Ledovskoy's prior analyses (Jan 2021 simulation;
     July-2021 FNAL test-beam data analysed Nov–Dec 2021).

Status keys: ☐ TODO   ◐ IN PROGRESS   ☑ DONE   ⊘ TESTED-NOT-WORTH-IT

────────────────────────────────────────────────────────────────────────────
## Already closed this project (context)
────────────────────────────────────────────────────────────────────────────
☑ Unified plot style; Nature-review presentation fixes
☑ DRS4 time-base verified (no cell-width calib; stop cell recovered) + stop-cell
   correction for combos (M2/M3), validated split-half, self-protected
☑ Histogram-binning artefact fixed (combo σ 96→78 ps at 150 GeV)
☑ Amplitude-dependent MCP jitter σ_MCP(A); MCP subtraction on combos (Page 3c)
☑ Saturation/spike flags in ntuple; full re-process of all 6 energies
☑ Systematic-uncertainty framework (cut variations) — partial (see #G2)
☑ Plot review: rendering bug + 9 code bugs + headline polish + enhancements
⊘ Cable-delay alignment — measured negligible/harmful, deliberately not applied

────────────────────────────────────────────────────────────────────────────
## Ledovskoy-comparison gaps (ranked) — THIS PLAN
────────────────────────────────────────────────────────────────────────────

### #1 ☑ Asymmetry / charge-sharing position reconstruction + position-binned
###      timing correction                          [DONE — validated]
RESULT (positionCorrection.C, split-half out-of-sample, all 6 energies):
  - VALIDATED: charge-sharing asymmetry tracks wire-chamber position,
    rho = -0.72 (x), +0.76 (y). Our 8-capillary geometry supports it.
  - A^2-combo (MCP-referenced): correction HELPS, largest at low E
    (25 GeV 135.6->121.9 ps, 50 GeV 79.5->72.3 ps), neutral at high E.
    Never hurts. Pattern RMS 53 ps @25 -> ~13-25 ps @high E.
  - Wire-chamber-based correction is WORSE than charge-sharing: the
    calorimeter's own charge sharing is a better predictor of the timing walk
    than the external tracker (sees the actual shower, not just track entry).
  - HEADLINE (DW-UP)/2 MCP-free corner estimator: NO benefit (neutral/slightly
    worse). It already cancels the position walk by the same common-mode
    cancellation that makes it MCP- and cell-width-free. Consistent with the
    whole story: the corner difference is robust by construction.
  Output: Summary/position_correction.{pdf,root}; runAll.sh step.
  OPTIONAL follow-on: integrate the charge-sharing correction into the
  MCP-referenced combos (M2/M3) in timingResolution (guarded, like stop-cell).
  Low urgency: combos are MCP-limited and not the headline.
Ledovskoy (Dec 2021) reconstructs beam impact position from calorimeter charge
sharing — asym = (A_a−A_b)/(A_a+A_b) — and selects a narrow position spot to
remove position-dependent time-walk: 154 ps → 113 ps.  His Jan-2021 sim predicted
this "position correction" is worth ~25 ps.  We MEASURE the non-uniformity
(uniformityScan: σ_t vs x,y) but never CORRECT for it.  Our 8-capillary
(NW/NE/SE/SW × U/D) layout is better suited than his 4 channels.
PLAN:
  - New macro positionCorrection.C
  - Build asym_x, asym_y from the 4 corner LG amplitudes
  - VALIDATE: asym vs wire-chamber x_trk/y_trk correlation (does charge sharing
    track position for our geometry? — empirical, no assumptions)
  - Map combo timing across (asym_x, asym_y); train a 2D position correction on
    EVEN events, apply to ODD, measure σ_t before/after OUT-OF-SAMPLE
  - Compare asym-based vs wire-chamber-based position correction
  - If validated → integrate into timingResolution / timingEnergyBins
ACCEPTANCE: out-of-sample σ_t improvement on the combination estimator.

### #2 ☐ GEANT4 + photostatistics / light-yield model
We have NO simulation and no photostatistics floor.  Ledovskoy's Jan-2021 sim is
the template: GEANT4 depositions → SLitrani optical transport → detector/amp
response → σ_t ∝ 1/√LY, predicting ~30 ps.  Answers "electronics- or
light-limited?" and is required for journal publication.
PLAN: GEANT4 e⁻ at 25–150 GeV into the W/LYSO geometry; compare E spectrum &
sampling term to data; estimate N_pe and the photostatistics timing floor.

### #3 ⊘ Amplifier non-linearity / saturation recovery — NOT APPLICABLE
FINDING: at 150 GeV, 0.00% of fiducial events have ANY saturated HG channel
(hg_saturated never trips; HG stays far below the 950 mV rail).  Energy uses LG
(3000 mV range, never saturates).  There is nothing to recover — Ledovskoy's
>350 mV amplifier non-linearity does not occur in our data.  The template-FIT
idea for timing he himself found "not understood" run-to-run; our CFD + Crystal
Ball already give robust timing.  => Closed with evidence, no code change.

### #4 ☑ Inter-channel calibration  (MIP-absolute: N/A, documented)
FINDING: no clean MIP peak at 150 OR 25 GeV (pure EM shower) -> MIP-based
ABSOLUTE deposited-energy calibration not available from electron-only data.
DONE: channelCalibration.C — inter-channel RELATIVE gains from LG amplitudes.
Equalisation reduces per-event 8-channel spread 14-16% -> ~8% at all energies
(matches Ledovskoy's "equal within 10%"; residual ~8% IS shower-position charge
sharing). 150 GeV gains: SW-U=0.74 (weak, MCP2 ref), SE-D=1.18, others ~1.0.
Output Summary/channel_calibration.{pdf,root}. Absolute scale via the existing
beam-energy linearity (27.5 mV/GeV).

### #5 ⊘ DRS empty-channel common-mode subtraction — TESTED, NEGLIGIBLE
The 9th (trigger) channel ch8 of D0G0 is the common-mode reference: its DC
baseline correlates rho=0.95-0.97 with every signal channel.  BUT our
per-channel pedestal subtraction (samples 3-52) already removes that DC
common-mode.  The PER-SAMPLE residual correlation (after DC removal) is only
rho=0.17, and subtracting ch8 sample-by-sample makes the noise floor slightly
WORSE (1.292 -> 1.300 mV) by adding ch8's own noise.  => Not applied.
(Ledovskoy needed it because his pipeline likely lacked per-channel pedestal
subtraction.)  Same outcome class as cable-delay alignment.

### #6 ◐ Cross-check Ledovskoy's unresolved "mysteries"
  - #6b BIMODAL timing: our SE-D hg_cfd is NOT single-Gaussian (chi2/ndf
    9.6 -> 1.1 with a 2nd component) — but this is the non-Gaussian TAIL we
    already handle with the Crystal Ball fit, not necessarily his SEPARATED
    two-peak bifurcation (which needs a waveform-level rising-edge study).
    Practical impact already covered by CB. ☑ checked.
  - #6a OPTIMAL CFD ≈ 3%: ☑ DONE & CONFIRMED. Added hg_cfd03/hg_cfd05 branches
    (processRun, reprocessed all 6). CFD-3% is much better PER-CHANNEL,
    especially Down channels (SE-D @150 GeV: 342 ps@CFD-20% -> 173 ps@CFD-3%);
    Up channels modest. HEADLINE (DW-UP)/2 4-corner: ~5 ps better @150 GeV
    (61->56), neutral @50 (per-channel gains wash out in the corner average,
    robust-by-construction). Branches available for production adoption.
    RECOMMEND: switch combination-method base CFD fraction to 3% (re-run timing
    chain on hg_cfd03) — modest headline gain, large per-channel/combination gain.

────────────────────────────────────────────────────────────────────────────
## Earlier physics gaps (from the 6-layer CMS assessment)
────────────────────────────────────────────────────────────────────────────
### #G1 ☑ Beam-momentum-spread — RESOLVED (literature value; negligible impact)
Replaced the arbitrary decreasing placeholders with the LITERATURE-cited H2
value: the H2 momentum bite is ~energy-independent, max Delta p/p = +-2%, ~1%
typical (CERN North Area H2 handbook).  Now kBeamSpread = 1.0% at all energies,
cited + flagged "literature-typical, not run-specific".  KEY FINDING: impact is
NEGLIGIBLE — our sigma_E/E is 11-19% (leakage-dominated in the 14 mm module), so
subtracting 1% in quadrature changes the constant term by <0.1% (detector-only
≈ measured).  Run-specific tune-sheet refinement would only matter for a larger,
better-contained module.  No longer "blocked" — resolved with evidence.

### #G2 ☑ Complete the systematic budget (fit-model term added)
DONE: systematicUncertainties.C now fits the nominal distribution with BOTH
Gauss and Crystal Ball and folds |sigma_Gauss - sigma_CB| into the total
systematic (quadrature) — the fit-shape ambiguity from non-Gaussian tails.
Cut-variation terms were already present.  (MCP-subtraction term: the headline
(DW-UP)/2 is MCP-free so it carries none; the MCP-limited combos already show
the σ_MCP sensitivity on mcpJitter Page 3c.)

### #G3 ☑ Energy from integrated charge (branches added)
DONE: hg_charge[8]/lg_charge[8] integrated-charge branches added to processRun
(WaveformUtils ExtractPulse/Multi), reprocessed all 6 energies.  Validated:
LG charge-vs-peak rho=0.996 — for our data (stable pulse shape, 0% saturation)
charge ~ peak, so the energy results are essentially unchanged, but the more
rigorous estimator is now available for the energy sum.

### #G4 ☑ Walk-correction cross-validation (done)
DONE: timingEnergyBins.C now does 5-fold CV of the M7 walk correction (train
slopes on 4/5, evaluate Method C on held-out 1/5, pooled over best bins).
Out-of-sample overfit is small: 0.3-2.6 ps at most energies (5.7 ps at one
low-stat point).  HEADLINE 150 GeV: in-sample 38.7 -> 5-fold OOS 38.9 ps
(overfit 0.3 ps) => the ~38 ps headline is cross-validation-robust.

────────────────────────────────────────────────────────────────────────────
## Notes
────────────────────────────────────────────────────────────────────────────
- Methodology rule (hard-won): ALWAYS validate a correction OUT-OF-SAMPLE
  (train-even / test-odd) before claiming a gain, and never ship a correction
  that can harm the result without a self-protection guard.
- ALWAYS open the regenerated PNG before claiming a plot is fixed.
- Headline result remains σ_t = 38 ps at 150 GeV via the MCP-free (DW−UP)/2
  estimator; combination methods are MCP-reference-limited.

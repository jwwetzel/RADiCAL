# Method-Gain Post-Fix Audit (pre-analysis) — §5.3 recompute

Date: 2026-06-09. Written BEFORE the script.

## 1. Old §5.3 claims and why they are retired
- Old numbers: "improves the resolution by **3.1 ps** at 150 GeV", "lowers the floor from
  **22.1 ± 1.8 to 19.5 ± 1.1 ps**", "a = 205→201", "best 25.4 ps". Source: pre-fix
  `methodCompare.C` runs (also embedded in `papers/timing/tab_methods.tex`: 28.5/28.1/25.3 at
  150 GeV; floors 22.1/20.0/19.5).
- Why retired: computed with the OLD `tebSigma` (fit window set by the tail-inflated 5σ RMS — the
  GATE-3-era estimator bug) and WITHOUT the in-event channel-consistency veto; also no paired
  uncertainty on the gain (panel flagged the missing CI). Currently quarantined behind TODO-P2 in
  the manuscript; prose softened to "a few picoseconds".

## 2. Corrected comparison (pre-registered design)
- Build: **DSB1** (the saturation-recovery showcase: clip fraction 4%→94% across 25→150 GeV);
  energies 25–150 (6 points).
- **Identical events**: an event enters the comparison set only if, for EACH of the three sources
  {srCFD (`hg_lgcfd`), cfd05 (`hg_cfd05`), LED (`hg_led`)}, all 8 timing ends are valid
  (peak ≥ 20 mV, valid crossing) AND the event passes that source's in-event median veto
  (|t−median| < 2 ns) — the intersection set. The brightest-1000 by `sum_lg` is then selected ONCE
  from this intersection, so all three sources are evaluated on the SAME 1000 events at each energy.
- Same cuts as production: wc_ok, MCP1 200–750 mV, TimingFiducialR(E) (continuous), same channels.
- Width: the PRODUCTION post-fix `tebSigma` per source (the published-convention width).
- Gain uncertainty: paired Poisson bootstrap (2000 replicas) of the fixed-window truncated-RMS
  difference Δσ = σ(cfd05) − σ(srCFD) (windows frozen per source from the full selected sample);
  the paired design cancels shower-population fluctuations. Same for LED − srCFD.
- Floors: a/√E ⊕ b fit per source over the identical-event σ(E) points (stat errors σ/√2000,
  PDG √(χ²/ndf) inflation).
- Saturation locality check: at 150 GeV, split the intersection set by the number of saturated ends
  (`hg_saturated` sum: ≥7 "fully clipped" vs ≤5 "less clipped") and report the gain in each — the
  recovery mechanism predicts the gain concentrates in the clipped class. At 25 GeV (4% clipped)
  srCFD ≈ cfd05 is the closure expectation.

## 3. Headline anchor
DSB1 at 150 GeV (most clipped, the published headline energy), with the 25 GeV closure point and
the per-energy Δσ(E) curve.

## 4. Physical question
How much does the saturation-recovered CFD improve identical-event timing over the clipped-peak CFD
(and over LED), is the improvement statistically real (paired CI), is it confined to
clipped/saturated events, and does it change the four-build narrative (expected: no — the method is
an estimator-bias correction, not a detector effect; the four-build table already uses the adopted
sources).

## 5. Pre-registered language
- ALLOWED: "srCFD improves the clipped-pulse timing width by X ps (Y%) relative to cfd05 in an
  identical-event comparison [paired-bootstrap CI]"; "the gain concentrates in saturated events";
  "method choice is a correction of estimator bias, not a new detector effect"; "on unclipped
  low-energy pulses srCFD and cfd05 agree (closure)".
- FORBIDDEN: the old 3.1 ps / 22.1→19.5 / 25.4 unless independently reproduced post-fix; "the
  method changes the physics conclusions".

## Outputs
`methodGainPostfix.C` (this dir) → `papers/figures/method_gain_postfix/method_gain_postfix.png`,
`papers/tables/method_gain_postfix_2026-06-09.md`, regenerated `papers/timing/tab_methods.tex`
(post-fix), §5.3 patch, log `method_gain_postfix_result.txt` (this dir).

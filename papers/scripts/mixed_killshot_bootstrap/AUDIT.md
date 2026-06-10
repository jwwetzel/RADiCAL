# MIXED Kill-Shot Bootstrap Audit (pre-analysis) — GATE 6

Date: 2026-06-09. Written BEFORE the bootstrap script.

## 1. What exactly is the "0.99 ratio"? (source read in full: `analyze/studies/mixedHeadToHead.C`)

- **Numerator / denominator:** mean over the five MIXED energies (50–150 GeV) of
  σ_group(DSB1)/σ_group(LuAG), where σ_group = √(mean of squares) of the PER-CAPILLARY intrinsic
  σ_t values in that material group at that energy. Convention: **DSB1/LuAG** (≈0.99 ⇒ DSB1 corners
  marginally better).
- **Per-capillary σ_t:** obtained from a same-layer pairwise solve: for every same-layer pair (i,j)
  of timing ends, the width σ_ij of the IN-EVENT difference t_i − t_j satisfies σ²_ij = s²_i + s²_j;
  the linear system over all pairs is solved per layer (4 Down ends; 3 Up ends).
- **Timing estimator:** CFD-5% from the LEGACY raw-slot branches `s_cfd05[36]` (slots 0–6), NOT the
  canonical `hg_cfd05[8]`. Slots 0–6 = NW-D, NE-D, SE-D, SW-D, NW-U, NE-U, SE-U.
  **SW-U (slot 7 region) is EXCLUDED — it lives outside DRS0 G0, so its pairwise differences would
  not cancel the group timebase.** The Up layer therefore has only ONE DSB1 end (NE-U).
- **Energies:** 50, 75, 100, 125, 150 GeV, equally weighted in the final mean.
- **Event-paired?** Yes at the pair level: every entry is an in-event difference (arrival time, group
  timebase, MCP all cancel). The widths are computed over events per pair; validity was per-pair
  (not a common event set).
- **Material labels:** **brightness-inferred in-script** (isD = mean amplitude > 0.65·max) — exactly
  the circularity GATE 1 addressed. GATE 1's pulse-shape result (8/8, capillary-pair consistent)
  fixes the map independently: **DSB1 = NE,SW → slots 1,3 (Down) and 5 (Up); LuAG = NW,SE →
  slots 0,2 (Down) and 4,6 (Up).**
- **The ratio is of σ (timing widths)** — not stochastic terms, not means, not amplitudes.
- Fiducial: r < 3.0 mm around the per-file WC centroid; slot validity `s_peak > 20 mV`.
- χ²/ndf = 0.4/5 was also reported for D−L consistency — suspiciously small; the per-point error
  bars (cap-to-cap spread / √n with a 3% floor) are likely conservative. To be re-examined here.

## 2. Pre-registered claim under test

**Claim:** the MIXED module shows no statistically meaningful timing penalty for LuAG:Ce relative to
DSB1 once the comparison is made within the same module / same showers / same DRS group / same
estimator; the large cross-build difference (a = 455 vs 201) is therefore dominated by detected
light yield and build/module effects, not intrinsic WLS re-emission kinetics.

**What this test CANNOT prove (pre-registered):**
- It does not prove LuAG and DSB1 are universally identical (different geometry/light levels could
  re-expose kinetics).
- It does not prove WLS kinetics are irrelevant in all configurations — only that they are not the
  dominant term at THESE light levels in THIS geometry under THESE estimators.
- It does not replace the 4-build energy-scan comparison (the a-term ordering stands on its own).
- It does not remove the obligation to report the light-yield and clipping differences between the
  corner materials (DSB1 corners are ×2–3 brighter and clip; LuAG corners do not).

## 3. Planned measurements (pre-registered)

**M1 — Direct event-paired width comparison (NEW, primary for the CI).** The Down layer contains a
same-material pair of each kind: DSB1 (slots 1,3) and LuAG (slots 0,2). On a COMMON event set (all
seven ends valid, identical cuts), per event compute Δ_D = t₁−t₃ and Δ_L = t₀−t₂; the material
widths are σ(Δ)/√2 (RMS-mean of the two capillaries) and the paired ratio is
R = σ(Δ_D)/σ(Δ_L) = √[(s₁²+s₃²)/(s₀²+s₂²)]. Strictly same events, same layer, pure same-material
in-event differences. Width definition: fixed-window truncated RMS (window from the full-sample
2.5σ core, frozen per energy — part of the estimator definition, enabling Poisson bootstrap).

**M2 — Pairwise-solve reproduction** (the original quantity) with GATE-1 labels replacing the
brightness threshold, on the canonical `hg_*` branches.

**M3 — Event bootstrap:** Poisson(1) weights per event (equivalent to resampling with replacement),
≥1000 replicas, per energy; report median, 68% and 95% CI of the per-energy ratio and of the
5-energy mean ratio (replicas combined energy-wise).

**M4 — Run jackknife:** leave-one-run-out (runs with ≥700 selected events), recompute the affected
energy's ratio and the 5-E mean; report spread and the worst-case excluded run; flag single-run
dependence.

**M5 — Method dependence:** repeat M1 for srCFD (`hg_lgcfd`, primary per current standards), cfd05
(`hg_cfd05`, the original estimator), LED (`hg_led`). Note physical appropriateness: cfd05 walks on
the clipped DSB1 corners; LED is appropriate for the dim LuAG corners; srCFD is the
clip-independent choice. In-pair differencing shares method systematics within each pair, so the
RATIO is expected more robust than absolute widths.

**M6 — Wrong-map diagnostics:** (a) SWAPPED map (NE,SW↔NW,SE): the ratio must invert (R→1/R);
(b) SCRAMBLED map (one bright + one dim per "material": pairs (1,0) vs (3,2) — cross-material
differences): the per-"material" widths must EQUALIZE and the amplitude ratio (the label-sensitivity
control) must collapse from ×2–3 to ≈1. Together these show the assignment carries real information
and the result is not label-insensitive bookkeeping. NOTE (pre-registered honesty): because the
kill-shot CONCLUSION is "σ ratio ≈ 1", the σ ratio alone is weakly map-sensitive by construction;
the amplitude ratio is the sharp label control.

## 4. Outcome definitions (pre-registered)
- **PASS:** the GATE-1-labeled, common-event paired ratio is consistent with unity within the 95%
  bootstrap CI under the primary estimator; jackknife shows no single-run dependence; the swapped
  map inverts and the scrambled map collapses the amplitude control as expected.
- **CONDITIONAL:** ratio consistent with unity but method- or run-dependent beyond the CI.
- **FAILED:** ratio significantly ≠ 1 under the confirmed map (then the kill-shot statement must be
  rewritten as a measured material difference, with its own CI).

## 5. Outputs
Script `mixedKillshotBootstrap.C` (this dir) → `papers/figures/mixed_killshot_bootstrap/`
(bootstrap distribution; jackknife plot; width comparison vs E; wrong-map/amplitude control;
method table printed + figure). Log: `killshot_bootstrap_result.txt` (this dir).

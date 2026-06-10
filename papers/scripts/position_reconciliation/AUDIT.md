# Position Reconciliation Audit (pre-analysis) — GATE 2

Date: 2026-06-09. Written BEFORE the reconciliation script, per gated-workflow protocol.
Goal: determine exactly what "1.5 mm" was, whether it is defensible, and what millimeter claim
Paper 3 is allowed. NOT to make the number look good.

## 1. The original "1.5 mm" — located and characterized

Source: `analyze/studies/showerLocalization.C` (commit lineage; quoted in
`papers/energy_position/README.md` as "1.52 mm x, 1.44 mm y @150 GeV").
Mechanism (read from code, lines 36–57):
- Estimator: `TLinearFitter` per coordinate, model `hyp3` = const + 3 linear terms in the corner
  light FRACTIONS f_i = A_i/ΣA, where A_i = lg_peak[i] + lg_peak[i+4] (corner = D+U ends summed;
  i = NW,NE,SE,SW). The 4th fraction is omitted (fixed by Σf = 1).
- Target: wire-chamber position relative to a HARDCODED module centre (6.6, 4.7) mm; window
  |x_w|,|y_w| < 6 mm; selection wc_ok + MCP window + ΣA > 500 mV.
- The number: `resx = sqrt(lx.GetChisquare()/n)` — i.e. the **RMS of unbinned event-level residuals
  of the linear fit, evaluated on the SAME events used to fit it (in-sample / training residual)**.
  It is NOT a Gaussian σ, NOT robust, NOT a binned-mean closure, NOT a standard error.
- So: event-level residual width — YES; train/test separated — NO (4 fitted parameters on ~10⁵
  events ⇒ overfitting bias is expected to be negligible, BUT this must be demonstrated, not assumed).
- ⚠ Schema rot found: the script reads branch `mcp_peak`, which DOES NOT EXIST in the current
  39-branch schema (now `mcp1_peak`). Against current files the MCP cut operates on an unset
  variable. The historical 1.5 mm was produced against an earlier reduction. The reconciliation
  script must use the current schema and re-derive everything.

## 2. The "3.6 mm" comparator — located and characterized

Source: `analyze/studies/wireChamberResolution.C` (data-driven, no second tracker; outputs in
`output/Summary/wire_chamber_resolution.{pdf,root}` when run).
Mechanism (read from code):
- Delay-line chamber; x = kWC_Scale·(t_R − t_L), kWC_Scale = 7/36 mm/ns; y likewise from (t_D − t_U).
- The 3.6 mm is **σ_x(UB) = kWC_Scale · σ(t_R + t_L) with PEAK-SAMPLE timing** — the spread of the
  END-TIME SUM. The script itself labels it an **UPPER BOUND**: the sum carries the common
  beam-arrival/trigger jitter t₀, "which CANCELS in the position difference t_R − t_L".
  It is therefore per-coordinate, peak-sample-timing, t₀-INFLATED — and **not the event-level noise
  of x_trk**. The t₀-cancelling cross-check (sumX−sumY)/√2 is flagged UNRELIABLE in-code
  (anti-correlation + secondary lobe). A CFD-50% sub-sample variant is also computed (expected smaller).
- Canonical DWCs are ~200 µm devices (Spanggaard); the binding limit here is the 1 GS/s digitizer
  (position quantization ≈ kWC_Scale × 1 ns ≈ 0.19 mm/sample), so the true difference-mode noise
  plausibly lies between ~0.3 and ~3.6 mm — it has never been measured directly in this setup.

## 3. What is mathematically allowed (stated before plotting)

- An event-level residual against a truth device with genuine event-level noise σ_WC CANNOT have
  width below σ_WC: σ²_resid = σ²_det + σ²_WC(diff-mode). If a measured residual of 1.5 mm coexists
  with a claimed 3.6 mm comparator, then EITHER the quantity is not an event-level residual width,
  OR the comparator number does not apply. **Reading the code resolves it: the quantity IS an
  event-level residual RMS, and the 3.6 mm does NOT apply (it is the t₀-inflated sum-side bound,
  not the difference-mode noise).** Corollary: the measured residual itself bounds the WC
  difference-mode noise: σ_WC(diff) ≤ σ_resid.
- A binned mean closure, centroid linearity, or fit-parameter uncertainty can legitimately be far
  smaller than any comparator resolution and must never be called event-level spatial resolution.
- Intrinsic unfolding σ_det = √(σ²_resid − σ²_WC) is permitted ONLY with a measured, applicable
  σ_WC(diff) and σ_resid > σ_WC; subtracting the inapplicable 3.6 mm from 1.5 mm gives an imaginary
  number — the script must print this check explicitly and refuse the unfolding.

## 4. Data and branches (current schema, verified 2026-06-09)

- WC truth: `x_trk`, `y_trk` (mm, WC frame), `wc_ok`, `wc_peak[4]` (plane amplitudes).
- Light division: `lg_peak[8]` (LG amplitudes; never clip). Corner sums A_i = lg_peak[i]+lg_peak[i+4].
- Quality: `mcp1_peak` (200–750 mV window — CORRECT current branch name), ΣA > 500 mV.
- Energy: per-file (`beam_energy`); primary DSB1 150 GeV (1.83 M entries), cross-checks 50 GeV
  (different beam tune) and 25 GeV (lowest light); fiducial variants: |x_w|,|y_w| < 6 (original
  window) vs r < 2.5 mm (timing fiducial).
- Module centre: derived per file from the LG-weighted beam centroid (RadView::beamCenter), NOT the
  hardcoded (6.6, 4.7).

## 5. Pre-registered outcomes

- **Outcome A — event-level position resolution established**: only if the held-out unbinned
  residual width exceeds a MEASURED, applicable comparator term and unfolding yields real σ_det.
  (Not expected: no measured σ_WC(diff) exists yet.)
- **Outcome B — comparator-limited event-level residual**: residual width consistent with the
  comparator contribution; language "comparator-limited transverse localization" only.
  Variant B′ (expected from §3): the residual **jointly bounds both terms**: σ_det ≤ σ_resid AND
  σ_WC(diff) ≤ σ_resid; language "the held-out event-level residual of X mm is an upper bound on
  both the capillary localization and the tracker term".
- **Outcome C — binned centroid closure only**: if the small number turns out to be a binned-mean
  artifact; language "sub-cell centroid response/closure" only. (Not expected — code shows unbinned.)
- **Outcome D — arithmetic/branch-definition error**: PARTIALLY REALIZED ALREADY: (i) the memory
  characterization "1.5 mm residual vs 3.6 mm event-level comparator = impossible" mischaracterized
  the 3.6 mm (sum-side UB, not difference-mode noise); (ii) the original script is schema-rotted
  (`mcp_peak`). Both go into the memory updates regardless of the rest of the result.

## Script + outputs
`positionReconcile.C` (this directory) → unbinned residuals (train AND held-out), binned closure,
comparator decomposition with the imaginary-unfolding check, fiducial + energy variants.
Figures → `papers/figures/position_reconciliation/`; log → this directory.

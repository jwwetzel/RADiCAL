# Paper 2 — Outline & report-derived inputs (energy + position + 4D)

Scope and intended figures: see `README.md`. This file collects the concrete analyses already in hand for
Paper 2, mined from the DSB1 detector report (`site/report/index.html`) and the docs drafts. **Every number
below is DSB1 single-build** unless noted — the per-build energy comparison (configResolutionFull.C) is a
separate config set; use these as the DSB1 anchor / position-infrastructure, not multi-build claims.

## §X Detector & selection inputs (the shared reference frame)
The position, energy, and 4D analyses all sit on one calibrated frame:
- **Module center = (6.58, 4.66) mm** in WC coordinates, square ~15.5×15.5 mm footprint, from the Σ_HG
  half-max edges of the 8 capillaries; cross-checked with Σ_LG to <0.05 mm; energy-independent to ±0.15 mm
  (25–150 GeV) → geometric, not a beam centroid; sets `kCalo_x0/y0 = (6.6, 4.7)`. *(moduleCenter.C; module_center.png)*
  — this is the **position origin** the localization needs.
- **MCP–RADiCAL alignment = 0.79 mm** (MCP (6.05, 4.08) vs RADiCAL (6.58, 4.66)); beam steered 2.48 mm
  off-center (centroid (8.93, 3.87)), both small vs the ~14–16 mm footprint → beam stays well inside fiducial.
  *(alignmentAnalysis.C; alignment.png)* — the common frame for the 4D demo.
- **Wire-chamber spatial resolution = σ_x ≈ 3.7, σ_y ≈ 3.6 mm** (peak-time; ~3.3 mm with sub-sample CFD-50%),
  measured single-tracker via the position-blind end-time-sum spread (σ_x = 7/36·σ(sum)); energy-independent;
  set by the slow ~1.3 ns/sample (1 GS/s) digitizer, not wire pitch. *(wireChamberResolution.C; wire_chamber_resolution-1.png)*
- **WC position formula:** x = (7/36 mm/ns)(t_R − t_L), y = (7/36)(t_D − t_U), CFD-50% delay-line peak times,
  `kWC_Scale = 0.194 mm/ns`.
- **Binning guardrail:** position quantizes in ~0.25 mm steps (one sample); hit/beam maps are fine at 1–2 mm,
  but **position-dependent analysis must use bins ≥ resolution (~3.5–4 mm)** — finer bins are smearing-correlated.

## §3 Energy at shower maximum
- **σ_E/E = 13.9 % @150 GeV** from the localized shower-max LG sum (~3 LYSO tiles ≈ 1 Molière radius; module is
  25 X₀ deep). The large constant term is **restricted transverse sampling / leakage, NOT longitudinal leakage.**
  *(layer5_energy.png)* — the DSB1 anchor for the README's "≈14 %, leakage-limited" claim.
- **Energy is read on the LG chain:** LG stays linear 25–150 GeV (hugs y∝E, ~55 mV/GeV); HG compresses above
  ~50 GeV by design (HG = timing). *(layer1_linearity_lg.png)* — the energy-observable justification.
- **Energy observable = SumLG Gaussian fit on timing-fiducial events**, per energy (6 runs). *({E}GeV_energy_distribution)*

## §4 Transverse position from the corner capillaries
- **Per-channel transverse maps:** mean pulse amplitude vs track impact for the 8 capillaries + 1×1 trigger,
  before/after the 1×1 trigger cut (ch15 > 60 mV, removes ~24 % no-hit pedestal), binned 0.25 mm x / 0.5 mm y.
  *(transverseMaps.C; transverse_maps_150GeV)* — the method/binning behind the localization.
- **Response is position-independent within the fiducial:** HG/LG amplitude ratio flat across the beam spot for
  every capillary (1 mm bins); no shadow/hot/dead region; SE-U/SW-U run ~10–15 % low (longer upstream WLS path).
  *(chargeProfiles.C; hg_lg_ratio_maps)*
- **⚠ Position-residual guardrail:** the README's "≈1.5 mm residual vs WC" is **below the WC's own ≈3.6 mm
  resolution** — present the 1.5 mm as **limited by the tracker comparator** (cite wireChamberResolution.C), not
  as a calibrated 1.5 mm capillary resolution, or a referee will object.

## §5 The E-type (full-length) energy capillary
- (No new TENERGY numbers in the DSB1 report; see README + `etypeEnergy.C` / `etypeChar.C`.)

## §6 Containment / event selection (shared with Paper 1)
- **4-population containment classification:** Pop A (good EM, passes 30 % cut), B (hadronic punch-through/edge,
  6.2 % @150), C (beam halo), D (off-axis); the SumPb < 0.30·SumLG cut separates A from B; containment falls
  toward r=3 mm (edge leakage), rises with energy; efficiency 95.3 % @150. *(pbglass_investigation; cross_energy_containment)*
  — the energy-selection cleanliness input + a useful Σ-scatter QC figure (contained-EM horizontal vs halo vertical bands).

## §6′ Simultaneous time, energy, position (4D, TENERGY)
- The **module center, MCP alignment (0.79 mm), and WC resolution (3.6 mm)** above are the shared-frame inputs
  `fourDdemo.C` needs to co-register the three observables. The 4D numbers (σ_t ≈57/best-bin ~35 ps, σ_E/E
  ≈11.9 %, σ_x ≈0.9 mm) are in the README.

## Scope flag
- Every report number here is **DSB1-only**. Keep σ_E/E = 13.9 % as the **DSB1 point**; the README's per-build
  13.5–15 % spread is the multi-build configResolutionFull.C result.

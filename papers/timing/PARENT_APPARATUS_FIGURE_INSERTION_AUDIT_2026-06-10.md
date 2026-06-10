# Parent-Apparatus Figure Insertion Audit — 2026-06-10

Insertion #1 of `PARENT_FIGURE_REUSE_REPORT_2026-06-10.md`: the three-panel apparatus/geometry
composite for Sec. 2.1. GEANT4 longitudinal-profile insertion (#2) remains DEFERRED to the
pre-submission tier, pending whether Sec. 7 needs more support after coauthor review.

## Sources

| Panel | Source file (`papers/timing/NIM_A_Figures/`) | Parent figure | Reuse mode |
|---|---|---|---|
| (a) module schematic | `lyso_w_cell_15.pdf` | Fig. 1 of Ref. [1] (NIM A 1068 (2024) 169737) | **as-is + minimal overlay labels** (beam arrow, UP/DW end tags); vector PDF embedded unmodified |
| (b) T-type capillary | `lyso_w_cell_16.pdf` | Fig. 6 of Ref. [1] | **as-is + minimal overlay labels** (quartz capillary, WLS filament, to-SiPM arrows); vector PDF embedded unmodified |
| (c) beam's-eye corner map | `lyso_plate_12.pdf` | Fig. 2 of Ref. [1] | **REDRAWN** (TikZ) — the engineering drawing is the dimensional source only; nothing is copied graphically |

License/credit: the parent paper is CC BY-NC-ND (Elsevier open access); the figures are
collaboration-own. Caption carries "adapted from Ref. [1]". Formal credit-line wording
("Adapted from [full citation], © 2024 The Authors, CC BY-NC-ND 4.0") is a submission-time
formality, tracked in the reuse report.

## Orientation and corner-map verification (the non-negotiable check)

What the authoritative record PINS (verified, not inferred):
1. **Corner names + channel pairing** — the author-confirmed 2023 logbook channel map
   (`dataset_2023.md`): kCap order NW-D/NE-D/SE-D/SW-D/NW-U/NE-U/SE-U/SW-U with the full
   [drs,grp,ch] HG and LG maps. The names NE/NW/SE/SW used in this figure are exactly the
   channel-map names used by every analysis script.
2. **MIXED material assignment** — NE+SW = DSB1, NW+SE = LuAG:Ce: data-inferred
   (`radcore/inferMixed.C`, per-corner brightness @150 GeV, recorded in `MIXED.json`) and
   independently CONFIRMED by the GATE-1 amplitude-independent pulse-shape discriminants
   (low-gain charge-to-peak and CFD rise time, validated on the pure builds). This is the
   post-fix map used by the same-shower control (Sec. 6) and by
   `makeMixedKillshotFigure.C` (DSB1 pair = NE−SW, LuAG pair = NW−SE).
3. **Hole geometry** — 14.0 mm square tile, 4 × ⌀1.3 mm corner through-holes at ±3.50 mm,
   1 × ⌀0.9 mm central monitoring-fiber hole, from the parent Fig. 2 engineering drawing;
   identical to the canonical module record (`detector_invariant.md`).

What the record does NOT pin (and how the figure handles it honestly):
4. **Absolute compass orientation in the beam's-eye view** (whether NE is upper-right or
   upper-left when looking downstream) is not recorded in the repo: the channel map is
   electronics-only, and the position estimator is a TRAINED linear fit
   (`positionReconcile.C`), so its coefficients absorb any orientation without recording it.
   THEREFORE panel (c) DECLARES its convention explicitly — "beam's-eye view, looking
   downstream; N up, E right (drawing convention)" — rather than asserting a physical
   orientation. The same-shower control of Sec. 6 depends ONLY on the diagonal pairing
   (opposite corners share a material), which IS pinned (items 1–2); no claim in the paper
   depends on the absolute orientation. Logbook confirmation of the physical orientation is
   folded into the existing coauthor question Q6 (corner-map logbook check).
5. **UP/DW end assignment in panels (a,b)**: the parent drawings carry one physical anchor —
   the T-type WLS filaments sit toward the LEFT end of Fig. 6, and shower maximum is in the
   upstream half of the 25 X₀ stack — so the composite marks beam entering from the left,
   left card = UP, right card = DW, consistent with that anchor. This too is stated as the
   drawing's convention; the timing observable (t_DW − t_UP)/2 is symmetric under relabeling
   up to a sign that is fixed by the measured negative depth drift (Sec. 7), not by this figure.

## Newly annotated vs copied

- Copied from Ref. [1] unmodified: the two vector schematics (a,b), including their original
  internal labels (W 2.5 mm / LYSO 1.5 mm / quartz capillary / monitoring fiber / SiPM /
  readout card; the WLS material-family label in (b)).
- Newly added by this paper: the beam arrow + UP/DW tags (a); the capillary/WLS-filament/
  to-SiPM callouts (b); the entire panel (c) drawing (corner labels, MIXED diagonal markers,
  convention note, dimension callouts from the parent's own engineering numbers).
- Panel (c) colors: neutral grayscale only (filled dark-gray = DSB1 corners, open = LuAG:Ce
  corners), so no color implies a performance or material claim.

## What this figure supports / must not support

SUPPORTS: the geometry definitions used by the analysis (module cross-section and length,
capillary positions, dual-end readout, the corner naming convention); the reader's ability to
follow the (DW−UP)/2 construction (Sec. 4) and the MIXED same-shower pairing (Sec. 6).

MUST NOT be used to support: the MIXED material assignment itself (the caption states the
assignment comes from the pulse-shape audit, NOT from this schematic); any radiation-hardness,
light-yield, or timing-performance claim; the absolute transverse orientation of the module in
the beam line; monitoring-fiber instrumentation (the central fiber was not in the 2023 readout
channel map — the label is geometric only).

## Build chain

`papers/scripts/apparatus_composite/apparatus_composite.tex` (standalone/TikZ; includes the two
parent PDFs + the redrawn panel (c)) → tectonic → `papers/timing/figs/radical_apparatus_composite.pdf`
(vector, used by the manuscript) + `.png` (raster preview copy). Inserted as the new Fig. 1
(figure*, top of Sec. 2.1), label `fig:apparatus`; subsequent figures renumber automatically.

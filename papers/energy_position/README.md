# Paper 2 — Shower-max energy & transverse position (and a 4D demo)

**Working title:** *Shower-maximum energy and position measurement with wavelength-shifting
capillaries, and a simultaneous time–energy–position demonstration, in a RADiCAL calorimeter
module.* (See the prose draft `docs/PAPER2_DRAFT.md`, sections 1–8.)

**Status:** skeleton. The localization analysis and final per-build figures still need to be
built/finalized; `radical_energy_position.tex` is a structural stub. This folder is the home for
the formal manuscript; working notes/draft prose stay in `docs/`.

## Thesis & scope
Energy and transverse position come from the **same shower-maximum measurement** — the four
corner-capillary amplitudes give a transverse centre-of-gravity (position, validated vs the wire
chamber; the module-ID concept), and their absolute sum gives energy. The E-type (full-length)
capillary is characterized in-beam, and the **4D demo** shows time + energy + position
simultaneously in one **TENERGY** module — the culminating demonstration of the published
§7 / Fig. 28 enhanced module. Covers the published roadmap **§7.1** (EM shower localization),
the energy / E-type studies, and **§7.3/7.4** (enhanced modular design).

Note (the physics link to Paper 1): timing and position are the **same observable** read two ways —
(DW−UP)/2 ∝ longitudinal shower depth, which is the timing *floor* and the position *signal*;
the depth-corrected estimator bridges them.

**Configs:** TENERGY (3×DSB1 T-type + 1 E-type) for the 4D demo and E-type; per-build energy for
the energy-resolution comparison.

## Intended figures (macros under `analyze/studies/`; current outputs in `site/assets/`)
| figure | macro | shows |
|---|---|---|
| `config_sigmaE_vs_E`, `config_sigmat_vs_E` | `configResolutionFull.C` | per-build σ_E/E and σ_t vs E |
| `shower_localization` | `showerLocalization.C` | capillary-light position vs wire-chamber truth |
| `transverse_maps_150GeV` | `transverseMaps.C` | per-capillary transverse response maps |
| `etype_energy_resolution`, `etype_vs_ttype` | `etypeEnergy.C`, `etypeChar.C` | E-type (energy) vs T-type response & linearity |
| `fourD_demo` | `fourDdemo.C` | one TENERGY module: time + energy + position together |
| `sigmat_vs_lightyield` (supporting) | `lightYieldTiming.C` | σ_t vs light yield (consistency view) |

Authoritative figure index (status LIVE-PAPER2 entries): memory `figure_catalog.md` → "Paper 2".

## Results so far (from the memory; to be finalized)
- **Energy at shower max:** σ_E/E ≈ **14 % @150 GeV**, essentially WLS-species-independent (13.5–15 %
  across DSB1/LuAG/MIXED); a localized, leakage-limited measurement consistent with the published
  52 %/√E stochastic term.
- **Transverse position:** residual ≈ **1.5 mm** (1.52 mm x, 1.44 mm y) @150 GeV vs the wire chamber,
  below the ~2.9 mm beam-spot spread.
- **E-type capillary:** ~4–5× dimmer at shower max (1.6 vs 6.7 mV/GeV) but linear in energy; excluded
  from the corner (DW−UP)/2 timing estimator.
- **4D demo (TENERGY, 150 GeV):** σ_t ≈ 57 ps all-fiducial (~35 ps best-bin), σ_E/E ≈ 11.9 %,
  σ_x ≈ 0.9 mm — from one event sample.

## Build
`tectonic radical_energy_position.tex` (first run fetches the elsarticle class). `\graphicspath{{figs/}}`;
copy the chosen figures into `figs/` before building.

## Pointers
- **Report-derived inputs already in hand:** `OUTLINE.md` — module center (6.58, 4.66 mm), MCP alignment
  (0.79 mm), wire-chamber resolution (~3.6 mm, the comparator that limits the ~1.5 mm position result),
  transverse-map method, the DSB1 energy backbone (σ_E/E = 13.9 % @150), and the 4-population containment.
- Prose draft: `docs/PAPER2_DRAFT.md`; outline: `docs/PAPER_OUTLINES.md`.
- Plan & physics: memory `radical_apparatus_conclusions.md` (★★ PAPER PLAN, §7 mapping).
- Apparatus / E-type / channel map: memory `experimental_setup.md`.

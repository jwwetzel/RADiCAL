> **Imported into the repository 2026-07-21 (workspace audit, Phase A4).** This file was
> maintained as analysis-assistant working memory during the 2026 publication campaign and is
> imported verbatim so the published workspace is self-contained: a reader needs no external
> notes. Line references cite repo files as of the import commit. Source: memory/detector_invariant.md.

# RADiCAL Detector — Tier 1 INVARIANT (the timeless concept)

> Tier 1 INVARIANT — reusable across all campaigns. Per-campaign instance: `dataset_<year>.md`. Methods: `methods_invariant.md`.

**What this is.** The conceptual specification of the RADiCAL ultracompact W/LYSO:Ce shashlik
EM-calorimeter module — the design that is true of *any* RADiCAL campaign, independent of which year,
which builds, or which numbers were measured. It encodes the INVARIANT once; everything that is a
per-campaign instance (channel-map indices, run numbers, SiPM bias voltage, clip values, measured
jitter/light yields, the data manifest) is **parameterized** in `dataset_<year>.md` — do not restate
those here.

**Parent / reference paper.** Perez-Lara, J. Wetzel, U. Akgun, et al. (RADiCAL Collaboration),
**arXiv:2401.01747 = NIM A 1068 (2024) 169737** — "Study of time resolution measurements and prospects
for energy resolution ...". The authoritative apparatus reference for the module concept.

---

## 1. The ultracompact W/LYSO:Ce shashlik design

An **ultracompact W/LYSO:Ce sampling EM calorimeter**: interleaved tungsten (W) absorber plates and
LYSO:Ce scintillator plates in a shashlik geometry, with thin reflective (Tyvek) inter-layer spacers to
avoid optical absorption at the LYSO/W surfaces. The module is **~25 radiation lengths (X₀) deep**
(≈0.9 nuclear interaction lengths) but only **~1 Molière radius wide** — a "finger-sized" volume.

**Rationale.** Contain an EM shower in a finger-sized volume and read its light *at shower maximum*. The
deliberately narrow (~1 R_M) transverse cross-section means the instrumented corner capillaries sample
the few LYSO:Ce tiles around shower max — a *localized shower-max* measurement, not full lateral
containment. This is what makes the same single physical structure yield three observables at once:
- **Timing** — a bright, localized, fast pulse near shower max;
- **Position** — dual-ended readout of the corner capillaries;
- **Energy** — longitudinal sampling of the shower.

LYSO:Ce is inorganic and radiation-hard; tungsten gives the short X₀ that makes the package compact.

---

## 2. The WLS capillary is the only variable (the standard-candle principle)

**INVARIANT WORDING (never violate): LYSO:Ce scintillator + W absorber are common to ALL builds. The
variable is the WLS capillary material** (DSB1 organic vs LuAG:Ce ceramic). Never say "crystal differs"
/ "DSB1 (LYSO) build" / "LYSO vs LuAG" — say "WLS differs", "DSB1 (organic WLS)", "DSB1 vs LuAG:Ce WLS".

This is the core experimental philosophy: a **one-variable-at-a-time / standard-candle** design. Because
the scintillator, absorber, geometry, readout, reference, and DAQ are all held fixed, changing only the
WLS-capillary fill gives an apples-to-apples comparison between candidate wavelength shifters
(e.g. a high-light-yield organic reference vs a radiation-hard ceramic candidate).

Four capillaries penetrate the full length of the module at its transverse corners. A separate central
monitoring fiber is distinct from the four instrumented WLS capillaries.

Filament-material roles (concept):
- **Organic WLS (e.g. DSB1):** high light yield, the reference / standard candle.
- **Ceramic WLS (e.g. LuAG:Ce):** radiation-hard candidate, lower light ("dim").

The timing filament is deliberately thin/small so that timing is set by **collected light**, not by
scintillator speed; the WLS material is varied for radiation hardness.

---

## 3. T-type vs E-type capillaries (the two capillary roles)

Two capillary geometries serve the two physics goals:

- **T-type ("timing"):** the WLS filament sits *only* in the shower-max region; the rest of the capillary
  core is clear quartz (waveguide). This gives a sharply localized, bright, fast pulse — optimized for
  **timing**.
- **E-type ("energy"):** the WLS runs the *full* capillary length, giving a uniform longitudinal
  response — optimized for **energy**. The E-type is several× dimmer at shower max than the T-type but is
  linear in energy. It is **excluded** from the corner timing estimator because including it degrades the
  four-corner combination.

(Build recipes — how many T-type vs E-type, and which corner is which material in a given campaign —
are a per-campaign instance; see `dataset_<year>.md`.)

---

## 4. Dual-gain readout: HG (timing) / LG (energy), with deliberate HG saturation

Each SiPM end is read out through **two parallel chains**:
- **High gain (HG) → timing.** Differentially amplified for a steep, fast rising edge. The timing chain
  runs at *very* high gain on purpose, so the rising edge is as steep as possible — **deliberately
  accepting full saturation / clipping of the pulse top** (digitizer-channel saturation, not SiPM-pixel
  saturation). Timing is taken from the edge, so clipping the top is an acceptable trade.
- **Low gain (LG) → energy.** Stays linear across the full energy range; carries the localized
  shower-max energy.

The design choice "saturate the HG on purpose for edge steepness, carry energy on the LG" is invariant.
A method consequence (parameterized per campaign): the clipped HG edge can be **recovered** by predicting
the true unclipped peak from the LG and timing on the steep edge below the clip.

HG (timing) pulses are negative-going; pedestal-subtracted, inverted amplitudes are returned positive.

---

## 5. Geometry: 4 corners × both ends → 8 SiPM timing channels

Four corner capillaries, each read out by a SiPM at **both** ends — its **upstream (UP)** and
**downstream (DW)** end — giving **8 SiPM timing channels for 4 capillaries** (4 corners × 2 ends).
The four corners are labeled by position (NW, NE, SE, SW). Each capillary thus produces a (DW, UP) pair.

This both-ends, four-corner layout is what enables the differential timing idea below and the
transverse center-of-gravity (position) and amplitude-sum (energy) reconstructions.

---

## 6. MCP-free differential timing (the concept)

The headline timing estimator is the **corner differential (DW−UP)/2** — mean downstream times minus
mean upstream times, divided by two (called *BestMinus* in the parent paper).

Why it works (invariant reasoning):
- **(DW−UP)/2 is a longitudinal-position estimator.** t_down − t_up is proportional to *where along the
  fiber* the light was produced, so the spread of (DW−UP)/2 is the **shower-depth fluctuation** — shower
  physics, near the ~1 X₀ depth-fluctuation limit, not an instrumental limit.
- **The difference cancels the references.** Taking down minus up at the same corner cancels the timing
  reference, the digitizer timebase error, and beam-arrival jitter — so the estimator is **MCP-free**: it
  does not depend on an external MCP/clock to reach its resolution.
- **The √2 subtraction penalty becomes a √2 averaging gain** when combined over the four corners
  (σ = σ_ch/√2).

The complementary observable, the **absolute shower time** (mean of the SiPM edge times referenced to an
external MCP), is the *other* timing output the same data supports; the MCP serves as an absolute
reference there but is **not required** for the (DW−UP)/2 depth resolution. Position and timing are the
same underlying measurement read two ways: (DW−UP)/2 is both the timing *floor* (depth fluctuation) and
the position *signal* (depth).

---

## 7. Program targets (goals, not any one module's measured values)

The RADiCAL program design targets — radiation-hard, ultracompact package:
- **~10 ps** EM-shower timing,
- **~1 mm** position,
- **~10%/√E** stochastic energy resolution.

These are program goals that motivate the design; measured per-campaign performance lives in
`dataset_<year>.md`.

---

## Cross-references (do not duplicate here)

- **Per-campaign instance** (channel-map indices, run numbers, SiPM bias voltage, the HG clip value,
  measured jitter / light yields / resolutions, the build recipes, the data/run manifest):
  `dataset_<year>.md` (2023 → `dataset_2023.md`).
- **Methods** (reduction pipeline, extraction, calibration, selection-cut machinery, analysis
  primitives): `methods_invariant.md`.
- **Per-paper thesis + figure index:** `paper_timing.md` (Paper 1) · `paper_energy_position.md` (Paper 2).
- **Global index:** `MEMORY.md`.

## Canonical module realization & WLS materials (design dimensions + material physics)
The realized RADiCAL module (the 2023 build; a future campaign notes any deviation):
- **Stack:** 29 × LYSO:Ce plates (1.5 mm) interleaved with 28 × tungsten plates (2.5 mm), Tyvek
  reflective sheet between layers. Envelope **14 × 14 × 135 mm³**, **~25 X₀** (~0.9 λ),
  **X₀ = 5.4 mm**, Molière radius **R_M = 13.7 mm** (~1 R_M wide → a *localized* shower-max
  measurement, not full lateral containment). [arXiv:2401.01747 §2; 135 mm length author-confirmed —
  the plate sum 28×2.5 + 29×1.5 = 113.5 mm is thinner than the envelope; the ~21.5 mm balance is
  Tyvek wrapping + optical interfaces + end/support structure. The CALOR-talk "~12 mm" R_M is rounded;
  use 13.7 mm.]
- **WLS capillaries:** quartz (fused-silica) waveguides, **183 mm long, 1150 µm OD / 950 µm ID**;
  inside each T-type capillary a WLS filament **900 µm × 15 mm** sits at the EM-shower-max region,
  the rest clear quartz. **DSB1** = organic WLS (absorption 425 nm, emission 495 nm, fluorescence
  decay τ = 3.5 ns); **LuAG:Ce** = radiation-hard ceramic WLS. **E-type** = WLS the full capillary
  length (uniform longitudinal response → energy). [arXiv:2401.01747 §2.]
- **SiPMs:** Hamamatsu **HDR2** (8 channels/module). **MCP-PMT reference:** Hamamatsu **R3809U-50**,
  intrinsic 10–20 ps. [arXiv:2401.01747 §4.]
- **Collaboration / parent paper:** Perez-Lara, J. Wetzel, U. Akgun, et al. (RADiCAL Collaboration),
  **NIM A 1068 (2024) 169737** = arXiv:2401.01747. **39 authors; 11 affiliations** — Adiyaman U.,
  Caltech, Coe College, Hofstra U., U. Iowa, Istanbul U., Istanbul U.-Cerrahpasa, Istanbul Technical U.,
  U. Notre Dame, U. Virginia, Yildiz Technical U. The collaboration sets author order; the per-paper
  author block lives in each `papers/<topic>/radical_*.tex`.

> VERIFY: the capillary physical/optical specs (183 mm, 1150/950 µm, 900 µm×15 mm filament, 425/495 nm,
> τ=3.5 ns, Tyvek) trace ONLY to the arXiv:2401.01747 quotation — no in-repo file to re-check against.
> (HDR2, R3809U-50, and the σ_E/E coefficients are independently corroborated in `docs/FOLLOWON_PLAN.md`.)

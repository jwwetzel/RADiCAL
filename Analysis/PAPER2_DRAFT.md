# Shower-maximum energy and position measurement with wavelength-shifting capillaries, and a simultaneous time–energy–position demonstration, in a RADiCAL calorimeter module

**RADiCAL Collaboration** — *[author list to be completed]*

*Draft v0.1 — prose first draft. Numbers from `configResolutionFull.C`,
`showerLocalization.C`, `etypeChar.C`, `fourDdemo.C`. Figures from
`Analysis/capillary_figs/`. Placeholders marked [...]. Convert to elsarticle (NIM A)
when settled. Companion to the timing paper (`PAPER1_DRAFT.md`) and to [1].*

---

## Abstract

The localized scintillation light produced at electromagnetic shower maximum, sampled
by wavelength-shifting (WLS) capillaries in an ultra-compact RADiCAL calorimeter module,
encodes both the shower energy and — through the relative light in the four corner
capillaries — its transverse position. Using electron beams of 25–150 GeV at the CERN
SPS H2 line, we report (i) the shower-max energy resolution, which is ~14 % and
essentially independent of the capillary scintillator (DSB1, LuAG:Ce, or mixtures); (ii)
a transverse position resolution of ≈1.5 mm obtained from the four-corner light, with no
assumption on the module orientation, validated against a wire chamber; and (iii) the
in-beam behaviour of a full-length ("E-type") energy capillary, which collects ~4–5×
less light at shower maximum than the shower-max-localized ("T-type") timing capillaries
yet responds linearly in energy — a uniform energy element. Finally, we demonstrate a
single module reading **time, energy and position simultaneously**: a module instrumented
with three T-type timing capillaries and one E-type energy capillary measures σ_t ≈ 35 ps
(energy-binned), σ_E/E ≈ 12 % and σ_x ≈ 1 mm in the same events — a proof of concept for
the four-dimensional RADiCAL module. These results realize the shower-localization and
energy/E-type directions of the RADiCAL R&D programme [1].

*Keywords:* electromagnetic calorimetry; shower-maximum sampling; position
reconstruction; 4D calorimetry; wavelength shifter; LuAG:Ce; SiPM; FCC.

---

## 1. Introduction

Future-collider electromagnetic calorimetry aims to measure, simultaneously and with
high granularity, the energy, time and position of each shower — "4D" calorimetry. The
RADiCAL concept samples the scintillation light at electromagnetic shower maximum, where
the energy deposition is greatest and most localized, and transports it to silicon
photomultipliers through wavelength-shifting (WLS) capillaries. A companion measurement
[1] established the timing capability (27 ps at 150 GeV), and a separate study
[companion timing paper] shows the timing is set by detected light yield, not by the
scintillator species.

This paper addresses the *energy* and *position* halves of the shower-max measurement,
which are physically inseparable: the same four corner capillaries that sample the local
energy also localize the shower transversely, because the relative light among them is a
center-of-gravity of the impact point. As [1] notes, the key value of the shower-max
energy measurement "is for the determination of EM shower position localization" and as
"a pattern-recognition constraint" to identify which module a shower struck. We quantify
both, characterize a dedicated full-length energy capillary, and close with a single
module measuring time, energy and position at once. Together these realize the
shower-localization and energy/E-type items of the RADiCAL programme [1].

---

## 2. Module and capillary types

The module is the ultra-compact 25 X₀ W/LYSO:Ce sampling calorimeter of [1] (29 LYSO:Ce +
28 W plates, 14 × 14 × 135 mm³, R_M = 13.7 mm), instrumented with four capillaries at the
transverse corners, each read out at both ends by SiPMs in dual gain.

Two capillary types are used, differing in where the WLS filament sits [1]:
- **T-type** — the WLS filament is positioned only at shower maximum; the rest of the
  capillary is clear quartz. The response is sharply localized at shower max — bright and
  fast, optimized for **timing**.
- **E-type** — the WLS runs the **full length** of the capillary, giving a uniform
  longitudinal response — optimized for **energy**.

The builds analyzed are DSB1, LuAG:Ce, and 2×DSB1+2×LuAG (all four T-type), and a
"TENERGY" build of three T-type DSB1 capillaries plus **one E-type energy capillary** at
the NW corner. The E-type identity of the NW capillary is confirmed by the run logbook
and directly by the data (it is ~4–5× dimmer than the three DSB1 corners). [Table 1:
builds, capillary types, energies.]

---

## 3. Energy at shower maximum

The energy deposited around shower maximum is estimated from the sum of the low-gain
(unsaturated) capillary signals. As detailed in [1], because the low-energy tail of the
measured-energy distribution would otherwise bias the timing, events are analyzed in bins
of measured energy. Figure 1 [`config_sigmaE_vs_E.png`] shows the shower-max energy
resolution σ_E/E as a function of beam energy for each build; it improves with energy to
**≈14 %** at 150 GeV and is **essentially independent of the capillary scintillator**
(13.5–15 % across DSB1, LuAG:Ce and the mixed build). For DSB1 the energy dependence is
consistent with the published parametrization (stochastic term ≈52 %/√E) [1].

We emphasize, following [1], that this is a **localized shower-max** energy resolution —
the capillaries sample only the ~3 LYSO:Ce tiles around shower maximum, in a module ~1
Molière radius wide — not a full-containment measurement; the latter requires a multi-
module array. Its value here is as the energy estimate underlying the energy-binned
timing and the position reconstruction of §4, and as the energy element of the 4D
demonstration of §6.

---

## 4. Transverse position from the corner capillaries

The transverse impact point of the shower follows from the **relative light in the four
corner capillaries**: a shower landing toward one corner deposits more light there, so the
light center-of-gravity tracks the position. Rather than assume the module-to-beamline
orientation, we determine the best linear estimator of the wire-chamber position from the
four corner light fractions and take its residual as the resolution (geometry-agnostic).

Figure 2 [`shower_localization.png`] shows the resulting capillary position estimator
versus the wire-chamber truth for 150 GeV electrons: a clean monotonic correlation in
both x and y, with the characteristic edge compression of a four-corner center-of-gravity.
The residual is **≈1.5 mm** (1.52 mm in x, 1.44 mm in y), well below the 2.9 mm spread of
the beam spot — i.e. the four-corner light carries real position information. As this
residual is the capillary estimator and the wire-chamber truth in quadrature, both are
≤1.5 mm here. This localized position, available per shower from the same signals that
measure the energy, also provides the module-identification / pattern-recognition handle
discussed in [1] for disentangling overlapping showers in a modular array. [Extension to
all energies and builds: §; a finer position truth would sharpen the capillary-only
resolution.]

---

## 5. The E-type (full-length) energy capillary

The TENERGY build carries one E-type capillary, allowing an in-beam comparison of the two
capillary types — the in-beam analog of the bench measurement of [1] (their Fig. 4).
Figure 3 [`etype_vs_ttype.png`] compares the E-type (NW) to the three T-type DSB1 caps.
Two findings stand out. First, the E-type collects **~4–5× less light at shower maximum**
(1.6 vs 6.7 mV/GeV) — expected, since its WLS is spread over the full length rather than
concentrated at shower max — yet its response is **linear in beam energy**, confirming it
as a valid, uniform energy element. Second, and contrary to the naive expectation that a
full-length cap must time poorly, its **single-channel timing is comparable** to the
T-type (~200 ps at 150 GeV): the prompt shower-max light still sets the leading edge. The
division of labour between the two capillary types is therefore one of light-collection
geometry and intended role — localized/bright for timing, uniform/linear for energy — not
of raw timing capability. (Consistently, the E-type corner is nonetheless excluded from
the (DW−UP)/2 timing estimator, because including it degrades the four-corner combination;
see the companion timing paper.)

Crucially, the E-type capillary delivers a **better per-capillary energy resolution** than
the T-type. Figure 4 [`etype_energy_resolution.png`] compares σ_E/E for the single E-type
cap, a single T-type cap, and the full eight-capillary low-gain sum: the E-type reaches
**~14% — better than a single T-type cap (~17–18%) and comparable to the entire 8-cap sum
(~14–16%)**. Despite being ~4–5× dimmer, the full-length E-type wins because the resolution
here is limited by *sampling* fluctuations, not photostatistics, and its uniform
longitudinal sampling captures the shower energy far more completely than a localized
shower-max sample. This directly validates the E-type as the energy element of the design.
A σ_E/E = a/√E ⊕ c fit gives a ≈ 50 %·√GeV (consistent with the published stochastic term
[1]) over a dominant constant term c ≈ 14 %; the near-flat energy dependence over 50–150 GeV
therefore reflects **constant-term (sampling / transverse-position) dominance** in this
localized shower-max, single-module measurement — not the absence of a stochastic term.
Exposing the steep a/√E behaviour requires full containment (a 3 × 3 array).

We note an important limitation: the 2023 run contains **no four-E-type configuration** (the
maximum is the single E-type cap of TENERGY), so a true energy-capillary-only measurement —
four full-length capillaries summed — is not yet available. The single-cap result above
already implies it would be excellent; a dedicated 4×E-type run is the clean next
measurement.

---

## 6. Simultaneous time, energy and position

The TENERGY module — three T-type timing capillaries plus one E-type energy capillary —
measures all three observables in the **same events**. Figure 4 [`fourD_demo.png`] shows,
for 150 GeV electrons: the (DW−UP)/2 timing of the three DSB1 capillaries (σ_t ≈ 57 ps
all-fiducial, ≈35 ps in the energy-binned best bin); the summed low-gain energy
(σ_E/E = 11.9 %); and the four-corner position residual against the wire chamber
(σ_x = 0.9 mm) — time, energy and position from one module, one event sample. This is a
direct proof of concept for the four-dimensional RADiCAL module envisaged in [1] (their
Fig. 28), in which dedicated T-type and E-type capillaries are deployed together for
simultaneous timing, energy and position measurement.

---

## 7. Discussion and outlook

The shower-max sampling delivers a coherent 4D measurement from a single ultra-compact
module: few-tens-of-ps timing, a ~14 % localized energy estimate, and ~1 mm position.
Each has a clear path to improvement. The energy resolution is limited by the localized
sampling and the single-module transverse leakage; a 3 × 3 array, by containing the full
shower, would reduce the stochastic and constant terms toward the calorimetric design
goals [1]. The position resolution is presently limited by the coarse wire-chamber truth
and the four-point sampling; a finer reference and additional capillaries would sharpen
it. The E-type capillary establishes the uniform energy element needed alongside the
T-type timing element. Combining all three in the enhanced modular design [1] is the
route to a full 4D RADiCAL calorimeter.

---

## 8. Conclusions

We have shown that the shower-maximum light sampled by WLS capillaries in a RADiCAL module
measures both energy and transverse position, and that a single module can measure time,
energy and position simultaneously. The shower-max energy resolution is ~14 % and
independent of the capillary scintillator; the four-corner light localizes the shower to
≈1.5 mm against a wire chamber; the full-length E-type capillary is a dim but linear,
uniform energy element with timing comparable to the T-type. In one TENERGY module we
measure σ_t ≈ 35 ps, σ_E/E ≈ 12 % and σ_x ≈ 1 mm in the same events — a proof of concept
for the four-dimensional RADiCAL calorimeter. These results realize the shower-
localization and energy/E-type directions of the RADiCAL programme [1] and motivate the
multi-module, mixed-capillary enhanced design.

---

## Acknowledgements
*[funding, CERN SPS/H2 operations, collaboration — to be completed]*

## References
1. RADiCAL Collaboration, *Study of time and energy resolution of an ultra-compact
   sampling calorimeter (RADiCAL) module at EM shower maximum over 25–150 GeV*,
   arXiv:2401.01747.
2. C. Hu et al., *LuAG ceramic scintillators for future HEP experiments*, Nucl. Instrum.
   Meth. A 954 (2020) 161723.
3. [companion RADiCAL timing/material paper — this work.]

# Comparison of DSB1 and LuAG:Ce wavelength-shifting capillaries for precision timing in an ultra-compact RADiCAL calorimeter module

**RADiCAL Collaboration** — *[author list to be completed]*

*Draft v0.1 — prose first draft. Numbers are harvested from the analysis pipeline
(`configResolutionFull.C`, `mixedHeadToHead.C`, `lightYieldTiming.C`,
`configSystematics.C`). Figures referenced from `Analysis/capillary_figs/`.
Placeholders marked [...]. Convert to elsarticle (NIM A) when content is settled.*

---

## Abstract

The RADiCAL programme develops ultra-compact, radiation-hard electromagnetic
calorimeter modules that measure both energy and time at electromagnetic shower
maximum. A previous test-beam study [1] established a timing resolution of 27 ps at
150 GeV using DSB1 organic wavelength-shifting (WLS) filaments in shower-max
("T-type") capillaries. Because the physics motivation for these modules is
operation in high-radiation environments, the wavelength shifter must ultimately be
radiation hard, and the ceramic scintillator LuAG:Ce is a leading candidate. Here we
compare DSB1 and LuAG:Ce T-type capillaries in the same module, beam and analysis at
the CERN SPS H2 line (electrons, 25–150 GeV). Using the MCP-free (DW−UP)/2 corner
estimator with a constant-fraction discriminator and run-folded out-of-sample
validation, we find that the per-capillary timing is governed by the **detected light
yield and is independent of the scintillator species**: a confound-free, in-event
comparison gives DSB1 ≈ LuAG:Ce, and the per-build differences are accounted for by
light yield. The uniformly bright DSB1 build
reaches 29.0 ± [..] ps at 150 GeV; builds containing the dimmer LuAG:Ce cluster about
8 ps higher, an offset fully accounted for by light yield. The total systematic
uncertainty is 1.1–1.7 ps. These results show that a radiation-hard wavelength shifter
delivers the same timing performance as the organic reference once it is made
sufficiently bright, validating the radiation-hard path for RADiCAL precision timing.

*Keywords:* electromagnetic calorimetry; fast timing; wavelength shifter; radiation
hardness; LuAG:Ce; SiPM; 4D calorimetry; FCC; CMS.

---

## 1. Introduction

Precision timing at the few-tens-of-picoseconds level is now a central requirement
for calorimetry at the energy frontier. At the High-Luminosity LHC and at proposed
future colliders (FCC-ee, FCC-hh), associating each electromagnetic shower with its
production time to ~30 ps suppresses pile-up and enables four-dimensional (4D) event
reconstruction. The RADiCAL programme addresses this with ultra-compact, radiation-
hard sampling modules in which scintillation light produced at electromagnetic
shower maximum is transported to silicon photomultipliers (SiPMs) through
wavelength-shifting (WLS) capillaries, providing a fast, spatially localized timing
signal [1].

A companion measurement [1] demonstrated a timing resolution of 27 ps at 150 GeV with
a RADiCAL module read out by four "T-type" capillaries — capillaries whose WLS
filament is positioned only in the shower-maximum region — filled with the organic
wavelength shifter DSB1. The motivation for the entire programme, however, is
operation under high irradiation, where organic shifters degrade. The radiation-hard
ceramic scintillator LuAG:Ce [2] is therefore a leading candidate to replace DSB1,
provided it can deliver comparable timing. Whether the *choice of scintillator*
affects the timing — or whether timing is instead set by the amount of detected light,
independent of the material — is the question this paper answers, realizing the
material-comparison item of the RADiCAL R&D programme outlined in [1].

We compare DSB1 and LuAG:Ce T-type capillaries directly, in the same module, the same
electron beam (25–150 GeV, CERN SPS H2), and the same analysis. Three builds are used:
a uniform DSB1 build, a uniform LuAG:Ce build, and a mixed build in which DSB1 and
LuAG:Ce capillaries are read out simultaneously in the same showers. The last permits
a confound-free, in-event comparison. We confirm that timing is set by light yield and
not by species, and we quantify the consequence for the choice of radiation-hard
shifter.

---

## 2. Experimental setup and datasets

**Module.** The RADiCAL module under test is an ultra-compact W/LYSO:Ce sampling
calorimeter, 25 radiation lengths (X₀) deep, described in detail in [1]: 29 LYSO:Ce
plates (1.5 mm) interleaved with 28 tungsten plates (2.5 mm), 14 × 14 × 135 mm³,
radiation length X₀ = 5.4 mm and Molière radius R_M = 13.7 mm. Four T-type capillaries
penetrate the full length of the module at its transverse corners; each contains a WLS
filament in the shower-maximum region and is read out by a SiPM at both its upstream
and downstream ends, giving eight timing channels. Each SiPM signal is amplified at
high gain (fast, for timing) and low gain (for the localized shower-max energy).

**Capillary builds.** Three configurations are compared, differing only in the WLS
material of the corner capillaries:
- **DSB1** — four DSB1 (organic) T-type capillaries;
- **LuAG** — four LuAG:Ce (ceramic) T-type capillaries;
- **MIXED** — two DSB1 and two LuAG:Ce T-type capillaries, read out in the same showers.

All three share an identical readout geometry, so a single analysis serves them all; a
config-agnostic reduction stores, per event, the (peak, constant-fraction time, charge)
of every digitizer channel together with the wire-chamber track and the timing
reference, with no channel map required at reduction time. [Table 1: run counts and
energies per build.]

**Reference and beam instrumentation.** A micro-channel-plate (MCP) photodetector
upstream provides a precision (~10–20 ps) timing reference; a delay-line wire chamber
measures the transverse impact point; a lead-glass array behind the module tags shower
leakage and hadronic contamination; scintillation counters trigger on single beam
particles. All channels are digitized by CAEN DT5742 (DRS4) samplers at 5 GS/s. The
H2 electron beam was taken in 25 GeV steps from 25 to 150 GeV, with at least 10⁶
triggers per step per build.

---

## 3. Timing reconstruction method

The time of each SiPM channel is extracted from the high-gain pulse and referenced to
the MCP. The module time is formed with the MCP-free **(DW−UP)/2 "corner" estimator**
("BestMinus" in [1]): the mean time of the four downstream channels minus the mean of
the four upstream channels, divided by two. Because each channel is referenced to the
same MCP, the difference cancels the reference time, the inter-group sampler jitter,
and the DRS4 cell-width error, exposing the intrinsic module resolution. As in [1], the
resolution is quoted in the best-contained measured-energy bin at each beam energy.

Two methodological choices differ from [1]. First, each pulse is time-stamped with a
**constant-fraction discriminator at 5 % of pulse height** rather than a fixed
threshold. The fixed-threshold analysis of [1] documented a ~15 % satellite population
shifted by ~0.2 ns, arising from the one-sample digitization quantization of the
leading-edge crossing; the constant-fraction crossing on the sharp low edge removes
this artifact. Second, the best-bin selection is validated **out-of-sample**: the bin
is chosen on one set of runs (5-fold, folded by run) and the resolution is measured on
the held-out runs, so the quoted value cannot be biased by selecting on a downward
fluctuation. Applied to the DSB1 build, this procedure reproduces the published
27.4 ps headline to within ~1.6 ps, validating the pipeline.

---

## 4. Results

### 4.1 Per-build timing resolution

Figure 1 [`config_sigmat_vs_E.png`] shows the out-of-sample best-bin timing resolution
as a function of beam energy for the three builds, with a/√E ⊕ b fits. The uniformly
bright **DSB1 build is lowest at every energy**, reaching **29.0 ps at 150 GeV**
(consistent with [1]); the **LuAG and MIXED builds cluster about 8 ps higher** at high
energy (37.4 and 38.2 ps at 150 GeV, respectively). The ordering is monotonic and
stable across the full energy range. [Table 2: σ_t(E) per build; fit coefficients.]

That the LuAG-containing builds are uniformly worse is, on its own, ambiguous — it
could reflect the crystal species or simply the lower light yield of the LuAG:Ce
capillaries. The next two subsections separate these.

### 4.2 In-event material comparison

The MIXED build reads DSB1 and LuAG:Ce capillaries in the **same showers**. Because all
timing channels share a common digitizer group, the pairwise difference of any two
capillary times cancels the event time, the reference, and the shared shower time,
leaving each capillary's intrinsic jitter — a confound-free comparison independent of
the reference and of any global offset. Figure 2 [`mixed_h2h.png`] shows the resulting
per-capillary intrinsic resolution for DSB1 and LuAG:Ce capillaries as a function of
energy. **Within uncertainties the two materials are identical** (~200 ps per
capillary, no systematic difference, equal at 150 GeV). With everything else held fixed
in the same showers, the scintillator species does not change the per-capillary timing.

### 4.3 The per-build differences track light yield

The in-event result of §4.2 is the clean, confound-free demonstration that species does
not matter; the remaining question is whether the *per-build* ordering of §4.1 is then
explained by light yield. Figure 3 [`sigmat_vs_lightyield.png`] plots the
reference-subtracted single-channel timing resolution of every capillary, at every
energy, against its light yield (the mean low-gain amplitude, an unsaturated relative
measure of detected light), with DSB1 and LuAG:Ce distinguished by colour. Both
materials occupy the **same σ_t–light-yield band**, falling together with light yield;
the dimmer LuAG:Ce capillaries lie at the low-light-yield end. We do not over-interpret
the small residual offset between the binned material means: this cross-build view is
**confounded**, because at equal light yield a (brighter) DSB1 capillary corresponds to
a lower beam energy than a LuAG:Ce one, and beam energy independently affects the
single-channel resolution. The clean statement is the in-event head-to-head of §4.2;
Fig. 3 is consistent with it. Together they imply that the per-build differences of §4.1
are an expression of light yield: the corner estimator is limited by its dimmest
capillary, so the uniformly bright DSB1 build wins, while a build containing dimmer
LuAG:Ce capillaries — of identical per-capillary quality — sits higher.

### 4.4 Systematic uncertainties

Systematic uncertainties are evaluated by varying the analysis cuts — fiducial radius,
shower-containment threshold and MCP amplitude window — and summing the resulting
shifts of the (stable, all-fiducial) resolution in quadrature. The total is **1.1–1.7
ps** across builds (DSB1 1.3, LuAG 1.1, MIXED 1.7 ps), comparable to [1] and dominated
by the fiducial-radius variation. The MCP-window term does not apply to the MCP-free
corner estimator. [Statistical uncertainties: §; the MIXED build has lower statistics
at the highest energies.]

---

## 5. Discussion

The picture is consistent across three independent views: the per-build ordering
(§4.1), the in-event head-to-head (§4.2), and the light-yield dependence (§4.3) all say
the same thing — **timing precision is governed by the detected light, not by the
scintillator species.** This follows from the underlying physics: the timing resolution
scales as the ratio of electronic noise to signal slew rate, σ_t ≈ noise/(dV/dt), and
what sharpens the leading edge relative to the noise floor is more detected
photoelectrons. A faster or slower intrinsic scintillator matters far less, at a thin
shower-max sampling element, than how much light reaches the SiPM.

The consequence for the programme is direct and favourable. LuAG:Ce is chosen for its
radiation hardness [2]; the concern is whether the ceramic shifter, which is dimmer in
the present capillaries, sacrifices timing. The data show it does not — *per capillary*,
LuAG:Ce times exactly as well as DSB1; the only penalty is its lower light yield, which
is an engineering parameter (filament geometry, coupling, photodetector efficiency),
not a fundamental limit of the material. A radiation-hard wavelength shifter made
sufficiently bright will therefore reach the same timing as the organic reference. This
removes a key risk from the radiation-hard module design.

[The deliberate high-gain saturation of the timing channels is timing-neutral at the
CFD-5 % operating point, since the constant-fraction crossing samples the intact low
edge below the clip; this is established in [Appendix / ref] and does not affect the
material comparison.]

---

## 6. Conclusions

We have compared DSB1 and LuAG:Ce wavelength-shifting capillaries for precision timing
in an ultra-compact RADiCAL calorimeter module, in the same beam and analysis as the
27 ps reference measurement [1]. By three independent methods — per-build resolution,
a confound-free in-event head-to-head, and the dependence on light yield — we find that
the per-capillary timing is **independent of the scintillator species and is set by the
detected light yield**. The uniformly bright DSB1 build reaches 29 ps at 150 GeV; builds
with the dimmer LuAG:Ce cluster ~8 ps higher, an offset fully explained by light yield,
with a 1.1–1.7 ps systematic uncertainty. Radiation-hard LuAG:Ce thus meets the timing
requirement once made bright enough, validating the radiation-hard path for RADiCAL
precision-timing modules. Increasing the LuAG:Ce light yield — through filament and
coupling optimization — is the clear next step toward a radiation-hard module matching
the organic-shifter timing.

---

## Acknowledgements
*[funding, CERN SPS/H2 operations, collaboration — to be completed]*

## References
1. RADiCAL Collaboration, *Study of time and energy resolution of an ultra-compact
   sampling calorimeter (RADiCAL) module at EM shower maximum over 25–150 GeV*,
   arXiv:2401.01747.
2. C. Hu et al., *LuAG ceramic scintillators for future HEP experiments*, Nucl. Instrum.
   Meth. A 954 (2020) 161723.
3. [CMS MTD TDR — timing-layer benchmark.]
4. [FCC-ee / FCC-hh CDR — physics motivation.]

# RADiCAL follow-on papers — outlines

Two papers from the May-2023 CERN-SPS data, both realizing the Section-7 roadmap of
**arXiv:2401.01747** (RADiCAL Collaboration; Perez-Lara, Wetzel, Akgun et al.).
Split by observable: **Paper 1 = timing**, **Paper 2 = energy + position (+ 4D)**.
Figure paths are in `Analysis/capillary_figs/` unless noted. Status tags:
✅ done · ◐ first-pass · ☐ TODO.

References for both: [1] arXiv:2401.01747 (self, the parent paper);
[2] Hu et al., *LuAG ceramic scintillators for future HEP experiments*, NIM A 954
(2020) 161723 (LuAG:Ce radiation hardness); [3] CMS MTD TDR (timing-layer benchmark).

---

# PAPER 1 — TIMING

**Title (working):** *Comparison of DSB1 and LuAG:Ce wavelength-shifting capillaries
for precision timing in an ultra-compact RADiCAL calorimeter module.*

**Thesis:** per capillary, DSB1 ≈ LuAG:Ce — **light yield, not crystal species, sets
the timing** — so a radiation-hard scintillator that is bright enough times just as
well. (Realizes ref [1] §7 item 2.)

**Abstract sketch:** Building on the 27 ps RADiCAL timing result [1], we compare two
WLS capillary fillings — DSB1 organic and the radiation-hard LuAG:Ce ceramic — in the
same module, beam, and analysis (CERN-SPS H2, electrons 25–150 GeV). Using the
MCP-free (DW−UP)/2 corner estimator [1] with a constant-fraction (5 %) discriminator,
and run-folded out-of-sample validation, we find the per-capillary timing is set by
the detected light yield and is **independent of the scintillator species**: an
in-event head-to-head gives DSB1 ≈ LuAG:Ce, and σ_t vs light yield places both
materials on a single curve. The uniformly-bright DSB1 build reaches 29 ps at 150 GeV;
builds containing the dimmer LuAG:Ce cluster ~8 ps higher, fully explained by light
yield. Radiation-hard LuAG:Ce therefore meets the timing goal once made bright enough.

**Sections & figures:**
1. **Introduction** — 4D calorimetry & HL-LHC/FCC timing; radiation hardness as the
   reason to test materials; the [1] result; this paper = its §7 item 2.
2. **Module, datasets & reduction** — the 25 X₀ W/LYSO module; 4 **T-type** capillaries
   read at both ends; the three builds **DSB1 / LuAG / MIXED** (all T-type); MCP, WC,
   DT5742; the config-agnostic reduction (`reduceRaw.C`). *Table:* run counts.
3. **Timing method** — CFD-5% (and that it removes the **0.2 ns / 15 % satellite** the
   fixed-threshold analysis of [1] documented); the (DW−UP)/2 "BestMinus" estimator
   [1]; energy-binned best bin; 5-fold run-folded OOS. *Optional fig:* CFD-fraction
   scan (`report/Summary/layer4_cfd_fraction_*.png`).
4. **Results**
   - σ_t(E) per build, OOS best-bin + a/√E⊕b fits — ✅ `config_sigmat_vs_E.png`
     (DSB1 29.0, LuAG 37.4, MIXED 38.2 ps @150; DSB1 reproduces the official 27.4
     within ~1.6 ps).
   - In-event, confound-free DSB1 vs LuAG per-capillary — ✅ `mixed_h2h.png`.
   - **σ_t vs light yield, DSB1 & LuAG one curve** (the thesis figure) — ◐
     `sigmat_vs_lightyield.png` (refine to reference-subtracted intrinsic σ_t).
   - Systematics — ✅ 1.1–1.7 ps (`configSystematics.C`).
5. **Discussion** — light yield ⇒ slew ⇒ timing; LuAG:Ce "bright enough"; the
   deliberate-saturation choice (CFD-5% timing-neutral, `desaturateCFD.C`, appendix).
6. **Conclusion** — species-independent; radiation-hard path validated.

**Remaining before submission:** ◐ intrinsic (reference-subtracted) σ_t-vs-LY via the
pairwise method; ☐ confirm LuAG config NW capillary type (T vs E); MIXED stats thin
@high-E (reduce more from Argon if needest). **Status: strong, ~80 % ready.**

---

# PAPER 2 — ENERGY + POSITION (+ the 4D combination)

**Title (working):** *Shower-maximum energy and position measurement with
wavelength-shifting capillaries, and a simultaneous time-energy-position
demonstration, in a RADiCAL module.*

**Thesis:** the shower-max capillary light yields **both** the local energy **and** the
transverse position; a module instrumented with timing (T-type) **and** energy (E-type)
capillaries measures **time, energy and position simultaneously** — a 4D
proof-of-concept. (Realizes ref [1] §7 item 1 + the energy/E-type direction.)

**Abstract sketch:** The localized energy deposited at EM shower maximum, sampled by
WLS capillaries in a RADiCAL module, determines both the shower energy and — through
the relative light in the four corner capillaries — its transverse position. We report
the shower-max energy resolution across capillary builds, a transverse position
resolution of ≈1.5 mm from the four-corner light center-of-gravity (validated against a
wire chamber), and the in-beam behaviour of a full-length (E-type) energy capillary.
A single module instrumented with three timing (T-type) and one energy (E-type)
capillary measures time (~35 ps), energy and position simultaneously — a proof of
concept for the enhanced 4D RADiCAL module.

**Sections & figures:**
1. **Introduction** — shower-max sampling → energy *and* position (quote [1]: the
   shower-max measurement's "key value is for … shower position localization"); the 4D
   goal; this paper = §7 item 1 + energy.
2. **Module & capillary types** — **T-type** (WLS at shower max, timing) vs **E-type**
   (full-length WLS, energy); the builds; **TENERGY = 3×T-type DSB1 + 1×E-type**.
   *Fig:* per-cap T/E identification (`configCapDiag.C` — the dim E-type corner). ◐
3. **Energy at shower max** — per-build σ_E/E(E) ✅ `config_sigmaE_vs_E.png`
   (build-independent ~13.5–15 %; extends the [1] DSB1-only result); the energy-binning.
4. **Position localization** — 4-corner light center-of-gravity → (x,y); resolution
   **≈1.5 mm** vs the wire chamber, below the 2.9 mm spot spread ✅
   `shower_localization.png`; the module-ID / pattern-recognition concept.
   *Optional:* transverse maps (`report/report_images/transverse_maps_150GeV-1.png`).
5. **The E-type energy capillary** — in-beam longitudinal uniformity & energy response
   of the full-length WLS cap (the in-beam analog of [1] Fig 4); TENERGY's NW cap. ☐
   (new analysis of the E-type signal.)
6. **Simultaneous time + energy + position (4D)** — the TENERGY module: 3 T-type DSB1
   give ~35 ps (≈ DSB1, `tenergyClean.C`) while the E-type cap measures energy and the
   four-corner light gives position — all in one event. Proof of concept for the
   enhanced module (ref [1] Fig 28). ◐
7. **Discussion / outlook** — toward the 3×3 array and full containment; full 4D.
8. **Conclusion.**

**Remaining before submission:** ☐ the E-type capillary in-beam characterization
(§5, the main new analysis); ◐ localization across energies & builds (currently DSB1
@150); ◐ position resolution with a better truth than the coarse WC; quantify the 4D
demonstration. **Status: substantial after adding localization; E-type analysis is the
main gap.**

---

## Shared status / next actions
- ✅ Corrections to depth (25 X₀) & estimator attribution (BestMinus = [1]) — live.
- ✅ Per-build σ_t(E), σ_E(E), systematics, light-yield, localization, T/E cap ID.
- ☐ Paper 2's E-type capillary in-beam analysis (longitudinal uniformity / energy).
- ☐ Confirm LuAG-NW capillary type (gates Paper 1's LuAG cleanliness).
- ◐ Refine the σ_t-vs-LY plot to intrinsic (reference-subtracted) σ_t.
- Draft order: **Paper 1 first** (nearly ready), then **Paper 2** once §5 (E-type) is done.

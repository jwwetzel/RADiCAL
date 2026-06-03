# RADiCAL capillary-materials comparison — CERN SPS, May 2023

Multi-configuration analysis of the four capillary builds taken in the May-2023
run: **DSB1**, **LuAG**, **2×DSB1+2×LuAG** ("MIXED"), **3×DSB1+1×Energy** ("TENERGY").
First-look, same-method numbers. (The validated, OOS-checked **27 ps DSB1 headline**
is the separate, polished analysis in the main report; these are exploratory
cross-config results.)

---

## TL;DR

1. **All four configs share the same 8+8 readout layout** — the labels describe the
   *crystal composition* of the capillaries, not the wiring. So one analysis geometry
   serves all of them.
2. A **config-agnostic reduction** (`reduceRaw.C`) turns the ~0.5 TB raw set into small
   per-config ntuples on Argon (≈12 MB/run, 160× smaller).
3. **In the same showers, DSB1 and LuAG capillaries time identically** (~200 ps per
   capillary, no systematic difference). **Crystal type is not the timing driver — light
   yield is.** This is the cleanest result and it confirms the detector's design premise.
4. The **high-gain / deliberate-saturation** strategy is roughly **timing-neutral**: CFD-5%
   lives on the sharp *low* edge, below the clip; the saturated peak is not used for
   timing. It neither breaks timing nor "mimics more light" (see §6).

---

## 1. Dataset & pipeline

Run logbook (`Analysis/hpc/logbook.csv`, 1278 rows). Electron, DONE, nominal-bias subsets:

| config (logbook label) | total runs | analyzed | energies |
|---|---|---|---|
| `DSB1` | 330 | 252 | 25–150 GeV |
| `LUAG` | 485 | 212 | 50–150 GeV |
| `2xDSB1, 2xLuAG` | 252 | 190 | 50–150 GeV |
| `3xDSB1, 1xEnergy` | 210 | 162 | 50–150 GeV |

Only DSB1 was run at 25 GeV. Reduction is config-agnostic (`reduceRaw.C`): for all 36 DRS
slots it stores `(peak, CFD-5% time, charge)` plus the config-invariant wire-chamber track
and both MCPs — no channel map needed at reduction. Argon recipe in
`Analysis/hpc/PRESCRIPTION.md`.

## 2. Channel layout (all configs identical)

`discoverReduced.C` on each config recovers the same structure:

- **8 HG timing capillaries** — DRS0 G0 ch0–6 + DRS0 G1 ch0
- **MCP1 / MCP2** — DRS0 G0 ch7 / DRS0 G1 ch7
- **8 LG energy capillaries** — DRS1 G0 ch0–7
- **Wire chambers** — DRS1 G1

In MIXED/TENERGY the DSB1 vs LuAG capillaries are told apart by amplitude (bright
~630 mV = DSB1, dim ~300 mV = LuAG).

## 3. Timing & energy resolution (same method, all-fiducial)

(DW−UP)/2 corner estimator, CFD-5%, r<3 mm fiducial, truncated-core σ — **all fiducial
events** (not the best-bin headline selection):

**σ_t [ps]:**

| E (GeV) | DSB1 | LuAG | MIXED | TENERGY |
|---|---|---|---|---|
| 25  | 70.2 | —    | —    | —    |
| 50  | 60.6 | 72.4 | 66.0 | 66.6 |
| 75  | 59.1 | 63.4 | 60.8 | 62.3 |
| 100 | 55.4 | 61.8 | 58.8 | 63.8 |
| 125 | 56.2 | 55.5 | 56.6 | 59.4 |
| 150 | 58.4 | 49.8 | 54.5 | 57.9 |

**σ_E/E ≈ 14–17 % for all configs, flat** — the short 4 X₀ prototype's leakage constant
term dominates regardless of crystal.

> **Selection-dependence (important).** All-fiducial favours LuAG at high E (49.8 < 58.4),
> but the **best-bin** selection (the one that sets the headline) flips it: DSB1 28 ps vs
> LuAG 36 ps at 150 GeV. Both best-bin numbers are tiny-efficiency (0.5–2.6 %) **in-sample
> overfits** — defensible per-config headlines require the run-folded **out-of-sample**
> treatment (TODO). The DSB1 best-bin (28 ps) reproduces the official 27.4 ps headline,
> validating the reduced pipeline. Figure: `compareConfigsPlot.C`.

## 4. The clean crystal test — in-event DSB1 vs LuAG (MIXED)

The MIXED module reads **both** crystals in the **same showers**. Timing slots 0–6 share
DRS0 Group 0, so the pairwise difference of any two CFD-5% times cancels event time,
group/run timebase offset, MCP, and (same-layer) shower-time — leaving exactly
√(s_i²+s_j²). Solving the pairwise system per layer gives each capillary's intrinsic σ_t,
**offset/MCP/proxy-free**. (`mixedHeadToHead.C`)

| E (GeV) | DSB1 caps σ_t | LuAG caps σ_t |
|---|---|---|
| 50  | 240 | 219 |
| 75  | 210 | 214 |
| 100 | 159 | 182 |
| 125 | 190 | 187 |
| 150 | 200 | 200 |

**Per capillary, DSB1 ≈ LuAG — no systematic difference.** This is the decisive,
confound-free comparison: with everything else held identical, the crystal choice does not
change the per-capillary timing.

## 5. Design concept — light yield drives timing (confirmed)

The timing capillary sits at shower-max and is thin and small, so its timing is set by the
**light it collects**, not by the scintillator's intrinsic speed — exactly what the in-event
test (§4) shows. The reason to test different materials is **radiation tolerance**: prove a
**bright** scintillator times well, then prove a **radiation-hard** one can be made bright
*enough*. The data supports the strategy — LuAG (more radiation-hard, dimmer) times
comparably per capillary, i.e. "bright enough works."

## 6. On the high-gain / deliberate-saturation strategy

The build deliberately uses very high HG gain so the rising edge is steep, accepting full
saturation, on the premise that *only the rising edge matters* and a "blown-up" edge gives
precise timing regardless of light yield. The data lets us check this:

**(a) Where on the edge is the timing best?** CFD-fraction scan on the saturated DSB1
channels (`cfdFractionDSB1.C`), (DW−UP)/2 σ_t [ps]:

| E | 5% | 10% | 20% | 30% | 50% |
|---|---|---|---|---|---|
| 50  | 60.7 | 62.0 | 62.1 | 64.2 | 76.3 |
| 100 | 55.5 | 56.3 | 57.5 | 58.6 | 61.4 |
| 150 | 58.4 | 61.6 | 63.9 | 64.5 | 65.4 |

Timing is **best at 5%** (the sharp low edge, well below the clip) and **degrades** as the
fraction climbs toward the saturated peak. So the *low* sharp edge is the useful part; the
clipped peak region is *not* usable for timing.

**(b) Does saturation hurt?** The in-event test (§4) is the clean answer: the **saturated**
DSB1 caps and the **unsaturated** LuAG caps time the **same** per capillary. So saturation
is roughly **neutral** for CFD-5% timing — it neither breaks it (DSB1 still reaches the
27 ps headline) nor helps it. *(An earlier slew-vs-amplitude plot suggested saturation
"hurt"; that plot was confounded by DSB1's clipping and by differing LG light yields, and is
superseded by the in-event result.)*

**(c) Why "high gain mimics more light" doesn't hold for timing.** Timing precision is
σ_t ≈ (noise)/(slew rate). A *linear* amplifier multiplies the slew **and the noise** by the
same gain, so σ_t = (G·noise)/(G·slew) is **unchanged** — you cannot improve timing by
turning up the gain. What genuinely sharpens the edge *relative to the noise floor* is **more
photoelectrons** (light yield + photostatistics), which gain cannot manufacture. Gain *is*
useful up to the point of lifting the signal clear of the **digitiser's fixed noise/
quantisation floor**; beyond that it only saturates.

**Net:** the deliberate-saturation choice is ~timing-neutral. It didn't cost the headline
(CFD-5% sits on the intact low edge), but it doesn't buy timing either — and it has real
costs: at high energy a clipped edge can't get steeper in its usable range (so you forfeit
the timing benefit of extra light there — consistent with DSB1's σ_t stalling above
~100 GeV while unsaturated LuAG keeps improving), and the HG carries no usable energy
information (energy must come from the LG chain).

## 6b. Apples-to-apples OOS headlines & the HG-vs-LG and de-saturation tests

**OOS best-bin (DW−UP)/2 σ_t @150 GeV, identical method across all four** (run-folded,
`configBestBinHGLG.C` for reduced configs, `dsb1OOSandCorr.C` for DSB1):

| config | σ_t (OOS best-bin) |
|---|---|
| **DSB1** | **29.4 ps** (official Gaussian-core pipeline: 27.4) |
| TENERGY | 35.7 |
| LuAG | 36.6 |
| MIXED | 36.7 |

DSB1 (uniformly brightest capillaries) wins by ~7 ps; configs containing dim LuAG caps
cluster at ~36 — the corner estimator is limited by its dimmest capillaries. OOS ≈ in-sample
throughout (defensible, not overfit).

**HG vs LG chain** (your proposed control), OOS best-bin σ_t @150 GeV: **HG ~36 ps vs LG
~64–88 ps.** The clean, unsaturated LG "perfect pulses" time ~2× *worse* than the saturated
HG — because LG is **slow** (low slew). Slew (a fast edge), not pulse cleanliness, sets
timing: the fast HG edge wins decisively even when its peak is clipped.

**LG→HG de-saturation idea** (reconstruct the true, unsaturated HG peak from the linear LG to
restore the CFD reference) — **now tested directly on the raw DSB1 waveforms**
(`desaturateCFD.C`). Per channel we fit the dual-gain slope `HG_true = b·LG` on unsaturated
events (b≈3–6), then on the saturated 150 GeV set we recompute the CFD crossing at *frac × the
reconstructed true peak* and compare (DW−UP)/2 σ_t (best-bin) head-to-head against standard
*frac × clipped peak*, in identical events:

| CFD point | standard (frac of clip) | de-saturated (frac of true) |
|---|---|---|
| **5 %**  | 32.4 ps | 32.6 ps |
| **20 %** | 35.3 ps | **31.1 ps** |

**At CFD-5% de-saturation buys nothing** (32.4→32.6): the 5 % threshold already sits on the sharp
*low* edge, far below the ~820 mV clip, so there is essentially no clipped-peak time-walk to
remove — reconstructing the true peak only adds its own LG-reconstruction jitter. **At CFD-20% it
works exactly as predicted** (35.3→31.1, a 4.2 ps best-bin gain): 20 %-of-clip lands high on the
edge near the saturation knee (a near-fixed voltage → real time-walk), and referencing 20 % of the
*true* peak removes that walk, pulling 20 % timing nearly back to 5 % performance. **Conclusion:**
the idea is sound and demonstrably effective *where saturation actually bites* (high CFD fraction),
but the detector's **CFD-5% operating point already sidesteps the clip by timing on the intact low
edge**, so de-saturation is moot for the headline — best achievable stays ~32 ps, set by light/slew,
not by the clipped peak. (The earlier correlation-only estimate of "marginal benefit" is superseded
by this direct measurement.)

## 7. Caveats & next steps

- **First-look only.** OOS-validated, best-bin per-config headlines + a/√E ⊕ b fits are TODO.
- MIXED's borderline-amplitude slot (SE-D) flips DSB1/LuAG classification across energies —
  mild noise on the §4 grouping; the DSB1≈LuAG conclusion is robust to it.
- TENERGY's "1×Energy" capillary is treated as a timing slot (all 8 HG slots are live);
  worth confirming its intended role.

## 8. Tools (all reproducible, committed)

`reduceRaw.C` · `discoverReduced.C` · `configResolution.C` · `configBestBin.C` ·
`configResolutionDSB1.C` · `configBestBinDSB1.C` · `mixedHeadToHead.C` · `cfdFractionDSB1.C` ·
`compareConfigsPlot.C` · `slewTest.C` (exploratory, confounded — kept for provenance) ·
`desaturateCFD.C` (LG-referenced de-saturation test on raw waveforms, §6b).
HPC: `Analysis/hpc/{PRESCRIPTION.md, reduceRaw via sge_reduce.sh, submit_reduce.sh, merge_reduced.sh}`.

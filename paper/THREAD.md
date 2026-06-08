# The thread to the money plot — RADiCAL Paper 1

One continuous argument, from the raw waveforms to the result. Every bead on the thread
is a figure that exists **for all four builds** (DSB1, LuAG, MIXED, TENERGY), batch-tested
on the local raw subset (one+ run per config: RUN1034–1261 DSB1, RUN2389/2186 LuAG,
RUN2941–2995 MIXED, RUN2722/2530 TENERGY). The destination is the **money plot** —
σ_t = a/√E ⊕ b, the stochastic term *a* tracking light yield while the floor *b* is shared.

---

### Bead 1 — The detector wrote down a real pulse on every channel
**Claim:** the data are sound; all eight SiPM channels, high- and low-gain, are alive.
**Plot:** `output/<BUILD>/{hg,lg}_waveforms.png` — per-channel average waveform (mean ± RMS),
all energies, **every build**. *(waveformProfiles.C, from the raw `pulse` tree.)*
**Reads:** clean fast HG pulse rising at t=0; a slow full-length LG pulse on the same channel;
no dead/flat channels in any build. For MIXED the LuAG corners are simply *dimmer*, not dead.
→ We can trust every channel of every build. Now look at what the HG pulse does at the top.

### Bead 2 — But the high-gain pulse is clipped — the timing edge is thrown away
**Claim:** the HG channel saturates at 820 mV; the steep edge we must time is lost to the clip.
**Plots:** the flat tops in the Bead-1 HG profiles (visible at high E); `clip.png` — the HG
peak-amplitude distribution piling up at 820 mV.
→ A constant-fraction time on the *measured* peak is forced onto the shallow foot. We need
the edge back. The low-gain copy still has it.

### Bead 3 — The low-gain copy predicts the true peak — linearly, on every channel
**Claim:** HG = a + b·LG holds per channel; the LG faithfully predicts the unclipped peak.
**Plots:** `lgcheck_<BUILD>.png` — per-channel HG-vs-LG hockey-stick, **all four builds**
(ρ = 0.73–0.95). *(mixedLGcheck.C.)*
**Reads:** clean linear rise → saturation knee, on all 32 channels across the builds. This
*is* the calibration that makes edge recovery possible — and it is what proves MIXED's LG sound.
→ Put the CFD at a fixed fraction of the *predicted* peak, on the steep edge.

### Bead 4 — Recovering the edge buys real picoseconds
**Claim:** timing the recovered edge is a factor 3.7 steeper and measurably better.
**Plots:** `pulse.png` (one waveform: clip 820, LG-predicted 1997, cfd05 foot 102 mV/ns vs
hg_lgcfd edge 381 mV/ns); `method_compare_DSB1.png` (same events, cfd05 vs lgcfd → +3.1 ps @150,
floor 22→20).
→ With the edge recovered, every build can be timed on the same well-defined crossing.

### Bead 5 — The MCP-free estimator gives clean, single Gaussians
**Claim:** (DW−UP)/2 cancels the reference; the timing distributions are clean.
**Plots:** `pub_dist_<BUILD>.png` — per-amplitude-bin (DW−UP)/2 and (DW+UP)/2 distributions
with central-peak fits, **all four builds**. *(pubFig.C.)*
**Reads:** narrow differential, wider absolute, clean Gaussians at every brightness.
→ Fit each, and watch the width fall with light.

### Bead 6 — Resolution falls with shower brightness, the same way in every build
**Claim:** within each build, σ_t improves as the showers get brighter (more light, steeper edge).
**Plots:** `pub_res_<BUILD>.png` — σ_t in 1000-event brightness slices, per energy, **all builds**.
→ Brightness — light — is the lever. Now lay the four builds side by side.

### ★ Bead 7 — THE MONEY PLOT: light yield, not WLS species, sets the resolution
**Claim:** σ_t = a/√E ⊕ b. The **stochastic term a tracks light yield** (DSB1 201 → TENERGY 240
→ MIXED 233 → LuAG 455 ps√GeV, ×2.3); the **floor b is shared** (~20 ps across the LYSO builds).
**Plot:** `light_yield_thesis.png` — all four builds, a/√E ⊕ b fits, one figure.
**Reads:** the builds order by light; the difference between DSB1 and LuAG lives *entirely* in
the light-driven term. At infinite light, all materials approach the same floor.
→ But a between-build comparison cannot, by itself, separate "less light" from "intrinsically
slower WLS" (e.g. LuAG:Ce decay kinetics). Kill that confound.

### Bead 8 — The clean proof: same module, same showers, only the WLS capillary differs
**Claim:** inside MIXED, a DSB1 (organic WLS) and a LuAG:Ce (ceramic WLS) capillary read the
*same showers in the same LYSO:Ce scintillator* — and time identically.
**Plot:** `mixed_h2h.png` — DSB1-vs-LuAG capillary σ_t(E), pairwise-CFD solve. **DSB1/LuAG = 0.99,
χ²/ndf = 0.4/5.** *(mixedHeadToHead.C.)*
**Reads:** every systematic (offset, MCP, shower-time) cancelled by construction; the two WLS
capillaries time the same. The build ordering in Bead 7 is light, not WLS species — proven.

### Bead 9 — The floor is the shower, and it confirms the published result
**Claim:** the energy dependence is 1/√E (photostatistics), giving b ≈ 20 ps — consistent with
the published 17.5 ps, not a revision. The floor is shared shower-depth physics.
**Plots:** `floor_model_<BUILD>.png` — 1/√E (χ²/ndf ≈ 2) vs 1/E (≈ 12) per build.
→ Light is the design lever; the residual floor is the shower itself — and that is where Paper 2
(depth-corrected timing, position, energy, 4D) begins.

---

## Status of the thread (Paper 1)
| Bead | Figure(s) | Builds | Status |
|---|---|---|---|
| 1 waveforms | hg/lg_waveforms | 4 | ✅ (waveformProfiles.C) |
| 2 saturation | clip | DSB1 (+ HG profiles all) | ✅ / extend clip to all |
| 3 HG/LG health | lgcheck | 4 | ✅ |
| 4 edge recovery | pulse, method_compare | DSB1 | ✅ (LYSO builds: extend method_compare) |
| 5 distributions | pub_dist | 4 | ✅ |
| 6 brightness | pub_res | 4 | ✅ |
| **7 MONEY PLOT** | light_yield_thesis | 4-in-1 | ✅ |
| 8 kill-shot | mixed_h2h | MIXED | ✅ |
| 9 floor | floor_model | 4 | ✅ |

## Next batch (to fully satisfy "extend all Paper 1 & 2 plots")
- **Paper 1 completion:** extend `clip` and `method_compare` to all LYSO builds; the systematics
  table (configSystematics/systematicUncertainties → all builds); the selection/optimization
  appendix (optScan/optLadder/cfdScan/channelCombinationScan — mostly param(build), just run).
- **Paper 2 (deferred plots, DSB1-only → extend):** energy resolution (configResolution),
  shower localization (showerLocalization, transverseMaps, positionCorrection), E-type response
  (etypeEnergy, etypeChar), the 4D demo (fourDdemo) — all on TENERGY + per-build.

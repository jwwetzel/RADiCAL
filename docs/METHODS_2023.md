> **Imported into the repository 2026-07-21 (workspace audit, Phase A4).** This file was
> maintained as analysis-assistant working memory during the 2026 publication campaign and is
> imported verbatim so the published workspace is self-contained: a reader needs no external
> notes. Line references cite repo files as of the import commit. Source: memory/methods_invariant.md.

# RADiCAL Methods — Tier 1 INVARIANT (the analysis playbook)

**Tier 1 INVARIANT — the analysis playbook, reusable across all campaigns. Numbers live in
`dataset_<year>.md`; paper-specific framing in `paper_*.md`.**

This file encodes HOW we analyze ANY RADiCAL campaign — the estimators, selection, fits, floor
physics, DRS4 quirks, and code abstractions that do not change when the beam/year/build changes.
The one physical invariant: **LYSO:Ce scintillator + W absorber are common to ALL builds; the only
variable is the WLS capillary (DSB1 organic vs LuAG:Ce ceramic). Never say "crystal differs" — say
"WLS differs".** When a concrete number is unavoidable below it is marked "(e.g. 2023: …)" — the live
tables live in `dataset_<year>.md`; per-paper framing in `paper_*.md`.

*Companion memories: `dataset_<year>.md` (the campaign's numbers/manifest/results), `paper_*.md`
(paper framing + figure index), `detector_invariant.md` (the detector concept), `MEMORY.md` (index). Hardware concepts that a
method depends on (channel groups, clip, dual-gain) are summarized here only enough to justify the
method; the full spec is the dataset file.*

---

## 1. The headline estimator: (DW−UP)/2 "BestMinus" — and WHY it cancels everything

The flagship timing observable is **(DW−UP)/2** = (mean of the downstream-end times − mean of the
upstream-end times) ÷ 2, over the valid corner capillaries. Called **"BestMinus"** in the parent
paper (its bins 6–8 of (DW−UP)/2). It is an **MCP-free differential**:

- **What cancels.** Both ends of a capillary carry the SAME per-event additive offsets — the shared
  reference time, the beam-arrival-in-window jitter, the DRS4 group timebase phase, and the SiPM/
  scintillation common-mode leading-edge shape. Subtracting DW−UP removes all of them by construction.
- **Why ÷2 and not √2 penalty.** Averaging the (correlated-noise-cancelled) ends turns the naive √2
  subtraction penalty into a √2 **averaging gain**: σ[(DW−UP)/2] = σ_channel/√2.
- **Reference referencing is per-group on purpose.** Each end is stored already referenced to the MCP
  copy *in its own DRS4 group* (see §6). So the same-group differential cancels the inter-group domino
  jitter EXACTLY — this is a deliberate channel-map design choice, not luck. (e.g. 2023: corners NW/NE/SE
  live wholly in DRS0 group 0 with MCP1; SW is split, SW-U in group 1 with MCP2; dataset_2023.md.)
- **What (DW−UP)/2 physically IS.** t_down − t_up ∝ the longitudinal depth where the light was produced.
  So (DW−UP)/2 is simultaneously the **timing floor** (its spread = shower-depth fluctuation, §4) and the
  **position/depth SIGNAL** for the energy+position paper — one observable read two ways. The depth-
  corrected estimator bridges them.

**Companion absolute-time observable.** The mean of the (MCP-referenced) leading-edge times,
`T_abs = mean_c(hg_led[c])`, is the absolute shower time (SiPM − MCP). It does NOT cancel the
common-mode edge shape, so it floors higher than (DW−UP)/2. Use `hg_led` (fixed-voltage edge, §3),
never cfd05-of-clipped, for anything absolute. (e.g. 2023: ~50 ps absolute vs ~25–28 ps depth.)

---

## 2. Selection: brightest-K, fiducial, containment — DETERMINISTIC cuts only

**Brightest-K slice (the "best achievable" reading).** Rank events by `sum_lg` (the LG light proxy)
and keep the brightest K (default **K=1000**; `rad::timingBrightestK(v,E,src,K)`). This gives identical
statistical tightness at every energy (same N) so energies can be compared cleanly, and quotes the
best-measured/best-contained σ_t. **A deterministic brightest-fraction (or quantile) cut does NOT
overfit** — it sits on the smooth σ(brightness) curve (proven by OOS run-folding, §3).

**Lock the selection — do not drift it silently.** THE headline = brightest-1000 + the best timing
source (§3) + (DW−UP)/2. Swapping to "brightest-2%" or a different K silently changes the number and
re-introduces non-monotonicity; James was burned by this repeatedly. If you must vary K, do it as a
labeled systematic (§5), not as a silent re-definition.

**Energy-bin best-bin (the published-style reading).** Slice `sum_lg` into bins and take the lowest-σ
bin with N≥500 (`rad::timingBestBin`). **Bins MUST be equal-POPULATION (quantile), never equal-WIDTH.**
σ_t falls monotonically with brightness, so the winner is the brightest bin; a fixed-width μ+2σ bin
sits in the Gaussian tail that EMPTIES as energy rises (the brightest bin can drop to a handful of
events → rejected → that energy capped at a dimmer bin → a false "higher E looks worse"). Equal-
population bins fix it and restore monotonicity. (This argmin-over-fixed-bins is also the ONE selection
that genuinely overfits — see §3; deterministic cuts avoid it.)

**Fiducial + centroid.** Center fiducial cuts on the per-run **signal-weighted beam centroid** (ΣHG of
the corner caps; ΣLG agrees as a cross-check), NOT the nominal geometric center. Timing fiducial radius
is OOS-optimized per energy (e.g. 2023: ~2.5 mm ≤100 GeV, 3.0 mm at 125/150; plateau at r≈3.0). Energy
uses a tighter radius (e.g. 2023: 2.0 mm).

**Containment.** Reject low-amplitude spikes from under-measured high-E tails: require ≥2 of 4 ends each
side AND a per-energy light floor (e.g. sum_lg > 0.40·median_E). A Pb-glass leakage veto rejects events
with sum_pb > 0.30·sum_lg (when sum_lg is above a threshold)). [STUDY-ERA NOTE, 2026-07: the Pb-glass containment veto is applied by the energy/uniformity STUDIES only — the production timing chain does not apply it; see the which-chain-applies-what summary in lib/physics/SelectionCuts.h.]

---

## 3. The timing SOURCE and robust width — and the HG/LG edge-recovery method

**The source enum (`RadView`: kCFD03, kCFD05, kCFD10, kCFD20, kCFD30, kCFD50, kLED, kLGCFD).** Every
event stores multiple per-end timing crossings; `RadView::timeOf(i,src)` / `hasSrc(src)` / `srcName(src)`
expose them, and the timing estimators take a `src` argument. Always run the per-build RANKING (one tree
pass over all available sources, fit each, rank) and pick the most ROBUST source per build — do NOT
hard-code one source for all builds.

- **cfd05** (CFD at 5% of the measured peak) = the published method. Fine on unclipped pulses; FRAGILE
  on small/slow/low-light pulses where 5% sits in noise on a shallow edge (crossing can jitter by ns →
  garbage (DW−UP)/2). So dim builds prefer **led** (fixed-voltage edge) or **lgcfd**.
- **led** (leading edge): time the pedestal-subtracted rising edge as it first crosses a FIXED ABSOLUTE
  threshold (e.g. 2023: 20 mV ≈ 4× HG noise), linear-interpolated between bracketing samples. Fixed-
  voltage, NOT a fraction of peak → immune to the clip that corrupts cfd05 (5% of a clipped 820 mV ≈
  41 mV of a bad peak). Best for dim/low-light builds and for ANY absolute-time work.
- **lgcfd = HG/LG edge recovery (paper name: srCFD; the key NEW method for high-light, heavily-clipped builds).** The HG
  timing chain hard-clips, so we are forced to time at the foot; lgcfd predicts the UNCLIPPED true peak
  from the linear LG and times on the steep edge below the clip:
  - Pre-pass fits `HG_peak = a + b·LG_peak` per end on UNCLIPPED events (the linear part of the HG-vs-LG
    "hockey stick"); per event `HG_true = a + b·LG_peak`; threshold = `lgcfd_frac·HG_true` (e.g. 2023:
    frac=0.15, the steep-but-safe sweet spot — higher fracs re-hit the clip on the brightest events and
    blow up). Crossing interpolated as a fraction of the measured peak, stored MCP-referenced as
    `hg_lgcfd[8]`.
  - **CALIBRATION IS CRITICAL: fit a,b on UNCLIPPED low-E data.** A per-file fit COLLAPSES at the top
    energy (all clipped → no LG lever-arm → slope→0). Fit per-channel from clean low-E points, take the
    MEDIAN of SANE slopes (auto-reject the over-clipped high-E AND any corrupted low-E energies), with a
    pulse-shape cut (drop charge/peak outliers = direct-hit/Cherenkov spikes), 3σ-trimmed. Write a per-
    build `<build>.hglg` sidecar (calibHGLG.C); `BuildConfig` reads it (`has_lgcal`); the Reducer applies
    it at EVERY energy, else falls back to a per-file fit WITH A WARNING.
  - The HG−LG transfer slope b is **energy-independent (drift ≤2%)** ⇒ ONE calibration per channel is
    correct. (A mild ~5% fit-range nonlinearity exists but is the same at every energy.) For energy work
    needing max rigor, fit a nonlinear HG(LG); not needed for timing.
  - **lgcfd is BUILD-DEPENDENT.** It wins where the pulse is high-light and heavily clipped; for dim,
    never-clipped, slower builds the foot (led/cfd03/05) wins and lgcfd is mid-pack. It buys a few ps at
    the headline + a cleaner floor; it does NOT improve the stochastic (1/√E) term. lgcfd's secondary
    benefit: it removes the ~0.2 ns fixed-threshold satellite population (a checkable claim — demonstrate
    in-data, do not assert).

**Robust Gaussian-core width `tebSigma` (never returns garbage).** To extract σ_t from a (DW−UP)/2
sample:
1. Compute a ROBUST core width: iterative 2.5σ truncated RMS, truncation-debiased (RMS_trunc = 0.9546·σ
   for a Gaussian, so divide by 0.9546).
2. Histogram (e.g. 120 bins over μ±4·RMS after a 5σ reject) and fit a 2σ-core Gaussian (`FitGaussCore`).
3. Use the Gaussian-core fit ONLY if it lands within 0.5×..2× of the robust core; otherwise fall back to
   the robust width. This keeps clean builds at their Gaussian value while preventing pathological low-
   light bins from returning 0.1 ps / 18000 ps garbage. Fit only the central peak (exclude the satellite
   shoulder) on the distribution figures.

---

## 4. FLOOR PHYSICS — a measured floor CONFIRMS the published value, never "revises" it

The (DW−UP)/2 σ_t falls with energy then PLATEAUS at high E. The plateau is **physics, not a failure to
extract the signal.** Decomposition (general form; plug the campaign's numbers from dataset_<year>.md):

| term | scaling | size (e.g. 2023) |
|---|---|---|
| slew = noise/(dV/dt) | improves with light (∝1/light) | ~2 ps for (DW−UP)/2 (already tiny → more light can't move the headline) |
| DRS4 cell-width | constant | <1 ps (negligible) |
| CFD/threshold optimization | one-time | ~2–3 ps gain only |
| **(DW−UP)/2 measured floor** | does NOT scale with light | the plateau itself |

- **The floor is shower-depth fluctuation.** (DW−UP)/2 is a longitudinal-depth estimator; its irreducible
  spread is the shower-to-shower fluctuation in EM shower-max depth, which is shower-development PHYSICS
  and does NOT scale with light yield. That resolves the paradox "more light, no improvement."
- **Self-consistency check: σ_z = v_fiber·σ_t.** With v_fiber ≈ 200 mm/ns, convert the measured σ_t to a
  depth spread and compare to **~1 X₀** of shower-max wander (X₀ ≈ 5.4 mm). They match (e.g. 2023: 25–28 ps
  → σ_z ≈ 4–5.6 mm ≈ 1 X₀). ⇒ the floor IS the ~1 X₀ shower-depth fluctuation — near the physics limit for
  an MCP-free corner-timing estimator, and material-INDEPENDENT by geometry (a prediction, not a coincidence).

**Photostat 1/√E beats slew 1/light — and the published floor STANDS.** Fit the brightest-K σ_t(E) with
both forms and compare χ²:
- **Standard photostat** σ_t = a/√E ⊕ b fits the energy dependence MUCH better than the amplitude/slew
  1/light form (e.g. 2023: χ²/ndf 2.1 vs 12.1; all builds agree photostat wins). The 1/light form OVER-
  extrapolates (it was the AMPLITUDE parametrization). ⇒ **adopt 1/√E ⊕ b.**
- A measured b consistent with the published floor **CONFIRMS** it; do NOT headline "we revised 17.5→25"
  or fight the published number. (e.g. 2023: photostat b ≈ 18–20 ps across LYSO builds ⇒ the published
  ~17.5 ps stands; the parent paper's brightest amplitude bins already reach 24–27 ps, so a measured ~25 ps
  is a confirmation/refinement, not a new result.)
- **Floor honesty rules.** (1) The floor b is an EXTRAPOLATION beyond the top beam energy — headline the
  MEASURED σ_t at the top energy, quote b with proper ±. (2) b is PER-TECHNIQUE: different sources (cfd05
  vs lgcfd) have genuinely DISTINCT floors; NEVER force a shared floor across techniques. (3) Fit b FREELY
  with per-point errors σ/√(2N), then PDG-scale by √(χ²/ndf) (point-to-point systematics, not stats,
  dominate the floor error); a↔b are strongly anti-correlated. (4) State plainly where a dim build's
  brightest top-energy point never actually reaches the extrapolated floor in data.

---

## 5. OOS run-folding / overfit avoidance + systematics

**OOS validation = MATCHED-FRACTION fold-by-run.** Split runs into folds; compare σ_t at a MATCHED
brightest-FRACTION across folds (NOT a matched K — a fixed K across half-size folds is a 2× looser
fraction and falsely looks unstable; this exact mistake was caught and corrected). At matched fraction,
fold-A ≈ fold-B ≈ full ⇒ a deterministic brightest-fraction/quantile cut is OOS-stable across the whole
selection continuum. The ONLY thing that genuinely overfits is the equal-WIDTH best-bin ARGMIN (picking
the luckiest of N fixed bins): in-sample optimistic, OOS worse — which is exactly why a published best-bin
headline IS already the honest OOS value and STANDS. Deterministic cuts avoid this; the argmin does not.

**Systematics budget — IN-MEMORY MIRROR of the locked pipeline.** Build a faithful mirror of the locked
`rad::timingBrightestK` that exposes the cut knobs the public API hides (fiducial r, MCP window, HG
threshold) and VALIDATE mirror==locked (Δ=0.00 ps every build) before trusting it. Then sweep the
selection: brightest-N (e.g. K=500/2000), fiducial r, MCP window, HG threshold, fit range. Report the
total selection systematic on σ_t at the top energy, and the dominant single variation (e.g. 2023: K=2000
adds +2 ps — it admits dimmer showers = the same light-resolution dependence, not an artifact). **EXCLUDE
the source choice (cfd05 vs lgcfd) from the systematic — that is the THESIS (a method gain), not a
selection systematic**, and it degenerates the low-light floor fit. The floor must SURVIVE the systematic
sweep for the "confirms published" claim to hold.

**Optimization landscapes (honesty appendix).** Show the σ_t surface vs each knob (brightest-N knee,
fiducial plateau, HG-threshold flat region) and show the adopted recovered-edge source lies BELOW the
entire clipped-peak CFD family for every build. Per-build best-bin landscapes diagnose the equal-width
bin bias (the brightest bin emptying at high E).

---

## 6. DRS4 / DT5742 quirks — and which cancel in same-group differentials

The DAQ is two CAEN DT5742 switched-capacitor digitizers (DRS4 chip), each = 2 inner groups × 9 channels ×
1024 samples. Quirks and their treatment:

- **Free-running domino + asynchronous trigger → the STOP CELL.** The trigger is asynchronous to the domino
  wave, so the stop cell rotates uniformly over [0,1023] (RMS ≈ 1024/√12). Recover it as the single zero-
  width step in the nominal time axis; store `stopcell[4]` per group. The time axis already rotates from
  the stop cell, so integer stop-cell correction barely helps the absolute inter-chip jitter (the residual
  is sub-cell trigger phase, removable only with a recorded common clock).
- **Same-group differentials cancel the domino/cell-width phase.** (DW−UP)/2 within one group, and
  MCP1−MCP2 (which share a group's stop cell with their channels) both cancel the cell-width timebase error.
  A cross-group stamp (t_HG − t_MCP with crossings far apart) does NOT — correct it with a stop-cell-indexed
  StopCellCorrection. Per-group MCP referencing (§1) is what makes the headline differential clean.
- **The "~100 ps" inter-group number is NOT the MCP.** σ(MCP1−MCP2)/√2 is FLAT with energy and amplitude-
  INDEPENDENT ⇒ a pure clock/digitization effect, not slew. The single MCP pulse is passively SPLIT into
  both groups, so the MCP's own intrinsic jitter (~10 ps) is COMMON-MODE and cancels; σ(MCP1−MCP2) measures
  the RELATIVE jitter of the two independent, free-running, not-phase-locked DRS4 domino oscillators (e.g.
  2023: ~71–74 ps inter-group). The published MCP intrinsic (~10–20 ps) is consistent — it is simply never
  directly measurable here (common-mode) and is sub-dominant when the MCP is used as a reference. This
  inter-group floor NEVER touches the headline because the corners share one group/domino/MCP copy.
- **Cell-width calibration is negligible for the timing board.** Data is written on a nominal uniform axis;
  measured cell-width RMS is float32-precision-limited (<1 ps) on the fast board → no per-cell width cal
  needed for timing. (A slow auxiliary board, if present, shows larger cell-width but is not used for timing.)
- **The HG clip is a DT5742 channel saturation, not SiPM-pixel saturation** (e.g. 2023: ~820 mV, lowered
  from nominal by a DC offset to capture the negative afterpulse). The legacy fixed saturation flag set for
  a different campaign may never fire — verify the clip level per campaign from the data, not the flag.

**Dual sampling rates (campaign-dependent — confirm per year).** The two boards may run at different rates
(e.g. 2023: DRS0 timing board @ 5 GS/s = 0.2 ns/cell carries all fine-timing channels incl. MCP copies,
Pb-glass, trigger scints; DRS1 @ 1 GS/s = 1 ns/cell carries LG energy + wire chambers). Confirm the rate of
each board against the campaign's logbook channel map before quoting per-cell time.

---

## 7. Absolute-time / TR0 method lessons (hard-won — don't repeat the bugs)

- **The absolute shower time is ALREADY in the reduced data.** The Reducer stores every HG time already
  MCP-referenced (`hg_led[i] = ledTime − ref`, ref = the same-group MCP copy). So `T_abs = mean_c(hg_led[c])`
  over the timing channels — there is NO second MCP subtraction. (DW−UP)/2 cancels the MCP (both ends carry
  the same −ref) → depth; mean(hg_led) keeps the −ref → absolute shower time.
- **Bugs that produced fake "ns-scale" absolute-time floors (ALL were analysis errors, not physics):**
  (a) DOUBLE-subtracting the MCP ((hg−mcp) − mcp_time); (b) a loose MCP peak cut letting junk events in
  (use the proper MCP window, e.g. 2023: min 200, max 750 mV); (c) a bad multi-run wildcard chain mixing
  files (verify per-run means agree before chaining). With clean leading-edge timing on RAW waveforms a
  single run already gives a sensible ~tens-of-ps absolute time; the reduced-data mean gives the headline.
- **Use the MCP as the absolute reference, TR0 only for inter-group sync.** The MCP is the better absolute
  reference; the discriminated trigger scintillator (TR0, split into ch8 of both groups) times worse but
  its `(t_TR0a − t_TR0b)` syncs the two free-running groups (~100 ps) as a cross-check. cfd05-of-the-clipped-
  peak injects large pulse-shape jitter into absolute SiPM times (it cancels in DW−UP so the headline never
  noticed) — fatal to absolute use; use `hg_led` instead.

---

## 8. Code abstractions — config-driven, year-parameterized, build-parameterized

The pipeline is **multi-dataset**: the campaign year and the build are parameters, not hard-coded paths.

- **`RAD_YEAR` + `DataPaths.h`.** `radYear()` reads the `RAD_YEAR` env (default = current campaign); the
  year flows as a parameter through `radReduced(build,E[,year])`, `radConfig(build[,year])` (build JSON),
  `radHglg(build[,year])` (the lgcfd sidecar), `radRaw(...)`. Switch campaigns with `export RAD_YEAR=...`
  — no recompile. Data layout: `$RAD_DATA/data/<year>/{raw,reduced/<BUILD>/<E>GeV.root,configs,metadata}`.
- **`FigPaths.h` year-namespacing.** Figures are campaign-namespaced so a new dataset never overwrites
  another's: `radFigP("figures/<sub>/x.png")` rewrites to `figures/<year>/<sub>/x.png` (idempotent, mkdir's
  the dir); `radFig(name, build, sub="narrative")` builds `figures/<year>/<sub>/<name>[_<build>].png`.
- **`BuildConfig` (config-driven reducer).** Per-build JSON (+ optional `.hglg` sidecar) carries the channel
  fill, clip level, lgcfd_frac, calibration, selection geometry. The Reducer is config-AGNOSTIC: it stores
  per-event observables for all DRS slots with no channel map needed at reduction; the map/material is
  applied at analysis time. `has_lgcal` flags a present lgcfd calibration.
- **`RadView`** = the per-event read interface over the reduced `rad` tree: `timeOf(i,src)` / `hasSrc(src)`
  / `srcName(src)` (the source enum kCFD03..kLED, kLGCFD), `hg_peak`, `is_timing`, LG-weighted `beamCenter`,
  fiducial helpers. Build-agnostic — the same analysis code runs on any build/year.
- **`RadTiming`** = the encapsulated headline pipeline: `tebSigma` (robust width §3), `timingBestBin(v,E,src)`
  (quantile-binned best-bin), `timingBrightestK(v,E,src,K)` (brightest-K slice). Both take a source arg.
- **Macros** take `(build[, year])` and call `radConfig(build)` etc. — e.g.
  `root -l -b -q 'analyze/timingHeadline.C+("DSB1")'`. `timingHeadline.C` runs the best-bin for EVERY
  available source in one tree pass, fits σ=a/√E⊕b, and RANKS them (pick the best/most-robust source per
  build). One-off studies live under `analyze/studies/`.
- **Run it:** `source setup.sh` (sets the lib include paths) then the macro. Heavy re-reduction runs on the
  HPC cluster from `reduce/hpc/` (compile → calibHGLG per build → submit_reduce → merge → rsync home); the
  cluster needs a `git pull` to pick up layout changes.

---

## 9. House plot style + ROOT/Cling pitfalls

- **House style:** `lib/viz/RADiCALStyle.h` (the RADiCAL ROOT style) + `lib/viz/PlotUtils.h` (`FitGaussCore`
  and shared plotting helpers). Use them for every figure so the look is consistent; write figures through
  the `FigPaths.h` helpers so they land year-namespaced.
- **Distribution figures:** fit the CENTRAL peak only (exclude the satellite shoulder), share x-range and use
  absolute shared-y across panels so per-energy/per-build distributions are visually comparable; draw the
  best-measured per-energy points as open circles when overlaying the floor fit.
- **ROOT/Cling pitfalls:** headers are bare-`#include`d (no library link) — `source setup.sh` must set the
  lib subdir include paths first or compilation fails. Compile macros with `+` (`macro.C+`). The no-crossing
  sentinel is a large negative value (kNoTime = −1e6) — always guard `if (t > −1e5)` before using a stored
  time, or sentinels poison means. `tebSigma` already guards against degenerate fits (§3); don't replace it
  with a bare Gaussian fit on low-light bins. A truncated/Zombie input file in a manifest silently corrupts a
  chained read — verify per-run means before chaining multiple runs/energies.

---

## 10. The invariant wording rule (repeat, because it is easy to get wrong)

LYSO:Ce scintillator and the W absorber are COMMON to all builds. The build label describes ONLY the WLS
capillary fill. The variable across builds is the **WLS capillary material** (DSB1 organic vs LuAG:Ce
ceramic) and the resulting **light level**. The timing ordering across builds is driven by LIGHT YIELD, not
WLS species (the in-module mixed-build comparison — same scintillator, same showers, same MCP, two
capillary materials side by side — is the confound-free proof). Never write "crystal differs", "DSB1 (LYSO)
build", or "LYSO vs LuAG"; write "WLS differs", "DSB1 (organic WLS)", "DSB1 vs LuAG:Ce WLS".

---

## 11. Reduced ntuple schema & pulse extraction (the invariant pipeline)
**Reduced `rad` tree.** The config-agnostic `reduce/reduceRaw.C` stores all 36 DRS slots; the named-format
`processRun.C` stores per-end arrays. Branches: `run, event, beam_energy, wc_ok, x_trk, y_trk, wc_peak,
sum_lg, sum_pb, pb_peak, mcp1_peak/time, mcp2_peak/time, tr0a/b_peak/time, stopcell[4], hg_peak[8],
hg_cfd03/05/10/20/30/50[8], hg_led/tot/charge/ped_rms/saturated/spike, hg_lgcfd[8], lg_peak[8], lg_charge`;
the config-agnostic form adds `s_peak[36], s_cfd05[36], s_charge[36]`. The `hg_*` crossing times are stored
**already MCP-referenced**. Output: ZSTD-5 at `data/<year>/reduced/<BUILD>/<E>GeV.root`.
**Pulse extraction** (`lib/waveform/WaveformUtils.h`): HG is negative-going; **pedestal = mean of samples
3–52**; peak = pedestal − min over samples 3–1023; **charge integrated over `[imin−15, imin+200]`** (~−3 to
+40 ns), positive contributions only. `ExtractPulseMulti`: LED/absolute threshold **20 mV**, minPeak 5 mV,
sat threshold = `hg_sat_mV` (per-build config); returns CFD crossings at 3/5/10/20/30/50% + LED + TOT. Spike
flag fires if any pedestal sample deviates **>5× ped-RMS** from the pedestal mean. MCP extracted at 20% CFD /
30 mV min; wire chambers at 50% CFD. Stop cell recovered by `drs4::FindStopCell` (the single zero-width step
in the nominal axis) → `stopcell[4]` = [D0G0, D0G1, D1G0, D1G1]; the cross-group `t_HG−t_MCP` (~130 cells
apart) is correctable via `StopCellCorrection`, while same-group differentials cancel it for free.

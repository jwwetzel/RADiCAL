# Statistical conventions — the timing analysis (one page, written once, reused everywhere)

Added 2026-07-21 (code audit, quick-win #6). Every statistical convention the published numbers
rely on, with the code location that implements it. If a macro appears to deviate, the macro's
in-file KERNEL-CLONE note states its deliberate delta; otherwise treat the deviation as a bug.

## 1. Width conventions — there are exactly TWO

- **Gaussian-core width (the convention of EVERY published resolution).** `rad::tebSigma`
  (`lib/physics/RadTiming.h`): iterative 2.5σ robust window → 2σ Gaussian-core fit inside it →
  tail guard (falls back to the truncation-debiased robust RMS only for pathological samples).
  All σ_t(E), all Table 1/2 numbers, the 25.7 headline and the ≈50 ps companion use this.
- **Tail-sensitive width (comparisons of estimator TAILS on identical events).** Fixed-window
  truncated RMS in a COMMON window frozen per sample (`cw`/`coreWindow` clones of the tebSigma
  window). Used ONLY in the method-gain/satellite studies (paper Sec. 5.3 core 1.3 ps vs
  tail-sensitive 4.0 ps distinction, App. B). The two conventions answer different questions and
  are never mixed in one number.

## 2. Statistical errors on σ points

σ/√(2N) with N = events in the selection (K=1000 → √2000) — `timingFitSummary.C`,
`methodGainPostfix.C`. This is the Gaussian σ-of-σ approximation; with N fixed by construction
it is uniform across builds/energies by design.

## 3. Fit-parameter errors

`σ_t(E) = a/√E ⊕ b` fits inflate parameter errors by the standard PDG scale factor
√(χ²/ndf) when χ²/ndf > 1 (`timingFitSummary.C`, `systematicsPostfix.C` floor block,
`methodGainPostfix.C` — the same convention at every fit; cited as [pdg] in the manuscript).

## 4. Selection-systematics total = RMS of signed shifts

Table 2's "total syst." is the RMS over the 8 single-knob signed shifts on identical events
(`systematicsPostfix.C`): ±1.0/±1.1/±0.9/±1.9 ps (DSB1/TENERGY/MIXED/LUAG). **Alternative
convention, for the referee who asks:** treating the 8 shifts as independent and summing in
quadrature gives ±2.8/±3.3/±2.4/±5.5 ps (arithmetic on the committed Table 2 shifts). The
paper's conclusions survive under either convention: the build ordering, the >2× stochastic-term
ratio (which comes from the fits, not this budget), and the DSB1-vs-LuAG separation
(25.7 vs 44.4 ps) are unaffected; the RMS convention is quoted because the knobs are
variations of the SAME selection, not independent error sources.

## 5. Bootstrap

Paired Poisson bootstrap (weights ~ Poisson(1) per event, same weights for both estimators —
`methodGainPostfix.C`), fixed seed `TRandom3 rng(20260610)`, 2000 replicas, 68% interval from
the 16/84 percentiles. The CI brackets the TAIL-SENSITIVE difference (convention 1b); its
central value is printed beside it in the generated table.

## 6. MIXED same-shower uncertainty

The 1.04 ± 0.05 uncertainty is the ENERGY-PERIOD SCATTER (SEM of the five per-energy ratios),
which dominates the per-point statistical errors — `makeMixedKillshotFigure.C`. Run-jackknife
spread < 0.01 (gate log).

## 7. Known small caveats (documented, no gated number affected)

- `tebSigma`'s robust FALLBACK debias constant 0.9546 is the single-pass 2.5σ truncation factor;
  the iterated fixed point is ~0.938, so the fallback is ~−2% biased for a pure Gaussian.
  Published resolutions use the Gaussian-core fit (fallback only trips for pathological
  samples). Correcting the constant requires a full gate rerun — post-submission item.
- `fullFiducialCheck` fixes r = 3.0 mm at all energies (systematics-nominal protocol);
  it matches production exactly at ≥125 GeV including the headline (see its AUDIT addendum).

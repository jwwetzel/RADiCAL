# MIXED Corner-Map Audit Note (pre-analysis) — GATE 1, data-only half

Date: 2026-06-09. Written BEFORE the discriminant script, per gated-workflow protocol.

## Question
The Paper-2 kill-shot (in-event DSB1/LuAG ratio ≈ 0.99) rests on the corner material map
**NE+SW = DSB1, NW+SE = LuAG:Ce**. That map was inferred FROM THE DATA via brightness/clipping
(NE/SW slam the 820 mV clip, NW/SE do not) — but brightness is the thesis variable, so a referee can
call the inference circular. This audit defines a brightness-independent pulse-SHAPE discriminant and
the confirmation/contradiction criteria. (The physical-logbook confirmation remains the other half of
GATE 1 and is a USER action; this analysis stands on its own either way.)

## Files used
- `data/2023/reduced/MIXED/{50,75,100,125,150}GeV.root` ('rad' tree; 39-branch schema verified
  2026-06-09 — see `papers/memory_dataset_inventory.md`).
- References (same schema, same reduction batch): `data/2023/reduced/DSB1/*.root` (8 known-DSB1 ends)
  and `data/2023/reduced/LUAG/*.root` (8 known-LuAG ends) — these provide LABELED template
  distributions of every discriminant, so MIXED corners are classified against measured references,
  not against assumptions.

## Channels present
8 capillary ends per build, kCap order: 0=NW-D 1=NE-D 2=SE-D 3=SW-D 4=NW-U 5=NE-U 6=SE-U 7=SW-U.
A corner's two ends (i, i+4) share the capillary — both must classify identically (built-in cross-check).

## Shape-sensitive branches available (verified to exist)
- `hg_charge[8]` (integrated HG pulse), `hg_peak[8]` (HG amplitude), `hg_tot[8]` (time over threshold)
- `lg_charge[8]`, `lg_peak[8]` (LG — never clips; 1 GS/s)
- Crossing times at multiple fractions: `hg_cfd05/hg_cfd(=20)/hg_cfd30/hg_cfd50[8]`, `hg_led[8]`
- `hg_saturated[8]` (clip flag — lets us select UNCLIPPED events for HG-shape work)

## Physics lever
DSB1 fluorescence decay τ ≈ 3.5 ns vs LuAG:Ce τ ≈ 70 ns (×20). Slower decay ⇒ for the same collected
charge, a LuAG pulse is much wider and lower-peaked, and its leading edge develops more slowly.

## Candidate discriminants (amplitude-normalized ⇒ brightness-independent to first order)
1. **D1 = lg_charge / lg_peak** (effective pulse width, LG). PRIMARY. LG never clips, so this works on
   ALL events including the bright DSB1 corners; the ×20 τ difference dwarfs the 1 GS/s sampling
   coarseness. Pure ratio ⇒ scale-invariant under brightness.
2. **D2 = hg_cfd50 − hg_cfd05 per end** (5%→50% rise time, HG @5 GS/s). CFD fractions are
   amplitude-independent by construction. Restricted to UNCLIPPED events (`!hg_saturated`) so the
   recorded peak is the true peak — note this biases the DSB1-corner sample toward dim events, which
   is acceptable because the discriminant is shape, not brightness.
3. **D3 = hg_charge / hg_peak** (HG width) on unclipped events — cross-check of D1 at 5 GS/s.
4. **D4 = hg_tot / hg_peak** — secondary; threshold-dependent, used only qualitatively.

## Circularity assessment
D1–D3 are ratios or fraction-to-fraction times: invariant under multiplying the pulse by a constant.
Brightness enters only through second-order effects (noise fraction on dim pulses, sampling phase).
The classification therefore does NOT reuse the brightness ordering that produced the original
inference. Residual coupling to acknowledge in print: clipping forces D2/D3 onto unclipped subsets;
D1 has no such restriction and is therefore the primary.

## Protocol (pre-registered)
1. Build labeled reference distributions of D1 (and D2, D3) per end from the pure DSB1 and pure LUAG
   builds at each energy (and combined). Quantify separation (means, widths, overlap).
2. Compute the same distributions for the 8 MIXED ends at each energy.
3. Classify each MIXED end against the references (nearest-template / likelihood-ratio of the D1
   distribution); record per-end classification margin.
4. Cross-checks: (a) the two ends of one capillary must agree; (b) classification must be
   energy-independent; (c) D2/D3 must agree with D1 where computable.

## Confirmation / contradiction criteria (pre-registered)
- **CONFIRMS** the map if: NE-D, NE-U, SW-D, SW-U classify DSB1-like AND NW-D, NW-U, SE-D, SE-U
  classify LuAG-like, at every energy, with per-end margins ≫ the reference widths, and D1/D2/D3 agree.
- **CONTRADICTS** if any end stably classifies opposite to the map (then: kill-shot re-derived with
  the corrected map — the analysis is map-parametric; rerun `mixedSeparate.C` + h2h with swapped indices).
- **INCONCLUSIVE** if reference distributions overlap too much to classify (not expected: τ ratio ×20),
  in which case the logbook becomes the only resolution and the kill-shot carries an explicit caveat.

## Outputs
Script: `papers/scripts/mixed_corner_map/cornerDiscriminant.C` (written AFTER this audit).
Figure: `papers/figures/mixed_corner_map/corner_discriminant.png` (reference distributions + the 8
MIXED ends overlaid + per-end classification table). Log: `corner_discriminant.log` next to the script.

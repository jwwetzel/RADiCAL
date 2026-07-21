# AUDIT — full-fiducial companion check (fullFiducialCheck.C)

> **Provenance honesty note: this audit was written POST HOC on 2026-07-21** (workspace
> world-class audit, Phase A1), documenting a verification script that was created and first
> run on 2026-06-10 during the prior-record audit (commit `164880b`). It is NOT a
> pre-registered gate audit like the others in `papers/scripts/*/AUDIT.md`; it records,
> after the fact, what the script does, why it exists, and what its committed log shows.
> Nothing here was backdated.

## Purpose

The claims law (`papers/memory_claims_and_forbidden_language.md`) requires the brightest-1000
headline to be quoted WITH its full-fiducial companion. The law's original companion value
("27–29 ps") was pre-fix-era and unverified; this script produced the verified replacement.

## What the script does (read-out only; no new selection or estimator)

`fullFiducialCheck.C` re-implements the locked production chain exactly as
`systematicsPostfix.C`'s nominal: DSB1 reduced data, wire-chamber fiducial r < 3.0 mm, MCP
window, HG ≥ 20 mV, in-event consistency veto W = 2.0 ns, srCFD (`hg_lgcfd`) times, post-fix
`tebSigma` width. It reports, per energy: (a) σ over ALL qualifying fiducial events
("sigma_full") and (b) σ of the brightest-1000 ("sigma_K1000") as a chain cross-check.

## Pass condition (implicit at creation, stated here)

sigma_K1000(150) must reproduce the gated headline 25.7 ps exactly — proving the chain is the
locked production chain. It does (25.7). sigma_full(150) is then the claims-law companion.

## Result (committed log: `full_fiducial_result.log`, rerun 2026-07-21, identical to 2026-06-10)

| E (GeV) | N_fid | sigma_full (ps) | sigma_K1000 (ps) |
|---|---|---|---|
| 25 | 41280 | 64.6 | 45.6 |
| 50 | 81488 | 53.7 | 33.5 |
| 75 | 76619 | 50.0 | 30.8 |
| 100 | 98994 | 50.1 | 26.4 |
| 125 | 77201 | 49.8 | 25.7 |
| 150 | 145423 | 50.5 | 25.7 |

**Companion number: σ_full(150 GeV) ≈ 50 ps (50.5).** Consumed by the manuscript (abstract,
Secs. 4/5.3, Conclusions) and by the amended claims-law Timing bullet; the retired "27–29 ps"
is documented in the law's own amendment note.

## Reproduce

```bash
source setup.sh
root -l -b -q 'papers/scripts/full_fiducial_check/fullFiducialCheck.C+' | tee full_fiducial_result.log
```
Deterministic (no RNG); requires `data/2023/reduced/DSB1/`.

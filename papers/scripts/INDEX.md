# papers/scripts/ — the gate index and runbook

One directory per **published result**: the generator macro, its pre-registered
`AUDIT.md`, and the committed result log are the evidence chain for every number in
the manuscript. This index records, for each gate, what it generates, what it clobbers,
and how to run it. Provenance overview: `../../ANALYSIS_GUIDE.md`; gate history and
retirements: `../memory_analysis_gates.md`.

## ⚠ The runbook rule

**Gate macros regenerate their committed outputs IN PLACE.** Running one is a
deliberate act of *verification*, never routine analysis:

1. Start from a **clean tree**, at the **repo root**, after `source setup.sh`.
2. Run the gate (safest: `tools/run_gate.sh <gate> --yes`, which refuses to run
   without data; or `tools/repro.sh` for the automated headline chain).
3. `git status --porcelain --untracked-files=no` must come back **clean** — byte-identical
   regeneration is the reproduction. The single allowed residual is
   `papers/timing/figs/thesis_postfix.pdf`, which differs only in its embedded PDF
   `/CreationDate`+`/ModDate` (environment drift; `tools/repro.sh` checks and restores it).
4. A dirty diff is a **FINDING** — investigate, never commit a silent "repair."

Never edit the macros' bytes, the committed `*_result.{log,txt}` evidence files, or an
existing `AUDIT.md` (append-only, dated addenda only). The committed logs are tracked
despite the `papers/**/*.log` gitignore pattern (added with `git add -f` — see the
comment in `.gitignore`).

All entry macros take **no arguments** and are deterministic (seeded bootstrap where
randomness is involved; no wall-clock randomness). Data-dependent gates need
`data/2023/reduced/` (~13 GB, not in git — `data/2023/README.md`; integrity:
`shasum -a 256 -c data/2023/MANIFEST.sha256`).

## The gates

| Gate | Entry macro(s) | Gates (paper item) | Regenerates | Data |
|---|---|---|---|---|
| `timing_fit_summary` | `timingFitSummary.C` | **money plot** Fig. 6 + Table 1: four-build σ_t(E) fits, headline **25.7 ± 0.6 ps**, floor 18.8 ± 0.8 | `papers/tables/timing_fit_summary_2026-06-09.md`, `papers/figures/timing_fit_summary/timing_fit_summary.png`, `papers/timing/figs/thesis_postfix.{png,pdf}` | all 4 builds |
| `full_fiducial_check` | `fullFiducialCheck.C` | ≈50 ps full-fiducial companion number (abstract, Secs. 4/5.3) | stdout only → committed `full_fiducial_result.log` (AUDIT honestly post-hoc; fixed r=3.0 mm at all E — see its addendum) | DSB1 |
| `systematics_postfix` | `systematicsPostfix.C` | Table 2 selection-systematics budget + Fig. A.11 | `papers/timing/tab_systematics.tex`, `papers/tables/systematics_postfix_2026-06-09.md`, `papers/figures/systematics_postfix/systematics_stability.png`, `papers/timing/figs/systematics.png` | all 4 builds |
| `method_gain_postfix` | `methodGainPostfix.C` | §5.3 method-gain recompute (1.3 ps core / 4.0 ps tail-sensitive; Fig. 7 material, currently quarantined behind TODO-P2 in the tex) | `papers/figures/method_gain_postfix/method_gain_postfix.png`, `papers/tables/method_gain_postfix_2026-06-09.md`, `papers/timing/tab_methods.tex` | DSB1 |
| `depth_dial` | `depthDial.C` (+ `depthDialDiag.C` diagnostics) | Fig. 9 depth dial: −33.6 ± 2.9 ps per e-fold (GATE 3) | `papers/figures/depth_dial/depth_dial.png`, `depth_dial_diag.png` | DSB1 (+3 cross-checks) |
| `mixed_killshot_bootstrap` | `mixedKillshotBootstrap.C`, then `makeMixedKillshotFigure.C` | Fig. 8 same-shower MIXED control: σ ratio **1.04 ± 0.05** (GATE 6) | `papers/figures/mixed_killshot_bootstrap/{killshot_bootstrap.png, mixed_h2h_corrected.png/.pdf, …CAPTION.txt}` | MIXED |
| `mixed_corner_map` | `cornerDiscriminant.C` | GATE 1 (data half): pulse-shape confirmation of the MIXED corner map (logbook half = open USER item) | `papers/figures/mixed_corner_map/corner_discriminant.png` | MIXED+refs |
| `position_reconciliation` | `positionReconcile.C` | GATE 2: what mm claim the companion energy/position paper is allowed | `papers/figures/position_reconciliation/position_reconcile.png` | DSB1 |
| `satellite_removal` | `satelliteRemoval.C` | Fig. B.12: satellite/tail suppression on identical events (4.60 → 3.30%) | `papers/figures/satellite_removal/satellite_removal.{png,pdf}` + caption | DSB1 |
| `apparatus_composite` | `apparatus_composite.tex` (tectonic, **data-free**) | Fig. 1 apparatus composite (provenance: `papers/timing/PARENT_APPARATUS_FIGURE_INSERTION_AUDIT_2026-06-10.md`) | `apparatus_composite.pdf` build artifact (untracked; canonical copy in `papers/timing/figs/`) | none |

Committed evidence log per gate: `papers/scripts/<gate>/*_result.{log,txt}` (stdout
captures). Figure bytes consumed by the manuscript are recorded in
`papers/timing/figs/MANIFEST.md` (drift check: `tools/check_figure_manifest.sh`).

## Canonical invocation

```sh
source setup.sh                     # once per shell, at the repo root
tools/run_gate.sh --list            # see gate names
tools/run_gate.sh timing_fit_summary --yes
git status --porcelain --untracked-files=no   # must be clean (or the PDF-timestamp residual)
```

Or the whole headline chain in one command with the verdict automated:

```sh
tools/repro.sh --check-data
```

Last full-chain byte-identical reproduction: 2026-07-21 (`CODE_AUDIT_2026-07-21.md`,
reproduction runner E4 — ~93 s of machine time, sole residual the PDF timestamp).

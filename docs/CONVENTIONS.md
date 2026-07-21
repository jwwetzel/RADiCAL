# CONVENTIONS — names, taxonomy, and house rules

(REORG_PLAN Phase B4, recorded 2026-07-21. Statistics conventions — width definitions,
error scaling, bootstrap, systematics totals — live in `STATS_CONVENTIONS.md`; the
claims wording law is `../papers/memory_claims_and_forbidden_language.md`.)

## Names that must never drift

- **srCFD ≡ `hg_lgcfd` ≡ `RadView::kLGCFD`.** Paper name ≡ branch name ≡ code enum —
  the saturation-recovered CFD (threshold = 0.15 × LG-predicted true peak). The branch
  is *never* renamed (binding REORG rule); the bridge is stated at all three definition
  sites and in `lib/physics/ESTIMATOR.md`.
- **Builds:** `DSB1`, `LUAG`, `MIXED`, `TENERGY` (dir + config-JSON spellings; prose
  uses DSB1, LuAG:Ce, MIXED, TENERGY). LYSO:Ce + W are common to all builds — **the WLS
  capillary is the variable** (claims law; never "crystal differs").
- **Canonical end order** (Schema arrays index 0–7): `NW-D, NE-D, SE-D, SW-D, NW-U,
  NE-U, SE-U, SW-U` — Down 0–3, Up 4–7. `eventDWUP` hardcodes this; `BuildConfig`
  guards it at load.
- **Energies:** 25 (DSB1 only), 50, 75, 100, 125, 150 GeV; files `<E>GeV.root`.
- **Timing sources** (per-end branches): `hg_cfd05`/`cfd03`/…/`cfd50` (fixed-fraction
  CFD), `hg_led` (leading edge), `hg_lgcfd` (srCFD), `hg_tot` (retired). Adopted
  per-regime: srCFD (DSB1, MIXED), LED (LuAG, TENERGY), `cfd05` diagnostic.

## Estimator taxonomy

- `eventDWUP` — the per-event reference-free kernel: (DW−UP)/2 + 2 ns in-event median
  veto. Cloned deliberately into gate macros (KERNEL-CLONE comments name each delta).
- `tebSigma` — production width: robust truncated-RMS seed → Gaussian-core fit, with
  tail/agreement guards. ("Robust fallback" = the seed used directly.)
- `timingBrightestK` (+ `kLGCFD`, K = 1000) — **the production headline chain**.
- `timingBestBin` — **RETIRED** (historical Method A); marked at its definition.
Full walk: `lib/physics/ESTIMATOR.md`. Every numeric threshold lives in
`SelectionCuts.h` with a WHY comment — new cuts go there, never inline.

## Where outputs go

- Studies write figures under `figures/` and scratch under `output/` — both
  regenerable, never load-bearing.
- Gate outputs are **committed evidence**: figures under `papers/figures/<gate>/`,
  tables under `papers/tables/*.md` (+ `papers/timing/tab_*.tex`), stdout logs beside
  the macro (tracked via `git add -f` against the `papers/**/*.log` ignore). Manuscript
  figure bytes: `papers/timing/figs/` + `MANIFEST.md`.

## House rules (the short list; binding source: WORKSPACE_REORG_PLAN §4)

1. **Append-only records:** AUDIT.md files, committed result logs, the campaign
   snapshot, memory files. Corrections are dated addenda, never edits; never backdate.
2. **Gate reruns are verification acts** — clean tree, repo root, `git diff` clean
   afterwards (sole residual: `thesis_postfix.pdf` /CreationDate). Dirty = finding.
3. **Never regenerate** current `papers/timing/figs/` bytes (incl. frozen `dist.png`,
   `optimization.png`) or renumber gated numbers pre-submission; never modify
   `papers/scripts/*/` macro bytes (guards live in `tools/run_gate.sh` instead).
4. **Committed PDF = tectonic build** (`cd papers/timing && tectonic radical_timing.tex`).
5. **Durable facts live in-repo** (docs/, data/2023/metadata/, papers/memory_*.md);
   assistant memory is working context only and is never committed.
6. **Never commit ROOT data files**; figures > 2000 px get a `sips -Z` /tmp copy
   before viewing.
7. DEPRECATED/HISTORICAL banners follow the house pattern (banner at top, pointer to
   the replacement); breadcrumbs are never deleted.

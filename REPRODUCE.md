# REPRODUCE — from a fresh clone to the published numbers

The reproduction claim this repo makes: every published number and figure regenerates
**byte-identically** from the committed macros over the reduced data (sole known
residual: one PDF's embedded timestamp). This file is the shortest honest path to
checking that yourself. The full provenance map is `ANALYSIS_GUIDE.md`.

## 1. Environment

- **ROOT** — verified byte-identical under **ROOT 6.40.00 (macosxarm64)** on
  2026-07-21 (`CODE_AUDIT_2026-07-21.md`, reproduction runner E4). ROOT's PNG/text
  output has proven byte-stable across sessions here; other ROOT versions or platforms
  may reproduce numerically but not byte-for-byte. **Honesty line: Linux is untested.**
- **tectonic** — only for rebuilding the manuscript PDF (`cd papers/timing && tectonic
  radical_timing.tex`; document class + orcidlink are vendored). Not needed for the
  numeric chain.
- No other dependencies: the analysis is plain ROOT macros; `source setup.sh` once per
  shell (sets `ROOT_INCLUDE_PATH` + `RAD_DATA`).

## 2. Data

The 21 reduced ntuples (~13 GB; 39-branch `rad` trees) are **not in git**:
`data/2023/reduced/<BUILD>/<E>GeV.root` for DSB1 (25–150 GeV), LUAG, MIXED, TENERGY
(50–150 GeV). See `data/2023/README.md` for the download source
*(public data deposit + fetch script pending — the venue decision (Zenodo vs CERN Open
Data) is an open pre-submission item; until then the share link in that README is the
source)*. Verify integrity before trusting any reproduction:

```sh
shasum -a 256 -c data/2023/MANIFEST.sha256      # 29 files, all must say OK
```

## 3. The one-command reproduction

```sh
source setup.sh
tools/repro.sh --check-data
```

This verifies the data manifest, runs the data-presence gate (`reduce/verify.C`), the
headline gate (`timing_fit_summary` → 25.7 ± 0.6 ps @ 150 GeV, DSB1, brightest-1000,
srCFD; a = 203 ± 6 → 440 ± 18 ps·√GeV; floor 18.8 ± 0.8 ps), and the full-fiducial
companion (≈50 ps), then applies the byte-identity verdict via `git diff`.

**Expected residual — exactly one:** `papers/timing/figs/thesis_postfix.pdf` differs
only in its embedded PDF `/CreationDate` + `/ModDate` (ROOT's PDF writer embeds the
wall clock; 18 bytes, characterized in `CODE_AUDIT_2026-07-21.md`). `repro.sh` verifies
the difference is *only* that and restores the file. Everything else — the table
markdown and every PNG — must be **byte-identical**. Any other dirty file is a finding.

## 4. The remaining gates

Each published figure/number has its own gate directory under `papers/scripts/<gate>/`
(macro + pre-registered AUDIT + committed log). Runbook and per-gate clobber lists:
`papers/scripts/INDEX.md`; guarded single-gate runs: `tools/run_gate.sh <gate> --yes`.
All gates are deterministic (seeded bootstrap; no wall-clock randomness); ~93 s of
machine time reproduces the full numeric chain (2026-07-21 measurement).

## 5. What reproduction does NOT cover

- **Raw → reduced**: the reduction ran on the Argon cluster (`reduce/hpc/`); raw data
  (~2 GB/run) is mirrored there. `reduce/validateReduce.C` is the committed
  bit-identity evidence for the reducer itself.
- **Frozen figures**: `dist.png` and `optimization.png` have no live generator
  (frozen-as-circulated; `papers/timing/figs/MANIFEST.md`) — re-verification of the
  optimization curves is a tracked pre-submission item.
- **The manuscript PDF**: rebuildable with tectonic; its committed bytes are the
  canonical build.

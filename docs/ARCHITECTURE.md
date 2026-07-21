# ARCHITECTURE — how the code is laid out and why

(REORG_PLAN Phase B4, recorded 2026-07-21. The *analysis* waterfall — raw → reduced →
gates → tables → manuscript — is `../ANALYSIS_GUIDE.md`; this file is the *code* map.)

## The four lib/ domains

`setup.sh` puts all four on `ROOT_INCLUDE_PATH`, so every macro uses bare `#include`s.

```
lib/waveform/   raw-signal layer (only needed for raw-waveform work)
  ChannelConfig.h    flat DRS4 array indexing + 2023 board defaults (DRS0 5 GS/s timing,
                     DRS1 1 GS/s energy). NOTE: the AUTHORITATIVE channel map is the
                     per-build JSON via BuildConfig.h; this header's arithmetic is
                     kept consistent with it (see header note).
  DRS4Calibration.h  DT5742 time-base recovery + stop-cell correction
  WaveformUtils.h    pulse extraction (ExtractPulse, Pulse, kNoTime)

lib/io/         paths, config, schema, reader
  DataPaths.h        THE path resolver ($RAD_DATA/data/<year>/...; no legacy fallbacks)
  FigPaths.h         campaign/year-namespaced figure output paths
  MiniJson.h         dependency-free JSON reader for build configs
  BuildConfig.h      build JSON → resolved flat DRS4 offsets + .hglg calibration
                     sidecar; enforces the D4U4 end-order invariant eventDWUP needs
  Schema.h           the ONE canonical 39-branch reduced-ntuple schema ('rad' tree)
  RadView.h          format-agnostic role-resolved reader over a 'rad' tree

lib/physics/    selection + timing science
  SelectionCuts.h    single source of truth for every cut threshold + TimingFiducialR(E)
  RadTiming.h        the timing kernel and estimators (see ESTIMATOR.md beside it)

lib/viz/        plotting
  RADiCALStyle.h     gStyle + shared color constants
  PlotUtils.h        FitGaussCore, DrawFitOverlay, ScanRunCenters
```

Intended layering: waveform → io → physics → viz. One **known inversion**: `RadTiming.h`
includes `PlotUtils.h` (viz) to reach `FitGaussCore` — documented in
`CODE_AUDIT_2026-07-21.md` (deep item D4); fixing it pre-submission was deliberately
deferred to keep gate-macro provenance frozen.

## The pipelines that use it

- **Reduction** (`reduce/`): `Reducer.C` is the one config-driven raw→reduced maker
  (driven by `data/2023/configs/<BUILD>.json`; HG↔LG calibration sidecar `<BUILD>.hglg`
  from `calibHGLG.C`). `reduceRun.C` is the per-run HPC driver (`reduce/hpc/` = Argon
  SGE array + merge scripts); `validateReduce.C` proved bit-identity against the legacy
  reducer; `verify.C` is the post-reduction cleanliness gate. Reduction normally runs on
  Argon, not locally.
- **Exploration** (`analyze/studies/`): ~120 single-purpose macros over the reduced
  ntuples — the preserved research trail (`analyze/README.md`, `studies/INDEX.md`).
- **Publication** (`papers/scripts/<gate>/`): one dir per published result — macro +
  pre-registered AUDIT.md + committed log (`papers/scripts/INDEX.md` is the runbook).
  Gate macros deliberately contain their own copies of the kernel (KERNEL-CLONE
  comments name each delta) so the published evidence chain has no hidden moving parts.
- **Verification** (`tools/`): `repro.sh` (headline chain, one command),
  `run_gate.sh` (guarded single-gate runs), `check_figure_manifest.sh` (figure drift),
  `pull_argon_data.sh` (raw mirror).

## Data resolution

Everything resolves through `DataPaths.h`: `$RAD_DATA` (default: repo root, set by
`setup.sh`) → `data/<year>/{raw,reduced,configs,metadata}`, year switched by
`$RAD_YEAR` (default 2023). Reduced files: `data/2023/reduced/<BUILD>/<E>GeV.root`,
39-branch `rad` tree (`lib/io/Schema.h`), ~13 GB total, not in git — integrity manifest
`data/2023/MANIFEST.sha256`. Configs + metadata are committed.

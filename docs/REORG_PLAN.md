# RADiCAL Workspace Reorganization Plan  —  ✅ COMPLETED

**Status:** DONE. Executed in four commits (data → code split → figures → README)
after all 26 reduced files verified clean. The end-state is documented in the
top-level `README.md`; this file is kept as a historical record of the plan.
(Some `datasets/`→`data/` strings below were auto-rewritten by the migration's
own path sweep — read them as the original intent.)

## Locked decisions
1. **Data home:** rename `data/` → **`data/`** (one home; absorbs the legacy
   `Data/` + top-level `reduced/` symlink shims).
2. **Disposition:** **archive everything** — nothing deleted. One-off macros →
   `analyze/studies/`; `legacy/` → `archive/`. All recoverable.
3. **Library layout:** **domain subfolders** under `lib/`.
4. **Scope:** **full reorg, phases A–F**, each a single compile-tested commit.

## Target tree
```
RADiCAL/
├── README.md            ← entry point: reduce in 3 cmds, analyze in 1
├── lib/                 ← shared library (all .h). ROOT_INCLUDE_PATH spans the 4 subdirs.
│   ├── waveform/  WaveformUtils.h  ChannelConfig.h  DRS4Calibration.h
│   ├── io/        MiniJson.h  BuildConfig.h  Schema.h  RadView.h  DataPaths.h
│   ├── physics/   SelectionCuts.h  RadTiming.h
│   └── viz/       PlotUtils.h  RADiCALStyle.h
├── reduce/              ← raw → reduced (CURRENT pipeline only)
│   ├── Reducer.C  reduceRun.C  calibHGLG.C  validateReduce.C
│   └── hpc/            (SGE scripts, env.sh, manifests — from reduce/hpc)
├── analyze/            ← reduced → results
│   ├── timingHeadline.C  sigmaT.C  timingLadder.C  slopeVsE.C   (curated)
│   └── studies/        (every other .C — kept, runnable, quarantined)
├── data/               ← ALL data, per year
│   └── 2023/{raw, reduced, configs, metadata}/  README.md  MANIFEST.csv
│       configs/  = <BUILD>.json + <BUILD>.hglg   (was data/2023/configs)
│       metadata/ = channel_map.yaml dataset.yaml runs.csv  (was data/2023/config)
├── figures/            ← all generated PNG/PDF (was radcore/figs, scattered)
├── report/             ← GitHub Pages site — UNTOUCHED
├── docs/               ← this plan, design notes; apparatus doc lives in memory
└── archive/            ← legacy/ + decommissioned scratch
```

## Key architecture rules (DRY / AI-proof)
- **One include root:** `ROOT_INCLUDE_PATH=lib/waveform:lib/io:lib/physics:lib/viz`.
  Bare `#include "X.h"` lines are unchanged — only files move + env updates.
- **One data resolver:** `lib/io/DataPaths.h` points only at `data/<year>/…`;
  delete the legacy `Data/`, `reduced/`, `output/` fallbacks.
- **Dependency direction:** `reduce/` and `analyze/` depend on `lib/`, never reverse.
  `lib/` knows nothing about specific studies.

## Macro routing (96 .C → destinations)
**reduce/** (current pipeline): Reducer.C, reduceRun.C, calibHGLG.C, validateReduce.C
**analyze/** (curated headline): timingHeadline.C, sigmaT.C, timingLadder.C, slopeVsE.C
**analyze/studies/** (everything else, archived-runnable):
  - from radcore/: absShowerTime, absTime, cfdInterp, cfdScan, checkCenter,
    depthCorr, edgeSlope, fixCFD, gateDSB1, headlineLG, hgLgClean, hgLgPlot,
    inferMixed, lgCFD, mcpJitterCanon, pbLeakage, perRun150, satLinearity,
    sigmaVsAmp, slopeVsThresh, testConfig, threshScan, trSync
  - from Analysis/ (all 65 legacy/study .C): alignmentAnalysis, analyzeResolution,
    averageWaveforms, cfdFractionDSB1, channelCalibration, channelCombinationScan,
    channelIntegrity, chargeProfiles, compareConfigsPlot, compareEnergies,
    configBestBin{,DSB1,HGLG}, configCapDiag, configResolution{,DSB1,Full},
    configSystematics, desaturateCFD, discoverChannels, discoverReduced,
    drs4Diagnostics, drs4TimeBase, dsb1OOSandCorr, edgeMechanism, elbow*,
    etype{Char,Energy}, fiducialTimingScan, fourDdemo, harvestResults,
    hgLgLinearity, idealUniform, investigatePbGlass, layer{1..6}Summary,
    lightYieldTiming, mcpJitter, mixedHeadToHead, moduleCenter, peaksHGLG125,
    positionCorrection, qualityPlots, showerLocalization, slewTest,
    systematicUncertainties, tenergyClean, timing*, transverseMaps,
    uniformityScan, walkCorrTest, wireChamberResolution
  - **legacy reducers** processRun.C, reduceRaw.C → archive/ (superseded by reduce/)
**lib/** headers: per domain table above.
**archive/**: legacy/ contents, Analysis/__pycache__, scratch outputs.

## Header → lib/ domain map
- waveform/: WaveformUtils.h, ChannelConfig.h, DRS4Calibration.h
- io/:       MiniJson.h, BuildConfig.h, Schema.h, RadView.h, DataPaths.h
- physics/:  SelectionCuts.h, RadTiming.h
- viz/:      PlotUtils.h, RADiCALStyle.h

## Path-coupling checklist (edit together, one commit per phase)
- [ ] `lib/io/DataPaths.h` — drop legacy fallbacks; canonical = `data/<year>/…`.
- [ ] `ROOT_INCLUDE_PATH` everywhere: reduce/hpc/env.sh, README cmds, macro
      header comments (`ROOT_INCLUDE_PATH=lib/waveform:lib/io:lib/physics:lib/viz` → the 4 lib subdirs).
- [ ] Hardcoded `data/2023/...` strings in macros (e.g. calibHGLG.C,
      slopeVsE.C, hgLgPlot.C, hgLgClean.C use `data/2023/configs/%s.json` and
      `data/2023/raw/...`) → `data/2023/...`. grep before moving.
- [ ] `radcore/figs/` output paths in macros → `figures/`.
- [ ] reduce/hpc scripts: REC_DIR / RAD_WORK / repo paths, manifest locations.
- [ ] rsync-home target: `data/2023/reduced/` → `data/2023/reduced/`.
- [ ] .gitignore: Data/, reduced/, data/*/raw, data/*/reduced,
      data/*/output, output, legacy/* → data/*/raw, data/*/reduced,
      figures cache, archive/.

## Execution order (each = compile-tested commit)
- **A — Data:** `git mv datasets data`; rename `data/2023/config`→`metadata`;
  drop empty `figures/`,`output/`; delete `Data/` + top-level `reduced/` shims;
  simplify DataPaths.h; fix .gitignore + rsync target. (Lowest risk — data already
  unified via symlinks.)
- **B — Library:** create `lib/{waveform,io,physics,viz}`; `git mv` headers; set
  `ROOT_INCLUDE_PATH` to the 4 subdirs; compile smoke-test (one reduce + one
  analyze macro). Include lines unchanged.
- **C — Code split:** create `reduce/`, `analyze/`, `analyze/studies/`; `git mv`
  per routing table; `git mv reduce/hpc reduce/hpc`; fix hardcoded data paths.
- **D — Figures:** `git mv radcore/figs/* figures/`; repoint macro output paths.
- **E — Legacy:** `git mv legacy archive/legacy`; move dead scratch in.
- **F — README:** top-level quickstart so a newcomer can drive (reduce + analyze
  one-liners, where data lives, how to add a year/build).

## Do-not-break
- `report/` GitHub Pages site (report/Summary/*.png, output/Summary/
  results.json) — keep building. If results.json path moves, update its writer.
- Argon access: port 40 (`ssh -p 40`), ROOT via `setup_root`. Cannot reach Argon
  from local Mac.
- Push to main only after each phase compiles; gate first push on data being home.

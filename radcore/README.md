# radcore — the RADiCAL processing hub

One centralized, config-driven pipeline for **any** RADiCAL module built from
**4 capillaries → 8 SiPM readouts → 2 DRS4 boards**. Drop in new data with a
different channel configuration, analyze it, and produce papers — all in one repo,
so previous results are always one path away.

## The idea

A RADiCAL build differs from another only in three declarative things:
the **channel map** (which DRS4 channel is each capillary end / MCP / wire-chamber
plane), the **per-capillary material & role** (timing vs energy), and the
**run list**. Everything downstream — waveform extraction, the reduced-ntuple
schema, and the analyses — is the same. So a build is fully described by one
JSON file, and the code is build- and year-agnostic.

```
datasets/
  registry.csv                 every build, every year, its config + status
  <year>/
    configs/<BUILD>.json        ← the only thing you write for new data
    raw/        RUN*.root        raw DRS4 waveforms (download from CERNBox)
    reduced/<BUILD>/<E>GeV.root  canonical reduced ntuple (one schema)
radcore/
  MiniJson.h        dependency-free JSON reader
  BuildConfig.h     loads a config JSON → resolved DRS4 offsets, geometry, runs
  Schema.h          the ONE canonical reduced-ntuple schema (RadEvent + branches)
  Reducer.C         (P2) config-driven raw → canonical reducer
  testConfig.C      gate: JSON config == legacy hardcoded ChannelConfig.h map
```

## Config JSON (one file = one build)

`channel_map.ends[]` lists the 8 readouts in canonical order
(NW-D, NE-D, SE-D, SW-D, NW-U, NE-U, SE-U, SW-U); each gives hardware
coordinates `[drs, grp, ch]` for its high-gain (timing) and low-gain (energy)
DRS4 channels and which MCP it references. `capillaries[]` carries the per-corner
material and role. `geometry`, `mcp_ref_jitter_ps`, and `runs` complete it.
See `datasets/2023/configs/DSB1.json`.

`BuildConfig::Load()` converts the `[drs,grp,ch]` triples to the same flat
offsets the legacy `Analysis/ChannelConfig.h` used:

```
chanOff(drs,grp,ch) = (1024*9*2)*drs + (1024*9)*grp + 1024*ch
timeOff(drs,grp)    = (1024*2)*drs   + 1024*grp
```

`radcore/testConfig.C` proves the loaded offsets equal the hardcoded `kCap[]`
map exactly (Phase-1 gate — passing).

## Canonical schema (one ntuple format)

`Schema.h` defines `rad::RadEvent`, a strict **superset** of the two legacy
formats: the rich per-capillary timing/energy features from `processRun.C`
(`hg_cfd03/05/10/20/30/50`, `hg_led`, `hg_tot`, `hg_charge`, `hg_ped_rms`,
`hg_saturated/spike`, `lg_peak/charge`, `sum_lg/sum_pb`, `pb_peak`, `wc_*`,
`mcp1/2`, `stopcell`) **plus** the generic per-slot arrays from `reduceRaw.C`
(`s_peak[36]`, `s_cfd05[36]`, `s_charge[36]`). One file feeds every analysis.

The per-capillary `[8]` arrays are in a fixed canonical end order, so index `i`
always means the same physical corner/end across builds and years.

## To add new data

1. Put raw files in `datasets/<year>/raw/`.
2. Write `datasets/<year>/configs/<BUILD>.json` (channel map + materials + runs).
3. Run the reducer → `datasets/<year>/reduced/<BUILD>/<E>GeV.root`.
4. Run analyses (build/year-agnostic) and harvest results into the report.

## Running the reducer

Locally, over the run list in the config (one or all energies):

```
ROOT_INCLUDE_PATH=radcore:Analysis \
  root -l -b -q 'radcore/Reducer.C+("datasets/2023/configs/DSB1.json", 150)'   # one energy
  root -l -b -q 'radcore/Reducer.C+("datasets/2023/configs/DSB1.json")'        # all energies
```

On the cluster (full statistics — all manifest runs), the SGE array job calls
the per-run entry point once per raw file, then `hadd` merges per energy:

```
root -l -b -q 'radcore/Reducer.C+'   # compile once
root -l -b -q 'ReduceFile("datasets/2023/configs/DSB1.json","<raw>.root",150,"<out>.root")'
hadd datasets/2023/reduced/DSB1/150GeV.root <byrun>/*.root
```

`ReduceFile` is config-driven and bit-identical to the legacy `processRun.C`
(`radcore/validateReduce.C` gate: worst |old−new| = 0 over a full run), so a
full re-reduction reproduces the published numbers by construction. Swap the
`reduceRaw.C` call in `Analysis/hpc/sge_reduce.sh` for this `ReduceFile` call to
upgrade the cluster pipeline to the canonical schema.

The canonical reader (`rad::RadEvent::ConnectBranches`) also reads the *existing*
legacy DSB1 merges transparently (it binds the old `mcp_peak`→`mcp1_peak`), so
analyses can move to the schema before the cluster re-reduction lands.

## Status

- **P1 (done):** schema + config loader + DSB1 config + gate.
- **P2:** unified config-driven reducer.
- **P3:** re-reduce DSB1 + validate it reproduces the published numbers.
- **P4:** analyses + report onto the schema.
- **P5:** roll out LUAG/MIXED/TENERGY + 2024–2026.
- **P6:** retire the legacy two-format split.

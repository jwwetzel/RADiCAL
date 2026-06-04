# RADiCAL 2023 — CERN SPS H2, May 2023 (reference campaign)

The original, fully-analysed RADiCAL campaign and the template for all later years.

- **Config:** `config/channel_map.yaml`, `config/dataset.yaml`, `config/runs.csv`
  (this campaign's wiring, metadata, and run list).
- **Live analysis & reports** currently sit at the repo root, not under here:
  - framework + 2023 channel map: `Analysis/` (`ChannelConfig.h`, `WaveformUtils.h`, …)
  - six-layer technical report: `report/`
  - publication-style paper: `paper.html`
  - raw data: repo-root `/Data`; reduced capillary-config files: repo-root `/reduced`
- **Headline:** σ_t = 27.4 ps @ 150 GeV; σ_E/E = 13.9 %; four capillary builds
  (DSB1 / LuAG / mixed / TENERGY) — see `Analysis/CAPILLARY_COMPARISON.md`.

This directory is kept as the **canonical worked example**: `config/channel_map.yaml`
is the accurate 2023 map, so new campaigns can diff against it and change only what
actually differs. 2023's bulk data has intentionally not been moved here (it would
break the live GitHub Pages site) — see `datasets/README.md` → *Status / migration*.

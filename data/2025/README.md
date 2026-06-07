# RADiCAL &lt;YEAR&gt; — &lt;campaign&gt;  (TEMPLATE)

Starting point for a new test-beam campaign. To use:

1. `cp -r data/_TEMPLATE data/<year>` (already done for 2024–2026).
2. Edit `config/dataset.yaml` — year, facility, beam, energies, capillary builds,
   geometry centre, any per-campaign cut overrides.
3. Edit `config/channel_map.yaml` — **verify every `[drs,grp,ch]`** against this
   campaign's wiring (the main thing that changes between years).
4. Fill `config/runs.csv` from the run logbook.
5. Drop a few sample raw runs into `raw/` and smoke-test the reduction locally
   before submitting the full set to the cluster.

`raw/`, `reduced/`, `output/` are gitignored (data payloads); `config/` and
`figures/` are tracked. See `data/README.md` for the full scheme and the
cluster workflow.

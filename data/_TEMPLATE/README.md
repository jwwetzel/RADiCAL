# RADiCAL &lt;YEAR&gt; — &lt;campaign&gt;  (TEMPLATE)

Starting point for a new test-beam campaign. Full playbook: `docs/PUBLICATION_ENGINE.md`
→ "New dataset". To use:

1. `cp -r data/_TEMPLATE data/<year>` (already done for 2024–2026).
2. Fill the **logbook** in `metadata/` — `dataset.yaml` (year, facility, beam, energies,
   builds, geometry centre, any per-campaign cut overrides), `channel_map.yaml`, and
   `runs.csv` from the run logbook. These are the human record (decorative — not read by code).
3. **Write the build config(s)** — the only files the code reads. For each build create
   `configs/<BUILD>.json` (copy a 2023 one, e.g. `data/2023/configs/DSB1.json`); **verify
   every `[drs,grp,ch]`** against this campaign's wiring and `hg_sat_mV` against its DC-offset
   positive rail (measure from the `hg_peak` pile-up). Optionally regenerate the `.hglg`
   calibration sidecar with `reduce/calibHGLG.C` on clean low-E data.
4. Drop a few sample raw runs into `raw/` and smoke-test locally:
   `RAD_YEAR=<year> root -l -b -q 'analyze/sigmaT.C+("DSB1",150)'` (no recompile to switch years).
5. Reduce on the cluster (`reduce/hpc/`, point it at this campaign's manifest).

`raw/`, `reduced/` are gitignored data payloads; `configs/` and `metadata/` are tracked.
The analysis switches campaign via `export RAD_YEAR=<year>`. See `data/README.md` for the
full scheme.

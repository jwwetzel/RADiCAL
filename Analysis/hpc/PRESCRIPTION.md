# Prescription — reduce ALL capillary configs on Argon, then analyze locally

Goal: turn the full ~0.5 TB raw test-beam dataset into **small reduced ntuples,
organized by capillary configuration**, that you copy home for the timing- and
energy-resolution analysis.

The four configs in the logbook (electron runs):

| config label (logbook col)   | ~runs | what it is                         |
|------------------------------|-------|------------------------------------|
| `DSB1`                       | ~330  | the already-validated config       |
| `LUAG`                       | ~485  | LuAG capillaries                   |
| `2xDSB1,2xLuAG`              | ~252  | mixed DSB1/LuAG                    |
| `3xDSB1,1xEnergy`            | ~210  | timing capillaries + 1 energy      |

**Key design choice — the reduction is config-AGNOSTIC.** `reduceRaw.C` does not
need a channel map: for every one of the 36 DRS slots it stores the pulse
features (`s_peak`, `s_cfd05`, `s_charge`) and it precomputes the
config-invariant wire-chamber track (`x_trk,y_trk,wc_ok`) and both MCPs. So the
*same* job reduces every config. Which slot is a timing vs energy capillary, the
corner/depth geometry, and the resolutions are all worked out **locally** on the
reduced files (where iterating is cheap), using `discoverChannels.C` +
the analysis chain. `discoverChannels.C` is validated on DSB1 — it re-derives the
known channel map from raw waveforms with no prior knowledge.

Reduced size: ~12 MB/run (160x smaller than raw); a merged per-energy file is a
few tens of MB; a whole config is ~1 GB. The full reduced set is ~15 GB.

---

## 0.  One-time setup (Argon)

```bash
git clone <repo> RADiCAL && cd RADiCAL
source Analysis/hpc/env.sh         # sets REC_DIR (raw), RAD_WORK (scratch), queue
#   check/override: REC_DIR=/Shared/lss_yonel/jwwetzel/RADiCAL_CERN_May2023/rec/rec
qsub Analysis/hpc/compile.sh       # prebuilds processRun + reduceRaw + discoverChannels .so
```

## 1.  Build a manifest per config (from the logbook)

`build_manifest.py` filters the logbook to electron runs of one capillary config.
Run it once per config (quote labels with commas):

```bash
LB=/path/to/logbook.csv
python3 Analysis/hpc/build_manifest.py "$LB" --capillary "DSB1"            -o Analysis/hpc/manifest_DSB1.csv
python3 Analysis/hpc/build_manifest.py "$LB" --capillary "LUAG"            -o Analysis/hpc/manifest_LUAG.csv
python3 Analysis/hpc/build_manifest.py "$LB" --capillary "2xDSB1,2xLuAG"   -o Analysis/hpc/manifest_MIXED.csv
python3 Analysis/hpc/build_manifest.py "$LB" --capillary "3xDSB1,1xEnergy" -o Analysis/hpc/manifest_TENERGY.csv
```
Each prints a per-energy run-count summary to stderr — sanity-check it.
(Add `--all-bias` to include the 41.25 V scan; default keeps the single 42.25 V gain.)

## 2.  Per config: resolve files → reduce (array) → merge

Do this block once per config (shown for LuAG):

```bash
export RAD_CONFIG=LUAG
export MANIFEST=$RAD_REPO/Analysis/hpc/manifest_LUAG.csv
export RAD_TASKLIST=$RAD_WORK/tasks_LUAG.txt

bash Analysis/hpc/discover_tasklist.sh                 # manifest runs -> raw file paths
N=$(wc -l < "$RAD_TASKLIST")
qsub -hold_jid rad_compile -t 1-$N Analysis/hpc/sge_reduce.sh   # embarrassingly parallel

# ...after the array finishes (qstat shows no rad_reduce):
bash Analysis/hpc/merge_reduced.sh                     # -> reduced/LUAG/<E>GeV.root
```
Repeat with `RAD_CONFIG`/`MANIFEST`/`RAD_TASKLIST` set to `DSB1`, `MIXED`,
`TENERGY`. (DSB1 here is the cross-check: its reduced file must reproduce the
known result.)

## 3.  (optional but recommended) confirm each config's slot layout

```bash
# pick any one raw run of the config and characterise its 36 slots
root -l -b -q 'Analysis/discoverChannels.C+("'"$REC_DIR"'/RUN<aLuAGrun>.root",3000)'
```
This prints, per slot, occupancy / amplitude / edge-speed / inferred role —
revealing the config's capillary count and timing-vs-energy split. (The same
role classification is repeated locally on the reduced files in step 5, using
the stored `s_charge/s_peak` ratio as the fast-vs-slow discriminator.)

## 4.  Copy reduced files home

```bash
rsync -av <argon>:.../RADiCAL/reduced/  ./reduced/
```
~1 GB/config. Layout: `reduced/<CONFIG>/<E>GeV.root`.

## 5.  Local analysis (per config) — done here

For each config, on the reduced ntuple:
1. **Channel map** — classify slots (timing / energy / MCP / WC) by pulse
   signature (`discoverChannels` logic ported to the reduced tree).
2. **Geometry** — assign corner/depth from energy-vs-WC-position correlation.
3. **Timing resolution** — (DW-UP)/2 corner estimator + per-channel, vs energy.
4. **Energy resolution** — Sigma(energy-capillary charge), sigma_E/E = a/sqrt(E) (+) b.
5. Cross-config comparison (DSB1 vs LuAG vs mixed vs timing+energy).

---

### Reduced ntuple schema (`rad` tree, per event)

| branch | meaning |
|---|---|
| `run, event, beam_energy` | provenance |
| `wc_ok, x_trk, y_trk` | wire-chamber track [mm] (config-invariant) |
| `mcp1_peak/time, mcp2_peak/time` | the two MCP references [mV, ns] |
| `stopcell[4]` | DRS4 stop cell per group (for timebase correction) |
| `s_peak[36]` | per-slot amplitude [mV] |
| `s_cfd05[36]` | per-slot CFD-5% crossing time [ns] |
| `s_charge[36]` | per-slot pulse integral |

slot index `s` → `drs=s/18, grp=(s/9)%2, ch=s%9`; `ch==8` is the DT5742
trigger/sync. DRS0 = fast/timing side, DRS1 G0 = slow/energy side, DRS1 G1 = WC.

> **Imported into the repository 2026-07-21 (workspace audit, Phase A4).** This file was
> maintained as analysis-assistant working memory during the 2026 publication campaign and is
> imported verbatim so the published workspace is self-contained: a reader needs no external
> notes. Line references cite repo files as of the import commit. Source: memory/dataset_2023.md.

# RADiCAL — 2023 Dataset Instance (CERN SPS H2, May 2023)

**Tier 2 PER-DATASET (CERN SPS, May 2023). Concept: `detector_invariant.md` · Methods: `methods_invariant.md`.**

This file is the source of every 2023 NUMBER: the four builds, the full authoritative channel map,
the SiPM bias, the measured detector characteristics, the selection-cut working points, and the
data/run manifest. The invariant detector CONCEPT (W/LYSO:Ce shashlik design, capillary roles,
(DW−UP)/2 estimator philosophy) lives in `detector_invariant.md`; the analysis METHODS (reduction
schema, pulse extraction, timing estimators, floor decomposition) live in `methods_invariant.md`.
Do not duplicate those here — point to the sibling. `RAD_YEAR=2023` selects this campaign in
`DataPaths.h`; figures land under `figures/2023/...` via `FigPaths.h`.

*Global index: `MEMORY.md`. Concept: `detector_invariant.md` · Methods: `methods_invariant.md` · Paper thesis +
figure index: `paper_timing.md` / `paper_energy_position.md`. Chronological log: `log_2023.md`.*

---

## Experiment one-liner

RADiCAL ultra-compact W/LYSO:Ce shashlik EM-calorimeter module tested at the **CERN SPS H2 line in
May 2023**, electron beam **25/50/75/100/125/150 GeV** (25 GeV steps), four WLS-capillary builds
(**DSB1, LuAG, MIXED, TENERGY**) — LYSO:Ce scintillator + W absorber common to all; the variable is
the WLS capillary fill. Parent/reference paper: **arXiv:2401.01747 = Perez-Lara, Wetzel, Akgun et al.
(RADiCAL), NIM A 1068 (2024) 169737** (DSB1 build validated against it).

---

## The four builds

All four builds share an **identical 8+8 DRS4 wiring / channel map**; the build label describes only
the capillary fill, not the wiring [CAPILLARY_COMPARISON:13-15,43-44; LUAG/MIXED/TENERGY.json `_note`].
Run-range → build map is verified against the manifests.

| Build | WLS capillary fill | Corner material map | Light level | Energies (GeV) | Run range / N | Bias |
|---|---|---|---|---|---|---|
| **DSB1** (baseline) | 4 × DSB1 organic, all T-type | all 4 corners DSB1 | bright | **25**,50,75,100,125,150 | 1034–1333 / 252 | 42.25 V |
| **LUAG** | 4 × LuAG:Ce ceramic, T-type | all 4 corners LuAG | dim (~3× less) | 50,75,100,125,150 | 2186–2427 / 212 | 42.25 V |
| **MIXED** | 2 × DSB1 + 2 × LuAG, all T-type | NE+SW = DSB1 (bright), NW+SE = LuAG (dim) | mixed in-module | 50,75,100,125,150 | 2787–3002 / 190 | 42.25 V |
| **TENERGY** | 3 × DSB1 T-type + 1 × E-type energy | NE,SE,SW = DSB1; **NW = E-type** | 3 bright + 1 dim NW | 50,75,100,125,150 | 2530–2754 / 162 | 42.25 V |

Sources: DSB1.json/LUAG.json/MIXED.json/TENERGY.json:3-5; manifest_*.csv; CAPILLARY_COMPARISON:31-36;
PAPER1_DRAFT:85-87; PAPER2_DRAFT:72-74; _superseded/figure_catalog.md:87-90.

**Per-energy run counts (manifests):**
- DSB1 = 34/37/36/35/35/75 (25/50/75/100/125/150 GeV)
- LUAG = 35/35/33/41/68 (50/75/100/125/150 GeV)
- MIXED = 34/24/43/45/44 (50/75/100/125/150 GeV)
- TENERGY = 34/36/33/28/31 (50/75/100/125/150 GeV)

Notes:
- **Only DSB1** was run at 25 GeV; LUAG/MIXED/TENERGY start at 50 GeV [CAPILLARY_COMPARISON:33-38].
- MIXED corner material is data-inferred (radcore/inferMixed.C, matching per-corner brightness @150 GeV)
  and flagged CONFIRM-against-logbook [MIXED.json:3-5,33-38].
- MIXED valid analysis era is restricted to runs **2912–3044** (early runs had bugged wire chambers);
  100/125/150 restricted to the valid 2975–3044 era (237k/283k/264k ev) [_superseded/figure_catalog.md:90; AC §0].
- Larger logbook counts exist (electron/DONE/nominal-bias rows differ from analyzed): DSB1 330/252,
  LUAG 485/212, MIXED 252/190, TENERGY 210/162; full logbook = 1278 rows [CAPILLARY_COMPARISON:29-36].
- A one-off LUAG run **2247** was 3 LuAG timing + 1 LuAG energy cap; **no four-E-type configuration
  exists** in the 2023 run [PAPER2_DRAFT:154-157].
- LuAG build NW capillary type (T-type vs full-length E-type) is unconfirmed (logbook "Added full LuAG
  Cap (Energy Cap)" motivates the check) — see VERIFY [CAPILLARY_COMPARISON:231-233, confidence=medium].

---

## SiPM bias

| Property | Value | Source |
|---|---|---|
| Bias voltage (nominal) | **v_up = v_down = 42.25 V**, identical upstream/downstream | manifest_{dsb1,LUAG,MIXED,TENERGY}.csv (816 rows, all 42.25/42.25) |

The 42.25 V nominal holds uniformly across **all 252 DSB1, 212 LuAG, 190 MIXED, 162 TENERGY runs**
[manifest_*.csv]. Settled after a SiPM bias scan aimed at finding the maximum bias before saturation
at 125 GeV; an earlier/lower test point of **41.25 V** (runs 100/1001–1033) was used before raising to
42.25 V and is dropped from the nominal manifest so all energies share one gain [logbook.csv lines
119,144-148,186; log_2023.md:370]. (SiPM model = Hamamatsu HDR2; active area / pixel
size not stated in the parent paper — see `detector_invariant.md`.)

---

## Full channel map — AUTHORITATIVE (logbook spreadsheet, J. Wetzel 2026-06)

Two CAEN **DT5742** switched-capacitor digitizers (DRS4-based), **DRS0** and **DRS1**, ran at
**different rates** (author-confirmed via the logbook channel map): **DRS0 @ 5 GS/s** = the *timing*
board (8 HG capillary channels + both MCP copies + **all 4 Pb-glass bars** + the 1×1/2×2 trigger
scintillators); **DRS1 @ 1 GS/s** = energy + tracking (8 LG capillary channels + the wire chambers).
DRS0 → **0.200 ns/cell**; DRS1 → **1 ns/cell**. 1024 cells/channel; each board = 2 inner groups × 9
channels (8 signal + 1 trigger). [ChannelConfig.h; DRS4Calibration.h; DSB1.json:6-15.]

`idx` = flat "Index in data" = `chanOff(drs,grp,ch) = 18432·drs + 9216·grp + 1024·ch`; validated
exactly against `ChannelConfig.h`. Time-axis offsets: `time_offset(drs,grp) = 2048·drs + 1024·grp`
→ kT_D0G0=0, D0G1=1024, D1G0=2048, D1G1=3072. Detector labels: `D/U` = Down/Up end (= DW/UP),
`H/L` = High/Low gain; e.g. "NE UH1" = NE-Up high-gain.

**DRS0 — 5 GS/s (timing board)**
| grp | idx | detector | role |
|---|---|---|---|
| G0 | 0 | SW DH1 | HG SW-Down |
| G0 | 1024 | NW DH2 | HG NW-Down |
| G0 | 2048 | NE DH3 | HG NE-Down |
| G0 | 3072 | SE DH4 | HG SE-Down |
| G0 | 4096 | NE UH1 | HG NE-Up |
| G0 | 5120 | NW UH2 | HG NW-Up |
| G0 | 6144 | SE UH3 | HG SE-Up |
| G0 | 7168 | MCP | **MCP1** |
| G0 | 8192 | Trigger | TRS0 (DRS trigger ch) |
| G1 | 9216 | SW UH4 | HG SW-Up (only HG in G1) |
| G1 | 10240 | PB2 SW | Pb-glass bar |
| G1 | 11264 | PB1 NW | Pb-glass bar |
| G1 | 12288 | PB3 NE | Pb-glass bar |
| G1 | 13312 | PB4 SE | Pb-glass bar |
| G1 | 14336 | 1×1 cm | trigger scint A1 |
| G1 | 15360 | 2×2 cm | trigger scint A2 |
| G1 | 16384 | MCP | **MCP2** |
| G1 | 17408 | Trigger | TRS0 |

**DRS1 — 1 GS/s (energy + tracking board)**
| grp | idx | detector | role |
|---|---|---|---|
| G0 | 18432 | SW DL5 | LG SW-Down |
| G0 | 19456 | NW DL6 | LG NW-Down |
| G0 | 20480 | NE DL7 | LG NE-Down |
| G0 | 21504 | SE DL8 | LG SE-Down |
| G0 | 22528 | NE UL5 | LG NE-Up |
| G0 | 23552 | NW UL6 | LG NW-Up |
| G0 | 24576 | SE UL7 | LG SE-Up |
| G0 | 25600 | SW UL8 | LG SW-Up |
| G0 | 26624 | Trigger | TRS0 |
| G1 | 27648 | WC-X-Anode | wire chamber X anode |
| G1 | 28672 | WC-X-Right | WC X right (delay line) |
| G1 | 29696 | WC-X-Left | WC X left (delay line) |
| G1 | 30720 | WC-Y-Down | WC Y down (delay line) |
| G1 | 31744 | WC-Y-Anode | wire chamber Y anode |
| G1 | 32768 | WC-Y-Up | WC Y up (delay line) |
| G1 | 33792 | 1×1 cm | scintillator |
| G1 | 34816 | Empty | — |
| G1 | 35840 | Trigger | TRS0 |

**kCap channel order (8 ends):** 0 = NW-D, 1 = NE-D, 2 = SE-D, 3 = SW-D, 4 = NW-U, 5 = NE-U,
6 = SE-U, 7 = SW-U [ChannelConfig.h]. NW-Up (kCap index 4) is the consistently weakest timing channel;
dropping it improves a channel combination by ~10 ps [log_2023.md:84,188].

**HG capillary map [drs,grp,ch]:** NW-D=[0,0,1], NE-D=[0,0,2], SE-D=[0,0,3], SW-D=[0,0,0],
NW-U=[0,0,5], NE-U=[0,0,4], SE-U=[0,0,6], **SW-U=[0,1,0]** (the only HG channel on DRS0 G1).
**LG map:** all eight on DRS1 G0 — NW-D=[1,0,1], NE-D=[1,0,2], SE-D=[1,0,3], SW-D=[1,0,0],
NW-U=[1,0,5], NE-U=[1,0,4], SE-U=[1,0,6], SW-U=[1,0,7] [ChannelConfig.h:81-98; DSB1.json:20-27].
**Consequence:** corners **NW, NE, SE are entirely in DRS0 G0** (both ends + MCP1); **SW is split**
(SW-D in G0, SW-U in G1 with MCP2). Each channel references the MCP copy in its **own** group
(caps 0–6 → MCP1; SW-U → MCP2) so same-group differential estimators cancel inter-group jitter
exactly — a deliberate design choice. (Why the (DW−UP)/2 floor sits far below the 71 ps inter-chip
floor — see `methods_invariant.md`.)

**MCP reference copies.** Single MCP-PMT (Hamamatsu R3809U-50) passively **split into both DT5742
groups** → MCP1 (DRS0 G0 ch7, `kMCP1=chanOff(0,0,7)=7168`) and MCP2 (DRS0 G1 ch7,
`kMCP2=chanOff(0,1,7)=16384`); corr(MCP1,MCP2)=1.000, equal amplitude ~236–237 mV. Configs quote
`mcp_ref_jitter_ps = 73` (BuildConfig default 73.0) [DSB1.json:54; BuildConfig.h:83].

**TR0 trigger.** Digitized DT5742 trigger input (TR0 / "TRS0") on the 9th channel of every group —
**TR0a=[0,0,8] (idx 8192), TR0b=[0,1,8] (idx 17408)** on the timing board; `(t_TR0a − t_TR0b)` gives
the inter-group mezzanine timing offset [Reducer.C:133-140; DSB1.json:31-32].

**Wire chambers.** Delay-line WC on DRS1 G1: X = Anode[1,1,0](27648)+Right[1,1,1]+Left[1,1,2];
Y = Down[1,1,3]+Anode[1,1,4](31744)+Up[1,1,5]. `x = kWC_Scale·(t_R−t_L)`, `y = kWC_Scale·(t_D−t_U)`;
**kWC_Scale = 7/36 mm/ns (~0.1944)** (`scale_mm_per_ns_frac [7,36]`); 50% CFD; effective WC position
res ~1 mm (kWC_resBin = 1.0 mm/bin; self-cal σ_x ≈ 3.6 mm). `wc_ok` requires all four planes,
peakTime>0, ≥20 mV (kWC_minPeak) [ChannelConfig.h:50-59; Reducer.C:110-121; SelectionCuts.h:196-200].

**Pb-glass leakage tagger.** 4 Pb-glass bars stacked 2×2 (each ~4×4×40 cm³), on DRS0 G1 ch1–4:
PB2-SW=[0,1,1](10240), PB1-NW=[0,1,2](11264), PB3-NE=[0,1,3](12288), PB4-SE=[0,1,4](13312)
[ChannelConfig.h:101-110; DSB1.json:34].

**Stop cell.** Trigger asynchronous to the free-running DRS4 domino wave → stop cell rotates uniformly
over [0,1023] (RMS 296–298 ≈ 1024/√12). Recovered by `drs4::FindStopCell`, stored `stopcell[4]` =
[D0G0, D0G1, D1G0, D1G1]. Beam lands at a fixed window position: MCP1 crossing ~91.3 ns ± 3 ns,
readout cell ~456 [DRS4Calibration.h:14-17,57-77; Reducer.C:104-107; log_2023.md:256-258].

---

## Measured detector characteristics (2023 instance)

**HG clipping & saturation.**
- HG hard-clips at **~820 mV** positive (clip onset ~818 mV, modal pile-up ~836 mV; `hg_sat_mV = 820`
  in the 2023 configs) — **DT5742 channel saturation, not SiPM-pixel saturation** — because the
  acquisition window was DC-offset down to capture the negative afterpulse, lowering the effective
  positive clip from the nominal ~950 mV [AC:341-344; DSB1.json:13].
- HG saturated (clipped) fraction by energy: 25 GeV=**4%**, 50=**71%**, 75=**90%**, 100–150=**93–94%**
  [AC:344].
- The legacy `hg_saturated` flag (950 mV threshold) **never fires (0.000%)** for 2023 data — wrong for
  this campaign [AC:344].
- Mean HG peak vs energy (DSB1): ~498/760/802/811/808/800 mV at 25/50/75/100/125/150 GeV — HG
  compresses above ~50 GeV [log_2023.md:56].

**Gains, edge slopes, noise.**
- LG is linear to 150 GeV (~4.8× rise 25→150 GeV); HG peak only plateaus (~1.6×, e.g. 465→788 mV)
  [AC:345-347].
- HG rising-edge dV/dt keeps rising with energy: **102→170→214→244→261→279 mV/ns** (25→150 GeV) — no
  SiPM saturation [AC:348].
- Edge slope at the 20 mV foot = 59→109 mV/ns; at the 400 mV steep level = 172→612 mV/ns (25→150).
  Single-channel slew σ = ped_rms/(dV/dt): foot ~22→12 ps, steep ~7.5→2.1 ps [AC:477-478].
- `hg_lgcfd` (LG-recovered) edge ~3.7× steeper than cfd05-on-clipped (e.g. 381 vs 102 mV/ns at
  ~125 GeV) [AC:230].
- HG pedestal/noise RMS `hg_ped_rms ≈ 1.31 mV` (slew budget); pedestal-RMS floor target ~5 mV/channel
  [AC:361; log_2023.md:499].

**Capillary light levels (amplitude discriminators).**
- MIXED/TENERGY brightness tags: bright ~630 mV = DSB1, dim ~300 mV = LuAG [CAPILLARY_COMPARISON:52-54].
- MIXED LG: LuAG corners (NW,SE) collect ~3× less than DSB1 corners (NE,SW): ~66–85 mV vs ~171–231 mV
  at 150 GeV [AC:92-93].
- TENERGY: NW E-type ~4–5× dimmer than DSB1 caps (LG ~110 mV vs ~480 mV at 150 GeV)
  [CAPILLARY_COMPARISON:215-216].
- E-type vs T-type LG amplitude slope: 1.6 vs 6.7 mV/GeV [PAPER2_DRAFT:127-129; CALOR2026 p.12].
- HG clip level on DSB1 ~820 mV [CAPILLARY_COMPARISON:198].
- Per-channel HG–LG correlation healthy across all builds: ρ=0.73–0.95 over all 32 channels
  (DSB1 0.81–0.95; MIXED 0.73–0.90 — MIXED LG is **not** "toast") [AC:89-95].

**DRS4 / reference jitter & timebase.**
- σ(MCP1−MCP2)/√2 ≈ **71–74 ps**, flat with energy and amplitude-independent — *not* the MCP intrinsic
  jitter (~10 ps common-mode, cancels); it is the relative **inter-group DRS4 domino-phase / stop-cell
  jitter** between the two free-running mezzanines [AC:378-392]. (Mechanism & why it never touches the
  headline: see `methods_invariant.md`.)
- DRS0 cell-width RMS < 1 ps (**0.84–0.90 ps** = float32 precision on 0.2 ns) → negligible for timing.
  (DRS1 LG/WC groups show ~4.8 ps but are not used for timing.) [AC:351; log_2023.md:309.]

**Optical / shower.**
- WLS-fiber light velocity **v_fiber ≈ 200 mm/ns** (σ_z = v_fiber·σ_t) [AC:372].
- LG energy sum droops ~9–20% by 150 GeV (longitudinal leakage out the back of the 25 X₀ module);
  PbGlass linearizes residual non-linearity from 9.2% to 2.2% [AC:349-350].
- Shashlik LG chain is AC-coupled: prompt pulse → ~35% balancing undershoot (~150–180 ns) → clean
  baseline by ~500 ns, no ringing (a few channels e.g. SE-U show a larger ~350 ns secondary lobe)
  [log_2023.md:491-493].

**Beam-spot / momentum.**
- Beam-spot transverse spread ~**2.9 mm** [PAPER2_DRAFT:111-112; _superseded/figure_catalog.md:195].
- Momentum spread: CERN H2 literature value ~1.0% typical (max Δp/p ±2%, energy-independent); not
  measured run-specific — see VERIFY [log_2023.md:344].

**Per-build headline σ_t (2023 numbers; method/interpretation → `methods_invariant.md` +
`_superseded/radical_apparatus_conclusions.md`).**
- Reference timing floor [arXiv:2401.01747]: **σ_t = 255.58/√E ⊕ 17.52 ps (~27 ps @150 GeV)**, cfd05
  on the clipped pulse; documents a ~0.2 ns satellite [Eq.2]. This work reproduces it within ~1.6 ps.
- Per-build (DW−UP)/2 σ_t @150 GeV (OOS best-bin): **DSB1 29.0 ps** (27.4 official Gaussian-core),
  **LuAG 37.4, MIXED 38.2, TENERGY 35.7** [PAPER1_DRAFT:134-136; CAPILLARY_COMPARISON:171].
- Brightest-1000 robust @150: **DSB1 25.3–25.4 ps** (lgcfd), **TENERGY ~28, MIXED ~39–40, LUAG ~41–42**
  [AC:176-184]. Monotonic ladder DSB1 (150<125<100<75<50<25): 25.3/25.6/28.9/30.7/35.6/43.4.
- Photostat-form (1/√E) floors per build: **DSB1 19.5±1.1, LuAG 19.8±5.6, MIXED 34.6±2.4,
  TENERGY 21.6±2.7 ps** (CONSISTENT with the published 17.5; do not revise it) [AC:40].
- In-event per-capillary intrinsic timing (MIXED pairwise): DSB1 ≈ LuAG ≈ ~200 ps/cap @150 GeV — no
  material difference [PAPER1_DRAFT:149-153].
- Energy resolution σ_E/E ≈ 52.04%/√E ⊕ 31.62%/E ⊕ 9.31% [arXiv:2401.01747 Eq.1]; shower-max σ_E/E
  ~14% @150 GeV, WLS-species-independent (13.5–15%) [PAPER2_DRAFT:85-89].
- Transverse position from 4-corner light: residual ~1.5 mm (1.52 x, 1.44 y) @150 GeV, validated vs
  wire chamber [PAPER2_DRAFT:110-112].
- Per-config systematics 1.1–1.7 ps (DSB1 1.3, LuAG 1.1, MIXED 1.7), dominated by fiducial-radius
  variation [PAPER1_DRAFT:176-180].

---

## Selection-cut working points (2023 instance)

Single source of truth: `lib/physics/SelectionCuts.h`. (Cut rationale / how they enter the estimators
→ `methods_invariant.md`.)

| Cut | Value | Source |
|---|---|---|
| MCP timing window | 200 mV (min) – 750 mV (max sat cut) | SelectionCuts.h:151-155 |
| MCP energy-only min | 50 mV | SelectionCuts.h:155 |
| HG per-channel timing min | **20 mV** (= LED thresh, ~4× the ~5 mV HG noise) | SelectionCuts.h:166; AC:590 |
| HG saturation | hg_sat_mV = **820 mV** (config); WaveformUtils/SelectionCuts default 950 mV | DSB1.json:13; WaveformUtils.h:124; SelectionCuts.h:168 |
| LG energy min | 5 mV | SelectionCuts.h:173 |
| PbGlass per-PMT min | 5 mV | SelectionCuts.h:178 |
| Containment cut | reject if sum_pb > 0.30·sum_lg (only when sum_lg > 300 mV) | SelectionCuts.h:160,187 |
| Module center (WC coords) | **kCalo_x0=6.6, y0=4.7 mm** (data-derived from shashlik edges; was 6.5/4.5) | DSB1.json:51; SelectionCuts.h:116-120 |
| Energy fiducial radius | **r = 2.0 mm** (tight) | SelectionCuts.h:125-146 |
| Timing fiducial radius | **r = 3.0 mm** nominal; OOS-optimized `TimingFiducialR(E)` = 2.5 mm for ≤100 GeV, 3.0 mm for 125 & 150 GeV | SelectionCuts.h:125-146; log_2023.md:813-814 |

Fiducial cuts are centered on the per-run data-derived signal-weighted beam centroid (ΣHG of 8 caps),
not the nominal center; ΣLG cross-check agrees <0.05 mm [SelectionCuts.h:116-120;
log_2023.md:656-658].

**Beam reference runs (DSB1, per energy):** 25 GeV=RUN1211, 50=RUN1148, 75=RUN1112, 100=RUN1075,
125=RUN1034, 150=RUN1258–1261 chained (~8.5 GB, ~120k events) [ChannelConfig.h:128-140; DSB1.json:56-63].

---

## Data / run manifest

- **Statistics:** ≥10⁶ triggers per energy step (parent paper); analyzed per-build totals DSB1=252,
  LUAG=212, MIXED=190, TENERGY=162 [arXiv:2401.01747 §5; manifest_*.csv].
- **Reduced files:** all 4 builds re-reduced + verified clean (21 files, `hg_lgcfd` + tr0a/b_time
  present); output `data/2023/reduced/<BUILD>/<energy>GeV.root`, ZSTD level 5 [Reducer.C:44-46,240-248;
  AC §0]. (Reduction schema + the `rad` tree branch list → `methods_invariant.md`.)
- **Manifests / logbook:** `reduce/hpc/logbook.csv` (1278 rows), `reduce/hpc/manifest_{dsb1,LUAG,MIXED,
  TENERGY}.csv` (816 rows total at 42.25/42.25). Build run ranges per the four-builds table above.
- **MIXED restricted era:** runs 2912–3044 (early runs had bugged wire chambers); 100/125/150 use the
  valid 2975–3044 subset (237k/283k/264k ev) [AC §0; _superseded/figure_catalog.md:90].
- **Local-vs-Argon raw layout.** Argon access via `p -40` (personal ssh wrapper, NOT
  `ssh hawkid@argon...`); login nodes argon-login-5 / argon-itf-login-3. Non-DSB1 raw:
  `/Shared/lss_yonel/jwwetzel/RADiCAL_CERN_May2023/rec/rec/RUN####.root` (no energy suffix);
  **DSB1 local raw:** `RUN####_<E>_GeV.root` (energy suffix). `RAD_WORK =
  /nfsscratch/jwwetzel/RADiCAL_CERN_May2023/analysis_work` [AC:540-543].
- Reduce pipeline runs on Argon from `reduce/hpc/` (compile.sh → submit_reduce.sh → merge); Argon needs
  a `git pull` to get the post-2026 layout. TENERGY: exclude `RUN2721.root` (truncated/Zombie) before
  reducing [AC:496].

---

## VERIFY (2023-specific low-confidence / unsourced / inconsistent — do not present as established)

- **MIXED corner material map** (NE+SW=DSB1, NW+SE=LuAG) is data-inferred (inferMixed.C; `MIXED.json`
  carries `_material_map_source` = "data-inferred… CONFIRM against logbook"); the inferred NW match is
  the weakest [MIXED.json:3-5].
- **LuAG build NW capillary type** (T-type vs full-length E-type) is unconfirmed; logbook line "Added
  full LuAG Cap (Energy Cap)" motivates the check [CAPILLARY_COMPARISON:231-233, confidence=medium].
- **Beam momentum spread ~1.0%** is the CERN H2 literature value, not measured run-specific
  [log_2023.md:344, confidence=medium].
- **MCP single-tube passive split into both DT5742 groups** (MCP1/MCP2, corr=1.000) is from local
  memory; not spelled out in the parent paper [confidence=medium].
- **[RESOLVED 2026-06 — author-confirmed]** The **timing (HG) board ran at 5 GS/s** (0.2 ns/cell); the
  parent paper states the DT5742 model but not the rate. The offline indices map cleanly to the two
  physical boards: DRS0 (idx 0–17408) = 5 GS/s (8 HG caps + both MCP copies + 4 Pb-glass bars + 1×1/2×2
  trigger scints) and DRS1 (idx 18432–35840) = 1 GS/s (8 LG caps + wire chambers). `chanOff` index math
  validated row-by-row against the logbook map (channel-map table above). No residual.

(Invariant / parent-paper-only specs — capillary physical dimensions, DSB1 optical peaks, SiPM active
area, Molière radius reconciliation — are NOT 2023-instance items; they live in `detector_invariant.md`.)

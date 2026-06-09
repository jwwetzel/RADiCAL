<!-- Synthesized 2026-06 from a 5-way repo audit (existing plans, multi-year data, code reuse,
figure pipeline, memory architecture). The organizing reference for scaling this repo into a
repeatable multi-dataset publication engine. A recommendation to execute incrementally; see the
migration path (§5) for ordered, low-risk steps. -->

# RADiCAL Publication Engine — Organizing Architecture

A single, concrete plan for turning the 2023-shaped repo into a repeatable multi-dataset publication engine. This **builds on decisions already made** (the config-driven reducer, `BuildConfig`/`RadView`/`RadTiming`/`DataPaths`, the per-year `data/<year>/` skeleton, the layered HTML report, the 2-paper split, the reorg's dependency direction) and only adds the missing connective tissue. It is a recommendation to approve, not an execution.

---

## 1. Organizing principle

**Encode the INVARIANT once; parameterize the VARIABLE everywhere.**

There are exactly two kinds of knowledge in this program:

- **INVARIANT** — true of *any* RADiCAL W/LYSO:Ce shashlik module in *any* beam-time: the detector concept (WLS capillary is the variable, not the crystal), the T-type/E-type roles, dual-gain HG-timing/LG-energy with deliberate HG saturation, the timing method `(DW−UP)/2` BestMinus, brightest-K, Gaussian-core σ, OOS run-folding, the floor-physics decomposition (floor = ~1 X₀ shower-depth fluctuation, σ_z = v_fiber·σ_t), the DRS4/DT5742 quirks and which cancel in same-group differentials, the code abstractions, the house plot style, and the parent paper arXiv:2401.01747. This is the **reusable core** — written, hardened, and read-only thereafter.
- **VARIABLE** — specific to one campaign or one paper: the builds table (DSB1/LUAG/MIXED/TENERGY), run ranges, per-energy counts, the literal channel-map indices, the measured numbers (820 mV clip, 71–74 ps jitter, 25.3/27.4 ps headlines), every figure PNG, and each paper's thesis/narrative/systematics state.

The repo **already tiers the VARIABLE** for data (`data/2023..2026/`, `_TEMPLATE/`, `registry.csv`) and papers (`papers/timing/`, `papers/energy_position/`). The failure mode today is that the INVARIANT is *not* separated out — it is smeared across 2023-specific code (macros hardcode `data/2023/...`), 2023-specific figure paths (no year namespace), and 2023-specific memory (`experimental_setup.md` carries run numbers; `radical_apparatus_conclusions.md` is 95% Paper-1 state).

**The whole architecture is one move: lift the invariant out of the 2023 instance, and make the year/build/observable a parameter that flows through data → code → figures → memory → paper without a source edit.**

---

## 2. Target layout

### 2a. Code (`lib/`, `reduce/`, `analyze/`)

```
lib/                         # INVARIANT — already ~90% portable, harden it
  io/    BuildConfig.h  DataPaths.h  RadView.h   FigPaths.h(NEW)
  waveform/ ChannelConfig.h(legacy 2023 lock — quarantine)
  physics/ RadTiming.h  SelectionCuts.h
  viz/   RADiCALStyle.h
reduce/                      # INVARIANT — gold standard, leave alone
  Reducer.C  reduceRun.C  hpc/{env.sh,submit_reduce.sh}   # already RAD_YEAR-clean
analyze/
  *.C                        # PORTABLE macros (build+year parameterized)
  studies/*.C                # one-off / legacy (quarantined, runnable)
  makeReport.py              # parameterize: --build --year, template "May 2023"
```

**Changes vs today:**
- `DataPaths.h:26` `kRadYear=2023` becomes **`RAD_YEAR` env-driven** (read once, like the existing `RAD_DATA`), so switching campaigns needs no recompile. The year-aware overloads already exist — wire the default to the env.
- Add **`lib/io/FigPaths.h`** mirroring `DataPaths.h`: `radFig("floor_model", build)` → `figures/<year>/narrative/floor_model_<build>.png`. This is the single highest-leverage fix — it makes 2023 and 2024 figures coexist instead of overwriting.
- Modern macros switch `Form("data/2023/configs/%s.json", build)` → **`radConfig(build)`** (already imported in 25 of them, just bypassed).
- Add **`radEnergies(build)`** helper (reads `cfg.runs` keys) to kill the 41 inline `{25,…,150}` literals.
- The ~46 `ChannelConfig.h`/`kRuns`-locked legacy macros stay in `analyze/studies/` (runnable archive); migrate only the most-cited (`timingEnergyBins.C`, `qualityPlots.C`, `compareEnergies.C`) onto `RadView`/`BuildConfig` when a paper needs them.

### 2b. Data (already correct — finish it)

```
data/
  registry.csv               # canonical build×energy×year index (make code read it)
  TEMPLATE.json              # fix stale "datasets/"→"data/" note
  _TEMPLATE/                 # clone-source for a new year
  2023/{raw,reduced,configs/<BUILD>.json,metadata/}   # real
  2024..2026/                # placeholders — fill configs/<BUILD>.json per build
```

**Changes vs today:** the JSON is authoritative; `metadata/*.yaml` are unread mirrors — keep them as the human logbook but stop treating the duplicated channel map as a second source of truth (or generate the yaml from the JSON). Fix the README path-rot (`config/`→`configs/`+`metadata/`, dead `Analysis/`, `datasets/`).

### 2c. Figures (year-namespaced, with a sync rule)

```
figures/<year>/{narrative,capillary,...}/<name>_<build>.png   # via FigPaths.h
output/<year>/{Summary,report_images,report.html}             # report, year-scoped
papers/<topic>/figs/                # synced from figures/<year>/ via a manifest
papers/<topic>/figs.map             # NEW: declarative "paperfig ← source" map
```

**Changes vs today:** every `Print()`/`SaveAs` path gains a `<year>` segment (today **zero** do); fix the `Analysis/Output` vs `output` case bug and the dead `Analysis/capillary_figs/` target (8 macros write to a dir that doesn't exist; outputs were hand-moved to `site/assets/`); replace the hand-copy + prose rename map in `papers/timing/README.md` with `figs.map` consumed by a `make figs` target.

### 2d. Memory (re-tier to mirror `data/<year>` + `papers/<topic>`)

| File | Tier | Single job |
|---|---|---|
| `MEMORY.md` | 0 INDEX | route to the right tier + scope map; never store facts; **sole home of the wording invariant** |
| `detector_invariant.md` | 1 INVARIANT | the RADiCAL detector concept, timeless (design, WLS-is-variable, T/E roles, dual-gain, parent paper, program targets) |
| `methods_invariant.md` | 1 INVARIANT | how we analyze ANY campaign (BestMinus, brightest-K, Gaussian-core σ, OOS folding, floor-physics decomposition, DRS4 quirks, the abstractions, house style, ROOT pitfalls) |
| `dataset_2023.md` | 2 PER-DATASET | the 2023 beam-time spec + measured results (builds table, run ranges, channel-map indices, 42.25 V, 820 mV/jitter/light numbers, manifest) |
| `paper_timing.md` | 3 PER-PAPER | Paper-1 thesis, kill-shot, floor-honesty, systematics, LOCKED HEADLINE METHOD, + P1 figure index |
| `paper_energy_position.md` | 3 PER-PAPER | Paper-2 plan, E-type/4D/TENERGY, localization, + P2 figure index |
| `log_2023.md` | X LOG | dated provenance trail (write-mostly; **never a lookup source**) |
| `dataset_TEMPLATE.md`, `log_TEMPLATE.md` | template | clone-source for a new campaign |

**Changes vs today:** split `experimental_setup.md` (≈70% invariant → Tier 1, ≈30% 2023 run-bookkeeping → `dataset_2023.md`); rename+split the misnamed `radical_apparatus_conclusions.md` (§6 floor physics + §10 primitives → `methods_invariant.md`; everything else → the two `paper_*.md`); dissolve `figure_catalog.md` into per-paper `figs:` sections; rename `project_radical_testbeam.md` → `log_2023.md`. De-duplicate the kCap order, MCP-split story, floor numbers, and HG-clip explanation to ONE home each; delete the stale 6.5/4.5 kCalo center from the log.

---

## 3. The reusable core (what exists, harden it; small gaps to close)

**Already built and genuinely dataset-agnostic — these are the asset:**

- **`reduce/Reducer.C` + `reduceRun.C`** — the gold standard. Zero literals; derives every path from `cfg.year`; stamps build/year into output. HPC side (`reduce/hpc/env.sh`, `submit_reduce.sh`) already `RAD_YEAR`-parameterized. **Reduction of a new year works today.**
- **`lib/io/BuildConfig.h`** — the keystone. One JSON fully describes a build (DRS4 `[drs,grp,ch]` map, per-corner material/role, geometry, run list, HG/LG sidecar). New dataset = new JSON, no recompile.
- **`lib/io/RadView.h`** — role-resolved, format-agnostic accessors; works on canonical and legacy slot files; `beamCenter()` works on every build.
- **`lib/physics/RadTiming.h`** — THE one timing method; build-agnostic via RadView; `is_timing(c)` reads roles from config, so MIXED/TENERGY heterogeneous corners just work.
- **`lib/io/DataPaths.h`** — year-aware overloads already present (`radReduced(build,E,year)`); `$RAD_DATA` override.
- **`lib/physics/SelectionCuts.h`** — all cuts as named constants, config-independent.
- **`analyze/makeReport.py`** — the layered HTML report; its `PlotEntry.png_stem` collision-avoider (lines 255–269) is exactly the right instinct, just needs lifting one level to dataset.

**The small gaps to close so a new dataset "just runs" (all cheap, the abstraction exists and is merely bypassed):**

1. **`RAD_YEAR` on the analysis path.** `DataPaths.h:26` hardcodes `kRadYear=2023` and every macro uses the year-less overload. Wire the default to the env. *(This is the one real engineering blocker; everything else is data-entry or mechanical sed.)*
2. **`radConfig(build)` instead of literal `data/2023/configs/`.** Mechanical sed across 25 macros that already import `DataPaths.h`.
3. **`radEnergies(build)`** reading `cfg.runs` keys → kills 41 inline energy literals.
4. **`FigPaths.h`** so figures are year-namespaced (today they collide across years).
5. **`makeReport.py --build --year`** + template the "CERN SPS May 2023" strings; fix the `Analysis/Output` case bug.
6. **`registry.csv` actually read by code** (it is the build×energy source of truth but no `.C` reads it).

---

## 4. Two playbooks (the anti-reinvent-the-wheel mechanisms)

### Playbook A — NEW DATASET (e.g. onboard 2024)

1. **Clone the skeleton.** `cp -r data/_TEMPLATE data/2024` (done). Fill `data/2024/metadata/{dataset.yaml,runs.csv}` from the beam logbook (replace the `# TODO`s).
2. **Write the build JSON(s).** For each build, copy `data/2023/configs/<BUILD>.json` → `data/2024/configs/<BUILD>.json`; set `year:2024`; **verify every `[drs,grp,ch]` triple against the 2024 wiring** and `hg_sat_mV` against the 2024 DC offset. This JSON is the only file that matters — yaml is decorative.
3. **(Optional) regenerate the HG/LG sidecar** via `reduce/calibHGLG.C` from clean low-E 2024 data.
4. **Reduce.** Build the manifest (`reduce/hpc/build_manifest.py`), then `RAD_YEAR=2024 RAD_CONFIG=<BUILD> RAD_MANIFEST=... reduce/hpc/submit_reduce.sh`. Lands in `data/2024/reduced/`. **Works today** — reduce is year-clean.
5. **Validate (out-of-sample, loudly).** Run the standard validation set with `RAD_YEAR=2024`: `analyze/sigmaT.C` (timing sanity), an energy-linearity check, `analyze/studies/lgcheck` per build (MIXED-LG health), and confirm `RadView::beamCenter()` lands on the fiducial. **Validate any correction train-even/test-odd before claiming a gain** (hard-won rule from `GAP_CLOSING_PLAN.md`).
6. **Standard figure set.** Run the portable narrative macros with `RAD_YEAR=2024`; outputs land under `figures/2024/...` via `FigPaths.h` (no collision with 2023). **Open every regenerated PNG before declaring it done.**
7. **Report.** `python analyze/makeReport.py --year 2024 --build <BUILD>` → `output/2024/report.html`; `analyze/deployReport.sh` → `site/`.
8. **Dataset memory entry.** `cp dataset_TEMPLATE.md dataset_2024.md`; fill builds table, run ranges, channel-map indices, measured numbers, VERIFY list. Inherit **all of Tier 1 unchanged** (write nothing already in `detector_invariant.md` / `methods_invariant.md`). Add one routing line to `MEMORY.md`. Start `log_2024.md`.

### Playbook B — NEW PAPER (e.g. a 2025 timing paper, or Paper 3)

1. **Pick configs + observable.** State the build set and the single observable (timing / energy / position). Confirm the two open author gates if relevant (LuAG NW capillary T-vs-E; exclude run 2247 from the LUAG timing set) — carry forward, don't re-derive.
2. **Reuse macros, don't write new ones.** The portable layer covers it: `RadTiming.h` for timing, `radReduced(build,E)` for data, the existing narrative macros for figures. Only write a macro if the observable genuinely has no precedent.
3. **Generate figures** into `figures/<year>/...` via `FigPaths.h`. Lead with the clean head-to-head where cross-build means are confounded (the `mixed_h2h.png` in-event proof — do **not** try to "fix" the LY overlap; it's understood and confounded by energy).
4. **Outline from template.** Copy the structure of `papers/timing/` (or `papers/energy_position/`). The 2-paper split and outlines in `docs/PAPER_OUTLINES.md` / `docs/FOLLOWON_PLAN.md` are the rationale of record — reuse them.
5. **Tex skeleton + figure sync.** Create `papers/<topic>/figs.map` (`thesis.png ← figures/<year>/narrative/light_yield_thesis.png`, …); `make figs` syncs sources into `papers/<topic>/figs/`. No hand-copy, no prose rename map.
6. **Systematics + catalog.** Reuse the systematics budget approach from `docs/` (already DONE for P1). Record each figure's LIVE/SUPERSEDED status **inline next to the claim it supports** in the per-paper memory file, with its generating macro — not in a global catalog.

---

## 5. Migration path (ordered; quick wins first)

**Quick wins (hours, mechanical, low-risk):**

1. **Wire `RAD_YEAR` into `DataPaths.h`** default (replace the `kRadYear=2023` constant with an env read). *Unblocks the entire analysis path for new years.* — the single most important change.
2. **Sed `Form("data/2023/configs/%s.json",build)` → `radConfig(build)`** across the 25 modern macros that already import `DataPaths.h`.
3. **Fix the `makeReport.py` case bug** (`Analysis/Output` → `output`) and the dead `Analysis/capillary_figs/` target in the 8 capillary macros. Unblocks case-sensitive HPC.
4. **Fix README/template path-rot** (`config/`→`configs/`+`metadata/`, dead `Analysis/`/`datasets/`/`organize_data.sh` references, fix the `datasets/` note in `data/TEMPLATE.json`).
5. **Rename memory files** in place: `project_radical_testbeam.md` → `log_2023.md` (+ "write-mostly, not a lookup source" header). Remove the copy-pasted companion-header block from 4 files; point all to `MEMORY.md`.

**Medium (a day each, higher payoff):**

6. **Add `lib/io/FigPaths.h`** and route narrative-macro `Print()` paths through it (year-namespaced). Ends all cross-year figure collisions; fixes the `opt_scan.png` fixed-name duplicate.
7. **Add `radEnergies(build)`** + replace the 41 inline energy literals.
8. **Split `experimental_setup.md`** → `detector_invariant.md` + `methods_invariant.md` (invariant) and `dataset_2023.md` (2023 bookkeeping). Create `dataset_TEMPLATE.md` + `log_TEMPLATE.md`.
9. **Rename + split `radical_apparatus_conclusions.md`** → §6/§10 to `methods_invariant.md`; the rest to `paper_timing.md` / `paper_energy_position.md`. Dissolve `figure_catalog.md` into per-paper `figs:` sections.
10. **`papers/<topic>/figs.map` + `make figs`** sync rule; delete the prose rename map (and its stale `narrative_hglg` entry).
11. **`makeReport.py --build --year`** argparse + template the "May 2023 / H2" strings.

**Larger (only when a paper demands it):**

12. **Migrate the heavy legacy macros off `ChannelConfig.h`** (`timingEnergyBins.C`, `qualityPlots.C`, `compareEnergies.C`) onto `RadView`/`BuildConfig`+`radReduced`. Highest effort; do it lazily, per-paper, not all at once.
13. **Make `registry.csv` the read source** for builds/energies (so a macro can enumerate a campaign).

The order is deliberate: items 1–5 make 2024 *analyzable*; 6–11 make it *publishable side-by-side with 2023*; 12–13 retire the last 2023 lock-in.

---

## 6. What NOT to do

- **Don't duplicate the invariant per dataset.** `dataset_2024.md` must inherit `detector_invariant.md` / `methods_invariant.md` verbatim — never re-state the BestMinus method or the floor physics per campaign. The whole point is to re-validate predictions, not re-derive them. If you find yourself copying a method paragraph into a dataset file, it belongs in Tier 1.
- **Don't over-template.** Templates only for what genuinely recurs (`data/_TEMPLATE/`, `dataset_TEMPLATE.md`, `log_TEMPLATE.md`, `figs.map`). Do **not** template per-paper *narrative* — the thesis of each paper is bespoke. A template that needs heavy per-use surgery is worse than a copied example.
- **Don't let figures lose their generating macro.** Every figure must be regenerable from a named macro to a cataloged path. The current `Analysis/capillary_figs/`→`site/assets/` hand-move and the byte-identical paper hand-copies are exactly the anti-pattern. `FigPaths.h` + `figs.map` exist to make the macro→path→paper chain mechanical and unbreakable.
- **Don't let the channel map drift between two homes.** The JSON is authoritative; the `metadata/channel_map.yaml` is decorative. Either generate the yaml from the JSON or mark it read-only — never hand-edit both.
- **Don't reintroduce path fallbacks.** The reorg deliberately deleted legacy `Data/`, `reduced/`, `output/` fallbacks so `DataPaths.h` is the single resolver. Keep it that way; add `RAD_YEAR`, not a search list.
- **Don't grow the log as a knowledge store.** `log_2023.md` is a dated trail (it records the same MIXED-LG question resolved three times). Distillations go UP to Tier 1/2; the log is provenance only. Cap it (~1,000 lines) and archive resolved threads.
- **Don't "fix" understood negatives.** The cross-build σ_t-vs-LY binned-mean overlap is confounded by energy and is *understood*; the in-event head-to-head is the clean proof. Don't burn effort re-litigating it each campaign — carry the conclusion forward.
- **Don't recompile to switch campaigns.** If onboarding 2024 ever requires editing `kRadYear` or `ChannelConfig.h` and rebuilding, the parameterization has failed. Year/build/observable flow as parameters, never as source edits.
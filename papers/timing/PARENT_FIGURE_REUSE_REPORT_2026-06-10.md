# Parent-Figure Reuse Report — timing manuscript (2026-06-10)

Scope: `papers/timing/NIM_A_Figures/` (88 files) mapped against the parent paper
NIM A 1068 (2024) 169737 (Figs. 1–28, read in full from the published PDF), cross-checked
against the Agent-5 figure audit (`PRIOR_RECORD_AUDIT_2026-06-09.md`) and Sec. 6 of
`PRIOR_RECORD_ACTION_PLAN_2026-06-09.md`, and against the current manuscript
`radical_timing.tex` (Figs. 1–8, A.9–A.10, B.11).

Naming note: the repo (`papers/timing/README.md`) labels this manuscript **Paper 1 — Timing**;
the task brief called it "Paper 2". This report follows the repo convention ("the timing paper").

Claims-law constraints applied throughout: no parent σ_t(E)/result plots as primary content
(published result lives only as the Table-1 control row + the separately planned grey reference
curve); canonical module numbers only (14×14×135 mm³, 29 LYSO:Ce + 28 W, X₀ = 5.4 mm,
R_M = 13.7 mm, ≈25 X₀); LYSO:Ce + W common to all builds, the WLS capillary is the variable;
every reused figure captioned "adapted from Ref. [1]"; GEANT4 figures only as the *published*
expectation of Ref. [1], never as new simulation work.

---

## 1. Inventory map (folder file → parent figure)

Apparatus / concept (all verified by direct viewing):

| File | Parent figure | Content |
|---|---|---|
| `lyso_w_cell_15.pdf` | **Fig. 1** | Module schematic: W (2.5 mm)/LYSO (1.5 mm) shashlik, 14 mm × 135 mm, quartz capillaries, monitoring fiber, SiPM, readout cards. Vector. |
| `lyso_plate_12.pdf` | **Fig. 2** | Tile cross-section: 14.0 mm square, 4 × ⌀1.3 mm corner through-holes at ±3.50 mm, 1 × ⌀0.9 mm central hole, 1.50 mm thickness. Vector. |
| `Caps.pdf` | **Fig. 3** | Photo: WLS capillaries under 365 nm UV; E-type (full-length WLS) upper/lower, T-type (shower-max filament) middle; DW/UP ends labelled. |
| `Caps2.pdf` | **Fig. 4** | LED scan: signal amplitude vs position for E-type (flat, blue) vs T-type (localised, orange). |
| `lyso_w_cell_17.pdf` | **Fig. 5** | Module schematic, E-type capillaries (full-length red WLS). Vector. |
| `lyso_w_cell_16.pdf` | **Fig. 6** | Module schematic, T-type capillaries (short red WLS segments at shower max). Vector. |
| `sim_v0.pdf` / `sim_v0.png` | **Fig. 7 (published)** | GEANT4 Edep/E0 per LYSO layer, 20/50/125 GeV, fixed WLS-filament band marked below axis ("Layer Index"). |
| `EDep.pdf` | Fig. 7 **pre-pub variant** | Same content, earlier render ("Tile Index", different palette). Use `sim_v0.pdf` instead. |
| `TimingSasha.pdf` | **Fig. 8** | GEANT4 simulated σ_t vs detected LY (log–log), single-end readout study. |
| `moliere.pdf` | **Fig. 9** | GEANT4 transverse energy map at shower max, 50 GeV, 14×14 mm² boundary. |
| `Fig10.png` | **Fig. 10(a)** | Module under UV photo. Low-res (543×356) — screen grab, not print quality. |
| `RADICAL_pic.pdf` | **Fig. 10(b)** | Assembly photos: tile stacking (top), POM frame with FE cards (bottom). |
| `3DPrint.pdf` | **Fig. 10(c)** | 3D-printed base plate with MCP enclosure + module, workshop photo. |
| `TandW.jpeg` | **Fig. 11(a)** | Installation on the CMS/HF lift table, H2 beam line (photo). |
| `Setup.pdf` | **Fig. 11(b)** | Overhead closeup: MCP, trigger counter, module enclosure, Pb-glass under cloth (unlabelled photo). |
| `SetupSchematic.pdf` | **Fig. 11(c) (published)** | Beam's-eye beamline photo with A1/A2 (triggers), B (MCP), C (module), D (Pb-glass), BEAM labels. |
| `Picturebeamview1.png` | Fig. 11(c) **pre-pub variant** | Same photo, yellow draft labels. Use `SetupSchematic.pdf` instead. |
| `profiles.png` | **Fig. 12 (published)** | Average normalised waveforms: Energy Channels (LG) / Timing Channels (HG, flat-top = deliberate saturation), downstream red / upstream blue. |
| `Figure12.png`, `figure12new.png` | Fig. 12 **pre-pub variants** | Per-channel waveform overlays (CH1–CH14, DRS0/DRS1) with visible ringing/reflections. Talk-grade. |
| `ninecap.png`, `ninecap2.png` | **Fig. 28** (two panels) | Enhanced-design capillary layouts (5 T-type + 4 E-type, future 18×18 mm² module). |

Results-class (published or preliminary — see Sec. 3):

| File | Parent figure | Note |
|---|---|---|
| `placeholder.001.png` | **Fig. 16 (published)** | E_meas (DW+UP) distributions, six energies, σ annotations. |
| `placeholder.png` | Fig. 16 variant | Same without σ annotations. |
| `placeholder2.png` | **Fig. 20** | (DW−UP)/2 distributions, six energies. |
| `ENE_1-2.png` | **Fig. 17 Left (published)** | Linearity, **29.4 mV/GeV**. |
| `ENE_1.png` | Fig. 17 Left **preliminary** | **29.8 mV/GeV — obsolete number.** |
| `ENE_2-2.png` | **Fig. 17 Right (published)** | σ_E/E = 9.31% ⊕ 52.04%/√E ⊕ 31.62%/E. |
| `ENE_2.png` | Fig. 17 Right **preliminary** | 9.50/50.76/31.62 — **obsolete numbers.** |
| `TIM_1-2.png` | **Fig. 27 Left (published)** | Six-method σ_t(E) survey (Downstream…BestMinus). |
| `TIM_1.png`, `TIM_1-3.png` | Fig. 27 Left variants | Draft palettes/legend orders. |
| `TIM_2-2.png` | **Fig. 27 Right (published)** | BestMinus σ_t(E) with 17.52 ⊕ 255.58/√E fit, error bars. |
| `TIM_2.png` | Fig. 27 Right variant | No error bars, draft legend. |

Per-run analysis class (group identification; ~50 files `DSB1_<E>.root_*` for E = 25…150 GeV,
`-2`/`-3` suffixes are alternate renders):

| File class | Parent figures | Content |
|---|---|---|
| `DSB1_<E>.root_EVE_1*` | **Fig. 13** class | Per-SiPM (x,y) beam maps from the wire chamber. |
| `DSB1_<E>.root_ENE_1*` | **Fig. 14** class | Per-SiPM LG amplitude histograms (DWNW…UPSE) with Gaussian fits. |
| `DSB1_<E>.root_ENE_3*` | **Fig. 15** class | Summed LG amplitudes DW / UP / DW+UP. |
| `DSB1_<E>.root_TIM_3*` | **Fig. 18** class | Per-SiPM t_SiPM − t_MCP distributions. |
| `DSB1_<E>.root_TIM_4*` | **Fig. 19** class | DW, UP, (DW−UP)/2, (DW+UP)/2 averages (satellite visible in the sum). |
| `DSB1_<E>.root_TIM_6*` | **Figs. 21–26** top/middle | Nine-bin timing distributions per measured-energy bin. |
| `DSB1_<E>.root_TIM_7*` | **Figs. 21–26** lower | σ_t vs measured amplitude, (DW−UP)/2 and (DW+UP)/2. |

Every file in the folder is accounted for; no unidentified files.

---

## 2. INSERT recommendations (ranked)

### #1 — Apparatus + corner-map composite → new Fig. 1, Sec. 2.1  (MUST; highest value)
**Sources:** `lyso_w_cell_15.pdf` (parent Fig. 1) + `lyso_w_cell_16.pdf` (parent Fig. 6) +
`lyso_plate_12.pdf` (parent Fig. 2, as the geometric basis for a redrawn beam's-eye corner map).
**Target:** Sec. 2.1 "Module and builds", placed at the top of Sec. 2 before the current
`fig:clip` (current Figs. 1–8 renumber by one). **Format: `figure*`** (three panels across two
columns; (a) wide, (b) and (c) narrower).
**Why:** the manuscript currently has *no* apparatus figure (audit gap 1a), and the NE/NW/SE/SW
corner names that carry the entire Sec. 6 MIXED argument are never drawn (gap 1b). One composite
closes both and gives the Fig. 7 caption an internal pointer instead of a prose digression.
**Compositing/redraw:**
- (a) overlay **DW / UP** end labels (beam enters from the upper left in the parent drawing —
  state the beam direction with an arrow); consider removing or re-captioning the
  "Monitoring fiber" arrow, since the central hole was uninstrumented in these tests
  (parent Sec. 2).
- (b) drop in as-is (vector); optionally crop the legend text.
- (c) **redraw** parent Fig. 2 as a simple beam's-eye square: four corner circles at
  ±3.50 mm offsets, central hole open, corners labelled NE/NW/SE/SW, and (either here or in a
  small twin) the MIXED colouring — DSB1 on NE–SW, LuAG:Ce on NW–SE. Simple TikZ/matplotlib job.
**DRAFT caption (claims-law compliant):**
> **Fig. 1.** The RADiCAL module and its readout geometry (adapted from Ref. [1]).
> (a) Schematic of the module: 29 LYSO:Ce scintillator plates (1.5 mm) interleaved with
> 28 tungsten absorber plates (2.5 mm), 14×14×135 mm³ overall (X₀ = 5.4 mm, R_M = 13.7 mm,
> ≈25 X₀ deep). Four quartz capillaries traverse the stack and are read out by SiPMs at both
> ends; the four downstream ends are labelled DW, the four upstream ends UP.
> (b) T-type capillary: a short WLS filament is positioned in the shower-maximum region, with
> fused quartz waveguides transporting the wave-shifted light to the two ends.
> (c) Beam's-eye view of the 14×14 mm² face: the four capillaries occupy corner positions
> offset ±3.50 mm from the module axes (central hole uninstrumented). In the MIXED build,
> DSB1 occupies the NE–SW diagonal and LuAG:Ce the NW–SE diagonal (Sec. 6).
**Effort:** light edit for (a)/(b), small redraw for (c) — about half a day total.
Also execute the action-plan caption pin for Fig. 7 (`mixed_h2h_corrected`): point it at panel (c).

### #2 — GEANT4 longitudinal shower profile → Sec. 7 (SHOULD; optional per action plan)
**Source:** `sim_v0.pdf` (the *published* parent Fig. 7 — NOT `EDep.pdf`, which is the pre-pub
"Tile Index" variant).
**Target:** Sec. 7 "Longitudinal timing asymmetry", single column, beside/before
`fig:depthdial`; also referenced once from the Discussion floor paragraph.
**Why:** makes the depth-drift mechanism visible — shower max migrating past the fixed filament
band is exactly the −33.6 ps/e-fold story — and converts "no simulations" (referee risk M2) into
"our measurement matches the programme's published expectation." Claims-law safe because it is
reproduced *as the published expectation of Ref. [1]*, with the paper's own G4 campaign still
honestly deferred.
**DRAFT caption:**
> **Fig. X.** Published GEANT4 expectation for the longitudinal shower profile in this module
> geometry: fraction of incident energy deposited per LYSO:Ce layer for 20, 50 and 125 GeV
> electrons, with the fixed longitudinal position of the T-type WLS filament indicated (band).
> Shower maximum migrates downstream with energy past the fixed filament position — the
> mechanism probed by the mean asymmetry drift of Fig. 8. Adapted from Ref. [1] (its Fig. 7).
**Effort:** drop-in (vector; no edits needed). The cite-only fallback of the action plan
(item 7) remains acceptable if figure count is a concern; if not inserted, add the citation
sentence to the Fig. 8 caption instead.

### #3 — Beamline arrangement → Sec. 2.2 (OPTIONAL; default remains cite-only)
**Source:** `SetupSchematic.pdf` (published parent Fig. 11(c), labelled A1/A2/B/C/D/BEAM).
**Target:** Sec. 2.2 "Readout, beam, and reference", single column.
**Why:** answers "where was the MCP / what triggered / what vetoed leakage" visually; directly
supports the Pb-glass veto and MCP-reference text. The audit (1e) judged a citation sufficient
("the H2 beamline arrangement of Ref. [1], its Fig. 11"), and that remains the default; insert
only if coauthors want the visual or a referee asks for the setup.
**DRAFT caption (if inserted):**
> **Fig. X.** Experimental arrangement in the CERN H2 beam line, viewed along the beam:
> scintillation trigger counters A1 (1×1 cm²) and A2 (2×2 cm²), the MCP-PMT timing reference (B),
> the RADiCAL module in its light-tight enclosure (C), and the Pb-glass backing array (D).
> Adapted from Ref. [1] (its Fig. 11).
**Effort:** drop-in.

### Caption-level reuse (no new figures; ride on existing edits)
- **Deliberate-saturation provenance** (action-plan item 4): append to the current `fig:clip`
  caption — "The high-gain chain is deliberately driven into saturation to maximise the slope of
  the recorded leading edge (Ref. [1], Sec. 4); Sec. 3 recovers the clipped peak event by event."
  The folder's `profiles.png` (published Fig. 12, flat-top timing waveforms) is the visual
  antecedent; cite it ("its Fig. 12") rather than reproducing — the manuscript's own Fig. 3
  (`pulse.png`) already shows a real clipped waveform better.
- **T-type concept anchoring:** if composite panel (b) is dropped for space, fall back to citing
  parent Figs. 3–6 in Sec. 2.1 text. `Caps.pdf` (Fig. 3 UV photo) is the photogenic alternative
  for panel (b) if a photo is preferred over the schematic; both carry the same provenance line.

---

## 3. DO-NOT-REUSE list

| File(s) | Reason |
|---|---|
| `TIM_1.png`, `TIM_1-2.png`, `TIM_1-3.png` (Fig. 27 Left + variants) | Published timing *results* — forbidden as primary content (claims law). Its method names (BestMinus/BestPlus) also conflict with the srCFD/LED/cfd05 taxonomy; if ever referenced, the one-time mapping sentence "(DW−UP)/2, called BestMinus in Ref. [1]" is mandatory. |
| `TIM_2.png`, `TIM_2-2.png` (Fig. 27 Right + variant) | The published σ_t(E) curve. Lives in the manuscript ONLY as the Table-1 control row plus the planned grey 256/√E ⊕ 17.5 ps reference curve on the money plot (Fig. 5 restyle, action-plan item 2). Re-inserting the parent plot would present superseded primary content alongside the unified-chain results. |
| `ENE_1.png` (29.8 mV/GeV), `ENE_2.png` (9.50/50.76 fit) | **Preliminary fit numbers that differ from the published record** (29.4 mV/GeV; 9.31 ⊕ 52.04/√E ⊕ 31.62/E). Never propagate. |
| `ENE_1-2.png`, `ENE_2-2.png` (published Fig. 17) | Correct numbers, wrong paper: energy response is the companion energy/position paper's territory (audit "not recommended for Paper 1"). |
| `placeholder.png`, `placeholder.001.png`, `placeholder2.png` (Figs. 16, 20) | Parent per-energy E_meas and (DW−UP)/2 distributions — superseded by the unified-chain `dist.png` (current Fig. 4); re-inserting would duplicate with legacy-estimator content. |
| `DSB1_<E>.root_TIM/ENE/EVE_*` class (~50 files; parent Figs. 13–26 source panels) | Per-run legacy-analysis output (fixed-threshold timing, measured-energy bins 6–8). Results-class, superseded by the unified chain; useful only as provenance archive for coauthor questions. |
| `Figure12.png`, `figure12new.png` | Pre-publication waveform variants (per-channel overlays with DRS ringing/reflections) — talk-grade; the published `profiles.png` exists, and the manuscript's `pulse.png`/`hglg.png` already cover waveform behaviour. |
| `Picturebeamview1.png` | Pre-pub variant of Fig. 11(c) (draft yellow labels). If the beamline photo goes in at all, use the published `SetupSchematic.pdf`. |
| `EDep.pdf` | Pre-pub variant of Fig. 7 ("Tile Index" axis, different palette). Use `sim_v0.pdf`. |
| `ninecap.png`, `ninecap2.png` (Fig. 28) | Future enhanced-design layout, **18×18 mm² module** — non-canonical dimensions; inserting risks implying new-design claims and dimension contamination. Cite Ref. [1] Sec. 7 in the Discussion only (already done for item 3). |
| `TimingSasha.pdf` (Fig. 8, σ_t vs LY simulation) | Audit 1d is explicit: **cite, do not reproduce** (single-end readout, different conditions). The one-sentence "anticipated by the programme's GEANT4 studies [1]" edit (E4) is the right vehicle. |
| `moliere.pdf` (Fig. 9) | Transverse-compactness simulation — supports the R_M narrative but is not load-bearing for any timing section; cite only if needed. |
| `TandW.jpeg`, `Setup.pdf`, `3DPrint.pdf`, `RADICAL_pic.pdf`, `Fig10.png`, `Caps2.pdf` | Assembly/installation photos and the LED scan: fine provenance, low information density for this paper; the timing paper needs geometry, not photographs. `Fig10.png` is additionally print-unusable at 543×356. |

---

## 4. Open questions / VERIFY

1. **Corner-label orientation for composite panel (c).** Parent Fig. 2 carries no compass
   labels; NE/NW/SE/SW must be drawn from the authoritative 2023 channel map
   (memory `dataset_2023.md`) with the viewing convention stated explicitly (beam's-eye =
   looking downstream, as in parent Fig. 13's wire-chamber maps). The DSB1 = NE–SW /
   LuAG:Ce = NW–SE assignment is from manuscript Sec. 6 (pulse-shape-confirmed); cross-check
   the drawn panel against that gate record before labelling. This is referee risk #7
   (corner-map provenance) — the panel must not introduce a left/right flip.
2. **T-type vs E-type schematic identity.** `lyso_w_cell_16.pdf` = T-type (short red WLS
   segments, parent Fig. 6) and `lyso_w_cell_17.pdf` = E-type (full-length red, Fig. 5) —
   identified visually with high confidence; one-line confirmation against the parent's
   figure source files would close it.
3. **Monitoring-fiber arrow in panel (a).** Parent Fig. 1 labels a central "Monitoring fiber",
   but the parent text states the central hole was not used in these tests. Decide: remove the
   arrow in the adapted version, or keep it with the "uninstrumented" note in the caption
   (current draft does the latter via panel (c)).
4. **Reuse permission formality.** Parent is CC BY-NC-ND 4.0 and collaboration-own (same author
   pool); "adapted from Ref. [1]" suffices scientifically. At submission, check Elsevier's
   own-work reuse boilerplate (NIM A → NIM A) in case the production team wants a credit line
   format ("Reprinted/adapted from [full citation], © 2024 The Authors, CC BY-NC-ND").
5. **Vector fidelity.** `lyso_w_cell_15/16/17.pdf` and `lyso_plate_12.pdf` are true vector PDFs
   (good); if labels are added, do it in the PDF/TikZ layer, not by rasterising.
6. **Paper numbering discrepancy.** Task brief says "Paper 2 (timing)"; repo and memory say
   Paper 1 = timing, Paper 2 = energy/position. No action in the manuscript; flagged so the
   circulation note doesn't inherit the wrong label.

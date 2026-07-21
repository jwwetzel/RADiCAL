# Campaign Snapshot — 2026-06-09

Frozen state of the publication campaign at the moment GATE 3 (depth dial) provisionally passed.
Purpose: reproducibility anchor for internal review; everything below was verified by direct command,
not recalled from memory.

## Repository state
- Branch: `main` · HEAD: `e7ac826` "Add methodDist.C — per-point (DW-UP)/2 distribution + Gaussian-fit diagnostic"
- **All campaign work is UNCOMMITTED** as of this snapshot.
- Modified (tracked): `analyze/studies/{methodCompare,methodDist,waveformProfiles}.C`,
  `lib/physics/{RadTiming,SelectionCuts}.h` (the estimator fix + continuous fiducial + veto constant),
  `lib/viz/PlotUtils.h` (DrawSuperTitle/DrawEnergyLegend/GridWithTitle), 6 regenerated narrative PNGs.
- Untracked (new): 9 analysis macros under `analyze/studies/` (adaptiveTiming, build_chain_html.py,
  mixedSeparate, monotonicityEvidence, monotonicityFix, outlierPeek, reductionQA, sigmaProbe,
  timingAllMethods, timingRegression), `chain_of_evidence.html`, 16 new narrative figures,
  `papers/` additions: EXPERT_PANEL_BLUEPRINT.md, two panel JSONs, six memory_*.md,
  `papers/scripts/depth_dial/` (AUDIT.md, depthDial.C, depth_dial_result.log),
  `papers/figures/depth_dial/` (depth_dial.png, depth_dial_diag.png).
- Suggested commit message (NOT yet committed):
  `Add RADiCAL paper memory infrastructure and depth-dial gate analysis`

## Environment
- ROOT 6.40.00 (Homebrew), macOS 15.7.7, arm64 (Apple Silicon).
- Env via `source setup.sh` (sets ROOT_INCLUDE_PATH to lib/{waveform,io,physics,viz}, RAD_DATA=repo root).
- RAD_YEAR unset → defaults to 2023 (DataPaths.h `radYear()`).

## Depth-dial run (GATE 3)
- Exact command:
  ```
  source setup.sh
  root -l -b -q -e '.L papers/scripts/depth_dial/depthDial.C+' -e 'depthDial()'
  ```
  (output teed to `papers/scripts/depth_dial/depth_dial_result.log`)
- Pre-registered method: `papers/scripts/depth_dial/AUDIT.md` (written BEFORE plotting).
- Inputs (reduced 'rad' trees, reduced on Argon Jun 7 2026):
  | file | bytes | mtime |
  |---|---|---|
  | data/2023/reduced/DSB1/25GeV.root | 501,953,492 | Jun 7 00:54 |
  | data/2023/reduced/DSB1/50GeV.root | 808,201,197 | Jun 7 00:55 |
  | data/2023/reduced/DSB1/75GeV.root | 596,740,047 | Jun 7 00:55 |
  | data/2023/reduced/DSB1/100GeV.root | 771,506,031 | Jun 7 00:51 |
  | data/2023/reduced/DSB1/125GeV.root | 723,530,216 | Jun 7 00:52 |
  | data/2023/reduced/DSB1/150GeV.root | 1,380,189,183 | Jun 7 00:53 |
  | + LUAG/TENERGY 50–150 GeV (cross-checks) | — | same reduction batch |
- Outputs: `papers/figures/depth_dial/depth_dial.png` (publication),
  `papers/figures/depth_dial/depth_dial_diag.png` (diagnostics),
  log `papers/scripts/depth_dial/depth_dial_result.log`.

## Scientific interpretation (one paragraph)
The mean dual-end timing asymmetry ⟨(t_DW − t_UP)/2⟩ of the DSB1 module, measured with the
amplitude-independent srCFD estimator on the full fiducial sample at fixed r = 2.5 mm and fixed
channel composition, drifts monotonically from −195.7 ps to −263.0 ps as the beam energy rises from
25 to 150 GeV — negative-going exactly as expected when shower maximum migrates downstream (closer to
the DW ends) by ≈1 X₀ per e-fold of energy. The fitted slope, −33.6 ± 2.9 ps per e-fold, agrees in
sign and order of magnitude with the geometric expectation −X₀/v_g ≈ −26 ps per e-fold; the ~28%
excess is consistent with modal dispersion lowering the effective light-propagation speed in the
quartz capillary, to be pinned by simulation. The dial is read by constant-fraction estimators only
(cfd05 mid-range slope ≈ −23 ps/e-fold concurs) while fixed-threshold LED crossings suppress it
(−5.5 ps/e-fold), establishing that the depth observable must be CFD-based. This is the first
measured, population-level demonstration that the RADiCAL dual-end timing carries longitudinal
shower-depth information — the physics hinge between the timing paper (the floor is depth) and the
energy/position paper (the fifth coordinate).

## Claim language in force at snapshot time
(authoritative copy: `papers/memory_claims_and_forbidden_language.md`)
- ALLOWED: "the measured mean of (t_DW−t_UP)/2 tracks the ln(E) shower-max migration, slope
  −33.6 ± 2.9 ps per e-fold, consistent with −X₀/v_g"; "toward 5D shower-maximum calorimetry"
  (population-level scoping mandatory; CFD-based caveat sentence required).
- FORBIDDEN (unchanged): per-event depth resolution / calibrated z; "5D" as achieved; t_diff as
  clock-referenced time-of-arrival; mm numbers as calibrated position resolution; any revision wording
  against the published 17.5 ps / 27 ps; E_SM-untagged energy-resolution numbers.

## UPDATE (later 2026-06-09): committed + GATE 2 closed
- Campaign committed: `45f15dc` "Add RADiCAL paper memory infrastructure and depth-dial gate analysis"
  (60 files). GATE 2 committed: `c8280f9` "GATE 2 passed: position arithmetic reconciled" (8 files).
  Working tree clean after both. Branch `main` is ahead of `origin/main` (not pushed).
- **GATE 2 outcome: PASSED (RECONCILED).** The 1.5 mm is a legitimate unbinned event-level residual
  RMS that survives a train/test split (Δ = 3.5 µm); the "3.6 mm comparator" was the t₀-inflated
  end-time-SUM upper bound and never applied to the position difference (memory mischaracterization,
  now corrected everywhere). Held-out residuals: 1.54/1.45 mm (x/y, ±6 mm window, 150 GeV),
  0.91/0.88 mm in the r<2.5 beam core, ~energy-independent. New finding: closure slope 0.698 —
  linear light-division estimator saturates beyond |x|≈3 mm. Joint-upper-bound language only; no
  intrinsic resolution; unfolding vs 3.6 demonstrated imaginary and refused.
  Products: `papers/scripts/position_reconciliation/`, `papers/figures/position_reconciliation/`.

## UPDATE 2 (later 2026-06-09): GATE 6 CONDITIONAL — the 0.99 kill-shot number RETIRED
- Pre-state verified: working tree clean at `c8280f9`; commits 45f15dc + c8280f9 confirmed.
- **GATE 6 outcome: CONDITIONAL.** The recorded MIXED ratio 0.99 (`mixedHeadToHead.C`) used
  brightness-threshold corner labels (amp > 0.65·max) that misassign SE-D (LuAG by GATE-1
  pulse-shape) into the DSB1 group at mid/high E — the 0.99 and χ²/ndf=0.4 are a label-mixing
  artifact and are retired from all claims. Replacement (common-event, GATE-1 labels, Down-layer
  same-material diagonal pairs, srCFD): per-capillary width parity within ~10–20%; 5-E mean ratio
  1.038, bootstrap 95% [1.033,1.043], honest scatter-based 1.04 ± 0.05; methods srCFD/cfd05/led =
  1.04/1.16/0.80; jackknife spread 0.003 (76 variants, no single-run dependence); swapped map
  inverts exactly; scrambled map collapses the amplitude control. Kill-shot claim rewritten
  (see claims memory): within-~20% same-shower parity vs ×2.3 cross-build → kinetics disfavored as
  dominant, with the position-coupling dilution caveat. `mixed_h2h.png` must be regenerated before
  any draft circulates. Products: `papers/scripts/mixed_killshot_bootstrap/`,
  `papers/figures/mixed_killshot_bootstrap/`.

## UPDATE 3 (later 2026-06-09): stale-claim purge + corrected kill-shot figure + post-fix (a,b) table
- Working tree was clean at `4c49888` before this pass.
- Stale 0.99/χ² language removed from all ACTIVE Paper-2 files (`papers/timing/THREAD.md`,
  `papers/memory_paper2_timing.md`); the panel blueprint + JSONs marked HISTORICAL with a
  supersession banner; `analyze/studies/mixedHeadToHead.C` carries a DEPRECATED header; the old
  figure renamed `figures/2023/narrative/mixed_h2h_DEPRECATED.png` + warning README.
- Official corrected figure: `papers/figures/mixed_killshot_bootstrap/mixed_h2h_corrected.{png,pdf}`
  (3 panels: widths / ratio with 1.04±0.05 scatter band / method dependence; approved caption in
  `makeMixedKillshotFigure_CAPTION.txt`).
- **POST-FIX authoritative four-build table:** `papers/tables/timing_fit_summary_2026-06-09.md` —
  DSB1 srCFD a=203±6, b=18.8±0.8, σ(150)=25.7±0.6 (HOLDS vs old 201/19.5; headline becomes 25.7);
  LuAG LED a=440±18, b=24.6±3.3 (a HOLDS; floor up from 19.8); TENERGY LED a=198±21, b=26.0±1.8
  (√(4/3) penalty); MIXED module-wide 253±29/34.3±2.4 (ill-posed ref). a-ratio 2.16 ("more than a
  factor of two"). **Floor narrative SOFTENED:** floors span 18.8–26.0 ps; "shared ~20 ps" retired
  without qualifiers; DSB1's 18.8±0.8 still CONFIRMS the published 17.5.

## Open flags at snapshot time
- 125→150 GeV local re-steepening of the dial (−56 ps/e-fold locally) — run-period check pending
  (addressed in DEPTH_DIAL_REVIEW.md diagnostics).
- GATE 1 (MIXED corner map): logbook confirmation pending (user action); data-only pulse-shape
  discriminant queued (`papers/scripts/mixed_corner_map/`).
- GATE 2 (position arithmetic), GATE 4 (T_abs), GATE 6 (bootstrap CI) open — see
  `papers/memory_analysis_gates.md`.

## UPDATE 4 (later 2026-06-09): manuscript reframe (structural pass)
- Working tree was clean at `56dda6d` before this pass. `radical_timing.tex` reframed around the
  gated story (audit: `papers/timing/MANUSCRIPT_REFRAME_AUDIT_2026-06-09.md`): new title
  ("Detected light yield governs…"), GATE-compliant abstract/conclusions, post-fix numbers
  (203±6/440±18, floors 19–26 qualified, 25.7±0.6, "more than a factor of two"), corrected MIXED
  mechanism, NEW same-shower-control section (§sec:mixed, fig mixed_h2h_corrected) and NEW depth-dial
  section (§sec:depth, acceptance-conditional), srCFD named at definition. Pre-fix method-gain and
  systematics numbers softened/quarantined behind TODO-P2 markers pending post-fix recompute.
  Builds clean under tectonic (PDF regenerated). 11 TODO-P2 markers tracked in memory_paper2_timing.


## UPDATE 5 (later 2026-06-09): pre-fix method/systematics derivations cleared
- Working tree was clean at `0a9fd36` before this pass.
- §5.3 method gain recomputed post-fix on identical events (audit + script + table + figure under
  `papers/scripts/method_gain_postfix/`): the old 3.1 ps / 22.1→19.5 are RETIRED and replaced by the
  two-convention split — core-width gain 1.3 ps @150 (26.9→25.6), tail-sensitive gain 4.0 ps
  68% CI [3.4,4.7] (paired bootstrap; the satellite-removal effect); gain localized in saturated
  events (27.4→25.3 fully-clipped subset); floors 21.8±2.7 (cfd05, poor fit) → 19.9±0.8 (srCFD).
- tab_systematics.tex regenerated with the full production chain on identical stored events
  (`papers/scripts/systematics_postfix/`): nominals 25.7/30.3/39.6/44.4 exactly reproduce the
  authoritative table; totals ±1.0/±1.1/±0.9/±1.9; new veto-window variation rows; DSB1 floor block
  18.8±0.8 (stat) with +0.2 fit-range shift. tab_methods.tex regenerated post-fix.
- tectonic build clean. Every active number in the manuscript now traces to a post-fix/gated product.
  TODO-P2 count 10 → 7 (remaining are prose-pass items, no pre-fix numbers behind any of them).


## UPDATE 6 (later 2026-06-09): satellite-removal demonstration — last analysis TODO cleared
- Working tree was clean at `ae850ee` before this pass.
- `papers/scripts/satellite_removal/` (pre-registered AUDIT + satelliteRemoval.C): on the §5.3
  identical-event sample @150 GeV, the tail fraction in a common 66 ps window drops 4.60%→3.30%
  (46→33 of 1000; Gaussian expectation 1.24%) while robust IQR widths are nearly identical —
  the cfd05 excess is non-Gaussian tail, supporting the estimator-bias interpretation. 75 GeV
  full-fiducial split: excess concentrates in clipped events (3.39 vs 1.95%); honest REVERSAL on
  dim unclipped pulses (srCFD 7.15%, noisy LG anchor) — the empirical basis of the per-regime
  adopted-source rule, reported in figure, caption, and §5.3. Manuscript builds clean.
  All remaining TODO-P2 items (7) are writing/citation/placement; no pre-fix or ungated number
  remains anywhere in the active manuscript.


## UPDATE 7 (later 2026-06-09): full prose + bibliography pass
- Working tree was clean at `b51eef1` before this pass.
- Title decided: "A same-shower comparison of wavelength shifters and light-yield limits in a
  compact RADiCAL electromagnetic-shower timing module". Abstract, introduction (five-way
  disentanglement framing), and conclusions (three-point landing) rewritten claims-law-compliant.
  Satellite figure moved to Appendix B; depth section kept in main text (decisions recorded
  in-file). Engineering extrapolation (labeled) + per-track-vs-per-shower benchmark sentence added.
  Bibliography: 24 entries, 20 verified, 4 TODO-P2-CITE placeholders (spacal, crilin, gundacker,
  lucchini-verify). Builds clean (557 KB). Compliance greps: no retired language or pre-fix numbers
  anywhere in active prose. Remaining work is citation verification, one figure restyle, the
  cite-blocked benchmark table, and collaboration items (title/authors/funding).


## UPDATE 8 (later 2026-06-09): coauthor-circulation prep complete
- Working tree was clean at `7cd8e1f` before this pass.
- All 4 TODO-P2-CITE placeholders RESOLVED with web-verified bibliography (An NIM A 1045 (2023)
  167629; Ceravolo JINST 17 (2022) P09033 + Cantone Front. Phys. 11 (2023) 1223183; Gundacker PMB
  64 (2019) 055012; Lucchini JINST 8 (2013) P10017, byline corrected). Benchmark table built
  (tab_benchmarks.tex, object-type column, no-ranking caption) and \input in Discussion with a
  per-shower comparator sentence. Money plot restyled (comparator bands, MIXED grey-dashed
  reference) → papers/timing/figs/thesis_postfix.{png,pdf}; numbers unchanged. Circulation note +
  referee-risk memo written. Manuscript builds clean; zero citation placeholders; compliance greps
  clean. STATUS: ready for internal coauthor circulation. Outstanding for external submission:
  author/funding blocks, title sign-off, GATE-1 logbook half (recommended), optional G4 campaign.


## UPDATE 9 (2026-06-09/10): production-formatting pass (layout only)
- Working tree was clean at `8a44b35` before this pass. NO claim, number, selection, or
  interpretation changed — layout/readability only (pre-registered inventory + dispositions in
  `papers/timing/FORMAT_AUDIT_2026-06-09.md`).
- Tables: benchmark table* was 472 pt over the page edge (Scale/Comparability columns invisible)
  → six wrapped p{} columns, footnotesize, all caveats retained; Table 1 (tab:builds) 63 pt
  column overflow (σ_t^150 clipped) and Table 2 (tab_systematics) 33 pt text-overprint → both fit
  via footnotesize + tabcolsep + @{} trim. The systematics GENERATOR (systematicsPostfix.C)
  patched identically so regeneration preserves the format. `\usepackage{array}` added.
- Figures: paper convention adopted — no internal ROOT super-titles on paper-bound canvases; the
  LaTeX caption carries the title. All five campaign scripts edited + re-run; values reproduced
  identically (MIXED 1.038±0.049 & 1.038/1.161/0.795; depth −33.6±2.9 PASS; method floors
  19.9/21.8/20.3). method_postfix bottom-panel y-title/label collision fixed; depth-dial &
  method-gain 90/100 log-label collisions fixed; satellite panel-C label shortened.
  Legacy hglg.png + dist.png garbled title strips cropped (panels untouched); clip/pulse/
  optimization/systematics recorded as acceptable this pass.
- Text: "Appendix Appendix A" (elsarticle \ref duplication) fixed; mixed_h2h figure promoted to
  figure* (3 panels were illegible at single-column width).
- Build: tectonic exit 0; overfulls 4 → 1 (the 2.22 pt page-head box, cosmetic). Compliance greps
  clean. Independent page-pair verification fan-out + table-numbers integrity check run on the
  rebuilt PDF before commit.
- Round 2: the verification fan-out (5 page-pair reviewers + integrity checker) confirmed ALL
  authoritative numbers unchanged; its 2 serious findings (Fig. 7 panel-title clip; appendix
  figures illegible at column width) + annotation collisions were fixed and re-verified (PDF now
  11 pages; B.11 on its own float page, clean). DISCOVERY (flagged, not fixed — out of scope):
  Fig. A.10's baked annotations are PRE-FIX (LuAG b=19.8±5.6 vs post-fix 24.6±3.3; MIXED syst
  ±1.5 vs ±0.9) — the stability figure escaped the stale-number purge. Table 2 is authoritative;
  regeneration from systematicsPostfix.C is the first pre-submission item (follow-up task spawned;
  flagged in the circulation note).


## UPDATE 10 (2026-06-10): appendix-figure consistency fix (pre-circulation)
- Working tree was clean at `9081d76` before this pass. Consistency correction only — the gated
  numbers did not change anywhere; the appendix figure was brought into agreement with them.
- Fig. A.10 REGENERATED POST-FIX: `systematicsPostfix.C` extended to emit the stability figure
  from the same nom/shift/tot arrays as tab_systematics.tex (agreement by construction; the rerun
  reproduced nominals 25.7/30.3/39.6/44.4, totals ±1.0/±1.1/±0.9/±1.9, DSB1 floor 18.8±0.8
  identically). New figure: per-build pads, nominal star + RMS band + 9 labeled variants (veto
  rows included), only the computed total annotated, NO baked floor numbers, no super-title.
  Caption rewritten (points to Tables 1–2; drops the no-longer-true "all within band except
  K=2000"). The stale pre-fix annotations (LuAG 19.8±5.6, MIXED ±1.5, DSB1 19.5±1.1) are gone;
  paperSystematics.C superseded for paper use.
- Same pass: hglg.png regenerated from hgLgPlot.C on local RUN1075 raw (identical slopes/counts;
  panel headers into the top-margin strip, y-title inside the pad, stats stamp below the 820 mV
  band, super-title removed — the 06-09 crop workaround is obsolete); depth_dial legend fully
  inside the frame; method-gain legend nudged left.
- Build exit 0; overfulls still 1 (2.22 pt page head, cosmetic); stale-annotation grep clean
  (hits only in audit/history files); pages 3/5/7/10 re-inspected. The circulation-note warning
  about Fig. A.10 is replaced by the consistency note. PDF IS CLEAN FOR COAUTHOR CIRCULATION.
- Remaining pre-submission figure item: verify Fig. A.9 (optimization.png) scan curves against
  the post-fix chain (working points consistent; curves from the older macro).


## UPDATE 11 (2026-06-10): prior-record continuity audit (6 auditors + synthesis)
- Working tree was clean at `cb19f51`. Sources: parent NIM A 1068 (2024) 169737; Instruments 6(3)
  27 (2022); IEEE TNS 70(7) 1296 (2023); CPAD/SLAC 2023 talk; FCC/MIT 2024 talk. Multi-agent
  workflow (Agents 1–6 parallel + coordinating editor). Outputs:
  `papers/timing/PRIOR_RECORD_{AUDIT,CONTINUITY_MATRIX,ACTION_PLAN}_2026-06-09.md`.
- VERDICT: "circulate after minor edits." Thesis, MIXED control, headline numbers, and the
  confirms-never-revises posture survived all six audits. Duplicate-publication risk LOW: the
  parent's Sec. 7 item 2 pre-registered both the WLS comparison and the DSB1 dataset reuse.
- KEY DISCOVERY: the claims-law "27–29 ps full fiducial" was pre-fix-era/unverified. Verified on
  the locked production chain (`papers/scripts/full_fiducial_check/fullFiducialCheck.C`): the
  full-fiducial σ_t(150) is ≈50 ps (50.5 ps, 1.45×10⁵ events; K=1000 reproduces 25.7 exactly).
  Law amended (27–29 RETIRED); companion now quoted in abstract, §4, §5.3, Conclusions.
- MUST-FIX edits landed (wording-level; no gated number changed): cfd05 → "clipped-peak
  discriminator class of Ref. [1]" + Table-1 caption rewrite (Q7 delegated); satellite
  sum-vs-differential/digitiser-vs-clipping qualification; LG digitised at 1 GS/s (DRS1) split
  from the 5 GS/s HG statement; "16 SiPMs" → "eight SiPMs (sixteen readout channels)" (+benchmark
  table); light-deficit ≈3× measured (not the a-ratio 2.2) with a∝1/√N → a/√2≈310 arithmetic and
  the Section-7-item-3 pointer fix. Pb-glass veto verified NOT in the timing chain → documented,
  no manuscript edit (Q8). Cheap SHOULD edits folded in: same-dataset disclosure (S1), GEANT4
  expectation citation (S2), per-regime "attempted uniformly" fix (S4), "designed to achieve"
  benchmark verb, Section-7-item-2 pointer, header [2]→[1].
- Claims law: +5 prior-record rules (estimator attribution, satellite provenance, talk quarantine,
  light-deficit arithmetic, apparatus facts); stale in-law numbers fixed (42→44.4, 19.8±5.6→
  24.6±3.3, 25.3→25.7). Circulation note: talk-facing mapping section + don't-reinstate bullets +
  Q7/Q8. Referee memo: risks 11–13 added, defended in-text. Stale figs/thesis.png + method.png
  deleted (pre-fix, unreferenced).
- Build exit 0; overfulls still 1 (2.22 pt cosmetic); compliance grep clean. SHOULD tier deferred
  to pre-submission: apparatus figure (parent Fig. 1 adaptation), money-plot published-fit
  reference curve + per-point errors, TNS lineage, T/E-type vocabulary, depth Tier-2 wording,
  HGTD/cantone2023/OSTI/byline verifications, A.9 post-fix curve check.


## UPDATE 12 (2026-06-10): apparatus composite inserted (reuse-report insertion #1)
- Working tree was clean at `6ee4634` (only the user-added NIM_A_Figures/ untracked).
- NEW FIG. 1 (`fig:apparatus`, figure*, Sec. 2.1): (a) parent Fig. 1 module schematic +
  (b) parent Fig. 6 T-type schematic, both embedded as-is with minimal overlays (beam arrow,
  UP/DW tags, capillary/WLS-filament/SiPM callouts); (c) beam's-eye corner map REDRAWN from the
  parent Fig. 2 engineering numbers with the author-confirmed channel-map names and the GATE-1
  MIXED diagonals (NE+SW=DSB1 filled, NW+SE=LuAG open; grayscale only). Caption: defines
  geometry, previews the (DW−UP)/2 construction and the same-shower logic, carries "adapted
  from Ref. [1]" + the the-schematic-is-not-the-source caveat for the material map.
- ORIENTATION HONESTY: the absolute beam's-eye compass orientation is NOT pinned by the repo
  record (channel map = electronics-only; position estimator = trained fit) → drawn as a
  DECLARED convention (N up, E right, looking downstream), stated in panel, caption, and audit;
  no claim depends on it (diagonal pairing is what carries Sec. 6); physical confirmation folded
  into coauthor Q6. UP/DW left-right assignment anchored by the filament position (shower max
  is upstream) in the parent drawing.
- Sec. 2.1 text: 4 sentences added (parent remains the detailed apparatus reference; figure
  defines geometry + corner convention; MIXED logic rests on the pulse-shape-confirmed map).
  No numerical claims added.
- Generator: papers/scripts/apparatus_composite/apparatus_composite.tex (standalone TikZ →
  vector PDF; PNG preview via sips). Assets: papers/timing/figs/radical_apparatus_composite.{pdf,png}.
  Audit: papers/timing/PARENT_APPARATUS_FIGURE_INSERTION_AUDIT_2026-06-10.md. Build exit 0,
  12 pages, overfulls still 1 (2.22 pt cosmetic); figures auto-renumbered (clip→Fig 2,
  satellite→B.12); compliance grep clean.
- GEANT4 longitudinal-profile insertion (#2) DEFERRED to pre-submission (pending Sec.-7 support
  needs after coauthor review). Earlier this session: author block completed from the parent
  roster (39 authors, 11 ORCIDs verified against the parent PDF link annotations; Wetzel first,
  Perez-Lara alphabetical at P) and the parent-figure reuse report committed.


## UPDATE 13 (2026-07-21): workspace world-class audit — Phase A executed
- Ten-agent panel audit committed (`WORKSPACE_AUDIT_2026-06-10.md` + REORG_PLAN; verdict
  C+ trending B-: interior record A-, stranger-facing shell D). User decisions: transparency
  perimeter KEPT (memos/panel docs publish as-is); Phase A executed in full.
- A1 full-fiducial evidence gap CLOSED: fullFiducialCheck.C rerun (identical: 50.5 ps @150
  full-fiducial, K=1000 reproduces 25.7), log committed, AUDIT.md added honestly labeled
  post-hoc. A2 root README rewritten to the claims-law headline + trust path + complete map;
  ANALYSIS_GUIDE.md added (waterfall + full provenance table, srCFD<->hg_lgcfd bridge,
  exploratory-vs-load-bearing split). A3 chain_of_evidence.html fixed (8 waveform PNGs
  vendored into figures/2023/narrative/chain/; renders complete on fresh clone). A4 external
  working memory imported verbatim with provenance headers: data/2023/metadata/DATASET_NOTES.md,
  docs/APPARATUS_2023.md, docs/METHODS_2023.md. A5 papers/README.md (master index + series-
  numbering declaration: timing = Paper 1; "Paper 2" alias explained) and papers/timing/README.md
  (actual 12-figure map + document index) rewritten. A6 papers/timing/figs/MANIFEST.md
  (md5 + generator + copy method per figure; dist/optimization frozen-as-circulated). A7
  DEPRECATED banners on paperSystematics.C + methodCompare.C + breadcrumb beside the stale
  narrative systematics.png. A8 HISTORICAL-PAGE banners on site/paper.html, timing_story.html,
  report/index.html. A10 manuscript rebuilt from committed sources; apparatus_composite build
  artifact + NIM_A_Figures (rights pending) + .DS_Store gitignored; push + tag
  `circulation-2026-06`.
- Phase B (pre-submission) remains: LICENSE/CITATION/Zenodo data DOI + fetch script +
  REPRODUCE.md, studies/gates indexes, metadata reconciliation, A.10-era optimization-curve
  re-verification, ARCHITECTURE/CONVENTIONS docs. Follow-up experts warranted: open-data/
  publishing-rights specialist; hardware coauthor (Q8); user logbook scan (GATE 1 half).


## UPDATE 14 (2026-07-21): code-audit quick-win list executed
- Six-expert code/analysis audit (CODE_AUDIT_2026-07-21.md: A-, SoC ~125 min measured, zero
  code-vs-paper mismatches, byte-identical reproduction) followed by the full quick-win list:
- (1) RadTiming.h header rewritten — timingBrightestK+kLGCFD (srCFD) declared THE production
  headline (25.7±0.6); timingBestBin marked RETIRED at its definition; stale tebSigma lead-in
  and dead paths fixed; sigmaT.C header corrected. (2) KERNEL-CLONE cross-reference comments on
  all 7 gate-macro reimplementations of eventDWUP/tebSigma, each naming its deliberate deltas
  and the sigma_K1000(150)=25.7 equivalence proof. (3) SelectionCuts.h summary rewritten as
  WHICH-CHAIN-APPLIES-WHAT (production timing = 1,2,3,5,6+veto+brightest-K; Pb/M5/M7 tagged
  study-era/retired); kMCP2*/kHG_maxPeak annotated not-applied-in-production. (4) Reproduction
  protocol: run-from-repo-root + expected PDF-timestamp residual documented in ANALYSIS_GUIDE +
  timing_fit_summary AUDIT addendum (numeric-record vs environment-noise split of the log).
  (5) srCFD alias at all three code definition sites + METHODS_2023. (6) docs/STATS_CONVENTIONS.md
  (two width conventions, sqrt(2N), PDG scaling, RMS-of-shifts + quadrature alternative
  2.8/3.3/2.4/5.5 with conclusions-survive paragraph, bootstrap, MIXED scatter); methodGainPostfix
  now prints the un-resampled tail-sensitive centrals — gate RERUN, diff = intended lines only,
  ALL numbers identical, and the paper's 4.0 ps tail-sensitive gain now sits in the committed
  table (evidentiary hole closed; LED tail-sensitive 4.3 ps also recorded); sentinel -1.0 -> n/a.
  (7) lgcfd_frac=0.15 WHY at BuildConfig.h + Reducer.C incl. the 20/780 mV guard rationale.
  (8) BuildConfig D4U4-order guard (protects eventDWUP's index invariant); systematicsPostfix
  FAIL-print insurance; fullFiducialCheck r=3.0-protocol AUDIT addendum; depth-dial fit-form
  addendum; tebSigma 0.9546-vs-iterated-0.938 fallback-debias caveat (comment only, post-
  submission item); ChannelConfig authority note; stale Analysis/ path tags purged.
- MANUSCRIPT (disclosure, no numbers changed): Sec. 4 now defines the in-event 2 ns median
  veto and discloses the per-energy fiducial ramp (2.5 mm ≤100 GeV -> 3.0 mm ≥125); App. A
  "chosen 3.0 mm" reworded. Rebuilt clean (one 2.22 pt cosmetic overfull unchanged).
- All 6 edited gate macros compile clean in isolation; expected post-quick-win SoC ~80-90 min.


## UPDATE 15 (2026-07-21): Phase B mechanical batch — reproduction tooling + indexes
- Documentation/tooling only; no numeric change, no gate-macro bytes touched. Commits
  `a7c400e` (batch) + `3c81e7f` (repro.sh verdict fix) on this machine (OWC checkout,
  ROOT 6.36.04/macosxarm64 — NOT the 6.40.00 machine of UPDATE 14's E4 run).
- SHIPPED: tools/repro.sh (D1: one-command headline-chain reproduction, PDF-timestamp
  exemption built in); tools/run_gate.sh (B3 wrapper: refuses without data, clobber
  warning, bash-3.2-safe); tools/check_figure_manifest.sh (B8 drift checker — run twice
  this pass, 14/14 OK both times); papers/scripts/INDEX.md (gate index + runbook);
  analyze/README.md + analyze/studies/INDEX.md (~120 macros, one line each, from their
  own headers; only hgLgPlot.C + narrativeFigs.C load-bearing; 3 DEPRECATED banners
  confirmed m–z, none a–l); docs/ARCHITECTURE.md + docs/CONVENTIONS.md (B4);
  lib/physics/ESTIMATOR.md (D3: kernel walk + naming bridge + SW-U/mcp_ref-2
  cancellation); data/2023/MANIFEST.sha256 (B2: sha256+bytes+entries, 21 ntuples +
  8 sidecars, generated programmatically, `shasum -c` verified 29/29 OK); REPRODUCE.md
  (env pinning, venue-pending placeholders); README/ANALYSIS_GUIDE pointers;
  .gitignore tracked-log-exception comment.
- VERIFICATION RUN (new evidence): `tools/repro.sh --check-data` executed twice on this
  machine. Final: STEP 0 manifest 29/29 OK; STEP 1 verify.C PASS; STEP 2
  timingFitSummary exit 0; STEP 3 fullFiducialCheck stdout IDENTICAL to committed log
  (diff -B); STEP 4 = NUMERIC PASS: papers/tables/timing_fit_summary_2026-06-09.md
  byte-identical, but the three figure renders (timing_fit_summary.png,
  thesis_postfix.{png,pdf}) came out byte-different under ROOT 6.36.04 (~0.3% size;
  encoder/writer version). FINDING RECORDED in REPRODUCE.md: image byte-identity is
  toolchain-pinned (6.40.00); numeric records are version-robust. Figures restored from
  HEAD (`git checkout --`), tree clean, figs/MANIFEST unchanged. repro.sh exit code 3
  now names this outcome; first-run false-positive (blank-line capture noise in the
  full_fiducial log compare) fixed with diff -B. Total 229 s.
- FINDING (deferred to B5 metadata reconciliation): data/2023/README.md still describes
  two ntuple formats; measured today all 21 reduced files carry the single 39-branch
  'rad' schema. Also: tectonic absent on this machine (manuscript builds happen on the
  other machine or after install).
- Phase B remaining: blocked-on-decisions B1 (license/CITATION), B2 venue (Zenodo vs
  CERN Open Data) + fetch script, B5 (metadata reconciliation), B6 (banner batch),
  B7 (GATE-1 logbook half, Q8), B8 tab_methods/NIM_A dispositions, B9 tag; plus the
  A.10 optimization-curve re-verification.

## UPDATE 2026-07-21 (late) — report pipeline migrated to the production chain (banner retired)

- WHAT: the full report pipeline (harvester + ~30 analyze/studies macros + makeReport.py +
  deployReport.sh) migrated off the retired best-bin/cfd05 chain onto the gated production
  chain (brightest-1000 (DW−UP)/2, srCFD/LED per regime, TimingFiducialR ramp, eventDWUP veto,
  robust tebSigma). site/report/ redeployed banner-free with current numbers.
- NEW PRODUCER: analyze/studies/timingProduction.C — writes Summary/timing_energy_bins.root
  under the legacy graph names + per-build TParameters + full-fiducial/cfd05 companions. Carries
  a KERNEL CLONE equivalence assertion vs rad::timingBrightestK and a GATED CHECK printout.
  Run evidence (output/logs/rerun_timingProduction.log): "GATED CHECK: s150 PASS (25.7 vs 25.7)
  | FF150 PASS (50.5 vs 50.5) | a PASS | b PASS"; all four builds reproduce the gated table
  exactly (DSB1 a=203±6 b=18.8±0.8 s150=25.7 light150=5185; LUAG 440±18/24.6±3.3/44.4/1800;
  MIXED 253±29/34.3±2.4/39.6/3230; TENERGY 198±21/26.0±1.8/30.3/4045).
- DATA: reduced ntuples symlinked as output/<E>GeV/ntuple.root (identical event sample,
  n=1,826,512 @150). Comprehensive branch-bind audit: mcp_time/mcp_peak → dynamic
  mcp1_* fallback in lib/viz/PlotUtils.h (ScanRunCenters ×2), alignmentAnalysis.C, + 37-macro
  mechanical sweep (silent-zero SetBranchAddress hazard eliminated).
- RERUNS: Batch A (18 QA macros) + Batch B (12 consumers/heroes) all exit=0 errors=0
  (output/logs/rerun_*.log). Rewritten data-driven: timingFloorComparison.C (was hard-coded
  "22 vs 17.5" story → 18.8±0.8 CONFIRMS 17.5), idealUniform.C (was hard-coded 27.4 ladder →
  reads production curve; projection now 25.7→23.7 ps @150). Label fixes: elbowFractionTrend,
  edgeMechanism, layer4Summary, layer5Summary (srCFD is the operating point, not "CFD-5% adopted").
- HARVESTER: harvestResults.C emits fit errors (PDG inflation), dsb1_sigma_ff[], per-build
  block (from timingProduction TParameters), mixed_sameshower_ratio 1.04±0.05, depth_slope
  −33.6±2.9, syst_sel (gate-log provenance). results.json sha256[:12]=bdb57d97668e.
- makeReport.py: 24 pre-fix passages rewritten under the claims law (exec summary wholesale,
  l5-floor section now "confirming the published 17.5 ps", fiducial-opt reframed as the
  historical radius derivation, per-regime estimator rule throughout, WLS-capillary-is-the-
  variable statement, per-build table, generated-on dateline). Claims-gate greps on
  output/report.html: 27.4=0, "22 vs 17.5"=0, 0.99=0, "identical timing"=0, "crystal differs"=0;
  required: 25.7×16, brightest-1000×27, srCFD×23, "confirming the published 17.5"×2.
  report.html sha256[:12]=77b5cdffe926.
- DEPLOY: bash analyze/deployReport.sh → 219 files, 36 MB; deployed index sha256[:12]=77b5cdffe926;
  verified in-browser (KPIs 25.7 / ≈50 / 18.8; 244/244 figures load; hero + floor figures
  visually inspected). Hub card (site/index.html) updated to "production chain · 2026-07".
- COMMANDS: root -l -b -q 'analyze/studies/timingProduction.C+' ; batch scripts (nohup) ;
  root -l -b -q 'analyze/studies/harvestResults.C+' ; python3 analyze/makeReport.py ;
  bash analyze/deployReport.sh.

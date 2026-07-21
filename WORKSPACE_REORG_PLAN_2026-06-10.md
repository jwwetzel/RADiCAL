# Workspace Reorganization Plan — 2026-06-10

Synthesis of the ten-agent panel (full reports: `WORKSPACE_AUDIT_2026-06-10.md`).
NOTHING in this plan has been executed; it awaits sign-off on the target tree.

## 1. Verdict — how far from world-class

**Distance to world-class: about two weeks of purely additive work — none of it touching a gated byte — but today the workspace would fail its own stated goal.** Overall grade: **C+ trending B-** (interior analysis record: A-; stranger-facing shell: D). The panel is unanimous on the diagnosis: this is a two-layer object whose inner layer (gates, claims law, pre-registered audits) is better than normal HEP collaboration practice, wrapped in an outer layer (front doors, data access, licensing, publication state) that is stale, contradictory, or simply not built. Nothing scientific needs rework; everything a stranger touches first does.

**Strongest assets (verified by multiple experts independently):**
1. The gate layer — papers/scripts/*/ with genuinely pre-registered AUDIT.md files, committed result logs, and papers/tables/timing_fit_summary_2026-06-09.md as a true single source of truth; the headline 25.7±0.6 ps traces script→log→table→figure→tex with byte-identical figure copies (E1 md5-verified).
2. A clean, verified dependency spine — lib ← reduce ← analyze ← papers/scripts, one schema, one cuts file (SelectionCuts.h with WHY-comments — "the best file in the repo," E8), one path resolver, bit-identity reduction gate, and a manuscript that rebuilds to within 4 bytes of the committed PDF (E6 ran it).
3. Honest falsification in public — the 0.99 kill-shot retired via commit + DEPRECATED figure + breadcrumb README; claims law and append-only snapshot governance; seeded bootstrap.
4. Primary provenance in-repo — the 1278-row beam logbook.csv, per-run manifests with bias, config-driven reduction with committed .hglg calibration sidecars.

**Widest gaps:**
1. Every entry surface states retired physics: root README headlines the forbidden "≈27 ps," the public site says 27.4/25.3 ps (zero mentions of 25.7), curated analyze/ macros advertise the retired estimator, and the README repo map omits papers/ — the layer where the entire scientific record lives.
2. The release shell is 0% built: no LICENSE, no CITATION.cff, no tags, no data DOI, a literal `<PASTE_CERNBOX_SHARE_LINK_HERE>` placeholder, and 31 unpushed commits — the public repo (github.com/jwwetzel/RADiCAL) ends at the pre-gate era.
3. Load-bearing knowledge lives outside the repo in an orphaned ~/.claude memory dir (authoritative 36-slot channel map, run counts, valid eras, jitter provenance), while in-repo metadata contradicts code on a method-relevant fact (channel_map.yaml: 5 GS/s flat vs the author-confirmed DRS1@1GS/s).
4. Provenance holes at the paper's own convention: the mandatory 50.5 ps companion has no AUDIT.md and no committed log; Fig. 4 (dist.png) is byte-orphaned; Fig. 1's inputs are 88 untracked Elsevier-figure files; chain_of_evidence.html renders broken on any clone.
5. Live hazards: two unmarked superseded macros still write the gated manuscript tables, and a data-less rerun of the money-table script silently overwrites the authoritative table with zeros (E6 verified this live).

## 2. Scoreboard

| Dimension | Expert | Grade | One-line reason |
|---|---|---|---|
| Data-analysis strategy | E1 | B- | Gate layer exemplary and headline fully traced; but every entry surface contradicts it and two figures + the mandatory companion number have provenance gaps |
| Calorimetry / physics organization | E2 | B- | Cuts/estimators/gated production world-class and timing-primary; three physics eras coexist unlabeled, estimator has 3 names and no in-repo home, energy has no publication-grade home |
| Instrumentation | E3 | B- | Consumable layer (code+configs+manifests) is A-grade; authoritative 36-slot map lives in orphaned external memory and tracked metadata contradicts code on sampling rates |
| Paper review / figure custody | E4 | B- | Numbers chain passes; figure chain fails the "seconds test" — undocumented renames, one byte-orphaned figure, untracked Fig. 1 inputs, three numbering schemes |
| Project strategy / release readiness | E5 | C+ | Interior governance world-class; release shell (LICENSE, CITATION, DOI, tags, data path, availability statement) literally not started; 31 commits unpushed |
| Reproducibility | E6 | C- | Fresh-clone experiment: paper rebuilds to 4 bytes, but zero data access, zero env pinning, and the documented money-table command verifiably clobbers the committed evidence with zeros |
| Software architecture | E7 | B- | Include-graph spine clean and verified; write-side back-edges from unmarked superseded macros into gated tables, manual figure custody with a verified stale-name collision |
| Fresh eyes / onboarding | E8 | C- | Trusts the 25.7 ps after one hour — but only by grepping past five front doors that each stated retired numbers, methods, or dead paths |
| **Overall (publishable-workspace goal)** | **Panel** | **C+ → B-** | **A-grade audit spine, D-grade shell; the goal weights the stranger's path, and the stranger's path is currently broken at every door** |

## 3. Target structure

```
RADiCAL/
├── README.md                    REWRITE: claims-law headline (25.7±0.6 ps brightest-1000 srCFD + ≈50 ps
│                                full-fiducial companion); repo map gains papers/, tools/, output/,
│                                chain_of_evidence.html; 10-line "How to trust this result" path;
│                                correct the "reduced data in the repo" fiction
├── ANALYSIS_GUIDE.md            NEW (the missing spine): raw → reduce/ (REREDUCE.md) → data/2023/reduced
│                                (39-branch 'rad') → lib/physics/{SelectionCuts,RadTiming}.h →
│                                papers/scripts gates (AUDIT.md discipline) → papers/tables/ →
│                                radical_timing.tex; plus figure/number → generator → log → commit table
├── REPRODUCE.md                 NEW: 5-command recipe (setup → fetch → verify.C → timingFitSummary.C →
│                                tectonic), env pins (ROOT 6.40.00 tested/macOS; tectonic), runtimes,
│                                expected outputs, and the WARNING that gate reruns overwrite tracked
│                                products in place (verify via git diff — never commit over the record)
├── LICENSE / CITATION.cff / .zenodo.json   NEW (Phase B; license + authorship are user/collab decisions)
├── chain_of_evidence.html       REPOINT ONLY: 8 <img> → tracked copies of the existing waveform PNGs,
│                                vendored (byte-copied, never regenerated) into figures/2023/narrative/chain/
├── docs/
│   ├── ARCHITECTURE.md          NEW: layer diagram + verified no-back-edge claim; intra-lib waiver;
│   │                            the THREE figure-output systems as declared policy; cwd=repo-root rule;
│   │                            srCFD ↔ hg_lgcfd ↔ RadView kLGCFD name map; frozen-gate-code doctrine
│   ├── CONVENTIONS.md           NEW: corner labels + viewpoint; capillary D/U (downstream/upstream) vs
│   │                            wire-chamber Down/Up collision; WC scales; kCap order; HG/LG apparatus
│   │                            wording pending coauthor Q8
│   ├── APPARATUS_2023.md        NEW: imported external-memory compendium (MCP model, DAQ-settings
│   │                            narrative → 820 mV clip, measured 71–74 ps jitter provenance)
│   └── PAPER1_DRAFT.md etc.     SUPERSESSION BANNERS (the EXPERT_PANEL_BLUEPRINT.md house pattern)
├── data/2023/
│   ├── README.md                REWRITE: single 39-branch schema (retire the two-format text),
│   │                            DOI-backed access replacing <PASTE_CERNBOX_SHARE_LINK_HERE>, real sizes
│   ├── MANIFEST.sha256          NEW: path, bytes, sha256, tree, entry count — 21 reduced files,
│   │                            8 configs/.hglg, 20 local raw
│   ├── metadata/channel_map.yaml  FIX: full 36-slot map (WC anodes, trigger scints, logbook labels)
│   │                            + per-board sampling DRS0@5GS/s / DRS1@1GS/s (kill the flat 5.0)
│   ├── metadata/dataset.yaml    REGENERATE: gated headline, live paths, .hglg mention
│   ├── metadata/DATASET_NOTES.md  NEW: run counts, valid eras (MIXED 2975–3044), light levels —
│   │                            imported from the orphaned ~/.claude memory (repo becomes the authority)
│   └── configs/MIXED.json       ANNOTATION-ONLY amend: _material_map_source cites the GATE-1 8/8
│                                pulse-shape gate, "logbook half pending" retained; values untouched
├── reduce/hpc/README.md         BANNER: HISTORICAL — REREDUCE.md is the sole current recipe
├── analyze/
│   ├── README.md + studies/INDEX.md  NEW: all 118 macros classified LOAD-BEARING (the named ~10-macro
│   │                            set with what each feeds) / SUPERSEDED-BY (successor path) / EXPLORATORY
│   ├── studies/paperSystematics.C, methodCompare.C  DEPRECATION BANNERS (house mixedHeadToHead.C
│   │                            pattern) — both still write papers/timing/tab_*.tex and must be defused
│   └── sigmaT.C, timingHeadline.C  COMMENT-ONLY header fixes, post-tag, "documentation only" commits
├── figures/2023/narrative/      ADD breadcrumb README beside the stale pre-fix systematics.png
│                                (rename to *_DEPRECATED optional in Phase C — verified uncited)
├── papers/
│   ├── README.md                REWRITE as master index: governing record (claims law, gates, snapshot),
│   │                            canonical table, 12-figure map, review-history list, paper-numbering
│   │                            declaration (series scheme: parent=I, timing=II, energy/position=III),
│   │                            fix the 3 dead external-memory references
│   ├── scripts/INDEX.md + README.md  NEW: gate ledger (status/script/AUDIT/log/figure/consumer) +
│   │                            runbook (order, exact commands, frozen-code rule, clobber warning)
│   ├── scripts/full_fiducial_check/  ADD honest post-hoc AUDIT.md (dated today, labeled retrospective)
│   │                            + rerun teeing full_fiducial_result.log (verified stdout-only — safe)
│   ├── scripts/apparatus_composite/  ADD post-hoc AUDIT.md; track output PDF; vendor the 2 parent-paper
│   │                            PDFs it includes (pending rights ruling) so Fig. 1 builds from a clone
│   ├── figures/ + tables/       UNTOUCHED — the gated record
│   ├── timing/figs/MANIFEST.md  NEW: 12 includes → generator → source → md5 → copy-method; dist.png and
│   │                            hglg.png documented as FROZEN-AS-CIRCULATED (no regeneration pre-submission)
│   ├── timing/README.md         REWRITE against the actual tex (real figure map; drop stale placeholder note)
│   ├── timing/review_history/   OPTIONAL Phase C git mv of the 9 dated process docs + OUTLINE/THREAD,
│   │                            recorded as a snapshot UPDATE with an old→new path table
│   └── memory_*.md, CAMPAIGN_SNAPSHOT, tables/, all AUDIT.md   APPEND-ONLY, never edited; date-drift
│                                fixed by ADDED header lines ("opened 06-09; updates through 06-10 inline")
├── site/                        BANNER or regenerate paper.html / timing_story.html / report/ BEFORE any
│                                push — the Pages workflow auto-deploys whatever lands on main
├── tools/fetch_data.sh          NEW: DOI-backed download + sha256 verification (stranger route);
│                                pull_argon_data.sh relabeled maintainer-only
└── git                          push all 31 commits AS-IS (no rebase) after site fix + perimeter decision;
                                 commit the rebuilt radical_timing.pdf; tag circulation-2026-06, later
                                 submission-v1 (each snapshotted by a Zenodo version DOI)
```

## 4. Migration plan (phased, risk-ranked)

**PHASE A — pre-circulation (days; coauthors are the first strangers).**
A1. SAFE — Close the full_fiducial gap: rerun `papers/scripts/full_fiducial_check/fullFiducialCheck.C` teeing stdout to a committed `full_fiducial_result.log` (verified stdout-only; touches no gated product), add an AUDIT.md honestly labeled "written post hoc 2026-07, documenting the verification run" — never backdated. Confirm via `git diff` that only new files appear.
A2. SAFE — Rewrite root `README.md` (claims-law headline, complete repo map incl. papers/, trust-path 10-liner, corrected data statement). Write `ANALYSIS_GUIDE.md` (the waterfall + provenance table).
A3. SAFE — Vendor-copy the 8 existing `output/<BUILD>/{hg,lg}_waveforms.png` into `figures/2023/narrative/chain/` and repoint `chain_of_evidence.html` (byte copies of existing images; regenerate nothing).
A4. SAFE — Import external memory: distill `~/.claude/.../memory/dataset_2023.md` (+ detector/methods invariants) into `data/2023/metadata/DATASET_NOTES.md` and `docs/APPARATUS_2023.md`; fix the 3 dead memory references in `papers/README.md`.
A5. SAFE — Rewrite `papers/README.md` (master index + series-numbering declaration) and `papers/timing/README.md` (actual 12-figure map).
A6. SAFE — Write `papers/timing/figs/MANIFEST.md` recording CURRENT bytes (md5), generators, copy-methods; dist.png and hglg.png entered as "frozen as circulated, manual provenance per FORMAT_AUDIT F7."
A7. SAFE — Defuse the clobber macros: house-pattern deprecation banners on `analyze/studies/paperSystematics.C` and `methodCompare.C`; breadcrumb README beside the stale `figures/2023/narrative/systematics.png`.
A8. CAREFUL — Site: banner (minimum) or regenerate `site/paper.html`, `timing_story.html`, `site/report/` — MUST land in the same push as everything else (Pages auto-deploys).
A9. USER DECISION — Transparency perimeter: REFEREE_RISK_MEMO, THREAD.md, panel JSONs go public on push (remote is public GitHub). Decide deliberately before A10.
A10. CAREFUL — Publish the record: rebuild `radical_timing.pdf` once from committed sources and commit it; decide NIM_A_Figures interim disposition (default: keep untracked + pointer note pending rights answer); commit or ignore apparatus_composite.pdf; push all 31 commits AS-IS (no rebase — hashes are cited in audits); tag `circulation-2026-06`.

**PHASE B — pre-submission (1–2 weeks).**
B1. LICENSE + CITATION.cff + .zenodo.json (after user/collaboration decisions on license and code-release authorship).
B2. Data tier: Zenodo deposit of the ~13 GB reduced set (+configs/.hglg), `data/2023/MANIFEST.sha256`, `tools/fetch_data.sh`, `REPRODUCE.md`, replace the CERNBox placeholder, add Data & Code Availability to the tex.
B3. Indexes: `analyze/studies/INDEX.md` + `analyze/README.md`; `papers/scripts/INDEX.md` + gate runbook (with the explicit rerun-clobber warning; optionally a `tools/run_gate.sh` wrapper that refuses to run without data — the gate macros' bytes stay untouched).
B4. `docs/ARCHITECTURE.md` + `docs/CONVENTIONS.md` (fold the estimator taxonomy in; srCFD↔hg_lgcfd bridge).
B5. CAREFUL — Metadata reconciliation: channel_map.yaml (36-slot, per-board rates), dataset.yaml regeneration, registry.csv/MANIFEST.csv `datasets/`→`data/` prefixes, configs' runs_manifest paths; MIXED.json annotation-only amend citing the GATE-1 pulse-shape gate with "logbook half pending" retained.
B6. Banners + comment-only header fixes: reduce/hpc/README, docs drafts; analyze/sigmaT.C, timingHeadline.C, lib/physics/RadTiming.h:6 — post-tag, each commit labeled "documentation only, no numeric change."
B7. USER ACTION — Close GATE 1's logbook half (scan the assembly page); then amend MIXED.json annotation and the gates file by dated append. Reconcile tex:224 "per-cell calibration" vs DRS4Calibration.h through the claims-law amendment process. Resolve coauthor Q8 (HG/LG architecture) and update apparatus prose everywhere at once.
B8. tab_methods.tex decision (input it or move to papers/tables/ with a "generated, intentionally not included" header); NIM_A_Figures final disposition per rights ruling; `tools/check_figure_manifest.sh` drift checker; `.gitignore` comment documenting the tracked-log exception.
B9. Tag `submission-v1`; mint the Zenodo version DOI against it.

**PHASE C — nice-to-have / post-publication.**
C1. Optional `git mv` of the 9 dated process docs into `papers/timing/review_history/` (snapshot UPDATE with path table) — riskiest mechanical move proposed; cosmetic value only.
C2. HEPData record for the sigma_t(E) tables; Dockerfile or environment.yml.
C3. Full site regeneration replacing banners.
C4. Energy backbone gate home (open GATE 7): `papers/scripts/energy_backbone/` with pre-registered AUDIT for Paper III; populate `papers/energy_position/figs/`.
C5. Post-`v1-paper` refactors ONLY: lib/core/ split, RadView::hg_saturated() accessor, shared gatherProduction() for 2024+ campaigns, StampBuild provenance extension for future reductions.
C6. Optional regeneration of dist.png from a committed script as a NEW, audit-logged act — a post-circulation decision, never a silent repair.
C7. Optional rename of the stale narrative systematics.png to the *_DEPRECATED pattern.

**WHAT MUST NOT BE DONE (union of all eight, binding):** never edit existing AUDIT.md files, committed result logs, memory files, or the snapshot (append-only); never regenerate or renumber any gated number or the current papers/timing/figs bytes (including dist.png) before submission; never modify code bytes of papers/scripts/*/ macros (input guards ruled out — wrapper + runbook instead); never rename/move analyze/studies macros or papers/scripts dirs; no rebase/squash of the 31 commits; never rename the hg_lgcfd/lgcfd branch or Schema layout; never delete the DEPRECATED breadcrumbs, the tracked evidence .log files, or the MIXED.json caveat (amend by dated append only); never backdate a new AUDIT.md.

## 5. Publication checklist

**USER / COLLABORATION DECISIONS (blocking, cannot be done mechanically):**
- [ ] License split: code (MIT or BSD-3) vs figures/prose (CC-BY-4.0) — and collaboration sign-off.
- [ ] Code-release authorship for CITATION.cff: Wetzel + contributors vs the full 39-author roster.
- [ ] Transparency perimeter: publish REFEREE_RISK_MEMO, THREAD.md, and panel JSONs deliberately (declared in README) or move to a private attic before the public push.
- [ ] NIM_A_Figures (88 files, 33 MB Elsevier parent-paper figures): redistribute with LICENSE_NOTE, vendor only the 2 files Fig. 1 needs, or exclude with a citation-only manifest — pending rights expert.
- [ ] GATE 1 logbook half: scan the MIXED assembly logbook page; then dated-append closure to gates file + MIXED.json annotation.
- [ ] Coauthor Q8: HG/LG = gain-split SiPM copies vs separate WLS optical paths — gates every apparatus-doc rewrite.
- [ ] Data venue choice: Zenodo vs CERN Open Data for the ~13 GB reduced ntuples (+ HEPData for sigma_t(E) tables).
- [ ] DT5742 run-control settings recovery from Argon/EOS (if they exist) for the apparatus compendium.

**MECHANICAL (executable once decisions land):**
- [ ] LICENSE, CITATION.cff, .zenodo.json at root.
- [ ] Zenodo deposit of reduced ntuples + configs/.hglg + logbook.csv; concept DOI + version DOI per tag.
- [ ] data/2023/MANIFEST.sha256 (path, bytes, sha256, tree, entries) + tools/fetch_data.sh with hash verification.
- [ ] Replace `<PASTE_CERNBOX_SHARE_LINK_HERE>`; rewrite data/2023/README.md to the single 39-branch schema.
- [ ] Data & Code Availability section in radical_timing.tex (currently zero matches for availab/zenodo/github).
- [ ] REPRODUCE.md with env pinning (ROOT 6.40.00 tested / macOS; tectonic; Linux-untested honesty line).
- [ ] README trust path + ANALYSIS_GUIDE.md; papers/README master index; figs/MANIFEST.md.
- [ ] Post-hoc AUDIT.md + committed log for full_fiducial_check; AUDIT.md for apparatus_composite.
- [ ] chain_of_evidence.html image vendoring (renders complete on a clone).
- [ ] Import external ~/.claude memory content → data/2023/metadata/ + docs/ (repo becomes the authority).
- [ ] Site banner/regeneration before push; push 31 commits; tags `circulation-2026-06` then `submission-v1`.
- [ ] Commit rebuilt radical_timing.pdf; .gitignore comment documenting the tracked-evidence-log exception.
- [ ] INDEX files (analyze/studies, papers/scripts), supersession banners, path-rot batch (registry.csv, dataset.yaml, channel_map.yaml, runs_manifest fields, stale headers).

## 6. Where experts disagreed — rulings

**1. Input guards inside gate macros (E6) vs frozen gate bytes (E7).** E6 verified live that a data-less rerun of timingFitSummary.C silently zeros the committed authoritative table; E7 forbids any byte change to papers/scripts macros. RULING: E7 wins — committed-code≡committed-numbers is the workspace's core evidentiary property. The hazard is mitigated without byte edits: a loud warning in REPRODUCE.md and the gate runbook ("reruns overwrite tracked products; verify by git diff; never run without data") plus an optional `tools/run_gate.sh` wrapper that refuses to invoke root when `data/2023/reduced` is absent.

**2. Regenerate dist.png now (E1) vs freeze the circulated bytes (E4, E7).** RULING: freeze. The figs/ bytes are what the circulating PDF embeds; regenerating before submission decouples the reviewed PDF from the tree. The MANIFEST documents the orphaned provenance honestly ("manual crop, FORMAT_AUDIT F7"); regeneration becomes an explicit, audit-logged post-circulation option (Phase C6).

**3. Comment-only banners on committed superseded macros (E7) vs strictly-additive-only (E7's own requested-expert question).** RULING: banners allowed. The house precedent exists — mixedHeadToHead.C's banner was itself a post-hoc comment edit and is universally praised by the panel; git preserves the as-run bytes; and the live clobber risk of paperSystematics.C/methodCompare.C (both still write papers/timing/tab_*.tex) justifies in-file warnings. INDEX.md is added regardless. lib/ header comment fixes are deferred to post-tag Phase B with "documentation-only" commit labels, honoring E7's ACLiC caution.

**4. Rename the stale narrative systematics.png (E7) vs never-move-cited-paths (E1, E8).** My spot-check found zero references to it in any tracked md/html/py, so the rename is safe in principle. RULING: breadcrumb-first (additive README beside it) in Phase A; the git mv rename is optional Phase C — conservative because only tracked files were checked.

**5. chain_of_evidence images: vendor existing bytes (E1, E5) vs regenerate via waveformProfiles.C (E6 option).** RULING: vendor-copy only. Nothing under the trust narrative gets regenerated, even non-gated QA images — copies of the existing bytes are sufficient and provably identical.

**6. MIXED.json flag: amend now to cite the gate (E3) vs keep until the logbook closes (E1, E2).** RULING: split the difference — the annotation field may be amended to cite the GATE-1 8/8 pulse-shape confirmation (so the config stops contradicting the manuscript's "pulse-shape-confirmed" claim) while explicitly retaining "logbook half pending"; full resolution only by dated append after the logbook scan. Never delete; channel-map values untouched.

**7. Paper numbering.** RULING: adopt E4's series scheme (parent = I, timing = II, energy/position = III), declared once in papers/README.md; memory filenames (memory_paper2_timing.md etc.) are never renamed; folder labels corrected at the next manuscript editing pass.

**8. review_history/ move (E4) vs move-nothing conservatism (implicit in E1/E5/E8).** RULING: demoted to optional Phase C. It is the riskiest mechanical move proposed (9+ files cross-reference those paths) for purely cosmetic gain; if done, git mv only, recorded as a snapshot UPDATE with an old→new path table.

**9. Push sequencing.** E4/E5 want the 31 commits pushed before/at circulation; E5 flags that Pages auto-deploys the stale site on push, and my spot-check confirmed the remote is public GitHub with zero tags. RULING: push in Phase A, but only after the site banner fix AND the transparency-perimeter decision, in one push, immediately tagged.

**10. Grade spread (B- from the code-adjacent experts vs C-range from the stranger-facing ones).** Not a factual conflict — the same verified facts weighted differently. Since the user's stated goal is a stranger-trustable published workspace, the stranger-facing weighting governs the overall grade: C+ today, B- after Phase A, A-range after Phase B.

## 7. Additional expertise warranted

Three follow-ups are genuinely warranted (five of the eight requests collapse into the first):

1. **Open-data / publishing-rights specialist (HEP: Zenodo, CERN Open Data, HEPData, Elsevier reuse)** — consolidates the requests of E1, E4, E5, E6, and E8 into two exact questions: (a) What is the right DOI vehicle and license for the ~13 GB reduced 'rad' ntuples + configs/.hglg sidecars + the beam logbook.csv so the paper can cite an immutable dataset, and does NIM A's data-availability expectation prefer Zenodo, CERN Open Data, or a Zenodo+HEPData split (sigma_t(E) tables to HEPData)? (b) Given the authors' retained rights for NIM A 1068 (2024) 169737 and its arXiv posting, may the 88-file papers/timing/NIM_A_Figures/ inventory — and the two panels adapted into the new Fig. 1 — be redistributed in a public GitHub repo, and under what license/attribution text?

2. **Detector-hardware coauthor / parent-paper hardware lead (E2's request)** — settle claims-law open question Q8: are HG/LG gain-split copies of eight SiPMs (the manuscript's "eight SiPMs, sixteen readout channels") or separate shower-max vs full-length WLS optical paths (ChannelConfig.h:75-76 and README framing)? Every apparatus-doc rewrite in Phase B depends on the answer.

3. **Not an expert but a required user action (E3's request)** — the DAQ operator/logbook holder (the user): scan the physical logbook page with the MIXED module assembly to close GATE 1's open logbook half, and check whether the DT5742 run-control settings files still exist on Argon/EOS for the apparatus compendium.

E7's requested FAIR/RO-Crate question (are post-hoc comment banners on committed macros acceptable?) is resolved by ruling 3 in Conflicts — the house precedent (mixedHeadToHead.C) plus git-preserved as-run bytes settle it; no further expert needed.

## 8. Final recommendation

The honest distance: this workspace is roughly two weeks of purely additive, non-scientific work away from world-class — the interior record (pre-registered gates, claims law, committed evidence, verified 25.7±0.6 ps chain) is already at analysis-preservation quality exceeding normal HEP practice, but today a stranger meets five consecutive front doors that state retired physics, cannot obtain a single event of data, finds the designated trust artifact broken on clone, and discovers that the public repo does not even contain the gated era (31 commits unpushed, zero tags, no license). The single highest-leverage next move is the Phase-A front-door package executed as one push before coauthor circulation: rewrite the root README to the claims-law headline with a complete repo map, add ANALYSIS_GUIDE.md as the waterfall spine, vendor the chain-of-evidence images, import the orphaned external-memory dataset knowledge into the repo, close the full_fiducial_check evidence gap with an honest post-hoc audit + log, banner the stale site, and push + tag `circulation-2026-06` — because coauthors are the first strangers, and every one of those steps is mechanical, reversible, and touches no gated byte.

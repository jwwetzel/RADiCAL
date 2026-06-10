# Prior-Record Audit — RADiCAL timing paper (2026-06-09)

Continuity audit of the Paper-2 draft (commit `cb19f51`) against the prior-public RADiCAL record:
parent NIM A 1068 (2024) 169737; Instruments 6(3) 27 (2022); IEEE TNS 70(7) 1296 (2023);
CPAD/SLAC 2023 talk (38 slides); FCC/MIT 2024 talk (49 slides). Six independent auditors + a
coordinating-editor synthesis (multi-agent workflow). Companion files:
`PRIOR_RECORD_CONTINUITY_MATRIX_2026-06-09.md`, `PRIOR_RECORD_ACTION_PLAN_2026-06-09.md`.

## Executive summaries

### Agent 1 — Prior-publication continuity

Continuity audit of the Paper-2 timing draft (commit cb19f51) against the parent NIM A 1068 (2024) 169737, Instruments 6(3) 27 (2022), IEEE TNS 70(7) 1296 (2023), and the CPAD-2023/FCC-MIT-2024 talks. The draft is cleanly differentiated: its core deliverables (four-build scan, LuAG/TENERGY/MIXED timing, srCFD recovery, same-shower control, systematics budget, depth drift) are all new, and it directly fulfills the parent's Sec. 7 item 2 promise (the DSB1-vs-LuAG:Ce comparison), which the intro already acknowledges. Duplicate-publication risk is LOW-MODERATE and confined to the DSB1 energy scan, which reuses the parent's own May-2023 dataset — a reuse the parent itself pre-registered ("compared with the data presented here"); one explicit same-dataset sentence would close it. Three substantive continuity defects: (1) the Discussion cites "Ref. [radical] (its Section 7.3)" but the parent has no §7.3 — it is item 3 of a numbered list (MUST fix); (2) the draft attributes "cfd05" to the parent, while the parent describes a per-channel fixed threshold on the HG leading edge — needs a reconciling footnote or relabel (SHOULD, verify internally); (3) the ~0.2 ns satellite is cited as "also noted in Ref. [radical]", but the parent observed it only in the MCP-referenced SUM and attributed it to a digitizer sampling artifact, explicitly absent from the difference — the draft sees a clipping-induced shoulder in the DIFFERENCE, so the conflation needs one qualifying clause (SHOULD). Also flagged: the headline 25.7±0.6 ps appears without the full-fiducial companion number required by the claims law (MUST per the law); missing lineage citations (TNS 28-GeV precursor, parent Fig. 8 light-yield simulation, BestMinus identification, T-type filament architecture). The talks contain only preliminary parent-era DSB1 results plus claims that violate the current claims law ("radiation hard demonstrated", "<10 ps potential") — obsolete content the draft correctly does not propagate.

### Agent 2 — Presentation history

Audited both public talks (CPAD/SLAC Nov 2023, 38 slides; FCC-US/MIT Mar 2024, 49 slides) against the current manuscript (papers/timing/radical_timing.tex), the coauthor circulation note, and the claims law. Findings: (1) The talks' headline numbers (17.52 ps floor / 255.58 ps·sqrt(GeV) preliminary fit, "25 ps @ 150 GeV, ~18 ps limiting") are the cfd05-chain ancestors of the published 17.5/256/27 and are CONFIRMED, not contradicted, by the new srCFD numbers (203±6, 18.8±0.8, 25.7±0.6) — but the a: 256→203 change needs one explicit sentence in the coauthor note or coauthors will read it as a revision. (2) Four talk claims are now forbidden or deliberately absent from the paper: "RADiCAL module is radiation hard / demonstrated to meet high-radiation needs," "<10 ps potential at >150 GeV," "meets needs of FCC EndCap," and "measures position with high precision" — the note should pre-empt their reinstatement. (3) Crucially, NO MIXED/LuAG comparison and NO srCFD method ever appeared in either public talk, so the retired 0.99 ratio lives only in internal drafts; the same-shower control and the saturation-recovery method are brand-new to talk-only audiences and need framing as such. (4) The manuscript has no module/setup schematic figure — the talks' standard module diagram (CPAD slide 13 / FCC slide 12) and the timing-vs-energy capillary scan are the obvious candidates to adapt. Concrete note/manuscript additions with exact text are provided.

### Agent 3 — Citations & prior art

Citation audit of radical_timing.tex (cb19f51) against the parent NIM A paper, Instruments 6:27, IEEE TNS 70:1296, and the two talks. All four published RADiCAL outputs are cited — no missing RADiCAL publications. Two must-fix items: (1) DSB1, the paper's central variable, has zero provenance citations (parent cites Eljen + OSTI waveshifter report; LuAG side has two refs) and (2) bibitem ruchti1996 is an orphan — present in the bibliography, never cited. One must-verify benchmark: the HGTD "35–50 ps per hit" (text + tab_benchmarks) looks like the TDR per-TRACK range; per-hit is plausibly 35–70 ps. Useful adds: MCP-PMT reference (R3809U-50 / Bortfeldt NIM A 960), cates2018 at the front-end amplifier description (its role in parent/TNS), FCC-ee CDR in the intro collider list, naming HDR2 SiPM, and the parent's arXiv:2401.01747. Bibliographic inconsistencies: the draft appends "(RADiCAL Collaboration)" to the parent citation though the published byline has no collaboration designation; internal comments call the parent "Ref. [2]" but it is [1]; cantone2023's title "Crilin beam-test results" reads as an unverified paraphrase. Intro anchoring verified strong: the parent's §7 item 2 is precisely the promised DSB1-vs-LuAG comparison, and the ≈0.2 ns satellite attribution and a=256/b=17.5/27 ps numbers all check against the parent PDF.

### Agent 4 — Technical consistency

Audited radical_timing.tex (commit cb19f51) against all five prior-record sources plus the authoritative 2023 channel map. Three HIGH-severity active defects: (1) the manuscript says HG and LG copies are "both digitised … at 5 GS/s," but the author-confirmed channel map puts all eight LG channels on DRS1 at 1 GS/s — method-relevant, since srCFD's amplitude anchor comes from the LG waveform; (2) the published estimator of Ref. [1] is attributed as "cfd05 … as in Ref. [radical]" (text + Table 1), but the parent paper states fixed thresholds on the HG leading edge (its Sec. 5.2 / Fig. 12) — the parent's lead author will flag this; (3) the Discussion's "16 SiPMs in a 14×14 mm² cell" contradicts both the manuscript's own "4×2 = 8 SiPM channels" and the parent's 2×4 HDR2 cards (8 SiPMs, 16 gain channels) — and the claims-law "16-SiPM module" line carries the same conflation. Medium findings: the ~0.2 ns satellite is cross-referenced to Ref. [1], which explicitly says its satellite does NOT appear in (DW−UP)/2; TENERGY is described as "multi-segment with an added longitudinal energy-readout section" instead of 3 T-type + 1 full-length E-type capillary; the Pb-glass containment veto and trigger counters are absent from the setup; HDR2/R3809U-50 models unnamed; the LuAG "factor ~2.2 light deficit" conflicts with the measured ~3× ratio. Prior-record-internal X₀/length disagreements (4.5–4.7 mm/114 mm vs 5.4 mm/135 mm) are historical-only; the current draft implicitly uses the parent's X₀ = 5.4 mm and should say so. Talk-only claims (rad-hard demonstrated, <10 ps potential) are obsolete and correctly not propagated.

### Agent 5 — Figures & tables

Audited all 11 figures/3 tables of the current manuscript against the five prior records. (1) Verified: Section 2 has NO apparatus figure — every prior RADiCAL publication and talk leads with the module schematic (parent Fig. 1) and a capillary-placement cross-section (parent Fig. 2 / Instruments Fig. 2); the MIXED corner map (NE/SW vs NW/SE) is described only in words, which weakens the paper's unique asset. Recommend adapting parent Fig. 1 plus a small beam's-eye corner-map schematic, and citing (not reusing) parent Figs. 7, 8, 11, 12. (2) Retirement list: the "Preliminary"-watermarked sigma_t(E) plot from both talks, the FCC-2024 "25 ps @ 150 GeV / ~18 ps limit" bullet, all talk radiation-hardness/"<10 ps potential"/"meets FCC EndCap" claims, the 29.8 mV/GeV preliminary linearity, and the stale pre-fix thesis.png/method.png still sitting in figs/. (3) The money plot diverges from the established visual story (parent Fig. 27 right has per-point error bars + single fit curve); the existing TODO-P2 restyling closes this, and a grey published-fit reference curve would make "confirms 17.5 ps" visible. (4) One real continuity defect: the manuscript twice attributes "cfd05, as in Ref. [1]" while the parent paper's own text describes an optimized fixed-threshold leading-edge time — exact hedged caption wording supplied. Provenance: keep "post-fix/corrected" labels OUT of public captions (no published antecedent exists); continuity is correctly carried by the Table 1 published row and "confirms" verbs; add "deliberately driven into saturation (Ref. [1] Sec. 4)" to defuse the saturated-digitizer objection.

### Agent 6 — Referee risk

Referee-risk audit vs the full prior record (parent NIM A 169737, Instruments 2022, TNS 2023, CPAD 2023 + FCC/MIT 2024 talks) finds the existing 10-risk memo solid but incomplete. Three NEW high-priority exposures: (1) the manuscript never states explicitly that the DSB1 energy scan re-analyzes the SAME dataset as the parent paper — the strongest "already published" attack surface; (2) the manuscript attributes the parent's timing method to "cfd05 (a CFD on the clipped measured peak, as in Ref. [1])", but the parent's published text describes a per-channel FIXED THRESHOLD on the leading edge — a method-attribution contradiction the parent's corresponding author will catch; (3) the manuscript violates its own claims law (forbidden #9): the full-fiducial companion number (27–29 ps) never appears beside the brightest-1000 headline 25.7 ps. Two MED findings: the ~0.2 ns satellite cross-reference contradicts the parent (parent says it appears ONLY in the MCP-referenced sum and is a digitizer sampling effect, absent from (DW−UP)/2; the manuscript reports it in the differential and blames clipping), and the absence of simulations is under-defended when the programme's own G4 (parent Fig. 8) predicted exactly the light-yield scaling the paper measures — a free supporting citation. Talk continuity is good (FCC 2024 already previewed 25 ps / ~18 ps); the talks' "<10 ps potential" and "radiation hard demonstrated" claims need only stock responses, not edits. Old-vs-new numbers (27/17.5 vs 25.7/18.8) and the retired MIXED 0.99 are well defended. Six exact text edits recommended; specific update bullets for REFEREE_RISK_MEMO and the coauthor note provided, including one new coauthor question (parent-chain estimator verification, for Perez-Lara).


---

# Full auditor reports



---

# Agent 1 — Prior-publication continuity

# Prior-Publication Continuity Audit — RADiCAL Timing Paper (commit cb19f51)

Sources read in full: [PARENT] NIM A 1068 (2024) 169737 (all 15 pp.); [INSTR] Instruments 6(3) 27 (2022) (6 pp.); [TNS] IEEE TNS 70(7) 1296–1300 (2023) (5 pp.); manuscript `papers/timing/radical_timing.tex`; circulation note; referee-risk memo; claims law. Talks skimmed for superseded-claim content: [CPAD] SLAC Nov 2023 (38 slides), [FCC] FCC/MIT Mar 2024 (49 slides).

---

## 1. Novelty map — result by result

| Draft result | Prior record | Verdict |
|---|---|---|
| Module concept: W/LYSO:Ce shashlik, 14×14×135 mm³, 25 X₀, WLS capillaries, dual-gain HG/LG, MCP reference | Fully published: INSTR §2, TNS §II, PARENT §2/§4 | Background — correctly treated as such; two missing first-use citations (HG/LG split, T-type filament; see §3 items c, d) |
| DSB1 σ_t(E), 25–150 GeV, (DW−UP)/2, a=256, b=17.5, 27 ps @150 | **Published** (PARENT Eq. 2, Figs. 20–27). Same module (M2925W2815L-415DSB1-T), same May 2023 H2 dataset | **Re-analysis, not new data.** Draft handles it correctly in spirit (published row in Table 1, floor verb "confirms", new numbers labelled srCFD/brightest-1000) but never states outright that it is the *same dataset* — see §4 and edit S2 |
| σ_t improves with shower brightness/amplitude | Published: TNS Fig. 10 (28 GeV), PARENT Figs. 21–26 (amplitude-binned σ_t) | Known effect, re-presented (Fig. dist bottom). Draft says "as expected from the slew argument" — should add the lineage citation (edit O7) |
| HG clips at 820 mV; LG-predicted-peak recovery (**srCFD**); ×3.7 edge steepening; identical-event method gain; tail decomposition; per-regime source rule | Nowhere. PARENT *deliberately* drove HG into saturation and timed with a fixed threshold; no recovery method exists in any prior document | **New** — and it is the methodological heart of the paper |
| LuAG:Ce build timing (a = 440±18; 44.4 ps @150) | Never measured. PARENT Sec. 7 item 2 promised exactly this comparison | **New** (promise delivered) |
| TENERGY multi-segment build timing | Not in any prior publication | **New** |
| MIXED same-shower control, σ_DSB1/σ_LuAG = 1.04±0.05, estimator spread 0.80–1.16 | Nowhere; PARENT envisioned only a *cross-module* comparison | **New, and exceeds the promise** — the same-shower cancellation design is the paper's unique asset |
| Four-build floors 19–26 ps, mutually consistent with caveats; DSB1 floor 18.8±0.8 confirms 17.5 | 17.5 ps published (PARENT) | New measurement framed correctly as **confirmation** (claims-law compliant) |
| Depth drift −33.6±2.9 ps per e-fold (acceptance-conditional consistency reading) | No precedent in any prior document or talk | **New** |
| Selection-systematics budget (Table 2), OOS-style optimization scans | PARENT quotes no selection systematics | **New** |
| Benchmark table (SPACAL, Crilin, BTL, HGTD) | Not in prior RADiCAL papers | New contextualization |
| ~0.2 ns satellite | PARENT noted a 0.2 ns satellite — but in the **MCP-referenced sum**, ~15% of signals, attributed to a digitizer one-sample-segment shift, and **explicitly absent from (DW−UP)/2** ("It does not appear in Fig. 19 Lower Left") | Draft's "(the ∼0.2 ns satellite also noted in Ref. [radical])" attached to the **difference** distribution conflates two distinct observations/mechanisms — see edit S4 |
| Estimator label for the published chain: "cfd05 … as in Ref. [radical]" | PARENT §5.2 describes "fixed thresholds … set for the high-gain signals … on the leading edge" — a fixed-threshold LED, not a constant fraction of the clipped peak | **Attribution tension** — a referee reading both will catch it; see edit S3 |

Talks check: CPAD-2023 and FCC-MIT-2024 contain only the preliminary parent-era DSB1 plots (17.52/255.58 fit, "Preliminary" watermark) and **no** LuAG/MIXED/TENERGY numbers, no srCFD, no depth dial. They also contain claims now forbidden ("The RADiCAL module is radiation hard" as a demonstrated result; "potential to reach <10 ps"; "25 ps @ 150 GeV with limiting resolution ~18 ps"). The draft propagates none of these — correct. Nothing in the talks pre-publishes any draft result.

---

## 2. "Promises kept" map — parent Sec. 7 (and earlier outlooks) vs this draft

The parent's Sec. 7 is a **numbered list (items 1–4), not subsections** — relevant to edit M1.

| Promise | Source | Status in draft |
|---|---|---|
| **Item 1**: EM shower localization study, 25–150 GeV, vs wire chamber (also §5.1: "will be discussed in a future paper") | PARENT | Correctly **deferred** to the companion energy/position paper (§7 and Conclusions). Clean division of labor; no leakage of position claims |
| **Item 2**: "Comparison of DSB1 and LuAG:Ce Wavelength Shifters … compared **with the data presented here** … in 25 GeV steps over 25 ≤ E ≤ 150 GeV" | PARENT | **Delivered — this is the paper.** The intro already says so ("its Section 7 … This paper presents that comparison"). Two refinements: (a) point to item 2 specifically; (b) the realized scan is 50–150 GeV for the non-DSB1 builds (25 GeV exists for DSB1 only) — one clause acknowledging the deviation from the promised 25 GeV steps would preempt a pedantic referee (edit O6). Note the phrase "with the data presented here" **pre-registers the DSB1 dataset reuse** — quote-worthy (see §4) |
| **Item 3**: Enhanced modular design — double the LYSO:Ce thickness at shower max → ×2 light → σ_t improved ~1/√2; FCC-hh energy goals; 3×3 array | PARENT | Not a beam result here (correct); used as the basis of the LuAG light-recovery **engineering extrapolation** in the Discussion, cited as "its Section 7.3" — **broken pointer** (no §7.3 exists); fix to "Section 7, item 3" (edit M1) |
| **Item 4**: New 18×18 mm² module, 5 T-type + 4 E-type capillaries, LuAG:Ce E-type, T-type DSB1-vs-LuAG, 25–200 GeV, then 3×3 array | PARENT | Future work, untouched — correct. If TENERGY is in fact the first step toward the simultaneous timing+energy readout of items 3–4, one clause in §2.1 saying so would strengthen the continuity narrative (edit O9; only if factually right) |
| "Future plans include measurements of the timing performance over an extended range 25 < E < 200 GeV" | TNS Conclusion | Delivered to 150 GeV (parent: one build; this draft: four builds). No statement needed |
| GEANT4 expectation: σ_t set by detected light level (Fig. 8 / Fig. 5b); 30–50 ps achievable; "By reading out … both upstream and downstream ends and using several timing capillaries … the aim has been to reach or exceed this expectation" | PARENT §3, INSTR §3.2 | The draft's central thesis (stochastic term tracks detected light) is the **measured confirmation of this published simulation expectation — and the draft never cites it.** Adding the citation both strengthens continuity and defuses the "your conclusion is obvious" objection (referee-risk #1) by reframing it as the quantitative verification (edit S5) |
| <50 ps timing verification | INSTR objective | Already delivered by TNS (42–49.5 ps @ ~28 GeV) and PARENT; draft's intro skips the TNS step entirely (edit S6) |

**Answer to the explicit question:** Yes — the draft *should* say it follows the parent's Sec. 7, and it **already does** (Intro: "That work also identified, as a priority for further study (its Section 7), a direct comparison… This paper presents that comparison."). This is the right framing; tighten it to the item level and consider quoting the promise verbatim (edits M1, O6).

---

## 3. Places the draft should cite/quote prior RADiCAL work

a) **§1, second paragraph (must-adjacent accuracy):** "A first RADiCAL test-beam measurement…" skips the actual first timing beam test (TNS, FTBF 2022). Suggested replacement: *"An initial timing beam test at the Fermilab Test Beam Facility demonstrated σ_t ≈ 42–50 ps for E ≈ 28 GeV electrons~\cite{wetzel2023}. A subsequent measurement at the CERN SPS, on a single high-light build, reported a time resolution of ≈27 ps at 150 GeV…"*

b) **§1 same paragraph:** point to the parent's item 2 specifically: *"…identified, as a priority for further study (its Section~7, item~2), a direct comparison of the timing performance of different WLS materials in otherwise-identical modules."*

c) **§2.1 (T-type architecture; also needed because §7 invokes "the fixed-position wavelength-shifter filament" without prior introduction):** after "read out longitudinally by WLS capillaries", add: *"Each capillary is of the `T-type' of Refs.~\cite{instruments2022,radical}: a short WLS filament positioned at the depth of shower maximum, fused between solid quartz waveguides that pipe the shifted light to SiPMs at both ends."*

d) **§2.2 (dual-gain readout):** *"…both digitised by DRS4 … as in Refs.~\cite{radical,wetzel2023}."*

e) **§4 (estimator lineage — defuses "method shopping"):** after Eq. (2): *"Equation~(\ref{eq:dwup}) is the `BestMinus' estimator identified as optimal in Ref.~\cite{radical}; the brightest-$K$ selection replaces the central energy-bin selection used there with a deterministic, fixed-size criterion."*

f) **§4 (brightness dependence):** *"Within each build, the resolution improves with shower brightness, as expected from the slew argument and as reported in Refs.~\cite{wetzel2023,radical}."*

g) **§5.2 interpretation paragraph (thesis anchored to the published simulation):** *"This is the behaviour anticipated by the GEANT4 study of Ref.~\cite{radical} (its Fig.~8), which predicted the timing resolution of this geometry to be governed by the detected light level; the four-build scan is its quantitative, measured confirmation."*

h) **§2.2 or Table 1 caption (same-dataset transparency — see §4):** *"The DSB1 exposure analysed here is the dataset of Ref.~\cite{radical}, whose Section~7 explicitly foresaw its reuse for this comparison; it is re-timed under the uniform four-build chain of Sec.~\ref{sec:method}, and the published values are quoted unchanged for comparison in Table~\ref{tab:builds}."*

i) **§4 (satellite reconciliation):** see edit S4 below.

---

## 4. Duplicate-publication risk assessment

**Overall: LOW–MODERATE, and reducible to LOW with one sentence.**

- **The only overlap with the published record is the DSB1 energy scan** — same physical module, same May 2023 H2 exposure as the parent. Mitigations already in place: the published (a, b, σ150) appear as a clearly-labelled comparison row in Table 1; the floor verb is locked to "confirms"; the new headline (25.7±0.6) is tied to a *different estimator and selection* and never presented as a revision; the paper's abstract, title, and conclusions center on content absent from every prior publication (four-build comparison, same-shower control, srCFD, depth drift). What is **missing** is the explicit statement that the DSB1 data are the *same dataset* re-analysed. Without it, a careful referee may read Fig. thesis/Table 1 as implying an independent DSB1 measurement. The parent itself pre-registered the reuse ("compared with the data presented here") — say so (edit S2). With that sentence the residual risk is negligible.
- **No figure is reproduced from any prior publication** (all eleven figures are new analysis products) — no copyright/permission exposure.
- **No text recycling observed** — the setup section is a fresh condensation, not copied prose.
- **The estimator-label tension (S3) is the one place the draft could be accused of misdescribing prior work**: it attributes "cfd05 (a CFD on the clipped measured peak)" to the parent, while the parent describes a per-channel fixed threshold on the HG leading edge. This is an accuracy-about-prior-work issue, not duplication, but it sits exactly where a continuity-minded referee will look (the published-row comparison).
- **Conference talks pose no duplication issue** (slides, preliminary watermarks, parent-era content only). Conversely, none of their stronger claims may be imported: "radiation hard demonstrated", "<10 ps potential", "meets needs of FCC EndCap", and the talk-rounded "25 ps @150 / limiting ~18 ps" all violate the claims law or the floor-verb rule. The draft currently imports none — keep it that way.

---

## 5. Recommended manuscript edits

| # | Location | Severity | Edit (exact suggested text) |
|---|---|---|---|
| M1 | §7 Discussion, "as already foreseen in Ref.~\cite{radical} (its Section~7.3)" | **MUST** | The parent has no subsection 7.3; Sec. 7 is a numbered list. Replace with: *"as already foreseen in Ref.~\cite{radical} (its Section~7, item~3)"* |
| M2 | Abstract, §5.3, Conclusions — every appearance of the 25.7±0.6 ps headline | **MUST** (per the claims law's own non-negotiable rule: "full-fiducial always quoted beside any selection") | The full-fiducial companion number appears nowhere in the manuscript. At first headline use in §5.3 and in the Conclusions add: *"($27$--$29$~ps over the full fiducial sample)"*. Reconcile with the 2026-06-09 reframe audit, which presumably should have caught this. This also *strengthens* continuity: the full-sample number visibly brackets the published 27 ps |
| S2 | §2.2, end of "Readout, beam, and reference" (or Table 1 caption) | SHOULD | Add the same-dataset sentence of §3(h) above. Closes the duplicate-publication gap and pre-empts "why do your DSB1 numbers differ from the published ones" |
| S3 | §5.3 "cfd05 (a CFD on the clipped measured peak, as in Ref.~\cite{radical})" and Table 1 published row "src = cfd05" | SHOULD (verify against the production chain before editing) | The parent describes a fixed threshold, not a clipped-peak CFD. Either relabel the published row ("LE thr., Ref.~[1]") or add a reconciling footnote: *"Ref.~\cite{radical} times the high-gain leading edge with a per-channel threshold referenced to the recorded (clipped) pulse; cfd05 — a constant fraction of the clipped peak — is the closest member of that discriminator family in the present chain and reproduces the published values (Table~\ref{tab:builds})."* If internal documentation shows the parent's chain was literally cfd05, keep the label but still add the reconciling clause, because the parent's *text* says "fixed threshold" |
| S4 | §4, "(the ∼0.2~ns satellite also noted in Ref.~\cite{radical})" | SHOULD | The parent observed the satellite in the MCP-referenced **sum**, attributed it to a digitiser one-sample-segment shift, and stated it does **not** appear in the difference. Replace the parenthetical with: *"(a structure at the ∼0.2~ns scale was also noted in Ref.~\cite{radical}, there in the MCP-referenced sum and attributed to a digitiser sampling artefact; the shoulder discussed here appears in the difference under clipped-peak estimators and is shown in Sec.~\ref{sec:method-gain} to be clipping-induced)"* |
| S5 | §5.2 interpretation paragraph | SHOULD | Add the parent-Fig.-8 simulation-confirmation sentence of §3(g). Directly serves referee-risk #1 |
| S6 | §1, "A first RADiCAL test-beam measurement…" | SHOULD | Insert the TNS precursor sentence of §3(a) — establishes the full lineage FTBF→parent→this paper and avoids over-claiming "first" |
| S7 | §4, after Eq. (2) | SHOULD | BestMinus identification sentence of §3(e) — ties the estimator to the published one (anti-"method-shopping") |
| S8 | §2.1 | SHOULD | T-type filament sentence of §3(c) — also repairs the forward dependency of §7's "fixed-position wavelength-shifter filament" |
| O1 | §1, "(its Section~7)" | OPTIONAL | Tighten to "(its Section~7, item~2)"; optionally quote: *"a comparison `in 25~GeV steps in electron beam energy over the range 25~GeV $\le$ E $\le$ 150~GeV'"* |
| O2 | §2.2, "the other builds at 50–150 GeV" | OPTIONAL | Append: *"(the 25~GeV point of the comparison foreseen in Ref.~\cite{radical} was recorded for DSB1 only)"* |
| O3 | §4 brightness sentence | OPTIONAL | Citation addition of §3(f) |
| O4 | §2.2 dual-gain sentence | OPTIONAL | Citation addition of §3(d) |
| O5 | Header comments, lines 3–4 and 25–26 ("Ref.~[2], NIM A 1068 (2024) 169737") | OPTIONAL | The parent is Ref. [1] in the current bibliography order; fix the comment to "[1]" so the roster-completion step doesn't grab the wrong reference |
| O6 | Keywords | OPTIONAL | "4D calorimetry" is Paper-3 territory; consider "shower-maximum timing" instead. Low risk either way (a keyword is not a claim) |
| O7 | §2.1 TENERGY description | OPTIONAL (only if factually correct) | If TENERGY realizes the simultaneous timing+energy direction of parent Sec. 7 items 3–4, add: *"a first step toward the combined timing--energy readout foreseen in Ref.~\cite{radical} (its Section~7, items~3--4)"* |

**Compliance note on these recommendations:** none revives the retired 0.99/χ²=0.4 ratio, uses "statistically indistinguishable"/"consistent with unity", claims universal WLS irrelevance, claims achieved 5D, calls the depth dial a calibrated z, or calls the position result an intrinsic resolution. Table 1, Table 2, the corrected MIXED figure, and the claims law remain authoritative over all talk-era numbers.


---

# Agent 2 — Presentation history

# Presentation-History Audit — RADiCAL timing paper vs CPAD 2023 & FCC-MIT 2024 talks

Sources compared: [CPAD] `/Users/jameswetzel/Downloads/CPAD_SLAC_2023_WETZEL_RADiCAL.pdf` (38 slides, Nov 7 2023); [FCC] `/Users/jameswetzel/Downloads/FCC_MIT_2024 Wide.pdf` (49 slides, Mar 26 2024); manuscript `/Users/jameswetzel/Dropbox/Research/CERN_TEST_BEAMS/CERN_May_2023/RADiCAL/RADiCAL/papers/timing/radical_timing.tex` (commit cb19f51); `COAUTHOR_CIRCULATION_NOTE_2026-06-09.md`; `REFEREE_RISK_MEMO_2026-06-09.md`; `papers/memory_claims_and_forbidden_language.md`.

Headline structural finding first, because it conditions everything else: **neither public talk contains any MIXED-module result, any LuAG timing number, any per-build comparison, or any mention of srCFD/HG-LG saturation recovery.** The talks present a single-build (DSB1) story with the preliminary cfd05-class "Best Estimator" chain. Therefore (a) the retired "0.99 / chi2=0.4" MIXED ratio was never public — its only audience is internal-draft readers, already addressed by circulation-note item 1; and (b) the paper's two central assets (the same-shower control and the srCFD estimator taxonomy) are *entirely new* to anyone whose memory is anchored on the talks. The circulation note currently frames "what changed" only against *older internal drafts*; it should also frame what is new/changed against the *talks*, which is what most coauthors actually saw.

## 1. Old claims now superseded

| # | Old claim | Where it appeared | What the current record says |
|---|---|---|---|
| 1 | sigma_t(E) fit "17.52 (const), 255.58 (sqrtE)", labeled Preliminary, "Best Estimator"; 150 GeV point reads ~26.4 ps | CPAD slide 30; FCC slide 38 (same plot, now tagged arXiv:2401.01747) | This is the ancestor of the **published** cfd05 values (a=256, b=17.5, sigma=27 ps, Ref. [1]), which the paper keeps as the comparison row in Table 1 and never revises. The new primary chain (srCFD, brightest-1000) gives a=203±6, b=18.8±0.8, 25.7±0.6 ps. The a: 256→203 difference is estimator + selection, not a re-measurement — this is stated in Table 1's caption but not in the coauthor note (see §4-A). |
| 2 | "Measured timing resolution of 25 ps @ 150 GeV, with limiting resolution of ~18 ps" | FCC slide 39 | Superseded by the precise post-fix statements: 25.7±0.6 ps (150 GeV, brightest-1000, srCFD) with 27–29 ps full-fiducial, floor 18.8±0.8 ps **confirming** the published 17.5. Continuity is good; only the unqualified rounding is obsolete. Claims law requires the selection to be named and full-fiducial quoted alongside. |
| 3 | "The RADiCAL module is radiation hard." / "The RADiCAL concept has been demonstrated to meet the needs of high radiation … environments" | CPAD slides 31, 33; FCC slides 40, 43 | **Forbidden in the paper** (claims law item 8: zero irradiation data in this dataset; the best-timing build is organic DSB1). The manuscript makes only component-level hardness statements by citation (Ref. [luag]) plus the LuAG timing-viability + engineering-extrapolation framing. The FCC backup slides 46–47 (Caltech RIAC plots, Co-60 capillary study) are exactly the component-level evidence the paper points to by citation — they are not module-level beam demonstrations. |
| 4 | "RADiCAL has potential to reach <10 ps timing resolution at >150 GeV" + "GEANT4 simulations predict time resolution of ~10 ps" | CPAD slides 17–18, 31; FCC slides 18–19, 40 | Implicitly superseded by the paper's floor physics: the 19–26 ps floors are read as longitudinal shower fluctuation (depth drift −33.6±2.9 ps/e-fold consistency reading), so raw light/energy scaling does **not** reach 10 ps; the stated route below the floor is a depth-corrected estimator (Discussion, §8). The 10 ps figure survives only as the long-term *program target* (memory: detector_invariant), not a claim of this dataset. |
| 5 | "<30 ps for 150 GeV electrons — meets needs of FCC EndCap" | CPAD slide 31; FCC slide 40 | The paper deliberately replaces collider-specific "meets needs" language with the comparator bands: CMS BTL 30–60 ps per track, HGTD 35–50 ps per hit, explicitly labeled per-MIP/track comparators, plus the per-shower benchmark table (SPACAL, Crilin) "without ranking them". Note the requirement-band overlay on the money plot is still an open TODO-P2 at radical_timing.tex line ~292. |
| 6 | "The module is potentially capable of measuring shower energy, shower time, and shower position with high precision" / "Shower position can be determined" / measured sigma_E/E fit 9.31% ⊕ 52.04%/sqrtE ⊕ 31.62%/E | FCC slides 16, 34–36, 41 | Energy/position/4D content is reserved for the companion paper (Paper 3) under its own gated language (E_SM tag, 1.5/0.9 mm as upper bounds, never "high precision position"). The timing manuscript correctly contains none of it; coauthors primed by the FCC talk may expect it (see §3). |
| 7 | Fermilab "45 ps @ 28 GeV" summary | CPAD slides 25–26; FCC slides 29–30, 39 | Not superseded — that is the published TNS result (Ref. [4]) from a different campaign. The paper cites Ref. [4] without re-quoting the number; no conflict. |
| 8 | Module WLS list "DSB1 plastic, LuAG:Ce, and QD glass or polysiloxane" | CPAD slide 13; FCC slide 12 | The 2023 campaign and the paper cover DSB1 and LuAG:Ce only (plus MIXED, TENERGY). QD/polysiloxane are program-level future options; no QD build exists in this dataset. |

Also checked: nothing in either talk asserts the old MIXED story, "statistically indistinguishable," "consistent with unity," a calibrated z coordinate, achieved 5D, or a revision of 17.5/27 ps — the talks are clean of the specifically forbidden formulations, mostly because they predate the analyses entirely.

## 2. Useful graphics/framing from the talks worth adapting

1. **Module schematic (CPAD slide 13 / FCC slide 12).** The W(2.5 mm)/LYSO(1.5 mm)/quartz-capillary exploded diagram with 14 mm × 114 mm dimensions and the WLS callout. **The manuscript currently has no setup or module figure at all** — §2 is text-only, and the first figure a reader meets is the clip histogram. Referees of an instrumentation paper will expect a geometry figure; this one (originating in the parent paper, so provenance is easy) is the natural Fig. 1. Adding the DW/UP end labels and the four-capillary corner naming (NE/SW/NW/SE) would also pre-empt the §6 corner-map geometry questions (referee risk #7).
2. **Timing-vs-energy capillary scan (K. Ford; CPAD slide 12 / FCC slide 11).** The single plot that explains *why* a short WLS segment at shower max is a timing device and a full-length filament is an energy device. Strong candidate for the setup figure's panel (b), or at minimum for the companion paper.
3. **"Shower max radius ~4.5 mm vs Molière 13.7 mm" framing (CPAD slides 14–16 / FCC slides 15–17).** The compactness argument in one number pair. The Introduction's "where the signal is largest and most localised" could carry the 4.5 mm vs 13.7 mm contrast in one parenthesis (values from the parent paper's GEANT4 study, A. Ledovskoy).
4. **Zhu fast-scintillator table (FCC slide 13).** Quantitative anchor for the §6 kinetics expectation: LuAG:Ce light yield 35–48% of LYSO:Ce, decay 820/50 ns vs 40 ns, light in first ns 240 vs 740 photons/MeV. One sentence with these numbers (cited to Zhu/Hu) would sharpen "this is the configuration in which kinetics would be expected to show" — it quantifies both the light deficit (~×2.2, matching the measured stochastic ratio) and the slowness. Use as text + citation, not a reproduced table.
5. **Estimator-ladder plot (FCC slide 49, backup).** Downstream / Upstream / BestPlus / BestPlus-MCP / (DW-UP)/2 / BestMinus vs energy — the visual proof that the differential estimator beats every MCP-referenced single-end combination. The paper asserts this in §4 in two sentences; the figure is the talk-era ancestor of the t_diff/t_sum taxonomy in the claims law. Optional appendix figure if a referee pushes on "(DW−UP)/2 cancels arrival time" (referee-defense pair already anticipates this).
6. **GEANT4 depth-of-shower-max profile (FCC slide 20, "Depth of shower max is slightly energy dependent," C. Perez Lara).** Direct conceptual precursor of §7's depth drift. Recommend citing the *concept* only (the paper already anchors on Longo/PDG); do **not** import the figure as supporting evidence, since the claims law records the G4 campaign as not-yet-run for floor decomposition — a talk-era illustrative sim should not be promoted to validation.
7. **The dime photo (CPAD slide 6 / FCC slide 6).** Not for the paper, but keep for the eventual seminar/press framing of "sub-30 ps from a 14×14 mm², 16-SiPM module."

## 3. Likely coauthor-confusion points

1. **"Why did a go from ~256 to 203?"** Coauthors remember 255.58/17.52 from both talks (and 256/17.5 from the parent paper). Table 1 explains it (different estimator: srCFD vs cfd05; brightest-1000 selection), but the circulation note's "what changed" list compares only against internal drafts (201→203, 25.3→25.7). A talk-anchored reader will perceive a 20% stochastic-term *revision* unless told it is an estimator/selection change with the published cfd05 row intact and unchanged.
2. **"We said the module is radiation hard — why does the paper retreat?"** The strongest talk claims (CPAD 31/33, FCC 40/43) are now forbidden. The note's limitations section says "no irradiated-module beam data in this dataset" but never connects this to the talks; a coauthor may try to reinstate talk language during review. This is the most likely wording fight.
3. **"Where did the <10 ps and FCC-EndCap claims go?"** Both talks end on them. The paper's floor-physics reading (floor = shower-depth fluctuation, common to any scintillator in this geometry) quietly *explains why* raw <10 ps is not reachable without depth correction — but nobody has said this to the coauthors explicitly.
4. **"What is srCFD, and is it method shopping?"** Talk-only readers have never seen the HG/LG recovery; their last memory is a simple best-estimator CFD chain. The note's item 3 describes the method gain but presumes context. One sentence of "this did not exist at the time of the talks; the published cfd05 numbers are reproduced unchanged as the control" defuses both confusion and the method-shopping reflex (referee risk #3).
5. **"Where are the MIXED results we never saw?"** Inverse confusion: because no talk showed MIXED/LuAG anything, §6 arrives with zero public history. Expect "when was this measured? same May 2023 run?" — worth one sentence (all four builds were exposed in the same May 2023 campaign; the comparison is new analysis, not new data).
6. **"Where is the energy resolution?"** The FCC talk title was "time **and energy** resolution" and showed the 52%/sqrtE fit (slide 36). The note's limitations mention the companion paper but a pointer tied explicitly to that talk content would help.
7. **Minor numerical memory hazards:** the ~26.4 ps point readable off the talks' 150 GeV plot (vs published 27, vs new 25.7); "45 ps @ 28 GeV" Fermilab (different campaign, still valid); "200 GeV this summer" (FCC slide 39 — whatever became of that run, it is out of scope for this paper and the note could say so in one clause).

## 4. Recommended additions (location · severity · exact suggested text)

**A. Coauthor note — new subsection after "What changed scientifically vs. older internal drafts" · HIGH**
Title: "What changed vs. the public talks (CPAD/SLAC 2023, FCC-US/MIT 2024)". Suggested text:

> If your memory of these results is anchored on the CPAD 2023 or FCC-MIT 2024 talks, four things are new or different. (1) The talks' preliminary fit (a≈255.6, b≈17.5 ps, "Best Estimator") became the published cfd05 result of Ref. [1] (a=256, b=17.5, 27 ps at 150 GeV); the paper reproduces those values *unchanged* as the control row of Table 1. The new primary numbers (a=203±6, b=18.8±0.8, 25.7±0.6 ps) come from a different estimator (srCFD) and a named selection (brightest-1000) — this is not a revision of the published result, and the floor verb remains "confirms". (2) The srCFD saturation-recovery method did not exist at the time of the talks; Sec. 5.3 and Appendix B isolate its gain on identical events. (3) The MIXED same-shower comparison (Sec. 6) was never shown publicly; it is new analysis of the same May 2023 campaign data. (4) The four-build comparison (LuAG, MIXED, TENERGY) is likewise public for the first time.

**B. Coauthor note — add to "Intentionally conservative claims (please don't strengthen)" · HIGH**
Suggested text:

> Three claims from our own talks are deliberately absent and must not be reintroduced: "the RADiCAL module is radiation hard / demonstrated to meet high-radiation needs" (zero irradiated-module beam data in this dataset; hardness is component-level by citation, per the wording law); "potential to reach <10 ps at >150 GeV" (the measured floors and their depth-fluctuation reading imply the route below ~19 ps is a depth-corrected estimator, not more light alone — Sec. 8); and "meets needs of FCC EndCap" (replaced by the labeled per-track/per-hit comparator bands and the per-shower benchmark table, which compare without ranking).

**C. Manuscript §2 (Experimental setup) — add a module/setup figure · MEDIUM**
Location: after the first paragraph of Sec. 2.1 in `/Users/jameswetzel/Dropbox/Research/CERN_TEST_BEAMS/CERN_May_2023/RADiCAL/RADiCAL/papers/timing/radical_timing.tex` (currently no geometry figure exists anywhere in the paper). Adapt the parent-paper module schematic used on CPAD slide 13 / FCC slide 12, adding DW/UP end labels and the NE/SW/NW/SE capillary corner names used in Sec. 6. Suggested caption skeleton:

> Figure 1: The RADiCAL module: alternating tungsten (2.5 mm) and LYSO:Ce (1.5 mm) plates, 14×14×114 mm³, traversed by four WLS capillaries read out at both ends (DW = downstream, UP = upstream; corners labeled NE/NW/SE/SW as used in Sec. 6). Adapted from Ref. [1].

**D. Manuscript §6, kinetics-expectation sentence · LOW–MEDIUM**
Location: the paragraph beginning "This is the configuration in which kinetics would be expected to show" (radical_timing.tex ~line 468). Suggested insertion after the first clause:

> LuAG:Ce emits roughly a third to a half of the light of LYSO:Ce, with only ≈240 photons/MeV in the first nanosecond against ≈740 for LYSO:Ce~\cite{luag}, so both the light deficit and the slower kinetics point the same way in this comparison.

(Verify the exact numbers against the Hu/Zhu source before insertion; they are read here from the R.-Y. Zhu CPAD-2019 table reproduced on FCC slide 13.)

**E. Manuscript §1 (Introduction), compactness hook · LOW**
Location: the sentence "...time the shower at its maximum, where the signal is largest and most localised" (~line 82). Suggested extension:

> (the energy-weighted shower radius at shower maximum is ≈4.5 mm, against a full Molière radius of 13.7 mm~\cite{radical})

**F. Coauthor note — "Known limitations" addition · LOW**
Suggested clause appended to the limitations paragraph:

> Energy and position results shown at FCC-MIT 2024 (sigma_E/E ≈ 52%/sqrtE at shower max; SiPM occupancy maps) belong to the companion energy/position paper and are intentionally absent here; the 2024-and-later beam periods (e.g. the planned 200 GeV run) are likewise out of scope for this 2023-dataset paper.

**G. Manuscript money plot TODO — close it before circulation · MEDIUM (already tracked)**
The TODO-P2 at radical_timing.tex ~line 292 (requirement bands 30–60 ps BTL / 35–50 ps HGTD on `thesis_postfix.png`) is precisely the defensible replacement for the talks' "meets needs of FCC EndCap" slide; closing it restores the talk's strongest rhetorical beat in claims-law-compliant form. No new analysis required (restyling of a gated product).

**H. Optional appendix figure — estimator ladder · LOW**
If referees press on the (DW−UP)/2 algebra (anticipated defense pair), the FCC slide-49 plot (Downstream/Upstream/BestPlus/BestPlus-MCP/(DW-UP)/2/BestMinus vs energy) can be regenerated from the production chain as an appendix panel demonstrating that the differential estimator beats every MCP-referenced combination at all energies. Defer unless requested.

## Compliance self-check on these recommendations
No recommendation revives the 0.99 ratio, uses "statistically indistinguishable"/"consistent with unity", claims universal WLS irrelevance, claims achieved 5D, calls the depth dial a calibrated z, or calls the position result an intrinsic resolution. Item D requires source verification before insertion (numbers currently sourced from a slide reproduction of the Zhu 2019 table). All quoted current-record numbers trace to Table 1/Table 2/claims-law authoritative values (25.7±0.6; 203±6/440±18; 18.8±0.8 confirming 17.5; 1.04±0.05; −33.6±2.9 ps/e-fold).


---

# Agent 3 — Citations & prior art

# Citation & Prior-Art Audit — RADiCAL timing manuscript (commit cb19f51)

Sources read: full `radical_timing.tex` (26 bibitems; 25 distinct keys cited) + `tab_benchmarks.tex` + `tab_systematics.tex`; complete reference lists of [PARENT] (24 refs), [INSTR] (9 refs), [TNS] (14 refs); text-extracted citations from both talks ([CPAD], [FCC]); the claims law, referee risk memo, and circulation note. Read-only; no files modified.

**Headline:** the prior RADiCAL record is fully cited — all four published RADiCAL outputs (parent NIM A 1068 169737; arXiv:2203.12806; Instruments 6(3) 27; IEEE TNS 70(7) 1296–1300) appear and are correctly keyed. The gaps are in *supporting* prior art, one orphan bibitem, and one suspect benchmark number.

---

## 1. Missing-citation list (priority-ranked)

### MUST-ADD / MUST-FIX

| # | Item | What the prior record has |
|---|------|---------------------------|
| M1 | **DSB1 provenance/properties reference — currently ZERO citations for the paper's central high-light material.** The LuAG side carries two refs (`luag`, `lucchini2013`); DSB1 carries none, yet the intro asserts "LuAG:Ce is intrinsically slower than the organic DSB1" and the whole thesis turns on the DSB1/LuAG contrast. | Parent [20]: DSB1 wavelength shifter manufactured by Eljen Technology, Sweetwater TX (URL eljentechnology.com) — also TNS [10]. Parent [21]: "Waveshifters and Scintillators for Ionizing Radiation Detection", https://www.osti.gov/servlets/purl/924745 (the Notre Dame WLS-development report; the parent's source for DSB1 λ_abs=425 nm, λ_em=495 nm, τ=3.5 ns, developed by Notre Dame with Eljen). |
| M2 | **Orphan bibitem `ruchti1996`** (R.C. Ruchti, Annu. Rev. Nucl. Part. Sci. 46 (1996) 281) is in the bibliography but never `\cite`d anywhere (verified by grep over tex + both table files). With `thebibliography` it will typeset as a numbered, never-referenced entry — a referee-visible defect. Either cite it or delete it. | Not cited in parent/INSTR/TNS either; it is the natural lineage reference for fiber/WLS optical readout (and Ruchti is a programme PI), so citing rather than deleting is defensible. |
| M3 | **HGTD benchmark number likely mislabeled** — "≈35–50 ps per hit (design)" appears in both the Discussion text and `tab_benchmarks.tex`. My strong recollection of CERN-LHCC-2020-007 is **35→70 ps per hit** (start→end of life) and **30→50 ps per track**; the draft's range looks like the per-track numbers attached to a per-hit label. Given the table's whole purpose is object-type honesty (per-MIP vs per-hit vs per-shower — a claims-law referee-defense pair), this is the most referee-visible benchmark risk. VERIFY against the TDR and fix either the label or the range. | Benchmarks (BTL/HGTD/SPACAL/Crilin) are new in this draft — no prior RADiCAL material cites them, so the TDRs themselves are the only check. |

### USEFUL

| # | Item | What the prior record has |
|---|------|---------------------------|
| U1 | **MCP-PMT reference.** The draft says only "A microchannel-plate photomultiplier (MCP-PMT) provides a common start" — no model, no timing spec, no citation. The (DW+UP)/2 absolute widths of Fig. 4 carry the MCP term, so its ~10–20 ps scale is load-bearing for interpreting the differential-vs-absolute gap. | Parent §4 (in text, verified): "MCP tube (Hamamatsu R3809U-50) with a timing resolution of 10 ps < σ_t < 20 ps"; parent §6 uses σ_t,MCP ≈ 10 ps. TNS [11]: J. Bortfeldt et al., "Timing performance of a micro-channel-plate photomultiplier tube", NIM A 960 (2020) 163592 — the prior record's formal MCP citation. |
| U2 | **`cates2018` at the front-end description (§2.2).** Both parent [22] and TNS [8] cite Cates et al. PMB 63 (2018) 185022 as the design basis of the high-gain differential amplifier ("modeled after a CERN design [22]"). The draft cites cates2018 only for photostatistics scaling in §6 — a different role. Add it where the HG/LG split is described; optionally let `gundacker` carry the §6 photon-count/emission-profile scaling alone (Cates 2018 is an electronics/SPTR paper; a TOF-PET referee may quibble at it as "photostatistics canon"). | Parent [22] = TNS [8] = draft `cates2018` (full authors per TNS: Cates, Gundacker, Auffray, Lecoq, Levin). |
| U3 | **FCC-ee motivation reference.** The intro's collider list cites cmsmtd/hgtd/hgcal/fcchh/aleksa — FCC-hh only on the FCC side. The parent's stated motivation is "FCC-ee, FCC-hh, fixed target and forward" and it cites the FCC-ee CDR first. | Parent [1] FCC physics opportunities, EPJ C 79(6) (2019) 474; parent [2]/[3] FCC-ee: The lepton collider (CDR vol 2), EPJ ST 228(2) (2019) 261–623. (NB parent [2] and [3] are duplicate citations of the same document — import once.) |
| U4 | **Name the SiPM (Hamamatsu HDR2) in §2.2**, citing the parent. The draft never names the photosensor; the parent does (in text, no separate bibitem needed). | Parent §4 text, verified in PDF. |
| U5 | **Append arXiv:2401.01747 to the `radical` bibitem.** Both FCC-MIT talk slides label the CERN 25–150 GeV results "arXiv:2401.01747, Submitted to NIM A" — i.e., the parent's preprint. Standard Elsevier practice permits "NIM A 1068 (2024) 169737; arXiv:2401.01747" and it serves readers without journal access. VERIFY the mapping before adding. | [FCC] slides 35 and 37. |

### OPTIONAL

| # | Item | Prior-record source |
|---|------|---------------------|
| O1 | Programmatic anchors: DOE BRN report (Fleming) — cited by ALL THREE prior publications (parent [6], INSTR [1], TNS [1]); ECFA roadmap (parent [7], CDS 2784893); CPAD/RDC links (parent [8–10]). The draft's "RADiCAL R&D programme" sentence (§1) could carry one of these. | Parent [6–10]. |
| O2 | Hu et al. LYSO/SiPM radiation-damage trio (INSTR [7–9]: IEEE TNS 65 (2018) 1018; 67 (2020) 1086; 68 (2021) 1244) — only needed if a component-level radiation-tolerance sentence is added; claims-law item 8 anticipates exactly this "tolerance BY CITATION" pattern. Currently no such sentence exists, so no action required. | INSTR [7–9]. |
| O3 | CAEN DT5742 product reference (parent [23]) — the draft already cites drs4 + ritt2014, which is stronger (peer-reviewed); adding the vendor URL is a style choice. | Parent [23]. |
| O4 | ROOT (Brun & Rademakers, TNS [13]) — common NIM A courtesy citation. | TNS [13]. |
| O5 | Ledovskoy simulation reports (parent [12,13], notredame.box.com) — gray literature behind the ≤20 ps goal and shower-max method; adequately inherited through the `radical`/`instruments2022` chain. Leave out unless a goals sentence is added. | Parent [12,13], INSTR [4,5]. |
| O6 | A shashlik-origin citation for '"shashlik" type' (§2.1). No precedent: the parent says "Shashlik/Kebab-style" uncited, so omission is continuity-consistent. If added, the reference must be independently verified — none exists in the prior record to copy. | None (VERIFY if pursued). |
| O7 | GEANT4 trio (TNS [4–6]) — NOT needed: this draft contains no simulation. Flag only so nobody "fixes" it; if a G4 floor-decomposition campaign lands later, these are the inherited refs. | TNS [4–6]. |

## 2. Exact placement + one-line reason

- **M1 (DSB1 refs)** → §2.1, first DSB1 mention ("\textbf{DSB1}, high-light organic WLS", ~line 134), and/or §1 at "the organic DSB1~\cite{luag,lucchini2013}" (~line 98–99). Reason: the paper's central variable needs the same provenance footing as its comparator, and the "intrinsically slower" kinetics claim needs the DSB1 τ=3.5 ns source, which currently rests on two LuAG papers.
- **M2 (`ruchti1996`)** → if kept: §2.1 at "read out longitudinally by WLS capillaries" (~line 125). Reason: fiber/WLS-readout lineage; otherwise delete the bibitem to avoid an unreferenced numbered entry.
- **M3 (HGTD)** → Discussion ~line 549 AND `tab_benchmarks.tex` line 22–23. Reason: per-hit vs per-track mislabeling in the one table whose stated purpose is object-type honesty.
- **U1 (MCP)** → §2.2, "(MCP-PMT) provides a common start" (~line 146–147). Reason: the absolute-timing panels (Fig. 4 middle) include the MCP term; name R3809U-50 and cite `radical` (its in-text 10–20 ps) and/or add Bortfeldt NIM A 960 (2020) 163592.
- **U2 (cates2018 at FE)** → §2.2, "split into a high-gain (HG) copy... and a low-gain (LG) copy" (~line 141–143). Reason: continuity — parent and TNS both cite Cates 2018 as the amplifier design basis.
- **U3 (FCC-ee)** → §1, first-sentence citation cluster `\cite{cmsmtd,hgtd,hgcal,fcchh,aleksa}` (~line 78). Reason: parent's primary motivation lists FCC-ee first; the draft currently motivates only via FCC-hh.
- **U4 (HDR2)** → §2.1/§2.2 at "silicon-photomultiplier (SiPM) channels". Reason: device anchoring, cite `radical`.
- **U5 (arXiv number)** → bibitem `radical`. Reason: open-access companion of the single most-cited reference (8 cites).
- **O1** → §1, "The RADiCAL R\&D programme develops..." sentence (~line 81). Reason: programme framing used by every prior RADiCAL publication.

## 3. Bibliographic inconsistencies (draft vs prior record)

1. **`radical` byline**: draft writes "C.~Perez-Lara, J.~Wetzel, U.~Akgun, et al. (RADiCAL Collaboration)". The published paper's byline (verified, p.1) is 39 named authors with **no collaboration designation**. Drop the parenthetical or VERIFY the collaboration credit against INSPIRE/journal record before keeping it.
2. **Internal comment defect**: lines 3–4 and 25–26 of the tex call the parent "Ref.~[2]" — it is **[1]** in the current bibliography. Comments only (not rendered), but a circulation-confusion hazard for coauthors completing the author block.
3. **`radhard2022` title variants**: draft (matching parent [14]) has "Precision timing, ultracompact, radiation-hard electromagnetic calorimetry"; TNS [3] cites it as "RADiCAL: Precision-timing, ..." and INSTR [2] as "RADiCAL—Precision-timing, Ultracompact, Radiation-hard Electromagnetic Calorimetry". The prior record disagrees with itself — VERIFY the actual arXiv:2203.12806 title and adopt it exactly.
4. **`instruments2022`**: published title uses an em-dash ("RADiCAL—Precision Timing, ...", verified on the PDF); draft uses a colon ("RADiCAL: Precision timing, ..."). Trivial; align.
5. **`wetzel2023`**: published title is "...A Radiation Hard Innovative EM Calorimeter" (no hyphen, verified); draft hyphenates "radiation-hard". Trivial; align.
6. **`fcchh`**: parent [4] includes the issue number "228 (4)"; draft omits "(4)". Trivial.
7. **`aleksa`**: the rest of the prior record (TNS [2], both talks) cites it by report number CERN-FCC-PHYS-2019-0003; the draft uses arXiv:1912.09962 only. Suggest carrying both identifiers.
8. **`lucchini2013`**: the pre-resolution draft (commit 7cd8e1f) had "M.T.~Lucchini" with a marked verify-note; the current entry says "M.~Lucchini" and the verify-note was removed. VERIFY initials/title (M.T. is the likely correct form).
9. **`cantone2023`**: title "Crilin beam-test results" reads as a paraphrase, not a journal title — Front. Phys. 11 (2023) 1223183 / arXiv:2308.01148 metadata are plausible, but the exact published title and author list were never verified in-repo (this was one of the four TODO-P2-CITE placeholders; `spacal`, `crilin`, `gundacker` were resolved with full, plausible metadata — `gundacker` PMB 64 (2019) 055012 is correct — but `cantone2023`'s title still looks unfinished).
10. **Caution when importing parent refs**: parent [2] and [3] are duplicate citations of the same FCC-ee CDR volume 2 — import only once (relevant to U3).

### Benchmark-usage check (BTL / HGTD / SPACAL / Crilin)
- **CMS BTL**: correctly cited to the MTD TDR (CERN-LHCC-2019-003) and correctly labeled "design, begin–end of life" in the table. One wording nit: intro line 79 says LYSO timing layers "**reach** 30–60 ps" while citing only a TDR — "target / are designed to achieve" is safer unless a measured beam-test reference (verified metadata only) is added.
- **ATLAS HGTD**: report number CERN-LHCC-2020-007 correct; the 35–50 ps **per hit** range is the suspect item (M3).
- **SPACAL**: An et al., NIM A 1045 (2023) 167629 / arXiv:2205.02500 — metadata consistent with the known publication; the 18.5±0.2 ps @ 5 GeV and ~30 ps light-guide variant match that paper's W-GAGG results; "situates... without ranking" complies with the claims-law "never fastest/best" rule. Sound.
- **Crilin**: Ceravolo JINST 17 (2022) P09033 correct; the "<100 ps for >1 GeV (Proto-0); <25 ps for >3 GeV" split is plausible but should be verified against the actual cantone2023 text when its title is fixed (item 9).
- None of these four systems were cited in any prior RADiCAL material (greps of parent/INSTR/TNS reference lists and both talks confirm) — they are new, deliberate additions, and apart from M3 they are handled *more* carefully than typical prior art.

### Intro anchoring verdict
**Sufficiently anchored, and the load-bearing attributions check out against the parent PDF:**
- All four prior RADiCAL publications cited in the opening programme sentence (§1).
- "That work also identified, as a priority for further study (its Section 7), a direct comparison of... WLS materials" — **verified**: parent §7 item 2 is verbatim "Data Analysis: Comparison of DSB1 and LuAG:Ce Wavelength Shifters... in 25 GeV steps... over the range 25 GeV ≤ E ≤ 150 GeV". This is the paper's strongest continuity asset.
- "a=256 ps√GeV, b=17.5 ps, ≈27 ps at 150 GeV" — verified against the parent's Eq. (2) and abstract; floor verbs ("confirms") comply with the claims law throughout.
- "the ~0.2 ns satellite also noted in Ref.~\cite{radical}" — **verified**: the parent discusses a small satellite peak shifted by ≈0.2 ns.
- "as already foreseen in Ref.~\cite{radical} (its Section 7.3)" — the content is real (parent §7 item 3: doubling LYSO:Ce thickness at shower max), but the parent formats §7 as a numbered list, not subsections; "Section 7, item 3" would be the precise pointer. Minor.

## 4. VERIFY — could not be confirmed from the materials at hand

| Item | What is known / what to check |
|------|-------------------------------|
| HGTD per-hit vs per-track ranges (M3) | High-confidence recollection: 35→70 ps per hit, 30→50 ps per track (CERN-LHCC-2020-007). Check the TDR directly before editing. |
| Exact arXiv:2203.12806 title (`radhard2022`) | Prior record self-inconsistent (with/without "RADiCAL" prefix). Pull the arXiv abstract page. |
| `cantone2023` exact title + authors | "Crilin beam-test results" is a paraphrase; Front. Phys. 11 (2023) 1223183 / arXiv:2308.01148 plausible but unverified in-repo. |
| arXiv:2401.01747 = parent preprint (U5) | Strong evidence from [FCC] slides 35/37 ("Submitted to NIM A"); confirm the arXiv page matches NIM A 1068 169737 before adding. |
| "LYSO-based calorimetric devices have demonstrated ~30 ps shower timing" (`anderson2015`, intro line 80) | D. Anderson et al. NIM A 794 (2015) 7–14 is the right paper (also TNS [14]); confirm the ~30 ps figure is what it reports at high amplitude. |
| `ritt2014` page range 3607–3617 | Volume/issue (TNS 61(6), 2014) believed correct; pages unverified. |
| Parent [21] OSTI purl/924745 (for M1) | Parent gives only the title-as-URL "Waveshifters and Scintillators for Ionizing Radiation Detection"; pull the OSTI record for author/year before creating the bibitem — do not invent metadata. |
| "FTBF / A. Bornheim et al." ([CPAD] slide 37) | Informal slide credit, no formal reference anywhere in the prior record's bibliographies; most plausibly subsumed by `anderson2015` (Bornheim is Caltech/CMS-timing) and `wetzel2023`. No action unless a coauthor identifies a citable publication. |
| Vendor-URL bibitems (Eljen, CAEN) vs footnotes | NIM A accepts both; parent used numbered bibitems — following the parent is the low-risk choice. Style decision for the corresponding author. |

**Claims-law compliance of these recommendations:** none of the above revives the retired 0.99/χ²=0.4 ratio, uses forbidden parity language, claims universal WLS equivalence, 5D, calibrated z, or intrinsic position resolution; all proposed numbers defer to Table 1/Table 2/the corrected MIXED figure, and every benchmark suggestion preserves the "situate, never rank" framing.


---

# Agent 4 — Technical consistency

# Technical-consistency audit — RADiCAL timing manuscript (commit cb19f51) vs the prior record

Sources checked: [PARENT] NIM A 1068 (2024) 169737 (all 15 pp); [INSTR] Instruments 6(3) 27 (2022) (all 6 pp); [TNS] IEEE TNS 70(7) 1296–1300 (2023) (all 5 pp); [CPAD] CPAD/SLAC Nov 2023 talk (38 slides); [FCC] FCC-US/MIT Mar 2024 talk (49 slides); current manuscript `papers/timing/radical_timing.tex`; claims law `papers/memory_claims_and_forbidden_language.md`; authoritative 2023 channel map / working points (memory `dataset_2023.md`, author-confirmed 2026-06).

## 1. Table of technical inconsistencies

| # | Item | Prior source + value | Current manuscript value | Severity | Active / historical | Recommended action |
|---|---|---|---|---|---|---|
| 1 | **LG digitisation rate** | Authoritative 2023 channel map (author-confirmed): DRS0 @ **5 GS/s** carries the 8 HG channels + MCP copies; DRS1 @ **1 GS/s** carries all 8 LG channels + wire chambers. ([PARENT] loosely says "All signals were digitized with 5GS/s", which the channel map supersedes.) | Sec. 2.2: HG and LG copies "**both** digitised by DRS4 … (CAEN DT5742) **at 5 GS/s** over a 1024-cell window" | **High** | Active | Split the sentence: HG timing channels digitised at 5 GS/s (0.2 ns/cell); the LG copies on a second DT5742 at 1 GS/s (used for amplitudes, not edge timing). This is method-relevant — srCFD's amplitude anchor (Eq. 1) is read from the 1 GS/s LG waveform, and a referee/coauthor checking LG peak sampling will ask. Note the parent's "all at 5 GS/s" wording should not be repeated. |
| 2 | **Estimator attributed to Ref. [1]** | [PARENT] Sec. 5.2 + Fig. 12 caption: "**fixed thresholds** were set for the high-gain signals … fast timing using the **fixed threshold on the rise of the leading edge**"; best estimator = "BestMinus" (Δt_DW−Δt_UP)/2 on energy bins 6–8. [TNS]: per-channel-optimized fixed thresholds. | Sec. 5.3: "cfd05 (a CFD on the clipped measured peak, **as in Ref.~\cite{radical}**)"; Table 1 published-DSB1 row lists **src = cfd05**. (Internal memory also asserts "cfd05 on the clipped pulse" for the parent — same conflation.) | **High** | Active | Verify the actual published chain with C. Perez-Lara, but as published the parent used a fixed-threshold leading-edge time, not a 5% CFD. Either relabel the Table 1 published row ("LED/fixed-threshold, publ.") and reword "as in Ref. [1]" → "cfd05, our clipped-peak CFD diagnostic, which reproduces the published chain to ~1.6 ps", or add a footnote stating the published estimator in the parent's own terms. Also note the published 27 ps used the parent's bins-6–8 BestMinus selection, not brightest-1000 — the comparison row mixes selections and should say so. The locked numbers (27/17.5/256) are untouched; only the estimator label changes. |
| 3 | **SiPM count** | [PARENT]: front-end cards "each housed **four** Hamamatsu HDR2 SiPM" × 2 cards = **8 SiPMs**; "eight high-gain readout channels", "all eight low-gain SiPM signals". [TNS]: "Each SiPM had two readout circuits" (LG + HG). | Sec. 2.1: "$4\times2 = 8$ SiPM channels" — but Discussion: "on **16 SiPMs** in a $14\times14$ mm² cell" | **High** | Active (internal contradiction + prior conflict) | The module has 8 SiPMs read out in 16 gain channels (8 HG + 8 LG). Fix the Discussion sentence ("on eight SiPMs (sixteen readout channels) in a 14×14 mm² cell") — and flag that the claims-law line "14×14 mm², 16-SiPM module" encodes the same conflation and needs the same correction after hardware confirmation. |
| 4 | **~0.2 ns satellite provenance** | [PARENT] Sec. 5.2.1/Fig. 19: satellite = DRS4 one-sample (≈0.2 ns) threshold-crossing shift, ~15 % of channels, visible in **MCP-referenced** averages; "It does **not** appear … in the timing difference between downstream and upstream SiPM channels" (i.e. not in (DW−UP)/2). | Sec. 4: "a small non-Gaussian shoulder (the ∼0.2 ns satellite **also noted in Ref.~\cite{radical}**)" excluded from the **t_DW−UP** fit; Secs. 5.3/App. B attribute the satellite to cfd05 clipping, removed by srCFD. | **Med-High** | Active | The two satellites are arguably different objects (parent: MCP-referenced sample-shift; here: clipping-induced cfd05 tail in the differential). A coauthor will quote the parent's "does not appear in DW−UP" sentence back at you. Qualify the cross-reference ("a structure analogous to the ≈0.2 ns satellite of Ref. [1], there seen in MCP-referenced distributions") or drop it and let the digitizer-forensics result stand on its own. |
| 5 | **TENERGY description** | Dataset record: TENERGY = **3 × DSB1 T-type + 1 × full-length E-type (energy) capillary at NW**; prior papers define the T-type/E-type capillary vocabulary ([PARENT] Sec. 2, [INSTR] Sec. 2, [TNS] Sec. II). | Sec. 2.1: "a **multi-segment** high-light build with an added **longitudinal energy-readout section**"; Table 1 WLS = "organic, seg". | **Med** | Active | "Multi-segment/added section" suggests extra module hardware; it is a capillary swap. Rewrite using the priors' vocabulary: "TENERGY: three DSB1 timing (T-type) capillaries plus one full-length energy (E-type) capillary; it therefore times with three capillaries" — which also makes the √(4/3) penalty self-explanatory. |
| 6 | **T-type capillary construction absent** | [PARENT]: T-type = quartz capillary 183 mm, OD 1150 µm / bore 950 µm, **15 mm DSB1 filament (900 µm) at shower max**, quartz rods elsewhere; DSB1 = Notre Dame/Eljen dye, abs 425 nm / em 495 nm, τ = 3.5 ns. | Setup says only "WLS capillaries sample shower maximum"; Sec. 7 later invokes "the fixed-position wavelength-shifter **filament**" without it ever being defined. | **Med** | Active (omission) | Add 1–2 sentences in Sec. 2.1 defining the T-type construction (localized 15 mm WLS filament at shower max in a quartz waveguide, cite [1,3,4]) and DSB1/LuAG:Ce optical basics. This grounds both the depth section's "filament" and the MIXED corner discussion. |
| 7 | **Pb-glass backing + trigger counters omitted** | [PARENT]: A1 (1×1 cm²)/A2 (2×2 cm²) trigger counters, **A2 = primary trigger**; 2×2 Pb-glass array (4×4×40 cm³) as leakage/non-EM tagger. Production selection (SelectionCuts.h) includes a containment veto (reject if ΣPb > 0.30·ΣLG). | Sec. 2.2/Sec. 4: only DWC fiducial + valid MCP mentioned; no trigger counters, no Pb glass, no containment cut. | **Med** | Active (omission) | If the containment veto is in the production timing chain, it must appear in the Sec. 4 selection list and the Pb-glass array in Sec. 2.2 (one sentence each, cite [1]). If it is disabled for the timing chain, confirm and leave out — but verify before circulation. |
| 8 | **Photosensor / MCP models unnamed** | [PARENT]/[INSTR]/[TNS]: **Hamamatsu HDR2** SiPMs; MCP-PMT **Hamamatsu R3809U-50** (10 ps < σ < 20 ps [PARENT]; ≈10 ps [TNS]). | "silicon-photomultiplier (SiPM)" and "microchannel-plate photomultiplier (MCP-PMT)", no models. | **Med** | Active (omission) | Name both (HDR2; R3809U-50 with its 10–20 ps class, noting the differential estimator cancels it). Coauthors Ruchti/Zhu will expect the hardware identity; it also harmonizes with all three prior papers. SiPM bias (42.25 V both ends) is also worth one clause for reproducibility. |
| 9 | **LuAG light-deficit factor** | Dataset record: LUAG build "dim (~**3×** less)"; MIXED LuAG corners collect ~2.6–3× less LG light than DSB1 corners. | Discussion: "recovering the factor **∼2.2** light deficit … would bring a → ≈300 ps√GeV". | **Med** | Active | 2.2 is the stochastic-term *ratio* (440/203), not a measured light ratio; under photostatistics (a ∝ 1/√N) a ×2.2 in `a` implies ×4.7 in light, while the measured LG ratio is ~3×. Restate precisely: quote the measured ~3× LG amplitude ratio, and present the projection as "doubling the detected light gives a → a/√2 ≈ 310 ps√GeV" (which is what the ≈300 actually assumes). As written, a numerate coauthor will catch the unit mismatch. |
| 10 | **Beam line not named** | [PARENT]: CERN SPS **H2** beam line, CMS/HF motion table; ≥10⁶ triggers per energy step; wire chamber ≈3 m upstream. ([INSTR]: H4 for the 2015 4×4 array; [TNS]: FNAL MT6.) | "CERN SPS (May 2023)". | Low | Active (omission) | Say "H2 beam line of the CERN SPS"; optionally add the trigger statistics sentence. Avoids any confusion with the earlier H4 (array) campaign. |
| 11 | **MCP role wording** | [PARENT]: A2 counter was the primary trigger; the MCP is an independent precision **timing reference**. | "A MCP-PMT provides a common **start**". | Low | Active | "Common start" reads as trigger. Prefer "provides a common time reference; the A2 scintillator triggered the readout" (matches parent). |
| 12 | **"its Section 7.3"** | [PARENT] Section 7 is a numbered list (items 1–4); there is no subsection 7.3. Item 3 contains the ×2 LYSO-thickness idea. | Discussion cites "Ref.~\cite{radical} (its Section~7.3)". | Low | Active | Change to "its Section 7" or "Section 7, item 3". (The intro's plain "its Section 7" for the WLS-comparison priority — item 2 — is correct.) |
| 13 | **Module geometry numbers across priors** | [PARENT] (canonical): 14×14×135 mm³, 29 LYSO (1.5 mm) + 28 W (2.5 mm), ≈25 X₀, **X₀ = 5.4 mm**, R_M = 13.7 mm. [TNS] text: 13 cm long, **X₀ = 4.7 mm**; [INSTR]: depth **114 mm** = 25 X₀ ("X₀∼4.5 mm"); the 114 mm schematic is recycled in [TNS] Fig. 1, [CPAD] s13, [FCC] s12; both talks state "shower-max radius ~4.5 mm = X₀". Capillary length: 18 cm [TNS] vs 183 mm [PARENT]. | Manuscript states only "≈25 radiation lengths deep" and "14×14 mm²" — no conflicting number. | Low (Med if dims get added) | **Historical-only** (priors disagree among themselves) | Current draft is safe by silence, but (a) if dimensions are added, use only the parent's canon (14×14×135, 29+28, X₀ = 5.4, R_M = 13.7); (b) the depth section's "−X₀/v_eff ≈ −26 ps per e-fold" implicitly uses X₀ = 5.4 mm and v_eff ≈ 200 mm/ns — state both numbers so the −26 ps is checkable and unambiguously tied to the parent's X₀, not the talks' 4.5 mm. |
| 14 | **Corner labels undefined** | [PARENT] channel naming DWNW/DWNE/DWSW/DWSE, UPNW/… (Figs. 14, 18); dataset map uses NE/NW/SE/SW × D/U. | NE/NW/SE/SW first appear in Sec. 6 (t_NE − t_SW) with no prior definition; Sec. 2 defines only DW/UP. | Low | Active | Add the corner convention where DW/UP are defined in Sec. 2.1 ("the four capillary positions are labelled by compass corner NE/NW/SE/SW, viewed from upstream"), matching the parent's channel names. |
| 15 | **Header comments cite "Ref. [2]"** | The parent is \bibitem #1 (`\cite{radical}`) in the manuscript's own bibliography. | Lines 4 and 25 (comments): "as in Ref.~[2], NIM A 1068 (2024) 169737". | Low | Active (comment-only, not rendered) | Fix to Ref. [1] to avoid confusing coauthors reading the source. |
| 16 | **Keywords** | Claims law: 4D demonstration language is Paper-3 scope; "toward 5D" only where already used. | Keywords include bare "**4D calorimetry**"; also "LYSO" (house style elsewhere is LYSO:Ce). | Low | Active | Consider "toward 4D calorimetry" or drop; "LYSO:Ce" for the keyword. |
| 17 | **Saturation as deliberate choice** | [PARENT]: "the high gain signals were **driven further into saturation**, to maximize the signal rise of the leading edge". | Sec. 3 presents clipping purely as an obstacle that srCFD overcomes. | Low | Active (framing) | One clause acknowledging the deliberate operating choice ("the HG gain was intentionally set so that showers clip, maximizing edge slew [1]; srCFD recovers the amplitude information this sacrifices") pre-empts a "why did you let it clip?" referee/coauthor question and matches the parent's framing. |
| 18 | **Talk-only claims** | [CPAD] s31 / [FCC] s40,43: "the RADiCAL module **is radiation hard**", "concept **demonstrated** to meet the needs of high radiation … environments", "potential to reach <10 ps", "45 ps @ 28 GeV" (TNS published 42–43 ± 4.3/4.4, 49.5 ± 5). [FCC] s39: "25 ps @ 150 GeV, limiting ~18 ps" (preliminary). | Manuscript does **not** propagate any of these (rad-hardness is by citation only; 25.7 ± 0.6 / 18.8 ± 0.8 are the gated numbers). | — | **Historical-only / obsolete** | No action — record that these talk claims are superseded and must never re-enter the draft (claims-law items 6, 8). The preliminary 17.52/255.58 fit in both talks is the same parent fit, correctly cited as 17.5/256 in the manuscript. |
| 19 | **LuAG:Ce emission wavelength** | [FCC] s14: LuAG:Ce **510 nm**; [FCC] s13 (Zhu table): λ_peak **520 nm**. | Not quoted. | Low | Historical-only | If optical numbers are added per item 6, pick one sourced value (Hu et al. [luag]) rather than the slides. |

## 2. Terminology harmonization notes

- **srCFD / hg_lgcfd / cfd05 / LED.** The manuscript's taxonomy (srCFD primary, cfd05 + LED diagnostic; branch name `hg_lgcfd` given once) is internally clean and matches the claims-law definition verbatim. What is missing is the **bridge to the prior record's vocabulary**: the parent speaks of "fixed thresholds on the leading edge" and of estimator combinations named *Downstream / Upstream / BestPlus / BestPlus-MCP / (DW−UP)/2 / BestMinus*. One sentence in Sec. 3 or 4 mapping "LED ≙ the fixed-threshold leading-edge times of Ref. [1]" and "Eq. (2) ≙ the (DW−UP)/2 / BestMinus family of Ref. [1]" would let every coauthor of the parent locate themselves immediately — and is required anyway once item 2 (cfd05 misattribution) is fixed.
- **T-type / E-type.** All three prior papers and both talks consistently use "timing (T-type) capillary" and "energy (E-type) capillary". The manuscript never uses either term, yet needs exactly this distinction for TENERGY (item 5) and the Sec. 7 "filament". Adopt the prior vocabulary once in Sec. 2.1.
- **DW/UP and corners.** Parent channel names DWNW…UPSE decompose into the manuscript's DW/UP × NE/NW/SE/SW — consistent, but define the corner half in Sec. 2 (item 14). Internal D/U ("Down/Up") never leaks into the draft — good.
- **Invariant wording rule** ("LYSO:Ce + W common; the WLS capillary is the variable") is correctly observed in the abstract, Sec. 2.1, and the Table 1 caption. No instance of "crystal differs". The MIXED kill-shot language, floor verbs ("confirms"), the 10–20 % parity form, the scatter-based 1.04 ± 0.05, the method-dependence sentence, the position-coupling caveat, and the corner-map provenance sentence all match the claims law exactly.
- **"Delay-wire chamber"** (with Spanggaard citation) is a fine harmonization of the parent's generic "beam/wire chamber"; keep.
- **"Shashlik"** spelling matches the priors (parent uses "Shashlik/Kebap-style" once; no action).
- **Adopted-source / per-regime rule** wording is consistent across Secs. 3, 5.3, Appendix A/B. Note that Table 1's "src" column (srCFD vs LED per build) is the surviving statement of the rule — after fixing item 2, make sure the published row's source label doesn't silently look like a fifth "regime".

## 3. What a prior-record-familiar coauthor would flag as wrong or unexplained

1. **"cfd05 … as in Ref. [1]"** (item 2). Perez-Lara and Ruchti wrote "fixed threshold on the rise of the leading edge" in the parent; seeing their estimator renamed to a 5 % clipped-peak CFD — in the comparison table that anchors the paper's continuity claim — is the single most likely circulation objection. The internal memory's own "cfd05 on the clipped pulse [arXiv:2401.01747]" assertion should be re-verified against the parent text at the same time.
2. **"16 SiPMs"** (item 3). Anyone who built the readout cards knows there are four HDR2s per card. The same person will note the claims-law sound-bite has the same error.
3. **"Both digitised at 5 GS/s"** (item 1). Whoever configured the DT5742s knows the LG/WC board ran at 1 GS/s; this is also the one place where the published parent record itself is loose, so the new paper is the place to state it correctly, not to inherit it.
4. **The satellite cross-reference** (item 4). The parent explicitly states the 0.2 ns satellite is absent from (DW−UP)/2; citing it as "also noted in Ref. [1]" inside a t_DW−UP discussion looks like a misreading of their own paper unless qualified.
5. **TENERGY as "multi-segment"** (item 5) — the builder knows it is a one-capillary swap (3 T-type + 1 E-type), and the priors' T/E vocabulary exists precisely for this.
6. **Where is the Pb glass?** (item 7). The parent devotes a paragraph to the backing calorimeter and the non-EM rejection it provides; its silent disappearance — while a Pb-glass-based containment veto sits in the production cut chain — reads as an unexplained selection difference.
7. **The −26 ps/e-fold expectation** (item 13): with the talks having broadcast "X₀ ≈ 4.5 mm" and the parent stating 5.4 mm, the depth section must pin its X₀ and v_eff numerically, or two coauthors will derive two different expectations (−22 vs −27 ps) and ask which one the dashed band is.
8. **The 2.2× vs 3× light deficit** (item 9): Ledovskoy/Zhang-class readers will square the stochastic ratio and notice 2.2 ≠ 4.7 ≠ 3; one precise sentence dissolves it.
9. **Missing hardware identities** (HDR2, R3809U-50, H2, bias, A2 trigger; items 8, 10, 11): not errors, but the parent's coauthors will perceive the setup section as under-specified relative to the standard their own paper set, and every one of these is a one-clause fix with a citation to [1].
10. **Lineage nuance** (low): the intro calls the parent "a first RADiCAL test-beam measurement"; strictly the FTBF June 2022 run ([TNS], 42 ± 4.3 ps at ≈28 GeV) was the first timing beam test of this module. A half-sentence acknowledging the FTBF demonstration before the parent's full energy scan respects the published lineage and costs nothing.

**Compliance confirmation (things checked and found correct):** a = 203±6 → 440±18 (×2.2) and 25.7±0.6 / 18.8±0.8 / 1.04±0.05 match the post-fix claims law; published 256 / 17.5 / 27 ps quoted exactly per the parent (Fig. 27: 255.58/17.52); "confirms, never revises" verbs respected; floors stated with the TENERGY √(4/3) and LuAG-extrapolation (44.4 ps measured) qualifiers; no "statistically indistinguishable"/"consistent with unity"; no 0.99/χ²=0.4; depth dial framed as acceptance-conditional consistency, not a calibrated z; position result absent (correctly deferred); 820 mV clip correctly attributed to the digitiser; energies 25–150 GeV with DSB1-only at 25 GeV matches the run manifest; 39-author/11-affiliation roster note matches the parent's masthead.


---

# Agent 5 — Figures & tables

# Agent 5 — Figures/tables reuse & divergence audit (RADiCAL timing paper, commit cb19f51)

Sources examined in full: parent NIM A 1068 (2024) 169737 [PARENT, 15 pp, Figs. 1–28]; Instruments 6(3) 27 (2022) [INSTR, Figs. 1–5]; IEEE TNS 70(7) 1296 (2023) [TNS, Figs. 1–10]; CPAD/SLAC 2023 talk (38 slides); FCC/MIT 2024 talk (49 slides); current manuscript PDF (11 pp; Figs. 1–8, A.9, A.10, B.11; Tables 1–3) and its `.tex` captions; the claims law; both 2026-06-09 memos; `papers/timing/figs/` inventory.

---

## 1. Figure reuse recommendations (prior figure → where it helps the draft)

**1a. Module schematic — VERIFIED GAP, significant.** The current paper has **no apparatus figure at all**: Figs. 1–3 are method figures (clip histogram, HG/LG calibration, single waveform), and Sec. 2 ("Experimental setup") is text-only. Every prior RADiCAL record leads with the module drawing: PARENT Fig. 1 (W/LYSO shashlik with quartz capillaries, monitoring fibers, SiPM readout card), INSTR Fig. 1a/b, TNS Fig. 1, CPAD slides 9/13, FCC slides 8/12. For a NIM A instrumentation paper this absence is anomalous, and it leaves the reader of Sec. 2.1 ("four capillaries traverse the stack... DW and UP") with no picture of the geometry the whole paper depends on. **Recommendation:** adapt PARENT Fig. 1 (with the T-type/WLS-filament detail of PARENT Fig. 6 if a single composite is preferred) as a new Fig. 1 in Sec. 2.1, captioned "adapted from Ref. [1]". It is the collaboration's own figure, so reuse is unproblematic; redrawing with the 2023 dimensions is even better.

**1b. Beam's-eye capillary/corner map — second gap, directly load-bearing for Sec. 6.** The MIXED result hinges on "DSB1 (NE, SW) and LuAG:Ce (NW, SE) on opposite diagonals", yet the corner labels NE/SW/NW/SE are never drawn. The visual ancestor exists: PARENT Fig. 2 (tile cross-section with the 4+1 hole pattern, 3.50 mm offsets) and especially INSTR Fig. 2 (beam's-eye schematic with two capillary types on opposite diagonals — exactly the MIXED topology). **Recommendation:** a small beam's-eye square with four labeled corners (colored by WLS material, central unused hole shown) either as an inset in Fig. 7(A) or as panel (b) of the new apparatus figure. This also visually supports the "same showers, same module, same DRS group" cancellation argument and the corner-map provenance defense (referee risk #7).

**1c. PARENT Fig. 7 (GEANT4 longitudinal shower profile) → cite in Sec. 7 (depth).** The depth-drift section currently leans on Longo–Sestili and the PDG. The program's own GEANT4 figure shows shower max moving from LYSO layers ~8–10 to ~11–13 across the energy range, with the fixed WLS-filament position marked — the exact mechanism behind both the −33.6 ps/e-fold drift and its acceptance-conditional caveat. **Recommendation:** cite it ("the GEANT4 study of Ref. [1] (its Fig. 7) shows shower maximum migrating past the fixed filament position over this energy range"), no need to reproduce.

**1d. PARENT Fig. 8 / INSTR Fig. 5b (simulated sigma_t vs detected light yield) → cite in Sec. 5.2 or Discussion.** The paper's thesis — stochastic term tracks detected light — has a simulation ancestor inside the program's own record (the log-log sigma_t vs LY line, "more light = better resolution" in both talks). One sentence, e.g. "consistent with the detected-light scaling anticipated in the simulation studies of Refs. [1,3]", roots the thesis in prior RADiCAL expectation and preempts "post-hoc interpretation" objections. Do not reuse the figure itself (simulation, single-end readout, different conditions).

**1e. PARENT Fig. 11 (beamline arrangement: A1/A2 triggers, MCP, module, Pb-glass) → cite in Sec. 2.2.** Sec. 2.2 describes the MCP reference but not the beamline; one clause — "the H2 beamline arrangement of Ref. [1] (its Fig. 11)" — saves a figure and answers setup questions by reference.

**1f. PARENT Sec. 4 / Fig. 12 (deliberate saturation) → fold into Sec. 2.2 or Fig. 1 caption.** The parent states the HG signals "were driven further into saturation, to maximize the signal rise of the leading edge". The current paper presents clipping as a found condition; stating it was a *deliberate design choice of Ref. [1]* both strengthens continuity and defuses the "you saturated your SCA digitizer" objection (referee risk row 3). Exact text in Sec. 3 below.

**1g. Benchmark Table 3** is fine as built; the talks' CMS-BTL/HGTD comparator framing (CPAD slide 8, FCC slide 48 FCC-hh dose table) is already superseded by the table's object-type column. No change.

Not recommended for Paper 1 (Paper 3 material; do not pull in): PARENT Figs. 13–17 (per-capillary (x,y) maps, energy linearity, sigma_E/E), INSTR Figs. 3–4 (4×4 array, 15.7%/sqrtE energy resolution), TNS Figs. 7–9, FCC slides 35–36.

---

## 2. Figure retirement list (must NOT return)

1. **The "Preliminary" sigma_t(E) plot (CPAD slide 30, FCC slide 38)** — same data as PARENT Fig. 27 right but watermarked `Preliminary` with build-label strip `M2925W2815L-415DSB1-T`. Superseded by the published figure; any reuse must be of the published version, cited as Ref. [1].
2. **FCC slide 39 bullet: "25 ps @ 150 GeV, with limiting resolution of ~18 ps."** Preliminary numbers that match neither the published record (27 / 17.5 ps, cfd05) nor exactly the current post-fix values (25.7±0.6 / 18.8±0.8, brightest-1000 srCFD). Their numerical closeness to the new headline is coincidental (different selection/estimator); quoting the talk would blur the locked "published 27/17.5 unchanged; new numbers confirm" line. Never cite the talk for numbers.
3. **CPAD slide 31 / FCC slide 40–41 / CPAD slide 33 claims**: "The RADiCAL module is radiation hard" (asserted from this data), "meets needs of FCC EndCap", "potential to reach <10 ps at >150 GeV", "demonstrated to meet the needs of high radiation and high luminosity environments". All forbidden under the claims law (items 1, 8); none may migrate into captions, conclusions, or cover-letter text. Hardness remains component-level by citation (FCC backup slides 46–47 = the Caltech RIAC and Co-60 studies → cite Hu et al. and parent refs, do not reproduce in Paper 1).
4. **CPAD slide 28 energy linearity "29.8 mV/GeV (Preliminary)"** — superseded by the published 29.4 mV/GeV (PARENT Fig. 17); Paper 3 territory anyway.
5. **TNS Fig. 10 (42±4.3 ps at ~28 GeV)** — not wrong, but it is an MCP-referenced sigma(t_RAD − t_MCP), i.e. t_sum-class. If ever mentioned, it must carry the taxonomy label; never place it beside the reference-free differential numbers as if comparable (claims law: both numbers always labeled by taxonomy).
6. **PARENT Fig. 27 left ("BestMinus/BestPlus" six-method plot)** — useful prior art but its method names conflict with the current estimator taxonomy (srCFD/LED/cfd05); if referenced, add the mapping "(DW−UP)/2, called BestMinus in Ref. [1]" once. Current Sec. 1 cites Ref. [1] generically — adequate.
7. **Local stale files `figs/thesis.png` and `figs/method.png` (pre-fix, June 7–8 timestamps)** — no longer referenced by the tex (which uses `thesis_postfix.png`, `method_postfix.png`) but still present; delete or archive so a future edit cannot silently re-include pre-fix numbers (the exact failure mode the Fig. A.10 episode already demonstrated).
8. **Any internal MIXED ratio plot from the retired 0.99/chi2=0.4 era.** No *published* figure ever showed it (verified: it appears in none of the five public sources), so retirement is an internal-discipline item: `mixed_h2h_corrected.png` is the only admissible MIXED figure.

---

## 3. Caption improvements for current figures (exact text)

**3a. Fig. 5 (money plot) — visual continuity with the published sigma_t(E) style.** PARENT Fig. 27 right is the established look: data points *with error bars*, one smooth fit curve, "Time Resolution (ps)" vs "Beam Energy (GeV)". Current Fig. 5 shows four curves with no visible per-point statistical errors — the one place the new paper looks *weaker* than its predecessor. The in-file TODO-P2 (error bars + requirement bands + dashed MIXED) already prescribes the fix; executing it is the single highest-value figure change. Additionally append to the caption:
> "The published DSB1 parametrisation of Ref. [1] ($256/\sqrt{E} \oplus 17.5$ ps, clipped-peak discriminator) is shown for reference (grey dashed); the present DSB1 fit is consistent with it (Table 1)."
(and draw that one grey curve — it makes "confirms, never revises" *visible* rather than only tabular).

**3b. Table 1 caption — the cfd05 attribution needs a continuity clause (real divergence found).** The parent paper's own text describes its timing as an **optimized fixed-threshold leading-edge time** ("fixed thresholds were set for the high-gain signals…", Sec. 5.2; Fig. 12 caption: "fast timing using the fixed threshold on the rise of the leading edge"). The current manuscript twice writes "cfd05 … as in Ref. [1]" (Sec. 5.3) and labels the published row "cfd05" (Table 1). A referee who reads both papers will flag this as a misattribution. Without relitigating the internal equivalence (Table 1 is authoritative), the wording should be hedged precisely. Table 1 caption, replace "The published DSB1 values (cfd05) are shown for comparison." with:
> "The final row reproduces the published DSB1 result of Ref.~[1] ($a=256$~ps$\sqrt{\mathrm{GeV}}$, $b=17.5$~ps, $\sigma_t^{150}\approx 27$~ps), obtained there with a clipped-pulse threshold discriminator; cfd05 is the present chain's implementation of that estimator class, and reproduces the published values on this dataset."
And in Sec. 5.3, change "cfd05 (a CFD on the clipped measured peak, as in Ref. [1])" to "cfd05 (a CFD on the clipped measured peak, the estimator class of Ref.~[1])". (If the gate record documents exact equivalence differently, use its wording — but "as in [1]" unqualified is a defect.)

**3c. Fig. 1 (clip.png) caption — add the deliberate-saturation provenance.** Append:
> "The high-gain chain is deliberately driven into saturation to maximise the slope of the recorded leading edge (Ref.~[1], Sec.~4); Sec.~3 recovers the clipped peak event by event."
This converts an apparent defect into a documented design choice with a published antecedent.

**3d. Fig. 7 (mixed_h2h_corrected) caption — pin the corner geometry.** After "(corner map confirmed by pulse-shape discriminants)" insert:
> "the four capillaries sit at the corners of the $14\times14$~mm$^2$ face (capillary-placement pattern as in Ref.~[1], its Fig.~2), DSB1 on the NE--SW diagonal and LuAG:Ce on the NW--SE diagonal"
— or, better, add the small beam's-eye inset of Sec. 1b and reference it.

**3e. Fig. B.11 caption — satellite continuity.** Change "(the $\sim$0.2 ns satellite of Sec. 4)" linkage to carry the prior-record anchor once in the appendix too:
> "…the common tail window … is marked. The $\sim\!0.2$~ns satellite structure was first noted in Ref.~[1] (its Fig.~19); here it is shown to be clipping-induced and removed by srCFD."
(Sec. 4 already says "also noted in Ref. [1]" — good; the appendix figure should be self-contained since referees read appendices in isolation.)

**3f. Fig. 8 (depth_dial) caption** — already carries "acceptance-conditional; see text" (compliant). Optionally add the parent anchor: "Shower-maximum migration with energy in this stack is shown in the GEANT4 study of Ref.~[1] (its Fig.~7)." No other change; do not add anything resembling a calibrated-z reading.

---

## 4. Provenance clarifications needed (post-fix / corrected / supersedes)

1. **Do NOT add "post-fix" or "corrected" language to any public caption.** Verified: none of the pre-fix numbers (25.3 ps, a=201/455, the 0.99 MIXED ratio) ever appeared in the five public sources — the only public prior numbers are the parent's 27/17.5/256 (and the talks' preliminary variants). Public continuity is therefore fully handled by the Table 1 published row + "confirms" verbs; in-paper "corrected/post-fix" labels would advertise an internal bug with no public antecedent and invite questions the paper need not answer. Keep provenance in the circulation note (where it is already excellent: items 1–2 of "What changed") and in `FORMAT_AUDIT`/gate files.
2. **Internal filenames vs captions are correctly decoupled** (`thesis_postfix.png`, `method_postfix.png`, `mixed_h2h_corrected.png` render with neutral captions). Maintain this; also remove the stale pre-fix `thesis.png`/`method.png` from `figs/` (Sec. 2 item 7) so the postfix/legacy distinction cannot regress.
3. **Fig. A.10 (systematics)**: the 2026-06-10 regeneration from the post-fix chain (figure and Table 2 now agree by construction; caption defers to Tables 1–2 as authoritative) is the right pattern — figure annotations should never be an independent source of numbers. Confirmed the rendered PDF carries the regenerated version (totals ±1.0/±1.1/±0.9/±1.9 ps).
4. **For coauthors who saw the FCC-2024 talk**: the circulation note should add one explicit mapping line so nobody "remembers" the slide-39 numbers into review comments, e.g. "The FCC-MIT 2024 slide quoted preliminary values (25 ps @ 150 GeV, ~18 ps limit); the published record remains 27/17.5 ps (cfd05, Ref. [1]), and the present brightest-1000 srCFD values are 25.7±0.6 / 18.8±0.8 ps — different selection and estimator; do not interchange." This is a note-level clarification, not manuscript text.
5. **MIXED module-wide row** (Table 1) is already labeled reference-only in both caption and text — no further provenance flag needed.

### Quick punch list (priority order)
1. Add apparatus figure (parent Fig. 1 adaptation) + beam's-eye corner map (Sec. 1a/1b).
2. Execute the Fig. 5 TODO-P2 restyle (error bars, bands) + grey published-fit reference curve (Sec. 3a).
3. Fix the "cfd05 … as in Ref. [1]" attribution in Sec. 5.3 and Table 1 caption (Sec. 3b) — the only true text-vs-prior-record contradiction found.
4. Add the deliberate-saturation sentence with Ref. [1] anchor (Sec. 3c).
5. Cite parent Figs. 7/8/11 at the three marked points (Secs. 1c–1e).
6. Delete stale `figs/thesis.png`, `figs/method.png`; add the FCC-talk number-mapping line to the circulation note.
No current figure was found to contain an obsolete number; the post-fix regeneration chain held up under this audit.


---

# Agent 6 — Referee risk

# Referee-Risk Audit — RADiCAL timing manuscript vs the prior record
Auditor: Agent 6 (read-only). Sources read: parent paper [PARENT] in full (15 pp), [INSTR] in full (6 pp), [TNS] in full (5 pp), [CPAD] and [FCC] talks sampled at the results/conclusions slides, plus `radical_timing.tex` (commit cb19f51), `REFEREE_RISK_MEMO_2026-06-09.md`, `COAUTHOR_CIRCULATION_NOTE_2026-06-09.md`, and the claims law.

## Key prior-record facts established
- **[PARENT]** (NIM A 1068 (2024) 169737): a = 255.58 ≈ 256 ps√GeV, b = 17.52 ≈ 17.5 ps, σ_t(150) = 27 ps, BestMinus = (Δt_DW − Δt_UP)/2 on **measured-energy bins 6–8** (itself a brightness selection). Timing extraction: *"fixed thresholds were set for the high-gain signals… a timestamp was set the moment a high-gain signal crossed a specific threshold. The level of this threshold was optimized"* (§5.2), and *"the timing algorithm, which uses a fixed threshold"* (§5.2.1 satellite paragraph). HG was *"driven further into saturation, to maximize the signal rise of the leading edge"*. The ~0.2 ns satellite (~15% area) appears **only in the MCP-referenced (DW+UP)/2** and is attributed to a one-sample-segment DRS shift relative to the MCP; the parent states explicitly *"It does not appear in Fig. 19 Lower Left"* — i.e., **absent from (DW−UP)/2**. Section 7 item 2 explicitly promises the DSB1-vs-LuAG:Ce comparison (the current paper); item 3 (a numbered list item, not a subsection "7.3") foresees doubling the LYSO thickness at shower max. The parent leaned on GEANT4 (Figs 7–9), including **Fig. 8: simulated σ_t vs detected light yield** — the exact scaling the new paper measures.
- **[TNS]** (FTBF June 2022, 28 GeV): MCP-referenced absolute σ(t_RAD − t_MCP) = 49.5±5 ps (E>15), 43±4.4 near peak, 42±4.3 after MCP subtraction. Different observable (absolute, not differential); no conflict.
- **[INSTR]**: G4 expectation 30 < σ_t < 50 ps; 4×4 array σ_E/E = 15.7%/√E ⊕ 0.1/E ⊕ 1% (full-module containment, distinct from E_SM); position goal "a few mm".
- **Talks**: [CPAD] (Nov 2023) and [FCC] (Mar 2024) show the *same* May-2023 CERN dataset, preliminary, with the 17.52 ⊕ 255.58/√E fit, "<30 ps at 150 GeV — meets needs of FCC EndCap", "RADiCAL has potential to reach <10 ps at >150 GeV", "the RADiCAL module is radiation hard", "demonstrated to meet the needs of high radiation… environments". [FCC] slide 39 already previewed "25 ps @ 150 GeV, with limiting resolution of ~18 ps" — **good** continuity with the current 25.7±0.6 / 18.8±0.8.

---

## (1) Top referee objections, ranked

### HIGH
**H1 — "This was already published" / same-dataset opacity (NEW; not in the 10-risk memo).**
The parent's 25–150 GeV CERN H2 dataset *is* the May-2023 DSB1 dataset of this paper. The manuscript says "DSB1, the brightest, is the build of Ref. [1]" (§2.1) but never states that the DSB1 **data** are the same as Ref. [1]'s, re-analyzed. A referee who realizes this independently will read silence as concealment — the worst version of this objection. (The defense material exists: Ref. [1] §7 promised this study; the published values are reproduced and shown in Table 1; the new content list in §1 is genuinely new — 3 new builds, unified chain, srCFD, MIXED control, systematics budget, depth drift.)

**H2 — "This contradicts your previous paper": method attribution of the published result (NEW).**
§5.3 reads "cfd05 (a CFD on the clipped measured peak, **as in Ref. [1]**)" and Table 1's published row lists src = "cfd05". The parent's published text describes a **per-channel optimized fixed threshold on the leading edge** — in the new taxonomy, LED-like, not a clipped-peak CFD. The parent's corresponding author (Perez-Lara, a coauthor here) and any careful referee will catch this. Either the production legacy chain genuinely is cfd05 and reproduces 27/17.5 (in which case say so and reconcile with the parent's text), or the attribution must be reworded. **Must be verified before circulation** — this also touches the satellite story (H/M-coupled, see M1).

**H3 — Claims-law self-violation: full-fiducial number missing (NEW; internal-compliance defect).**
The claims law (allowed-language block + forbidden #9) requires "27–29 ps over the full fiducial sample" to be quoted beside *any* brightest-1000 number. Nowhere in the manuscript (abstract, §5, conclusions) does a full-fiducial σ_t(150) appear. The external form of this objection: "your headline is the best 1000 events — what does the module do for everyone?" The K-scan (App. A, panel b) and the K=2000 systematics row partially defend, but the law's own fix is cheap and mandatory.

### MED
**M1 — Satellite provenance contradicts the parent (NEW).**
§4: "a small non-Gaussian shoulder (the ∼0.2 ns satellite **also noted in Ref. [1]**)" — in the context of the **(DW−UP)/2** distribution, attributed in §5.3/App. B to clipped-pulse cfd05 bias. The parent says its satellite lives **only** in the MCP-referenced sum and is a **digitizer sampling-segment** effect, explicitly absent from the differential. As written, the cross-reference invites "Ref. [1] says the opposite about your own estimator." The two observations may both be true (different estimator, different mechanism, different distribution), but the citation must be qualified or dropped.

**M2 — "Where are the simulations?" (under-defended; memo has no row).**
The parent built its case on G4 (Figs 7–9); this paper has none and uses Longo + X₀/v_eff parametrics. Currently the only acknowledgment is in the circulation note's limitations, not in the manuscript. Free win available: the parent's **Fig. 8 predicted σ_t ∝ (detected LY)^(−1/2)** — citing it converts "no simulation" into "our measurement confirms the programme's own simulation expectation."

**M3 — Estimator heterogeneity inside the four-build comparison (partially covered by memo risk 3, distinct from it).**
Table 1 mixes srCFD (DSB1, MIXED) with LED (LUAG, TENERGY); §3 says "applied uniformly to every build" and two sentences later "…LED is used instead" — internal tension a hostile reader will quote. The defense exists (App. A panel a: the adopted source lies *below* the whole clipped-peak family for every build, dramatically for LUAG — i.e., LED makes LuAG look *better*, so the >×2 stochastic gap is conservative) but the conservative-direction argument is never stated explicitly.

**M4 — carried-over memo risks that remain MED and adequately defended:** kinetics-confound framing (risk 1), MIXED ratio not unity (risk 2), floors not all 20 ps (risk 4), geometry-specificity of the same-shower control (risk 7). No change to their dispositions.

### LOW
**L1 — Talk continuity ("you changed your story").** "<10 ps potential at >150 GeV" (CPAD s31, FCC s40) vs the 18.8 ps confirmed floor; "the RADiCAL module is radiation hard" / "demonstrated to meet the needs of high-radiation environments" (both talks) vs the paper's citation-only hardness stance; "meets needs of FCC EndCap" vs the paper's no-requirements-met framing. All are conference-preliminary claims the paper deliberately narrows; the [FCC] talk's "25 ps / ~18 ps" preview actually *helps*. Prepare stock responses; no manuscript edits. (Stock line for <10 ps: the measured floor is shower-depth physics common to any scintillator in this geometry; the route below it is the depth-corrected estimator named in §8, plus the §7-item-3 light-recovery design; the 10 ps program target is unchanged.)
**L2 — Citation precision:** "its Section 7.3" (Discussion, line ~544) — the parent's §7 is a numbered list; item 3 is not subsection 7.3.
**L3 — Depth-dial wording nuance:** §7 "consistent in sign and magnitude" is a notch stronger than Tier-2's mandated "order-of-magnitude level", though the 30% excess and acceptance qualifier are present. Borderline-compliant; tighten if cheap.
**L4 — TNS absolute-vs-differential:** no conflict; taxonomy already in §4. If a referee juxtaposes 42 ps @ 28 GeV (absolute) with the differential scan, the answer is the two-estimator taxonomy.
**L5 — Hygiene:** author-roster placeholder, funding TODO (line 664), TODO-P2 money-plot comment (lines 292–297) — known circulation items.

## (2) Current defense per objection + adequacy
| # | Current defense | Adequate? |
|---|---|---|
| H1 | §1 frames the paper as the comparison promised by Ref. [1] §7; Table 1 carries the published row; §5.3 reproduces 26.9/27.4 ps under the legacy-style estimator on identical events | **Partial** — the same-dataset fact must be explicit; silence is the risk |
| H2 | None — the attribution is asserted | **Inadequate**; needs verification + rewording (edit E1) |
| H3 | App. A K-scan + K=2000 systematics row | **Inadequate vs the paper's own law**; add the 27–29 ps companion (edit E3) |
| M1 | None — "also noted in Ref. [1]" asserted without qualification | **Inadequate**; qualify or drop (edit E2) |
| M2 | Floor argument framed as consistency reading; limitations live only in the circulation note | **Thin**; one cross-cite sentence fixes it (edit E4) |
| M3 | App. A panel (a); per-regime rule of §5.3; "method choice… does not alter the four-build conclusions" | **Mostly adequate**; fix the "uniformly/instead" tension + state the conservative direction (edit E5) |
| M4 (memo 1,2,4,7) | As per memo | Adequate; unchanged |
| L1 | Paper is strictly narrower than the talks; FCC talk previewed the new headline | Adequate; stock responses only |
| L2–L5 | — | Trivial fixes |
Old-numbers objections ("27 vs 25.7", "17.5 vs 18.8") remain **well defended**: published row in Table 1, "confirms"-only verbs, identical-event reproduction of the legacy values, fit-range stability (+0.2 ps) in §6. The retired 0.99 MIXED ratio never appeared in any external talk or paper I examined (CPAD/FCC show no MIXED result at all), so the "changed your story" exposure on MIXED is **internal-only** and fully handled by the circulation note.

## (3) Exact text edits recommended
- **E1 (H2) — §5.3, line ~350.** After verifying against the parent production chain: if the parent used a fixed threshold (as its text says), change to: *"…computed with cfd05, a CFD referenced to the clipped measured peak --- the representative clipped-peak discriminator of this analysis; Ref.~\cite{radical} itself timed a per-channel fixed threshold on the saturated leading edge, whose behaviour the LED diagnostic reproduces ---"* and change Table 1's published-row "src" from `cfd05` to `thr.~\cite{radical}` (or `LED`). If the legacy chain *is* cfd05 and reproduces 27/17.5, instead add: *"(the cfd05 chain reproduces the published values on this dataset; we use it as the legacy reference)"* and reconcile the wording with the parent's fixed-threshold description in a footnote.
- **E2 (M1) — §4, line ~231.** Replace "(the $\sim\!0.2$~ns satellite also noted in Ref.~\cite{radical})" with: *"(Ref.~\cite{radical} noted a satellite of similar $\sim\!0.2$~ns offset in its MCP-referenced sum distribution and attributed it to a digitiser sampling effect; the shoulder discussed here appears in the differential under the clipped-peak discriminator and is suppressed by the recovered edge, Sec.~\ref{sec:method-gain})"*.
- **E3 (H3) — §5.3 end (line ~388) and Conclusions (line ~574), optionally abstract.** After "…is $25.7\pm0.6$~ps." add: *"Over the full fiducial sample, without the brightness selection, the same estimator gives $27$--$29$~ps at $150$~GeV."* (Use the exact audited full-fiducial value from the gated products.)
- **E4 (H1 + M2) — §2.1, line ~137.** Extend: *"DSB1, the brightest, is the build --- and the dataset --- of Ref.~\cite{radical}; its timing is re-derived here under the unified analysis chain so that all four builds are compared on an equal footing, and the published results are reproduced for reference (Table~\ref{tab:builds}, Sec.~\ref{sec:method-gain})."* And in §5.2 (after "the stochastic term tracks it", line ~278) or Discussion: *"This light-yield scaling is the behaviour anticipated by the programme's GEANT4 studies, which predicted the single-module time resolution to be set by the detected light level~\cite{radical,instruments2022}."*
- **E5 (M3) — §3, lines ~210–214.** Replace "The method is applied uniformly to every build. For the low-light builds…, a robust fixed-threshold… (LED) is used instead;" with: *"The recovery is attempted uniformly for every build; where the small pulses of the low-light builds (LuAG, TENERGY) defeat the fractional-threshold reconstruction, the per-regime rule of Sec.~\ref{sec:method-gain} selects the robust fixed-threshold LED time instead (Table~\ref{tab:builds})."* Optionally add in §5.2: *"This assignment is conservative for the cross-build comparison: the adopted source lies below the entire clipped-peak family for LUAG (Fig.~\ref{fig:opt}a), so the more-than-factor-of-two stochastic-term difference is not an artefact of the per-build estimator choice."*
- **E6 (L2) — Discussion, line ~544.** "(its Section~7.3)" → "(its Section~7, item~3)".

## (4) Memo and circulation-note updates
**REFEREE_RISK_MEMO_2026-06-09.md — add/modify rows:**
- ADD **Risk 11 (HIGH→LOW after E4):** "DSB1 scan re-analyzes the published dataset — say so explicitly" (defense: §1 promised-study framing + Table 1 published row + §5.3 reproduction; blocker until E4 lands).
- ADD **Risk 12 (HIGH, verification gate):** "Parent method attribution — manuscript calls Ref. [1]'s estimator cfd05; Ref. [1]'s text says per-channel fixed threshold on the leading edge. Verify the legacy production chain; apply E1; coauthor question to C. Perez-Lara."
- ADD **Risk 13 (MED):** "Satellite cross-reference contradicts parent (sum-only, digitizer mechanism vs differential, clipping mechanism) — qualify per E2."
- ADD **Risk 14 (MED→LOW after E4):** "No simulations — cite parent Fig. 8 G4 LY-scaling as supporting expectation; G4 campaign honestly deferred."
- ADD **Risk 15 (MED→LOW after E5):** "Per-build estimator heterogeneity in tab:builds — make the conservative-direction argument explicit; fix the 'uniformly/instead' sentence."
- ADD **Risk 16 (compliance, blocker):** "Full-fiducial companion number absent beside every brightest-1000 quote — violates the wording law's forbidden #9; fix via E3 before circulation."
- UPDATE **Risk 9:** the four citation placeholders appear resolved in the current bibliography (verify once); remaining hygiene = author roster, funding text (line 664), TODO-P2 figure restyling comment.
- UPDATE **Risk 1:** the Gundacker citation is now complete; drop the placeholder caveat.
- ADD a **talk-continuity appendix** (LOW): stock responses for the CPAD/FCC "<10 ps potential", "radiation hard demonstrated", and "meets FCC EndCap needs" claims; note the FCC 2024 "25 ps / ~18 ps" preview as favorable continuity.

**COAUTHOR_CIRCULATION_NOTE_2026-06-09.md:**
- ADD to "Specific questions for coauthors" (as Q7): *"Ref. [1] coauthors (esp. C. Perez-Lara): the manuscript currently labels the published timing method 'cfd05'; Ref. [1]'s text describes a per-channel fixed threshold on the HG leading edge. Which estimator produced the published 27 ps / 17.5 ps, and does our cfd05 row correctly represent it?"*
- ADD to "What changed scientifically" or "Known limitations": one line stating the DSB1 sample is the Ref. [1] dataset re-analyzed under the unified post-fix chain (mirrors edit E4), so no coauthor is surprised by the same-data relationship.
- ADD to "Known limitations": the satellite-attribution difference vs Ref. [1] (sum/digitizer vs differential/clipping) and its planned qualification (E2).
- ADD to the pre-circulation checklist: the full-fiducial companion number (E3) as a wording-law blocker alongside the existing citation-verification task.

**Bottom line:** the existing memo's risks 1–8 hold up against the actual prior record. The genuinely new exposures are continuity facts a parent-paper coauthor or diligent referee will find in an afternoon: same dataset (H1), method attribution (H2), satellite provenance (M1) — plus one self-inflicted claims-law gap (H3). All four are fixable with the six edits above; none threatens the thesis, the MIXED control, or the headline numbers.

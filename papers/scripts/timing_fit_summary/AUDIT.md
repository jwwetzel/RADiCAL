# Timing Fit Summary Audit (pre-analysis) — post-fix four-build (a, b) recompute

Date: 2026-06-09. Written BEFORE the recompute script.

1. **Builds:** DSB1, LUAG, MIXED, TENERGY (all four 2023 builds with reduced data).
   MIXED is included as a module-wide single-estimator row WITH AN EXPLICIT CAVEAT: a single
   module-wide method is ill-posed for MIXED (half the corners clip, half do not — established this
   campaign); its row exists for completeness, and the paper treats MIXED per material (GATE 6).
2. **Primary estimator per build (the adopted-source rule, stated ex ante):** clipped/bright builds
   → srCFD (`hg_lgcfd`): DSB1, MIXED; dim/unclipped builds → LED (`hg_led`): LUAG, TENERGY.
3. **Energies:** DSB1 25–150 GeV (6 points); LUAG/MIXED/TENERGY 50–150 GeV (5 points). Fit ranges
   restricted to each build's measured span.
4. **Cuts (production):** `wc_ok`; MCP1 200–750 mV; TimingFiducialR(E) (2.5 mm ≤100 GeV → 3.0 mm
   ≥125 GeV, continuous ramp between — the post-fix continuous fiducial); per-end validity
   hg_peak ≥ 20 mV; the in-event channel-consistency veto (|t − median| < 2 ns, `eventDWUP`);
   brightest-1000 selection by sum_lg (the locked headline selection).
5. **Functional form:** σ_t(E) = √( (a/√E)² + b² ), fit over the measured span only.
6. **Units:** a in ps·√GeV; b in ps; E in GeV.
7. **Historical vs recomputed:** the circulating numbers (DSB1 a≈201 b=19.5±1.1; LuAG a≈455
   b=19.8±5.6; MIXED 233/34.6; TENERGY 240/22±3) are PRE-FIX (old tebSigma window, no veto) except
   the partial post-fix regression of 2026-06-09 (DSB1 lgcfd b=18.8±0.8, σ(150)=25.7). THIS RUN
   recomputes ALL four builds with the production post-fix estimator (`RadTiming.h` robust window +
   tail guard + `eventDWUP` veto) and becomes the single authoritative table.
8. **Uncertainties:** the fit uses per-point statistical errors σ/√(2K), K=1000; scatter/systematic
   beyond statistics is folded by inflating parameter errors by √(χ²/ndf) when χ²/ndf > 1 (PDG
   scale-factor convention). Selection systematics are NOT re-derived here (they live in
   `papers/timing/tab_systematics.tex`, ≤1–1.5 ps on σ(150)).

**Narrative checks to perform (pass/fail printed by the script):** a_DSB1 ≈ 200; a_LuAG ≈ 455;
a-ratio ≈ 2.3; b shared ≈ 20 ps among LYSO builds; best measured σ(150) ≈ 25–26 ps (DSB1, srCFD,
brightest-1000); floor verb remains "confirms 17.5 ps". If recomputed values differ, THE NARRATIVE
IS UPDATED to the new numbers — not the reverse.

Outputs: `papers/tables/timing_fit_summary_2026-06-09.md`,
`papers/figures/timing_fit_summary/timing_fit_summary.png`, log in this directory.

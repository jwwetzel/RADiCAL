# CLAUDE.md — session bootstrap for the RADiCAL 2023 analysis

This file auto-loads for any Claude Code session in this repo, on any machine. It replaces
nothing — it points. The repo is deliberately self-contained (workspace audit, Phase A).

## Read first, in order

1. `README.md` — headline + trust path.
2. `ANALYSIS_GUIDE.md` — the waterfall + full provenance table (every figure/number → generator
   → committed evidence).
3. `papers/CAMPAIGN_SNAPSHOT_2026-06-09.md` — append-only campaign record; the LAST update is
   the current state.
4. `papers/memory_claims_and_forbidden_language.md` — the claims law. **Authoritative over any
   other phrasing, including older docs/talks.** Companion: `papers/memory_analysis_gates.md`.

## Standing rules (non-negotiable, learned the hard way)

- **Claims law governs every public number.** Headline: σ_t = 25.7 ± 0.6 ps (brightest-1000,
  150 GeV, DSB1, srCFD), ≈50 ps full fiducial beside it; floor 18.8 ± 0.8 **confirms** (never
  revises) the published 17.5 ps; LYSO:Ce + W are common to all builds — **the WLS capillary is
  the variable** (never "crystal differs").
- **Gate discipline:** published results live in `papers/scripts/<gate>/` (macro + AUDIT.md +
  committed log). Gate macros REGENERATE their committed outputs in place — run them from the
  **repo root**, to verify only, and treat any dirty `git diff` as a finding. Expected
  reproduction residual: `thesis_postfix.pdf` /CreationDate only. Never edit committed logs or
  generated tables by hand; never backdate an audit.
- **Snapshot is append-only** — record each work pass as a new UPDATE with commands + hashes.
- **The committed PDF is always the tectonic build** (`cd papers/timing && tectonic
  radical_timing.tex`); Texifier may be used for editing (class + orcidlink vendored).
- **Never commit ROOT data files**; `papers/**/*.log` is gitignored — result logs need
  `git add -f`. Figures >2000 px: `sips -Z` a /tmp copy before viewing.
- **Memory duplication rule:** durable project facts belong IN-REPO (`docs/`,
  `data/2023/metadata/`, `papers/memory_*.md`); the local `~/.claude` memory holds working
  context only. When a fact changes, update the repo copy.
- The repo is **public** (github.com/jwwetzel/RADiCAL) and the whole workspace publishes with
  the paper — anything committed here is world-readable; assistant memory must never be
  committed into it.

## Environment

`source setup.sh` once per shell (ROOT include path + RAD_DATA). Reduced ntuples (~13 GB) live
in `data/2023/reduced/` (Dropbox-synced, not in git); raw mirror on Argon
(`tools/pull_argon_data.sh`) / local subset in `data/2023/raw/`.

## Where things stand / what's next

See the last CAMPAIGN_SNAPSHOT update, plus `WORKSPACE_REORG_PLAN_2026-06-10.md` §4–5
(Phase B pre-submission checklist: LICENSE, CITATION, data DOI, repro script) and
`CODE_AUDIT_2026-07-21.md` (deep items). Open human items: GATE-1 logbook scan, coauthor
Q7/Q8, NIM_A_Figures rights decision.

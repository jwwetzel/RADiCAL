# Manuscript figure manifest — papers/timing/figs/

Recorded 2026-07-21 (workspace audit, Phase A6). Every file the manuscript
includes, its md5 (first 12 hex) AT THE TIME OF RECORDING, its generator, and how the copy
got here. If a regeneration changes a checksum, update this manifest in the same commit —
a mismatch between this table and the tree means drift and must be treated as a finding.
Full provenance chain: ../../../ANALYSIS_GUIDE.md.

| File | md5(12) | Generator | Copy method |
|---|---|---|---|
| `clip.png` | `c803f0a28f6a` | analyze/studies/narrativeFigs.C (fig 1) | cp figures/2023/narrative/narrative_clip.png |
| `depth_dial.png` | `50d0842f8297` | papers/scripts/depth_dial/depthDial.C | cp papers/figures/depth_dial/depth_dial.png |
| `dist.png` | `2accb5798139` | FROZEN AS CIRCULATED: legacy distribution study; top title strip manually cropped 2026-06-09 (FORMAT_AUDIT F7) — regeneration deferred to protect content identity | manual (documented) |
| `hglg.png` | `089aa487c817` | analyze/studies/hgLgPlot.C("DSB1",RUN1075,100) | cp figures/2023/hglg_panel_DSB1.png |
| `method_postfix.png` | `a549e8037ab4` | papers/scripts/method_gain_postfix/methodGainPostfix.C | cp papers/figures/method_gain_postfix/method_gain_postfix.png |
| `mixed_h2h_corrected.png` | `569e9b0d09e9` | papers/scripts/mixed_killshot_bootstrap/makeMixedKillshotFigure.C | cp papers/figures/mixed_killshot_bootstrap/mixed_h2h_corrected.png |
| `optimization.png` | `7c6484709243` | legacy selection-scan macro (optAppendix.C era); working points verified post-fix, curve re-verification is a pre-submission item | historical copy |
| `pulse.png` | `30dc775a5920` | analyze/studies/narrativeFigs.C (fig 2) | cp figures/2023/narrative/narrative_pulse.png |
| `radical_apparatus_composite.pdf` | `0f6ae9469b94` | papers/scripts/apparatus_composite/apparatus_composite.tex (tectonic) | direct output copy (cp from generator dir) |
| `radical_apparatus_composite.png` | `553204c1f95a` | same, rasterized | sips -s format png --resampleWidth 3200 |
| `satellite_removal.png` | `ef08efddc685` | papers/scripts/satellite_removal/satelliteRemoval.C | cp papers/figures/satellite_removal/satellite_removal.png |
| `systematics.png` | `83a033e33911` | papers/scripts/systematics_postfix/systematicsPostfix.C (post-fix regeneration 2026-06-10; annotations = Table 2 arrays by construction) | written directly by the gate script |
| `thesis_postfix.pdf` | `8db4962d0889` | same | written directly by the gate script |
| `thesis_postfix.png` | `6d457c2a317d` | papers/scripts/timing_fit_summary/timingFitSummary.C | written directly by the gate script |

Notes: `dist.png` and `optimization.png` are the two files without a live gate generator —
both are frozen-as-circulated with their status recorded above and in
`FORMAT_AUDIT_2026-06-09.md` / `WORKSPACE_REORG_PLAN_2026-06-10.md` (Phase B).

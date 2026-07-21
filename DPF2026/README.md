# RADiCAL @ DPF 2026 — slide deck

HTML slide deck for **"RADiCAL — Research and Development of Optical Materials for
Ultra-compact, Fast-Timing EM Calorimetry"**, J. Wetzel (Iowa) for the RADiCAL
Collaboration. DPF 2026, Fermilab, Instrumentation (Calorimetry & Scintillators).
14 slides + speaker notes, ~15 minutes.

## Files

| File | Use |
|---|---|
| `index.html` | The deck. Edit this. References `assets/` (keep the folder together). |
| `standalone.html` | **Travel copy** — one self-contained file, images baked in, works offline anywhere. Email it / USB it. Rebuild after edits: `python3 build_standalone.py`. |
| `assets/` | Figures (the 8 used + a few spares for backup slides). |
| `SPEAKER_NOTES.md` | Per-slide talking points + a 12-minute timing budget. |

## Present

Open `index.html` (or `standalone.html`) in any browser, press **F** for fullscreen.

| Key | Action |
|---|---|
| **→ / Space / click-right** | next slide |
| **← / click-left third** | previous |
| **F** | fullscreen |
| **S** | toggle speaker notes (per-slide) |
| **O** | jump to slide number |
| **Home / End** | first / last |

Also: swipe on a touchpad/iPad, and `#7` in the URL deep-links to slide 7.

## Export to PDF (backup / handout)

Open in **Chrome** → Print → **Landscape**, margins **None**, "Background graphics" **on**,
scale to fit → Save as PDF. Each slide becomes one landscape page. (Have this PDF on a USB
as a projector fallback.)

## Rebuild the standalone after editing

```bash
python3 build_standalone.py     # index.html + assets/  ->  standalone.html
```

## Story arc (14 slides)

1. Title · 2. Why fast-timing calorimetry (targets: 10 ps / 1 mm / 10%/√E) ·
3. The module (apparatus) · 4. Four builds — WLS capillary is the variable ·
5. srCFD saturation-recovery method · 6. **25.7 ps** headline · 7. Light yield, not
species (MIXED same-shower) · 8. The depth dial · 9. Preliminary energy & position ·
10. Repositioning: a calorimeter **test bed** · 11. The three knobs · 12. Now: LuO:Yb +
electronics · 13. Roadmap → CERN T10 (Aug 2026) · 14. Summary.

## Provenance

Every figure comes from the committed analysis chain (see `../ANALYSIS_GUIDE.md`):
money plot, method gain, MIXED, depth dial from `papers/scripts/*/`; energy/position from
the companion analysis; apparatus composite from `papers/scripts/apparatus_composite/`.
All quoted numbers obey `../papers/memory_claims_and_forbidden_language.md` (the claims law):
timing 25.7 ps brightest / ≈50 ps full-fiducial with the floor **confirming** (not revising)
the published 17.5 ps; energy/position labeled **preliminary**; radiation tolerance stated at
component level; no 4D/5D or calibrated-z claims. 2023 tiles were LYSO:Ce; **LuO:Yb** is the
new ceramic scintillator for the T10 run.

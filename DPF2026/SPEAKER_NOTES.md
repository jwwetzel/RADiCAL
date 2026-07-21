# Speaker notes — RADiCAL @ DPF 2026

~12 min talk + ~3 min questions. The same notes are embedded live in the deck — press **S**.
Target ~50 s/slide; the three "landing" slides (6, 7, 8) deserve a beat longer.

| min | slide | beat |
|----|-------|------|
| 0:00 | 1 Title | who / what: R&D on the *optical materials* for fast-timing calorimetry. |
| 0:40 | 2 Motivation | why timing; the three program targets (10 ps / 1 mm / 10%/√E). |
| 1:30 | 3 Module | index-finger W/LYSO shashlik; WLS at shower max; **(DW−UP)/2 needs no clock**. |
| 2:40 | 4 Builds | the experiment's design: **LYSO+W fixed, WLS capillary is the one variable**. The question: light or kinetics? |
| 3:40 | 5 srCFD | HG clips at 820 mV → recover the true peak from the linear LG → time the steep edge (×3.7). A method contribution. |
| 4:45 | 6 **Headline** | **25.7 ± 0.6 ps** brightest / ≈50 ps full fiducial (quote BOTH); floor 18.8 **confirms** published 17.5. |
| 6:00 | 7 **Thesis** | a: 203→440 (>2×) AND same-shower 1.04±0.05 ⇒ **light, not species** — within this geometry & light levels. |
| 7:20 | 8 **Depth dial** | −33.6 ps/e-fold, consistent with shower-max migration. NOT a calibrated z. ⇒ *depth is a knob.* |
| 8:30 | 9 Energy+position | **preliminary**: σ_E/E ≈14% (species-independent); position ~1.5 mm (upper bound), ~0.9 mm core. Toward (x,y,E,t). |
| 9:30 | 10 Pivot | light dominates + depth is a knob ⇒ build a **test bed**, not one more fixed prototype. |
| 10:10 | 11 Knobs | tile (LYSO→**LuO:Yb**) · WLS species **+ depth** · readout (SiPM adapter cards). Swap each independently. |
| 11:00 | 12 Now | LuO:Yb tiles **in hand, drilling**; new cards; adjustable depth mechanics. |
| 11:40 | 13 Roadmap | 2022 → 2023 → now → **CERN T10, end of August** → targets. |
| 12:10 | 14 Summary | 25.7 ps; light≫species; every knob. LuO:Yb drilling → T10. Thanks / questions. |

## Careful-wording reminders (claims law)

- Always pair the headline with the full-fiducial companion: "25.7 ps for the brightest showers, ≈50 ps full fiducial."
- Floor **confirms** the published 17.5 ps — never "improves" or "revises."
- "Light yield, not WLS species, sets the timing — **within this geometry and at these light levels.**"
- Depth dial is an **ensemble consistency** reading, **not** a per-event / calibrated z.
- Energy & position are **preliminary** (companion paper in prep); position numbers are **upper bounds** (tracker-limited).
- Radiation tolerance is **component-level** (LuAG:Ce, W by citation) — not demonstrated in this beam data.
- "Toward simultaneous (x, y, E_SM, t)" — a demonstration goal, not an achieved 4D/5D calorimeter.
- 2023 scintillator tiles were **LYSO:Ce**; **LuO:Yb** is the *new* ceramic scintillator for T10.

## Likely questions

- *Why 25.7 and not the program's 10 ps?* 10 ps is the target for the optimized detector; 25.7 ps is this
  prototype at 150 GeV. The test bed is how we close the gap — more light, tuned depth.
- *Is LuAG really as good as DSB1?* Same-shower widths agree at 10–20% (1.04±0.05); the cross-build gap is the
  light each collects, not the shifter kinetics — in this geometry.
- *Radiation hardness?* Component-level for LuAG:Ce/W by citation; irradiated-module beam validation is future work.
- *Why adjustable capillary depth?* The depth dial: timing sensitivity tracks shower-max depth, which moves with energy.

# site/ — the published RADiCAL website

Self-contained GitHub Pages site. Everything it references lives **inside this
folder**, so it deploys as one unit.

```
site/
├── index.html            redirect → going_radical.html (the landing)
├── going_radical.html    the field guide (main hub)
├── paper.html            the paper
├── calor2026_radical.html / calor2026_slides.html   CALOR-2026 talk + slides
├── data_lab.html / geant4_lab.html                  student lab exercises
├── report/               generated analysis report (built by analyze/makeReport.py,
│                         published here by analyze/deployReport.sh)
├── assets/               every image the pages embed (capillary plots + lab figs)
└── .nojekyll             serve raw (no Jekyll)
```

## How it is served

Published via GitHub Actions (`.github/workflows/pages.yml`), which uploads this
`site/` folder as the Pages artifact. GitHub maps the artifact to the **site
root**, so `site/going_radical.html` is served at `<pages-url>/going_radical.html`
— i.e. moving the pages here did **not** change any public URL.

**One-time setup** (repo Settings → Pages → Build and deployment):
set **Source = "GitHub Actions"** (previously "Deploy from a branch / root").
After that, any push touching `site/**` auto-deploys.

## Editing

- The HTML pages are hand-authored — edit them here directly.
- `report/` is generated: `python3 analyze/makeReport.py` then
  `bash analyze/deployReport.sh` (writes into `site/report/`).
- `assets/` holds published copies of embedded plots. The capillary plots live
  here; the four lab figures are copies of `student/figs/expected_*.png`.

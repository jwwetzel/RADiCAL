# legacy/ — superseded early analysis (do NOT use)

This folder holds the original, monolithic first-pass analysis. It has been
**fully superseded** by the modular pipeline in [`../Analysis/`](../Analysis),
which is the current, maintained codebase. These files are kept only for
provenance and are not referenced by anything in `Analysis/`.

| File | What it was |
|------|-------------|
| `analyzeRad.C`, `analyzeRadRDF.C` | original single-file analysis macros |
| `rootlogon.C` | early ROOT style logon (now `Analysis/RADiCALStyle.h`) |
| `Output_old/`, `*.root`, `*.pdf` | outputs from the old workflow (gitignored) |

➡  **To run the analysis, see the top-level [`README.md`](../README.md) Quick Start.**

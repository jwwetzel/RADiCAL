# RADiCAL collaboration / authorship model

Authorship is **time-varying** — people join and leave, affiliations change, and not every
member authors every paper. So the collaboration record is organized **by test-beam era**,
matching the `RAD_YEAR` era-container model the rest of the workspace uses.

## Where things live

| File | Role | Edited by |
|------|------|-----------|
| `data/<era>/metadata/authors.yaml` | **Source of truth** — the ordered author list, affiliations, and roles for that campaign's paper(s). | **Hand** (with coauthor sign-off) |
| `data/_TEMPLATE/metadata/authors.yaml` | Schema + template for a new era. | Hand (copy it) |
| `collaboration/people.yaml` | Aggregate roster across all eras; `years_active` is **computed**. | **Generated** — do not edit |
| `data/<era>/metadata/generated/authors.tex` | elsarticle `\author`+`\affiliation` block, a drop-in `\input`. | **Generated** — do not edit |
| `data/<era>/metadata/generated/zenodo_creators.json` | `creators` for that era's Zenodo/data-DOI deposit. | **Generated** — do not edit |
| `CITATION.cff` (repo root) | GitHub "Cite this repository". | **Generated** — do not edit |

## The one rule: ORCID is the identity join key

The same person keeps the **same ORCID** across every era they appear in. That is what lets the
aggregate roster and each person's `years_active` be *computed* rather than maintained by hand —
so there is no second list to drift. The only fields that repeat across era files are the stable
ones (name, ORCID); the things that legitimately change per era — **affiliation, role, and whether
the person is on that paper** — live in the era file, where they belong.

`tools/gen_authors.py` cross-checks this: if one ORCID carries two different names across eras it
**errors** (free typo / disambiguation detection). Authors without an ORCID are listed as a warning
so the gaps can be filled over time.

## Workflow

Add or update a campaign:

1. Copy `data/_TEMPLATE/metadata/authors.yaml` to `data/<era>/metadata/authors.yaml` (or edit the
   existing one). **Reuse each returning person's exact ORCID from 2023** so the join works.
2. Regenerate every artifact from the repo root:
   ```bash
   python3 tools/gen_authors.py
   ```
   (`--check` validates without writing; `--era 2023` regenerates one era; `--citation-era 2023`
   picks which era populates the root `CITATION.cff`.)

The generated files are **committed and regenerated in place**, exactly like the gate outputs:
after editing a source `authors.yaml` and running the tool, a dirty `git diff` on the generated
files is expected; a dirty diff *without* a source change is a finding (someone hand-edited a
generated file).

## Notes

- **Governance.** The author list and collaboration membership are collaboration-governed. These
  files should reconcile with the collaboration's official roster (see the open coauthor items,
  `papers/timing/COAUTHOR_CIRCULATION_NOTE_2026-06-09.md`); they are a structured mirror, not a new
  authority.
- **Compound surnames.** The Zenodo/CITATION "Family, Given" split defaults to *last token =
  family*. For compound surnames (e.g. "Sunar Cerci", "Karasu Uysal") set explicit `family:` /
  `given:` in the era file — the tool warns on any multi-token name lacking them.
- **The manuscript.** `data/2023/metadata/generated/authors.tex` is a drop-in `\input` for
  `papers/timing/radical_timing.tex`; the manuscript currently carries the block inline (identical
  roster). Wire in the `\input` when convenient — the generated block uses appearance-order
  affiliation letters, cosmetically different from the hand-written one but equivalent.
- **Zero dependencies.** `gen_authors.py` is stdlib-only; it uses PyYAML if importable, otherwise a
  vendored strict-subset YAML loader that covers the `authors.yaml` grammar.

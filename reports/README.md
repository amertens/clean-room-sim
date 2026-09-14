# reports/

The written outputs of the cleanTMLE project. One canonical source per
document; rendered copies sit next to their `.qmd`. Superseded drafts
(the JSS and methods-journal forks of the manuscript, the applied-analyst
tutorial, the May review notes, and the standalone talk-figure sets) were
removed in September 2026 and live in git history if ever needed.

## Current documents

- `manuscript_outcome_blind_dq.qmd` (+ `.html`, `.docx`): the manuscript.
  Two-failures framing, the estimability methods section, and the
  multi-arm Rescue.Co case study, with every case-study table computed
  from `rescueCo/results/multiarm/*.csv` at render time.
- `cleanTMLE_presentation.qmd` (+ `.pptx`): the talk. Same result files,
  five data-drawn figures; renders with the `templates/lab_template.pptx`
  reference doc.
- `references.bib`: the shared bibliography.

## Companion analyses the manuscript cites

- `workflow_contrast_c.qmd` (+ `.html`): does plasmode candidate
  selection beat a fixed-library TMLE across DGPs (the workflow-contrast-C
  question)? Source updated to the 0.2.0 API in September 2026; the
  rendered `.html` is the May run under the pre-consolidation calls,
  which fit the same estimators. Re-render to refresh.
- `cleanTMLE_vs_causalRisk_actg.qmd` (+ `.html`, `.docx`): side-by-side
  reproduction of causalRisk's ACTG 320 analysis with the cumulative-risk
  grammar. Grammar verbs only, current API.

## Support files

- `make_figures.R` -> `figures/`: the two base-R appendix flowcharts
  (roadmap, governance) used where mermaid cannot render (docx).
- `templates/`: the pptx reference document.

The plain-language tutorial formerly at
`cleanTMLE_for_applied_analysts.qmd` was consolidated into the package
vignettes: the staged-analysis vignette's first half is the applied
walkthrough, and the function reference carries the dictionary.

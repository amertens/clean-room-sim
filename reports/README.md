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
  from `rescueCo/results/multiarm/*.csv` at render time. The prose and
  code listings describe cleanTMLE 0.3.0 (the sixteen-verb surface); the
  simulation sections read `results_new/`, which the revision protocol
  is regenerating under 0.3.0, so those sections render only once the
  reruns exist (see `revision/CHANGELOG.md`).
- `cleanTMLE_presentation.qmd` (+ `.pptx`): the talk. Same result files,
  data-drawn figures; renders with the `templates/lab_template.pptx`
  reference doc.
- `references.bib`: the shared bibliography, also used by the tutorial in
  `tutorials/`.

## Companion analyses the manuscript cites

- `workflow_contrast_c.qmd` (+ `.html`): does plasmode candidate
  selection beat a fixed-library TMLE across DGPs (the workflow-contrast-C
  question)? The rendered `.html` is the May run under the
  pre-consolidation calls, which fit the same estimators; the source
  predates the 0.3.0 verbs and is on the rerun list.
- `cleanTMLE_vs_causalRisk_actg.qmd` (+ `.html`, `.docx`): side-by-side
  reproduction of causalRisk's ACTG 320 analysis with the cumulative-risk
  grammar. Grammar verbs only, current API.

## Support files

- `make_figures.R` -> `figures/`: the two base-R appendix flowcharts
  (roadmap, governance) used where mermaid cannot render (docx).
- `templates/`: the pptx reference document.

The plain-language walkthrough that used to live here as
`cleanTMLE_for_applied_analysts.qmd` now has two homes: the package's
Get-started vignette (`vignette("cleanTMLE")`) and Full workflow article
cover the software surface, and `tutorials/clean_room_targeted_learning_tutorial.qmd`
covers the reasoning on a simulated cohort with a known truth.

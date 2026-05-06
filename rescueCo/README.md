# rescueCo/ — RescueCo Kenya Trauma Case Study

A **self-contained** clean-room TMLE workflow built on the
[cleanTMLE](https://github.com/amertens/cleanTMLE) package. This folder
can be lifted into any project (or made the root of a new project) and
the analysis will reproduce without any external paths.

For the original purpose / methods description see
[`README_legacy.md`](README_legacy.md).

## What's in this folder

```
rescueCo/
├── README.md                         (this file)
├── README_legacy.md                  original detailed README
├── config/clean_room_config.yml      analysis configuration
├── data/                             VENDORED Stata .dta files
│   └── Final Data Cut 03.27.2026/
├── R/                                helper functions sourced by scripts
│   ├── utils.R, ps_matching.R, plotting.R, simulation.R,
│   ├── negative_controls.R, tmle_binary.R, tmle_survival.R
├── scripts/                          numbered run scripts (00..06)
├── results/                          pipeline outputs (CSV / RDS / PNG)
│   └── manuscript_artifacts/         tables/figures for the methods paper
├── reports/                          rendered HTML reports
└── logs/                             pipeline run log
```

## Reproducibility

```r
# From the project root containing rescueCo/, OR from inside rescueCo/:
Rscript rescueCo/scripts/00_run_all.R
```

The pipeline auto-detects whether it's running from the project root or
from inside `rescueCo/` and adjusts paths accordingly. Everything it
reads (input data, helper R files) and everything it writes (results,
reports, logs) lives **inside `rescueCo/`** — no dependencies on
sibling folders.

## Required R packages

`cleanTMLE` (≥ 0.1.0, ideally 0.1.1+), `tmle`, `survtmle`,
`SuperLearner`, `glmnet`, `MatchIt`, `caret`, `haven`, `data.table`,
`janitor`, `ggplot2`, `yaml`, `survival`, `MASS` (for ordinal GOSE),
`jsonlite`, `here`, `kableExtra`.

Install cleanTMLE from source if not on CRAN:

```r
remotes::install_github("amertens/cleanTMLE")
```

## Stage map

| Stage | Script | What it does |
|-------|--------|--------------|
| 1 | `01_stage1_build_cohort.R` | Build cohort, compute GOSE, lock spec, register negative controls |
| 1b | `01b_stage1b_eda.qmd` | Optional EDA render (skipped if quarto not on PATH) |
| 2 | `02_stage2_design_diagnostics.R` | PS estimation, balance, IPTW, matching diagnostics |
| 2b | `02b_stage2b_negative_controls.R` | Negative-control TMLE (5 NCs); residual-bias gate |
| 3 | `03_stage3_outcome_blind_simulation.R` | Plasmode + DQ stress, candidate selection, gate |
| 4 | `04_stage4_binary_outcome_analysis.R` | Locked GOSE > 4 analysis + ordinal GOSE PO sensitivity |
| 5 | `05_stage5_survival_outcome_analysis.R` | survtmle (primary), Cox, discrete-time TMLE (sensitivity) |
| 6 | `06_stage6_render_report.R` | Render `rescueCo/reports/rescueco_clean_room_case_study.qmd` |

## Manuscript artifacts

`rescueCo/results/manuscript_artifacts/` holds the tables and figures
the cleanTMLE methods paper §8 case-study cites:

- `table_effects.csv` — one row per estimator
- `table_dq_coverage.csv` — DQ stress curves
- `fig_love_plot.png`, `fig_dq_gradient.png`, `fig_ic_plot.png`,
  `fig_clever_covariate.png`
- `audit_summary.csv`, `decision_summary.csv` — audit trail and decision
  log exported via `cleanTMLE::export_audit_trail()` /
  `export_decision_log()`
- `case_study_metadata.json` — package version, lock hash, primary
  effect estimate

Read in the methods Quarto via:

```r
meta <- jsonlite::fromJSON(file.path(rescueco_path,
  "rescueCo/results/manuscript_artifacts/case_study_metadata.json"))
te   <- read.csv(file.path(rescueco_path,
  "rescueCo/results/manuscript_artifacts/table_effects.csv"))
```

## Lifting this folder to a new project

```bash
cp -r rescueCo/ /path/to/new/project/
cd /path/to/new/project
Rscript rescueCo/scripts/00_run_all.R
```

Everything runs out of the box. The only dependency that must be
installed externally is the `cleanTMLE` R package (and its CRAN
dependencies).

## Negative-control reconciliation

The reconciliation memo for the disagreement between the legacy NC
analysis (significant for SES NCs) and the cleanTMLE NC TMLE
(non-significant) is at `../reports/nc_reconciliation_memo.md`. Short
version: SES NCs flag in the **full ambulance cohort** (which includes
interfacility transfers) but pass in the **primary cohort** (transfers
excluded). Excluding transfers is the right primary-cohort decision;
the cleanTMLE result is correct for that cohort.

Three of the five NCs (`chronic_diabetes_insulin`, `chronic_hiv_art`,
`fuel_wood`) are now re-attached to the lock data after the NZV filter
so the cleanTMLE NC TMLE evaluates all five rather than only the two
that survived NZV. See `01_stage1_build_cohort.R` after the W-matrix
construction.

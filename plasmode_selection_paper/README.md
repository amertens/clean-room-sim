# Plasmode-selection paper

A paper-quality simulation study comparing **fixed-library TMLE** (the traditional approach) against **plasmode-selected-library TMLE** (the cleanTMLE workflow) across five data-generating processes that span the kinds of outcome-model mis-specification an applied analyst is likely to face.

This folder is self-contained and is **not** required to build the main `cleanTMLE` manuscript or the package. The deeper question it answers is whether plasmode-based candidate selection improves TMLE operating characteristics relative to a fixed-library default — and if so, under which DGPs.

## Folder layout

```
plasmode_selection_paper/
├── R/
│   ├── dgps.R              # five data-generating processes
│   ├── candidates.R        # the cleanTMLE candidate grid (6 candidates)
│   ├── workflows.R         # fixed-library and plasmode-selected fitters
│   └── run_simulation.R    # the driver function
├── scripts/
│   ├── run_pilot.R         # ~20-30 min smoke-test budget
│   └── run_full.R          # paper-quality budget (~6-24 h)
├── reports/
│   ├── manuscript.qmd      # draft paper with [PLACEHOLDER] tags
│   └── references.bib
├── results/                # generated output (gitignored)
└── README.md (this file)
```

## How to run

```r
# Pilot (fast, for development and to populate the manuscript draft)
setwd("plasmode_selection_paper")
source("scripts/run_pilot.R")

# Full study (paper-quality replicate budget)
source("scripts/run_full.R")

# Render the manuscript draft. It auto-loads the full results if present,
# otherwise falls back to the pilot results, and prints a warning banner.
quarto::quarto_render("reports/manuscript.qmd", output_format = "html")
```

### Interruption-resistant runs

`run_simulation()` writes incremental checkpoints to its `out_path` every `checkpoint_every` replicates (default 5) **and** at the end of each DGP. Writes are atomic: results are saved to `<out_path>.tmp` and renamed, so the final rds file is never half-written even if the process is killed during a save.

If `out_path` already exists when `run_simulation()` starts and `resume = TRUE` (the default), the function reads the existing file, identifies all `(dgp, rep)` cells whose three workflow rows are already present without NA estimates, and skips them. So you can:

```r
# Start the full study
source("scripts/run_full.R")

# Kill it (Ctrl-C / close R / power outage)
# ... later, in a fresh R session ...
source("scripts/run_full.R")   # picks up where it stopped
```

Each replicate uses a deterministic seed derived from `seed_base + 1000 * dgp_index + rep_i`, so resumed runs reproduce the same sequence of datasets and the same Monte Carlo numbers as an uninterrupted run.

To force a fresh start (e.g. after editing a DGP), delete the rds file at `out_path` or pass `resume = FALSE` to `run_simulation()`.

## What's prespecified vs configurable

**Prespecified (DGP and contrasts).** The five DGPs, the propensity-score model, the marginal RD (-0.05), the candidate grid, and the three workflows compared are all fixed in `R/`. The headline contrast is fixed-parametric vs fixed-rich vs plasmode-selected.

**Configurable (replicate budget and inner plasmode).** `n_per_rep`, `n_reps`, and `inner_reps` are arguments to `run_simulation()`. The pilot uses 40 outer replicates × 20 inner; the full study uses 500 outer × 50 inner.

## Five DGPs

All share the same baseline-covariate distribution (six variables: age, sex, biomarker, comorbidity, BMI z-score, smoke) and the same propensity-score model. They differ only in the outcome generation:

| DGP | Outcome model deviates from linear-in-baseline how |
|---|---|
| `linear` | Not at all — fixed-parametric library is correctly specified. |
| `nonlinear_smooth` | Quadratic biomarker + sinusoidal BMI. |
| `interactions` | Two strong two-way interactions. |
| `sparse` | Only biomarker matters; analyst includes all six. |
| `high_dim_noise` | Linear in the six real covariates plus 20 noise covariates in the adjustment set. |

## Three workflows

1. **Fixed-parametric TMLE.** SuperLearner library = `(SL.glm, SL.mean)`, truncation = 0.01.
2. **Fixed-rich TMLE.** SuperLearner library = `(SL.glm, SL.glmnet, SL.gam, SL.ranger, SL.mean)`, truncation = 0.01.
3. **Plasmode-selected TMLE.** cleanTMLE candidate grid (6 candidates spanning parametric / parametric+glmnet / rich / rich+screener × two truncation levels). Inner plasmode selects by min-RMSE. The FIORD two-stage rule is reported as a sensitivity.

## Performance metrics

Bias, RMSE, 95% CI coverage (with MC SE), empirical SD, mean SE, SE/SD ratio, and — for the plasmode workflow — the empirical candidate-selection cross-tab by DGP.

## Replicate budget and Monte Carlo standard errors

- **Pilot (`scripts/run_pilot.R`):** 10 outer replicates per DGP, 10 inner, n = 500. Coverage MC SE ≈ 0.07. Used to populate placeholders and end-to-end sanity-check the pipeline; takes roughly 5–7 hours on a single workstation. **Not inferentially defensible.**
- **Full (`scripts/run_full.R`):** 500 outer replicates per DGP, 50 inner. Coverage MC SE ≈ 0.010. Paper-quality.

## Dependencies

- `cleanTMLE` (this repository, `../cleanTMLE/`).
- `SuperLearner` and at least: `glmnet`, `gam`, `ranger`.
- `quarto` for rendering the manuscript.
- `dplyr`, `tidyr`, `ggplot2` for the manuscript chunks.

## Status

Pilot results exist at `results/sim_pilot.rds`. The manuscript draft auto-detects the pilot file and renders with a warning banner. The full study has not yet been run; placeholders in the manuscript (`[PLACEHOLDER]`) mark sections to be written or edited once `results/sim_full.rds` is generated.

## Open questions for collaborators

1. Should the DGP set include a sixth, positivity-strained DGP? The current five span outcome-model mis-specification but not positivity.
2. Should we vary the propensity-score model alongside the outcome model, or keep it fixed?
3. What is the right primary metric for the comparison — bias, RMSE, coverage, or a composite? The Phillips workshop materials lean coverage; the FIORD selector leans the same way.
4. Should plasmode selection be evaluated under multiple sample sizes (n = 500, 1000, 2000) to characterise the n at which the selection adds value? Cost is linear in n_per_rep and roughly quadratic if we also vary the inner plasmode.

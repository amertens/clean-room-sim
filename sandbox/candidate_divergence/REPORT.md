# Candidate divergence under data-quality stress: a constructed demonstration

## Purpose

This study demonstrates that the rule used to select among prespecified TMLE
candidates can change the candidate that is locked once a data-quality threat
sweep is taken into account. The main simulation in `run_simulation.R` does not
exhibit this divergence, because there the `min_rmse` and `min_max_rmse` rules
select the same candidate. The study here is a constructed existence
demonstration: the severities were chosen to exhibit the mechanism, so it shows
that divergence occurs under realistic threats rather than that divergence is a
generic property of the workflow.

## Data-generating process

The design reuses `generate_data()` and `compute_truth()` from
`run_simulation.R` without modification: five baseline covariates, a logistic
treatment model at marginal overlap (`overlap_strength = 1.6`), and a 24-month
binary outcome carrying a true risk difference of -0.0297 (`effect_size`
-0.05). A single reference cohort of 2,000 patients is drawn under a fixed
top-level seed (20260530) and held constant across the study.

## Candidate grid

The three candidates differ only in propensity-score truncation, with a
main-effects `SL.glm` propensity model in every case:

  - `aggressive`  truncation 0.001 (efficiency-first)
  - `middle`      truncation 0.025 (intermediate)
  - `robust`      truncation 0.20  (robustness-first)

## Threat sweep

Three prespecified threats at increasing severity, run through the Stage 2b
plasmode and scored by `select_tmle_candidate`:

  - `near_positivity`: the covariate-to-treatment association is amplified
    (centred propensity log-odds scaled by 2, 3, 4) so a subgroup approaches
    deterministic treatment and the estimated propensity score reaches the
    boundary. This threat is built in the sandbox because the package menu has
    no positivity stress; candidate fits reuse the package's own
    `.plasmode_fit_one_candidate`, and the engineered rows are merged into the
    package `plasmode_dq_results` object before selection.
  - `unmeasured_confounding`: a latent binary factor (prevalence 0.20) shifts
    treatment and outcome on the odds scale, odds ratio sweep 3, 4, 5, 6, 7, 8.
  - `covariate_missingness`: missing completely at random with median
    imputation, fractions 0.10 and 0.20.

## Locked thresholds

Declared before the stress test: maximum absolute bias 0.02, minimum coverage
0.88, and a maximum RMSE ratio of 2.0 relative to a candidate's own baseline.

## Monte Carlo design

200 synthetic replicates per cell, five independent batches (lock seeds spaced
by 1,000,000 so no synthetic draw is shared across batches). Monte Carlo
standard errors are the between-batch standard deviation divided by the square
root of the batch count.

## Result

| candidate  | truncation | baseline RMSE (MC SE) | worst-case RMSE (MC SE) | worst threat              |
|------------|-----------:|----------------------:|------------------------:|---------------------------|
| aggressive |      0.001 |      0.0214 (0.0004)   |      0.1481 (0.0011)    | near_positivity slope x4  |
| middle     |      0.025 |      0.0211 (0.0004)   |      0.0922 (0.0005)    | unmeasured_U OR 8         |
| robust     |      0.20  |      0.0184 (0.0003)   |      0.0977 (0.0006)    | unmeasured_U OR 8         |

`select_tmle_candidate(rule = "min_rmse")` selects `robust` (lowest baseline
RMSE). `select_tmle_candidate(rule = "min_max_rmse", dq_results = ...)` selects
`middle` (lowest worst-case RMSE). The two rules disagree, and the disagreement
holds in all five batches. The worst-case separation between the two selected
candidates is 0.0055 with a combined Monte Carlo standard error of 0.0008, a
ratio of about 6.7, so the change in decision is a property of the design
rather than simulation noise.

## Interpretation

Propensity-score truncation trades variance against robustness, and the two
selection rules reward opposite ends of that trade. At baseline the heaviest
truncation has the lowest RMSE, because under double robustness clipping the
smallest propensity scores removes weight variance at almost no cost in bias, so
`min_rmse` selects `robust`. The same heavy truncation is the most exposed to
strong unmeasured confounding, since the observations it discards carry the
information needed to offset the latent factor, so `robust` has the highest
worst-case RMSE among the moderate candidates (Panel B). The lightest
truncation has the opposite weakness and is destroyed by the positivity
violation, where its inverse-probability weights reach a mean maximum near 1,000
and its RMSE rises past seven times baseline (Panel A). The intermediate
truncation has neither failure at full severity, so it has the lowest worst-case
RMSE and `min_max_rmse` selects it. The minimax rule encodes a preference for
protection against the most damaging plausible threat over the best showing on
undisturbed data, and in this design those preferences select different
candidates.

## Files

  - `divergence_study.R`      self-contained study (modes: smoke, direction, full)
  - `results/candidate_divergence_full.rds`   feasibility and DQ tables, winners, config
  - `figures/degradation_gradient.png`        RMSE-ratio gradient, two threats, threshold line
  - `figures/selection_table.png`             baseline and worst-case RMSE with MC SE and winners
  - `results/selection_table.csv`             the selection table as data
  - `results/tuning_log.md`                   the path to the final configuration

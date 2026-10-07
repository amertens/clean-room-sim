# Pilot report: steps 1 and 2 (calibration, DGP diagnostics, smoke runs)

Date: 2026-10-07. Branch `simulation-pilot` in clean-room-sim. Stopped after step 2, as BUILD-BRIEF.md asks; the contrast pilot (step 3) has not been started.

cleanTMLE: version 0.4.1, installed into `simulation/lib` from the submodule's `main` branch as it stood at the time of install (commit 6770bbd, exported with `git archive` because the submodule's working tree was on its `vignettes` branch with local edits). The superproject's submodule pointer is 72bf62c; between 72bf62c and 6770bbd the package R code differs only in two note strings in `estimate_effect()` and the scope message of `simulate_design()`, neither of which the simulation reads.

Runtimes below were measured while another project's R job (`analysis/99_run_all.R`, with three parallel worker processes) was running on the same machine, so they may be inflated. A profile of five map repetitions run after the smoke step took 3.4 s per repetition at n = 500, against about 1.6 s per repetition inside the smoke runs, so the timing varies by roughly a factor of two with machine load.

## Step 1: calibration and DGP diagnostics

Population draw: n = 1e6, seed 20261007. Overlap is the share of the population with true pi outside 0.05 to 0.95. b0 sets the mean control-arm risk mu(0, X) to 0.15. The share outside depends on gamma only. b0 depends on gamma only for the inside model, and only in the fourth decimal.

| gamma | share outside | b0 inside | b0 near | b0 outside |
|---|---|---|---|---|
| 0.25 to 1.00 | 0.000000 | -2.061 | -3.007 | -2.030 |
| 1.25 | 0.000235 | -2.060 | -3.007 | -2.030 |
| 1.50 | 0.041898 | -2.060 | -3.007 | -2.030 |
| 1.75 | 0.158404 | -2.060 | -3.007 | -2.030 |
| 2.00 | 0.263758 | -2.060 | -3.007 | -2.030 |
| 2.25 | 0.345756 | -2.060 | -3.007 | -2.030 |
| 2.50 | 0.410962 | -2.060 | -3.007 | -2.030 |
| 2.75 | 0.464430 | -2.060 | -3.007 | -2.030 |
| 3.00 | 0.509422 | -2.060 | -3.007 | -2.030 |
| 3.25 to 4.00 | 0.547 to 0.632 | -2.060 | -3.007 | -2.030 |

Chosen gamma (the same for every model, because the treatment model does not depend on the outcome model). The rule picks the grid value closest to 0.02 among shares under 0.05, closest to 0.15 among shares in 0.10 to 0.20, and closest to 0.40 among shares over 0.30:

| Overlap | gamma | Share outside |
|---|---|---|
| good | 1.25 | 0.0002 |
| moderate | 1.75 | 0.158 |
| severe | 2.50 | 0.411 |

All three levels can be reached. The full grid is in `calibration_grid.csv`.

Inside model scale: the PLAN leaves the constant c in `c (0.5 z(logit pi) + 0.5 z(X beta))` unspecified. I set c = 1 / SD(0.5 z(logit pi) + 0.5 z(X beta)) on the population draw, so the prognostic index has SD 1, as the package's linear surface does with `prognostic = 1`. That gives c = 1.066, 1.058 and 1.053 at the three chosen gammas. The direction beta is `rnorm(4)` under seed 7, applied to X1 to X4 standardized by their exact population moments.

### Population truths and DGP diagnostics (Table S2)

Risk-difference scale. Prevalence is the overall mean of Y with treatment drawn from the true pi. Variance explained is Var(mu(A, X)) / Var(Y). The c-statistic of the true pi and the treated share depend on overlap only: 0.827 and 0.492 (good), 0.888 and 0.493 (moderate), 0.938 and 0.495 (severe). Every null scenario has all four truths equal to 0 and prevalence 0.150.

| Overlap | Model | ATE | trimmed ATE | ATT | ATO | Prevalence (effect) | Var. explained (effect / null) |
|---|---|---|---|---|---|---|---|
| good | inside | 0.1032 | 0.1032 | 0.1638 | 0.0942 | 0.231 | 0.262 / 0.104 |
| good | near | 0.0918 | 0.0918 | 0.1569 | 0.0808 | 0.227 | 0.433 / 0.244 |
| good | outside | 0.0985 | 0.0986 | 0.1362 | 0.0903 | 0.217 | 0.253 / 0.096 |
| moderate | inside | 0.1033 | 0.0969 | 0.1753 | 0.0880 | 0.236 | 0.269 / 0.104 |
| moderate | near | 0.0918 | 0.0902 | 0.1694 | 0.0715 | 0.234 | 0.439 / 0.244 |
| moderate | outside | 0.0985 | 0.0893 | 0.1424 | 0.0855 | 0.220 | 0.263 / 0.096 |
| severe | inside | 0.1033 | 0.0822 | 0.1845 | 0.0801 | 0.241 | 0.274 / 0.105 |
| severe | near | 0.0918 | 0.0641 | 0.1793 | 0.0585 | 0.239 | 0.441 / 0.244 |
| severe | outside | 0.0985 | 0.0781 | 0.1469 | 0.0810 | 0.223 | 0.269 / 0.096 |

Full precision is in `dgp_diagnostics.csv` and `calibration.rds`. The control-arm risk is 0.150 in every scenario, by construction. Step 1 took about 4.5 minutes.

## Step 2: smoke runs

Two datasets (idx 1 and 2) for each of moderate / inside / effect at n = 500 and severe / outside / effect at n = 2,000. For each dataset: the map with generate-treatment, the map with sample-treatment (reps = max_reps = 200), and the answer-key fit. Everything ran without error, through `run_stage()` with caching. A repeated call on a cached cell returned identical rows in 0.01 s. The summary functions `answer_key()` and `agreement()` also ran end to end on these rows; with two datasets their numbers are only a code check.

### Runtime per call (one dataset, mean of 2)

| n | Map, generate | Map, sample | Answer key |
|---|---|---|---|
| 500 | 320 s | 400 s | 0.6 s |
| 2,000 | 730 s | 941 s | 1.4 s |

Individual map runs ranged from 318 to 430 s at n = 500 and from 644 to 1,178 s at n = 2,000. assess_design() takes about 1 s; the rest is the 200 repetitions. In a profile of one map repetition, 88% of the time is in `tmle::tmle()`: 24 calls per repetition (4 surfaces x ATE, trimmed ATE and ATT x 2 truncations).

### Example map verdicts

Moderate / inside / effect, n = 500, dataset 1:

| Estimand | Generate: status | Generate: half-width | Sample: status | Sample: half-width | Unresolved (gen / sample) |
|---|---|---|---|---|---|
| ATE | infeasible | | infeasible | | yes / no |
| trimmed ATE | feasible (glm_t0.01) | 0.086 | infeasible | | yes / yes |
| ATT | infeasible | | infeasible | | no / no |
| ATO | feasible (glm_t0.01) | 0.076 | feasible (glm_t0.01) | 0.075 | yes / yes |

Design diagnostics on this dataset: share outside the band 0.208, ESS fraction 0.21 (ATE), 0.44 (trimmed ATE), 0.52 (ATT), 0.52 (ATO), largest weight 68 (ATE), 16 (trimmed ATE), 39 (ATT), 0.99 (ATO), preference-score equipoise share 0.24.

Severe / outside / effect, n = 2,000, dataset 1:

| Estimand | Generate: status | Generate: half-width | Sample: status | Sample: half-width | Unresolved (gen / sample) |
|---|---|---|---|---|---|
| ATE | infeasible | | infeasible | | no / yes |
| trimmed ATE | feasible (glm_t0.01) | 0.051 | feasible (glm_t0.01) | 0.050 | yes / yes |
| ATT | infeasible | | infeasible | | no / no |
| ATO | feasible (glm_t0.01) | 0.043 | feasible (glm_t0.01) | 0.043 | yes / yes |

Across all 8 map runs: ATE and ATT were infeasible every time, ATO was feasible every time, and the trimmed ATE was split (generate: feasible in 2 of 4 runs; sample: feasible in 3 of 4). The two variants disagreed on the trimmed ATE in both n = 500 datasets, in opposite directions. Many verdicts were marked unresolved, as expected with max_reps = reps. All rows are in `smoke.rds` and `step2_output.txt`.

### Answer-key check (40 datasets, not cached)

The two smoke datasets of moderate / inside at n = 500 both gave estimates about two SEs above the truth. To rule out a bug, I fitted the key on datasets 1 to 40 of three scenarios (`key_check.R`, output in `key_check_output.txt`). At n = 20,000 (one dataset per model, moderate overlap) every estimate was within 2.3 SEs of its truth.

| Scenario | ATE | trimmed ATE | ATT | ATO |
|---|---|---|---|---|
| good / inside, n = 500 | estimable (bias -0.005, cov 0.95) | estimable | estimable (bias -0.015, cov 0.90 to 0.925) | estimable |
| moderate / inside, n = 500 | not (bias 0.022, cov 0.85) | estimable (0.017, 0.925) | not (0.035 to 0.039, 0.68 to 0.70) | estimable (0.013, 0.925) |
| severe / outside, n = 2,000 | not (0.010 to 0.013, 0.83 to 0.85) | estimable (0.012, 0.925) | not (0.020 to 0.023, 0.63 to 0.68) | estimable (0.008, 0.925) |

So the answer key changes with overlap and estimand, and on these two scenarios the generate map's calls for ATE, ATT and ATO matched it. The trimmed ATE, which the key calls estimable in both, is where the map's verdict was split.

## Projected compute (sequential, from the smoke timings)

Full grid as specified: 18 scenarios x 2 n, R = 100 map runs per variant per scenario and M = 500 key datasets per scenario.

| n | Map hours (both variants) | Key hours | Total |
|---|---|---|---|
| 500 | 360 | 1.6 | 362 |
| 2,000 | 835 | 3.4 | 839 |
| **Total** | | | **1,201 h** |

Contrast pilot, per scenario: 100 map runs of both variants plus 200 key datasets take 20.0 h at n = 500 and 46.4 h at n = 2,000 (the key part is under 5 minutes). For 100 generate runs alone: 8.9 h at n = 500 and 20.3 h at n = 2,000.

## Concerns and proposed changes (PLAN.md is unchanged)

1. **Compute is the binding constraint.** 1,201 sequential hours is far beyond budget. Proposals, in order of value:
   - *Share map datasets across the outcome models.* The map never sees the outcome, so its verdict on a dataset is the same whatever the true outcome model is. If X and A were seeded from (overlap, n, idx) and Y from (scenario, n, idx), one map run would serve all six model x effect scenarios of an overlap level. The comparison by model relation would then also be paired on the same map runs. Map hours drop from 1,196 to 199. The brief currently seeds the whole dataset from the scenario id, so this needs your approval.
   - *Parallel workers.* `run_stage(workers = )` is implemented but untested, because of the one-R-job rule. The machine has 16 cores and 31 GB. With 8 workers the shared design would take about 25 h.
   - *reps = 100 instead of 200* halves map time but makes more verdicts unresolved. *One truncation (0.01)* would cut the tmle calls roughly in half: the trimmed ATE estimates were identical under both truncations in every key fit, and glm_t0.01 was the selected candidate in every feasible smoke verdict.
2. **The inside model is not exactly inside the declared family.** The plan uses the package's default surfaces, whose `baseline_risk` is 0.10, while the DGP's control-arm risk is 0.15. The map's surfaces also use their own seeded X beta direction, which is not the DGP's. The inside model therefore has the family's form, not one of its four surfaces. Proposal: set `surfaces = list(baseline_risk = 0.15)` in the plan (an analyst would know a registry's baseline risk), and describe "inside" as "of the family's form".
3. **Outcome prevalence.** The control-arm risk is 0.15, as the brief specifies. The overall prevalence is 0.22 to 0.24 in the effect scenarios (0.15 under the null). If PLAN section 4's "baseline risk about 0.15" was meant as overall prevalence, b0 needs recalibrating.
4. **Easy cells may dominate Table 1.** On the smoke and key-check evidence, ATO is feasible and estimable everywhere, and ATE and ATT are infeasible and not estimable at moderate and severe overlap. Pooled sensitivity and specificity could then look good because of easy cells, while the informative cells (trimmed ATE, and ATE and ATT near the good-to-moderate boundary) are few. Proposal: report Table 1 by estimand as well as pooled, and consider a fourth overlap level at gamma = 1.5 (4.2% outside), where ATE and ATT should sit near the estimability boundary. The good level at gamma = 1.25 has almost no patients outside the band (0.02%).
5. **Map verdicts are noisy near the boundary.** With max_reps = reps = 200, many decision cells are unresolved, and the trimmed ATE verdict flipped between datasets and between variants. This is the tool as implemented, so it is fair, but it will lower the map's agreement in the borderline cells.
6. **Trimmed ATE definition.** The truth trims on the true pi, while the estimator trims on the cross-fitted estimate. At these sample sizes this added no visible bias (key check), but it is part of what the answer key measures and should be stated.
7. **Implementation choices to confirm.**
   - The data seed and the plan seed are hashed separately from ("data" or "plan", scenario id, n, idx). If they shared one seed, the fold assignment would reuse the random numbers that drew X1, so folds would be correlated with the covariates.
   - The sample-treatment copy of `.one_rep` still makes the generate-treatment draw and discards it. The resampled rows, folds and outcome draws therefore use the same random numbers as the generate variant, so the two variants are paired by repetition.
   - The answer key also applies the package's usable-cell rule (at most 5% failed fits) when it judges a candidate. The brief names only the bias and coverage tolerances.
   - Cache keys include the full scenario row, n, idx, reps, the cleanTMLE version and the text of the stage's code files. Editing a code file therefore invalidates that stage's cache.

## Files created

- `simulation/R/dgp.R`, `R/map.R`, `R/key.R`, `R/run.R`, `R/summarize.R`, `simulation/pilot.R`, `simulation/setup.R` (committed)
- `.gitignore`: added `simulation/cache/` and `simulation/lib/` (committed)
- `simulation/pilot/`: `calibration.rds`, `calibration_grid.csv`, `dgp_diagnostics.csv`, `step1_output.txt`, `smoke.rds`, `step2_output.txt`, `key_check.R`, `key_check_output.txt`, this report
- `simulation/cache/` (git-ignored): 12 smoke cells; `simulation/lib/` (git-ignored): cleanTMLE 0.4.1

## Commits (branch `simulation-pilot`, not pushed)

- 74488b2 Simulation pilot code: DGPs, support map (generate and sample treatment), answer key, cached runner, summaries
- 82ad787 Simulation pilot steps 1 and 2: calibration, DGP diagnostics, smoke runs, report

PLAN.md and BUILD-BRIEF.md are still untracked; I did not commit them. Unrelated changes (rescueCo, the cleanTMLE submodule pointer) are untouched and unstaged.

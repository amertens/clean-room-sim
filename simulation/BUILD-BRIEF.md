# Build brief: pilot code for the estimability simulation

Read first: `simulation/PLAN.md` (the agreed design; sections 1 to 4) and the skill files
`C:/Users/andre/.claude/plugins/cache/ctml-skills/stats-research/0.12.0/skills/design-and-report-simulations/references/implementation.md`,
`.../references/piloting.md`, and `.../assets/simulation-scaffold.R`, `.../assets/pilot_check.R`, `.../assets/performance_measures.R` (read and adapt; do not copy blindly).

## Where

All code in `C:/Users/andre/OneDrive/Documents/clean-room-sim/simulation/R/`, run scripts in `simulation/`, outputs in `simulation/cache/` (add `simulation/cache/` to clean-room-sim's `.gitignore`) and pilot results in `simulation/pilot/`. The cleanTMLE package is used as installed from the submodule at its current `main` (0.4.1). Install it once into a project library, `simulation/lib/` (git-ignored), with `install.packages("../cleanTMLE", repos = NULL, type = "source", lib = "simulation/lib")`, and put that library first on `.libPaths()` in every script, so results do not depend on whatever cleanTMLE is installed in the user library. Use package internals through `cleanTMLE:::` where needed; do not edit the package.

## Code structure (one codebase for pilot and full run)

1. `R/dgp.R`
   - `gen_covariates(n, seed)`: X1 ~ Uniform(-2, 2), X2 ~ Bernoulli(0.5), X3, X4 ~ Uniform(-1, 1).
   - `true_pi(X, gamma)`: expit(-0.2 + gamma X1 + 0.3 X2 + 0.3 X3).
   - `true_mu(model, a, X, pi, b0)` for model in "inside", "near", "outside", and `null = TRUE` versions (A terms removed), exactly as in PLAN.md section 4. For "inside", z() standardizes over the population draw used for calibration (store the means and SDs with the scenario), and X beta uses a fixed seeded direction over X1 to X4.
   - `calibrate(model, gamma_grid)`: on one large draw (n = 1e6), choose b0 so the control-arm risk mu(0, X) averages 0.15, and report for each gamma the share of pi outside 0.05 to 0.95.
   - `population_truth(scenario)`: the true ATE, trimmed ATE (patients with TRUE pi in 0.05 to 0.95), ATT (weights pi), ATO (weights pi(1 - pi)) on the risk difference scale, from the 1e6 draw.
   - `draw_data(scenario, n, idx)`: X, A ~ Bernoulli(pi), Y ~ Bernoulli(mu(A, X)); seed derived from a hash of (scenario id, n, idx) only, never from loop position.
   - Scenario table: overlap level (good, moderate, severe, with gamma from calibration) x model (inside, near, outside) x effect (non-null, null), each with an id.
2. `R/map.R`
   - `run_map(scenario, n, idx, variant = c("generate", "sample"), reps)`: build the plan (outcome "y"; estimands ATE, trimmed_ATE, ATT, ATO; trim_band c(0.05, 0.95); ps_library "glm"; candidates list(library = list(glm = "glm"), truncation = c(0.01, 0.05)); default surfaces; K = 5; reps = max_reps = reps; seed from the dataset hash), lock the design data (no y column), assess_design(), then:
     - variant "generate": `simulate_design()` as is.
     - variant "sample": the same loop with each patient's OBSERVED treatment kept instead of drawing A* from the propensity score. Implement by a local copy of `cleanTMLE:::.one_rep` that sets the treatment to `A_obs[idx]` and is otherwise identical, then the package's `.sim_metrics()` and `.sim_verdict()` on its results (no extension batches, since max_reps = reps). Keep the copy short and comment the single difference.
   - Return one row per estimand: verdict (feasible or infeasible), selected candidate, ci_halfwidth, unresolved, and the design diagnostics for that estimand from the same assess_design() call: share outside the band, max_weight, effective sample size as a fraction of n (sum of both arms over n), and the share of patients with a preference score (Walker) in 0.3 to 0.7, computed from the cross-fitted propensity score. Also the runtime.
3. `R/key.R`
   - `run_key(scenario, n, idx)`: draw the dataset with its outcome, fit `cleanTMLE:::.fit_and_target()` once with the GLM candidate library for all four estimands and both truncations (same plan, folds from `cleanTMLE:::.make_folds(n, 5, seed)`), and return per (estimand, truncation): estimate, se, ci_lower, ci_upper, failed. Also the per-estimand smallest-p-value bookkeeping needed for claim 3 (p from estimate and se).
4. `R/run.R`
   - `run_stage(stage, scenarios, ns, idx)`: for each (scenario, n, idx) not already cached, compute and save `cache/<stage>/<scenario>_n<n>_i<idx>.rds`; a cache hit equals recomputation (key on everything that determines the result, including the package version and reps). Runs sequentially by default; optional `future.apply` with a set number of workers. Any subset runs alone and returns the same rows as the full run.
   - `collect(stage)`: bind cached rows into one data frame.
5. `R/summarize.R`
   - Answer key per (scenario, n, estimand): for each candidate, bias and coverage against population_truth over the key datasets, with Monte Carlo SEs; estimable = some candidate meets |bias| <= 0.02 and coverage >= 0.90 (point estimates, as the package does).
   - Agreement: for each verdict source (map generate, map sample, and each diagnostic with its best threshold by Youden's index chosen across all scenario-datasets) sensitivity and specificity by model relation, with Monte Carlo SEs.

## Pilot script `simulation/pilot.R` (run in this order; stop and report after each)

1. Calibration and diagnostics: gamma grid (for example 0.25 to 4 by 0.25) per model; pick gamma for good (< 5% outside), moderate (10 to 20%), severe (> 30%); b0; population truths; Table S2 diagnostics (truths, share outside, c-statistic of true pi, outcome prevalence, variance explained). Save `pilot/calibration.rds` and print the table.
2. Smoke: 2 datasets each of run_map (both variants) and run_key, for one scenario at n = 500 and one at n = 2000; record the runtime of each call. Print projected hours for: R = 100 map runs (both variants) and M = 500 key datasets, per scenario, for the full grid (3 overlap x 3 models x 2 effects x 2 n).
3. STOP after step 2 and report. Do not start the contrast pilot (about 100 map runs) until the controller confirms, because its cost depends on step 2's timing.

## Environment

- R: "/c/Program Files/R/R-4.4.2/bin/Rscript.exe". One R job at a time (this machine crashes under concurrent R jobs); long jobs with the Bash tool's run_in_background and wait. R processes you did not start are not yours.
- A shell guard blocks output redirection to paths built from variables; use literal paths.
- ASCII-only R code. Commit the code (not cache, lib or pilot outputs larger than a few hundred KB) in clean-room-sim on a new branch `simulation-pilot`, with messages ending "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>". Never push. Leave unrelated uncommitted files in clean-room-sim (rescueCo changes, top-level C1..C4/PRIMARY folders) untouched and unstaged.
- Do not dispatch subagents.

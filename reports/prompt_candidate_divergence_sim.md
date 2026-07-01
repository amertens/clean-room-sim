# Claude Code prompt: engineer a candidate-divergence scenario

Purpose. The cleanTMLE paper claims the data-quality stress test lets an analyst
select a candidate estimator that is robust across data-quality threats rather
than the one with the best baseline RMSE, via `select_tmle_candidate(rule =
"min_max_rmse")`. The current simulation never demonstrates this: the candidate
grid in `run_simulation.R` is two near-identical GLM truncation candidates
(`glm_t01`, `glm_t05`) that produce the same metrics, so `min_rmse` and
`min_max_rmse` always agree. The task is to construct a scenario in which the two
rules select different candidates, where the `min_rmse` winner is best at
baseline but collapses under a prespecified threat and the `min_max_rmse` winner
is slightly less efficient at baseline but stable across threats.

Paste everything below into a coding session with the `clean-room-sim`
repository open.

---

```
You are working in the clean-room-sim repository. The cleanTMLE R package source
is in cleanTMLE/. The existing simulation harness is run_simulation.R at the repo
root; reuse its DGP and patterns. Do NOT edit run_simulation.R or the package
source unless a bug forces it; build this study as a standalone, self-contained
script so it cannot disturb the main results.

GOAL
Produce a Monte Carlo demonstration in which two prespecified TMLE candidates
diverge under data-quality stress, such that:
  - under the baseline plasmode, candidate "aggressive" has the lowest RMSE and
    would be chosen by select_tmle_candidate(rule = "min_rmse");
  - under at least one prespecified data-quality threat, candidate "aggressive"
    degrades sharply (large RMSE / bias) while candidate "robust" stays within
    the locked thresholds;
  - select_tmle_candidate(rule = "min_max_rmse", dq_results = ...) therefore
    selects "robust", a DIFFERENT winner than min_rmse.
The deliverable proves that the minimax rule changes the decision in a realistic
setting, which the main simulation does not currently show.

SETUP
1. Create a sandbox folder: sandbox/candidate_divergence/ with subfolders
   results/ and figures/. Put all scripts there. Nothing outside this folder
   should be modified.
2. Load cleanTMLE with devtools::load_all("cleanTMLE") (or
   pkgload::load_all) so you run against the working tree, and record
   sessionInfo() and the package version into results/session_info.txt.
3. Confirm the API by reading the roxygen for these functions in
   cleanTMLE/R/cleanroom.R and cleanTMLE/R/plasmode_dq.R before calling them:
   create_analysis_lock, attach_estimand, tmle_candidate,
   run_plasmode_feasibility, run_plasmode_dq_stress, summarize_dq_degradation,
   select_tmle_candidate (rules: min_rmse, min_bias, max_coverage,
   min_max_rmse), lock_primary_tmle_spec. Note run_plasmode_dq_stress takes
   data_quality_scenarios as a named list and has a fit_timeout argument; use a
   finite fit_timeout so a runaway fit cannot hang the run.

DESIGN HYPOTHESIS (start here, then tune empirically)
- DGP: reuse the generate_data(n, overlap_strength, effect_size, ...) function
  from run_simulation.R (copy it into the sandbox script; do not source the
  whole harness). Use marginal overlap (overlap_strength around 1.5 to 2.0) so
  that propensity-score truncation actually bites. Keep a real, nonzero effect
  (effect_size around -0.05) and compute the truth with the same compute_truth
  routine.
- Candidates: define a grid where truncation is the lever that trades baseline
  efficiency against fragility under stress. Begin with:
    aggressive: tmle_candidate("trunc_min", g_library = "SL.glm",
                truncation = 0.001)   # tiny truncation -> low baseline RMSE
                                      #   but heavy-tailed weights under stress
    robust:     tmle_candidate("trunc_hi",  g_library = "SL.glm",
                truncation = 0.10)    # heavier truncation -> slightly higher
                                      #   baseline RMSE, stable under stress
  Optionally add a middle candidate (truncation = 0.025) to show a gradient.
  AVOID SL.glmnet in the PS library: run_simulation.R documents that glmnet runs
  away on near-positivity-violation synthetic designs. Stay with SL.glm-based
  candidates, or cap any flexible learner and set fit_timeout.
- Threat that should break the aggressive candidate: the most likely lever is
  unmeasured confounding (latent U shifting both A and Y) combined with marginal
  overlap, because tiny truncation lets a few near-zero propensity scores create
  extreme weights once U distorts the treatment mechanism; covariate missingness
  at 0.20 is a secondary candidate threat. Configure data_quality_scenarios with
  a severity sweep, e.g. unmeasured_confounding with U_treatment_OR and
  U_outcome_OR in c(1.5, 2.0, 3.0) and covariate_missingness fractions in
  c(0.05, 0.10, 0.20).

ITERATE IN THIS ORDER (do not jump to the full run)
Step 1 - Smoke test (minutes). Run the whole pipeline once at tiny scale:
  n = 1000, reps = 5, dq reps = 5, one effect size. Confirm it executes end to
  end with no error, that run_plasmode_feasibility returns a metrics table with
  DIFFERENT baseline RMSE across candidates, and that run_plasmode_dq_stress
  returns a degradation table. Print the baseline metrics and the per-scenario
  worst-case RMSE per candidate.
Step 2 - Direction check (low reps). At n = 2000, reps = 20, verify the
  qualitative pattern: aggressive has the lowest baseline RMSE, and under the
  strongest stress cell aggressive has a HIGHER worst-case RMSE than robust.
  Then check the decision:
    win_min   <- select_tmle_candidate(plas, rule = "min_rmse")
    win_minmax<- select_tmle_candidate(plas, rule = "min_max_rmse",
                                       dq_results = dq)
  The success condition is win_min$candidate_id != win_minmax$candidate_id, with
  win_min == aggressive and win_minmax == robust.
Step 3 - Tune if the rules agree. If they do not yet diverge, adjust ONE knob at
  a time and re-run Step 2, in this priority order:
    (a) widen the truncation gap (e.g. 0.001 vs 0.20);
    (b) increase overlap_strength so low truncation is more fragile;
    (c) increase the top U strength (OR 3.0 or 4.0) or the top missingness
        fraction;
    (d) lower max_abs_bias / tighten the locked thresholds so the fragile
        candidate trips them.
  Keep a short tuning log (results/tuning_log.md) recording each knob change and
  the resulting baseline RMSE, worst-case RMSE, and the two winners, so the path
  to the final configuration is reproducible.
Step 4 - Debug guidance. If a candidate errors or produces NaN: inspect the
  weight tails (max weight, ESS) of the aggressive candidate under the stress
  cell; confirm the plasmode generator's realized effect matches effect_size;
  set a finite fit_timeout; reduce flexible-learner use. If baseline RMSE is
  identical across candidates, truncation is not biting: raise overlap_strength
  until propensity scores reach the truncation bounds (print the PS range).
Step 5 - Full run (once Steps 2-3 show stable divergence). Scale to n = 2000,
  reps = 200, dq reps = 200, fixed seed. Report Monte Carlo standard errors for
  every RMSE, bias, and coverage figure so the divergence is demonstrably not
  noise. The divergence must hold with the min_rmse and min_max_rmse winners
  separated by more than their combined Monte Carlo error on worst-case RMSE.

OUTPUTS
1. results/candidate_divergence_results.rds: the plasmode feasibility table, the
   DQ degradation table, both rule winners, and the config used.
2. figures/degradation_gradient.png: a degradation-gradient plot with threat
   severity on the x axis, RMSE (or RMSE ratio to baseline) on the y axis, one
   line per candidate, and a horizontal line at the locked RMSE/bias threshold.
   The aggressive line must cross the threshold while the robust line stays
   under it.
3. figures/selection_table.png or a kable: baseline RMSE, worst-case RMSE across
   DQ scenarios, min_rmse winner, min_max_rmse winner, each with MC SE.
4. sandbox/candidate_divergence/REPORT.md: a one-page write-up stating the DGP,
   the candidate grid, the threat sweep, the locked thresholds, the two winners,
   and one paragraph interpreting why the minimax rule changed the decision.
   Write it in flowing prose suitable to drop into the manuscript's simulation
   section; follow the project style rules (no em-dashes, no "not X but Y", no
   editorializing adverbs, do not use "GO/FLAG/STOP" generically).

ACCEPTANCE CRITERIA (all must hold on the full run)
  - win_min == aggressive, win_minmax == robust, and they differ.
  - aggressive has the strictly lowest baseline RMSE.
  - aggressive has the strictly highest worst-case RMSE across the DQ sweep, and
    that worst case exceeds the locked threshold while robust's does not.
  - the separation in worst-case RMSE exceeds the combined Monte Carlo SE.
  - the degradation-gradient figure visibly shows the aggressive candidate
    crossing the threshold and the robust candidate staying under it.
Stop and report if, after the Step 3 tuning options are exhausted, the rules
still agree; in that case summarize what was tried and what the closest
divergence looked like, rather than forcing a result.

CONSTRAINTS
  - Self-contained: everything under sandbox/candidate_divergence/; do not touch
    run_simulation.R, results/, results_new/, or the package source except to
    fix a genuine bug (and if you must, isolate the fix and flag it in REPORT.md).
  - Reproducible: one fixed top-level seed, recorded in the config and the rds.
  - Show me the Step 1 smoke-test output and the Step 2 direction-check winners
    before launching the Step 5 full run.
```

---

## Why this design should work (reviewer's note, not part of the prompt)

Truncation is the cleanest lever for engineered divergence because it trades
variance against bias in a way that the stress threats amplify asymmetrically. A
tiny truncation threshold gives the lowest variance, hence the lowest RMSE, when
overlap is adequate, so it wins the baseline `min_rmse` comparison. Under
marginal overlap combined with unmeasured confounding or heavy covariate
missingness, the same tiny threshold lets a handful of near-zero propensity
scores produce extreme weights, and the candidate's RMSE and bias inflate. A
heavier truncation threshold sacrifices a little baseline efficiency but caps the
weights and stays inside the locked envelope under stress, so it wins the
worst-case `min_max_rmse` comparison. That asymmetry is exactly the
robust-versus-fragile contrast the minimax rule exists to resolve, and it is the
contrast the current two-candidate GLM grid cannot show because both candidates
truncate in the same benign range.

If truncation alone proves too weak a lever on this DGP, the next most reliable
source of divergence is library complexity: an unscreened flexible learner that
overfits the synthetic outcome surface at baseline against a parsimonious learner
that does not. That path is secondary because of the documented SL.glmnet
runaway-compute issue on near-positivity designs, so it should be attempted only
with a finite fit timeout and a capped learner.

# Run manifest — reviewer-revision simulations

All code below is integrated and smoke-validated. Run **one job at a time**
(this machine crashes long/concurrent R jobs). Each is launched from the repo
root with the bundled R: `"C:/Program Files/R/R-4.4.2/bin/Rscript.exe"`.

The cleanTMLE package has already been reinstalled with the MNAR feature; if you
edit the package again, reinstall first:
`Rscript -e "devtools::install('cleanTMLE', quick=TRUE)"`.

| # | Command | Purpose | Output | Approx. runtime |
|---|---------|---------|--------|-----------------|
| 1 | `SIM_REPS=1000 Rscript run_simulation.R` | Main sim incl. **Scenario D** (misspecified surface) and **MAR/MNAR** DQ threats, MC SE for all measures | `results_new/` | hours (4 scenarios; Scenario D uses gam) |
| 2 | `SIM_REPS=1000 SIM_NOBS=500  SIM_RESULTS=results_n500  Rscript run_simulation.R` | N-sweep point | `results_n500/` | hours |
| 3 | `SIM_REPS=1000 SIM_NOBS=5000 SIM_RESULTS=results_n5000 Rscript run_simulation.R` | N-sweep point | `results_n5000/` | hours |
| 4 | `Rscript _bootstrap_variance.R` | Variance study now covering **TMLE & TMLE_CF** (not just IPTW/Match_TMLE) | `results_new/bootstrap_variance.{rds,csv}` | ~1-2 h |
| 5 | `Rscript rescueCo/scripts/rescueco_dq_full.R 200 50` | Case-study DQ with all 5 threats + MAR/MNAR, bounded learner | `rescueCo/results/plasmode_dq_*.{rds,csv}` | ~5-6 h |
| 6 | `Rscript sandbox/validation/validate_vs_tmle.R 20000` | Tier-2 validation table | `sandbox/validation/validation_vs_tmle.csv` | ~2 min |
| 7 | `Rscript sandbox/candidate_divergence/_stability.R` | #3 selection-stability table (from existing divergence batches) | console | seconds |

Notes
- Rep count, sample size, DQ reps, and output dir are env-overridable in
  `run_simulation.R` (`SIM_REPS`, `SIM_NOBS`, `SIM_DQREPS`, `SIM_RESULTS`).
- If a long run dies partway (the observed failure mode), just relaunch it; none
  of these overwrite inputs, only outputs.
- For reliability, consider: add the R install dir and the repo to the antivirus
  exclusion list, and pause OneDrive sync on `clean-room-sim/` during runs.
- Smoke first if unsure: pass small args, e.g. `SIM_REPS=20 Rscript run_simulation.R`,
  or `Rscript rescueCo/scripts/rescueco_dq_full.R 3 3`.

## What each run establishes (reviewer items)
1. **#1 Misspecified DGP (Scenario D)** — under a nonlinear, effect-modified
   surface, GLM-based estimators are biased and undercover while the flexible
   cross-fitted TMLE recovers the truth. Smoke: crude bias +0.13 (wrong sign),
   GLM-TMLE bias +0.04 (cov 0.67), SL-TMLE bias +0.01 (cov ~1.0).
2. **#2 Variance** — TMLE and TMLE_CF are anti-conservative under marginal
   overlap (IF SE/SD ~0.85, cov ~0.87) and the bootstrap does **not** fully
   repair it (the gap is partly bias). Run #4 quantifies this at scale.
3. **#3 Selection stability** — `min_max_rmse` is stable (5/5 batches) when
   candidates are well separated (margin = 6.7x MC SE); it flips only when they
   are nearly tied (the case study). Run #7 prints the table.

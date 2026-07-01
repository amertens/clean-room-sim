# Package bug fixes found while building the Section 9.4 separating-grid study

These are real code fixes (not manuscript edits), discovered and applied with your
approval while running the separating-grid study. Recorded here so the revision
task has a complete account of package changes. None changes a reported numeric
result.

## 1. `run_plasmode_dq_stress()` crashed on cells with fewer than 2 converged reps

- **File / line:** `cleanTMLE/R/plasmode_dq.R:1027`.
- **Symptom:** `Error in if (emp_sd > 0) ... : missing value where TRUE/FALSE
  needed`, aborting the entire stress test.
- **Cause:** `se_cal = round(if (emp_sd > 0) mean_se / emp_sd else NA, 3)` uses a
  scalar `if` on `emp_sd`. When a degraded cell has fewer than two converged
  replicates, `emp_sd = sd(<2 values) = NA`, and `if (NA)` errors. This is
  reachable at any `reps` (a single degraded cell with one convergence triggers
  it), and certain at `reps = 1`. It bites hardest exactly where the stress test
  matters most: aggressive truncation under a positivity or confounding threat,
  where non-convergence is most likely.
- **Why it is a latent inconsistency:** the sibling function
  `run_plasmode_feasibility()` already computes the same quantity with the
  NA-safe vectorised form at `cleanTMLE/R/cleanroom.R:1468-1469`
  (`ifelse(metrics$emp_sd > 0, ...)`). The two functions diverged; the stress
  path was never NA-guarded. This is the same class of plumbing-robustness gap
  the manuscript reports the rescueCo run surfaced (NA-outcome guards, near-zero
  variance handling).
- **Fix applied:**
  ```r
  se_cal = round(if (!is.na(emp_sd) && emp_sd > 0)
                   mean_se / emp_sd else NA_real_, 3),
  ```
  `se_cal` is `NA` in the affected cells either way, so no previously reported
  number changes; the guard only prevents the abort.
- **Regression test:** `cleanTMLE/tests/testthat/test-plasmode_dq.R`, new test
  "run_plasmode_dq_stress is NA-safe when a cell has < 2 converged reps" (reps = 1
  must return a `plasmode_dq_results` with all-`NA` `se_cal` rather than erroring).
- **Verification:** package reinstalled; full `test-plasmode_dq.R` passes
  (all tests green); the previously-crashing reps = 1 repro now returns cleanly
  for all four anchored threats.

### Suggested follow-up for the revision task (not done here)
A reviewer-facing note: consider auditing the rest of `plasmode_dq.R` and
`cleanroom.R` for other scalar `if (<stat>)` guards on quantities that can be `NA`
under non-convergence (e.g. coverage or bias on empty `valid` subsets), and
factor the `se_cal` computation into one shared helper so the feasibility and
stress paths cannot drift again.

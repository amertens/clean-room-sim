# Simplifying cleanTMLE for a 0.5.x release

Proposal, 2026-10-22. Written from the read-only audit of 22 October 2026: `audit/2026-10-22/REPORT.md`, recommendations 4 and 7 to 17, with the detail in `01-parsimony.md`, `02-wrapper.md`, `06-check-guarantees.md`, `07-check-framing.md` and `08-literature.md`. Where an audit and the report differ, this spec follows the report. Nothing here is approved: every decision below awaits Andrew, and no package code is written until he has answered them. Line numbers refer to `cleanTMLE` at commit `2e96f98` on branch `vignette-limits`, read in the working tree on 22 October; the one uncommitted edit there touches only `vignettes/cleanTMLE.Rmd` (its new lines `290-301` and `367-378`). Calls made in writing the spec are marked **[call]** and listed at the end.

## Aim

Make the package smaller and easier to learn without changing any number it produces, and make the guarantee that the simulation tests the estimator applied after unblinding rest on a recorded build rather than on a version string. The cuts remove code that the case study, the vignettes and the paper do not use: about 400 lines of R and two packages from Imports (report, summary and size table); five packages leave Suggests, and cobalt moves there from Imports. Before any cut, tests pin the in-house numerical code to the established packages it stands in for, and a stored run of the design stage on 0.5.0 fixes the numbers that every later step must reproduce.

## Scope

In scope are recommendations 4 and 7 to 14 (code and tests), 15 (vignettes and `README.md`; the proposed opening of the main vignette is in `10-proposed-text.md`, section 1), 16 (housekeeping) and 17 (four small items from the literature).

Recommendations 3 and 6 enter only as a pointer. They belong to the learner-defaults item, which has its own brainstorm and spec (`docs/superpowers/specs/2026-10-18-outcome-blind-design-stage.md`, section "Package item: learner defaults"). This spec changes no learner, tuning grid, penalty, default library or targeting step; `R/learners.R` is read by the new tests and not edited. The rules recommendation 3 sets for that item (measure before fixing, write each decision rule into its spec before the run, judge every fix on simulated data with a known truth) govern it, and step 1 below gives it the tests it needs.

Out of scope: the Rescue.Co case study, which stays on 0.5.0 as installed; the simulation code in `simulation/R/`, which recommendation 26 moves to 0.5.x after pilot step 3a; the paper, apart from noting which sentences a step makes stale; and recommendations 31 to 36, which wait until after the paper (`save_lock()`, the `after_treatment` role, restrictions, jsonlite to Suggests, the principal-component highly adaptive lasso as a learner, a fitter named in the plan). No export is added, and none is removed.

## Decisions for Andrew

Each awaits approval and carries the report's recommended choice. D1 comes first, because four other recommendations depend on it (report, recommendation 11). D11 to D13 are raised by this spec, not by the report.

### D1. Comparators in the package (recommendation 11)

Status: awaiting approval.

Recommended: cut both comparators, `.estimate_weighting()` and `.estimate_matching()`. Write the comparators the simulation needs in `simulation/R/` with the same flexible propensity learner as the TMLE arm. If the paper keeps a comparator in the case study, call WeightIt in script 25 and label it as not evaluated by the simulation.

Reason: script 25 calls neither, the simulation calls neither, and the package's weighting comparator uses a main-terms logistic score without truncation (`R/estimate_effect.R:39-40`), so a difference from the TMLE mixes the learner, the truncation and the estimator. Matching estimates the ATT only and relies on a `lm_weightit(weights = )` argument that later WeightIt versions ignore, by audit 02's reading.

### D2. The fitter stays in-house (recommendation 4)

Status: awaiting approval.

Recommended: keep the fitter in-house and pin each component to the package it stands in for; a fitter named in the plan (recommendation 34) is the route to delegation after the paper.

Reason: one fit per fold yields a whole boosting or elastic-net path (`R/learners.R:123-147`), bootstrap copies share folds (`R/learners.R:117-121`), and SuperLearner's wrappers differ in fold draws, leaf sizes, net settings and an inner `cv.glmnet` (`02-wrapper.md:28-36`), so delegation would change every design-stage number.

### D3. A build that differs at estimation (recommendation 7)

Status: awaiting approval.

Recommended: `estimate_effect()` refuses a build that differs from the one that ran the simulation, unless a written override is given; elsewhere a difference warns.

Reason: the lock hashes the version string and not the code (`R/lock.R:17-23`), and a mismatch only warns (`R/lock.R:159-162`), so a rebuild that kept the version number passes silently.

### D4. Time-to-event path (recommendation 9)

Status: awaiting approval.

Recommended: cut it and list it in `WISHLIST.md`.

Reason: it is the one route on which the estimator run after unblinding is not the one the simulation evaluated (`R/estimate_effect.R:168`), and neither the case study, the vignettes nor the manuscript uses it.

### D5. E-value summary and implausibility flag (recommendation 10)

Status: awaiting approval.

Recommended: cut both.

Reason: only the print method reads the flag (`R/methods.R:91`), and a sign that differs from the crude difference is what adjustment for confounding can produce.

### D6. Interval rule only for negative controls (recommendation 12)

Status: awaiting approval.

Recommended: keep only the interval rule, and keep `rule = "interval"` accepted so that the case-study plans still validate.

Reason: the Methods and the case study use only that rule (`manuscript/methods.qmd:102-107`; `rescueCo/scripts/22_design_stage.R:138`).

### D7. `select` and `surfaces$q0` (recommendation 13)

Status: awaiting approval.

Recommended: remove both.

Reason: an analyst-written selection function is a degree of freedom in a package built to remove them, and `surfaces$q0`, which has no test, would let predictions fitted to the real outcome by a steward enter the design stage.

### D8. Vignettes: register, re-locking policy and terms (recommendation 15)

Status: awaiting approval.

Recommended: keep "you" and "we" (`style/voice.md`, tutorial register); allow re-locking with looser tolerances after a negative result only as a new lock, with the first lock's design report reported beside it; use "design report" for what `design_report()` returns, "proxy covariate" for a stand-in, and give the verdict labels `GO`, `STOP` and `FLAG` once each, after a description in words.

### D9. Delete `attic/` (recommendation 16)

Status: awaiting approval.

Recommended: yes; the `v0.3.0` tag keeps the history.

### D10. Equipoise share and the registration format (recommendation 17)

Status: awaiting approval.

Recommended: add the share of patients with a preference score between 0.3 and 0.7 to `assess_design()` as a reported field that gates nothing; leave a registration format for `export_design_log()` on the wishlist, since the demonstration has no registration.

### D11. Version number (raised by this spec)

Status: awaiting approval.

Recommended: 0.5.1, with the `v0.5.0` tag on the commit that `rescueCo/lib` was built from. Nothing numerical changes, and the arguments removed (`method`, `select`, the time-to-event arguments) had no user outside this project, since 0.5.0 is not yet public. The alternative is 0.6.0, which would signal the removed arguments by the usual convention. Either way, if this release is published before 0.5.0, the paper says that the case study ran on 0.5.0 and what the later release changed (report, recommendation 3).

### D12. One branch (raised by this spec)

Status: awaiting approval.

Recommended: a single branch, `simplify-0.5.1`, with one commit per step and the tests run after each step. The steps are ordered and every cut leans on the tests of step 1, so separate branches would add more branches to combine and more checks to run, without making any step independent.

### D13. Order relative to the learner-defaults item (raised by this spec)

Status: awaiting approval.

Recommended: step 1 (tests and pins) runs first, since the learner item needs the same tests; the rest of this spec follows on the same branch and goes into main before the learner item changes any code, so that the learner item's numerical changes are judged against a package whose numbers are pinned and unchanged. If Andrew decides that a confirmed learner defect holds the 0.5.0 release (report, recommendation 3), the learner item may go first after step 1, and this spec rebases onto it; the reference values of step 1 are then regenerated on the learner build before step 2 begins, with that regeneration recorded in `NEWS.md`.

## Guarantees that must not change

Every step keeps the following, and the tests named beside each must still pass after every step, with what they assert unchanged; step 9 rewrites only the call at `test-simulate-design.R:183`, which passes the arguments it removes, to use `pairs`.

The fitter is shared. The simulation and the final analysis fit and target through the same function, `.fit_and_target()` (`R/fit_candidate.R:164-191`), reached from `.one_rep()` (`R/simulate_design.R:112`) and from `fit_candidate()` (`R/fit_candidate.R:193-196`), which `estimate_effect()` (`R/estimate_effect.R:180`) and `negative_control_ladder()` (`R/nc_ladder.R:229`) call with the candidate, folds and seed the design report recorded. No step adds a second estimation path. Protected by `test-estimate-effect.R:4` and `test-nc-ladder.R:342` and `:366`.

Cross-fitting is unchanged: `.cross_fit()` (`R/fit_candidate.R:13-29`), `.make_folds()` (`R/utils.R:45-48`), `.inner_folds()` (`R/learners.R:117-121`) and the grouping of bootstrap copies in `.one_rep()` (`R/simulate_design.R:100-104`). Protected by `test-fit-candidate.R:11` and the new test of `.fit_Q()` (step 1).

The edge check is unchanged: `.cv_setting()` keeps each setting's selected point, path length and risk (`R/learners.R:123-147`), and `.edge_check()` reads them (`R/learners.R:216-253`). Protected by `test-learners.R:53` and the glmnet and nnet cases of step 1.

The hashes are unchanged. `.hashable()` and `.hash()` (`R/utils.R:3-32`), `.dossier_hash()` (`R/dossier.R:52-55`), `.stamp()` and `.check_stamp()` (`R/lock.R:187-201`) and `.unblind_hash()` (`R/unblind.R:13-15`) keep their recipes. `.lock_hash()` (`R/lock.R:17-23`) keeps its elements in this order: plan, data, id, treatment, covariates, seed, K, version, dictionary. Step 3 makes the dictionary unconditional in that last position, which gives the same hash for every lock made by 0.5.0. A plan field removed from the code stays as a harmless extra field in stored objects, because the stored plan is what is hashed. Protected by `test-lock.R:18` and the stored 0.5.0 objects of step 1.

The blinding in `unblind()` is unchanged: the lock refuses the plan's outcome columns (`R/lock.R:85-88`), `.blocked_reasons()` lists a design report with no feasible estimand or a `STOP` verdict (`R/unblind.R:7-11`), `unblind()` refuses those without a written override (`R/unblind.R:94-98`), and `.verify_unblinded()` rechecks the lock, the design report, the unblind hash and the gate before anything is estimated (`R/unblind.R:43-59`). Protected by `test-lock.R:1`, `test-blinding.R:1` and `test-dossier-unblind.R:56`, `:242` and `:260`.

A lock reproduces from its seed. Given the same design data, plan, dictionary, seed and build, `create_analysis_lock()` gives the same hash (the creation time is not hashed, `R/lock.R:141-143`), and every design-stage function gives the same numbers: seeded steps run under `.with_seed()` (`R/utils.R:40-43`), each repetition has its own seed (`R/simulate_design.R:91-92`, `:112`), the net is seeded (`R/learners.R:71-73`) and boosting uses no subsampling (`R/learners.R:77-79`). Protected by `test-simulate-design.R:174` and `:198` and by the reference-values test of step 1, which reruns a fixed design stage and compares it with values stored from 0.5.0.

The numbers do not move. No step of this spec changes an estimate, a standard error, a coverage, a verdict or a selected candidate. The reference-values test is the check; a step that changes a stored number stops until Andrew has seen why.

## Constraints already set

No R process runs while the Rescue.Co design stage runs: code and tests may be written then, but no test, `R CMD check` or render runs until it has finished. Only one R job runs at a time (`2026-10-18-outcome-blind-design-stage.md:133`), and the report's order puts this work after pilot step 3a on the 0.4.1 library, the rehearsal of script 25 on synthetic outcomes and the light outcome-free case-study additions, and before the learner-defaults experiments and the full simulation.

Nothing is installed into `rescueCo/lib` until script 25 has run, and nothing into `simulation/lib` until pilot step 3a has run. The case study stays on 0.5.0 as installed. Tests use `devtools::test()` and the check builds into its own temporary library; neither `rescueCo/lib` nor `simulation/lib` may be on `.libPaths()` when they run.

## Order of work

Step 0 sets the starting point. Step 1 writes the tests and pins (recommendation 4) and fixes the 0.5.0 numbers before any R file changes. Steps 2 and 3 change protected code (recommendations 7 and 8). Steps 4 to 9 are the cuts, in the order that touches each file least often: the time-to-event path (9), the comparators (11), the dead methods and the flag (10), the grading rules (12), `select` and `surfaces$q0` (13), then the merged tests and repeated code (14). Step 10 adds the small items from the literature (17), step 11 the vignettes and `README.md` (15), and step 12 the housekeeping and the release (16).

### Step 0. Starting point (no R)

Andrew records the fingerprint of `rescueCo/lib` in the ledger (report, recommendation 1) and tags the commit it was built from as `v0.5.0` (recommendation 16). The uncommitted edit to `vignettes/cleanTMLE.Rmd` on `vignette-limits` is committed and merged to main, or set aside; that is his call. D1 to D13 are answered and dated in this file. The branch starts from main after that.

Done when the tag exists, the branch exists, and every decision above has a status and a date.

### Step 1. Tests and pins (recommendation 4)

This step adds tests and changes two small things in the package. Sub-steps 1a to 1i run on the unchanged 0.5.0 code; 1j comes last.

1a. Reference values and stored objects. A script `data-raw/make-test-reference.R` **[call]** runs on the branch at its first commit, which is the 0.5.0 code, and saves two fixtures under `tests/testthat/fixtures/`. The first, `reference-0.5.0.rds`, holds for two fixed plans the cross-fitted and single-fit propensity scores, the weight, balance and population tables of `assess_design()`, the numeric columns of `simulate_design()`'s metrics and its verdict, the estimates and intervals of `negative_control_ladder()`, and the `fit_candidate()` rows for every estimand, together with the versions of R, tmle, SuperLearner, glmnet, xgboost, earth and nnet. The two plans are the vignette's plan on `cr_example` (`helper-fixtures.R:60-70`, interval rule) and a small default-library plan on the first 300 rows of `cr_example` with `K = 2`, `V = 2` and two repetitions **[call]**. The second fixture, `v050/`, holds the lock, design, simulation, ladder and design report of the vignette's run, saved with `saveRDS()`. Generating them is an R job and waits its turn in the queue above.

1b. `test-reference.R` checks that the stored 0.5.0 objects still pass `verify_lock()` (with its version warning), `.check_stamp()` for the design, simulation and ladder, and `.verify_dossier()`. It then reruns both plans and compares every stored number to a tolerance of `1e-10`, columns compared by name so that a column a later step removes drops out of the comparison and nothing else does. When the installed dependency versions differ from the recorded ones the comparison skips with a message naming the difference **[call]**, so that a dependency update shows up as a visible skip, not as a failure blamed on a cut.

1c. A direct test of `.fit_Q()` (`R/fit_candidate.R:38-46`) beside `test-fit-candidate.R:11`: with a `"glm"` library and two folds, the predictions for fold 1 under treatment and under control equal those of a logistic regression with the treatment as a column, fitted on fold 2 only.

1d. Cases for glmnet and nnet in the edge check (`R/learners.R:229-237`), beside `test-learners.R:53`: hand-built risk rows in which the penalty sits at the smallest value on the path, and in which the decay sits at the bottom and at the top of its grid, each flagged with the right share.

1e. One test that fits the default library (`.default_library`, `R/plan.R:3`) through `.super_learner()` (`R/learners.R:182-214`) on the 300-row data of `test-learners.R:1-7`: 17 risk rows, none failed, weights on the simplex.

1f. Pin `.nnloglik()` (`R/learners.R:159-180`) to `SuperLearner::method.CC_nloglik()` on the minimized log-loss, not on the weights, which are poorly determined among near-duplicate members. The sketch below uses the `computeCoef(Z, Y, libraryNames, obsWeights, control, verbose)` signature of the installed SuperLearner 2.0-29, with `trimLogit` set to the bound `.nnloglik()` uses.

```r
cc <- SuperLearner::method.CC_nloglik()$computeCoef(
  Z = Z, Y = y, libraryNames = colnames(Z), obsWeights = rep(1, length(y)),
  control = list(trimLogit = 1e-6), verbose = FALSE)
ll <- function(w) .logloss(y, stats::plogis(drop(stats::qlogis(.bound(Z, 1e-6)) %*% w)))
expect_lte(ll(.nnloglik(Z, y)), ll(cc$coef) + 1e-6)
expect_lte(ll(cc$coef), ll(.nnloglik(Z, y)) + 1e-6)
```

1g. Pin the glmnet and xgboost paths of `.cv_setting()` (`R/learners.R:123-147`) on the same inner folds, with a sample size that splits into folds of equal size so that fold-weighted and unweighted means agree. For each elastic-net mixing value, the selected point equals `which.min()` of `glmnet::cv.glmnet(foldid = inner, lambda = path, type.measure = "deviance")$cvm`, and the risk at it equals half of `cvm` there, to a tolerance that allows for glmnet's bounding of fitted probabilities. For one boosting setting, the selected round and its risk equal those of `xgboost::xgb.cv(folds = split(seq_len(n), inner), nrounds = 500)`.

1h. Pin the targeting of the ATE and the ATT (`.target_tmle()`, `R/fit_candidate.R:63-74`) to values stored in a fixture on two fixed datasets, one on which tmle refits the propensity score for the ATT (more than a tenth of the controls fall outside the overlap subset, `R/fit_candidate.R:49-52`) and one on which it does not. Pin the ATE also to a short independent fluctuation with the clever covariate, written in the test:

```r
ate_ref <- function(Y, A, g, Q1, Q0, trunc) {
  g <- pmin(pmax(g, trunc), 1 - trunc)
  H <- A / g - (1 - A) / (1 - g)
  eps <- stats::coef(stats::glm(Y ~ -1 + H + offset(stats::qlogis(ifelse(A == 1, Q1, Q0))),
                                family = stats::binomial()))
  Q1s <- stats::plogis(stats::qlogis(Q1) + eps / g)
  Q0s <- stats::plogis(stats::qlogis(Q0) - eps / (1 - g))
  psi <- mean(Q1s - Q0s)
  ic <- H * (Y - ifelse(A == 1, Q1s, Q0s)) + Q1s - Q0s - psi
  c(psi = psi, se = stats::sd(ic) / sqrt(length(Y)))
}
```

The tolerance against this reference is set after reading how tmle bounds the outcome regression; a gap larger than that is a finding for Andrew before step 2. The existing check at `test-fit-candidate.R:119` stays as the test of the speed-up it protects, not as a pin, because it calls tmle a second time with the same arguments.

1i. Coverage of the overlap-weighted TMLE (`.target_ato()`, `R/fit_candidate.R:86-114`) with the true nuisance functions, on 200 datasets of 2,000 patients from the process at `test-fit-candidate.R:39-55`: coverage of the 95% CI between 0.92 and 0.98. This test and the default-library reference run only when `CLEANTMLE_SLOW_TESTS` is set, and Andrew runs them once before the tag **[call]**.

1j. Require `tmle (>= 2.1.1)` (`DESCRIPTION:28`). In `.tmle_args()` (`R/fit_candidate.R:48-61`), pass the two arguments whose defaults tmle 2.1.1 changed (its `NEWS:1-4`), `target.gwt = FALSE` as the 2.1.1 usage shows it and `min.retain = NULL` **[call]**, so that a later change of default cannot alter results silently, and remove the `formals()` guard at `:59`. The `skip_if_not()` lines at `test-fit-candidate.R:109` and `:120` go, and `test-fit-candidate.R:108` gains an assertion that `target.gwt` is `FALSE`.

```r
.tmle_args <- function(Y, A, X, g, Q1, Q0, truncation, which) {
  args <- list(Y = Y, A = A, W = as.data.frame(X), Q = cbind(Q0, Q1), g1W = g,
               gbound = c(truncation, 1 - truncation), family = "binomial",
               g.SL.library = "SL.glm", target.gwt = FALSE, min.retain = NULL)
  if (which != "ATT") args$evalATT <- FALSE
  args
}
```

No test is removed in this step.

Done when every new test passes on the 0.5.0 code (1a to 1i) and again after 1j with the reference test unchanged, and any pin that fails has gone to Andrew as a finding before step 2 starts.

### Step 2. Record the build and refuse a different one at estimation (recommendation 7)

A new `.build_record()` in `R/lock.R`, after `.pkg_version()` (`:1`), returns the cleanTMLE version, a hash of every object in the package namespace made with the package's own `.hash()`, the source commit when the installer recorded one, the R version and the versions of tmle, SuperLearner, glmnet, xgboost, earth and nnet. `simulate_design()` stores it in its result (`R/simulate_design.R:388-393`), inside the simulation's stamp, so the build that ran the simulation is what the design report records. `design_report()` copies `simulation$build` into the design report as `build`, beside `version` (`R/dossier.R:167`), inside the design-report hash, and warns when the running build differs. `estimate_effect()` (`R/estimate_effect.R:141-195`) gains `build_override = NULL`; after `.verify_unblinded()` (`:143`) it compares the running build with the design report's and stops on any difference unless `build_override` is text, which then goes into the estimate's note. `export_design_log()` writes the build and any override into the JavaScript Object Notation (JSON) log (`R/design_log.R:153-165`). The lock object, `.lock_hash()` and the version warning of `verify_lock()` (`R/lock.R:159-162`) do not change, so every existing lock still passes `verify_lock()`.

```r
.build_record <- function() {
  ns <- asNamespace("cleanTMLE")
  nms <- sort(grep("^\\.__", ls(ns, all.names = TRUE), value = TRUE, invert = TRUE))
  pk <- c("tmle", "SuperLearner", "glmnet", "xgboost", "earth", "nnet")
  list(cleanTMLE = .pkg_version(),
       code = .hash(lapply(mget(nms, envir = ns), function(x)
         if (is.function(x)) deparse(x, control = NULL) else x)),
       commit = utils::packageDescription("cleanTMLE")$RemoteSha %||% NA_character_,
       R = as.character(getRversion()),
       packages = vapply(pk, function(p) as.character(utils::packageVersion(p)), ""))
}
```

The commit is filled only by an installer that records it, such as `remotes::install_github()`; for the case study the ledger record of recommendation 1 carries it. The implementer checks once that `devtools::load_all()` and an installed build give the same code hash; if they do not, tests compare records only within one mode.

Tests added in `test-build.R`: the record has its fields; a design report carries the simulation's build, and editing it fails `.verify_dossier()`; `estimate_effect()` refuses when `.build_record()` is mocked to report another tmle version, names the package in its message, and proceeds with `build_override` text, which the estimate's note and the JSON log carry; a design report from the 0.5.0 fixture, which has no build, is refused at estimation without an override.

Done when the new tests pass, the reference test is unchanged and the stored 0.5.0 lock still passes `verify_lock()`. The paper sentence at `manuscript/methods.qmd:71` can then carry the version condition the report gives under recommendation 5.

### Step 3. One refusal of locks made before 0.5.0 (recommendation 8)

Remove the compatibility guards: `.need_design_fields()` (`R/utils.R:50-57`) and its calls at `R/dossier.R:93`, `R/methods.R:29`, `R/methods.R:148` and `R/simulate_design.R:285`; `.need_dictionary()` (`R/dictionary.R:229-236`) and its calls at `R/assess_design.R:93`, `R/dossier.R:90`, `R/nc_ladder.R:152`, `R/simulate_design.R:284`, `R/unblind.R:47` and `R/unblind.R:81`; the checks for objects from earlier versions in `design_report()` (`R/dossier.R:94-96` and `:101-109`, keeping the stamp and hash checks at `:110-114`) and in `negative_control_ladder()` (`R/nc_ladder.R:155-157`); the branch for a design report without an identification table in `.identification_rows()` (`R/design_log.R:33`) and the branch without `nc_estimand` in `.muntner_log()` (`R/design_log.R:73-79`); the defaults for plans stored without `prognostic` or `max_reps` (`R/plan.R:86`, `R/simulate_design.R:28` and `:325`) and for a lock without a creation time (`R/dossier.R:150`); and the conditional dictionary in `.lock_hash()` (`R/lock.R:15-16`, `:21`).

Put one refusal in their place, at the top of `verify_lock()` (`R/lock.R:153-164`), before the hash is compared, so that a 0.4.x lock gets a clear message instead of a hash mismatch:

```r
.lock_hash <- function(lock) {
  .hash(list(plan = lock$plan, data = lock$data, id = lock$id,
             treatment = lock$treatment, covariates = lock$covariates,
             seed = lock$plan$seed, K = lock$plan$K, version = lock$version,
             dictionary = lock$dictionary))   # same elements, same order
}
# first lines of verify_lock(), after the class check
if (is.null(lock$dictionary) || is.null(lock$version) ||
    utils::compareVersion(lock$version, "0.5.0") < 0)
  stop("This lock was made by cleanTMLE ", lock$version %||% "before 0.5.0",
       ", which has no variable dictionary. Re-lock with this version.", call. = FALSE)
```

Tests removed: `test-lock.R:422` (a lock without a dictionary passes the hash check), `test-dossier-unblind.R:120`, `:327` and `:336`, `test-design-log.R:139` and `:189`, the earlier-version half of `test-nc-ladder.R:487` (keeping its first case, a ladder run on another simulation) and `test-nc-ladder.R:513`, and `test-simulate-design.R:469`. Test added: a lock with version `"0.4.1"` and no dictionary is refused, with the new message, by `verify_lock()`, `load_lock()` and every function that takes a lock.

Done when no `.need_dictionary`, `.need_design_fields` or "earlier cleanTMLE" remains in `R/`, the stored 0.5.0 objects verify, and the reference test is unchanged.

### Step 4. Cut the time-to-event path (recommendation 9, after D4)

Remove `.sl_equivalents()` and `.estimate_tte()` (`R/estimate_effect.R:78-121`), the time-to-event branch of `estimate_effect()` (`:160-170`) and its roxygen lines (`:128-129`); the time-to-event arguments, documentation, parsing and fields of `analysis_plan()` (`R/plan.R:95-96`, `:131-132`, `:175-176`, `:191-203`, `:209-210`, and `outcome_type`, `target_time` and `hazards` at `:305-306`); the time-to-event half of `.check_outcome_coding()` (`R/unblind.R:26-37`) and its documentation (`:70-73`); the time-to-event scope note (`R/simulate_design.R:385-387`); `hazards` and `outcome_type` in the design report (`R/dossier.R:10-11`, `:151`); and `outcome_type` in `print.cr_plan()` (`R/methods.R:5`). `analysis_plan()` refuses an outcome of more than one column, and `.check_outcome_coding()` refuses a stored plan that has one, each with a message that names the wishlist **[call]**. Drop concrete and data.table from Suggests; SuperLearner stays there for the pin of step 1 **[call]**. Add the time-to-event estimand to `WISHLIST.md`, in its section on parity with `causalRisk` (`:52-55`).

Tests removed: `test-estimate-effect.R:88`, `test-dossier-unblind.R:289`, and the time-to-event cases at `test-plan.R:18-21`. Test added: the two-column refusal.

Done when no `concrete`, `hazards`, `target_time` or `outcome_type` remains in `R/`, `analysis_plan()` has 20 arguments, and the reference test is unchanged.

### Step 5. Remove the comparators (recommendation 11, after D1)

Remove `.syntactic()`, `.estimate_weighting()` and `.estimate_matching()` (`R/estimate_effect.R:19-76`), the `method` argument and its branches (`:141`, `:144`, `:177`, `:186-194`) and their documentation (`:125-139`); `estimate_effect()` becomes `estimate_effect(unblinded, estimand = NULL, build_override = NULL)`. The estimate keeps its `method` (always `"tmle"`) and `stress_tested` fields so that `as.data.frame()` keeps its columns **[call]**. Change the advice in `.no_candidate_message()` (`:13-17`), which names a comparator, and remove the comparator note in the design report (`R/dossier.R:141`), the `globalVariables()` call (`R/cleanTMLE-package.R:5-7`), the comparator sentence in `DESCRIPTION` (`:12-13`) and the weighting chunk of the main vignette (`vignettes/cleanTMLE.Rmd:349-354`). Replace the one cobalt call (`R/assess_design.R:112-113`, `:121-135`) with a base-R standardized difference:

```r
.max_smd <- function(X, A, w) {
  m <- function(a) colSums(X[A == a, , drop = FALSE] * w[A == a]) / sum(w[A == a])
  v <- function(a) apply(X[A == a, , drop = FALSE], 2, stats::var)
  d <- abs(m(1) - m(0)) / sqrt((v(1) + v(0)) / 2)
  d <- d[is.finite(d)]
  if (length(d)) max(d) else NA_real_
}
```

The pin below decides the variance used for binary columns: it must match `cobalt::bal.tab(s.d.denom = "pooled", binary = "std")`. WeightIt leaves Imports, cobalt moves from Imports to Suggests, and MatchIt leaves Suggests.

Tests removed: `test-estimate-effect.R:32`, `:42`, `:113`, `:121` and `:144`; the weighting and matching calls inside `test-estimate-effect.R:49`. Test changed: `test-assess-design.R:43` triggers the warning by mocking `.max_smd()` instead of cobalt. Test added: `.max_smd()` equals the largest absolute `Diff.Adj` of `cobalt::bal.tab()` for each estimand's weights on `cr_example`, skipped when cobalt is absent.

Done when no `WeightIt`, `MatchIt` or `cobalt::` remains in `R/`, and the reference test is unchanged, including the balance table, to `1e-8`. Outside the package, `manuscript/methods.qmd:120` and the comparator question of recommendation 1 follow D1.

### Step 6. Cut the dead methods, the E-value summary and the flag (recommendation 10, after D5)

Remove the `IMPLAUSIBLE` line of `print.cr_estimate()` (`R/methods.R:91`), `.evalue()`, `summary.cr_estimate()` and its print method (`:95-123`), `tidy.cr_estimate()` (`:133-138`), the `as.data.frame()` methods for simulations and ladders (`:140-144`) and `plot.cr_simulation()` (`:157-166`); `.implausibility_check()` (`R/fit_candidate.R:131-148`), its call (`:175`) and the `guard` argument and the two flag columns of `.estimate_row()` (`:150-162`); and the flag fields of `.new_estimate()` (`R/estimate_effect.R:9`). Regenerate `NAMESPACE`. Drop EValue and generics from Suggests, and correct `WISHLIST.md:62-63`, which says the risk ratio is computed only for the E-value.

Tests removed: `test-estimate-effect.R:80` and `:134`, `test-fit-candidate.R:85`. The reference test compares columns by name, so the two flag columns drop out and nothing else changes.

Done when `NAMESPACE` lists 10 S3 methods and the reference test is unchanged. Recommendation 26 must check the column set of `.estimate_row()` when the simulation moves to 0.5.x, since `simulation/R/key.R:16` and `simulation/R/map.R:53` read it.

### Step 7. Keep only the interval rule (recommendation 12, after D6)

Change `analysis_plan()` to accept only `rule = "interval"` and refuse `"point"` and `"ci"` with a message (`R/plan.R:146-153`, `:266-270`); reduce `.nc_grade()` to the interval reading and drop the `in_band` column (`R/nc_ladder.R:7-26`); reduce `.nc_rule_text()` to one rule (`R/design_log.R:8-27`). The template states the rule a second time (`inst/templates/dossier.qmd:150-170`); `design_report()` stores the text of `.nc_rule_text()` in the design report and the template prints it, so the sentence has one source **[call]**.

Tests changed: the fixture's default `nc_rule = "point"` becomes `"interval"` (`helper-fixtures.R:32`), and the fixtures that need a fixed `GO` or `STOP` are rewritten with controls whose effects give the same verdict under the interval rule (about ten, `01-parsimony.md:277`); `test-nc-ladder.R:95` keeps its precedence cases only; `test-nc-ladder.R:246-247` and `:360` drop `in_band`; `test-design-log.R:121` and `:152` merge into one; `test-plan.R:80` checks that the interval rule is accepted and the others refused. Tests removed: `test-nc-ladder.R:276` and `:292`.

Done when the reference ladder (interval rule) is unchanged and a stored 0.5.0 plan with `rule = "interval"` still validates. This touches the grading that gates `unblind()`; the interval logic itself does not change.

### Step 8. Remove `select` and `surfaces$q0` (recommendation 13, after D7)

Remove the argument, documentation, check and field of `select` (`R/plan.R:156-157`, `:181`, `:293-294`, `:310`) and its uses in `.decision_candidates()` (`R/simulate_design.R:205-214`) and `.sim_verdict()` (`:238-242`); remove `surfaces$q0` from the documentation and `.n_surfaces()` (`R/plan.R:84-88`, `:126-127`), `.surfaces()` (`R/simulate_design.R:17-19`), the lock's length check (`R/lock.R:127-128`) and the design report's plan summary (`R/dossier.R:5`, `:7-8`), and its sentence in the second vignette (`vignettes/plan-and-dossier.Rmd:112-115`). Closures stay hashable for `surfaces$custom` (`R/utils.R:3-30`), so recommendation 34 stays open.

Tests changed: `test-lock.R:61`, which shows that a function in the plan changes the hash, uses `surfaces$custom` closures instead of `select`; the `select` case of `test-simulate-design.R:83` (at `:110`) goes.

Done when `analysis_plan()` has 19 arguments and the reference test is unchanged.

### Step 9. Merge repeated tests and give repeated code one source (recommendation 14)

In the code: replace `.is_text()` (`R/unblind.R:1-3`) with `.has_text()` (`R/plan.R:19`) at `R/unblind.R:55`, `:85` and `:95`; fold `.role()` (`R/design_log.R:1-4`) into its one caller (`:43`); replace `.deparse1()` (`R/dossier.R:1`) with base `deparse1()` at `:5-11`, which changes only the text of new design reports, while `.hashable()` (`R/utils.R:3-30`) keeps its own `deparse()` call because it feeds the hashes; turn the five stray-role blocks of `.check_dictionary()` (`R/dictionary.R:142-174`) into one loop with the same messages; keep one copy of the before, at and after vocabulary (`R/plan.R:54`, `R/dictionary.R:10`); drop the arguments `surfaces_run` and `libraries_run` of `.one_rep()` (`R/simulate_design.R:80-89`), keeping `pairs`; drop the outputs `coverage_mcse` and `rmse_mcse` (`R/simulate_design.R:159`, `:162`). Replace `.wilson()` (`R/simulate_design.R:128-136`, called at `:155`) with `stats::prop.test(x, n, correct = FALSE)$conf.int`, keeping the guard for zero repetitions, only if a pin over every count of successes, for every number of trials up to a thousand, agrees to `1e-12` **[call]**; otherwise `.wilson()` stays, with the pin as its test.

In the tests, follow `01-parsimony.md:146-238`: one render of a rich fixture in place of the six Quarto render tests (`test-dossier-unblind.R:149`, `:179`, `:200`, `:215`, `:372`, `:448`); one test per repeated check (the stamp checks at `test-assess-design.R:34`, `test-lock.R:50` and `:73`, `test-nc-ladder.R:197`; the control that is also a covariate at `test-lock.R:106` and `:140`, `test-nc-ladder.R:615`, `test-dossier-unblind.R:395`); cut the tests that only check a class or restate the code (`test-blinding.R:8`, `test-assess-design.R:159`, `test-dossier-unblind.R:235`, `test-utils.R:18`, `test-simulate-design.R:415`); give `test-design-log.R:59` a fixture that produces `STOP`, so that its `skip_if_not()` at `:70` goes; rewrite the fold test (`test-simulate-design.R:383`); update the mock of `.one_rep()` (`test-simulate-design.R:283-286`) and the case at `:183`; delete the empty `tests/testthat/_snaps/`.

Every test named under "Guarantees that must not change" stays, as do the new tests of steps 1 to 3.

Done when the tests pass, the reference test is unchanged, and the test lines are near the estimate in the size table below.

### Step 10. Four small items from the literature (recommendation 17, after D10)

Add the equipoise share to `assess_design()` (`R/assess_design.R:104-111`, result list `:138-142`), computed on the cross-fitted score as the simulation's comparator computes it (`simulation/R/map.R:87-93`), shown by `print.cr_design()` (`R/methods.R:28-36`), carried into the design report's design fields (`R/dossier.R:152-153`), the JSON log (`R/design_log.R:157-159`) and the template, and never read by a decision:

```r
pref <- stats::plogis(stats::qlogis(g) - stats::qlogis(mean(A)))
share_equipoise <- mean(pref >= 0.3 & pref <= 0.7)
```

Add one sentence on what withholding the outcome does not guarantee to the tamper-evidence paragraph of `README.md` (`:45-47`) and to the notes `design_report()` writes (`R/dossier.R:141-146`): choices made on the covariates can still bias the effect when the outcome is predictable from them, citing Fisher 2020 (`arXiv:2007.02514`) once sci-cite has added it. Add a short `README.md` paragraph naming the study diagnostics of the Observational Health Data Sciences and Informatics network and the REFINE2 tool, with one sentence each on what they fix and what they read. Add a worked example to the candidates section of `vignettes/plan-and-dossier.Rmd` (`:117-125`) of a truncation set from the sample size (Gruber and colleagues 2022, `doi:10.1093/aje/kwac087`), which `analysis_plan()` already accepts (`R/plan.R:227-230`), with the bound read from that paper first. The registration format stays on `WISHLIST.md`, which lists it (`:46`).

Test added: the equipoise share equals a direct computation on the reference design. The reference test compares the stored design fields by name, so the new field does not disturb it.

Done when the tests pass and the two sentences are in `README.md` and the design-report notes.

### Step 11. Vignettes and `README.md` (recommendation 15, after D8)

In the main vignette: the opening proposed in `10-proposed-text.md`, section 1, which says what a user supplies and gets back, defines time zero, states the boundary of the simulated family, says which parts are in-house and which go to tmle, and says that the design data's freedom from the outcome rests on the steward; pairs D and E of `05-writing.md`; two candidates in the example plan (`vignettes/cleanTMLE.Rmd:90`) so that selection is visible; one refusal of the lock shown with `tryCatch()`; the forecast wording at `:232-233`; and the term set of D8. In the second vignette: what to do when cells stay unresolved at the repetition cap or no estimand is feasible; the re-locking policy of D8; that `approved_by` is a label and not an authentication (`R/unblind.R:85-86`); and one wording for the bias tolerance (`vignettes/plan-and-dossier.Rmd:75-77`). In `README.md` (`:3-7`) and the `DESCRIPTION` description (`:8-13`): the boundary of the simulated family. Copy into `vignettes/references.bib` the entries for every key the new text cites.

Done when both vignettes have been rendered once, after the design stage has finished, with no error and with the rendered text read.

### Step 12. Housekeeping and the release (recommendation 16)

Delete `attic/` and its line in `.Rbuildignore` (`:10`, after D9). Bump `DESCRIPTION:3` to the version of D11. Write `NEWS.md` (draft below). Replace the rendered copies in `docs/vignettes-0.5.0/`, which were built from the uncommitted working tree (`05-writing.md:63`), with renders from the release commit under a folder named for the release. Andrew runs `R CMD check` once, after the Rescue.Co design stage has finished and never while it runs, with `CLEANTMLE_SLOW_TESTS` set for one run of the slow tests. Merge, then tag.

Done when `R CMD check` reports no error, warning or note, the reference test passes with matching dependency versions, the stored 0.5.0 objects verify, the counts in the size table have been recounted from `NAMESPACE`, `DESCRIPTION`, `R/` and `tests/testthat/`, and the supplement records the commit the case study ran on.

## Expected size afterwards

The first two columns are the report's size table (`REPORT.md`, "Projected size of the package"). The last column is this spec's estimate for the 0.5.x release: it keeps the export `save_lock()` and keeps jsonlite in Imports (both wait until after the paper), keeps SuperLearner in Suggests for the pin, moves cobalt to Suggests, and adds the code of steps 1, 2 and 10, which the report sized in days rather than lines. Its line and function counts are approximate.

| Measure | 0.5.0 today | After every package recommendation (report) | After this spec |
|---|---|---|---|
| User-facing functions | 12 exports | 11 | 12 |
| S3 methods | 16 | 10 | 10 |
| Internal functions | 102 | 88 | 90 |
| Arguments of `analysis_plan()` | 22 | 19 | 19 |
| Lines of R | 3,068 | 2,670 | 2,720 |
| Lines of tests | 3,457 | 3,030 | 3,250 |
| Imports entries | 13 (11 packages) | 11 (9 packages) | 11 (9 packages) |
| Suggests entries | 12 | 7 | 8 |

Imports after this spec: digest, earth, ggplot2, glmnet, jsonlite, nnet, stats, tmle, utils, withr and xgboost. Suggests: cobalt, future, future.apply, knitr, quarto, rmarkdown, SuperLearner and testthat. The conclusion of the report stands: by user-facing surface the package is about half the size of `causalRisk`, and most of what remains is code the outcome-blind design relies on.

## Release

One branch, `simplify-0.5.1` (D11, D12), one commit per step, each commit message naming its step and recommendation. The tests run after each step once the design stage has finished; `R CMD check` runs once at the end, by Andrew, after the Rescue.Co design stage has finished and never during it. No build of this branch is installed into `rescueCo/lib` before script 25 has run or into `simulation/lib` before pilot step 3a has run.

Draft `NEWS.md` entries, one per step:

```markdown
# cleanTMLE 0.5.1

* Numbers are unchanged from 0.5.0: a stored run of the design stage on 0.5.0
  is reproduced by the tests.
* New tests pin the in-house code to established packages: the super learner's
  weights to SuperLearner's convex log-likelihood combination, the elastic-net
  and boosting paths to cv.glmnet() and xgb.cv(), the targeting of the ATE and
  the ATT to stored values and an independent fluctuation, and the coverage of
  the overlap-weighted TMLE.
* tmle (>= 2.1.1) is required, and the targeting arguments whose defaults
  tmle 2.1.1 changed are passed explicitly.
* simulate_design() records the build (the code, R and the fitting packages),
  the design report carries it, and estimate_effect() refuses a different
  build unless `build_override` gives a written reason.
* Every function refuses a lock made before 0.5.0 with one message; the
  separate checks for objects from earlier versions are gone.
* The time-to-event path is removed (see WISHLIST.md).
* The weighting and matching comparators are removed, with the `method`
  argument of estimate_effect(); WeightIt and MatchIt are no longer used, and
  cobalt is needed only to test the balance table.
* summary() with the E-value, tidy(), the as.data.frame() methods for
  simulations and ladders, plot() for simulations and the implausibility flag
  are removed.
* Negative controls are graded by the interval rule only; "point" and "ci" are
  refused.
* `select` and `surfaces$q0` are removed from analysis_plan().
* assess_design() reports the share of patients in equipoise (preference score
  between 0.3 and 0.7); it gates nothing.
```

## Calls for Andrew

1. Reference values and stored 0.5.0 objects as test fixtures, generated by `data-raw/make-test-reference.R` on the branch's first commit, with the comparison skipped, visibly, when dependency versions differ from the recorded ones (step 1a, 1b).
2. The default-library reference runs on the first 300 rows of `cr_example`, with `K = 2`, `V = 2` and two repetitions, to keep it to about a minute (step 1a).
3. The coverage test of the overlap-weighted TMLE and the default-library reference run only with `CLEANTMLE_SLOW_TESTS` set, once before the tag (step 1i).
4. `min.retain = NULL` is passed explicitly along with `target.gwt = FALSE`; the report names only `target.gwt`, and audit 06 names both as defaults tmle 2.1.1 changed (step 1j).
5. SuperLearner stays in Suggests, because the pin of `.nnloglik()` calls it; recommendation 9 would have dropped it (step 4).
6. The build is recorded in the simulation object and copied into the design report, so that a design report built on one build from a simulation run on another is caught; the override argument is named `build_override` (step 2).
7. A plan with a two-column outcome is refused with a message naming the wishlist, in `analysis_plan()` and in `.check_outcome_coding()` (step 4).
8. The estimate keeps its `method` and `stress_tested` fields, constant after step 5, so that `as.data.frame()` output keeps its columns for code that reads it (step 5).
9. The grading-rule sentence is stored in the design report and printed by the template, which gives it one source (step 7).
10. `.wilson()` is replaced by `prop.test()` only if they agree to `1e-12` on every count for every number of trials up to a thousand (step 9).

## Writing check

The check was run from the project root in Git Bash, `python "C:/Users/andre/OneDrive/Documents/scientific-writing-skills/checks/check.py" "audit/2026-10-22/09-spec-package-simplification.md" --min-level warn`, and again after the fixes. The first run flagged two words Andrew does not use (a verb for combining branches, twice, and one for passing a hash check, three times), an algorithm name read as a dash code, a package name read as an identifier, two undefined abbreviations (the log format and the identifier of a paper), and typed numbers that are not results (line numbers of the vignette edit, test tolerances, a preprint identifier and a count of trials). Each was reworded or put in backticks.

What remains is expected in a spec. Every error flags code formatting in prose: a function, argument, file, path with line numbers, package or test name in backticks, which the spec needs to say what each step touches. Six warnings read the decision labels as undefined abbreviations; they follow the house format of the specs in `docs/superpowers/specs/`. One warning concerns the project configuration, not this file: the overlap source in the unblinded analysis's repository is not found, and I did not look for it, so overlap with that manuscript was not checked. The code sketches and the draft `NEWS.md` entries sit in fenced blocks, which the check does not read. No dash, banned word, bold lead-in or glossary label remains.

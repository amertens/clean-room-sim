# Parsimony audit of cleanTMLE 0.5.0

Audit of 2026-10-22 on the package source at commit `2e96f98` (the merge of release/0.5.0). I read every file in `cleanTMLE/R/`, the tests, `NAMESPACE`, `DESCRIPTION`, `.Rbuildignore`, the 0.5.0 and 0.4.x parts of `NEWS.md`, `WISHLIST.md`, and listed `attic/` and `data-raw/`. I ran no R code.

cleanTMLE 0.5.0 has 12 exported functions, 16 S3 methods and 102 internal functions in 17 R files, 13 reference pages, two vignettes, and test files that hold 13% more lines than the R code. By surface it is already smaller than `causalRisk`, whose function reference documents 93 topics (29 under five headings meant for users and 64 labelled as used internally) and 8 datasets; I could not size its code or dependencies from the PDFs. The weight of cleanTMLE sits behind its exports, in a hand-written cross-validated learner with an edge check, the dictionary and identification checks, the adaptive simulation and the hashing. Most of that is what the outcome-blind design stage relies on and should stay. The reductions I can support are about an eighth of the lines of R and of tests, one export, six S3 methods, 14 internal functions, two packages in Imports and a net five in Suggests; the size table near the end gives the counts. More than half of the saving in R comes from four features that the Rescue.Co scripts never call: the time-to-event path, the weighting and matching comparators, two of the three negative-control grading rules, and the E-value summary.

Every change below has to leave `verify_lock()`, `unblind()` and `estimate_effect()` able to read the lock and design-report files that the running design stage writes. A plan field removed from the code stays as a harmless extra field in stored objects, but the element list and order inside `.lock_hash()` and the hashing helpers must not change.

## Functions

The design report is the object `design_report()` returns, called `dossier` in the tables and in the code. Callers list the call sites in `R/` by file and line, then the number of lines that mention the name in the tests, in the vignettes and the template, in the Rescue.Co re-execution scripts (`rescueCo/scripts/2*.R`) and in the paper simulation (`simulation/R`). The counts come from a text search, so they include comments and messages. "Protected" in the Risk column marks code the outcome-blind design stage relies on: the shared fitter, cross-fitting, the edge check, the lock and design-report hashes with `verify_lock()`, and the blinding in `unblind()`.

| Function | File:line | Exported | Purpose | Callers | Recommendation | Reason | Risk |
|---|---|---|---|---|---|---|---|
| `.estimand_weights` | assess_design.R:1 | no | Inverse-probability weights for each estimand from the cross-fitted score | assess_design.R:115, :122; tests 1 | Keep | Feeds the weight and balance tables. `WeightIt::get_w_from_ps` gives the same weights for ATE, ATT and ATO, but only helps if WeightIt stays. | Low |
| `.ess` | assess_design.R:9 | no | Effective sample size of a weight vector | assess_design.R:117 | Keep | Two lines. `WeightIt::ESS` is the same formula, same condition as above. | Low |
| `.c_statistic` | assess_design.R:12 | no | Area under the ROC curve of the propensity score, by ranks | assess_design.R:140; tests 1 | Keep | Tested against a pairwise count. | Low |
| `.small_trim_group` | assess_design.R:27 | no | True when one to four patients lie on either side of the trim band | assess_design.R:46, :57, :64 | Keep | Small-cell rule that the dossier tables rely on. | Low |
| `.std_shift` | assess_design.R:30 | no | Standardized shift of a mean from the cohort mean | assess_design.R:49, :75 | Keep | Two callers. | Low |
| `.populations` | assess_design.R:35 | no | Covariate means and shifts in each estimand's target population | assess_design.R:136; tests 5; script 23 reads the result from the dossier | Keep | The case study applies its population gate by hand from this table (spec 2026-10-20, section 5). | Low |
| `.trimmed_profile` | assess_design.R:59 | no | Counts, and covariate means of removed and kept patients | assess_design.R:137; tests 6; script 23 | Keep, drop its `kept` mean | The kept mean repeats the `trimmed_ATE` column of `.populations` (:41-48), and both apply the same suppression rule. | Low |
| `assess_design` | assess_design.R:91 | yes | Cross-fits the propensity score, grades overlap, builds weights, balance, populations, the edge check and the single-fit score of the simulated world | tests 41; vignettes 1; scripts 2; simulation 4 | Keep | Core design-stage verb. | Protected: cross-fitting, edge check, shared fitter |
| `.role` | design_log.R:1 | no | Reads one role from the plan | design_log.R:43 | Merge into `.muntner_log` | One caller, three lines. | None |
| `.nc_rule_text` | design_log.R:8 | no | Sentence stating the negative-control grading rule | design_log.R:84; tests 3 | Keep, as the only copy | The same three sentences are written again in `inst/templates/dossier.qmd:154-168`. It shrinks to one variant if only the interval rule stays. | Low |
| `.identification_rows` | design_log.R:31 | no | Turns flagged identification rows into log rows | design_log.R:54 | Keep, remove the NULL branch (:30-33) | The branch serves dossiers built before 0.5.0. | Low |
| `.muntner_log` | design_log.R:41 | no | Builds the decision log in the Muntner columns | design_log.R:148 | Keep, remove the branch for a dossier without `nc_estimand` (:74-79) | Same reason. | Low |
| `export_design_log` | design_log.R:133 | yes | Writes the decision log as CSV (Muntner columns) or JSON | tests 22; vignettes 2; scripts 1 (25_estimate.R:38, Muntner format) | Keep, move jsonlite to Suggests | The case study uses only the Muntner format. The README leads with JSON (README.md:35). | Touches the record that carries both hashes |
| `.dictionary_stop` | dictionary.R:14 | no | `stop()` without the call | dictionary.R, 27 sites | Keep | Removes repetition. | None |
| `.q` | dictionary.R:16 | no | Back-quotes names for messages | dictionary.R, 5 sites | Keep | Tiny. | None |
| `.dictionary_coerce` | dictionary.R:21 | no | Puts the dictionary in a plain canonical form so the hash ignores its type and row order | dictionary.R:64 | Keep | The stored dictionary is hashed. | Protected: lock hash |
| `.kept_columns` | dictionary.R:52 | no | Lists the columns the lock keeps | dictionary.R:101; lock.R:134 | Keep | One definition shared by the check and the lock. | Low |
| `.check_dictionary` | dictionary.R:63 | no | Applies every dictionary rule | lock.R:130 | Keep, turn the five stray-role blocks (:142-174) into one loop | 165 lines, the longest function. The refusals are the contract and have tests. | Medium |
| `.need_dictionary` | dictionary.R:231 | no | Refuses a lock that has no dictionary | assess_design.R:93; dossier.R:90; nc_ladder.R:152; simulate_design.R:284; unblind.R:47, :81; tests 1 | Cut | Only a 0.4.x lock lacks a dictionary, and 0.5.0 is the first public version (commit 2e96f98 message). | Low |
| `.deparse1` | dossier.R:1 | no | One-line text of a function or formula | dossier.R:5, :6, :9, :11 | Replace with base `deparse1()` | Base R has it since 4.0.0 and the package needs R 4.1. | Low: line breaks in dossier text may change |
| `.plan_summary` | dossier.R:3 | no | Makes the plan printable for the dossier | dossier.R:148; tests 2 | Keep | Shrinks when `select`, `hazards` and `q0` go. | Low |
| `.flat` | dossier.R:23 | no | One-line text of a plan field | dossier.R:26, :49; tests 5 | Keep | Used for every row of the plan table. | Low |
| `.plan_table` | dossier.R:47 | no | The prespecified plan as a table | dossier.R:149; tests 2 | Keep | Shown in the dossier. | Low |
| `.dossier_hash` | dossier.R:52 | no | Hash of the dossier without its own hash | dossier.R:60, :169; tests 3 | Keep | Part of the approval chain. | Protected: dossier hash |
| `.verify_dossier` | dossier.R:57 | no | Recomputes and compares the dossier hash | design_log.R:145; unblind.R:48, :82; tests 1 | Keep | Part of the approval chain. | Protected: dossier hash |
| `design_report` | dossier.R:88 | yes | Checks every stamp, walks the ladder, builds and hashes the dossier, optionally renders it | tests 42; vignettes 2; scripts 1 | Keep, remove the 0.4.x checks (:93-96, :102-110) | Those checks cannot fire on a 0.5.0 object. | Protected: dossier hash (the cut touches the build, not the hash code) |
| `.flatten_text` | dossier.R:177 | no | Removes line breaks from table cells before Quarto | dossier.R:183, :199; tests 2 | Keep | Needed for pipe tables. | Low |
| `.render_dossier` | dossier.R:187 | no | Renders the template with Quarto | dossier.R:170 | Keep | The only route to the HTML or Word dossier. | Low |
| `.new_estimate` | estimate_effect.R:1 | no | Builds a `cr_estimate` | estimate_effect.R:166, :182, :191 | Keep | Loses the `stress_tested` and `method` fields if the comparators go. | Low |
| `.no_candidate_message` | estimate_effect.R:13 | no | Message for an estimand with no selected candidate | estimate_effect.R:164, :178 | Keep | Tested through its error. | None |
| `.syntactic` | estimate_effect.R:19 | no | Makes column names syntactic | estimate_effect.R:25, :54, :89 | Cut with the comparators and the time-to-event path | The TMLE path does not use it. | Low |
| `.estimate_weighting` | estimate_effect.R:24 | no | Weighting comparator through WeightIt | estimate_effect.R:186; tests 1 | Cut (decision for Andrew) | Not stress-tested, not called by the case study (25_estimate.R:36 uses the default), and it keeps WeightIt in Imports for 25 lines. Editing needed in `manuscript/methods.qmd:120` and `vignettes/cleanTMLE.Rmd:349-353`. | Low |
| `.estimate_matching` | estimate_effect.R:50 | no | ATT matching comparator through MatchIt | estimate_effect.R:187; tests 1 | Cut | ATT only, repeats the extraction code of `.estimate_weighting` (:40-47 and :68-75), and needs a `globalVariables("weights")` workaround (cleanTMLE-package.R:5-7). | Low |
| `.sl_equivalents` | estimate_effect.R:78 | no | Maps learner names to SuperLearner wrappers for concrete | estimate_effect.R:103 | Cut with the time-to-event path | Only the time-to-event path uses it. | Low |
| `.estimate_tte` | estimate_effect.R:83 | no | Time-to-event estimate through concrete | estimate_effect.R:165; tests 1 | Cut (decision for Andrew) | The one path where the estimator run after unblinding is not the one the simulation evaluated (note at :168). The case study, vignettes and paper do not use it, and the spec moved parity features to the wishlist (spec 2026-10-18, release scope; WISHLIST.md:54-55). | Low |
| `estimate_effect` | estimate_effect.R:141 | yes | Rechecks the unblinded object, then estimates the dossier's estimand with the dossier's candidate | tests 25; vignettes 5; scripts 2 | Keep, becomes TMLE only | The `method` argument and the time-to-event branch (:160-170) go with the cuts above. | Protected: shared fitter and the unblind checks |
| `.design_matrix` | fit_candidate.R:5 | no | Model matrix from the locked covariates | assess_design.R:95; estimate_effect.R:89, :172; nc_ladder.R:179; simulate_design.R:288; tests 11; simulation 2 | Keep | The paper simulation calls it (simulation/R/key.R:10). | Protected: shared fitter |
| `.cross_fit` | fit_candidate.R:13 | no | Runs K-fold cross-fitting | fit_candidate.R:32, :40 | Keep | One driver for both nuisances. | Protected: cross-fitting |
| `.fit_g` | fit_candidate.R:31 | no | Cross-fitted propensity score | 7 sites in R; tests 3; simulation 1 | Keep | Used by the design, the simulation and the analysis. | Protected: cross-fitting, shared fitter |
| `.fit_Q` | fit_candidate.R:38 | no | Cross-fitted outcome regression with the treatment as a column | fit_candidate.R:167; tests 0 | Keep, add a direct test | No test calls it. The open learner-defaults item changes it (spec 2026-10-18, lines 145-147). | Protected: shared fitter |
| `.tmle_args` | fit_candidate.R:48 | no | Argument list for `tmle::tmle` | fit_candidate.R:67; tests 2 | Keep | Thin wrapper. | Low |
| `.target_tmle` | fit_candidate.R:63 | no | Runs `tmle::tmle` on supplied nuisance estimates | fit_candidate.R:118, :119, :122; tests 1 | Keep | Thin wrapper. | Protected: shared fitter |
| `.ato_eif` | fit_candidate.R:76 | no | Influence function of the overlap-weighted effect | fit_candidate.R:107 | Keep | Part of the ATO estimator. | Low |
| `.target_ato` | fit_candidate.R:86 | no | Targeting loop for the overlap-weighted effect | fit_candidate.R:127; tests 1 | Keep | The only estimator written in this package. A search of the export lists of the installed packages for overlap-weighted targeting found none. | Medium: tested against the truth |
| `.target` | fit_candidate.R:116 | no | Dispatches on the estimand | fit_candidate.R:169; tests 5 | Keep | Four-way switch. | Low |
| `.implausibility_check` | fit_candidate.R:131 | no | Flags a sign that differs from the crude difference, a ratio above 5, or a difference above the largest crude risk | fit_candidate.R:175; estimate_effect.R:44, :72; tests 2 | Cut (decision for Andrew) | A sign change from the crude difference is what adjustment for confounding can produce, and no decision reads the flag. Only the print method shows it (methods.R:91). | Low |
| `.estimate_row` | fit_candidate.R:150 | no | Builds one result row | fit_candidate.R:171, :174, :175; simulate_design.R:77 | Keep | Shared by the simulation and the analysis. | Low |
| `.fit_and_target` | fit_candidate.R:164 | no | Fits both nuisances once and targets every estimand and truncation | fit_candidate.R:194; simulate_design.R:112; tests 4; simulation 2 | Keep | The shared fitter. | Protected: shared fitter |
| `fit_candidate` | fit_candidate.R:193 | no | Adapter from a candidate record to `.fit_and_target` | estimate_effect.R:180; nc_ladder.R:229; tests 11 | Keep | The only unexported function without a leading dot. Rename later. | Protected: shared fitter |
| `.identification_table` | identification.R:13 | no | Table of what the plan and dictionary state about time zero, eligibility and timing | dossier.R:149; tests 2; scripts 1 (22_design_stage.R:190, through `:::`) | Keep | The case-study check mode calls it. | Low |
| `.logloss` | learners.R:14 | no | Log-loss | learners.R:108, :114, :171; tests 1 | Keep | Basis of the cross-validated risk. | Protected: edge check |
| `.setting` | learners.R:19 | no | Record for one learner setting | learners.R:30-35 | Keep | | Protected: edge check |
| `.settings` | learners.R:27 | no | Expands learner names into settings (17 for the default library) | learners.R:185; tests 2 | Keep | The nnet decay grid in it is the subject of the open learner-defaults item. | Protected: edge check |
| `.prepare_setting` | learners.R:40 | no | Fixes the glmnet penalty path | learners.R:128 | Keep | | Protected: edge check |
| `.xgb_range` | learners.R:49 | no | Iteration range for xgboost 1.x versus 2.x | learners.R:98 | Keep | Both major versions are in use (the manuscript reports 1.7.8.1, methods.qmd:213). | Low |
| `.learner_fit` | learners.R:53 | no | Fits one setting | learners.R:133, :144, :187; tests 4 | Keep | | Protected: edge check |
| `.learner_predict` | learners.R:91 | no | Predicts from one fit at the chosen point of the path | learners.R:114, :142, :145, :190; tests 1 | Keep | | Protected: edge check |
| `.learner_path_loss` | learners.R:102 | no | Validation loss along the path | learners.R:135 | Keep | Gives the cross-validated risk of every setting that the edge check reads. | Protected: edge check |
| `.inner_folds` | learners.R:117 | no | Inner folds that keep copies of a row together | learners.R:193; tests 2 | Keep | | Protected: cross-fitting |
| `.cv_setting` | learners.R:123 | no | Cross-validates one setting and keeps its position, path length and risk | learners.R:195 | Keep | | Protected: edge check |
| `.risk_row` | learners.R:149 | no | One row of per-setting risk | learners.R:188, :202, :203 | Keep | | Protected: edge check |
| `.nnloglik` | learners.R:159 | no | Convex weights that minimize log-loss on the logit scale | learners.R:206; tests 2 | Keep | | Low |
| `.super_learner` | learners.R:182 | no | The cross-validated super learner | fit_candidate.R:22; tests 5 | Keep | I know of no package that returns every setting's cross-validated risk along a path, which the edge check needs. This is an inference from the packages' documented interfaces; I did not test SuperLearner wrappers. | Protected: edge check, shared fitter |
| `.edge_check` | learners.R:216 | no | Share of fits whose chosen setting sits on a grid edge | assess_design.R:142; simulate_design.R:390; tests 2 | Keep, add glmnet and nnet tests | The test covers xgboost only (test-learners.R:53-65), and the nnet bottom edge is the finding behind the open learner-defaults item. | Protected: edge check |
| `.pkg_version` | lock.R:1 | no | Version string | lock.R:139, :159, :161; dossier.R:167; tests 1 | Keep | | Low |
| `.check_nc_not_covariate` | lock.R:7 | no | Refuses a negative control that is also a covariate | lock.R:115; dossier.R:91; nc_ladder.R:153 | Keep | Andrew decided on 2026-10-19 that later verbs re-run it. | Low |
| `.lock_hash` | lock.R:17 | no | Hash of plan, data, columns, seed, K, version and dictionary | lock.R:140, :156, :193; tests 4 | Keep, make the dictionary an unconditional element in the same position | The conditional (:15-16, :21) serves 0.4.x locks. Hashes of 0.5.0 locks stay the same if the element order stays. | Protected: lock hash |
| `create_analysis_lock` | lock.R:74 | yes | Validates, refuses outcome columns, stores the dictionary, hashes | tests 42; vignettes 2; scripts 1; simulation 1 | Keep | | Protected: lock hash, blinding |
| `verify_lock` | lock.R:153 | yes | Recomputes the hash and warns on a version change | 8 sites in R (every verb); tests 12 | Keep | | Protected: lock hash |
| `save_lock` | lock.R:173 | yes | Verifies, then calls `saveRDS` | tests 2; scripts 1 (22_design_stage.R:233) | Cut, replace with `saveRDS()` | `load_lock()` and every verb verify the hash, so a changed lock is caught on read. | Low: the script needs one edit, after the run |
| `load_lock` | lock.R:181 | yes | Reads a lock and verifies it | tests 2; scripts 1 (25_estimate.R:31) | Keep | The verified reader. | Protected: lock hash |
| `.stamp` | lock.R:187 | no | Ties an object to the lock hash and its own hash | assess_design.R:143; nc_ladder.R:261; simulate_design.R:394; lock.R:198; tests 14 | Keep | | Protected: lock hash |
| `.check_stamp` | lock.R:192 | no | Checks the lock hash, the stamp and that the object is unmodified | dossier.R:97, :98, :110; nc_ladder.R:158, :159; simulate_design.R:286; tests 7 | Keep | lock.R:156-158 and :193-195 repeat the comparison and the message. They could share one helper, but both sites are protected, so I would leave them. | Protected: lock hash |
| `.short` | methods.R:1 | no | First 12 characters of a hash | methods.R:21, :65 | Keep | | None |
| `print.cr_plan` | methods.R:4 | S3 method | Prints the plan | tests 1 | Keep | | None |
| `print.cr_lock` | methods.R:19 | S3 method | Prints the lock | no vignette, script or test | Keep | Default display of the object. | None |
| `print.cr_design` | methods.R:28 | S3 method | Prints overlap, weights and edge check | vignettes 1 (cleanTMLE.Rmd:155); tests refuse an old object only | Keep, remove the version guard (:29) | The guard serves 0.4.0 objects. | None |
| `print.cr_simulation` | methods.R:39 | S3 method | Prints the verdict and the decision note | no vignette, script or test | Keep | Default display of the object. | None |
| `print.cr_nc_ladder` | methods.R:54 | S3 method | Prints the verdict by estimand | no vignette, script or test | Keep | Default display of the object. | None |
| `print.cr_dossier` | methods.R:63 | S3 method | Prints the decision and short hashes | vignettes 1 (cleanTMLE.Rmd:276) | Keep | | None |
| `print.cr_unblinded` | methods.R:77 | S3 method | Prints who approved and any override | no vignette, script or test | Keep | Default display of the object. | None |
| `print.cr_estimate` | methods.R:86 | S3 method | Prints one estimate | vignettes 1 (cleanTMLE.Rmd:335); tests 1 | Keep, drop the IMPLAUSIBLE line (:91) | Goes with the guard. | None |
| `.evalue` | methods.R:95 | no | E-value of a risk difference, through the risk ratio | methods.R:110; tests 0 | Cut | Only `summary()` calls it. EValue in Suggests brings metafor, dplyr and others for 12 lines. | Low |
| `summary.cr_estimate` | methods.R:109 | S3 method | Estimate plus an E-value | tests 3 | Cut | Not used by the case study, the vignettes or the manuscript. | Low |
| `print.summary.cr_estimate` | methods.R:114 | S3 method | Prints the summary | tests 1 | Cut | Goes with `summary()`. | Low |
| `as.data.frame.cr_estimate` | methods.R:126 | S3 method | One-row data frame of an estimate | scripts 1 (25_estimate.R:39); tests 1 | Keep | The case study uses it. | None |
| `tidy.cr_estimate` | methods.R:134 | S3 method | The same row with broom names | no caller anywhere | Cut | Duplicates `as.data.frame()`, and the only reason for `generics` in Suggests. | None |
| `as.data.frame.cr_simulation` | methods.R:141 | S3 method | Returns `x$metrics` | no caller anywhere | Cut | Users read `$metrics`. | None |
| `as.data.frame.cr_nc_ladder` | methods.R:144 | S3 method | Returns `x$table` | no caller anywhere | Cut | Users read `$table`. | None |
| `plot.cr_design` | methods.R:147 | S3 method | Histogram of the propensity score by arm | vignettes 1 (cleanTMLE.Rmd:176); tests 1 | Keep | The vignette shows it. | None |
| `plot.cr_simulation` | methods.R:158 | S3 method | Coverage by candidate and surface | no caller, no test | Cut | Dead code. The dossier tables carry the same numbers. | None |
| `.nc_grade` | nc_ladder.R:1 | no | Reads each control, grades domains and rungs | nc_ladder.R:250; tests 23 | Keep, drop the point and ci rules and the `in_band` column | The paper describes only the interval rule (methods.qmd:102-107), the case study uses it (22_design_stage.R:138), and the NEWS notes that a point-estimate band narrower than the sampling error fails a control with no effect most of the time (NEWS.md:106-109). `in_band` is "kept for compatibility" with no earlier release. | Touches the grading that gates unblinding; the interval logic is unchanged |
| `.nc_candidates` | nc_ladder.R:53 | no | Picks the candidate that each estimand's controls use | nc_ladder.R:184; tests 2 | Keep | | Protected: shared fitter |
| `negative_control_ladder` | nc_ladder.R:150 | yes | Estimates each control for every estimand on the full cohort and each restriction | tests 30; vignettes 1; scripts 1 | Keep | The case study passes no restrictions (22_design_stage.R:129-141); see the decision on restrictions below. | Protected: shared fitter |
| `.expand_library` | plan.R:5 | no | Expands "default" and checks learner names | plan.R:220, :225 | Keep | | Low |
| `.has_text` | plan.R:19 | no | True for one non-blank string | 11 sites (plan.R, dictionary.R) | Keep, as the only copy | `.is_text` (unblind.R:1) has the same body. | Low |
| `.norm_time_zero` | plan.R:26 | no | Validates and trims time zero | plan.R:217 | Keep | | Low |
| `.norm_eligibility` | plan.R:59 | no | Validates and trims the eligibility criteria | plan.R:218 | Keep | | Low |
| `.n_surfaces` | plan.R:84 | no | Counts the outcome surfaces | plan.R:257; tests 2 | Keep, drop its `q0` term and `%\|\|% 1` default | Both serve options or plans that no longer exist. | Low |
| `analysis_plan` | plan.R:168 | yes | Validates and stores every plan field | tests 19; vignettes 5; scripts 1; simulation 1 | Keep, drop `select`, the time-to-event arguments, `surfaces$q0` and the point and ci rules | The case study (22_design_stage.R:129-141) uses none of these options. It passes `rule = "interval"`, so keep that value accepted and refuse the others. The plan is hashed, so removing fields changes new hashes only. | Protected: lock hash |
| `.outcome_columns` | plan.R:316 | no | Names of the outcome columns | lock.R:85; unblind.R:87 | Merge (inline `plan$outcome`) | One column once the time-to-event form goes. | Low |
| `.candidates` | plan.R:318 | no | Crosses libraries with truncations | methods.R:7; nc_ladder.R:54; tests 7 | Keep | | Low |
| `.std` | simulate_design.R:1 | no | Standardizes a vector, with zero for a constant | simulate_design.R:13, :16, :18, :23, :25; tests 3 | Keep | | Low |
| `.surface_name` | simulate_design.R:6 | no | Name of an outcome surface | simulate_design.R:36 | Keep | | Low |
| `.surfaces` | simulate_design.R:10 | no | Builds the family of outcome surfaces | simulate_design.R:292; tests 14; simulation 1 | Keep, drop the `q0` branch (:17-19) | The `q0` path has no test and lets a steward's outcome-based predictions enter the design stage. | Protected: simulation (the paper simulation calls it, map.R:70) |
| `.true_values` | simulate_design.R:49 | no | True value of each estimand on a surface | simulate_design.R:294; tests 3; simulation 1 | Keep | | Protected: simulation |
| `.expected_events` | simulate_design.R:62 | no | Expected events among the rows each estimator uses | simulate_design.R:298; tests 1 | Keep | Shown in the dossier and used in the vignette (cleanTMLE.Rmd:310). Its role as a version marker (dossier.R:94, nc_ladder.R:155) goes. | Low |
| `.failed_rows` | simulate_design.R:73 | no | Rows that mark a failed repetition | simulate_design.R:111, :116; simulation 3 | Keep | The paper simulation calls it (key.R:16, map.R:53). | Low |
| `.one_rep` | simulate_design.R:85 | no | One simulation repetition | simulate_design.R:300, :311, :312; tests 9; simulation 1 | Keep, drop `surfaces_run` and `libraries_run` | Only tests pass them. `simulation/R/map.R:26-65` holds a forked copy of this function. | Protected: shared fitter |
| `.wilson` | simulate_design.R:129 | no | Wilson interval for a proportion | simulate_design.R:155 | Replace with `stats::prop.test(x, n, correct = FALSE)$conf.int` | I believe this is the same score interval (read from the documentation, not run). Keep the guard for zero repetitions. | Low |
| `.sim_metrics` | simulate_design.R:138 | no | Bias, coverage, RMSE, precision, pass and unresolved flags per cell | simulate_design.R:331; tests 15; simulation 1 | Keep, drop `coverage_mcse` and `rmse_mcse` | Computed at :158-159 and :162, read nowhere in the package or the template. | Protected: simulation |
| `.candidate_summary` | simulate_design.R:194 | no | One row per candidate of one estimand's metrics | nc_ladder.R:51, :66; simulate_design.R:223, :233; tests 3 | Keep | | Protected: simulation |
| `.decision_candidates` | simulate_design.R:210 | no | Candidates whose unresolved cells can change a decision | simulate_design.R:224, :245 | Keep, simplify when `select` goes | | Protected: simulation |
| `.decision_cells` | simulate_design.R:219 | no | Unresolved cells that can change a decision | simulate_design.R:320, :336; tests 3 | Keep | | Protected: simulation |
| `.sim_verdict` | simulate_design.R:230 | no | Feasible or infeasible, and the selected candidate | simulate_design.R:332; tests 6; simulation 1 | Keep | | Protected: simulation |
| `simulate_design` | simulate_design.R:282 | yes | Plasmode simulation, selection and targeted extension | tests 23; vignettes 2; scripts 2; simulation 2 | Keep | | Protected: shared fitter |
| `.is_text` | unblind.R:1 | no | True for one non-blank string | unblind.R:55, :85, :95 | Merge into `.has_text` | Same body as plan.R:19. | Low |
| `.blocked_reasons` | unblind.R:7 | no | Reasons a dossier blocks unblinding | unblind.R:54, :94; tests 4 | Keep | | Protected: blinding |
| `.unblind_hash` | unblind.R:13 | no | Hash of dossier hash, approval and outcomes | unblind.R:51, :114; tests 1 | Keep | | Protected: blinding |
| `.check_outcome_coding` | unblind.R:17 | no | Checks the outcome coding | unblind.R:93 | Keep the binary half, drop the time-to-event half (:26-37) | | Touches `unblind()` input checks, not the blinding |
| `.verify_unblinded` | unblind.R:43 | no | Rechecks lock, dossier, hash and gate | design_log.R:136; estimate_effect.R:143; unblind.R:6; tests 2 | Keep | | Protected: blinding |
| `unblind` | unblind.R:79 | yes | Admits the outcome after the dossier is approved | tests 39; vignettes 5; scripts 3 | Keep | | Protected: blinding |
| `%\|\|%` | utils.R:1 | no | Null default | 16 sites in R | Keep, replace with base `%\|\|%` when R 4.4 is the minimum | Base R has it from 4.4.0; the package needs 4.1. | None |
| `.hashable` | utils.R:3 | no | Turns plan content into something stable to hash, including closures | utils.R:16, :24, :32 | Keep | Needed for `custom` surfaces. It could drop the closure branch (:4-20) if `custom` also went. | Protected: lock hash |
| `.hash` | utils.R:32 | no | SHA-256 of a hashable object | dossier.R:54; lock.R:22, :189; unblind.R:14; tests 12 | Keep | | Protected: lock hash |
| `.bound` | utils.R:34 | no | Keeps probabilities inside (eps, 1 - eps) | 12 sites in R; tests 5 | Keep | | Low |
| `.with_seed` | utils.R:40 | no | Runs code under the default random generators and restores the caller's | learners.R:71; simulate_design.R:15, :92, :112; tests 1; simulation 3 | Keep | | Protected: reproducibility |
| `.make_folds` | utils.R:45 | no | Cross-fitting fold assignment | 8 sites in R; tests 17; simulation 2 | Keep | | Protected: cross-fitting |
| `.need_design_fields` | utils.R:53 | no | Refuses a design made by 0.4.0 | dossier.R:93; methods.R:29, :148; simulate_design.R:285; tests 0 direct | Cut | Unreachable in practice. A 0.4.x lock is refused by `.need_dictionary` first, and a 0.4.x design carries the hash of a 0.4.x lock. | Low |
| `.utc_now` | utils.R:60 | no | UTC time stamp | dossier.R:167; lock.R:143; unblind.R:107 | Keep | | None |

## Tests

A test block that only checks a class, restates the implementation, or repeats another file is marked cut or merge. The folder `_snaps/` is empty and no test uses a snapshot.

| Test file | What it protects | Recommendation | Reason |
|---|---|---|---|
| helper-fixtures.R:1-133 | Shared fixtures, and `contains_outcome()`, the blinding probe | Keep | `vignette_plan()` (:60-70) copies the vignette's plan by hand, so the two can drift. |
| test-assess-design.R:1-32 | Overlap grade, weights, balance, and that the simulated world's score depends on covariates only | Keep | The second test protects the design stage's central correction (0.4.0 NEWS). |
| test-assess-design.R:34-41 | Stamp refuses a modified lock | Merge into test-lock.R | `.check_stamp` is tested again at test-lock.R:50-59 and :73-79. |
| test-assess-design.R:43-57 | A balance failure warns and returns NA | Keep | Protects against a silent pass. It mocks cobalt only to trigger our own handler. |
| test-assess-design.R:59-142 | Population tables and small-cell suppression | Keep, merge the two suppression tests | :78-98 and :100-142 assert the same suppression rule twice. |
| test-assess-design.R:144-157 | Population tables with and without a trim band | Keep | |
| test-assess-design.R:159-163 | That the dossier carries the population tables | Cut | Checks `is.data.frame` and a list type only. The dossier tests build the dossier anyway. |
| test-assess-design.R:165-183 | C-statistic, largest and 99th-percentile weights | Keep | The first two assertions in :172-183 recompute from `.estimand_weights`, which is circular. |
| test-assess-design.R:185-191 | `plot()` on a design | Fold into one smoke test | Checks the class and two factor levels only. |
| test-blinding.R:1-6 | No design-stage object holds the outcome | Keep and strengthen | It finds only a field named `y` or an exact copy of the vector. The 0.3.0 test compared every verb's output under a permuted outcome (NEWS.md:187-191), which is stronger. |
| test-blinding.R:8-12 | Simulation runs without an outcome column | Cut | A class check plus a column check already made at test-lock.R:1-6. |
| test-design-log.R:1-20 | The JSON log carries both hashes, the decision, approval and estimates | Keep, merge the two tests | Protects the record. Stays only if the JSON format stays. |
| test-design-log.R:22-28 | `cr_example` shape | Move to test-example-data.R | Misfiled. |
| test-design-log.R:30-37 | Exactly the exported verbs | Keep | Protects the surface; update the list after any export is cut. |
| test-design-log.R:39-78 | Muntner columns, one row per step, roles and override | Keep | The second test calls `skip_if_not(... "STOP")` (:70), so it can skip without notice. Make the fixture produce a STOP. |
| test-design-log.R:80-104, :106-119 | Negative-control row names the estimand; by-rung verdicts in the JSON | Keep | |
| test-design-log.R:121-137, :152-169 | Wording of the grading rule, twice | Cut to one test | String matching on prose, tested at two places. Almost nothing is left once only the interval rule remains. |
| test-design-log.R:139-150, :189-195 | Dossiers built before 0.5.0 | Cut | The code they test is unreachable. |
| test-design-log.R:171-187 | Each identification flag is a log row | Keep | |
| test-dossier-unblind.R:1-30, :32-118 | Ladder walk, refusal of foreign objects, id alignment, STOP gate, gate follows the primary, no feasible estimand | Keep | Core contract of the dossier and `unblind()`. |
| test-dossier-unblind.R:120-129 | Ladder made before 0.4.1 | Cut | Duplicates test-nc-ladder.R:487-511 and tests unreachable code. |
| test-dossier-unblind.R:131-147 | Dossier needs `nc`, blank override and approver, no matching id | Keep | |
| test-dossier-unblind.R:149-233, :372-393, :448-457 | Six Quarto renders with about 35 phrase matches | Merge into one render of a rich fixture | Each render takes seconds and skips without Quarto. One render with an interval rule, an extension table, an exclusion and an identification flag covers the same phrases. |
| test-dossier-unblind.R:235-240 | `unblind_hash` equals `.hash(list(...))` | Cut | Restates the implementation. The tamper tests at :242-271 are the contract. |
| test-dossier-unblind.R:242-287 | Hand-built object past STOP is refused; edits to override or outcome are detected; binary coding | Keep | Protects the blinding. |
| test-dossier-unblind.R:289-315 | Time-to-event coding | Cut with the time-to-event path | |
| test-dossier-unblind.R:317-325 | Creation time and plan table | Keep | |
| test-dossier-unblind.R:327-334, :336-353 | Lock without a creation time; objects from 0.4.0 | Cut | Unreachable code. |
| test-dossier-unblind.R:355-361 | `.flat` formatting | Merge into the plan-table test | |
| test-dossier-unblind.R:363-370 | Times end in Z | Merge into test-lock.R:93-104 | The same pattern is tested there. |
| test-dossier-unblind.R:395-401 | A control that is a covariate is refused at the report | Keep | One of four tests of one check, see test-lock.R and test-nc-ladder.R. |
| test-dossier-unblind.R:403-446, :459-521 | Identification table content, and restriction and follow-up flags | Keep | Overlaps test-identification.R:6-89 and test-design-log.R:171-187; run them on one shared fixture. |
| test-estimate-effect.R:1-30 | The dossier's candidate, K and fold seed reach the shared fitter; a non-primary estimand is labelled | Keep | Protects the central contract. :23-30 can skip without notice. |
| test-estimate-effect.R:32-40 | Weighting equals a direct WeightIt call | Cut with the weighting comparator | Compares our wrapper with the call it makes. |
| test-estimate-effect.R:42-47 | Matching is ATT only | Cut with matching | |
| test-estimate-effect.R:49-78 | Factor and oddly named covariates; the message for an estimand with no candidate | Keep, drop the weighting calls | |
| test-estimate-effect.R:80-86, :134-142 | E-value and tidy | Cut | `is.numeric` and a one-row check only; re-tests EValue. |
| test-estimate-effect.R:88-111 | Time-to-event through concrete | Cut with the path | |
| test-estimate-effect.R:113-132, :144-150 | Weighting and TMLE target one trimmed population; comparator note | Cut with weighting | The TMLE side is tested in test-fit-candidate.R. |
| test-example-data.R:19-26 | On `cr_example` the ATE is infeasible and the trimmed ATE and ATO are feasible | Keep | Ties the vignette story to the data. It runs a full simulation. |
| test-example-data.R:28-47, :49-56 | Controls read on the primary estimand; the shipped dictionary locks | Keep | The first overlaps test-nc-ladder.R; the second absorbs test-design-log.R:22-28. |
| test-fit-candidate.R:1-37 | Design matrix, cross-fitted score for fold 1 from fold 2 only, truncation reaches ATE and ATT | Keep | Protects cross-fitting. Add the same check for the outcome fit. |
| test-fit-candidate.R:39-83, :91-106 | ATO against the truth and the augmented estimator; one row per estimand and truncation; failures become rows | Keep | |
| test-fit-candidate.R:85-89 | The implausibility flag | Cut with the guard | |
| test-fit-candidate.R:108-134 | The ATT step of tmle is skipped without changing the ATE | Keep :119-134, cut :108-117 | The first only checks an argument name. The second compares with `tmle::tmle` and protects the speed-up. |
| test-identification.R:6-89 | Integration test on a cohort with registry-style timing problems | Keep | |
| test-identification.R:91-111 | White-space trimming in time zero and eligibility | Keep, move to test-plan.R | |
| test-identification.R:113-130 | `.flatten_text`, two tests for one function | Merge into one | |
| test-identification.R:132-135 | The template text has no `cleanTMLE:::` | Cut | A string search over a file; a render test would fail if the template called the package. |
| test-learners.R:9-14 | The default library has 17 settings, one type per family | Keep | It will change with the learner-defaults item. |
| test-learners.R:16-51, :67-103 | xgboost round prediction, super learner contract, failing learner, single glm, simplex weights, grouped folds | Keep | Protects the learner and cross-fitting code. |
| test-learners.R:53-65 | Edge check, xgboost only | Keep, add glmnet and nnet cases | The nnet branch is the one behind the open learner-defaults item. |
| test-lock.R:1-48, :61-71, :81-91 | Outcome columns refused, hash breaks on any change, bad inputs, save and load, closure hashing, 0/1 controls | Keep | Protects the blinding and the lock hash. |
| test-lock.R:50-59, :73-79 | Stamps detect a foreign lock and a modified object | Merge into one | Same function, same setup, tested twice here and once in test-assess-design.R and test-nc-ladder.R:197-206. |
| test-lock.R:93-104 | Roles hashed, creation time not hashed, time format | Keep | |
| test-lock.R:106-112, :140-146 | A control that is a covariate is refused at the lock | Merge into one | Same check as test-dossier-unblind.R:395-401 and test-nc-ladder.R:615-623. |
| test-lock.R:114-138 | `nc_exclude` values, and its place in the hash | Keep | |
| test-lock.R:172-255, :276-284, :350-385 | Dictionary refusals, responses, treatment timing, one-variable agreement | Keep | The time-zero and dictionary refusals are the contract of 0.5.0. |
| test-lock.R:257-274 | Nodes after treatment are kept but never adjusted for | Keep | Option to cut with the role, see below. |
| test-lock.R:286-348 | Restriction variables | Keep | Option to cut with restrictions, see below. |
| test-lock.R:387-420 | Dictionary type and order do not change the hash; stored form; hashed | Keep, merge :401-412 into :387-399 | :401-412 checks structure only (class, names, row names). |
| test-lock.R:422-438 | A lock without a dictionary verifies but no verb accepts it | Cut | Tests the 0.4.x guard. |
| test-nc-ladder.R:14-93 | No controls, null and treatment-driven controls, every estimand per control, inestimable control, one-arm rung | Keep | |
| test-nc-ladder.R:54-68 | Verdict equals the full-cohort rung | Merge into :165-195 | Same assertion. |
| test-nc-ladder.R:95-163 | `.nc_grade` for the point and ci rules and precedence | Keep the precedence cases, cut the point and ci cases | About 40 of 69 lines are specific to the two rules. |
| test-nc-ladder.R:165-195, :208-274 | Gate is the full cohort's verdict; non-finite controls; interval rule | Keep | |
| test-nc-ladder.R:197-206 | A design stamped by another lock is rejected | Merge into test-lock.R | Stamp test again. |
| test-nc-ladder.R:276-303 | Point and ci mapping; interval versus point | Cut with the rules | |
| test-nc-ladder.R:305-485 | Each estimand's controls use the simulation's candidate; whole candidate reaches the fitter; refusal of foreign simulations | Keep | Protects the shared-fitter contract in the ladder. |
| test-nc-ladder.R:487-521 | Controls run on another simulation; ladders and simulations from earlier versions | Keep :487-493, cut the rest | The rest tests unreachable code. |
| test-nc-ladder.R:523-613 | `nc_exclude` removes the right columns, shares one propensity fit, runs on the real fitter | Keep, merge the three mock tests | |
| test-nc-ladder.R:615-623 | A control that is a covariate is refused at the ladder | Merge | Third copy of the check. |
| test-plan.R:1-87 | Defaults, refusals, whole numbers, `max_reps`, surface vectors, half-width tolerance, candidates | Keep | The rule test at :80-87 shrinks when only the interval rule stays. |
| test-plan.R:89-110, :121-147, :153-162 | `nc_exclude` checks, time zero and eligibility, eligibility in the plan table | Keep | |
| test-plan.R:112-119, :148-152 | The printed `nc_exclude` line; plan table names exist | Cut | A regular expression over printed text, and a name check that :153-162 covers. |
| test-simulate-design.R:3-142 | Surfaces and true values, metrics, unresolved cells, verdict | Keep | Core of the simulation contract. |
| test-simulate-design.R:144-217, :219-257 | Constant outcome, discriminating power, draws independent of the pair, RNG kind, extension to `max_reps` | Keep | :174-196 changes when `surfaces_run` goes. |
| test-simulate-design.R:259-381 | Extension reaches only the cells that can change the decision | Keep | It mocks `.one_rep` with a copy of its signature, so any signature change breaks it. |
| test-simulate-design.R:383-396 | The plan's K reaches the folds | Rewrite | It asserts every call of `.make_folds` uses K, which holds only for a glm library because inner folds use V. |
| test-simulate-design.R:398-413 | The world uses the single-fit score | Keep | Protects the design stage's central correction. |
| test-simulate-design.R:415-421 | Bootstrap copies never straddle folds | Cut | Indexes `.make_folds()` output itself, so it holds by construction. |
| test-simulate-design.R:423-467, :483-573 | NaN intervals, default surfaces pinned, duplicate names, precision columns, expected events, half-width failures | Keep | The default-surface pin restates the formula by design. |
| test-simulate-design.R:469-481 | A plan stored without `prognostic` | Cut | Tests a 0.4.0 plan. |
| test-utils.R:1-8, :22-31 | `.hash` ignores formula environments, sees function bodies and captured values | Merge into one | Two tests of one function. |
| test-utils.R:10-16 | Folds are balanced and seeded | Keep | |
| test-utils.R:18-20 | `.bound` | Cut | Restates `pmin(pmax(...))`. |
| `_snaps/` | Nothing | Delete the directory | Empty, and no test uses a snapshot. |

## Dead code and duplication

Code with no caller, or reachable only from tests, is `plot.cr_simulation`, `tidy.cr_estimate`, the two `as.data.frame` methods for simulations and ladders, `.evalue` with `summary()`, and `.estimate_matching`. The outputs `coverage_mcse` and `rmse_mcse` (`simulate_design.R:158-162`) are computed and read nowhere. The `in_band` column is "kept for compatibility" (`nc_ladder.R:7-9, :24-26`) with no earlier release to be compatible with. The arguments `surfaces_run` and `libraries_run` of `.one_rep` are passed only by tests. The fields `outcome_type` (`dossier.R:151`) and `component_hashes` (`dossier.R:165`) of the design report are written and never read; the second is an audit record and should stay.

Several guards cannot fire on a 0.5.0 object. Two are `.need_design_fields` and `.need_dictionary`, because a 0.4.x lock is refused first and a 0.4.x design carries a 0.4.x lock hash. The same holds for the version checks in `design_report` (`dossier.R:94-96, :102-110`), `negative_control_ladder` (`nc_ladder.R:155-157`), `.identification_rows` (`design_log.R:30-33`), `.muntner_log` (`design_log.R:74-79`), the `%||%` defaults for `prognostic` and `max_reps` (`plan.R:86`, `simulate_design.R:28, :325`) and the conditional dictionary in `.lock_hash` (`lock.R:15-16, :21`).

Two reachable options are unused. The `surfaces$q0` branch (`simulate_design.R:17-19`, `lock.R:127-128`, `dossier.R:7-8`) is untested, and it lets predictions that a steward fitted to the real outcome enter the design stage (`plan-and-dossier.Rmd:112-115`), which the blinding ledger would have to record. The `after_treatment` role is tested, but nothing reads those columns; the comment at `lock.R:132-133` ties them to the analysis of how the outcome is observed, which is on the wishlist.

Duplicated logic is easy to find. `.has_text` (`plan.R:19`) and `.is_text` (`unblind.R:1`) have the same body. The vocabulary before, at, after is defined twice (`plan.R:54` and `dictionary.R:10`). The three sentences that state the negative-control rule are written in `.nc_rule_text` (`design_log.R:8-26`) and again in the template (`dossier.qmd:154-168`), and a change to one will not reach the other. The lock-hash comparison and its message appear in `verify_lock` (`lock.R:156-158`) and `.check_stamp` (`lock.R:193-195`). Five stray-role blocks in `.check_dictionary` (`dictionary.R:142-174`) have one shape. `.estimate_weighting` and `.estimate_matching` copy the same extraction of the estimate, interval and risks (`estimate_effect.R:40-47` and `:68-75`). The mean of the kept population is computed in `.populations` (`assess_design.R:41-48`) and `.trimmed_profile` (`assess_design.R:59-77`). The trim-band test `g >= lo & g <= hi` is written seven times (`assess_design.R:4, :42`; `simulate_design.R:54, :66`; `fit_candidate.R:121`; `estimate_effect.R:36, :97`). Three helpers compute a hash over a stored object with slightly different recipes (`.stamp`, `.dossier_hash`, `.unblind_hash`); I would leave them alone because they are protected. Outside the package, `simulation/R/map.R:26-65` is a forked copy of `.one_rep`, so a change to the package function will not reach the paper simulation, and `simulation/R/summarize.R:36-38` recomputes metrics that `.sim_metrics` also computes.

Tests repeat each other in four places. The stamp check is tested in `test-assess-design.R:34-41`, `test-lock.R:50-59, :73-79` and `test-nc-ladder.R:197-206`. The check that a control is not a covariate is tested in `test-lock.R:106-112, :140-146`, `test-nc-ladder.R:615-623` and `test-dossier-unblind.R:395-401`. The wording of the grading rule is matched in `test-design-log.R:121-137, :152-169`. The same refusal of an old ladder is tested in `test-dossier-unblind.R:120-129` and `test-nc-ladder.R:487-511`. Four tests re-test another package or restate the implementation: `test-estimate-effect.R:32-40` (WeightIt), `test-estimate-effect.R:80-86` (EValue), `test-fit-candidate.R:119-134` (tmle, though it protects our speed-up) and `test-dossier-unblind.R:235-240` (our own hash). Two tests check a class and little else: `test-blinding.R:8-12` and `test-assess-design.R:159-163`. Six of the eight print methods have no test that prints a valid object.

## Size if every recommendation were taken

Line counts are estimated from the line ranges cited in the tables, not from an edit, so read them as about. The clear cuts need no decision: unreachable guards, dead methods, duplicate code, and tests that repeat or check nothing. The full set adds the decisions marked for Andrew (time-to-event path, both comparators, point and ci rules, `select`, `surfaces$q0`, the implausibility flag, the E-value summary, `save_lock`, and jsonlite moved to Suggests) and the new tests for the protected code.

| Measure | causalRisk | 0.5.0 today | After the clear cuts | After every recommendation |
|---|---|---|---|---|
| Exported functions (user-facing entries) | 29 entries, 26 names | 12 | 12 | 11 |
| S3 methods | 6 (plot, hist and forest plot methods) in its user-facing list | 16 | 12 | 10 |
| Internal functions | 64 reference entries labelled internal | 102 | 96 | 88 |
| Reference pages | 93 topics | 13 | 13 | 13 |
| Datasets | 8 | 1 | 1 | 1 |
| Lines of R | not available | 3,068 | about 2,950 | about 2,670 |
| Lines of tests | not available | 3,457 | about 3,250 | about 3,030 |
| Lines of the dossier template | not applicable | 231 | about 215 | about 215 |
| Imports entries (packages, plus stats and utils) | not in the PDFs | 13 (11 packages) | 13 | 11 (9 packages) |
| Suggests entries | not in the PDFs | 12 | 11 | 7 |

The next table shows where the saving comes from, in the order of the recommendations at the end. A negative figure is a saving. The first row adds tests, so its test figure is positive.

| Recommendation | Lines of R | Lines of tests | Imports | Suggests | Functions removed |
|---|---|---|---|---|---|
| 1. Add tests for protected code | 0 | +60 | 0 | 0 | 0 |
| 2. Drop compatibility code for pre-0.5.0 objects | -53 | -112 | 0 | 0 | 2 internal |
| 3. Cut the time-to-event path | -95 | -56 | 0 | -3 | 3 internal |
| 4. Cut dead methods, E-value summary, implausibility flag | -78 | -28 | 0 | -2 | 6 S3 methods, 2 internal |
| 5. Interval rule only | -22 | -80 | 0 | 0 | 0 |
| 6. Cut both comparators | -82 | -45 | -1 | -1 | 3 internal |
| 7. Cut `select` and `surfaces$q0` | -20 | -10 | 0 | 0 | 0 |
| 8. Consolidate tests | 0 | -158 | 0 | 0 | 0 |
| 9. One source for repeated text and code | -45 | 0 | 0 | 0 | 4 internal |
| 10. Move jsonlite to Suggests and keep the JSON log | +3 | 0 | -1 | +1 | 0 |
| 11. Cut `save_lock` (restrictions and the after-treatment role are not counted) | -8 | 0 | 0 | 0 | 1 export |

The clear cuts are rows 1, 2, 8 and 9, plus the dead methods in row 4. The Imports figure for row 6 is WeightIt. Installed size changes little, because `cobalt` already brings `ggplot2`, `chk` and `rlang`, and `tmle` depends on `glmnet` and `SuperLearner`. The saving is in code and tests, not in what a user installs.

Two further options are not in these totals because they change what the paper describes. Both are decisions for Andrew, and neither is a hidden cost of the cuts above.

| Option | Lines of R | Lines of tests | Why it waits |
|---|---|---|---|
| Remove `restrictions` from the plan, the ladder and the dictionary | -55 | -150 | `methods.qmd:102` and the vignettes describe restriction rungs, and `negative_control_ladder` would become a misnomer. |
| Remove the `after_treatment` role | -8 | -25 | Low value, but it changes the identification table the paper shows. |

## Top recommendations in priority order

1. Add the missing tests for protected code before any cut and before the learner-defaults change: a direct test of `.fit_Q` (`fit_candidate.R:38-46` has none), glmnet and nnet cases for the edge check (`test-learners.R:53-65` covers xgboost only), and a permuted-outcome test of the design verbs, which the 0.3.0 suite had (`NEWS.md:187-191`). It adds little code and protects what the learner-defaults item is about to change. Timing: 0.5.x release.

2. Remove the compatibility code for objects from before 0.5.0, with its tests. The merge commit calls 0.5.0 the first public version, and these guards cannot fire on a 0.5.0 object. It touches `.lock_hash`, so keep the element list and order, and every existing lock hash stays valid. Timing: 0.5.x release.

3. Cut the time-to-event path (decision for Andrew). It is the only route where the estimator run after unblinding is not the one the simulation evaluated (`estimate_effect.R:168`), nothing in the case study, the vignettes or the manuscript uses it, and the spec moved parity features to the wishlist (`2026-10-18-outcome-blind-design-stage.md:112`, `WISHLIST.md:54-55`). It is the largest single saving in R and removes three packages from Suggests. It touches the outcome checks in `unblind()` but not the blinding. Timing: 0.5.x release.

4. Cut the dead methods and the unused add-ons: `plot.cr_simulation`, `tidy`, and the two `as.data.frame` accessors have no caller at all. The E-value `summary()` and `.implausibility_check` are decisions for Andrew, because a sign change from the crude difference is what adjustment can produce and no decision reads the flag. EValue and generics leave Suggests. Timing: 0.5.x release.

5. Keep only the interval rule for negative controls (decision for Andrew), so that the code matches the paper (`methods.qmd:102-107`) and the case study (`22_design_stage.R:138`), which passes `rule = "interval"` and must keep running unchanged. About ten test fixtures use the point rule to get fixed verdicts and need rewriting. It touches the grading that gates `unblind()`, although the interval logic does not change. Timing: for the paper.

6. Decide on the comparators. Cut matching now (ATT only, unused). Weighting needs a decision because `methods.qmd:120` and `cleanTMLE.Rmd:349-353` describe it and the case-study scripts do not run it. My recommendation is to cut it, run WeightIt directly in the case-study script if the paper keeps a comparison, and update the Description field of `DESCRIPTION`, which names both comparators. It touches `estimate_effect()` but not the fitter it calls. Timing: for the paper.

7. Remove `select` and `surfaces$q0` from the plan (decision for Andrew) before the first release, while removal is cheap. An analyst-written selection function is a degree of freedom in a package built to remove them, and `q0` is untested and admits outcome-based predictions into the design stage. The vignette sentence at `plan-and-dossier.Rmd:112-115` goes with it. It touches the simulation's selection rule and surface builder but not the shared fitter. Timing: 0.5.x release.

8. Consolidate the tests: one Quarto render instead of six, one test per repeated check, no class-only or circular tests, and delete the empty `_snaps/`. This is the largest saving in tests after the compatibility tests, and it shortens the run time of the render tests. Timing: 0.5.x release.

9. Give repeated text and code one source: the grading-rule text in the template, `.is_text`, the timing vocabulary, a loop for the stray-role checks, base `deparse1`, `prop.test` for the Wilson interval, the unused arguments of `.one_rep`, and the unused outputs `coverage_mcse`, `rmse_mcse` and `in_band`. After the case study, point `simulation/R/map.R` at the package function instead of its forked copy. Timing: 0.5.x release, except the fork, which waits until after the paper.

10. Move jsonlite to Suggests, or cut the JavaScript Object Notation (JSON) log (decision for Andrew). The case study uses the Muntner format only, but `README.md` leads with JSON and the JSON log is the only one that carries the full design-report hash, so I would move jsonlite to Suggests and keep the format for now. Timing: after the paper.

11. After the case study and the paper, decide on `restrictions`, the `after_treatment` role and `save_lock`. The case study passes no restrictions and calls `save_lock` once (`22_design_stage.R:233`). Rename `negative_control_ladder` only if restrictions go. Timing: after the paper.

12. Housekeeping. Delete `attic/` (a rendered 0.3.0 vignette with its cache, a long functions document and an update script, all excluded from the build; the tag `v0.3.0` keeps the history). If any of these cuts land before the 0.5.0 tag, record the commit of the build that the case study ran on in the supplement. Timing: 0.5.x release.

## Writing check

I ran `check.py` on this file from the project root with `--min-level warn`. It reports errors of one kind and one warning. Every error is the rule that flags code formatting in prose: the backticked function names, file names with line numbers, and package names such as `causalRisk`. The instruction for this audit allows these, because every finding needs its evidence as `path:line` and the function names are the subject. No em dashes, banned words, number flags, undefined acronyms or other pattern flags remain. The one warning says that an overlap source named in `writing/config.yml` (a report in the earlier unblinded analysis folder) cannot be found; it does not concern this file, and I did not open that folder.

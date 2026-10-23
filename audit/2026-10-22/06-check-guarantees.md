# Check of the five priority recommendations against the code

Adversarial review, 22 October 2026, of the five highest-priority recommendations of the synthesis report, read against the cleanTMLE source, the Rescue.Co scripts, the simulation code and the installed dependency packages. I ran no R process. `audit/2026-10-22/REPORT.md` was not on disk when I checked the folder (it held `01` to `05`, `08-literature.md` and `lit/`), so I reviewed the five recommendations as they were given to me.

None of the five is refuted, and none survives unchanged. The common weak point is the guarantee that the estimator applied after unblinding is the one the design-stage simulation evaluated. In the code that guarantee rests on a version string: the lock and the design report record only cleanTMLE's version (`cleanTMLE/R/lock.R:139`, `cleanTMLE/R/dossier.R:167`), `verify_lock()` only warns when it differs (`lock.R:159-162`), nothing records which commit `rescueCo/lib` was built from, and the numbers also depend on tmle, glmnet, xgboost, earth and nnet, whose versions nothing records. tmle 2.1.1 itself changed its default targeting step (`C:/Users/andre/AppData/Local/R/win-library/4.4/tmle/NEWS:1-4`). Recommendation 1 needs that check widened. Its rehearsal of script 24 is not possible as written. Recommendation 2 would let the pilot drop the comparisons where claim 1 is weakest. The test bed that recommendation 3 names cannot exercise the defects it is meant to judge. In recommendation 4, the permuted-outcome test protects nothing in 0.5.0, and pins that call tmle a second time move with tmle. Recommendation 5 holds from the software side only if the hashes leave this machine before outcomes exist.

## 1. Freeze and rehearse every step after unblinding before script 24

### Verdict: survives with changes

### What I tried

I read `rescueCo/scripts/25_estimate.R`, `24_steward_outcomes.R` (with its uncommitted diff), `22_design_stage.R` and `23_registration_bundle.R`. On the package side I read `unblind.R`, `lock.R`, `estimate_effect.R` and `dossier.R`. I tried to show that the recommendation misreads the evidence, that its version check would catch what it needs to catch, and that its rehearsal could run without touching real outcomes.

The core of the recommendation holds. The loop has no handler (`25_estimate.R:29-41`), `unblind()` stops when no estimand is feasible or the negative-control verdict is `STOP` (`unblind.R:7-11, 94-98`), only the primary outcome is passed (`25_estimate.R:33`), no comparator is run (`25_estimate.R:36`), the longitudinal TMLE is a comment (`25_estimate.R:44-52`), script 22 builds no 180-day lock (`22_design_stage.R:130` fixes the outcome as `severe_disability`), and the uncommitted edit replaced the `APPROVAL.txt` guard with a check that nine design reports exist (`git diff -- rescueCo/scripts/24_steward_outcomes.R`; now `:17-21`). Seven parts need changing.

First, the refusal is stated as current fact ("Today a refused protocol comparison ... stops the loop"). The final run has finished only the comparison between operators (`rescueCo/logs/pipeline.log:1226`; `rescueCo/results/reexec_050/final/` holds only `L3_dossier.html`), so nobody yet knows whether the protocol comparison is refused. The risk is real but conditional: `list.files()` puts `L1_lock.rds` first, and any refusal, including a `STOP` verdict on a sensitivity lock, ends the run.

Second, the version check is too weak as proposed. The lock hashes the version string, not the code (`lock.R:17-23`), and a mismatch only warns (`lock.R:159-162`). `rescueCo/lib` was "built from the cleanTMLE submodule" (`22_design_stage.R:20-21`), the repository where the learner-defaults branch will be made. I found no record of the commit it was built from, and no `v0.5.0` tag exists (`git tag` in `cleanTMLE` lists only `v0.3.0`; the merge commit 2e96f98 says "not yet tagged or released"). The simulation library, by contrast, records its commit (`simulation/pilot/pilot-report.md:5`). A rebuild from a branch that kept the version string would pass `verify_lock()` silently. The fitted numbers also depend on packages that `rescueCo/lib` probably does not hold. I infer this from `22_design_stage.R:20-22`; under the task's rules I did not open `rescueCo/lib`. Those packages would come from the user library, where another project's update would change them. One example of such a change: tmle 2.1.1 changed its default targeting and the default `min.retain` used when it refits the treatment model (`tmle/NEWS:1-4`), and cleanTMLE passes neither argument (`cleanTMLE/R/fit_candidate.R:53-55`). The code itself treats xgboost's behaviour as version-dependent (`learners.R:50, 84`).

Third, script 24 cannot be rehearsed end to end on synthetic data. It reads the real registry and follow-up files and writes the real store (`24_steward_outcomes.R:24-29, 41-43`), and its arm-blind check compares the rescored indicator with the real one (`:33-35`). Script 25, too, writes its products beside the real locks (`25_estimate.R:35, 37`), so a rehearsal would leave synthetic unblinded objects and fits in `locks/final`.

Fourth, no exported function can produce "the other candidate". `estimate_effect()` takes the candidate from the design report (`estimate_effect.R:158`) and has no argument to choose another, and `fit_candidate()` is not exported (`cleanTMLE/NAMESPACE`). The script would have to call `cleanTMLE:::fit_candidate()` after `estimate_effect()` has verified the object. The recommendation also does not name which of the three other candidates is meant. If it means the super learner, that candidate failed the plan's coverage tolerance for the ATE: `final_operators_sim_super_learner_ate_coverage_min` is 0.88 against the tolerance of 0.90 at `22_design_stage.R:134`. It has to be labelled as an estimator the design stage did not accept.

Fifth, the longitudinal TMLE needs frozen inputs as well as frozen code. Every final lock holds complete cases only (`22_design_stage.R:164`), so the cohort the longitudinal analysis needs, with the patients lost at six months, sits in no lock and under no hash. lmtp's defaults differ from the package's conventions: it trims density ratios at their 0.999 quantile, uses `"SL.glm"` for both nuisances, and uses 10 inner folds (usage of `lmtp_control()` and `lmtp_tmle()` in `C:/Users/andre/AppData/Local/R/win-library/4.4/lmtp/help/lmtp.rdb`, read by decompressing the help database; `02-wrapper.md:28` reports the same from lmtp's source). Learners, folds, seed and trimming have to be written into the script before the store exists.

Sixth, `estimate_effect()` drops missing outcomes without a note (`estimate_effect.R:175-176`) and redraws the folds on the shorter vector (`:179`), so the analysis could use rows the simulation never described. For the primary outcome this should not happen, because "six-month status recorded" means the primary outcome is known (`rescueCo/BLINDING_LEDGER.md:42`). For a 180-day lock built on the same rows it may. `unblind()` already counts the gap (`unblind.R:111`), so the script can assert it is zero.

Seventh, approval needs a record outside the machine. Script 23 writes every lock and design-report hash (`23_registration_bundle.R:57-63, 76-85`), but `rescueCo/results/` is ignored by git (`.gitignore:10`) and this clone has no remote-tracking branch for `simulation-pilot`, which is 32 commits ahead of `origin/main`. Scripts 23 and 25 are already committed (55c5d81); the uncommitted ones are 22, `22_run_both.R` and 24.

### Evidence

| Point | Path and line |
|---|---|
| Final run has only the comparison between operators | `rescueCo/logs/pipeline.log:1226`; `rescueCo/results/reexec_050/final/` |
| Version string hashed, not code; mismatch warns | `cleanTMLE/R/lock.R:17-23, 139, 159-162`; `dossier.R:167` |
| No record of the commit behind `rescueCo/lib`; simulation records its commit | `22_design_stage.R:20-22`; `simulation/pilot/pilot-report.md:5`; `git tag` |
| tmle changed its default targeting in 2.1.1; cleanTMLE does not set it | `tmle/NEWS:1-4`; `fit_candidate.R:53-55` |
| Script 24 reads real files and writes the real store | `24_steward_outcomes.R:24-29, 33-35, 41-43` |
| No candidate argument; `fit_candidate()` internal | `estimate_effect.R:158`; `cleanTMLE/NAMESPACE` |
| Super learner below the coverage tolerance | `final_operators_sim_super_learner_ate_coverage_min`; `22_design_stage.R:134` |
| Locks hold complete cases only | `22_design_stage.R:164` |
| Missing outcomes dropped silently | `estimate_effect.R:175-179`; `unblind.R:111` |
| Hashes not tracked or pushed | `.gitignore:10`; `23_registration_bundle.R:76-85` |

### Proposed revision

> Before script 24 runs, make every step after unblinding fixed and checkable.
>
> - Rewrite script 25 so that it attempts every lock and records a refusal by `unblind()` as a result, with no override. As written, the script ends at its first refusal, and a protocol comparison sorts first.
> - Have script 25 stop unless `n_design_without_outcome` is zero for each lock.
> - Have script 25 estimate one named alternative candidate for the confirmatory comparison, written into the plan with the label "did not meet the design-stage tolerances" if it did not. This needs the internal `fit_candidate()` called on the verified object.
> - Pass the secondary outcome only if its lock exists, and build that lock on the rows whose secondary outcome is known.
> - Write the lmtp code in full, with its learners, folds, seed and `.trim` fixed, and save and hash its input cohort, including the patients lost at six months, before the store is built.
> - Have scripts 24 and 25 write to a directory given as an argument.
> - Once the design stage finishes, and before anything touches `rescueCo/lib`, record in the ledger the cleanTMLE commit, a hash of the installed namespace's functions, and the versions of R, tmle, SuperLearner, glmnet, xgboost, earth and nnet. Script 25 stops on any difference.
> - Replace the removed guard with a machine-readable approvals file (lock, design-report hash, approver, date) that scripts 24 and 25 check against each design report.
> - Commit the edits to scripts 22, 22_run_both and 24, the rewritten 25, the approvals file and a tracked copy of the registration bundle, and decide whether to push them before script 24 runs.
> - Decide whether the comparators at `methods.qmd:120` are run.
> - Rehearse script 25 on a synthetic outcome store drawn from the plasmode, written to a separate directory. Test the steward functions in `rescueCo/R/steward.R` on a synthetic follow-up file with the real file's columns. Never run script 24 itself before approval.

Decisions for Andrew: whether to push the branch (it is not on the remote today), and whether to tag the commit `rescueCo/lib` was built from.

## 2. Read the simulation pilot, then fix the claims, rules and scenarios

### Verdict: survives with changes

### What I tried

I checked whether the cached 0.4.1 results can stand in for 0.5.0, whether step 3a needs the package at all, and whether any part of the plan lets the pilot decide which claims survive.

The cached results can stand in. Between `release/0.4.1` and `main`, `R/learners.R`, `R/fit_candidate.R` and `R/utils.R` are identical, and `simulate_design.R` and `assess_design.R` gain only a `.need_dictionary()` call and one word (`git diff release/0.4.1 main -- R/` in `cleanTMLE`). Step 3a reads the cache through `collect()` (`simulation/R/run.R:87-94`) and the summary code calls no cleanTMLE function, so only the assertion at `simulation/pilot.R:15-16` ties it to 0.4.1. The registry evidence is right: every key in `writing/registry.json` begins `final_` or `first_four_months_` (my count), so none comes from the simulation study. `simulation/R/key.R:24-31` implements only the smallest-P rule, while `simulation/PLAN.md:119` names the narrowest CI.

Two things fail. The first is the instruction to "drop contrasts the budget cannot resolve". The sizing rule needs more repetitions the smaller the observed gap is relative to its spread (`simulation/R/pilot_check.R:7-10`). So dropping by budget removes exactly the comparisons in which the outcome-blind simulation's verdict and the best overlap diagnostic are closest, which is where claim 1 (`PLAN.md:9, 24`) would fail. That is the pilot setting the claim, which the recommendation says it wants to prevent, and it breaks the project rule that a contradicted claim is rescoped. `03-simulation.md:54` already lists dropping a truth after seeing labels as a move to rule out.

The second is the cache. A cached cell is reused when its key matches, and the key holds the cleanTMLE version string and the text of the simulation's own files, not cleanTMLE's code or tmle's version (`run.R:24-29`). Moving to 0.5.0 invalidates every cached cell in `simulation/cache/`, the answer-key fits and both map variants alike. A learner-defaults build that kept the version string would be reused silently.

### Evidence

`git diff release/0.4.1 main -- R/` in `cleanTMLE`; `simulation/pilot.R:15-16, 120-156`; `simulation/R/run.R:24-29, 87-94`; `simulation/R/pilot_check.R:7-10`; `simulation/R/key.R:24-31`; `simulation/PLAN.md:9, 24, 119`; `03-simulation.md:54`.

### Proposed revision

Keep the order the recommendation gives: rules first, then step 3a, then the full run. Replace its last bullet with:

> Use 3a only to size the full run. Keep every contrast in the planned tables. A contrast the budget cannot resolve is reported as unresolved, and the claim it serves is rescoped to match. Before the full run, add R's version and those of tmle, glmnet, xgboost, earth and nnet to the cache key at `run.R:24-29`, as `03-simulation.md` section 6.2 proposes.

## 3. Rules for the learner-defaults item, and the case study on 0.5.0 as installed

### Verdict: survives with changes

### What I tried

I checked whether staying on 0.5.0 keeps the shared fitter, whether the proposed test bed can show the defects, and whether keeping risk paths is cheap.

Staying on 0.5.0 keeps the shared fitter only if `rescueCo/lib` is never rebuilt from the repository where the item's branch will live (point two under recommendation 1). The guarantee listed as at stake is also overstated. A rerun of the design stage on a patched build before script 24 would read no outcome, so it would not break the outcome-blind claim. The reasons to stay on 0.5.0 are that the case study shows the release, and that choosing between two design reports would be a selection the paper must explain. Both are good reasons, and Andrew has made the call (`docs/superpowers/specs/2026-10-18-outcome-blind-design-stage.md:133`).

The test bed does not fit the defects. The recommendation judges every fix "on the simulation study's data-generating processes". Those use four covariates (`simulation/R/map.R:8`), a main-terms logistic propensity score and a main-terms outcome candidate only (`map.R:19-23`). The super-learner subset that `simulation/PLAN.md:101` and `manuscript/methods.qmd:158` describe is not in the code. In that setting the nnet decay grid and the elastic net's penalty on the treatment column never run. The in-sample world score is a logistic fit, not an in-sample fit of boosted trees with `min_child_weight = 0` (`learners.R:78`), so the mechanism the staged diagnostic is meant to test does not arise. I infer that boosted trees with that setting can fit their training rows closely; I did not measure it. The wrapper audit's own experiments use 4 and 75 covariates (`02-wrapper.md:187`), so the recommendation narrows the audit it cites.

Keeping "each learner setting's held-out risk path" is not free. The path for one fold of the default library has about 5,150 values: ten boosting settings of up to 500 rounds, three elastic-net paths of 50 penalties, three nets and one spline fit (`learners.R:6-12`; my arithmetic). Over five folds and four outcome surfaces that is about 0.8 megabytes of plain numbers per repetition, and the comparison between operators ran up to 300 repetitions on some cells (`final_operators_sim_datasets_ate`). A worker has already run out of memory on the larger comparisons (`rescueCo/scripts/22_run_both.R:5-7`). The plateau test needs each setting's risk and its standard error, not every round. The risks of the three decays are already computed per fit (`learners.R:201-204`) and discarded at `simulate_design.R:390`.

The report also leaves out the cheapest check in the wrapper audit. Comparing the c-statistic and the share outside the band of the world's score with those of the cross-fitted score needs no simulation (`02-wrapper.md:195`). The design object holds both scores (`assess_design.R:138`).

### Evidence

`docs/superpowers/specs/2026-10-18-outcome-blind-design-stage.md:133, 137-138, 145`; `simulation/R/map.R:8, 19-23`; `simulation/PLAN.md:101`; `manuscript/methods.qmd:158`; `cleanTMLE/R/learners.R:6-12, 78, 146, 201-204`; `simulate_design.R:16, 119, 238, 390`; `assess_design.R:103, 138`; `rescueCo/scripts/22_run_both.R:5-7`; `02-wrapper.md:187, 195`; registry `final_operators_sim_main_terms_ate_unresolved_cells`, `final_operators_sim_super_learner_ate_unresolved_cells`, `final_operators_sim_datasets_ate`.

### Proposed revision

> - Confirm that the re-execution stays on 0.5.0 as installed, and that `rescueCo/lib` is not rebuilt until script 25 has run.
> - First, compare the world's score with the cross-fitted score from the saved design objects.
> - Before any fix, keep each setting's cross-validated risk and its standard error, computed in place from the held-out losses. Keep path learners' risks at a few points, not every round, and only in the experiment code.
> - Run the staged known-truth diagnostic, with each decision rule written into the item's spec first.
> - Judge each fix on known-truth processes built without Rescue.Co data: the simulation study's processes with the default library added, and a synthetic process with about as many columns, binary and continuous, as the Rescue.Co design. Rescue.Co covariates serve only to check that a setting runs.
> - Hold the release only if the diagnostic finds a defect in the fitter or the targeting.

Decision for Andrew: whether a synthetic process at the Rescue.Co dimension is acceptable as the test bed. It is stricter than the spec, which forbids tuning to Rescue.Co results but not a plasmode on its covariates.

## 4. Test and pin the in-house numerical code

### Verdict: survives with changes

### What I tried

I looked for a reason to delegate the fitter instead, and asked whether each proposed test could catch the failure it is meant for.

Delegation would change results and cost time. One fit per fold gives the whole boosting and elastic-net path (`learners.R:123-147`), and bootstrap copies share folds through `groups` (`learners.R:117-121`; `simulate_design.R:100-104`). The wrapper audit's comparison with SuperLearner's wrappers shows different fold draws, leaf sizes, net settings and an inner `cv.glmnet` (`02-wrapper.md:28-36`). Keeping the fitter in-house survives.

Three of the proposed tests need changing. The ATE "pin" at `tests/testthat/test-fit-candidate.R:119-134` calls `tmle::tmle()` with the same arguments `.tmle_args()` builds, so it moves with tmle. It could not have caught the change of default targeting in tmle 2.1.1 (`tmle/NEWS:1-4`); the installed help shows `target.gwt = FALSE` in the usage (`tmle/help/tmle.rdb`, read by decompressing). The proposed ATT pin to "a direct `tmle::tmle()` call" has the same weakness. A pin needs either stored reference values or an independent implementation, such as a short logistic fluctuation with the clever covariate of the ATE. Requiring `tmle (>= 2.1.1)` is right, for a stronger reason than the `formals()` guard. With that floor, `.tmle_args()` can pass `target.gwt = FALSE` explicitly, so a later change of default cannot alter results silently.

The `.nnloglik()` pin should compare losses, not weights. `method.CC_nloglik` solves the same problem, a convex combination minimizing log-likelihood on the logit scale (its description in `SuperLearner/help/SuperLearner.rdb`). The default library has 17 members (`test-learners.R:9-11`), several of them near-duplicates, so the weights are not well determined. The tolerance of about 1e-3 on the weights in `02-wrapper.md:215` would fail on realistic inputs while the fits agree. This is an inference from the library's make-up; I did not run it.

The permuted-outcome invariance test protects nothing in 0.5.0. In 0.3.0 the lock held a masked outcome, so invariance had to be tested (`NEWS.md:189-191`). In 0.5.0 the lock refuses the plan's outcome columns (`lock.R:85-88`, tested at `test-lock.R:1-5`), no design-stage function takes outcome values, and `test-blinding.R:1-6` checks that no design-stage object holds the outcome. A permuted store has no path into the design stage, so the test would pass by construction.

The other proposals stand. No test calls `.fit_Q()` (`grep` of `tests/`). The edge check is tested only on boosted trees (`test-learners.R:53-65`). No test fits the net or the default library (fixtures use `ps_library = "glm"`, `helper-fixtures.R:19, 64`). The overlap-weighted TMLE is checked on one dataset with true nuisances (`test-fit-candidate.R:39-55`), not for coverage.

### Evidence

`cleanTMLE/R/fit_candidate.R:48-61, 86-114`; `learners.R:117-147, 159-180`; `tests/testthat/test-fit-candidate.R:39-55, 119-134`; `test-learners.R:9-11, 53-65`; `test-lock.R:1-5`; `test-blinding.R:1-12`; `helper-fixtures.R:19, 64`; `cleanTMLE/DESCRIPTION:28`; `tmle/NEWS:1-4`; `02-wrapper.md:28-36, 215`.

### Proposed revision

Replace the first two groups of bullets with:

> - Add a direct test of `.fit_Q()` (cross-fitting and the treatment column), glmnet and nnet cases for `.edge_check()`, and one test that fits the default library.
> - Pin `.nnloglik()` to `SuperLearner::method.CC_nloglik` on the minimized log-loss.
> - Pin the glmnet and xgboost paths to `cv.glmnet(foldid = )` and `xgb.cv(folds = )` on the same folds, comparing the selected path point.
> - Pin ATE and ATT targeting to stored reference values on a fixed dataset, and the ATE also to a short independent fluctuation, with tmle's propensity refit both triggered and not.
> - Check the coverage of the overlap-weighted TMLE with true nuisances.
> - Require `tmle (>= 2.1.1)`, pass `target.gwt` explicitly, and remove the `formals()` guard.

Drop the permuted-outcome test.

## 5. Restate the paper's claim as the joined workflow

### Verdict: survives with changes

### What I tried

My lens covers the parts of the claim the software must carry: "the estimator applied being the one evaluated" and "tamper-evident". The literature reading is outside it.

The first phrase holds only with the same code and the same dependency versions (recommendation 1, point two). It also needs one qualification the Methods already make. For the ATT, tmle can refit a logistic treatment model instead of using the cross-fitted score (`fit_candidate.R:49-52`; `manuscript/methods.qmd:112`). The simulation runs the same step, so the identity holds, but "cross-fitted nuisance models" is not true of every ATT fit.

"Tamper-evident" needs a copy of the hashes held where the analysts cannot change it, made before outcomes exist. In the case study the hashes now sit in a git-ignored folder on one machine (`.gitignore:10`). The governance decision moves them to the supplement (`2026-10-18-outcome-blind-design-stage.md:129`), which comes after unblinding. Unless they are committed and pushed before script 24 runs, a reader cannot tell whether a lock was made before or after the outcomes were scored.

The software also supports demoting the comparison of TMLE with weighting. The weighting comparator uses WeightIt with a main-terms logistic score and no truncation (`estimate_effect.R:39-40`), while the TMLE uses the plan's library with truncation. A difference between them mixes the learner, the truncation and the estimator.

### Evidence

`cleanTMLE/R/lock.R:139, 159-162`; `fit_candidate.R:49-55`; `estimate_effect.R:39-40`; `manuscript/methods.qmd:71, 112`; `.gitignore:10`; spec `2026-10-18-outcome-blind-design-stage.md:129`.

### Proposed revision

Add two conditions to the recommendation: the version record of recommendation 1 exists before the claim is written, and the hashes are pushed before script 24. Proposed text for the paper, to use only once both are true:

Before (`manuscript/methods.qmd:71`):

> The simulation fits every candidate with the same fitter, the same number of cross-fitting folds, the same library and truncation and the same targeting steps that the final analysis later uses, so that the candidate the simulation selects is the estimator applied to the outcome.

After:

> The simulation fits every candidate with the same fitter, the same number of cross-fitting folds, the same library and truncation and the same targeting steps that the final analysis later uses, under the same installed versions of cleanTMLE and of the packages that fit and target the models (tmle, glmnet, xgboost, earth and nnet), which we recorded before the outcomes were joined. The candidate the simulation selects is therefore the estimator applied to the outcome.

For the limits on the first page, a sentence for the Discussion:

> The hashes of the plan, the design data and the design report show whether any of them changed after it was fixed, but only against a copy held where the analysts cannot alter it; in the re-execution that copy was the project repository, pushed on [date of the push], before the outcomes were scored.

## Other problems

The synthesis report itself was not in `audit/2026-10-22/` when I looked, so I could check only the five recommendations as quoted.

Recommendation 1 says the old approval guard was removed. That is accurate, but the removed guard asked for a registration link (`git diff -- rescueCo/scripts/24_steward_outcomes.R`), and the governance decision of 20 October dropped registration. The edit followed a decision, so restoring "an approval check" is right and restoring the old one is not.

Recommendation 1 asks to "commit scripts 22 to 25". Scripts 23 and 25 are committed (55c5d81); the uncommitted files are 22, `22_run_both.R` and 24.

Recommendation 1 also says "every lock whose primary estimand is feasible". `unblind()` also refuses a feasible lock whose negative-control verdict is `STOP` (`unblind.R:9-10`), so the scope should be every lock `unblind()` admits.

The evidence for recommendation 4 calls `test-fit-candidate.R:119-134` an ATE pin to tmle. It checks that skipping tmle's ATT step leaves the ATE unchanged, against a second call to tmle with the same inputs. It pins nothing against a change in tmle.

I checked that a lock reproduces from its seed, given the same versions. Every repetition runs under its own seed (`simulate_design.R:91-92, 112`), folds and nets are seeded (`utils.R:293-301`, `learners.R:71-73`), and the boosted trees use no subsampling (`learners.R:77-79`). That holds, and it is one more reason the versions must be recorded.

## Writing check

Ran `python "C:/Users/andre/OneDrive/Documents/scientific-writing-skills/checks/check.py" audit/2026-10-22/06-check-guarantees.md --min-level warn` from the project root in Git Bash, several times. The first run flagged bold lead-ins, one banned word, two terms from the project glossary, a hash name read as a dash code, three undefined abbreviations and four numbers not in the registry. I fixed these in this file, dropping three of the four numbers.

The final run leaves two kinds of finding, both expected in this audit. Most remaining errors are flags for code formatting in prose; every one is a function name, file path, line reference or registry key in backticks, which this internal audit needs. One number is flagged, in recommendation 3: the count of values per fold, my arithmetic from the grid sizes at `cleanTMLE/R/learners.R:6-12`, so it is not a registry value. The one warning concerns the project's `writing/config.yml`, whose overlap source in `rescueCo_analysis` is not found; it is not about this file.

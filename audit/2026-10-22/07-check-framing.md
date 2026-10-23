# Adversarial check of the five top recommendations: science and framing

Written 22 October 2026 for Andrew Mertens. This is a review of proposals, and nothing in the package, the paper, the vignettes, the specs or any script was changed.

None of the five recommendations is refuted outright. Each protects something real, and most of the cited evidence says what the recommendation says it does. All five survive only with changes, and two of the changes come before anything else.

The first concerns recommendation 1. Script 24 cannot be rehearsed on synthetic data as it stands, because it opens the real follow-up file, and script 25 writes into the real lock folder and into the folder this audit treats as outcome-free. A rehearsal run as worded would either build the real outcome store early or put synthetic estimates next to the design-stage results.

The second concerns recommendation 5. Its headline, the joined workflow as the advance, is the claim the literature survey finds unsupported because each part has a published source (`audit/2026-10-22/08-literature.md:13`). Its proposed narrower gap points at the estimator, where outcome-blind simulation is already published [@nanceCausalRoadmapSimulations2024a], and not at the estimand, where the survey finds the gap.

Three recommendations also turn decisions that are Andrew's into directives without saying so. Recommendation 3's release rule narrows decision 10 and his instruction of 22 October (`docs/superpowers/specs/2026-10-18-outcome-blind-design-stage.md:125, 133, 141`). Recommendation 2 brings back the smallest-P comparator he replaced on 7 October (`simulation/PLAN.md:119`). Recommendation 5 demotes the comparison of TMLE with weighting, which he made a main point on 20 October.

`REPORT.md` was not on disk when I wrote this (I looked at 22:58 on 22 October), so I reviewed the five recommendations as quoted in my task and could not check the rest of the report.

## How I checked

I read both specs, audits 02, 03, 04 and 08 in full and the parts of audits 01 and 05 the recommendations cite, the blinding ledger, scripts 22 to 25 and the uncommitted diffs, `simulation/PLAN.md`, `pilot.R`, `key.R`, `map.R`, the pilot logs, the TMLE comparison proposal, the brief, the introduction and Methods drafts, `writing/references.bib`, the registry, and the cleanTMLE code each claim rests on. Git was used read-only. I started no R process and opened nothing under `rescueCo/data/`, `rescueCo/reports/`, `rescueCo/lib/` or directly in `rescueCo/results/`, nor `rescueCo/logs/pipeline.log`. From `rescueCo/deviations_inventory.csv` I printed only the id, departure, date, code, estimate-status and check-status columns of rows that mention scoring, which hold no outcome values.

## Recommendation 1: freeze and rehearse every step after unblinding

### Verdict

It survives with changes.

### What I tried

I tested whether the code says what the recommendation says. It does: script 25 loops over every lock file in alphabetical order with no handler, `unblind()` refuses a lock without an override, the protocol locks sort first, script 25 passes only the primary outcome, the longitudinal TMLE is a comment, the version check only warns, and script 22 builds no lock for death within 180 days. The approval guard in script 24 was replaced in the working tree by a check that nine design reports exist. These hold up.

The rule "estimate every lock whose primary estimand is feasible" misses one refusal. `unblind()` also refuses a lock whose negative-control verdict is `STOP` (`cleanTMLE/R/unblind.R:9-10`), so a lock with a feasible ATE and that verdict would stop the loop again. The rule has to be "every lock that `unblind()` admits without an override". The re-execution spec already says script 25 unblinds each lock (`docs/superpowers/specs/2026-10-20-rescueco-reexecution-0.5.0.md:91`), so this rule adds little new discretion, which is good, but the role of each estimate (one confirmatory, the rest descriptive or sensitivity) still has to be written into the plan. Audit 04 left which locks are estimated as Andrew's decision (`audit/2026-10-22/04-case-study.md:186`); the report states a rule in its place.

"Commit scripts 22 to 25" is inaccurate. All four are in commit 55c5d81. The working tree changes only `22_design_stage.R`, `22_run_both.R` and `24_steward_outcomes.R` (git status). The edit to script 22 changes the `reviewers` text of the plan (`rescueCo/scripts/22_design_stage.R:139-140`), and the plan is inside the lock hash (`cleanTMLE/R/lock.R:17-22`), so only the working-tree version reproduces the locks on disk. Audit 04's sentence that the script "is not itself in git" (`04-case-study.md:53`) is wrong in the same way.

"Stop if the installed cleanTMLE is not the build that made the locks" cannot be done from the lock. The lock stores only `packageVersion()` as a string (`cleanTMLE/R/lock.R:1, 139`), and a patched build that keeps the version 0.5.0 would pass. A fingerprint of the build has to be recorded now, before any learner-defaults patch exists.

The rehearsal, as worded, is unsafe. Script 24 reads the registry and the follow-up file directly (`rescueCo/scripts/24_steward_outcomes.R:24-29`) and has no synthetic input, so "rehearse 24 end to end" means building the real outcome store unless a path argument is added first. Script 25 writes `_unblinded.rds` and `_fit.rds` into the lock folder and the estimates and design logs into `rescueCo/results/reexec_050/` (`25_estimate.R:24-25, 35-42`), the folder this audit's rules treat as outcome-free and from which `writing/export_results.R:13` builds the registry inputs. A rehearsal must write to a scratch folder, and the real run should too.

The approval check is an addition to Andrew's governance decision, not a restoration. The removed guard required an approval file that named a registration link (git diff of script 24), which the governance decision dropped. That decision records approval through `unblind(approved_by = )` (`docs/superpowers/specs/2026-10-18-outcome-blind-design-stage.md:129`), which runs in script 25, after script 24 has scored the outcomes. The re-execution spec still says script 24 runs only after approval (`2026-10-20-rescueco-reexecution-0.5.0.md:13, 17`), so the check reconciles the two specs, but it is Andrew's call. The plan names him as analyst, steward and approver (`22_design_stage.R:139-140`), so the check is a dated record of order, not a separate control.

I looked for anything from the unblinded analysis that the freeze would carry in. Script 24 scores the outcome as the unblinded analysis does (the re-execution spec, line 91), and the approved inventory dates three parts of that definition to the unblinded analysis with estimates in view: the primary outcome (row 8, 5 May 2026, estimates cited), scoring of the Glasgow Outcome Scale-Extended and of impossible answers (row 18, 15 June 2026, cited), and the six-month horizon (row 30, 31 August 2026, available). This was approved on 20 October and is not a reason to reopen it now, but the freeze record and the paper have to say that the outcome definition was inherited, and audit 04's decision table (`04-case-study.md:99-108`) omits these rows. Every new choice the freeze makes, such as whether the other candidate is estimated, is also made by people who had seen the unblinded estimates (`rescueCo/BLINDING_LEDGER.md:7`), so the fewer new choices, the better.

The stake is overstated. "The case study's whole outcome-blind claim" is already limited to the software withholding the outcome (`writing/brief.md:32`; `manuscript/methods.qmd:200`). What the freeze protects is narrower: that the re-execution's plan and code did not change after its outcome store existed, and that the primary estimator applied is the one the simulation evaluated (`methods.qmd:71`).

### Evidence

| Claim in the recommendation | What the source says | Holds? |
|---|---|---|
| A refused protocol lock stops the loop first | `25_estimate.R:29-41`; `unblind.R:7-11, 94-98` | Yes |
| Rule: estimate every lock with a feasible primary estimand | `unblind.R:9-10` also refuses a `STOP` verdict | Needs "admitted without override" |
| Commit scripts 22 to 25 | 55c5d81 holds all four; three files differ in the working tree | Narrow to the three |
| Stop if the build differs | `lock.R:1, 139` stores the version string only | Needs a fingerprint |
| Rehearse 24 and 25 on synthetic outcomes | `24_steward_outcomes.R:24-29` reads real files; `25_estimate.R:24-25, 35-42` writes into real folders | Needs input and output paths |
| Restore an approval check | Spec 2026-10-18 line 129 records approval in `unblind()`; re-execution spec lines 13, 17 | Addition for Andrew |
| No lock for death within 180 days | `22_design_stage.R:130`; re-execution spec line 73 requires one; `methods.qmd:187, 197` promise it | Yes, and a spec departure |
| Comparators promised | `methods.qmd:120` | Yes |

### Proposed revision

> **Freeze and rehearse every step after unblinding, before script 24 builds the outcome store.** Write into the plan that script 25 unblinds every lock, estimates each lock that `unblind()` admits without an override with the candidate the simulation selected, reports each refusal (no feasible estimand, or a `STOP` negative-control verdict) as a result, and labels one estimate confirmatory and the others descriptive or sensitivity analyses. Andrew decides whether the candidate the simulation did not select is estimated for the confirmatory comparison; if it is, the plan gives the outcome-free reason (its cells never resolved) and says the selected candidate stays primary whatever the comparison shows.
>
> Decide the outcome of death within 180 days before script 24: build its lock, as the re-execution spec says (line 73), or remove it from the Methods, since without a lock its estimate would not come from an estimator the simulation evaluated. Decide whether the comparators at `methods.qmd:120` run, after the parsimony decision on the comparators. Write the longitudinal TMLE code into script 25.
>
> Record a fingerprint of the cleanTMLE build in `rescueCo/lib` (its source commit and a hash of the installed package folder) in the ledger now, and have script 25 stop if the installed build differs. Andrew decides whether script 24 checks for a ledger record of his approval of each design-report hash; this adds to his governance decision and is a dated record of order, since he holds every role. List in the freeze record each scoring rule script 24 applies, with its inventory row and whether estimates existed when it was set (rows 8, 18, 30).
>
> Have script 25 write its estimates and logs to a new folder, not `rescueCo/results/reexec_050/`. Commit the working-tree edits to `22_design_stage.R`, `22_run_both.R` and `24_steward_outcomes.R`. Rehearse script 25 end to end on a synthetic outcome store drawn from the plasmode, with every output in a scratch folder; rehearse script 24 only on a synthetic follow-up file with the same columns, through a path argument added for the purpose. What is at stake: that the re-execution's plan and code did not change after its outcome store existed, and that the primary estimator applied is the one evaluated.

## Recommendation 2: fix the simulation's claims and rules, then read the pilot

### Verdict

It survives with changes.

### What I tried

The evidence checks out. Every registry key I counted begins `final_` or `first_four_months_`, so no simulation result exists yet. `simulation/pilot.R:15-16` asserts version 0.4.1, and step 3a is at `:120-156`. The pilot log shows the label fits and both map variants cached (`simulation/pilot/step3_run_output.txt:1`; `step3_resume2_stdout.txt`), and `simulation/pilot/` holds no `step3.rds` or `step3_*.csv`. `pilot-report.md:3` still says step 3 has not started. `simulation/R/key.R:24-31` implements only the smallest-P rule.

The title and the body disagree. The title says to read step 3a and then fix the claims, which is audit 03's order (`03-simulation.md:270`). The body says to write the rules first and run step 3a afterward, and calls this stricter. The body's order is the one I support; the title has to change.

The recommendation reopens two approved decisions without saying so. On 7 October Andrew approved the narrowest CI as the claim 3 comparator, "a milder rule than the smallest P value" (`simulation/PLAN.md:114, 119`). Audit 03's reason for adding the smallest-P rule back (`03-simulation.md:38`) predates any claim 3 result, so it is not tuning, but it is his decision. The baseline-risk factor in audit 03 (`03-simulation.md:208`) varies the truth's baseline risk on three overlap levels, which replaces the approved family-baseline variants at 0.10 and 0.30 (`PLAN.md:116`). The nonlinear-treatment truth comes from the TMLE proposal (`simulation/PROPOSAL-tmle-comparison-2026-10-20.md:56`), whose three decisions are still open (`:66-70`; `writing/outline-results.md:13`).

"Use 3a only to size the full run and to drop contrasts the budget cannot resolve" fails in two ways. Step 3a prints the agreement table and every verdict-against-diagnostic contrast (`pilot.R:147-151`), so whoever runs it sees their direction; "only to size" cannot be enforced, and the safeguard is fixing the rules first. More seriously, a small pilot gap on the main contrast is the case in which claim 1 fails, and the plan already says that then "claim 1 fails and the table shows it" (`PLAN.md:51`). Dropping that contrast would hide a null result. It has to stay in the table as unresolved, or be judged against the equivalence margin fixed beforehand, with the claim rescoped. Sizing then needs the run length for the equivalence test, not only the gap-based size the pilot check computes (`03-simulation.md:99`).

The baseline-risk factor makes part of claim 1 true by construction. Audit 03 adds it because "the design gives the verdict no condition in which a propensity-only summary must fail" (`03-simulation.md:26`). Choosing a factor that the comparators cannot see is acceptable when it is fixed before results, but against propensity-only diagnostics on that factor the verdict wins by design. A fair comparator is an outcome-free events rule an analyst could apply, such as expected events per covariate column [@peduzziSimulationStudyNumber1996], which the package can already supply before the outcome (`cleanTMLE/R/simulate_design.R:60, 389`). Table 1 should report that factor separately, and the claim should say what the simulation adds (overlap and expected events together), not that it beats overlap diagnostics in general. Inside the family the agreement is also partly definitional, since the inside truth is built from the family and the label uses the verdict's own tolerances (`03-simulation.md:52, 70`).

"Fix the rules before the pilot can shape the claims" is already partly late. Pilot steps 1 and 2, the smoke runs and the 40-dataset label check have been read (`simulation/pilot/key_check_output.txt`; `pilot-report.md:91, 99-101`), and audit 03 uses them to argue for rescoping the TMLE claim (`03-simulation.md:158`). That is legitimate, but the simulation plan should record beside each rule which pilot outputs had been seen when it was written, as it already does for the added overlap levels (`PLAN.md:111`).

### Evidence

| Claim in the recommendation | What the source says | Holds? |
|---|---|---|
| No simulation result in the registry | Every key begins `final_` or `first_four_months_` (my count) | Yes |
| Step 3a asserts 0.4.1 and has not been run | `pilot.R:15-16, 120-156`; no `step3*` output files | Yes |
| The code has only the smallest-P rule, the plan names the narrowest CI | `key.R:24-31`; `PLAN.md:119` | Yes, and the plan's rule was approved |
| Fix the rules, then run 3a | Title says the reverse | Title must change |
| Use 3a only to size and to drop contrasts | `pilot.R:147-151` prints the contrasts; `PLAN.md:51` keeps a failed claim in the table | Revise |
| Add the baseline-risk factor | `PLAN.md:116` approves a narrower variant; `03-simulation.md:26` gives a design reason | Needs a fair comparator and Andrew's decision |

### Proposed revision

> **Fix the simulation's claims, decision rules and scenarios in `simulation/PLAN.md`, then run the pilot analysis to size the full run.** Before any R job, Andrew decides and the plan records: the pass rule and the equivalence margin for claim 1; whether claim 3 reports the smallest-P rule beside the narrowest-CI rule he approved on 7 October, with the claim then scoped to the finding that selection on the estimate inflates the type I error while selection on interval width does not; whether the truth's baseline risk becomes a factor, replacing the approved family-baseline variants, and if so with an outcome-free events-per-covariate rule among the comparators and that factor reported separately; whether the nonlinear-treatment truth enters, which settles part of the TMLE proposal; and the rule for dropping a truth. Beside each rule, the plan records which pilot outputs had been seen when it was written.
>
> Then, as the first R job after the design stage, run `pilot.R` step 3a on the 0.4.1 library. Use it to size the full run, including the size each equivalence test needs. A contrast the budget cannot resolve stays in its table as unresolved, or is judged against the margin fixed beforehand, and its claim is rescoped; it is not removed. Scenarios the pilot did not run get their own small sizing run under the same rule.

## Recommendation 3: rules for the learner-defaults item, and the case study stays on 0.5.0

### Verdict

It survives with changes.

### What I tried

I checked whether keeping the case study on 0.5.0 contradicts a decision. It does not: the spec says "the re-execution stays on 0.5.0 as installed" (`2026-10-18-outcome-blind-design-stage.md:133`), so closing audit 04's rerun proposal (`04-case-study.md:190`) follows Andrew. It is also the opposite of tuning, since no fix can now change the case study's selection. The registry keys hold (`final_operators_sim_main_terms_ate_unresolved_cells` and `final_operators_sim_super_learner_ate_unresolved_cells` are 8; `final_operators_sim_datasets_ate` is 300), as do `simulate_design.R:16, 238` and the in-sample world score (`assess_design.R:103`, where every row is in one fold, so `fit_candidate.R:19-20` trains and predicts on all rows). The statement that every outcome model loads half on the score's logit holds for the linear form only; the nonlinear form adds a square and an interaction (`simulate_design.R:23`), which dilutes that share.

The release rule contradicts the spec. "Hold the 0.5.0 release only if the diagnostic finds a defect in the fitter or targeting" narrows decision 10, under which a defect the design stage exposes "goes into a patch release first" (`2026-10-18-outcome-blind-design-stage.md:125`), and the item's own text, under which a confirmed defect goes into a patch release before 0.5.0 is published and whether that holds the release "is Andrew's call" (`:133`). It also sits badly with Andrew's instruction that the decay grid and the treatment penalty "are to be fixed, not only tested" (`:141`): under the report's rule a confirmed decay-grid defect would not count.

"Today only the chosen point's risk is kept" (`learners.R:146`) is true for the path learners (the elastic-net penalty and the boosting rounds, `learners.R:138-146`) but not for the decay grid, which is what question 1 is about. Each decay is its own setting (`learners.R:12, 33`), and its mean held-out risk and its ensemble weight are kept (`learners.R:201-207`). What is missing is the per-fold spread, since the fold losses are summed (`learners.R:135-137`), and, in the simulation, the whole risk table, which is reduced to edge counts (`simulate_design.R:388-393`). Because the three decays enter the ensemble together rather than one being chosen, whether the net carries weight decides whether its edge matters (inference; audit 04 makes the same point at `04-case-study.md:83`), so the weights should be kept in the simulation too.

Judging every fix "on the simulation study's data-generating processes" alone may miss the problem. Those processes have four covariates and samples of 500 and 2,000 (`03-simulation.md:48`), while the behavior that opened the item came from `final_count_covariate_columns` covariate columns and training folds drawn from `final_operators_n` patients (`02-wrapper.md:169, 175`). The decay penalty scales with training size, and the penalized treatment is one column among 76 (`02-wrapper.md:175, 183`). Audit 02 proposed the simulation processes plus a wide variant with about 75 covariates (`02-wrapper.md:171, 179, 187`), and the report drops the wide variant without a reason.

The recommendation is silent on what the paper says if a patch ships. The paper describes the released version (`2026-10-18-outcome-blind-design-stage.md:116`), and a patch before publication would leave the case study on an earlier build. The decay grid also enters the propensity score of the selected estimator, because both candidates use the default library for the score (`22_design_stage.R:132`; `04-case-study.md:77`), so a decay-grid patch touches the estimator the case study applies, not only the candidate it rejected.

"Staged known-truth diagnostic" is a coined name; audit 02 describes the experiment as supplying the true nuisances to the targeting step one at a time (`02-wrapper.md:195`).

### Evidence

| Claim in the recommendation | What the source says | Holds? |
|---|---|---|
| Case study stays on 0.5.0; no rerun | Spec 2026-10-18 line 133 | Yes |
| Hold the release only for a fitter or targeting defect | Spec lines 125, 133, 141 | Contradicts the spec |
| Only the chosen point's risk is kept | True for path learners; decay settings keep mean risk and weight (`learners.R:201-207`); no standard error; simulation keeps only edge counts | Partly |
| World score is an in-sample fit; outcome models load on it | `assess_design.R:103`; `fit_candidate.R:19-20`; `simulate_design.R:16, 23` | Yes, half only in the linear form |
| Judge fixes on the simulation study's processes | Four covariates (`03-simulation.md:48`); audit 02 adds a variant as wide as the application (`02-wrapper.md:171`) | Too narrow |

### Proposed revision

> **Set the rules for the learner-defaults item before it changes any code, and keep the case study on 0.5.0 as installed.** Confirm the spec's decision: the existing design report governs, and the comparison between operators is not rerun on a patched build. Before any fix, keep for every learner setting its per-fold held-out losses (so its risk has a standard error), its ensemble weight and, for path learners, the risk along the path, in both the design assessment and the simulation; and run the experiment that supplies the true nuisances to the targeting step one at a time. Write each decision rule into the item's spec before its run.
>
> Judge every fix on simulated data with a known truth: the simulation study's data-generating processes and a wide variant with about 75 covariates and training sizes from a few hundred to several thousand, never on Rescue.Co results, with the Rescue.Co covariates used only to check that a setting runs. Whether a confirmed defect holds the 0.5.0 release stays Andrew's call under decision 10, and his instruction to fix the decay grid and the treatment penalty stands whatever the bias experiment finds. If a patch changes the default library before 0.5.0 is published, the paper says the case study ran on the earlier build and what the patch changed, including the propensity score the selected estimator uses.

## Recommendation 4: test and pin the in-house numerical code

### Verdict

It survives with changes.

### What I tried

The evidence holds. No test calls `.fit_Q()` (git grep of `tests/`); the edge-check test uses boosting rows only (`cleanTMLE/tests/testthat/test-learners.R:53-65`); the neural net is never fitted in the tests (it appears only in the settings count at `test-learners.R:13` and the plan defaults at `test-plan.R:5-6`); the ATE is pinned to a direct `tmle::tmle()` call and the ATT is not (`test-fit-candidate.R:119-134`); and tmle has no version floor (`cleanTMLE/DESCRIPTION:28`), with the `formals()` guard at `fit_candidate.R:59`.

The permutation test "as in 0.3.0" would test nothing in 0.5.0. In 0.3.0 the outcome store sat beside the lock and the test checked that design verbs ignored it (`cleanTMLE/NEWS.md:184-191`). In 0.5.0 no design-stage function receives an outcome: the lock's design data must not contain the plan's outcome columns and the simulation runs without one (`lock.R:27-29`; `test-blinding.R:8-12`). The risk in 0.5.0 is a design column that carries outcome information, since the lock keeps the outcome-observation indicators (`lock.R:132-134`) and the ledger records what the six-month indicator says about the outcome (`rescueCo/BLINDING_LEDGER.md:26, 43`). A useful test plants a copy of the outcome as an undeclared design column and checks that the lock drops it and that every output is unchanged when it is permuted.

Keeping the fitter in-house is presented as settled, but it bears on Andrew's second goal, a thin wrapper around established TMLE packages. Audit 02 concludes that the numerical layer cannot become much thinner without giving up the one-fitter guarantee and asks Andrew to confirm (`02-wrapper.md:150, 211`). The recommendation should leave that decision to him and say what follows for the framing: the package delegates targeting of the ATE, trimmed ATE and ATT to tmle, and owns its learner library, its cross-fitting and its overlap-weighted TMLE, so calling it a thin wrapper in the paper or vignettes would overclaim. The pins show agreement with established code, not that the estimators are valid; only the coverage check with oracle nuisances speaks to validity.

The WeightIt pins depend on a decision the report does not order: the parsimony audit proposes cutting both comparators (`audit/2026-10-22/01-parsimony.md:278`), while recommendation 1 asks whether the comparators run and recommendation 5 keeps weighting rows in a procedure comparison.

### Evidence

| Claim in the recommendation | What the source says | Holds? |
|---|---|---|
| No direct test of `.fit_Q()` | No match in `tests/` | Yes |
| Edge check tested on boosting only; net never fitted | `test-learners.R:13, 53-65`; `test-plan.R:5-6` | Yes |
| Permutation test as in 0.3.0 | `NEWS.md:184-191`; 0.5.0 never passes an outcome to design functions (`test-blinding.R:8-12`) | Redesign |
| ATT not pinned to tmle | `test-fit-candidate.R:119-134` | Yes |
| No tmle floor | `DESCRIPTION:28` | Yes |
| Keep the fitter in-house | `02-wrapper.md:150, 211` leaves it to Andrew | Decision for Andrew |

### Proposed revision

> **Test and pin the in-house numerical code before any cut or learner change.** Add a direct test of `.fit_Q()` (cross-fitting and the treatment column), glmnet and nnet cases for `.edge_check()`, and a test that a copy of the outcome planted as an undeclared design column is dropped by the lock and leaves every design-stage output unchanged when permuted. Pin `.nnloglik()` to `SuperLearner::method.CC_nloglik`, the glmnet and xgboost paths to `cv.glmnet(foldid = )` and `xgb.cv(folds = )` on the same folds, and ATT targeting to a direct `tmle::tmle()` call with and without tmle's refit of the propensity score; check the coverage of the overlap-weighted TMLE with oracle nuisances. Require tmle (>= 2.1.1).
>
> Andrew decides whether the fitter stays in-house. If it does, the paper and vignettes say that cleanTMLE delegates targeting of the ATE, trimmed ATE and ATT to tmle and owns its learner library, cross-fitting and overlap-weighted TMLE, and they do not describe it as a thin wrapper. Settle whether the comparators stay in the package before writing the WeightIt pins.

## Recommendation 5: restate the paper's claim after reading the nearest prior work

### Verdict

It survives with changes. The literature survey refutes its headline framing; the reading, the limits and the fallback paragraphs survive.

### What I tried

I checked the proposed headline against the literature. The survey in `08-literature.md` finds a source for each part of the workflow: staging with a clean room [@muntner2024]; outcome-blind simulation, run after data collection and before effect estimation, to choose the nuisance and variance estimation of a TMLE [@nanceCausalRoadmapSimulations2024a] (its abstract in `writing/references.bib` says so); a targeted-learning analysis plan [@gruber2022]; and gating on prespecified diagnostics before unblinding [cite: Conover 2025 objective study validity diagnostics]. It concludes that the claim that running the causal roadmap with TMLE inside a clean room "is itself an advance is not supported by this literature", and that the defensible advance is extending outcome-blind simulation from the estimator to the estimand, with a prespecified order, evidence on how often the verdict is right, and a package that makes the sequence checkable (`08-literature.md:7-13`). The report's recommendation points the other way, and the report does not draw on audit 08.

The proposed narrower gap, "no pre-outcome check tests whether an estimator delivers the estimand", is contradicted by a paper already in the bibliography [@nanceCausalRoadmapSimulations2024a] and by two the survey found [cite: Balzer 2026 prespecified data-adaptive TMLE in trials tested by plasmode simulation] [cite: Meng REFINE2 plasmode comparison of estimators for a prespecified ATE] (`08-literature.md:55, 63-64`). The gap the survey supports concerns the estimand: "in the uses we found, the estimand is fixed before the simulation is run" (`08-literature.md:160`).

I checked what evidence the joined-workflow claim would rest on. The simulation's main grid has one candidate (`simulation/R/map.R:14, 21`), so the estimator-choice half has no simulation evidence, and the super learner subset cannot be afforded (`03-simulation.md:46, 162`). The clean-room half is a software property with no evidence of bias prevented (`05-writing.md:184`), and the case study's approver had seen the outcomes (spec 2026-10-18 line 129). So the workflow claim could only be argued, and the paper's evidence speaks to the estimand verdict and to choosing before or after the outcome.

The recommendation contradicts a decision without saying so. Andrew's project record says that on 20 October he made "the advantage of TMLE with outcome-blind learner choice" a main point of the paper, with the proposal's three decisions still open (`C:/Users/andre/.claude/projects/C--Users-andre-OneDrive-Documents-clean-room-sim/memory/project_rescueco_reexec.md:13`; `writing/outline-results.md:13`). Audit 05 left both the headline and this question to him (`05-writing.md:408`). The report turns both into directives.

Replacing the TMLE claim with a procedure comparison risks the move the project rule forbids. The procedure comparison (`03-simulation.md:172, 239-250`) was proposed after the pilot's TMLE fits failed the tolerances at limited overlap (`03-simulation.md:158`; `pilot-report.md:99-101`). Swapping a claim the pilot contradicted for a comparison in which the package's procedure is expected to win is rescoping only if the contradicted result is reported as a finding. Inside the simulated family the procedure's coverage is also partly guaranteed, because the verdict and the label share the same tolerances (`03-simulation.md:52`), so its informative cells are the true outcome models that the simulated family cannot reproduce.

The limits the recommendation puts on the first page are right but incomplete. Hiding the outcome does not prevent bias when the outcome is predictable from the covariates on which the design stage chooses the estimand and the estimator (`08-literature.md:11, 81`) [cite: Fisher 2020 blinding outcomes is neither necessary nor sufficient]. The simulation has four covariates and two sample sizes, against `final_count_covariate_columns` covariate columns and comparisons of `first_four_months_operators_n` and `final_operators_n` patients in the application (`03-simulation.md:48`). And only the primary estimator applied is the one evaluated; the longitudinal TMLE and the comparators are not (`methods.qmd:120, 207`).

The evidence on the bibliography holds: `writing/references.bib` has no entry for Barnard, Conover, Liu 2025 or Karim, and `song2026`, `dang2023`, `gruber2022`, `gattoStructuredProcessIdentify2022`, `desaiProcessGuideInferential2024` and `nanceCausalRoadmapSimulations2024a` are present. Barnard's paper is in `writing/lit/liu_estimand2.jsonl`. The introduction's placeholder for Song (`manuscript/introduction.qmd:8`) can become `song2026` once its figure is checked (`05-writing.md:179`). Reading only four papers is too few; the survey lists more near work with verified identifiers (`08-literature.md:193`).

### Evidence

| Claim in the recommendation | What the source says | Holds? |
|---|---|---|
| Frame the claim as the joined workflow | `08-literature.md:13`: each part has a source; the advance is the estimand extension | Refuted as the headline |
| The gap may narrow to estimator checks | Nance abstract in `references.bib`; `08-literature.md:63-64, 160` | Refuted; the gap is the estimand |
| Do not make TMLE against weighting a main claim | Memory note line 13; `outline-results.md:13`; `05-writing.md:408` | Decision for Andrew |
| Limits on the first page | `08-literature.md:11, 81`; `03-simulation.md:48` | Incomplete |
| Four nearest papers missing from the bibliography | `references.bib`; `writing/lit/` | Yes |

### Proposed revision

> **Restate the paper's claim as an extension of outcome-blind simulation from the estimator to the estimand, inside a workflow that makes the sequence checkable, after reading the nearest prior work.** First read and add, through sci-cite, the near work that `08-literature.md` lists with verified identifiers, at least Barnard, Conover, Liu, Karim, Balzer, Meng, Haynes and Fisher. Then Andrew chooses the headline. The one the literature supports: outcome-blind simulation has been used to choose the estimator before effect estimation [@nanceCausalRoadmapSimulations2024a], and staging with a clean room restricts access to the outcome [@muntner2024]; this paper uses the simulation to judge which estimand a design can deliver, turns the verdict into a decision through a prespecified order, tests how often the verdict is right, and provides software in which the primary estimator applied after the outcome is joined is the one the simulation evaluated. The joined workflow is the setting, not the claimed advance.
>
> Put the limits on the first page: the verdict holds over the simulated family; the lock makes changes detectable, not impossible; nothing removes unmeasured confounding; withholding the outcome does not prevent bias when the outcome is predictable from the covariates the design stage uses; the simulation has four covariates and two sample sizes; and the case study's analysts had seen the outcomes. Write the Discussion's conclusion and a fallback paragraph for each simulation claim now.
>
> Whether the comparison of TMLE with weighting stays a main claim is Andrew's decision, taken on 20 October with three questions open. If the pilot's result that TMLE also fails the ATE and ATT tolerances at limited overlap holds, the paper reports it as a finding; a procedure comparison may be added as a separate claim, with the note that, when the true outcome model is one of the simulated ones, its coverage is partly guaranteed by the shared tolerances.

## Other problems

`REPORT.md` was not on disk at 22:58 on 22 October, so I could check only the five recommendations quoted in my task.

The report is described as built from audits 01 to 05 and does not use `08-literature.md`, which changes recommendation 5 and supplies the Fisher caveat. It should be folded in before the report goes to Andrew.

Audit 02 records that its author printed the first 40 lines of `rescueCo/logs/pipeline.log`, a log of the unblinded analysis whose header includes outcome counts, and says none of it entered a recommendation (`02-wrapper.md:264`). Audit 03 cites timing lines from the same file (`03-simulation.md:162, 350`). The report should carry this disclosure, the blinding ledger should record it, and the instructions to future auditors should list that log among the files not to open. I did not open it.

Script 25 writes the real estimates and design logs into `rescueCo/results/reexec_050/` (`25_estimate.R:25, 38, 42`). After unblinding, the folder that the project's audit rules call outcome-free would hold estimates. The output folder should change before script 25 runs.

Decisions about the comparators are spread across four places and nobody orders them: the parsimony cut (`01-parsimony.md:278`), recommendation 1 (whether they run), recommendation 4 (WeightIt pins) and recommendation 5 and the TMLE proposal (weighting "as the package implements it", `PROPOSAL-tmle-comparison-2026-10-20.md:30`). The cut should be decided first.

Two Methods sentences conflict with decisions or with the plan, and nothing in the five recommendations fixes them. The design-stage paragraph gives a bias tolerance of 0.02 (`manuscript/methods.qmd:74`), while the application's plan uses 0.0275 (`22_design_stage.R:134`; `methods.qmd:196`), so the paragraph should say that 0.02 is the simulation study's value. And `methods.qmd:199` still promises an approver and a registration on the Open Science Framework, which the governance decision rules out (spec 2026-10-18 line 129).

## Writing check

The check was run from the project root in Git Bash as `python "C:/Users/andre/OneDrive/Documents/scientific-writing-skills/checks/check.py" "audit/2026-10-22/07-check-framing.md" --min-level warn`. The outcome is recorded below.

The check ran. After my fixes it reports no banned words, no dashes, no bold lead-ins, no untraced numbers and no undefined abbreviations; the fixes replaced the bold verdict lines with subheadings, one word from the banned list, the spec names the check read as internal codes, two phrases from the project glossary, one number that is not in the registry, and three abbreviations. What remains is expected. Nearly every error flags a file path, line reference, registry key, function name or file name in backticks, which this internal audit needs and the task allows. Four errors flag the citation placeholders for Conover, Balzer, Meng and Fisher, which the citation rule requires where the project bibliography lacks a source. The one warning concerns the project configuration, not this file: the overlap source named in the writing configuration is a file in the unblinded analysis's repository, which I was told not to open, so that part of the check could not run.

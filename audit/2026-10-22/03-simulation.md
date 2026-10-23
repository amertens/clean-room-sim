# Audit 03: the paper's simulation study

Prepared 2026-10-22. This is a read-only audit. No R process was started, no `.rds` file or cache entry was opened (the cache was only listed), and nothing outside this file was written. Numbers come from `writing/registry.json` (the key is quoted) or from files cited by path and line; where a number is my own arithmetic on an assumed input, the text says so. "The verdict" means the output of `simulate_design()` for one estimand; "the estimability label" is whether some candidate meets the tolerances on the true outcome model. A Monte Carlo standard error (MCSE, or SE in running text) measures the noise from a finite number of simulated datasets, and a data-generating process (DGP) is a rule that draws a simulated dataset.

## Overall judgment

The planned study is a sound validation of one feature of the package, and its machinery is good: seeds come from names, results are cached per cell, design-stage runs are shared across outcome models so comparisons are paired, and a pilot gate exists. As planned, it cannot make the strongest honest case, for five reasons. The pilot already holds the data for its main contrasts (500 cached design-stage runs per variant and 6,000 label fits at n = 500) but nobody has analyzed them, so no promised contrast has been checked for resolvability. The scenarios vary overlap, the one factor the comparator diagnostics read, so the claim that the verdict beats those diagnostics may be untestable, and the claim that agreement falls when the family cannot reproduce the true model can appear only through the labels because the verdict never sees the outcome. Every scenario gives the estimators a correct main-terms propensity model, so outcome models the family cannot reproduce may leave the labels of a doubly robust estimator unchanged, and nothing in the grid needs a flexible learner. Claim 3 as reworded (narrowest CI) is expected to show no harm, and the TMLE comparison pits TMLE against weighting rather than testing the clean room, with an expected ordering that the pilot's own fits already contradict. The super learner subset in `manuscript/methods.qmd:158` cannot be run. A smaller design with four main displays and five supplementary ones carries the claims that survive. Moving the code to 0.5.0 changes two call sites and no numerics.

## 1. Claims and aims

### 1.1 The plan against the skill's stages

| Stage | Status | Evidence |
|---|---|---|
| 1. Claims | Four claims, prioritized. Claim 1 holds two claims. No claim states what would count as met. | `simulation/PLAN.md:7-14` |
| 2. Aims | Three goals with subclaims by scenario; comparators pinned; the label is the oracle; a semi-oracle is optional. No pass criteria. | `PLAN.md:16-34` |
| 3. Mockups | Written before code, but stale after the 2026-10-07 revisions: Table 1 is still pooled, Figure 1 still has three overlap levels, and the rows promised for the baseline-risk variants are missing. | `PLAN.md:40`, `:57`, `:116` |
| 4. DGPs | Reasons given for forms and for the overlap levels; none for several coefficients. Diagnostics lack the nonlinearity measure and the working-model bias check. | `PLAN.md:84-103`, `simulation/R/dgp.R:176-190` |
| 5. Build | One pipeline, subsets run alone, cache keyed on version and code text. Per-repetition results are discarded, which the skill says to keep. | `simulation/R/run.R:14-51`, `simulation/R/map.R:109-111` |
| 6. Pilot | Steps 1 and 2 reported, and `pilot-report.md:3` still says step 3 has not started; it did run (cache listing, launch logs), but its analysis step never did (no `step3.rds` or `step3_*.csv` in `simulation/pilot/`). | `simulation/pilot.R:105-156`, `simulation/pilot/step3_resume2_stdout.txt` |

### 1.2 Claim 1: the verdict against overlap diagnostics, and its boundary

Claim 1 (`PLAN.md:9`) is two claims. Claim 1a says the verdict agrees with the label more often than overlap diagnostics when the truth lies inside or near the family. Claim 1b says agreement falls when the truth lies outside. `PLAN.md:24` asks the verdict to beat every best-threshold diagnostic, which is four diagnostics, four estimands, two sample sizes and two measures; two readers could disagree about whether winning most of 64 comparisons meets it. A rule should be written down before the full run. For each of the ATE, trimmed ATE and ATT, the verdict is better when its balanced accuracy exceeds the best tuned diagnostic's by at least 3 Monte Carlo SEs, no worse when the lower 95% bound of the difference exceeds -0.05 (a margin for Andrew to set), and unresolved otherwise, with a sentence prepared for each outcome.

Claim 1a may be untestable as designed. The scenarios cross overlap with outcome shape, but the outcome shape changes only the label, never the verdict, and the diagnostics are all monotone summaries of the one thing that varies. A threshold tuned in-sample by Youden's index (`simulation/R/summarize.R:52-64`, `:83-103`), separately for each estimand and sample size, can then track the label up to dataset-level noise. That noise is real at n = 500 (two datasets at one level had a share outside the band of 0.208 and 0.148, and an effective sample size fraction of 0.21 and 0.44, in `simulation/pilot/step2_output.txt`), and the verdict is noisy too (`simulation/pilot/pilot-report.md:126`), so I cannot say who wins. I can say that the design gives the verdict no condition in which a propensity-only summary must fail. Two changes supply one. First, add a factor that no propensity summary can see, the outcome's baseline risk, which sets how many events the cohort holds; the plan's family already takes a baseline risk as an input (`cleanTMLE/R/plan.R:239`), and the package already reports expected events before the outcome is seen (`cleanTMLE/NEWS.md`, 0.4.1). Second, compare against diagnostics as an analyst could apply them: at the package's own grade cut-offs (`cleanTMLE/R/assess_design.R:111`), and at a threshold tuned on the other sample size, not only at the best in-sample threshold. If the verdict still does not beat a tuned diagnostic, the honest rescoping is that it matches the tuned diagnostic without needing a threshold that no analyst can tune, because the true outcome model is unknown.

Claim 1b is testable only through the labels. The verdict depends on (overlap level, n, dataset) alone, since covariates and treatment are drawn per level and the verdict never sees the outcome (`simulation/R/dgp.R:192-210`, `simulation/R/summarize.R:68-78`). At a given level the inside, near and outside models therefore receive the same verdict, and agreement can differ across them only where their labels differ. For a label to differ, the same estimator must have different bias or coverage under different outcome models. I infer, without having seen the labels, that this is unlikely for the planned outside model. TMLE is consistent when either nuisance estimator is [cite: TMLE consistency when either the outcome regression or the propensity score is estimated consistently], and every planned propensity model is exactly the main-terms logistic model the candidate fits (`dgp.R:38`, `simulation/R/map.R:19-24`). The outside model's two features, effect heterogeneity in a covariate unrelated to treatment and a threshold in the treated-dominated region (`dgp.R:70-71`), act on the treated arm's outcome, where data are plentiful. The control-arm outcome, which the ATE must extrapolate into the region of thin support, is linear in the first covariate in both the inside and outside models and curved only in the near model (`dgp.R:65-71`, the term `0.6 * X1^2`). So the model most likely to move labels is the one the plan calls closest to the family. The second half of claim 1b has two defences. Read the labels first (section 7, recommendation 1). Then, if labels do not differ, add the planned nonlinear-treatment scenario from the TMLE proposal (`PROPOSAL-tmle-comparison-2026-10-20.md:56`), the only planned truth where a main-terms estimator is wrong on both nuisances, and state 1b as a comparison with an equivalence margin. A row with an unmeasured confounder (`PROPOSAL:57`) gives one certain disagreement, since the verdict cannot see the confounder, and it states the limit the Methods already name (`methods.qmd:45`).

If claim 1b stays a main claim, one ablation attributes any fall to the family and not to resampling noise or to the fitted propensity score: run the verdict again with the family extended by the true outcome model as a custom surface (`cleanTMLE/R/plan.R:253`; the four default surfaces cannot be dropped, so the family becomes those four plus the truth) for the near, outside and nonlinear truths at two levels. If agreement recovers, the family is the cause. This is the skill's semi-oracle for a mechanism claim, and it applies only to truths that are functions of the covariates alone, which excludes the inside model because its index uses the true propensity score. If 1b becomes a supplementary sentence, skip it.

### 1.3 Claim 2: new treatments against observed treatments

The mechanism is published [@shaw2026]: keeping each patient's observed treatment carries regions without one arm into the simulated data, so the variant should call estimands infeasible more often, which is a one-sided prediction. The pilot's smoke runs showed the two variants disagreeing on the trimmed ATE in both n = 500 datasets, in opposite directions (`pilot-report.md:91`), so the average gap may be near zero. The variant costs about 1.25 times the generate run (400 against 320 seconds at n = 500, `pilot-report.md:64-65`) and was the slowest pilot stage. Keep claim 2 only if the pilot's paired contrast on the trimmed ATE clears |z| of 2; otherwise reduce it to one sentence that cites the published mechanism and say the study could not separate the variants.

### 1.4 Claim 3: choosing the estimand before the outcome

The plan changed the outcome-informed comparator from the smallest P value to the narrowest CI after the Methods review (`PLAN.md:119`; `methods.qmd:157`), but the code still implements the smallest P value (`simulation/R/key.R:24-31`, `simulation/R/summarize.R:186-211`). The two rules differ in what they select on. The width of a Wald interval is close to independent of the point estimate, so selecting the narrowest interval should leave coverage and the type I error near nominal; the harm the claim describes arises from selecting on the estimate. I infer that the narrowest-CI rule will show no gap at any sample size, which would look like a failure of claim 3 and is a feature of the comparator. Under the rule that selects on the estimate, the gap is small because the four estimands are highly correlated estimators of nearly the same quantity (arithmetic in section 3). Changing back to the smallest P value after seeing a null with the narrowest-CI rule would be shopping for a comparator, so both rules should be fixed now and both reported, with the claim scoped to the finding that selection on the estimate inflates the rejection rate and selection on width does not. The contrast also exists only where two or more estimands are feasible: in the smoke runs at the moderate and severe levels the ATO was feasible every time, the trimmed ATE in some runs, and the ATE and ATT never (`pilot-report.md:91`), so at those levels the prespecified and outcome-informed choices coincide or have only two options.

### 1.5 Claim 4: forecast against realized half-width

The forecast is the largest, across the four surfaces, of the mean half-width of the selected candidate (`cleanTMLE/R/simulate_design.R:249-250`), so the ratio of realized to forecast should sit below 1 by construction. A correlation across scenarios (`PLAN.md:34`, `simulation/R/summarize.R:258`) is dominated by the spread of widths between overlap levels and sample sizes and will be high whatever the quality of the forecast. The display should show the ratio's distribution within each level and estimand, state the expected direction, and set a calibration target before the run. The Rescue.Co plan uses the half-width forecast as a feasibility gate (the ATT was infeasible at 0.0551 against a tolerance of 0.055, `writing/digest.md:23`), which is the strongest reason to keep this claim.

### 1.6 What the claims leave out

The paper's thesis is that the estimand and the estimator are chosen before the outcome. Claims 1 to 4 concern the estimand only. The main grid has one candidate (`simulation/R/map.R:14`, truncation 0.01 only, main-terms learners), so selecting a candidate is vacuous in every scenario except the super learner subset. The proposal supplies an estimator claim, discussed in section 4.

Three scope gaps should be stated in the paper whether or not they are closed. The simulation has 4 covariates and sample sizes of 500 and 2,000; the application has 75 columns (`final_count_covariate_columns` = 75) and comparisons of 226 (`first_four_months_operators_n`) and 1,204 (`final_operators_n`) patients, and the main-terms candidate's worst-case bias in the 226-patient comparison was 0.050 to 0.067 (`first_four_months_operators_sim_main_terms_worst_bias_min` and `_max`), against biases of up to 0.039 in the pilot's n = 500 fits (`simulation/pilot/key_check_output.txt:3-11`). The simulation uses `ps_library = "glm"` (`map.R:20`) where the package default and the application use the super learner. The simulation runs `max_reps` equal to `reps` (`map.R:23`), so no unresolved cell receives extra batches, while the application's verdicts used the extension (`final_operators_sim_datasets_ate` = 300 after two batches, and every cell still unresolved, `writing/digest.md:23`). Each is a limit on what the simulation's agreement figures say about the application.

### 1.7 Where the comparisons favour the package, where they do not, and the tuning risks

Four choices favour the package. The family is given the true baseline risk, which the Methods acknowledge (`methods.qmd:153`) and the misspecified variants are meant to offset. The inside model is built from the family. The verdict is estimator-aware, since it uses the same fitter, folds and plan seed as the label fits (`simulation/R/key.R:1-3`), while every diagnostic is estimator-blind, so part of any advantage is that the verdict knows the estimator. And the label is defined by the verdict's own tolerances, which is the claim but makes the agreement partly definitional. Three choices work against the package: the diagnostics' thresholds are tuned in-sample against the labels, and the best of the four is picked by agreement in the same cell (`simulation/R/summarize.R:152-154`), which is generous to the comparator; the verdict is a conjunction over four surfaces (section 3.3); and no unresolved cell receives extra batches. The first group does not cancel the second, so the paper should state both.

Four moves would amount to tuning the design until a claim survives, and each should be ruled out in `PLAN.md` before the full run. Adding overlap levels after seeing the contrast pilot's results: the two already added were disclosed, and no further change should follow. Changing the claim 3 comparator after seeing a null: fix both rules now. Redefining the outside model, or dropping a truth, after seeing the labels: state the rule in advance (section 5.2). Choosing the equivalence margin or the pass rule after the full run: choose them from the pilot's sizing record, before it.

## 2. Data-generating processes (DGPs)

### 2.1 What the skill finds sound

The covariates are bounded and few; the treatment model has one coefficient that sets overlap and has been calibrated on a population draw of 10^6; the true value of each estimand is computed by brute force; variance explained, the c-statistic, prevalence and the share outside the band are reported (`simulation/pilot/dgp_diagnostics.csv`); the control-arm risk is the same in every scenario (0.150); and datasets are seeded from names (`dgp.R:23-28`, `:198-217`). The treatment model leaves the skill's default range for the linear predictor of about -3 to 3 on purpose, since overlap is the factor under study, and the skill's default of separate training, validation and estimation draws in place of cross-fitting does not apply, since cross-fitting is part of what is evaluated. Adding two overlap levels after seeing the pilot was disclosed (`PLAN.md:111`), and a design chosen to make a mechanism visible is allowed by the skill if it is reported. The five levels should now be frozen.

### 2.2 What the skill finds missing or unsafe

| Item | Problem | Evidence | Fix |
|---|---|---|---|
| One manipulated factor | Overlap is the only factor that moves the label, and every comparator diagnostic reads it (section 1.2). | `dgp.R:112` | Add baseline risk as a second factor for the three non-trivial levels at n = 500 (section 5). |
| Correct propensity model everywhere | A main-terms propensity model is exactly right in every scenario, so a doubly robust estimator is robust to the outcome model's shape. The design-stage world also draws treatment from a main-terms fit, so a nonlinear true propensity score is invisible to the verdict. | `dgp.R:38`, `cleanTMLE/R/assess_design.R:103` | Add the nonlinear-treatment scenario; report that the verdict cannot see propensity misspecification when the plan's learners are main-terms. |
| Cancellation check missing | Nobody has computed the asymptotic bias of the main-terms TMLE on a very large draw for each truth. If it is below 0.02, no number of repetitions will show a failure. For a truth whose propensity model is correct I expect it to be about zero (inference from double robustness), which would settle that the first three truths can differ in label only through finite-sample behavior. | skill, "Cancellation"; `dgp.R:176-190` | One fit per truth at n = 10^6 before any replication; amplify the truth if it is below the tolerance and say so. |
| Nonlinearity not reported | The skill's fourth diagnostic is absent. | `dgp.R:176-190` | Add the share of Var(mu) captured by the best linear approximation. |
| "Inside" is built from the family | The inside model reproduces the package's linear surface with heterogeneity 0.5, except for the covariate direction and a baseline risk set to match (`dgp.R:54-58`, `:65-68`; `pilot-report.md:123`). It favours the verdict by construction; the near and outside models were chosen by the author. | `dgp.R:54-68` | State that "outside" means these specific models and that the baseline risk is given the true value. |
| The family's direction is seeded per dataset | The family's covariate direction and its squared covariate come from the plan seed (`cleanTMLE/R/simulate_design.R:15-23`), which here differs by dataset (`map.R:12`). Verdicts therefore vary with the seed as well as the data. | as cited | Report it as a source of verdict variation; it is realistic. |
| Null scenarios double-count | `run_table` merges all six model-by-effect scenarios of a level into Table 1, so each verdict is counted up to six times with possibly different labels. | `summarize.R:68-78` | Table 1 uses effect scenarios only; null scenarios serve claim 3. |

### 2.3 What the study must include for true models the family cannot reproduce (learner-defaults item, question 4)

The default family is logistic in a smooth index with one square and one interaction (`cleanTMLE/R/simulate_design.R:23`), and the design stage selects, among the candidates that pass, the one with the smallest worst-case root mean squared error across the family (`:238`). A family that main-terms regression fits well therefore tends to select main-terms regression whether or not the real surface is smooth. The Rescue.Co design stage selected main-terms logistic regression with truncation 0.05 for the ATE (`writing/digest.md:23`). To say anything about true models the family cannot reproduce, the study needs the following.

| Requirement | Why | Where it sits in the minimal design |
|---|---|---|
| A truth, for outcome and treatment, that main-terms logistic regression cannot represent and whose misspecification moves the estimate by more than the tolerance | Without both nuisances wrong, double robustness hides the difference; without a bias above 0.02 there is nothing to detect | The nonlinear-treatment scenario, after the large-draw bias check |
| Two or more candidates in the label, reported separately | "Some candidate passes" can hide that only the flexible one does | Key fits for main-terms, a reduced library (main-terms plus MARS), and the default super learner on a few hundred datasets |
| Truths on both sides of the family, crossed with the library | Any claim that the super learner helps or hurts depends on which side the DGP is on | Table S3: smooth in-family and strongly nonlinear, by library |
| The verdict's selected candidate on the strongly nonlinear truth | The point of the question is whether a family that favours main-terms still selects the flexible candidate when it is needed | A small number of design-stage runs with two candidates at one level and n = 500 (section 5.4) |
| The truth written as a custom surface | Lets the same function become the package's example of a truth outside the declared surfaces, which `cleanTMLE/WISHLIST.md` lists | `surfaces$custom` (`cleanTMLE/R/plan.R:253`) |
| DGPs fixed before running, justified by the claim | The spec bars tuning to Rescue.Co results | Freeze the list in `PLAN.md` before the full run |

Covariate dimension is the remaining gap: the application's main-terms candidate has 75 columns against 226 to 1,204 patients, and a scenario with about 30 covariates and 5 active ones would test the regime in which the family's preference matters most. It adds code (the DGP is written for four named covariates) and I place it after the paper, with the paper stating that the simulation has four covariates.

### 2.4 The code against the skill's build stage

The pipeline has the properties the skill asks for. One run function serves the pilot and the full run, with a subset of units, sample sizes and dataset indices (`simulation/R/run.R:58-84`). Each dataset's seed comes from its name and not from a loop position (`simulation/R/dgp.R:23-28`, `:198`). Each cell is written atomically under a key that includes the cleanTMLE version and the text of the code files (`run.R:24-51`), a cache hit was shown to equal recomputation, and parallel rows were shown to equal sequential rows (`simulation/pilot/pilot-report.md:58`, `simulation/pilot/parallel_test_output.txt`). The Monte Carlo SEs are formed from dataset indices, which is right when one design-stage run serves six scenarios (`summarize.R:1-5`).

Four things fall short. The raw per-repetition results are discarded (`map.R:109-111`), so no verdict can be recomputed at 100 repetitions, under a different tolerance or under the extension rule without rerunning it. The cache key leaves out R, `tmle`, WeightIt and MatchIt, although `tmle::tmle` does most of the work. The observed-treatment variant is a private copy of a package function (`map.R:26-62`) and ten other internals are called directly (section 6.2). And the summary functions that produce Table 1, claim 3 and claim 4 were run end to end only on two datasets (`pilot-report.md:58`), so the first run of the analysis step on 6,000 label fits and 1,000 verdicts may need debugging; that is a reason to run it early and a risk to the one-hour estimate in section 7.

## 3. Contrasts and Monte Carlo error

### 3.1 Method and what is missing

For each contrast the skill reduces the comparison to a per-dataset difference d and sizes the run at `n_sim >= (k sd(d) / mean(d))^2`, with k = 3 to see a gap and 5 to be comfortable, using a lower bound on the gap and treating a pilot with |z| below 2 as too small to size anything (`simulation/R/pilot_check.R:98-146`). Here the unit is the dataset index, because one design-stage run serves every outcome model at its level, so per-index means are formed first (`summarize.R:1-5`). The pilot's contrast data are in the cache (5 levels x 100 verdicts per variant, 30 scenarios x 200 label fits, all at n = 500) but the analysis `pilot.R 3a` (`pilot.R:120-156`) has not produced output. I therefore could not use the pilot's gap, `sd_d` and |z|. The inputs below are measured where a text file exists (the 40-dataset label check in `simulation/pilot/key_check_output.txt`) and otherwise assumed and labelled; running 3a replaces every assumed input.

### 3.2 Label precision

The standard deviation of an estimate is the bias MCSE reported at 40 datasets times the square root of 40 (`key_check_output.txt:3-11`, `:13-22`, `:24-33`). For example, moderate overlap, inside truth, n = 500: ATE 0.0564, ATT 0.0943, ATO 0.0413, trimmed ATE 0.0451.

| Cell (n = 500 unless stated) | Bias at 40 datasets | MCSE at M = 1,000 | Distance from 0.02 in MCSE, M = 1,000 | Same, M = 5,000 |
|---|---|---|---|---|
| Moderate, inside: ATE | 0.0218 | 0.0018 | 1.0 | 2.2 |
| Moderate, inside: trimmed ATE | 0.0171 | 0.0014 | 2.0 | 4.5 |
| Moderate, inside: ATT | 0.0352 | 0.0030 | 5.1 | 11 |
| Severe, outside, n = 2,000: ATT | 0.0202 | 0.0026 | 0.1 | 0.2 |
| Good, inside: ATT | -0.0149 | 0.0021 | 2.5 | 5.6 |

Coverage has the larger label noise. The MCSE of a coverage near 0.90 is 0.0095 at M = 1,000 and 0.0042 at M = 5,000 (formula in the skill's `performance-measures.md`), so cells with true coverage within about 0.019 of the 0.90 tolerance are unresolved at M = 1,000 and within 0.008 at M = 5,000. A label fit costs about 0.66 worker-seconds at n = 500 (6,000 cells took 664 seconds on 6 workers, `simulation/pilot/step3_run_output.txt`), so M = 5,000 for the 30 effect and null scenarios of the pilot is about 28 worker-hours. I recommend it, with the report of cells within two MCSE that `PLAN.md:117` already promises. Some cells will stay unlabellable at any affordable M (the severe ATT above); the sensitivity analysis excluding them is the right treatment.

### 3.3 Noise in the verdict itself

The verdict is a conjunction. A candidate counts only if it passes on all four surfaces (`cleanTMLE/R/simulate_design.R:197`), and at 200 repetitions each pass is a noisy event. With a true coverage of 0.90 in each surface, a single surface passes the point-estimate rule with probability 0.56 (binomial arithmetic, 200 repetitions); four independent surfaces all pass with probability 0.10, and perfectly dependent surfaces with 0.56. At a true coverage of 0.92 the figures are 0.59 and 0.88, and at 0.88 they are 0.00 and 0.23. For bias with an estimate SD of 0.05, a true bias of 0.017 passes with probability 0.80 and 0.023 with 0.20. So near a tolerance the verdict is conservative, which is a property of the tool at this budget, and sensitivity will be lower than specificity there. The paper should say so before a reader finds it.

### 3.4 The planned contrasts

The per-level paired difference in correct-call rate between two sources has variance q - g^2, where q is the share of datasets on which the two sources disagree and g the gap. The standard error at R verdicts is sqrt((q - g^2) / R). The q and g values below are assumed.

| Contrast (where) | Inputs | MCSE at the planned size | Difference the claim needs | Resolved? |
|---|---|---|---|---|
| Verdict against best diagnostic at one level (Table 1; R = 100 in the pilot, `pilot.R:106-107`) | q = 0.2 or 0.4; g = 0.10 | 0.044 or 0.062 (z = 2.3 or 1.6) | g at 3 SE: 0.13 or 0.19 | No at R = 100 for g below about 0.15; at R = 300, z = 4.0 or 2.8 |
| Same, g = 0.05 | q = 0.2 | 0.044 | 3 SE needs R = 712; 5 SE needs 1,975 | No at any planned R |
| Pooled over five levels, trimmed ATE, two informative levels (assumed q = 0.4, g = 0.10) | `sd_d` = 0.177, gap 0.04 | 0.0177 at R = 100; 0.0102 at 300; 0.0079 at 500 | R of 176 for 3 SE, 488 for 5 SE | No at R = 100 (z = 2.3); yes at 300 for k = 3 |
| Inside against outside (Table 1) | labels differ at few levels; q about 0.2, g about 0.05 | as above | R of 712 or more | Probably not |
| Generate against observed treatment, trimmed ATE | smoke: both datasets disagreed (`pilot-report.md:91`); q about 0.5 | 0.049 at R = 200 for g = 0.15 (z = 3.1) | g = 0.10 needs R = 441 (k = 3) | Only for g of 0.15 or more at R = 200 |
| Claim 3, rejection under the null, smallest P value among 4 feasible estimands | correlation among estimators 0.70, 0.85, 0.95 (assumed) | gap 0.081, 0.055, 0.030 | required repetitions 103, 155, 292 (k = 3); 284, 430, 810 (k = 5) | Yes at k = 3 with R of 300 if four estimands are feasible; k = 5 needs up to 810 |
| Claim 3, same with 2 feasible estimands | same correlations | gap 0.034, 0.025, 0.015 | required repetitions 259, 358, 603 (k = 3) | Marginal; at the severe levels the choice does not exist |
| Claim 3, narrowest CI | width nearly independent of the estimate (inference) | gap about 0 | none | Structurally null |
| Claim 4, ratio of realized to forecast half-width | assumed sd of the ratio 0.15 | 0.015 at R = 100 | departure from the target of 0.05 | Yes, for a ratio; no for a correlation across levels, which is high by construction |
| Proposal (b), coverage of weighting against TMLE (`PROPOSAL:64`) | M = 1,000, q about 0.25, g = 0.05 | 0.016 | 3 SE needs 0.05 | Yes with the main-terms TMLE arm |
| Proposal (b), same with the super learner chosen by the design stage per dataset | M about 20 to 100 (section 4) | 0.05 to 0.11 | gap of 0.15 to 0.33 | No |
| Proposal (c), bias of weighting against TMLE | paired; gap above 0.03, sd of the paired difference about 0.03 (assumed) | 0.003 at M = 100 | gap of 0.03 | Yes, if the gap is real (cancellation check first) |

### 3.5 Choosing R and M, and what it costs

The size of R is set by the largest requirement among the contrasts the paper keeps; contrasts that need more than the budget are dropped from the claims, not rescued by changing the design (the skill's rule, and the proposal's own, `PROPOSAL:64`). Costs below use the pilot's measured seconds per repetition with one truncation (2.12 at n = 500 and 3.63 at n = 2,000, `simulation/pilot/parallel_test_output.txt`), times 200 repetitions, with 1.25 times that for the observed-treatment variant (`pilot-report.md:64-65`), and the measured speed-up of 4.94 on 6 workers. They agree with the pilot's actual run, which took on the order of a day for 500 verdicts in each variant at n = 500 (`simulation/pilot/step3_resume_stdout.txt`, `step3_resume2_stdout.txt`).

| Stage | Seconds per unit | R = 100 | R = 300 |
|---|---|---|---|
| Generate verdicts, 5 levels, n = 500 and 2,000 | 424 and 726 | 32 hours wall | 97 hours (4.0 days) |
| Observed-treatment verdicts added, same grid | 530 and 910 | 41 hours more | 121 hours more |
| Label fits, M = 5,000, 40 scenarios at each sample size | 0.66 and about 1.4 (the second assumed) | n/a | 23 hours wall (114 worker-hours) |

The observed-treatment variant at both sample sizes is the least informative use of compute. I recommend running it at n = 500 on three levels (section 5.4).

## 4. The TMLE comparison proposal

### 4.1 Judgment

The proposal (`simulation/PROPOSAL-tmle-comparison-2026-10-20.md`) is well built as a plan, with comparators pinned, a stated boundary and a rule that a subclaim is dropped when its contrast cannot be resolved (`:64`). It has four problems.

It tests TMLE against weighting and matching, not the clean room. The justification (`:11`) is that flexible nuisance estimation is where an analyst has most freedom, and it is defensible only if made blind. Every comparator is prespecified, so the comparison cannot show what blindness buys. The arm that would, a candidate chosen after seeing the outcome, is absent.

Its expected ordering under limited overlap is already contradicted. Subclaim (b) (`:22`) expects TMLE to keep closer to the tolerances for the ATE and ATT at limited overlap. The pilot's fits of the main-terms TMLE at moderate overlap and n = 500 gave an ATE bias of 0.0218 with coverage 0.85 and an ATT bias of 0.035 to 0.039 with coverage 0.68 to 0.70, and at severe overlap and n = 2,000 an ATE bias of 0.0125 with coverage 0.83 (`pilot-report.md:99-101`, `key_check_output.txt:3-22`). TMLE fails these tolerances there. The result that follows is that at limited overlap no estimator meets the tolerances for the ATE or ATT, and what rescues inference is moving down the ladder, which the proposal's Table 2 hides by counting estimator-by-estimand combinations (dominated by the ATO, which every estimator passes).

Two of its rows coincide. With one candidate in the grid (`map.R:14`), "TMLE with the candidate selected before outcomes" is the main-terms candidate in every scenario except the nonlinear one, so rows 1 and 2 of Table 2 are identical in three of four columns.

Its headline arm is unaffordable. Choosing the candidate per dataset needs a super learner design-stage run on each dataset. The Rescue.Co pilot took 889 seconds per repetition at n = 1,204 with 75 covariates (`rescueCo/logs/pipeline.log:1204`), and 2,251 seconds at n = 6,494 (`:1205`). Even if four covariates and n = 500 make a repetition ten times cheaper, a run of 200 repetitions takes about 5 hours, so 100 datasets would take 500 worker-hours for one level and one sample size. Only about 20 to 30 datasets are affordable, which resolves nothing smaller than a gap of about 0.3 at 3 SE (paired standard deviation about 0.5, assumed).

The comparison also stacks the odds toward TMLE in two places: the weighting comparator uses a main-terms propensity model while the selected TMLE gets a super learner in scenario (c), and the matching comparator estimates the ATT of matched patients only. A super-learner-weighted comparator is needed to attribute any advantage to TMLE and not to flexibility. The package's weighting comparator is a main-terms WeightIt model with the M-estimation standard error (`cleanTMLE/R/estimate_effect.R:24-48`), which I infer accounts for fitting the weights (WeightIt's default; not checked against the installed version), so its coverage under limited overlap is an empirical question and "undercovers" should not be assumed.

### 4.2 Recommendations on the three decisions

#### Decision 1: is the claim stated as you want it, with its boundaries?

No. The boundaries (no advantage at good overlap with linear models; no fix for unmeasured confounding) are right and should stay, but the claim's core is wrong in two respects and should be rewritten as a claim about the procedure. I recommend:

> Across overlap levels, the procedure the package supports (the estimand chosen from a prespecified order by the outcome-blind simulation, the estimator fixed in the plan) reports an interval that covers the true value of the estimand it reports at close to the nominal rate, while an analysis that fixes the protocol's estimand and an estimator does not once overlap is limited; the price is that the reported estimand differs from the protocol's in a stated share of datasets. Where the protocol's estimand is estimable, targeted maximum likelihood estimation and weighting [agree | differ] by an amount the study will report. When an unmeasured confounder is present, every procedure is biased alike.

The bracket is left open because the pilot cannot yet say which holds. The reason is that the comparison's value is at the level of the decision (which estimand, with which estimator, with what coverage), and the estimator-only comparison cannot show that.

#### Decision 2: is the comparison of TMLE candidates (selected against main-terms) part of the main claim?

Supplement, and restricted. The selected arm cannot be run per dataset at any resolvable size, with the default super learner as the flexible candidate. Rows should be defined by candidate: main-terms, a reduced flexible library, and the default super learner on a few hundred datasets, all fitted once per dataset like the label fits (cheap relative to a design-stage run), plus the design stage's selection on 20 to 60 datasets of the nonlinear scenario, reported as a selection rate with its interval. This is also the evidence the open learner-defaults item needs (section 2.3), so it should be designed once for both.

#### Decision 3: Figure 3 and Table 2 as proposed in the main text?

No. Table 2 should be replaced by the procedure comparison in section 5.3 (coverage of the reported estimand, share of datasets reporting the ATE, share with no feasible estimand), where each strategy is judged by the estimand it reports and the price of moving is visible. Figure 3 (bias and coverage by estimator, estimand and scenario) goes to the supplement as the full-results display; matching, which targets a different population, belongs there too.

#### Additions the proposal needs

The proposal also needs an outcome-informed arm (the estimator or estimand chosen after the outcome); weighting with the same flexible propensity learner as the TMLE arm; the semi-oracle TMLE (true nuisances) in the nonlinear scenario only, as a check that the DGP is not the problem; the large-draw bias check for the main-terms TMLE before any replication (the skill's cancellation test); and an equivalence margin for the "the three agree" sentence (`:9`), because agreement is a claim of no difference and needs one (the skill's rule, `piloting.md`, "Claims of no difference").

## 5. The proposed minimal design

### 5.1 Claims the design carries

| Claim | Display | Condition under which it fails and is reported as such |
|---|---|---|
| 1a. The design-stage simulation's verdict matches the estimability label at least as well as a tuned overlap diagnostic, and better than diagnostics at conventional or transported thresholds | Figure 1, Table 1 | A tuned diagnostic matches it everywhere; then the claim becomes "without needing a tuned threshold" |
| 1b. Agreement does not fall (or falls by a stated amount) when the true model is one the family cannot reproduce | Figure 1, Table 1 outside block | Labels never differ across truths; then 1b is reported as untestable with these truths |
| 2. Keeping observed treatments changes the verdict in the direction published | Table 1 row, n = 500 | Contrast unresolved at the planned R; one sentence |
| 3 and the proposal. The procedure keeps coverage for the estimand it reports, at the price of reporting a different estimand in a stated share of datasets; selection on the estimate inflates the rejection rate | Table 2, Table S1 | Both outcome-informed rules show no gap; then the claim is "selection on width is harmless; selection on the estimate is not" or is dropped |
| 4. The half-width forecast is calibrated and conservative | Figure 2 | Calibration outside the target |
| Scope. Truths that need a flexible learner | Table S3 | None; it states a limit |

### 5.2 Scenarios

| Scenario | What varies | Tests | Cost |
|---|---|---|---|
| Five overlap levels (the plan's gammas 1.25, 1.5, 1.75, 2.0, 2.5; share outside 0.02%, 4.2%, 16%, 26%, 41%), frozen | Overlap | All claims; the first level is the base case in which every method should work | Design-stage runs at n = 500 and 2,000 |
| Truths: inside, near, outside (as built, after the large-draw bias check), and nonlinear treatment and outcome | Outcome and treatment shape | 1b, scope, procedure comparison | Label fits only; verdicts shared |
| Null versions of inside and outside | No effect | Table S1 | Label fits only |
| Baseline risk grid: family baseline in {0.03, 0.10, 0.30} against truth baseline in {0.03, 0.15, 0.30}, at the levels with 4.2%, 16% and 26% outside, n = 500 | Events; misspecified family | 1a (second factor), Figure 2 | Three extra verdict runs per dataset and level (family 0.15 exists) |
| Inside plus an unmeasured binary confounder | Identification | A certain disagreement; states the limit | Label fits only |
| Observed-treatment variant at the three middle levels, n = 500 | Treatment draw | 2 | 3 levels x R = 200 |
| Family extended by the true model (custom surface) for the near, outside and nonlinear truths, two levels, n = 500 | Family | 1b mechanism (skip if 1b is a sentence) | 3 truths x 2 levels x R = 100, generate only; about 1.25 times a standard run (assumed, one extra surface) |
| Two candidates (main-terms; main-terms plus MARS) on the nonlinear truth, one level, n = 500 | Library | Scope, selection | 60 runs at about 4x cost (assumed) |

I would drop the near model from the main display if, in the pilot, its labels equal the inside model's in every level, estimand and sample size. A rule stated now cannot be tuned later.

### 5.3 Displays and mockups

#### Figure 1 (main text; claims 1a, 1b)

*Verdict by overlap level against the estimability label.* Four rows (ATE, trimmed ATE, ATT, ATO) by two columns (n = 500; n = 2,000). x: overlap level, labelled by the share of patients outside 0.05 to 0.95. y: proportion of datasets on which the simulation judged the estimand feasible, 0 to 1, with 95% Wilson intervals over the R datasets. Four plotting symbols per level, one for each true outcome model, filled when the label says estimable and open when not. A grey line gives the feasible rate of the best tuned overlap diagnostic. Because the verdict never sees the outcome, the four symbols share one height at each level and differ only in fill, so a reader sees at once whether labels flip across truths where the verdict cannot. A good simulation puts filled symbols near 1 and open symbols near 0, with a transition at the level where the label changes; filled symbols low or open symbols high are misgrades. If the grey line coincides with the verdict's curve, the simulation adds nothing over a tuned diagnostic at that cell.

#### Table 1 (main text; claims 1a, 1b, 2)

*Agreement with estimability on the true outcome model.* Rows: estimand (4 groups) by verdict source (5 rows): design-stage simulation; same with observed treatment (n = 500 only); best single diagnostic at its tuned threshold; the same diagnostics at the package's grade cut-offs; best diagnostic with the threshold tuned at the other sample size. Columns: for each sample size, sensitivity and specificity for truths inside or near the family, and the same for truths outside (the outside model and the nonlinear scenario). Cells: estimate (Monte Carlo SE), effect scenarios only.

| Estimand and source | n = 500: inside or near, sens. | spec. | outside, sens. | spec. | n = 2,000: inside or near, sens. | spec. | outside, sens. | spec. |
|---|---|---|---|---|---|---|---|---|
| ATE: design-stage simulation | | | | | | | | |
| ATE: observed treatment | | | | | | | | |
| ATE: tuned diagnostic | | | | | | | | |
| ATE: grade cut-offs | | | | | | | | |
| ATE: threshold from the other n | | | | | | | | |
| (same five rows for the trimmed ATE, ATT and ATO) | | | | | | | | |

A reader should see the first row at least as high as the third in the first block, clearly above rows four and five if the claim is "without needing a tuned threshold", and the outside block either level with the first (claim 1b: no fall) or lower. A footnote gives the number of (level, truth) combinations in each label class and the share within two Monte Carlo SEs of a tolerance, with the sensitivity analysis that drops them.

#### Table 2 (main text; claim 3 and the proposal)

*Coverage and the price of moving.* Strategies as rows, overlap levels as columns, n = 2,000 (n = 500 in the supplement), inside truth with an effect.

| Strategy | Coverage of the reported estimand's true value, by level (5 columns) | Share reporting the ATE (5) | Share with no feasible estimand (5) |
|---|---|---|---|
| ATE, main-terms TMLE, no design stage | | | |
| ATE, weighting (main-terms) | | | |
| Prespecified order, feasibility from the design-stage simulation | | | |
| Prespecified order, feasibility from the tuned diagnostic | | | |
| Feasible estimand with the smallest P value | | | |
| Feasible estimand with the narrowest CI | | | |

The first two rows should fall below 0.90 as overlap worsens; the third should stay at or above 0.90 at every level while its ATE share falls, which is the price; rows five and six should show similar coverage, with the inflation visible only in Table S1. The table makes the protocol-estimand analysis and the package's procedure comparable on one scale, and the estimator comparison appears inside the first two rows.

#### Figure 2 (main text; claim 4)

*Calibration of the forecast half-width.* x: level, within estimand and sample size. y: ratio of realized to forecast half-width on a log scale, as a distribution of the datasets where the estimand was feasible, with a reference line at 1 and the target band marked. A second column repeats it for a family baseline of 0.10 and 0.30 against a truth of 0.15. A reader should see ratios at or below 1 with narrow spread when the family baseline matches the truth, and a shift of the right sign when it does not.

#### Supplement

| Display | Content | Claim |
|---|---|---|
| Table S1 | Rejection under the null and coverage for the prespecified order and both outcome-informed rules, by level, among datasets with two or more feasible estimands | 3 |
| Table S2 | DGP diagnostics (existing columns plus nonlinearity and the main-terms estimator's large-draw bias), the pilot sizing record (gap, SD, |z| and `n_sim_safe` per contrast), and the share of verdicts resting on unresolved cells | All |
| Table S3 | Bias, coverage, model SE over empirical SE and half-width by library (main-terms; main-terms plus MARS; default super learner on a few hundred datasets) for in-family and strongly nonlinear truths | Scope; learner item |
| Table S4 | Agreement and forecast calibration by family baseline against truth baseline | 1a, 4 |
| Table S5 | Bias, empirical SE, root mean squared error and coverage by scenario, estimand, estimator and sample size, with Monte Carlo SEs, including weighting, matching (ATT only) and the semi-oracle TMLE | The proposal |

Cut from the plan: the pooled Table 1; the super learner subset in the design-stage simulation; the proposal's Figure 3 and its Table 2 as drafted; the observed-treatment variant at n = 2,000; the correlation across scenarios in Figure 2; baseline variants at every level and sample size; matching in the main text. The receiver operating characteristic figure (`PLAN.md:59-61`) can stay as an optional supplementary figure because it is derived from stored rows and costs no compute.

### 5.4 Sizing and order of work

Order: run the analysis of the existing pilot first; read the labels and the paired contrasts; fix the decision rules and the list of scenarios in `PLAN.md`; implement the code changes; run the large-draw bias checks; then the full run in the order generate verdicts, label fits, baseline grid, observed-treatment variant, library comparison. Choose R as the maximum `n_sim_safe` over the kept contrasts, capped at what the budget allows (about 4 days of wall time at R = 300 for the generate verdicts across both sample sizes, section 3.5). The extra stages add: baseline grid, 3 families x 3 levels x R = 200 x 424 seconds = 212 worker-hours (43 hours wall); observed-treatment variant, 3 levels x 200 x 530 seconds = 88 worker-hours (18 hours wall); library comparison and design-stage selection, about 60 runs at an assumed 4 times the generate cost at one level, about 28 worker-hours (6 hours wall).

## 6. Moving the 0.4.1 pilot to 0.5.0

### 6.1 What changed between the pilot's install and 0.5.0

The pilot was installed from commit 6770bbd (`pilot-report.md:5`). `git diff 6770bbd..HEAD -- R` in `cleanTMLE` shows 12 changed files. The numerical code the simulation calls is unchanged: `R/fit_candidate.R`, `R/learners.R`, `R/utils.R` and `R/estimate_effect.R` have no differences, and `R/simulate_design.R` and `R/assess_design.R` each gain one line, `.need_dictionary(lock)` (`simulate_design.R:284`, `assess_design.R:93`). What changed is the interface: `analysis_plan()` takes a required `time_zero` (`R/plan.R:168`, `:217`), and `create_analysis_lock()` takes a required `dictionary` (`R/lock.R:74-85`, fields in `R/dictionary.R:9-11`). The main branch (2e96f98, ahead of origin by 2) is the 0.5.0 release; the submodule's working tree is on `vignette-limits` with an uncommitted vignette edit, so the install should be exported with `git archive main`, as `simulation/setup.R:5-8` already says.

### 6.2 Code changes

| File and line | Change |
|---|---|
| `simulation/R/map.R:18-24` (`make_plan`) | Add `time_zero = "treatment assignment (simulated)"` to `analysis_plan()`. `key.R:9` calls the same function. |
| `simulation/R/map.R:104-105` (`run_map`) | Pass `dictionary = sim_dictionary()` to `create_analysis_lock()`; define `sim_dictionary()` once (id, treatment fixed at time zero, four covariates fixed before; the template is `cr_example$dictionary`, `lock.R:81-84`). |
| `simulation/pilot.R:15-16` | Change the `stopifnot` version from "0.4.1" to "0.5.0". Run step 3a first, before the 0.5.0 install replaces the library, or the assertion stops it. |
| `simulation/R/map.R:26-62` | Keep the local copy of `.one_rep` (unchanged in 0.5.0) but add a check that, with the treatment draw swapped back, it reproduces `cleanTMLE:::.one_rep` on one repetition, so drift is caught. Update the comment. |
| `simulation/R/map.R:109-111` | Return `s$results` as well as the verdict and metrics, so the raw per-repetition rows are saved. |
| `simulation/R/key.R:24-31`, `simulation/R/summarize.R:186-211` | Implement both outcome-informed rules (smallest P value; narrowest CI among estimands the verdict calls feasible). |
| `simulation/R/map.R:21`, `simulation/R/dgp.R:13` | Make the family baseline risk and the truth's baseline risk arguments, for the baseline grid. |
| `simulation/R/summarize.R:68-78` | Restrict Table 1 to effect scenarios. |
| `simulation/R/run.R:24-29` | Add R, tmle, WeightIt and MatchIt versions to the cache key; `tmle::tmle` does 88% of the work (`pilot-report.md:67`), so its version belongs in the key. |
| New `simulation/R/adapter.R` | Collect the 10 internals now called through `cleanTMLE:::` (`.design_matrix`, `.failed_rows`, `.fit_and_target`, `.fit_g`, `.make_folds`, `.sim_metrics`, `.sim_verdict`, `.surfaces`, `.true_values`, `.with_seed`; `map.R:36-79`, `key.R:10-16`) behind one file with a test that each exists with the expected arguments. If the package is slimmed further, only this file changes. |
| `simulation/PLAN.md:3`, `simulation/BUILD-BRIEF.md:9`, `writing/brief.md` (last "What not to say" bullet) | Replace "0.4.1" with "0.5.0". |

### 6.3 What the cached results are good for

| Item | Valid on 0.5.0? | Action |
|---|---|---|
| `simulation/pilot/calibration.rds`, `calibration_grid.csv`, `dgp_diagnostics.csv` | Yes. They come from `dgp.R`, which calls no cleanTMLE code. The `.rds` records version 0.4.1 (`pilot.R:42-45`) but does not depend on it. | None |
| `simulation/cache/key` (6,000 cells, n = 500) | Numerically yes in expectation, since `.fit_and_target` and the fold code are unchanged. Stale by key: the version and the code hash differ. | Recompute under 0.5.0 (about 11 minutes) and compare with the old rows; no reason to carry them over |
| `simulation/cache/map_generate`, `map_sample` (500 cells each, n = 500) | Same reasoning; costly (about a day of wall time to recreate). They hold verdicts and metrics but not the 200 repetitions per run. | Use for the pilot analysis as they stand; verify three cells under 0.5.0 (bit-for-bit apart from the timing columns); the full run recomputes them, now with repetition rows saved |
| `simulation/pilot/smoke.rds`, `parallel_test_w6.rds`, `cache/_old_smoke_pre_revision`, `cache/_partest` | No. Pre-revision code or test scaffolding. | Andrew can delete them |
| Timings in `pilot-report.md` | Planning only; they vary about two-fold with machine load (`pilot-report.md:7`) | Re-time on the free machine |

The pilot analysis does not need the cache keys: `collect()` reads every file in a stage folder (`run.R:87-94`). Run step 3a under the existing 0.4.1 library, then install 0.5.0. The check that 0.5.0 reproduces 0.4.1 (recompute the label fits, compare, verify three design-stage cells) takes under an hour and is what lets the pilot's sizing carry over; if any cell differs, the pilot is a record of 0.4.1 only and the full run decides.

## 7. Top recommendations in priority order

| Priority | Recommendation | Effort | Risk | Decision for Andrew |
|---|---|---|---|---|
| 1 | Run the pilot analysis (`Rscript simulation/pilot.R 3a`) as soon as the Rescue.Co job has finished, on the 0.4.1 library, before any new compute. Read the labels across truths, the verdict against the best diagnostic, the paired contrasts and the sizes. | Minutes of machine time, 1 hour to read | None: it only reads the cache | Approve the first R job after the design stage ends (one R job at a time) |
| 2 | Rescope claim 1: write the pass rule, add baseline risk as a second factor, add the grade cut-off and transported-threshold comparators, state 1b with an equivalence margin | 1 to 2 days of code; about 43 hours of machine time for the baseline grid | The verdict may not beat a tuned diagnostic; the honest claim is then that it matches one without tuning | Choose the equivalence margin and the pass rule; accept that "more accurately than overlap diagnostics" in `writing/brief.md` may become "as accurately, without a tuned threshold" |
| 3 | Fix claim 3 and merge it with the proposal as a procedure comparison: both outcome-informed rules, a table that judges each strategy by the estimand it reports and shows the share reporting the ATE | 1 day | Both rules may show nothing at the levels with one feasible estimand | Approve reporting both rules and the share-reporting-ATE column |
| 4 | Answer the proposal's three decisions as in section 4.2: reword the claim, move the candidate comparison to the supplement, replace Table 2 and Figure 3; add the outcome-informed arm, the flexible-weighting comparator and the large-draw bias check | 1 day of design; compute small | The nonlinear scenario's bias may be below 0.02 and need amplifying (disclosed, per the skill) | Decisions 1 to 3 |
| 5 | Close question 4 of the learner-defaults item with the minimal library-by-truth table and the selection check on the nonlinear truth; write the truth as a custom surface | 1 day; about 28 worker-hours for the selection check, plus the label fits | The default family may select main-terms on a truth that needs the flexible candidate; that is a result for the Discussion, not a defect to hide | Which strongly nonlinear truth to commit, before any run |
| 6 | Move the code to 0.5.0 (section 6.2): two call sites, repetition rows saved, adapter file, drift test, cache key; verify numerics against the cache | 0.5 day | Low; numerics unchanged | Whether the simulation pins to the 0.5.0 release tag when it is cut |
| 7 | Remove the super learner subset from the design-stage simulation (`methods.qmd:158`) and say so; use the super learner only for the label fits | 1 sentence | The paper loses a demonstration it could not afford | Approve the rewording in section 8 |
| 8 | Freeze the five overlap levels, the truths, the displays and the decision rules in `PLAN.md` before the full run; record any later change | 2 hours | None | Sign-off on the frozen list |
| 9 | Size the full run from the pilot: M = 5,000 label fits (update `methods.qmd:148`, which says 1,000); R from the largest kept requirement; observed-treatment variant at n = 500 only | Compute about 5 days wall at R = 100 and about 8 days at R = 300 for the whole design | Budget | R and the compute window |

## 8. Proposed text changes (before and after)

None is applied. The paper text avoids internal codes.

#### `manuscript/methods.qmd:157`

Before: "To assess choosing the estimand after the outcome, we compared, on the same datasets with the outcome, the estimate of the primary estimand from the prespecified order with the estimate of whichever feasible estimand had the narrowest CI, each against the true value of the estimand it reports, on the datasets for which the simulation judged at least one estimand feasible."

After: "To assess choosing the estimand after the outcome, we compared, on the same datasets with the outcome, the primary estimand from the prespecified order with two ways of choosing among the feasible estimands once the outcome was seen: the estimand with the smallest P value, which selects on the estimate, and the estimand with the narrowest CI, which selects on a quantity that is nearly independent of it. Each was evaluated against the true value of the estimand it reports, and the rejection rate was computed under the null models, among datasets in which at least two estimands were judged feasible."

#### `manuscript/methods.qmd:158`

Before: "The candidates were the main-terms logistic regression for the outcome and the propensity score, with the propensity score truncated at 0.01 and five-fold cross-fitting, across the whole grid; the package's default super learner was used in the scenarios with the middle level of overlap at n = 2,000."

After: "The candidate in the design-stage simulations was main-terms logistic regression for the outcome and the propensity score, with the propensity score truncated at 0.01 and five-fold cross-fitting. Running the package's default super learner inside every repetition of every simulation was not feasible (one repetition took [R: seconds per repetition] seconds on the registry of the application), so the super learner and a reduced flexible library were evaluated on the true outcome models only, on [R: number] datasets per scenario."

#### `manuscript/methods.qmd:162`

Before: "with Monte Carlo standard errors that treat datasets as clustered within scenarios"

After: "with Monte Carlo standard errors computed from the spread, across dataset indices, of the dataset-level agreement, because one simulation per dataset serves every true outcome model at its overlap level"

#### `writing/brief.md`, key claim 1 (conditional on the pilot analysis)

Before: "agrees with whether the estimand is estimable on the true outcome model more often than overlap diagnostics at their best thresholds, when the true outcome model lies inside or near the simulated family; the agreement falls when it lies outside."

After, if the verdict does not beat a tuned diagnostic: "agrees with whether the estimand is estimable on the true outcome model about as often as the best overlap diagnostic, whose threshold was tuned on the true labels, and more often than diagnostics at conventional or transported thresholds; the agreement [does not fall | falls by an amount the study reports] when the true outcome model lies outside the simulated family."

## 9. Limits of this audit and check outcome

Three kinds of statement here are inferences and not readings, and I have marked them where they occur. The first is the argument from double robustness about which labels can differ across truths, to be settled by the labels themselves and by the large-draw bias check. The second is the independence of interval width from the estimate, to be settled by the claim 3 run. The third is the cost of a super learner repetition at four covariates and n = 500, which I scaled from the Rescue.Co timings in `rescueCo/logs/pipeline.log:1204-1205` and did not measure. Every q and g in section 3.4 and the standard deviations of the paired differences are assumed; the measured inputs are the 40-dataset label check, the timing files, the smoke runs and the cache listing. I did not open any cache entry, `.rds` file, Rescue.Co outcome, or the unblinded analysis, and nothing in the recommendations rests on them.

The writing check was run from the project root as `python "C:/Users/andre/OneDrive/Documents/scientific-writing-skills/checks/check.py" "audit/2026-10-22/03-simulation.md" --min-level warn`. Every finding that is not expected was fixed. What remains is of four kinds. Most of the errors flag code identifiers, file paths and line references in backticks, which an internal audit needs and which are expected. One error is the unresolved citation placeholder in section 1.2, which the audit rules require when a claim needs a source that is not in `writing/references.bib`. A few errors flag numbers that are not in the registry: counts read from the cache listing, the plan's number of label datasets, my recommended number of label datasets, and one arithmetic result for compute, each derived where it appears. One warning comes from the project's writing configuration (its overlap source file is not found), not from this file. The file has no em dashes, and the proposed paper text in section 8 contains no internal codes.

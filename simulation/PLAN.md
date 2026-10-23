# Simulation plan: estimability paper (Statistics in Medicine)

Backwards design (stats-research:design-and-report-simulations). Stages 1 to 3 agreed with Andrew on 2026-10-07; stage 4 is a draft that the pilot (stage 6) must confirm before anything is frozen. The simulation code lives in `simulation/` and calls the cleanTMLE package (0.4.1 or later); it does not modify the package.

## 1. Claims

| # | Claim | Evidence | Priority |
|---|---|---|---|
| 1 | The outcome-free support map's verdict for each estimand agrees with whether that estimand is in fact estimable more often than overlap diagnostics that do not simulate, when the true outcome model lies inside or near the declared family of outcome surfaces; the agreement falls when the true model lies outside it. | simulation | main |
| 2 | Drawing new treatments from the estimated propensity score (generate-treatment) grades estimability correctly; keeping each patient's observed treatment (sample-treatment) misgrades it. | simulation (one more comparator in claim 1's design) | main |
| 3 | Choosing the estimand from the ladder before outcomes keeps confidence-interval coverage nominal for the chosen estimand; choosing it after seeing outcomes, for example the estimand with the smallest p-value, undercovers and inflates the type I error. | simulation (and a short argument) | supplement |
| 4 | The pre-outcome expected confidence-interval half-width predicts the realized half-width. | simulation (same runs as claim 1) | secondary display |

The boundary in claim 1 (outside the family) answers the argument that a plasmode simulation can refute but cannot certify.

## 2. Simulation goals

### Goal 1 (claims 1 and 2)

- **Unit:** a scenario (a covariate and treatment model, a true outcome model, a sample size n) and an estimand (ATE, trimmed ATE with band 0.05 to 0.95, ATT, ATO).
- **Answer key (agreed):** an estimand is estimable in a scenario if at least one candidate meets the plan's tolerances on the true outcome model: absolute bias at most 0.02 on the risk difference scale and 95% CI coverage at least 0.90, over M datasets drawn from the true model with the real (simulated) outcome. This mirrors the map's own feasibility rule, applied to the truth instead of the declared surfaces. Population truth for each estimand comes from one very large draw.
- **Verdicts compared:** (a) the support map as implemented (generate-treatment), run outcome-free on each of R datasets per scenario; (b) the same map with sample-treatment; (c) overlap diagnostics computed on the same datasets, each with its best threshold chosen against the answer key (agreed): share of patients outside the band, largest weight for the estimand, effective sample size as a fraction of n for the estimand, and the share of patients with a preference score in 0.3 to 0.7.
- **Measures:** sensitivity and specificity of each verdict for "estimable", with Monte Carlo SEs, by the true outcome model's relation to the family (inside, near, outside); and the ROC curve of each continuous diagnostic with the map's operating point.
- **Sharpened:** across the five overlap levels crossed with the true outcome models at n = 500 and 2,000, the generate-treatment support map's balanced accuracy for estimability is higher than each comparator's, or no more than 0.05 below the best tuned diagnostic's, for each informative estimand when the true outcome model lies inside or near the family; and its agreement falls when the model lies outside it. The pass rules are in "Decision rules, fixed before the contrast pilot is summarized".
- **Oracle:** the answer key is itself the oracle (it uses the true outcome model). A semi-oracle map that uses the true propensity score instead of the estimated one separates errors in the map's simulated world from propensity-score estimation error; include it if the pilot shows the map's errors concentrate at small n.

### Goal 2 (claim 3)

- On the same R datasets per scenario, with the real outcome: the prespecified estimate (the ladder's primary from that dataset's map) versus two outcome-informed choices among the estimands the dataset's map judged feasible, the one with the smallest p-value and the one with the narrowest CI. Measures: coverage of the 95% CI for the chosen estimand's true value; under a null effect (separate null versions of each scenario), the rejection rate at 0.05; both on the datasets where at least two estimands are feasible.
- **Sharpened:** choosing by the smallest p-value rejects a true null more often than the prespecified choice, while choosing by the narrowest CI changes neither coverage nor the rejection rate by more than 0.03.

### Goal 3 (claim 4)

- Per scenario, estimand and dataset: the map's expected half-width for the selected candidate versus the realized half-width of that candidate's CI on the real outcome. Measure: their ratio and correlation across scenarios; calibration plot.

## 3. Mockups (empty cells; captions written first)

### Table 1 (main text; claims 1 and 2)

*Caption:* Agreement between each pre-outcome verdict and whether the estimand was in fact estimable (some candidate met |bias| <= 0.02 and coverage >= 0.90 on the true outcome model), by the true outcome model's relation to the declared outcome surfaces. Pooled over estimands, overlap levels and sample sizes; Monte Carlo SEs in parentheses. Overlap diagnostics use their best threshold, chosen against the answer key.

| Verdict | Inside: sensitivity | Inside: specificity | Near: sensitivity | Near: specificity | Outside: sensitivity | Outside: specificity |
|---|---|---|---|---|---|---|
| Support map, generate-treatment | | | | | | |
| Support map, sample-treatment | | | | | | |
| Share outside the band | | | | | | |
| Largest weight | | | | | | |
| Effective sample size fraction | | | | | | |
| Preference-score equipoise | | | | | | |

*Checks:* if the map wins inside and near and loses its edge outside, claims 1 and 2 hold with their boundary; if a diagnostic matches the map everywhere, claim 1 fails and the table shows it.

### Figure 1 (main text; claim 1)

*Caption:* How often the support map judged each estimand feasible, by overlap and true outcome model, against the answer key. Each point is a scenario at n = 2,000; filled points are estimable by the answer key, open points are not. Bars are Monte Carlo 95% intervals.

Facets: estimand (ATE, trimmed ATE, ATT, ATO). x: overlap (good, moderate, severe). y: share of map runs with verdict feasible (0 to 1). Color: inside, near, outside. A good map puts filled points near 1 and open points near 0.

### Figure S1 (supplement; claim 1 fairness)

*Caption:* ROC curves of each overlap diagnostic for estimability, with the support map's single operating point, by true outcome model.

### Table S1 (supplement; claim 3)

*Caption:* Coverage of the 95% CI for the chosen estimand, and rejection rate under a null effect, when the estimand is chosen from the ladder before outcomes versus after seeing outcomes (smallest p-value or narrowest CI), on the datasets where at least two estimands are feasible. By overlap; n = 2,000; Monte Carlo SEs in parentheses.

| Choice | Good: coverage | Good: rejection (null) | Moderate: coverage | Moderate: rejection | Severe: coverage | Severe: rejection |
|---|---|---|---|---|---|---|
| Prespecified (ladder) | | | | | | |
| Outcome-informed (smallest p) | | | | | | |
| Outcome-informed (narrowest CI) | | | | | | |

### Figure 2 (main text or supplement; claim 4)

*Caption:* Expected CI half-width from the outcome-free simulation versus the realized half-width on the outcome, for the selected candidate. One point per scenario and estimand; the line is equality.

### Table S2 (supplement)

DGP diagnostics: true value of each estimand, share of patients outside the band, c-statistic of the true propensity score, outcome prevalence, variance explained by the covariates.

**Main text:** Table 1, Figure 1, Figure 2. **Supplement:** Figure S1, Tables S1 and S2, the super learner subset.

## 4. DGPs and settings (draft; the pilot decides)

Covariates: X1 ~ Uniform(-2, 2) (the strong confounder, like severity), X2 ~ Bernoulli(0.5), X3 ~ Uniform(-1, 1), X4 ~ Uniform(-1, 1) (unrelated to treatment). Bounded so that overlap is controlled by the treatment model alone.

Treatment: logit pi(X) = -0.2 + gamma X1 + 0.3 X2 + 0.3 X3, with gamma set in the pilot so that the share of patients with pi outside 0.05 to 0.95 is under 5% (good), about 10 to 20% (moderate) and over 30% (severe). Reason: overlap is the factor the map is meant to read.

Outcome (binary, baseline risk about 0.15; b0 calibrated):

| Model | Relation to the declared family | logit mu(A, X) | Reason |
|---|---|---|---|
| Inside | inside | b0 + c (0.5 z(logit pi) + 0.5 z(X beta)) + A (0.5 + 0.5 z(logit pi)) | the family's own form, with the true pi |
| Near | near | b0 + 0.8 X1 + 0.6 X1^2 + 0.4 X3 + A (0.5 + 0.6 X1) | nonlinear prognosis and effect heterogeneity in the confounder; close to the family's nonlinear surface |
| Outside | outside | b0 + 0.8 X1 + 0.4 X3 + A (0.5 + 1.0 X4 + 0.8 1{X1 > 1.5}) | heterogeneity in a covariate unrelated to treatment and a threshold in the region of poor overlap, which the family's smooth heterogeneity in pi cannot mimic |

Null versions (claim 3): the same models with the A terms removed.

- **Sample sizes:** 500 and 2,000 (a small and a typical registry cohort).
- **Estimands:** all four; the main text pools them, the supplement reports each.
- **Methods:** cleanTMLE support map (generate-treatment, as implemented); sample-treatment variant (paper code); overlap diagnostics from assess_design(); answer key from fit_candidate() on the true model.
- **Learners:** main grid: main-terms GLM for the propensity score and outcome, truncation 0.01 and 0.05, K = 5 cross-fitting (agreed, given compute). Subset: the package's default super learner (one learner per family, set up per the supervised-learning skill) on the moderate-overlap scenarios at n = 2,000.
- **Map settings:** the package's default surface family (four surfaces); reps = 200 with max_reps = reps, so every map run has the same budget.
- **Repetitions:** R map runs per scenario and M answer-key datasets per scenario, set by the pilot; the R map datasets are the first R of the M, so the map, the key, claim 3 and claim 4 read the same datasets.

### Revisions after pilot steps 1 and 2 (2026-10-07, approved by Andrew)

Evidence: simulation/pilot/pilot-report.md.

- **Shared map datasets.** The map never sees the outcome, so its verdict on a dataset does not depend on the true outcome model. Covariates and treatment are seeded from (overlap, n, dataset index); each true outcome model draws its outcome on those same covariates and treatment. Map runs therefore depend only on (overlap, n, dataset), and every comparison across outcome models is paired.
- **Baseline risk of the simulated family set to 0.15**, the true control-arm risk, so the inside model is inside the family in both form and level; the near and outside models differ from the family only in form.
- **Overlap levels near the boundary of estimability.** Two levels are added between the existing ones, at gamma 1.5 (about 4% outside the band) and 2.0 (about 26%), because with good, moderate and severe overlap alone the ATO was feasible everywhere and the ATE and ATT infeasible at moderate and severe overlap, cases any diagnostic gets right. Table 1 is reported per estimand, not pooled. This choice was made to make the comparison informative, after seeing the pilot, and the paper says so.
- **Compute.** Main grid with one truncation level (0.01): the claims concern feasibility, not truncation choice, and the pilot's trimmed-ATE estimates did not differ between 0.01 and 0.05. Parallel workers (6 to 8) after a stability test. reps stays 200. The number of map runs per scenario R is set by the contrast pilot.

### Revisions after the Methods review (2026-10-07, approved by Andrew)

- **Misspecified baseline risk.** Superseded on 2026-10-23 by the baseline-risk factor in "Decision rules, fixed before the contrast pilot is summarized", which keeps these two runs, at the three middle overlap levels and n = 500, and reports them in their own table.
- **Answer key.** M = 1,000 key datasets per scenario. Report the share of scenario-estimand combinations whose bias or coverage lies within 2 Monte Carlo SEs of a tolerance, and a sensitivity analysis of Table 1 that excludes them.
- **To implement after the contrast pilot finishes** (the cache keys include the code's text, so the code must not change during a run): the narrowest-CI rule for claim 3, M = 1,000, the two baseline variants, and the ambiguity sensitivity analysis.
- **Claim 3 comparator.** Superseded on 2026-10-23: both the narrowest-CI rule and the smallest-P rule are reported, each among the feasible estimands, as set out in "Decision rules, fixed before the contrast pilot is summarized".
- **Coverage tolerance 0.90** is kept and described in the paper as somewhat arbitrary.

### Decision rules, fixed before the contrast pilot is summarized (2026-10-23)

**What had been seen when these rules were written.** Pilot steps 1 and 2 (calibration, DGP diagnostics, and smoke runs on two datasets in each of two scenarios), the 40-dataset check of the answer key in three scenarios (`pilot/key_check_output.txt`), and the timing runs. The contrast pilot's verdicts and label fits are cached, but its summary (`pilot.R` step 3a) has not been run, and no agreement, contrast, claim 3 or claim 4 result from it has been seen. Every rule below is fixed before step 3a is run. A rule changed after that is recorded here with its date and with the outputs seen before the change.

**Claim 1a: the verdict against overlap diagnostics.**

- *Measure.* Balanced accuracy, the mean of sensitivity and specificity, for each estimand and sample size, over the five overlap levels and the truths inside or near the family, computed per dataset index and paired across verdict sources on the same datasets. Its Monte Carlo SE is the SD of the per-index difference over the square root of R.
- *Comparators.* Three, in this order. First, the best tuned overlap diagnostic, whose threshold is chosen against the answer key in the same cell; it is generous to the comparator, since no analyst knows the answer key. Second, the same diagnostics with the threshold tuned at the other sample size, which an analyst could carry over from a published simulation. Third, the package's overlap grade (good below 5% outside the band, moderate below 20%), read as calling the ATE and ATT estimable when the grade is good and the trimmed ATE and ATO estimable when it is good or moderate.
- *Pass rule, per estimand and sample size.* The verdict is **better** than a comparator when the lower bound of the 95% interval for the difference in balanced accuracy is above 0. It is **no worse** than the best tuned diagnostic when that lower bound is above -0.05. Otherwise the comparison is **unresolved**, and is reported as such.
- *Claim 1a holds* in its stronger form if the verdict is better than all three comparators for every informative estimand and sample size, and in its weaker form if it is no worse than the best tuned diagnostic and better than the second and third comparators in every one. The weaker form is stated as "as accurate as a diagnostic tuned on the true outcome model, without a threshold that no analyst could tune". An estimand whose label does not vary across levels and truths at a sample size tests nothing; it is reported with its agreement rate and named as uninformative, and it does not count against the claim.
- *Margin.* 0.05 in balanced accuracy, one estimand in twenty misgraded on average over the two label classes. Differences smaller than this are below the misgrading that the verdict's own noise near a tolerance produces at 200 repetitions per verdict (audit 03, section 3.3, binomial arithmetic).

**Claim 1b: agreement outside the family.** The measure is the verdict's balanced accuracy for the truths outside the family (the outside model, and the nonlinear-treatment truth if it is added) minus that for the truths inside or near it, paired by dataset index. Agreement **falls** when the upper bound of the 95% interval for this difference is below 0, and is **not shown to fall** when the interval lies within ±0.05. Otherwise the comparison is unresolved. Because the verdict never sees the outcome, this difference comes only from labels that differ across truths, and the paper says so.

**Claim 3: the outcome-informed comparators.** Both rules are fixed now and both are reported, so that neither can be chosen after the other is seen. Each chooses among the estimands that the dataset's own verdict judged feasible, the set the prespecified order chooses from: the estimand with the smallest P value, which selects on the estimate, and the estimand with the narrowest CI, which selects on the width. The measures are coverage of the 95% CI for the chosen estimand's true value (effect scenarios) and the rejection rate at 0.05 (null scenarios), paired by dataset with the prespecified choice. They are computed on the datasets where the verdict judged at least two estimands feasible, since elsewhere the choices coincide, and the share of datasets with two, three or four feasible estimands is reported.

- Selection on the estimate **inflates the type I error** when the lower bound of the 95% interval for its rejection rate minus the prespecified choice's is above 0.
- Selection on the width is **reported as not changing** coverage or the rejection rate when the 95% intervals of both differences lie within ±0.03. Otherwise its estimates and intervals are reported without a claim.
- The claim is scoped to what these show. No claim is made that the gap grows with the number of feasible estimands unless a trend over two, three and four feasible estimands has a 95% interval that excludes zero.

**The baseline-risk factor.** Overlap is the only factor in the main grid that changes whether an estimand is estimable, and every comparator diagnostic reads overlap, so the design gives the verdict no condition that a summary of the propensity score cannot see. The outcome's baseline risk is such a condition, because it sets the number of events.

- *Truths.* The inside and outside models with control-arm risk 0.03, 0.15 (the main grid) and 0.30, b0 recalibrated for each on the population draw, at the three middle overlap levels (gamma 1.5, 1.75 and 2.0) and n = 500. Covariates and treatment are shared with the main grid, so these add label fits only.
- *Family.* The verdict is run with the family's baseline risk at 0.03, 0.10, 0.15 and 0.30 on the same datasets. The run at 0.15 is the main grid's. The runs at 0.10 and 0.30 are the misspecified variants approved on 7 October, here at these three levels and n = 500. The run at 0.03 is new. Each family setting is crossed with each truth, so a matched family (an analyst who knows the baseline risk) and a misspecified one are both shown.
- *A fair comparator.* An analyst who looks only at the propensity score cannot see the baseline risk, so against those diagnostics the verdict wins on this factor by construction. The comparator for this factor is therefore the best overlap diagnostic combined with an events rule an analyst could apply before the outcome: expected events in the arm with fewer, from the baseline risk the analyst gave the family, at least 10 per covariate column (`peduzziSimulationStudyNumber1996`). The package already reports these expected events before the outcome (`simulate_design()`, `expected_events`).
- *Reporting.* The factor is reported in its own table, not pooled into Table 1. The claim it can support is what the simulation adds over overlap and an events rule judged separately, using claim 1a's rule and margin.

**Dropping a truth.** A truth leaves the main display, and moves to the supplement, only if its estimability labels in the full run equal those of a truth kept in the display in every overlap level, estimand and sample size. Truths are kept in the order inside, outside, nonlinear treatment, near, so the later one moves. No truth is moved or redefined for any other reason, in particular not because the verdict agrees with it poorly.

**Contrasts, levels and sizing.** Step 3a is used to size the full run: R and M are set by the largest requirement among the contrasts above, at k = 3 for each "better" or "falls" comparison and by the equivalence bound (`pilot_equivalence()`) for each margin. A contrast that the affordable R cannot resolve stays in its table, marked unresolved, and the claim it serves is rescoped to match. It is not removed. No overlap level, truth or comparator is added or removed after step 3a is read, except by the rule for dropping a truth.

## 5. Build (to do after the pilot plan is agreed)

One function per stage, each runnable for any subset of scenarios, datasets and methods; each dataset seeded from (scenario, n, dataset index); per-dataset results cached under `simulation/cache/` (git-ignored).

## 6. Pilot (to do)

Checks: runs at 2 datasets; one large draw per scenario for Table S2's diagnostics and the gamma calibration; then about 100 map runs and 200 key datasets on the hardest contrasts (map versus best diagnostic, inside versus outside) to size R and M with pilot_check.

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
- **Sharpened:** across scenarios crossing overlap (good, moderate, severe) with the true outcome model (inside, near, outside the declared family) at n = 500 and 2,000, the generate-treatment support map has higher sensitivity and specificity for estimability than sample-treatment and than every best-threshold overlap diagnostic when the truth is inside or near the family, and its agreement falls when the truth is outside it.
- **Oracle:** the answer key is itself the oracle (it uses the true outcome model). A semi-oracle map that uses the true propensity score instead of the estimated one separates errors in the map's simulated world from propensity-score estimation error; include it if the pilot shows the map's errors concentrate at small n.

### Goal 2 (claim 3)

- On the same R datasets per scenario, with the real outcome: the prespecified estimate (the ladder's primary from that dataset's map) versus an outcome-informed choice (among the estimands, the one with the smallest p-value). Measures: coverage of the 95% CI for the chosen estimand's true value; under a null effect (separate null versions of each scenario), the rejection rate at 0.05.
- **Sharpened:** the prespecified choice covers at the nominal rate (within Monte Carlo error) where the map says the estimand is feasible, while the outcome-informed choice undercovers and rejects a true null more often than 5%, with the gap growing as more estimands are available.

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

*Caption:* Coverage of the 95% CI for the chosen estimand, and rejection rate under a null effect, when the estimand is chosen from the ladder before outcomes versus after seeing outcomes (smallest p-value). By overlap; n = 2,000; Monte Carlo SEs in parentheses.

| Choice | Good: coverage | Good: rejection (null) | Moderate: coverage | Moderate: rejection | Severe: coverage | Severe: rejection |
|---|---|---|---|---|---|---|
| Prespecified (ladder) | | | | | | |
| Outcome-informed (smallest p) | | | | | | |

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

## 5. Build (to do after the pilot plan is agreed)

One function per stage, each runnable for any subset of scenarios, datasets and methods; each dataset seeded from (scenario, n, dataset index); per-dataset results cached under `simulation/cache/` (git-ignored).

## 6. Pilot (to do)

Checks: runs at 2 datasets; one large draw per scenario for Table S2's diagnostics and the gamma calibration; then about 100 map runs and 200 key datasets on the hardest contrasts (map versus best diagnostic, inside versus outside) to size R and M with pilot_check.

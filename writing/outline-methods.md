# Outline: Methods (Statistics in Medicine estimability paper)

**Message (from the brief, proposed):** before any outcome is seen, an outcome-blind plasmode simulation can tell which causal estimands a design can deliver, and choosing the estimand at that stage keeps inference valid.

Digest items: the registry is empty. Case-study counts and simulation settings that are outputs (calibrated overlap levels, repetitions) are tagged [R: ...] and must come from the registry once the simulation and the Rescue.Co re-run exist. Design constants the authors chose (tolerances, trim band, truncation levels) are part of the method, not results, and are written as declared values.

Citations needed are tagged [cite: ...]; sci-cite resolves them against the project's Zotero collection.

## 2.1 Data structure, estimands and identification

1. Notation: O = (W, A, Y) with covariates W, binary treatment A, binary outcome Y; n independent observations; propensity score g(W) = P(A = 1 | W); outcome regression Q(a, W). A notation table goes in the supplement (audience: methods-stats). [cite: TMLE notation, van der Laan and Rose]
2. The four estimands as risk differences, written formally: ATE; trimmed ATE over {W: g(W) in [0.05, 0.95]}; ATT; ATO with weights g(1 - g). Each describes a different population; in ICH E9(R1) terms they differ in the population attribute. [cite: Crump 2009; Li, Morgan and Zaslavsky 2018; ICH E9(R1)]
3. Identifying assumptions in words: consistency, no unmeasured confounding, and positivity, with the point that positivity is the assumption the estimand choice can change: the ATE needs g bounded away from 0 and 1 everywhere, the ATT only away from 1, the ATO not at all. Practical positivity violations and why truncation does not remove extrapolation. [cite: Petersen 2012]

## 2.2 The staged workflow

4. Overview in one paragraph: the plan is declared and fingerprinted with the outcome-free design data; the design is assessed and simulated without the outcome; negative controls are checked; the SAP (the design report, with the decision) is reviewed and approved; only then is the outcome admitted and the selected estimator run. Why: separating design from outcome removes analytic flexibility, and the design stage answers estimability. [cite: Muntner 2024 staging and clean room; Rubin 2008 objective causal inference]
5. What the lock does and does not guarantee: tamper-evident (a hash of the plan and design data; any later change is detectable), not tamper-proof (it cannot stop an analyst from analysing the outcome elsewhere).

## 2.3 The outcome-blind plasmode simulation

6. Purpose and construction: keep the real covariates; draw new treatments from the estimated propensity score (generate-treatment) rather than keeping each patient's observed treatment, because keeping them builds a positivity violation into the simulation [cite: Shaw 2025]; draw outcomes from a declared family of outcome models fixed in the plan. [cite: plasmode simulation origins, Franklin 2014]
7. The family of simulated outcome models: logit Q_s(a, W) = b0 + gamma m(W) + a (tau + delta h(g(W))), with m a prognostic index (linear or nonlinear), tau the effect on the log-odds scale, delta effect heterogeneity in the propensity score (the variation that separates the ATE, ATT and ATO), b0 calibrated to a declared baseline risk; the grid crosses declared values. Why heterogeneity in the propensity score: it is what makes extrapolation costly. Scope: estimability is supported only over this family.
8. True values computed exactly from each simulated outcome model and the propensity score, so the simulation needs no outcome and has no Monte Carlo error in its truths.
9. Candidates and the shared fitter: each candidate is an outcome library and a propensity-score truncation level; the simulation and the final analysis use the same fitter and the same candidate, with K-fold cross-fitting. [cite: Super Learner, van der Laan 2007; cross-fitting / CV-TMLE, Zheng and van der Laan 2011]
10. The feasibility rule: an estimand is feasible when some candidate meets the declared tolerances on every simulated outcome model, |bias| at most 0.02 on the risk difference scale and 95% CI coverage at least 0.90, judged on point estimates; the candidate with the smallest worst-case RMSE is selected; an optional interval half-width tolerance. Why these tolerances (prespecified, on the decision scale).
11. Monte Carlo uncertainty: Monte Carlo SE for bias, Wilson interval for coverage; a cell is unresolved when its interval contains a tolerance; unresolved cells that could change the decision get further repetitions up to a cap; decisions resting on unresolved cells are flagged in the SAP. [cite: Morris, White and Crowther 2019 simulation studies]
12. Pre-outcome precision: expected CI half-width and expected events per arm from the simulated models and the observed treatment.

## 2.4 The prespecified order of estimands and the decision

13. The analyst declares the estimands in order of preference; the primary estimand is the first one the simulation judges feasible; the others are reported as secondary. Why ordering before outcomes matters (selection on outcomes invalidates inference; this is claim 3). The populations each estimand describes are reported (covariate means and standardized mean differences), with small-cell suppression.

## 2.5 Negative-control outcomes

14. Each negative-control outcome is estimated for every estimand, on the full cohort and under declared restrictions; graded by its 95% CI against a declared null band (inside, outside, inconclusive); verdict per estimand GO, STOP or FLAG; only the full-cohort verdict for the primary estimand can block unblinding. Why on the primary estimand: a control read on an estimand the design cannot deliver mixes extrapolation error with confounding. [cite: Lipsitch 2010 negative controls]

## 2.6 Estimation after approval

15. TMLE for the ATE, trimmed ATE and ATT; TMLE for the ATO with overlap weights; influence-function-based CIs; super learner with one learner per family (elastic net, MARS, gradient-boosted trees, a small neural network) and an edge check on tuned grids. Weighting and matching as comparators, not stress-tested. [cite: van der Laan and Rubin 2006; Gruber and van der Laan; overlap-weighted TMLE source]

## 2.7 Simulation study (ADEMP)

16. Aims: claims 1 to 4 (from the brief).
17. Data-generating processes: four bounded covariates; treatment model with one coefficient on the main confounder set to give good, moderate and severe overlap [R: calibrated coefficients and shares outside the band]; three true outcome models (inside, close to, and outside the simulated family), each with its reason; null versions; n = 500 and 2,000. Provisional until the pilot (simulation/PLAN.md section 4).
18. Estimands and the answer key: population truths from one large draw; an estimand is estimable in a scenario if some candidate meets the tolerances on the true outcome model over [R: M] datasets.
19. Methods compared: the outcome-blind plasmode simulation (generate-treatment), the sample-treatment variant, and overlap diagnostics at their best thresholds (share outside the band, largest weight, effective sample size fraction, preference-score equipoise [cite: Walker 2013]); for claim 3, the prespecified choice versus the estimand with the smallest p-value after outcomes.
20. Performance measures and repetitions: sensitivity and specificity with Monte Carlo SEs; ROC curves; coverage and rejection rate; half-width ratio; [R: R map runs and M datasets per scenario, set by the pilot]; learners: GLM across the grid, the default super learner on a subset.

## 2.8 Application: the Rescue.Co Kenya trauma registry

21. Setting, data source and dates; cohort and eligibility; exposure contrasts in words (from the glossary); outcome and its timing; covariates and why. [R: cohort counts]
22. The plan declared before outcome access, the order of estimands, the negative controls and restrictions, and that this is a re-execution under masking: the analyst had seen outcomes in a parallel unblinded analysis.
23. Missing outcomes: complete case. Because incomplete rows are dropped before estimation, the estimate targets the effect in patients with an observed outcome; it equals the effect in the full cohort only if outcome missingness is unrelated to the outcome and to the covariates. This is a change from the protocol's inverse probability of censoring weighting, which assumes missingness depends only on measured covariates and treatment; that analysis is reported as a sensitivity analysis. Deviations from the protocol listed. [R: number and share of missing outcomes by arm]

## 2.9 Software

24. R version, cleanTMLE version, and the packages it uses for estimation (tmle, WeightIt, MatchIt, xgboost, glmnet, earth, nnet); code and simulation scripts available [repository].

## Display plan (Methods)

- Table S1: notation.
- Figure 1 (could sit in Methods): the staged workflow as a diagram, with what is known at each stage. Message: the outcome enters only after approval.

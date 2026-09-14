# Target trial emulation and the R causal software landscape: what cleanTMLE should borrow

Research notes compiled 2026-09-13 for the cleanTMLE interface and visualization revision. Scope: TTE methods and reporting standards, the API design of the TTE and causal-pipeline R ecosystem, balance/overlap/weight diagnostics and their plot idioms, positivity-specific visualization, TMLE reporting conventions, and estimand-switching practice. Because the analysis machine has no CRAN access, every recommendation below is framed as an idiom to reimplement natively in base R plus ggplot2, not a dependency to add.

## 1. Target trial emulation core and reporting

### 1.1 The framework papers

Hernan MA, Wang W, Leaf DE. Target Trial Emulation: A Framework for Causal Inference From Observational Data. JAMA. 2022;328(24):2446-2447. (PubMed 36508210; https://pubmed.ncbi.nlm.nih.gov/36508210/). The JAMA Guide to Statistics and Methods piece fixes the seven protocol components that every emulation must specify: eligibility criteria, treatment strategies, treatment assignment, time zero and follow-up, outcome, causal contrast, and analysis plan. Its central claims for cleanTMLE: the protocol is written before outcomes are analyzed, and mishandling time zero is the dominant self-inflicted bias. The original methods paper is Hernan MA, Robins JM. Using Big Data to Emulate a Target Trial When a Randomized Trial Is Not Available. Am J Epidemiol. 2016;183(8):758-764 (https://academic.oup.com/aje/article-abstract/183/8/758/1739860). Its Table 1 is the canonical "emulation table" (see 1.4).

What cleanTMLE should borrow: the design-object fields should map one-to-one onto the seven protocol components, so the printed design spec reads as a target trial protocol rather than as a list of function arguments.

### 1.2 The TARGET reporting guideline (published, final)

Cashin AG, Hansford HJ, Hernan MA, et al. Transparent Reporting of Observational Studies Emulating a Target Trial: The TARGET Statement. JAMA. 2025;334(12):1084-1093. doi:10.1001/jama.2025.13350 (https://jamanetwork.com/journals/jama/fullarticle/2837724). Status: final and published September 2025, developed under EQUATOR (systematic review, two-round survey 2023-2024, consensus meeting June 2024, pilot with 108 stakeholders through February 2025). PLOS Medicine already requires it for new TTE submissions (https://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1004796), and other journals are following.

Structure: 21 items across Abstract, Introduction, Methods, Results, Discussion, Other Information. The load-bearing items for software are 6 and 7, which run in parallel columns: item 6 specifies the target trial protocol across the component domains (6a eligibility, 6b treatment strategies, 6c assignment, 6d follow-up, 6e outcomes, 6f causal contrasts, 6g identifying assumptions including a positive probability of receiving each strategy, 6h analysis plan), and item 7 describes how each component was emulated with observational data. Item 14 covers sensitivity and supplementary analyses including negative controls and quantitative bias analysis.

What cleanTMLE should borrow: an exporter that renders the design object as a TARGET-compliant items 6-7 table plus a checklist stub. A package that emits this table automatically removes the single largest reporting burden of a TTE manuscript and makes the clean-room protocol audit-ready by construction.

### 1.3 Applied and systematic reviews of TTE quality

Three reviews matter, and they agree that protocol reporting is the weak point.

Hansford HJ, et al. Reporting of Observational Studies Explicitly Aiming to Emulate Randomized Trials: A Systematic Review. JAMA Netw Open. 2023 (RTI copy: https://www.rtihs.org/sites/default/files/34062_Garcia_2024_Reporting%20of%20observational%20studies%20explicitly%20aiming%20to%20emulate%20randomized%20trials.pdf). Roughly 200 emulations reviewed; only about 58 percent fully described their protocol components. This review seeded TARGET.

Wang C, Tang D, von Dadelszen P, et al. Concordance between target trial emulation and randomised controlled trials: systematic review and meta-analysis. BMJ. 2026;393:e086810 (https://pmc.ncbi.nlm.nih.gov/articles/PMC13184834/). 106 emulation-RCT pairs from 49 studies. Overall Pearson correlation between emulation and trial estimates 0.58, rising to 0.83 for close emulations. Only 35 percent of emulation studies registered protocols and 33 percent followed a reporting guideline. Concordance improved with good outcome emulation, similar baseline age and sex distributions, and registry or linked data.

Desai RK, et al. Scoping review: Investigating target trial emulation approaches in oncology research. Cancer Epidemiology. 2026 (https://www.sciencedirect.com/science/article/abs/pii/S1877782126000317; PubMed https://pubmed.ncbi.nlm.nih.gov/41762536/). 94 oncology TTEs, 61 percent published 2024-2025. All reported treatment strategies, outcomes, and an analysis plan, but only 75.5 percent presented a clearly defined time zero. Note: the review the prompt tentatively attributed to "Scola" appears to be this Desai scoping review; no TTE quality review by a Scola was found. The "Matthews" thread resolves to Anthony Matthews' benchmarking work with Dahabreh and colleagues (Dahabreh IJ, Matthews A, Steingrimsson JA, Scharfstein DO, Stuart EA. Using Trial and Observational Data to Assess Effectiveness: Trial Emulation, Transportability, Benchmarking, and Joint Analysis. Epidemiologic Reviews. 2024) and the benchmarking framework in Dahabreh et al., arXiv:2203.14857 (https://arxiv.org/pdf/2203.14857).

What cleanTMLE should borrow: the failure modes these reviews document are exactly what the clean room design report should prevent by default: explicit time zero, registered protocol before outcome unlock, and a machine-generated protocol table.

### 1.4 The emulation table format

The convention from Hernan and Robins 2016, now codified by TARGET items 6-7, is a table with one row per protocol component and two or three columns: protocol component, target trial specification, emulation with the available data (plus, in the best examples, a column noting deviations and their expected direction of bias). Good TTEs also add rows for identifying assumptions and the analysis plan, which is where estimand and positivity reporting live: the causal contrast is named (ITT analog, per-protocol analog), the target population is named, and positivity is addressed either as an eligibility restriction (redefine the population so all strategies are plausible) or as a diagnostics-plus-response statement (propensity overlap assessed, trimming or weighting bounds prespecified). The Gruber 2023 RWE reanalysis (section 5.3) shows the current best practice: a table contrasting the original and refined covariate categorizations after positivity checks, and a propensity overlap figure in the main text.

What cleanTMLE should borrow: emit this table from the design object with a fourth column recording the cleanTMLE argument that operationalizes each row. The estimand row should name the target population explicitly, and the assumptions row should link to the design report artifacts (overlap figure, PoRT-style table, ESS accounting).

## 2. R packages for TTE and causal pipelines: API design

### 2.1 TrialEmulation (Su, Rezvani, Seaman, Starr, Gravestock)

CRAN: https://cran.r-project.org/package=TrialEmulation; repo: https://github.com/Causal-LDA/TrialEmulation; paper: Su L, Rezvani R, Seaman SR, Starr C, Gravestock I. TrialEmulation: An R Package to Emulate Target Trials for Causal Analysis of Observational Time-to-event Data. arXiv:2402.12083 (https://arxiv.org/abs/2402.12083). Sequential trial emulation for time-to-event data: expansion into a sequence of trials, IPW for switching and censoring, pooled logistic MSM, marginal ITT and per-protocol effects.

The instructive part is that it ships two generations of API. The legacy interface is a monolithic initiators() call plus data_preparation() and trial_msm(). The modern interface is estimand-first and staged: trial_sequence("ITT" | "PP" | "AT") creates an estimand-specific object, then set_data(), set_switch_weight_model(), set_censor_weight_model(), set_outcome_model(), set_expansion_options() configure it, then expand_trials(), calculate_weights(), fit_msm() execute, with show_weight_models() to inspect the weight fits and predict() for marginal cumulative incidence with bootstrap CIs. The maintainers moved away from the all-in-one function because staged configuration gives better validation, incremental error messages, and inspectable intermediate state.

Borrow: the estimand-named constructor (the estimand is the first thing the user declares, and it changes which downstream setters are required); the set_* configuration verbs on a spec object; show_weight_models() as a named inspection verb; and the general lesson that a monolith with 30 arguments is the design to avoid.

### 2.2 causalRisk (NoviSci / Target RWE, commercial; mimic the grammar only)

Docs: https://docs.novisci.com/causalRisk/articles/estimator_check.html and https://docs.novisci.com/causalRisk/reference/make_table1.html. The grammar is the cleanest spec-estimate-report separation in the ecosystem:

models <- specify_models(identify_treatment(Statin, ~W), identify_censoring(EndofEnrollment), identify_outcome(Death))
fit <- estimate_ipwrisk(data, models, times = ..., trim = ..., label = "Primary ITT")
make_table1(fit, DxRisk, Sex, smd = TRUE); make_table2(fit1, fit2, risk_time = 20); plot(fit)

Three ideas carry the design. First, variable roles are declared with identify_*() verbs inside a single specification object, so the model spec exists and prints before any estimation. Second, every estimation call takes a label argument, and labels flow into every table and plot, which makes multi-analysis comparison tables trivial. Third, the reporting layer is two canonical functions: make_table1 (weighted and unweighted baseline characteristics with SMDs) and make_table2 (effect estimates at named times), plus a plot method for risk curves. estimate_aipwrisk() swaps the estimator without touching the spec. There are also identify_competing_risk(), identify_subject(), identify_interval() for longitudinal structure.

Borrow: the identify_*() role-declaration verbs, the label argument threaded through all outputs, and the make_table1/make_table2 pairing as the canonical table API. For the clean room this maps naturally onto a stage split: specify and table1 run outcome-blind, table2 runs after unlock.

### 2.3 targeted (Klaus Holst)

CRAN: https://cran.r-project.org/package=targeted; site: https://kkholst.github.io/targeted/. AIPW estimators for ATE, CATE via cate(), risk regression (riskreg) for RD and RR, assumption-lean GLM inference. Its distinctive contribution is the learner R6 class: learner_glm(), learner_grf(), learner_xgboost(), learner_sl(), learner_expand_grid(), each with estimate(), predict(), cv() methods. Nuisance models are passed as learner objects through response.model and treatment.model arguments, and cross-fitting is a fold count. Returned estimate objects carry influence-function inference and can be merged and transformed.

Borrow: the idea that nuisance model specification is an object with a uniform contract, not a string or a formula convention. cleanTMLE's learner registry (whatever it can support offline) should present one constructor per learner plus a uniform predict contract, and the fit object should retain the fitted learners for diagnostics.

### 2.4 lmtp (Williams and Diaz)

CRAN: https://cran.r-project.org/package=lmtp; repo: https://github.com/nt-williams/lmtp; book-style docs: https://www.beyondtheate.com/. Four estimators with identical signatures (lmtp_tmle, lmtp_sdr, lmtp_sub, lmtp_ipw) over data, trt, outcome, baseline, time_vary, cens, shift, mtp, folds, learners. The shift argument takes a policy function, which is the cleanest abstraction in the ecosystem for "the intervention is a function, not a level." lmtp_contrast() compares fitted policies on chosen scales, output objects have print and tidy methods, and progressr gives progress reporting across folds.

Borrow: identical signatures across estimators so switching estimators is a one-token change; tidy() methods returning a one-row-per-estimand data frame; lmtp_contrast() as the pattern for post-fit contrasts; progress reporting for long SL fits. The shift-function abstraction is also the right way for cleanTMLE to express stochastic or trimmed interventions if it ever moves beyond static contrasts.

### 2.5 tmle3 / tlverse

Repo and handbook: https://tlverse.org/tlverse-handbook/tmle3.html. tmle3 is maximally general: spec objects (tmle_ATE(), tmle_TSM_all()) bundle the estimand; a nodes list (W, A, Y) declares the NPSEM; learner lists per nuisance feed sl3. The fit prints psi_transformed, se, Wald CI, and the initial (pre-targeting) estimate next to the targeted one, and it can estimate all treatment-specific means plus contrasts simultaneously with joint inference. Weaknesses worth learning from: the R6/OO surface is heavy for applied users, sl3 is a hard dependency, and the handbook chapter presents no positivity or weight diagnostics at all, which is a real gap given how central truncation is in TMLE practice.

Borrow: the spec-object idea in a lighter S3 form; printing the initial vs targeted estimate side by side (a useful honesty diagnostic showing how much targeting moved the estimate); simultaneous TSM plus contrast reporting. Avoid: making users learn a class hierarchy to run one ATE.

### 2.6 ltmle and the classic tmle package

ltmle (https://cran.r-project.org/package=ltmle) implements longitudinal TMLE with node-based specification (Anodes, Cnodes, Lnodes, Ynodes, abar), gbounds truncation defaulting to c(0.01, 1), and summaries that report TMLE and IPTW estimates side by side. Its interface shows the cost of positional node conventions: powerful but easy to misuse, with long argument lists.

tmle (Gruber and van der Laan; https://cran.r-project.org/web/packages/tmle/tmle.pdf; JSS 2012, tmle: An R Package for Targeted Maximum Likelihood Estimation) remains the reference point-treatment implementation. Two things matter for cleanTMLE. First, its print method reports the additive effect, the effect among the treated and controls, RR, and OR, each with IC-based CI and p value, in labeled plain-text blocks. Second, since version 1.5.0-1 its default lower propensity bound is the Gruber 2022 adaptive formula 5/(sqrt(n) ln n) (section 4.6), a defensible default that the package documents explicitly.

Borrow: multi-scale labeled effect blocks in print output, and the adaptive truncation default with the bound value echoed in output.

### 2.7 survtmle and concrete

survtmle (Benkeser and Hejazi; https://benkeser.github.io/survtmle/, https://github.com/benkeser/survtmle) computes covariate-adjusted cumulative incidence under right censoring and competing risks via survtmle(ftime, ftype, trt, adjustVars, SL.trt, SL.ftime, SL.ctime, method = "hazard" or "mean", t0), with timepoints() to extend fits over a grid and a plot method for adjusted cumulative incidence curves.

concrete (Chen, Rytgaard, Fong, Tarp, Petersen, van der Laan, Gerds; arXiv:2310.19197, https://arxiv.org/abs/2310.19197) implements continuous-time one-step TMLE for survival and competing risks and has an unusually disciplined three-verb API: formatArguments() validates and assembles all inputs into a checked argument object (and prints what it inferred), doConcrete() runs estimation, getOutput() produces estimates, risk curves, RD and RR with simultaneous confidence bands across time points and events.

Borrow: from concrete, the validate-then-run split where the argument checker is itself a first-class object the user can inspect (this is essentially a design report step); from survtmle, the timepoints pattern of extending a fit without refitting nuisances. Simultaneous confidence bands are the right default when cleanTMLE reports survival curves at multiple horizons.

### 2.8 riskRegression and adjustedCurves

riskRegression (Gerds et al.; https://cran.r-project.org/web/packages/riskRegression/refman/riskRegression.html) contributes two idioms. ate() estimates average treatment effects on risk scales for censored outcomes with G-formula, IPTW, or AIPW selected by an estimator argument on one function. Score() is a model-audit engine: it takes a list of candidate risk models and returns AUC, Brier score, calibration, with plotROC() and plotCalibration() companions. adjustedCurves (Denz et al.; https://github.com/RobinDenz1/adjustedCurves; arXiv:2402.15292, https://arxiv.org/abs/2402.15292) wraps 15 adjustment methods for survival curves and 7 for cumulative incidence behind one front door, adjustedsurv(data, variable, ev_time, event, method = ...), with uniform plotting, adjusted RMST, and curve-difference tests regardless of method.

Borrow: the single front door with a method or estimator argument and a method registry underneath, so gcomp vs IPW vs TMLE curves are one-token swaps with identical downstream reporting; Score() as the pattern for the outcome-blind propensity model audit (discrimination and calibration of g-hat, never of the outcome model).

### 2.9 AIPW and causaldrf

AIPW (Zhong, Kennedy, Bodnar, Naimi; https://yqzhong7.github.io/AIPW/; AJE 2021, https://academic.oup.com/aje/article/190/12/2690/6322284) is an R6 estimator with new() then fit() then summary(), k-fold cross-fitting, and two canonical diagnostic plots attached to the fitted object: plot.p_score() (propensity distributions by arm) and plot.ip_weights() (boxplots of truncated IP weights by arm). causaldrf (Galagate and Schafer; https://cran.r-project.org/package=causaldrf) estimates average dose-response functions for continuous treatments; relevant only as the reminder that estimand generality costs interface simplicity, and cleanTMLE is right to stay binary point-treatment.

Borrow: shipping the two positivity plots as methods on the fit object itself, so the diagnostic is one call away from any estimate and appears in the same namespace.

### 2.10 Packages that already produce design reports or diagnostics dossiers

The closest existing analogs to the cleanTMLE design report are in the OHDSI stack. CohortMethod (https://github.com/OHDSI/CohortMethod; diagnostics described in The Book of OHDSI, Method Validity chapter, https://ohdsi.github.io/TheBookOfOhdsi/MethodValidity.html and the MultipleAnalyses vignette https://ohdsi.github.io/CohortMethod/articles/MultipleAnalyses.html) computes a standard battery per analysis: preference score distributions and empirical equipoise, covariate balance before and after PS adjustment, attrition, power (minimum detectable RR), and negative-control based empirical calibration. Crucially, results tables carry an unblind column: effect estimates are withheld unless prespecified diagnostics pass, with published default thresholds (all absolute SMDs below 0.1; equipoise, the fraction of subjects with preference score in 0.3 to 0.7, above 0.2). This is the same governance idea as the cleanTMLE clean room, implemented at industrial scale. smdi (Weberpals et al.; https://janickweberpals.gitlab-pages.partners.org/smdi/; JAMIA Open 2024, https://academic.oup.com/jamiaopen/article/7/1/ooae008/7595634) produces a one-call diagnostics dossier for partially observed confounders, smdi_diagnose(), returning a compact table of three diagnostic domains. cobalt's bal.tab() is the de facto printable balance dossier for the matching/weighting world.

Borrow: the unblind gate as an explicit column in the decision summary (cleanTMLE already has decision_summary.csv; align its vocabulary with pass/fail per named diagnostic plus an overall unblind flag), and the one-verb dossier entry point.

## 3. Balance, overlap, and weight diagnostics packages and their visualizations

This section catalogs plot types for native ggplot2 reimplementation. None of these packages can be dependencies; the idioms are the deliverable.

### 3.1 cobalt (Greifer)

Site: https://ngreifer.github.io/cobalt/; vignette: https://cran.r-project.org/web/packages/cobalt/vignettes/cobalt.html; love.plot reference: https://ngreifer.github.io/cobalt/reference/love.plot.html.

- bal.tab(): the balance table. SMDs (with the option of population vs treated SDs), variance ratios, KS statistics; un/adjusted columns side by side; thresholds flag rows; handles clusters, multiple imputations, multiple weight sets in one table.
- bal.plot(): distributional balance for one covariate at a time; density or histogram by treatment group, faceted unadjusted vs adjusted; works for the propensity score itself, which makes it an overlap plot.
- love.plot(): the summary graphic. Design details worth copying exactly: one point per covariate per sample (unadjusted open circle, adjusted filled), dashed vertical threshold line (0.1 convention), var.order to sort covariates by unadjusted SMD (largest at top), abs = TRUE so everything reads rightward from zero, drop.distance to exclude the PS row, sample.names for arbitrary weight-set labels, optional lines connecting the same covariate across samples, and multi-panel display when several statistics (mean differences, KS, variance ratios) are requested at once. Love plot named for Thomas E. Love.

Borrow: the whole love.plot argument vocabulary for cleanTMLE's balance plot, and bal.tab's un/adjusted column pairing with threshold flags for the design report table.

### 3.2 halfmoon (r-causal; Barrett, D'Agostino McGowan)

Site: https://r-causal.github.io/halfmoon/; reference index: https://r-causal.github.io/halfmoon/reference/index.html; repo: https://github.com/r-causal/halfmoon; CRAN: https://cran.r-project.org/package=halfmoon; intro post: https://r-causal.github.io/r-causal-blog/posts/introducing-halfmoon/; design discussion: https://livefreeordichotomize.com/posts/2023-08-04-visual-diagnostic-tools-for-causal-inference/; companion book chapter: https://www.r-causal.org/chapters/09-evaluating-ps.

Complete function inventory from the pkgdown reference:

- Geoms: geom_mirror_histogram() and geom_mirror_density() (one arm plotted upward, the other reflected below the axis, with weighted overlays on the observed distributions); geom_ecdf() (weighted and unweighted ECDFs); geom_qq2() (two-sample QQ); geom_roc() (weighted ROC of treatment on PS); geom_calibration() (PS calibration with CIs); geom_love() (SMD dot-line plot).
- Plot functions: plot_mirror_distributions(), plot_qq(), plot_ess() (effective sample size bars per weight set), plot_balance() (multi-metric balance display), plot_model_roc_curve() and plot_model_auc() (weighted ROC and AUC as balance checks, where post-weighting AUC near 0.5 indicates balance), plot_model_calibration(), plot_stratified_residuals().
- Balance metrics: bal_smd(), bal_vr() (variance ratio), bal_corr(), bal_ks(), bal_energy() (energy distance, a whole-distribution multivariate metric), bal_qq(), bal_ess(), bal_prognostic_score() (balance on a prognostic score, an outcome-model-free proxy for confounding relevance), bal_model_roc_curve(), bal_model_auc().
- Check functions returning tidy tibbles across multiple weight sets at once: check_balance(), check_ess(), check_qq(), check_model_roc_curve(), check_model_auc(), check_model_calibration().
- Utilities: add_ess_header() (ESS instead of raw N in table headers), weighted_quantile(), stat_qq2(), stat_roc(); nhefs_weights demo data.

The architecture is the real lesson: check_* functions compute tidy data frames with one row per covariate per metric per weight set, plot_* functions and geoms render them, and a single call accepts multiple weight columns (.weights = c(w_ate, w_att, w_atm, w_ato)) so every diagnostic can compare candidate weighting schemes in one figure. Matching indicators are treated as just another weight column. The observed (unweighted) sample is always drawn alongside weighted versions.

Borrow: this is the primary template for cleanTMLE's native diagnostics layer. Reimplement geom_mirror_histogram, ECDF, love/SMD, ESS bars, weighted ROC-AUC balance, and PS calibration; keep the compute/plot split (a check table feeding a plot function) because the check tables are exactly what the clean-room audit trail wants to persist as CSV artifacts.

### 3.3 WeightIt (Greifer)

Reference: https://ngreifer.github.io/WeightIt/reference/summary.weightit.html and https://ngreifer.github.io/WeightIt/reference/weightit.html. summary() on a weightit object reports, per treatment group: weight ranges, the top 5 weights with unit ids, coefficient of variation, mean absolute deviation scaled by the mean, negative entropy, count of zero weights, and ESS before and after weighting. plot(summary(w)) draws weight distribution histograms with a dotted line at the mean weight, and for ATT only the non-focal group is drawn since focal weights are 1. The estimand argument ("ATE", "ATT", "ATC", "ATO", "ATM", "ATOS") is validated per method, and the documentation for each method states which estimands it can target, which operationalizes the estimand-to-method mapping.

Borrow: the weight summary statistic set verbatim (ranges, CV, top-k with ids, zero count, ESS by arm) as cleanTMLE's weight diagnostic table, and the validated estimand argument whose allowed values depend on the chosen scheme.

### 3.4 PSweight (Zhou, Matsouaka, Thomas, Li)

R Journal paper: https://journal.r-project.org/articles/RJ-2022-011/; CRAN: https://cran.r-project.org/web/packages/PSweight/PSweight.pdf; repo: https://github.com/thuizhou/PSweight. Two-object design: SumStat() is the outcome-free diagnostics object (ps.formula or externally supplied scores, weight scheme, trimming delta, estimation method glm/gbm/SuperLearner) and PSweight() is the estimator. plot.SumStat() gives type = "balance" (absolute standardized differences and population standardized differences across covariates with a threshold line, one panel per weighting scheme), type = "density" (PS densities by arm), type = "hist" (PS histograms, binary only). The paper's framework table maps tilting function h(x) to target population and estimand: h = 1 gives ATE, h = e(x) gives ATT, h = e(x)(1 - e(x)) gives ATO, h = min(e, 1 - e) gives ATM, entropy tilting gives ATEN.

For the ATO with augmentation, which is the detail requested: PSweight(ps.formula, yname, weight = "overlap", augmentation = TRUE, out.formula = ...) implements the augmented overlap-weighted estimator of Mao, Li, and Greene (2018), combining the overlap-weighted outcome means with outcome regression predictions. The returned object contains the estimated average potential outcomes per arm (muhat), their variance-covariance matrix from nuisance-adjusted sandwich estimation (accounting for estimation of both the propensity and outcome models) or bootstrap (default R = 50) when bootstrap = TRUE, plus the weights and propensity scores. summary.PSweight() then forms contrasts on DIF, RR (log), or OR (log) scales with delta-method CIs. Trimming is a separate delta argument, applied before weighting, and the package refits the PS after trimming.

Borrow: the outcome-free SumStat object as the model for cleanTMLE's stage-2/3 design object (diagnostics exist as a class with plot and print methods, before any outcome is touched); the h(x)-to-population table as documentation and as metadata on each estimand option; nuisance-adjusted variance reporting for the augmented ATO fallback; and the convention that trimming triggers PS refitting.

### 3.5 MatchIt diagnostics

Vignette: https://kosukeimai.github.io/MatchIt/articles/assessing-balance.html; plot reference: https://kosukeimai.github.io/MatchIt/reference/plot.matchit.html. summary.matchit() reports SMDs, variance ratios, and eQQ statistics (mean and max distance between empirical quantile functions). plot.matchit() types: "qq" (empirical QQ per covariate), "ecdf", "density", "histogram", and "jitter", which is the useful one for positivity: a strip plot of propensity scores in four horizontal bands (matched treated, matched control, unmatched treated, unmatched control), points jittered, sized by matching weight, so discarded units and the common-support boundary are directly visible.

Borrow: the jitter plot as a "who is in, who is out" display once cleanTMLE trims; the eQQ max statistic as a cheap whole-distribution balance number.

### 3.6 twang

Vignette: https://cran.r-project.org/web/packages/twang/vignettes/twang.pdf. plot(ps.object) types: "optimize" (balance criterion vs GBM iteration, the tuning trace), "boxplot" (PS boxplots by arm), "es" (standardized effect sizes before and after weighting, connected), "ks" (KS statistics or p values before and after). Its stop.method idea (es.mean, ks.max) selects tuning by explicit balance criteria.

Borrow: the tuning-trace plot idea if cleanTMLE ever tunes g-model complexity against balance rather than likelihood; otherwise superseded by cobalt/halfmoon idioms.

### 3.7 tipr and EValue (sensitivity)

tipr (D'Agostino McGowan; JOSS 2022, https://joss.theoj.org/papers/10.21105/joss.04495; repo https://github.com/r-causal/tipr) computes tipping-point analyses with a unified grammar of the form {action}_{effect}_with_{confounder}: fix the confounder-exposure association and solve for the confounder-outcome association that nullifies the estimate, or adjust an estimate for a hypothesized confounder. EValue (VanderWeele and Mathur; https://cran.r-project.org/package=EValue; site https://louisahsmith.github.io/evalue/) computes E-values for RR, OR, HR, and standardized differences via evalue(), and bias_plot() draws the contour of confounder association pairs sufficient to explain away the estimate.

Borrow: report an E-value column alongside every unlocked effect estimate (it needs only the point estimate and CI, so it costs nothing), and consider a tipr-style helper for the sensitivity section of the manuscript pipeline.

### 3.8 smdi (missingness diagnostics)

Site: https://janickweberpals.gitlab-pages.partners.org/smdi/; paper: Weberpals J, et al. smdi: an R package to perform structural missing data investigations on partially observed confounders in real-world evidence studies. JAMIA Open. 2024;7(1):ooae008 (https://academic.oup.com/jamiaopen/article/7/1/ooae008/7595634). Three diagnostic domains for each partially observed confounder: (1) do observed characteristics differ between complete and incomplete units (ASMD comparisons, Hotelling and Little tests), (2) can missingness be predicted from observed data (random forest AUC), (3) is missingness associated with the outcome. smdi_diagnose() runs all three and returns one compact table.

Borrow: the domain-structured one-table dossier pattern, and the specific three-question framing for cleanTMLE's data-quality stage (the package already emits table_dq_coverage.csv; the smdi triad is a principled upgrade path).

## 4. Positivity-specific visualization and reporting

### 4.1 Mirrored propensity histograms and densities by arm

Construction (halfmoon idiom): histogram of PS for the treated with positive counts, histogram for controls with negative counts (count * -1), shared bins, weighted overlays drawn semi-transparent on top of the observed distributions, horizontal reference at zero. Real examples: halfmoon docs and the r-causal book chapter 9 (URLs in 3.2); Gruber et al 2023 Figure 3 uses an overlap density plot as the main positivity exhibit in a regulatory-grade reanalysis (https://pmc.ncbi.nlm.nih.gov/articles/PMC10394864/). cobalt bal.plot with the distance variable and PSweight type = "density" are the single-panel equivalents. cleanTMLE should draw the mirrored version with the trimming bounds as vertical lines and the excluded mass shaded (see 4.5).

### 4.2 Weight distribution plots

AIPW's plot.ip_weights() is the reference: boxplots of truncated IP weights by arm (ggplot2::geom_boxplot), which makes arm asymmetry in weight tails obvious (https://yqzhong7.github.io/AIPW/). WeightIt's histogram with a dotted mean line plus its printed top-5 weights per arm is the complement. Recommended composite for cleanTMLE: boxplot or violin of weights by arm, annotated with max weight, CV, and ESS per arm, and a caption stating the truncation bound in force. The clever covariate H(A, W) in TMLE is the same quantity up to sign, so cleanTMLE's existing fig_clever_covariate.png can adopt this design directly.

### 4.3 ESS displays

halfmoon plot_ess() draws ESS as grouped bars per weight set and arm; bal_ess()/check_ess() return the numbers; add_ess_header() substitutes ESS for N in table headers, a small but excellent honesty device for weighted Table 1s. WeightIt prints ESS before and after weighting per arm (Kish formula). cleanTMLE should report ESS in every weighted artifact header and plot ESS bars per candidate estimand in the design report.

### 4.4 Who is trimmed: profile tables and plots

Practice in pharmacoepidemiology: report trimming rule and percentiles (for example the GLORIA-AF program reported asymmetric trimming at the 1.5th percentile of the exposed and 98.5th of the comparator PS distribution; https://cdn.clinicaltrials.gov/large-docs/07/NCT01671007/SAP_001.pdf), the number removed per arm, and a comparison of removed vs retained characteristics. Methods anchors: Sturmer asymmetric trimming (trim below the 5th percentile of the treated PS and above the 95th of the untreated); Crump optimal subset (commonly PS in 0.1 to 0.9; Crump, Hotz, Imbens, Mitnik, Biometrika 2009), which WeightIt exposes as the ATOS estimand; comparison study https://pmc.ncbi.nlm.nih.gov/articles/PMC11476304/ (alternative trimming approaches) and https://dx.doi.org/10.1093/aje/kwab041 (weighting and trimming strategies simulation). Visual: MatchIt's jitter plot showing discarded units. Recommended cleanTMLE artifact: a trimmed-vs-kept table with one row per covariate showing mean in kept, mean in trimmed, and SMD between them, plus counts by arm, so the design team can state in words who the trimmed estimand no longer describes.

### 4.5 Common-support region shading

Construction: on the mirrored histogram, draw vertical lines at the enforced bounds (truncation or trimming) and shade the excluded region (annotate rect with alpha fill) on both half-axes; add the percentage of each arm excluded as text. The equipoise variant from OHDSI: transform PS to the preference score (Walker et al. A tool for assessing the feasibility of comparative effectiveness research. Comparative Effectiveness Research. 2013; https://www.dovepress.com/a-tool-for-assessing-the-feasibility-of-comparative-effectiveness-rese-peer-reviewed-fulltext-article-CER), shade the 0.3 to 0.7 equipoise band, and report the fraction of subjects inside it (OHDSI passes the diagnostic when that fraction exceeds 0.2; https://ohdsi.github.io/TheBookOfOhdsi/MethodValidity.html). Multigroup extension: Yoshida et al., Pharmacoepidemiol Drug Saf 2019 (https://onlinelibrary.wiley.com/doi/10.1002/pds.4767).

### 4.6 Bias vs truncation curves (hockey sticks)

Gruber S, Phillips RV, Lee H, van der Laan MJ. Data-Adaptive Selection of the Propensity Score Truncation Level for Inverse-Probability-Weighted and Targeted Maximum Likelihood Estimators of Marginal Point Treatment Effects. Am J Epidemiol. 2022;191(9):1640-1651. doi:10.1093/aje/kwac087 (https://academic.oup.com/aje/article/191/9/1640/6580570). The simulation figures plot bias, variance, and MSE against the truncation level for n in {100, 1000, 10000}: bias rises as bounds tighten, variance falls, and MSE traces the hockey-stick shape whose elbow the adaptive bound 5/(sqrt(n) ln n) tracks. The bound shrinks toward zero as n grows, and it became the default in the tmle package (1.5.0-1). Recommended cleanTMLE plot: estimate with CI (and, in simulations, bias and coverage) on the y axis against a grid of truncation levels, with the adaptive bound marked; in the applied setting this doubles as a stability-vs-truncation sensitivity display. Related: Leger et al., Causal inference in case of near-violation of positivity: comparison of methods. Biometrical Journal. 2022 (https://onlinelibrary.wiley.com/doi/abs/10.1002/bimj.202000323), which found truncation limits bias under near-violations but introduces bias when positivity holds, an argument for reporting the whole curve rather than one point.

### 4.7 PoRT-style subgroup violation tables

Danelian G, Foucher Y, et al. Identification of in-sample positivity violations using regression trees: The PoRT algorithm. Journal of Causal Inference. 2023 (https://www.degruyterbrill.com/document/doi/10.1515/jci-2022-0032/html); implemented as port() in the RISCA package (https://rdrr.io/cran/RISCA/man/port.html). Mechanics: fit one shallow tree per covariate, flag leaves that are both large enough (subgroup at least alpha of the sample, alpha default around 0.05) and extreme enough (exposure prevalence below beta or above 1 - beta, beta for example 0.05); remove flagged covariates, repeat over pairs, then trios, up to depth gamma (default 2). Output is a plain-language table: subgroup definition (covariate cuts), subgroup size, proportion exposed. The sequential extension sPoRT handles longitudinal regimes (https://arxiv.org/html/2412.10245, PubMed https://pubmed.ncbi.nlm.nih.gov/40856329/). This is the single most clean-room-compatible positivity diagnostic because it needs only covariates and treatment, never outcomes, and its output is readable by clinicians. cleanTMLE can reimplement it with rpart-free recursive partitioning on its own or with a coarse exhaustive search over discretized covariates, since gamma is small.

### 4.8 Support maps: heatmaps over scenario grids

For the simulation side of the cleanTMLE manuscript, the reporting idioms are in rsimsum (Gasparini; https://ellessenne.github.io/rsimsum/, nested loop plot vignette https://ellessenne.github.io/rsimsum/articles/D-nlp.html): autoplot types include heat plots (performance metric as a fill over a scenario grid, for example bias or coverage by overlap severity times sample size), nested loop plots (all scenarios in lexicographic order on one axis; Rucker and Schwarzer's format, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4272778/), and zipper plots for coverage. Methodological anchor: Morris TP, White IR, Crowther MJ. Using simulation studies to evaluate statistical methods. Statistics in Medicine. 2019;38(11):2074-2102. cleanTMLE's scenario-grid figures (bias and gate operating characteristics across overlap scenarios) should adopt the heat plot with Monte Carlo SEs acknowledged, and zipper plots where coverage is the claim.

### 4.9 The numeric companion: parametric bootstrap positivity diagnostic

Petersen ML, Porter KE, Gruber S, Wang Y, van der Laan MJ. Diagnosing and responding to violations in the positivity assumption. Stat Methods Med Res. 2012;21(1):31-54 (https://journals.sagepub.com/doi/abs/10.1177/0962280210386207; open version https://biostats.bepress.com/ucbbiostat/paper269/). The parametric bootstrap estimates the finite-sample bias attributable to sparsity by simulating from an estimated data-generating distribution where the truth is known, then reporting estimator bias across bootstrap replicates. It remains the strongest single number to accompany the plots above, and it is fully outcome-blind if simulated outcomes are generated under a null or synthetic outcome model.

## 5. TMLE-specific reporting

### 5.1 tlverse / tmle3 presentation conventions

From the handbook chapter (https://tlverse.org/tlverse-handbook/tmle3.html): results tables show, per parameter, the initial substitution estimate and the targeted estimate side by side with SE and Wald CI (psi_transformed, lower, upper, init_est), and multiple parameters (all TSMs plus the ATE contrast) are estimated jointly. Super learner composition is reported through sl3 cross-validated risk tables: one row per candidate learner with its CV risk and ensemble coefficient, per nuisance (g and Q). Notable gap: the chapter presents no positivity, weight, or truncation diagnostics, so a package that pairs TMLE output with the section 3-4 diagnostics is ahead of the reference implementation. cleanTMLE should print initial vs targeted estimates together and emit SL weight tables (learner, CV risk, coefficient) for both g and Q as standard artifacts.

### 5.2 survtmle and the classic tmle package

survtmle reports per-arm cumulative incidence with IC-based SEs and plots adjusted incidence curves over timepoints() grids (https://benkeser.github.io/survtmle/articles/survtmle_intro.html). The tmle package prints labeled blocks per effect scale (additive ATE, ATT, ATC, RR, OR) with CIs and p values, plus the IC-based variance, and documents its adaptive truncation default (section 4.6). Both are IC-first: every estimate is accompanied by its influence curve, which is what makes the IC distribution plot (already in cleanTMLE as fig_ic_plot.png) the natural TMLE diagnostic: heavy IC tails signal near-positivity trouble and fragile inference, and the clever covariate distribution is the weight diagnostic in TMLE clothing.

### 5.3 The TL-SAP line of work

Gruber S, Lee H, Phillips R, Ho M, van der Laan M. Developing a Targeted Learning-Based Statistical Analysis Plan. Statistics in Biopharmaceutical Research. 2023;15(3):468-475 (https://www.tandfonline.com/doi/full/10.1080/19466315.2022.2116104). Principles: define the target parameter separately from the model and estimator; prespecify TMLE plus super learner entirely (library, folds, truncation) while remaining data-adaptive; prespecify a nonparametric sensitivity analysis for causal-gap violations; represent intercurrent events explicitly. The companion application, Gruber S, Phillips RV, Lee H, Concato J, van der Laan M. Evaluating and improving real-world evidence with Targeted Learning. BMC Med Res Methodol. 2023;23:178 (https://pmc.ncbi.nlm.nih.gov/articles/PMC10394864/), walks the roadmap on a ritodrine pulmonary-edema study and is the best published template for cleanTMLE reporting: age categories were coarsened after positivity checks (shown as an original-vs-refined table), the causal question was downgraded from dose-response to any-vs-none when identifiability failed, the SAP prespecified the SL libraries per nuisance (outcome: linear regression, BART, lasso; propensity: logistic regression, BART, GAM), V = 20 folds chosen for the event rate, PS truncation at 0.06 from the 5/(sqrt(n) ln n) rule, and a sensitivity analysis plotted estimate shifts across hypothetical causal-gap values. The umbrella prespecification framework is Dang LE, et al. A Causal Roadmap for Generating High-Quality Real-World Evidence. J Clin Transl Sci. 2023 (arXiv:2305.06850, https://arxiv.org/abs/2305.06850).

Borrow: cleanTMLE's clean_room_config.yml is already an SAP-shaped object; align its sections with the TL-SAP headings (estimand, identification, SL library per nuisance, folds, truncation rule, sensitivity analysis, estimand fallback rules) so the config file can be exported as SAP text.

### 5.4 CV-TMLE recommendations

Smith MJ, Phillips RV, Maringe C, Luque-Fernandez MA. Performance of Cross-Validated Targeted Maximum Likelihood Estimation. arXiv:2409.11265 (https://arxiv.org/abs/2409.11265; published version indexed at https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12270713/). Findings: CV-TMLE (cross-fitting the initial estimators) restores CI coverage under Donsker-class violations, small samples, and near-positivity violations, without adding bias; it is much less sensitive to the SL library than plain TMLE; adding regression trees to non-cross-validated TMLE inflated both bias and variance. Recommendation for software defaults: cross-fit whenever the library contains flexible learners, which for cleanTMLE means CV-TMLE should be the default estimator path (or at least the default whenever the library is not purely parametric), the number of folds should be reported in output, and the docs should warn against tree-based learners without cross-fitting.

## 6. Estimand switching and feasibility in practice

### 6.1 Greifer and Stuart's estimand-choice framework

Greifer N, Stuart EA. Choosing the Estimand When Matching or Weighting in Observational Studies. arXiv:2106.10577 (https://arxiv.org/abs/2106.10577v1); applied summary by the same author in the IQSS primer chapter Planning the Analysis (https://iqss.github.io/dss-ps/planning.html). The framework selects the estimand from substantive questions before any modeling: what is the policy question (would treatment be expanded to everyone, giving ATE; withheld from current recipients, giving ATT; extended to current non-recipients, giving ATC/ATU); who could realistically be intervened on; and only then, is the estimand feasible given overlap. When profiles present in one arm are absent in the other, the advice is to move to an overlap-defined population (ATO) or a trimmed sample rather than extrapolate, with the explicit caveat that the ATO population is defined statistically rather than substantively and must be characterized after the fact (describe the weighted covariate profile so readers know who the estimate is about). The chapter emphasizes documenting estimand, effect measure, time horizon, and target population before analysis. The ecosystem operationalization: WeightIt's estimand argument validates ATE/ATT/ATC/ATO/ATM/ATOS per method, and its docs state the allowed estimands per weighting scheme (https://ngreifer.github.io/WeightIt/reference/weightit.html).

### 6.2 The weights-to-target-population mapping

Li F, Morgan KL, Zaslavsky AM. Balancing covariates via propensity score weighting. JASA. 2018;113(521):390-400. The balancing-weights family: target population defined by tilting h(x), weights h(x)/e(x) and h(x)/(1 - e(x)); overlap weights h = e(1 - e) minimize asymptotic variance among all balancing weights and give exact finite-sample mean balance when the PS is logistic. Li F, Thomas LE, Li F. Addressing extreme propensity scores via the overlap weights. Am J Epidemiol. 2019;188(1):250-257 (https://academic.oup.com/aje/article/188/1/250/5090958) argues OW dominates trimming under poor overlap (no arbitrary cutpoint, no refitting, smooth downweighting). The PSweight paper's h(x) table (section 3.4) is the cleanest published mapping of scheme to population to estimand, and the figures in these papers showing the tilting functions over the PS axis (ATE flat, ATT proportional to e, ATO the inverted parabola peaking at 0.5) are the standard way to visualize what population each estimand describes. cleanTMLE should reproduce that one-panel figure (weight functions over the PS axis, shaded by the achieved sample density) in its estimand documentation and design report.

### 6.3 Prespecified fallback ladders

Petersen et al 2012 (section 4.9) gives the canonical ordered menu of responses to positivity failure: (1) restrict the covariate adjustment set, (2) redefine the parameter through a marginal structural working model projection, (3) restrict the sample (trimming, with the estimand becoming the effect in the retained population), (4) change the target intervention (dynamic or stochastic regimes, or an overlap-tilted population). The TL-SAP and Causal Roadmap papers (5.3) turn this into prespecification practice: the SAP names the primary estimand and the diagnostic thresholds that trigger each fallback before outcomes are seen. The Gruber 2023 BMC paper is a worked example of executing such a switch credibly. OHDSI's unblind gate (section 2.10) is the same governance made mechanical: named diagnostics with pass thresholds, and effect estimates withheld until they pass. Crump's optimal subset (ATOS in WeightIt) is the formal version of "trimmed ATE" as a first fallback; ATT is the natural second when the treated are well supported inside the comparator; ATO is the terminal fallback that always exists.

What cleanTMLE should borrow: encode the ladder in config as an ordered list of (estimand, trigger, threshold) entries, for example primary ATE; if min ESS fraction or equipoise fails, trimmed ATE with Crump or Sturmer rule named; if still failing, ATT; terminal ATO with augmentation. The design team's stage decision then reduces to reading the design report against prespecified triggers, and decision_summary.csv records which rung fired and why. Each rung must also swap the reported target-population description, reusing the 6.2 figure and a weighted covariate profile table of the achieved population.

## Synthesis: the 12 most concrete adoptions, ranked

1. Design report with an explicit unblind gate. From OHDSI CohortMethod and The Book of OHDSI Method Validity chapter. Give cleanTMLE one verb (design_report()) that computes the outcome-free battery: mirrored PS overlap, balance table with SMD threshold 0.1, equipoise fraction with the 0.3 to 0.7 preference band, ESS by arm, weight summary, PoRT table, attrition. Emit a pass/fail row per diagnostic and one overall unblind flag, and have the stage-5 estimation functions refuse to run (or loudly warn) unless the gate object says pass or a documented override is supplied. This is the OHDSI unblind column translated into the clean-room vocabulary cleanTMLE already uses in decision_summary.csv.

2. A prespecified estimand fallback ladder in config. From Petersen et al 2012, the TL-SAP papers, and Greifer and Stuart. Represent the fallback sequence as ordered (estimand, trigger, threshold) entries in clean_room_config.yml, evaluated by the design report; record the rung that fired and the resulting target population in the decision summary. Ship the PSweight h(x) table and the tilting-function figure as the documentation of what each rung estimates, and always output a weighted covariate profile of the achieved population after a switch.

3. The spec-then-estimate-then-report grammar with role verbs and labels. From causalRisk (specify_models plus identify_treatment/identify_outcome/identify_censoring, label arguments, make_table1/make_table2) and TrialEmulation's estimand-first trial_sequence() with set_* configuration verbs. cleanTMLE's user-facing API should read: spec <- clean_design(data, identify_treatment(...), identify_covariates(...), estimand = "ATE"); diagnose(spec); then, after unlock, fit <- estimate_tmle(spec, outcome = ...). Every constructor takes a label that flows into all tables and figures. The two-generation history of TrialEmulation is the argument to give reviewers for why cleanTMLE is not one function with 30 arguments.

4. halfmoon's diagnostics architecture, reimplemented natively. From halfmoon: check_* functions that return tidy one-row-per-covariate-per-weight-set tables, and plot_* functions plus geoms that render them; multiple candidate weight sets or estimands compared in a single call; the observed sample always drawn alongside weighted versions. Persist the check tables as CSV artifacts (they are the audit trail) and build the figures from them. Reimplement at minimum: mirrored histogram, weighted ECDF, love plot, ESS bars, weighted ROC-AUC balance check, and PS calibration.

5. cobalt's love plot design vocabulary. From love.plot(): absolute SMDs sorted by unadjusted magnitude, open vs filled points for unadjusted vs adjusted, dashed threshold line at 0.1, sample.names for arbitrary labels, optional multi-statistic panels (SMD plus KS), drop.distance handling for the PS row. cleanTMLE's existing fig_love_plot.png should adopt this argument vocabulary so users coming from cobalt find the same knobs.

6. The weight diagnostics summary, verbatim from WeightIt plus AIPW's plot. Report per arm: min and max weight, top 5 weights with row ids, coefficient of variation, count of zero weights, ESS before and after. Pair with a boxplot of weights by arm annotated with the max weight and the truncation bound in force (AIPW's plot.ip_weights design). In TMLE terms this is the clever covariate summary, so the same table should back fig_clever_covariate.png.

7. Bias-vs-truncation curve plus the adaptive bound as default. From Gruber et al 2022 AJE and the tmle package default. Default the lower PS bound to 5/(sqrt(n) ln n), echo the realized bound and the number of truncated units per arm in print output, and provide one function that re-estimates over a truncation grid and plots estimate with CI against the bound (the applied hockey stick). In the simulation reporting, pair with rsimsum-style heat maps and zipper plots over the scenario grid.

8. A PoRT-style subgroup violation table. From Danelian and Foucher's PoRT (RISCA::port). Reimplement the shallow recursive search (alpha, beta, gamma hyperparameters, defaults 0.05, 0.05, 2) over covariates, pairs, and trios, and print a plain-language table of subgroup definition, size, and exposure prevalence. This is the design team's most communicable positivity artifact and is fully outcome-blind.

9. CV-TMLE as the default path with SL transparency tables. From Smith et al arXiv:2409.11265 and tlverse presentation conventions. Cross-fit by default whenever the learner library is not purely parametric, report folds, and warn on tree learners without cross-fitting. Emit SL weight tables (learner, CV risk, ensemble coefficient) for g and Q as standard artifacts, and print the initial vs targeted estimate side by side the way tmle3 does.

10. TARGET and emulation-table exporters. From the TARGET statement (JAMA 2025) and Hernan and Robins 2016 Table 1. Add report_target_table(spec): a table with one row per protocol component and columns for target trial, emulation, and the cleanTMLE argument or config key that operationalizes it, plus a checklist stub for the 21 items. With PLOS Medicine already requiring TARGET, this converts a compliance chore into a package feature and makes the clean-room protocol externally legible.

11. A trimmed-vs-kept profile artifact. From Sturmer and Crump trimming practice, MatchIt's jitter plot, and the GLORIA-AF style of reporting trim percentiles. Whenever trimming or an optimal subset fires, emit: counts removed per arm, the trimming rule and realized cutpoints, a covariate table of kept vs removed with SMDs between them, and the mirrored histogram with the excluded region shaded. This states in data who the new estimand no longer describes, which is the honest core of estimand switching.

12. SAP-shaped config aligned with the TL-SAP headings. From Gruber et al, Statistics in Biopharmaceutical Research 2023, and the Causal Roadmap. Restructure clean_room_config.yml sections to mirror the TL-SAP: estimand and target population, identification assumptions and diagnostics with thresholds, SL library per nuisance, folds, truncation rule, sensitivity analyses (E-value by default, optional tipr-style tipping point), and the fallback ladder from item 2. Provide an exporter that renders the config as SAP prose so the design document, the code configuration, and the manuscript methods section are one artifact.

## URLs consulted

- https://pubmed.ncbi.nlm.nih.gov/36508210/
- https://academic.oup.com/aje/article-abstract/183/8/758/1739860
- https://jamanetwork.com/journals/jama/fullarticle/2837724
- https://pubmed.ncbi.nlm.nih.gov/40899949/
- https://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1004796
- https://www.rtihs.org/sites/default/files/34062_Garcia_2024_Reporting%20of%20observational%20studies%20explicitly%20aiming%20to%20emulate%20randomized%20trials.pdf
- https://pmc.ncbi.nlm.nih.gov/articles/PMC13184834/
- https://www.sciencedirect.com/science/article/abs/pii/S1877782126000317
- https://pubmed.ncbi.nlm.nih.gov/41762536/
- https://arxiv.org/pdf/2203.14857
- https://arxiv.org/abs/2402.12083
- https://cran.r-project.org/package=TrialEmulation
- https://cran.r-project.org/web/packages/TrialEmulation/refman/TrialEmulation.html
- https://github.com/Causal-LDA/TrialEmulation
- https://docs.novisci.com/causalRisk/articles/estimator_check.html
- https://docs.novisci.com/causalRisk/reference/make_table1.html
- https://kkholst.github.io/targeted/
- https://cran.r-project.org/package=targeted
- https://cran.r-project.org/package=lmtp
- https://github.com/nt-williams/lmtp
- https://www.beyondtheate.com/
- https://tlverse.org/tlverse-handbook/tmle3.html
- https://cran.r-project.org/package=ltmle
- https://cran.r-project.org/web/packages/tmle/tmle.pdf
- https://benkeser.github.io/survtmle/
- https://benkeser.github.io/survtmle/articles/survtmle_intro.html
- https://github.com/benkeser/survtmle
- https://arxiv.org/abs/2310.19197
- https://cran.r-project.org/web/packages/riskRegression/refman/riskRegression.html
- https://github.com/RobinDenz1/adjustedCurves
- https://arxiv.org/abs/2402.15292
- https://yqzhong7.github.io/AIPW/
- https://academic.oup.com/aje/article/190/12/2690/6322284
- https://cran.r-project.org/package=causaldrf
- https://github.com/OHDSI/CohortMethod
- https://ohdsi.github.io/TheBookOfOhdsi/MethodValidity.html
- https://ohdsi.github.io/CohortMethod/articles/MultipleAnalyses.html
- https://ngreifer.github.io/cobalt/
- https://cran.r-project.org/web/packages/cobalt/vignettes/cobalt.html
- https://ngreifer.github.io/cobalt/reference/love.plot.html
- https://r-causal.github.io/halfmoon/
- https://r-causal.github.io/halfmoon/reference/index.html
- https://github.com/r-causal/halfmoon
- https://cran.r-project.org/package=halfmoon
- https://r-causal.github.io/r-causal-blog/posts/introducing-halfmoon/
- https://livefreeordichotomize.com/posts/2023-08-04-visual-diagnostic-tools-for-causal-inference/
- https://www.r-causal.org/chapters/09-evaluating-ps
- https://ngreifer.github.io/WeightIt/reference/summary.weightit.html
- https://ngreifer.github.io/WeightIt/reference/weightit.html
- https://journal.r-project.org/articles/RJ-2022-011/
- https://cran.r-project.org/web/packages/PSweight/PSweight.pdf
- https://github.com/thuizhou/PSweight
- https://kosukeimai.github.io/MatchIt/articles/assessing-balance.html
- https://kosukeimai.github.io/MatchIt/reference/plot.matchit.html
- https://cran.r-project.org/web/packages/twang/vignettes/twang.pdf
- https://joss.theoj.org/papers/10.21105/joss.04495
- https://github.com/r-causal/tipr
- https://cran.r-project.org/package=EValue
- https://louisahsmith.github.io/evalue/
- https://janickweberpals.gitlab-pages.partners.org/smdi/
- https://academic.oup.com/jamiaopen/article/7/1/ooae008/7595634
- https://academic.oup.com/aje/article/191/9/1640/6580570
- https://onlinelibrary.wiley.com/doi/abs/10.1002/bimj.202000323
- https://www.degruyterbrill.com/document/doi/10.1515/jci-2022-0032/html
- https://rdrr.io/cran/RISCA/man/port.html
- https://arxiv.org/html/2412.10245
- https://pubmed.ncbi.nlm.nih.gov/40856329/
- https://journals.sagepub.com/doi/abs/10.1177/0962280210386207
- https://biostats.bepress.com/ucbbiostat/paper269/
- https://www.dovepress.com/a-tool-for-assessing-the-feasibility-of-comparative-effectiveness-rese-peer-reviewed-fulltext-article-CER
- https://onlinelibrary.wiley.com/doi/10.1002/pds.4767
- https://cdn.clinicaltrials.gov/large-docs/07/NCT01671007/SAP_001.pdf
- https://pmc.ncbi.nlm.nih.gov/articles/PMC11476304/
- https://dx.doi.org/10.1093/aje/kwab041
- https://ellessenne.github.io/rsimsum/
- https://ellessenne.github.io/rsimsum/articles/D-nlp.html
- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4272778/
- https://www.tandfonline.com/doi/full/10.1080/19466315.2022.2116104
- https://pmc.ncbi.nlm.nih.gov/articles/PMC10394864/
- https://arxiv.org/abs/2305.06850
- https://arxiv.org/abs/2409.11265
- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12270713/
- https://arxiv.org/abs/2106.10577v1
- https://iqss.github.io/dss-ps/planning.html
- https://academic.oup.com/aje/article/188/1/250/5090958

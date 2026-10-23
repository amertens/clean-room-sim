# Outline: Introduction (Statistics in Medicine estimability paper)

**Gap (one sentence):** Analysts are told to fix the estimand before seeing outcomes, but whether a design can deliver an estimand depends on positivity in a way that overlap diagnostics cannot judge and the existing simulation checks either need the outcome or keep the observed treatment; it is not known whether an outcome-blind plasmode simulation can tell, before outcome access, which estimand the design can deliver, or whether choosing the estimand at that stage keeps inference valid.

**Aim type:** methodological. The objects of study are performance measures (agreement with estimability on the true outcome model; coverage and type I error), not a causal effect. The Rescue.Co application has a causal estimand, stated in the Methods.

Four paragraphs, no results, one sentence of methods.

## P1. The problem: estimands are fixed before outcomes, but estimability is only learned from the data

Topic sentence: Observational comparisons of treatments now inform clinical and regulatory decisions, and their credibility rests on fixing the question and the analysis before the outcome is seen, as a trial protocol would.
- Fixing the analysis in advance limits choices made after seeing results [@kahanHowDesignPrespecified2020]; target trial emulation and the causal roadmap ask for the estimand first [@hernan2016; @dangStartTargetTrial2023]; staging and clean-room constructs allow design changes decided before outcome access [@muntner2024; @rubin2008].
- Quantified fact for the opening (I-01): [cite: how often observational comparative studies change the analysis or estimand after seeing data, or how often positivity problems are found, with a number].
- Ends on the open problem: the estimand the protocol asks for may not be estimable in the data, and that is often discovered only after the outcome has been seen.

## P2. What is known: the estimands need different positivity, and propensity-score diagnostics cannot say which is estimable

Topic sentence: The common causal estimands differ in the population they describe and in how much overlap they need.
- ATE needs overlap everywhere; trimmed ATE only in the band; ATT only away from 1; ATO none [@crump2009; @li2017]; in estimand terms the population attribute changes [@internationalcouncilforharmonisationoftechnicalrequirementsforpharmaceuticalsforhumanuse2019; @chenEstimandsRealWorldEvidence2023].
- Near-violations of positivity bias estimators by amounts that depend on the outcome regression, not on the propensity score alone [@petersenDiagnosingRespondingViolations2012; @wiederkehr2026].
- Overlap diagnostics (share outside a band, largest weight, effective sample size, preference-score equipoise) read only the propensity score [@walker2013]; the parametric bootstrap that does judge the bias fits an outcome model to the outcome [@petersenDiagnosingRespondingViolations2012].
- Ends: no pre-outcome check judges estimability itself.

## P3. What is known: outcome-blind simulation is used to choose estimators, but its design matters and its verdicts have not been checked against the truth

Topic sentence: Simulations built on a study's own covariates can evaluate estimators without the outcome.
- Plasmode simulation [@franklin2014]; outcome-blind simulations in the causal roadmap to prespecify the estimator [@nanceCausalRoadmapSimulations2024a; @williamsonApplicationCausalRoadmap2023].
- Keeping each patient's observed treatment builds a positivity violation into the simulated data [@shaw2026]; plasmode simulations can refute an estimator but not certify it [cite: Karim 2026].
- Concurrent work chooses the estimand from positivity diagnostics [cite: Liu et al. 2025, Stat Med] or gates unblinding on design diagnostics [cite: Conover et al. 2025, JAMIA].
- Ends: whether an outcome-blind simulation's verdict on which estimand is estimable agrees with estimability on the true outcome model, and whether choosing the estimand from that verdict keeps CI coverage, has not been evaluated (stated as untested, not as "no study", unless a logged search bounds it).

## P4. Objective

Topic sentence: We propose choosing the estimand at the design stage with an outcome-blind plasmode simulation, inside a staged workflow implemented in the R package cleanTMLE.
- One sentence of method: new treatments drawn from the estimated propensity score, outcomes from a prespecified family, the first feasible estimand in a prespecified order becomes primary, and the outcome is joined only after the SAP is approved.
- Objective (I-05): to estimate, in a simulation study, how often each pre-outcome verdict agrees with estimability on the true outcome model, by estimand and by the true model's relation to the simulated family, and the coverage and type I error of the estimate when the estimand is chosen before versus after outcome access; and to show the workflow in a re-execution of a trauma registry analysis with the outcome withheld.
- Aim type stated: methodological evaluation.

## Citations still needed (run sci-cite)
- P1 quantified fact.
- Karim 2026; Liu et al. 2025; Conover et al. 2025 (named by the earlier novelty review; not yet in the Zotero collection).

## Notes
- The brief's claim 5 still describes an old-run result ("found the protocol's ATE not estimable"); it should be restated without a result until the 0.4.1 re-execution exists. No result appears in the Introduction.

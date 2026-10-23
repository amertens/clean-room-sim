---
product: manuscript          # manuscript | results-report | supplement | sap | concept-note | slides | tutorial | response | scope-memo
audience: methods-stats      # Statistics in Medicine readers
venue: "Statistics in Medicine"
limits: {words: null, abstract_words: null, abstract: null}   # check the current author guidelines before drafting
design: methods              # a methods paper: simulation study plus a real-data application
causal_estimand: true        # causal estimands (ATE, trimmed ATE, ATT, ATO); identifying assumptions stated in Methods
length_target: {}
collaborator_sections: []
status: draft
---

# One-sentence message
(proposed) Outcome-blind simulation, which has been used to choose the estimator for an estimand fixed in advance, can also judge before any outcome is seen which estimands a design can deliver, and an order of estimands fixed in advance turns that judgment into the choice of the primary estimand. The simulation study measures how often the judgment is right, when the true outcome model resembles the simulated ones and when it does not, and what choosing the estimand before rather than after seeing outcomes does to coverage and type I error; cleanTMLE makes the sequence checkable and applies after unblinding the estimator that the simulation evaluated.

# Key claims (3 to 5)
1. (proposed) The outcome-blind plasmode simulation's verdict for each estimand (ATE, trimmed ATE, ATT, ATO) agrees with whether the estimand is estimable on the true outcome model more often than overlap diagnostics at their best thresholds, when the true outcome model lies inside or near the simulated family; the agreement falls when it lies outside. (Table 1, Figure 1)
2. (proposed) Drawing new treatments from the estimated propensity score grades estimability correctly; keeping each patient's observed treatment misgrades it. (Table 1)
3. (proposed) Choosing the estimand from the prespecified order before outcome access keeps 95% CI coverage nominal; choosing it after seeing outcomes by the smallest P value inflates the type I error, while choosing it by the narrowest CI changes neither coverage nor the type I error by more than 3 percentage points. (Table S1)
4. (proposed) The pre-outcome expected CI half-width predicts the realized half-width. (Figure 2)
5. (proposed) In a re-analysis of the Rescue.Co trauma registry in Nairobi, Kenya, with the outcome withheld from the software, the design stage, run on the first four months of enrollment and on the final cohort, recorded before any outcome was joined which estimands each comparison could deliver and which estimator the simulation selected for each, and the estimate for the confirmatory comparison came from the estimator the simulation selected. At four months it found no estimand that the comparison between ambulance operators could deliver, although no patient in that comparison lay outside the overlap band. The analysts wrote the decision rules after seeing the results of an earlier, unblinded analysis of the same registry, so the application shows that the rules can be applied without the outcome, not that they would have made the same choices without that knowledge. (Confirm the estimate clause once the estimates exist, and align this claim with the case-study claims of recommendation 19 of the 22 October audit.)

# Objective or estimand
Methods paper. The objects of study are (a) the agreement between pre-outcome verdicts and estimability on the true outcome model, by scenario and estimand, and (b) the coverage and type I error of the estimate for an estimand chosen before versus after outcome access. Estimands in the simulation and the case study: risk differences for the ATE, the trimmed ATE (band 0.05 to 0.95), the ATT and the ATO. Design and settings: simulation/PLAN.md.

# What not to say
- Nothing about the data-quality stress test, threat families, the GO/FLAG/STOP stress verdict, or a STOP result under unmeasured confounding (removed from cleanTMLE 0.4.0; the old manuscript's results do not carry over).
- No min_max_rmse over threats, variance-method selection, matched-cohort TMLE bootstrap, or IPCW for missing outcomes as package features.
- Not "sixteen verbs", "outcome store", "masking", or "tamper-proof": the lock makes changes detectable (tamper-evident) and keeps the outcome out of the design stage; it cannot stop anyone analysing the outcome elsewhere.
- Not that the simulation certifies estimability: it supports it over the declared family of simulated outcome models only.
- A FLAG negative-control verdict is not evidence of no residual confounding.
- The case study re-executes the blinded design with the outcome withheld from the software; the analyst had seen the outcomes in the parallel unblinded analysis. Say so, using the word blinding.
- The Rescue.Co case study's missing outcomes are handled by complete case, the estimator the design stage evaluated; say why it is weaker here than a longitudinal TMLE (early deaths are always observed, only survivors can be lost), and report the longitudinal TMLE as a sensitivity analysis not evaluated in the design stage.
- Not that the identification table shows identification holds: it lists what was declared.
- No results from the 0.1.x to 0.3.0 runs (archive/); every number comes from the 0.4.1 runs through the registry.

# Notes for this document
- Reuse within the project is allowed (reports/manuscript_outcome_blind_dq.qmd, tutorials/, reports/cleanTMLE_presentation.qmd), but the old manuscript's structure and most of its claims do not carry over.
- Cite only keys in writing/references.bib (the project's Zotero collection). The novelty review found the closest prior and concurrent work to position against: Liu et al. 2025 (Stat Med) on positivity-driven estimand choice, Karim 2026 on plasmode simulations refuting but not certifying, Conover et al. 2025 (JAMIA) on gating unblinding on diagnostics, Nance et al. 2024 (Epidemiology), Shaw et al. 2025, Petersen et al. 2012, Li, Morgan and Zaslavsky 2018, Crump et al. 2009. Add them through sci-cite.

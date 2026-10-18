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
(proposed) Before any outcome is seen, an outcome-blind plasmode simulation can tell which causal estimands a design can deliver, more accurately than overlap diagnostics, when the true outcome model resembles the simulated ones, and choosing the estimand at that stage keeps inference valid where choosing it after seeing outcomes does not.

# Key claims (3 to 5)
1. (proposed) The outcome-blind plasmode simulation's verdict for each estimand (ATE, trimmed ATE, ATT, ATO) agrees with whether the estimand is estimable on the true outcome model more often than overlap diagnostics at their best thresholds, when the true outcome model lies inside or near the simulated family; the agreement falls when it lies outside. (Table 1, Figure 1)
2. (proposed) Drawing new treatments from the estimated propensity score grades estimability correctly; keeping each patient's observed treatment misgrades it. (Table 1)
3. (proposed) Choosing the estimand from the prespecified order of estimands before outcome access keeps 95% CI coverage nominal; choosing it after seeing outcomes, as the estimand with the narrowest CI, undercovers and inflates the type I error. (Table S1)
4. (proposed) The pre-outcome expected CI half-width predicts the realized half-width. (Figure 2)
5. (proposed) In the Rescue.Co Kenya trauma registry, the design stage found the protocol's ATE not estimable and moved to estimands the data could support, before any outcome was joined.

# Objective or estimand
Methods paper. The objects of study are (a) the agreement between pre-outcome verdicts and estimability on the true outcome model, by scenario and estimand, and (b) the coverage and type I error of the estimate for an estimand chosen before versus after outcome access. Estimands in the simulation and the case study: risk differences for the ATE, the trimmed ATE (band 0.05 to 0.95), the ATT and the ATO. Design and settings: simulation/PLAN.md.

# What not to say
- Nothing about the data-quality stress test, threat families, the GO/FLAG/STOP stress verdict, or a STOP result under unmeasured confounding (removed from cleanTMLE 0.4.0; the old manuscript's results do not carry over).
- No min_max_rmse over threats, variance-method selection, matched-cohort TMLE bootstrap, or IPCW for missing outcomes as package features.
- Not "sixteen verbs", "outcome store", "masking", or "tamper-proof": the lock makes changes detectable (tamper-evident) and keeps the outcome out of the design stage; it cannot stop anyone analysing the outcome elsewhere.
- Not that the simulation certifies estimability: it supports it over the declared family of simulated outcome models only.
- A FLAG negative-control verdict is not evidence of no residual confounding.
- The case study re-executes the blinded design with the outcome withheld from the software; the analyst had seen the outcomes in the parallel unblinded analysis. Say so, using the word blinding.
- The Rescue.Co case study's missing outcomes are handled by complete case analysis, a change from the protocol's IPCW; say so and report IPCW as a sensitivity analysis.
- No results from the 0.1.x to 0.3.0 runs (archive/); every number comes from the 0.4.1 runs through the registry.

# Notes for this document
- Reuse within the project is allowed (reports/manuscript_outcome_blind_dq.qmd, tutorials/, reports/cleanTMLE_presentation.qmd), but the old manuscript's structure and most of its claims do not carry over.
- Cite only keys in writing/references.bib (the project's Zotero collection). The novelty review found the closest prior and concurrent work to position against: Liu et al. 2025 (Stat Med) on positivity-driven estimand choice, Karim 2026 on plasmode simulations refuting but not certifying, Conover et al. 2025 (JAMIA) on gating unblinding on diagnostics, Nance et al. 2024 (Epidemiology), Shaw et al. 2025, Petersen et al. 2012, Li, Morgan and Zaslavsky 2018, Crump et al. 2009. Add them through sci-cite.

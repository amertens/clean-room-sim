# Rescue.Co application: re-execution on cleanTMLE 0.4.1

Decisions approved by Andrew on 2026-10-08. The scripts are to be written and run after the simulation pilot finishes (one R job at a time on this machine). The existing clean-room scripts (rescueCo/scripts/10_multiarm_build.R, 11_multiarm_design.R, 13_multiarm_estimation.R) run on cleanTMLE 0.3.0 and differ from the paper's Methods on each point below.

## Decisions

1. **Primary outcome:** poor function at the six-month follow-up interview, Glasgow Outcome Scale-Extended 4 or less, with death scored 1 (the protocol's primary outcome, on which the study was powered). Death within 180 days of arrival is secondary. Reason: the analyst had seen outcomes in the unblinded analysis, so a case study about prespecification cannot switch the protocol's primary outcome.
2. **Negative controls:** keep the five (urban household; cooking with wood or charcoal; hypertension; insulin-treated diabetes; HIV on antiretroviral therapy). Fit each without the covariate it is derived from (cooking fuel; chronic illness), using the new `nc_exclude` plan argument. Describe them as characteristics measured before transport that test residual imbalance, not as negative-control outcomes in the strict sense. Urban household comes from a near-constant variable and will probably be too rare to grade; report it as not graded rather than replacing it.
3. **Adjustment set:** primary set restricted to covariates fixed before transport. Exclude payment method (recorded at discharge) and the distance from the injury to the receiving hospital (depends on where transport took the patient). Arrival physiology (heart rate, blood pressure, Glasgow Coma Scale components, palpable pulse) is added in a sensitivity analysis.
4. **Missing outcomes:** complete case primary, with the design stage run on the complete-case cohort. Inverse probability of censoring weighting as a sensitivity analysis, with the probability of an observed outcome estimated by the plan's propensity-score super learner from treatment and the primary adjustment set (in the application script, through tmle's missingness argument).
5. **Unconfirmed Rescue.Co patients:** the treated patients without a confirmed Rescue.Co ID (being checked by Rescue.Co) stay in the primary analysis and are excluded in a sensitivity analysis.
6. **Package:** `analysis_plan()` refuses a control that is also an adjustment covariate and takes `nc_exclude` (branch nc-exclude).

## Resolved 2026-10-08

- Protocol date: 9 May 2023 (Andrew).
- The five control characteristics are in the hospital registry file (enrolment form, page 1: household area, cooking fuel, chronic illness checkboxes), recorded for essentially every registry patient, not in the six-month follow-up file. They describe the patient before the injury and are not restricted to patients reached at follow-up.
- Causal diagram: report that no sufficient adjustment set exists while phone access, awareness of the service, wealth and injury severity not captured by the injury severity score are unmeasured, as a limitation, and report the E-value for the primary estimate (cleanTMLE's summary() already computes it through EValue).

- Missing-measurement indicators for arrival heart rate, blood pressure and pulse (colliders in the causal diagram; treatment predicts them): out of the primary set; included only with arrival physiology in the sensitivity analysis, noted there as colliders (Andrew, 2026-10-08).

## Open items

- The quantitative bias analysis for deaths before arrival (rescueCo_analysis/reports/time-zero-issue.md), if adopted for the protocol analysis.

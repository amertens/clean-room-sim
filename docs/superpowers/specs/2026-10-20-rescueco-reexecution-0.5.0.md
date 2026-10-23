# Rescue.Co re-execution on cleanTMLE 0.5.0

Plan, 2026-10-20. Implements the decisions in `2026-10-18-outcome-blind-design-stage.md` (D3, D6, D9, D10, and "Case study and paper on 0.5.0", all approved by Andrew on 2026-10-20). Supersedes `2026-10-08-rescueco-reexecution-0.4.1.md` where they differ. Calls made without Andrew are marked **[call]** and listed at the end for him.

## 1. Scripts (new; the 0.2.0-era scripts 01 to 14 stay as the legacy record)

| Script | Who runs it, when | Reads | Writes |
|---|---|---|---|
| `20_steward_indicators.R` | The data steward, before the design stage | Registry disposition; follow-up records (linkage and vital-status fields) | `data/design/observation_indicators.rds`: id, disposition recorded, six-month status recorded. Nothing else. |
| `21_build_design_data.R` | The analyst | Registry fields fixed before transport, transport, cohort flags, arrival state, negative-control sources (whitelisted columns only; no follow-up file) | `data/design/design_data.rds`, the variable dictionary, the column manifest |
| `22_design_stage.R` | The analyst | Design data, indicators | For each lock: plan, lock, design, simulation, negative-control ladder, dossier (HTML); checkpointed per lock |
| `23_registration_bundle.R` | The analyst | Dossiers | The outcome-free rule for the confirmatory comparison, the population gate, lock hashes and plan summaries for the OSF registration |
| `24_steward_outcomes.R` | The data steward, only after approval (D9) | Follow-up and disposition | `data/outcomes/outcome_store.rds`; arm-blind item checks logged |
| `25_estimate.R` | The analyst, after `unblind()` | Locks, dossiers, outcome store | Estimates; the longitudinal TMLE sensitivity analysis (lmtp, outside the package) |
| `26_reconcile.R` | The analyst, last | This study's estimates and the unblinded registry | Reconciliation table |

Scripts 24 to 26 are written now and run only after approval. No script before 24 derives, reads or logs an outcome value. Patient-level files live under `rescueCo/data/` (git-ignored); only aggregate tables go under `rescueCo/results/`.

## 2. Design data

**Linkage.** The steward links the follow-up file as the unblinded analysis does: exact id match after the accepted typo corrections and reconciliation actions, read from the unblinded analysis's two id-only audit files, then de-duplication preferring a row with vital status. The analyst's build never opens the follow-up file.

**Treatment and arms**, following the unblinded analysis (inventory row 6). Rescue.Co: `ems_ambulance` = 1, or a Rescue.Co id with the checkbox blank. Other operator: `ems_ambulance` = 0. An ambulance arrival with neither is "operator unknown" and has no treatment value (it is not filled as other operator). Arrival by other means: transport mode recorded and not ambulance; mode 99 (unknown) has no arm. Walk-ins: mode label "Walked in". Rescue.co's cross-reference tier is kept (1, 2, 3, 5); tier 5 (the 14 without a dispatch record) is excluded from every lock; tier 3 (probable matches) is kept in the primary analysis and excluded in a sensitivity analysis (D10).

**Cohort flags.** Interfacility transfer and prior care elsewhere (code 99 read as no, as in both analyses); direct arrival = neither; enrolment date = the date of hospital arrival.

**Adjustment set: characteristics fixed before transport**, as the unblinded analysis defines them (its 75 columns from about 30 variables): age, sex, marital status, education, occupation, home ownership, household assets, cooking fuel, habitual alcohol and tobacco use, alcohol at injury, chronic illness, medical history, prior surgery, sub-county of injury, injury mechanism, intent, place and activity, penetrating wound, injury severity score, number of regions injured, head injury severity, months since the start of enrolment, hour of injury, weekend, strike period. Encoding follows the unblinded analysis: covariates missing for more than 40% are dropped; continuous covariates median-imputed with an indicator when missing for at least 1%; nominal levels under 2% pooled; checkbox items under 1% pooled per family; zero-variance, duplicate and near-constant (under 2% or over 98%) columns dropped. Each design-stage run anchors these filters on its own primary cohort **[call 5]**.

**Excluded from the design data** (recorded after transport or outcome-related): payment, distance to the hospital reached, hospital reached, arrival time and arrival-time covariates, dispatch and scene times, signs of life, scene care, operator name, disposition values, treatments, every follow-up item.

**Arrival state** (sensitivity analysis only): systolic and diastolic blood pressure, heart rate, pulse examined, the three Glasgow Coma Scale components, and the missingness indicators of blood pressure and heart rate. Declared in the dictionary as stand-ins, recorded and fixed after time zero, with the reason (the closest record of severity at the scene) and the sensitivity analysis (the primary lock, which leaves them out).

**Negative controls**: urban household, cooking with wood or charcoal, hypertension, insulin-treated diabetes, HIV on antiretroviral therapy (unchecked box read as absent). `nc_exclude`: cooking fuel columns for the cooking control; chronic-illness columns for the three chronic controls, plus the medical-history items that duplicate diabetes or HIV if the data's value labels confirm it **[call 6]**.

**Outcome-observation indicators (D3)**, made by the steward (script 20), passed to the design data as two columns only: disposition recorded (disposition present and not 99), and six-month status recorded, defined as in the unblinded analysis's longitudinal TMLE: 1 for a death at disposition (whose outcome is known) and for a patient whose follow-up gives vital status or any GOSE item; 0 otherwise. A 0 among patients with disposition recorded therefore means a patient who left the emergency department alive and was not reached, as the blinding ledger says **[call 2]**.

**Timing fallbacks.** Months since start, weekend and strike period use the injury date, falling back to the arrival date when the injury date is missing (same day for nearly all). The hour of injury does not fall back to the arrival hour, which depends on transport: a missing injury hour is median-imputed like any continuous covariate **[call 4]**.

## 3. Variable dictionary

Time zero: the transport decision (dispatch, or departure from the scene), measured as the time of injury because arrivals by other means have no dispatch time; follow-up counted from arrival, with the response that the time from injury to arrival is short relative to six months (I2, row 30). Eligibility: arrived alive at a participating hospital, assessed after time zero, with the response citing D6's registered bias analysis and restriction. Every registry field is recorded after time zero (the registry is completed at hospital); `fixed` separates them.

| Variables | Fixed | Note |
|---|---|---|
| Demographic, socioeconomic, habits, medical history, chronic illness, prior surgery | before | enrolment form |
| Injury mechanism, intent, place, activity, sub-county, alcohol at injury, timing covariates | before | injury form or derived from the injury date |
| Injury severity score, regions injured, penetrating wound | before | coded at hospital but describe the injury; the reason says so **[call 3]** |
| Head injury severity | before | from the injury-form checkboxes; the reason records the open question to co-authors whether they were graded from the arrival Glasgow Coma Scale **[call 3]** |
| Treatment | at | recorded at arrival, confirmed against dispatch records |
| Direct arrival (transfer or prior care) | after | restriction variable; reason: route, and care sought elsewhere, may follow the transport decision **[call 3]** |
| Arrival state | after | stand-in (sensitivity locks only) |
| Negative controls | before | enrolment form |
| Disposition recorded, six-month status recorded | after | outcome-observation indicators; reason: D3, 2026-10-18 |

## 4. Locks

Each comparison has its own lock (decision 5); every lock is for one outcome.

| Lock | Comparison | Cohort | Covariates |
|---|---|---|---|
| L1 | Rescue.Co vs all other transport (protocol) | all arrivals with known arm | before transport |
| L2 | Rescue.Co vs all other transport (protocol) | direct arrivals | before transport |
| L3 | Rescue.Co vs other ambulance operators | direct arrivals | before transport |
| L4 | Rescue.Co vs arrival by other means | direct arrivals | before transport |
| L5 | Any ambulance vs arrival by other means | direct arrivals | before transport |
| S1 | as L3 | direct arrivals | before transport + arrival state (stand-in) |
| S2 | as L4 | direct arrivals, walk-ins excluded | before transport |
| S3 | as L5 | direct arrivals, walk-ins excluded | before transport |
| S4 | as L3 | direct arrivals, tier 3 excluded | before transport |

All locks are built on the complete-case rows: disposition recorded and six-month status recorded (decision 3; the design stage evaluates the estimator run after unblinding). The share lost by arm, from the indicators, is reported for each lock before locking.

**Run 1, the protocol's check point:** L1 to L5 for the primary outcome on patients enrolled in the first four months (19 February to 18 June 2024). This shows what the design stage would have said when the protocol's own check passed. In mid-2024 no six-month follow-up existed for these patients, so run 1's locks use every patient with a known arm, without the complete-case restriction or the observation indicators. They are built from the final data cut restricted by enrolment date, which approximates the data as entered in 2024 **[call 11]**. **Run 2, the final cohort:** L1 to L5 and S1 to S4 for the primary outcome, then the confirmatory comparison's lock for death within 180 days.

## 5. Plans

Estimands in order ATE, trimmed ATE (band 0.05 to 0.95), ATT, ATO. Propensity score: the package default super learner (elastic net, MARS, tuned boosted trees, small neural net; one learner per family). Candidates: main-terms logistic regression and the default super learner, each with truncation 0.01 and 0.05. Cross-fitting K = 5, V = 5. Outcome surfaces: package defaults (linear and nonlinear prognostic index, heterogeneity 0 and 0.5, log odds ratio 0.5), with the baseline risk declared: 0.22 for the primary outcome (the protocol's expected comparator risk of severe disability); 0.05 for death within 180 days **[call 7]**. Repetitions: set by a timing pilot (section 6) before any lock is final.

**Tolerances (decision 6).** Primary outcome: the protocol powered the study for a relative risk of 0.75 at a comparator risk of 22%, a difference of 0.055; bias tolerance 0.0275 (half of it); half-width tolerance 0.055 **[call 8]**. Death within 180 days: the protocol gives no power calculation; bias tolerance 0.01 and half-width tolerance 0.03 **[call 7]**. Coverage 0.90. Negative-control null band: ± 0.25 times the median prevalence of the five controls in the lock's cohort (computed from design data), at least ± 0.02 **[call 9]**.

**Population gate (decision 7), applied by hand from each dossier's populations table:** an estimand below the ATE counts only if it keeps at least 80% of each arm (trimmed ATE) and shifts no named covariate (age, injury severity score, head injury severity, assault intent, injury at home; road-traffic injury has no column of its own in the design) by more than 0.25 standardized mean differences; the overlap-weighted effect counts only for L3, S1 and S4 **[call 10]**.

**Confirmatory comparison rule (decision 5), applied by script 23 from run 2's dossiers:** among L1 to L5, the first in declared order whose overlap grade is good or moderate, whose primary estimand passes the population gate keeping at least 80% of each arm, and whose negative-control verdict is not STOP. Expected: L3. Whatever it picks is reported.

## 6. Compute

Before building every lock, a pilot times two repetitions of `simulate_design()` on run 2's L3 (n about 1,600) and L1 (n about 8,000), with the plan above. If the projected time for all locks exceeds about 72 hours of machine time, repetitions are reduced to 100 with `max_reps` 300, and then the super-learner candidate is dropped from L1 and L2 only. Each choice is recorded with its reason before any lock is final; none uses an outcome. Long runs are launched detached, one R job at a time, checkpointed per lock.

## 7. After approval

Script 24 builds the outcome store with the unblinded analysis's scoring: GOSE from the follow-up items with the shop-and-travel rule and the rule for impossible "cannot say any words" answers; GOSE 1 for every death at disposition whatever the interview; GOSE 4 or below as the primary outcome; death within 180 days of arrival as the secondary outcome, with deaths dated after 180 days, outside hospital, released. Script 25 unblinds each lock, estimates with the selected candidate (complete case), and runs the longitudinal TMLE with lmtp (disposition recorded, death at disposition, six-month status recorded, outcome) as the sensitivity analysis, labelled as not evaluated in the design stage.

## 8. Calls for Andrew

1. **Arms and the unknown operator** follow the unblinded analysis: the 12 ambulance arrivals with no recorded operator, and the 8 with transport mode unknown, have no arm (the clean-room build had put them in the comparator).
2. **Six-month status recorded** is made by the steward with the unblinded analysis's definition (a death at disposition counts as recorded), so the analyst's build reads no death value and the ledger sentence holds.
3. **Dictionary timing:** injury severity score, regions injured, penetrating wound and head injury severity declared fixed before transport, each with a reason that says it is coded at hospital (head injury severity: the co-author question is recorded); direct arrival declared fixed after time zero, so it is flagged as a restriction on a characteristic that may follow the transport decision.
4. **Hour of injury** is not filled from the arrival hour.
5. **Encoding filters** anchor on each run's own primary cohort, so run 1's columns may differ from run 2's, as they would have in 2024.
6. **Medical-history items** that duplicate diabetes or HIV go into `nc_exclude` if the value labels confirm it.
7. **Death within 180 days:** baseline risk 0.05, bias tolerance 0.01, half-width tolerance 0.03, all judgments in the absence of a protocol power calculation for mortality.
8. **Primary outcome half-width tolerance** equal to the smallest important difference, 0.055.
9. **Negative-control null band** ± 0.25 times the controls' median prevalence, at least ± 0.02.
10. **Population gate thresholds:** 80% of each arm kept; 0.25 standardized mean differences on five named covariates; overlap-weighted effect allowed only for the comparisons between ambulance operators.
11. **Run 1 data:** the final data cut restricted to the first four months approximates the data entered by mid-2024, and run 1 has no complete-case restriction because no follow-up existed then.

## 9. Decisions on the steps after the design stage (Andrew, 2026-10-23)

1. **Locks estimated.** Script 25 attempts every lock in the approvals file. A lock that `unblind()` refuses (no feasible estimand, or a STOP negative-control verdict) is reported as refused and never overridden. The comparison the rule picks carries the confirmatory claim; the other comparisons are descriptive and the S locks are sensitivity analyses.
2. **The rejected candidate.** The candidate the simulation did not select is not estimated after unblinding; the two libraries are compared from the design stage's own results, without the outcome.
3. **Approval record.** Script 23 writes an approvals template with each lock's role and hashes; Andrew completes it as `rescueCo/registration/approvals_final.csv`. Scripts 24 and 25 check every approved design report against it. This adds to the governance decision of 20 October; since Andrew holds every role, it records the order of events rather than an independent control.
4. **Push.** `rescueCo/registration/` is committed and pushed to GitHub before script 24 runs, and script 24 checks this, so a dated copy of the hashes exists off this machine before any outcome is scored.
5. **Death within 180 days.** The confirmatory comparison's lock for this outcome (section 4, run 2) is built by script 22b from that comparison's final lock, changing only the outcome, the baseline risk (0.05) and the tolerances (bias 0.01, half-width 0.03), because script 22 cannot change while the design stage runs. If the outcome is missing on any of the lock's rows, it is not estimated until Andrew decides.
6. **Comparators.** Weighting and matching are not run in the case study: the package's weighting comparator uses a main-terms propensity score without truncation, so a difference from the TMLE would mix learner, truncation and estimator, and matching targets the ATT. The Methods are changed to match.
7. **Cut-off.** A sensitivity lock or the 180-day lock without a design report 7 days after the last of L1 to L5 has one, or finished after that, is reported as not run, and script 24 waits only for the rest.

The build is fingerprinted (script 22c) and scripts 22b and 25 stop if it changes. Script 25 writes to `rescueCo/data/estimates/final/` and `rescueCo/results/estimates_050/`, rehearses on synthetic outcomes under `rescueCo/data/rehearsal/`, and runs the longitudinal TMLE with settings fixed in the script: the comparison's rows before the complete-case restriction, saved and hashed by script 23; nodes disposition recorded, death at disposition, six-month status recorded and the primary outcome; a super learner ensemble with one learner type per family (elastic net, MARS, boosted trees, a small neural net); 5 cross-fitting folds and 5 inner folds; density ratios trimmed at their 0.999 quantile.

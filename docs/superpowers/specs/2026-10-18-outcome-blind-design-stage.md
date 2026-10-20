# A fuller outcome-blind design stage for Rescue.Co, and cleanTMLE 0.5.0

Started 2026-10-18. This spec records Andrew's decisions, each with its date and reason. It amends `2026-10-08-rescueco-reexecution-0.4.1.md` where a decision below says so.

Evidence: `rescueCo_analysis/reports/cleantmle_rescueco_review_2026-10-18.md` and the approved inventory `rescueCo/deviations_inventory.csv` / `.md`.

## Aim

The case study shows what a fuller outcome-blind design stage would have found and planned for before any outcome was seen. cleanTMLE is extended where it cannot yet do that. The paper's message is that a staged, outcome-blind workflow gets closer to valid causal inference in a hard observational study, with a plain account of what it cannot do. A feasible verdict says an estimator can recover the estimand if identification holds. The identification table, the negative controls and the bias analysis speak to whether it does.

## Approvals so far

- 2026-10-18. Inventory of departures approved (34 departures, 2 identification problems; O split into estimates available and estimates cited). Three points resting on records outside both repositories were dropped because they cannot be backed up.
- 2026-10-18. Uncommitted work kept and committed on `simulation-pilot` (5d4a144 to 8c96548), not pushed.

## Decisions

Each is asked one at a time, recommendation first. Status "pending" until answered.

### D1. Scope
Status: decided 2026-10-18 (Andrew; recommended option).

cleanTMLE 0.5.0 gets P1 (time zero, eligibility, variable timing), P2 (how the outcome is observed), P3 (several comparisons, restricted populations, who is unsupported) and P4 (gates and tolerances). The paper rescopes simulation claim 1 to its conditions (Phase 4), shows the outcome-observation check in the case study with one supplementary simulation table, and discusses selection after treatment with its prespecified bias analysis. P5, P6 and P7 wait. Because P5 waits, the bias analysis for deaths before arrival (D6) is prespecified and registered outside the package, with its parameters in the registered plan.

Reason: these four items cover the departures coded D or S that the inventory shows a design stage could raise (P1 to P3) and the tolerances and population shift that decide the ladder (P4); the rest can follow without changing what the case study shows.

### D2. Missing outcomes
Status: decided 2026-10-18 (Andrew; recommended option). Reverses decision 4 of the 2026-10-08 spec.

The plan declares how the outcome is observed. Single-stage censoring weights are refused when an event that forces the outcome to be observed (death at emergency-department disposition) precedes censoring. `simulate_design()` evaluates complete case and a cross-fitted longitudinal TMLE, fitted with the same fitter, and the design stage chooses between them before any outcome is seen, so the estimator run after unblinding is still the one evaluated. `methods.qmd:188-191` is to be rewritten (Phase 6).

Reason: why the longitudinal TMLE is needed is knowable without outcomes (deaths at disposition are always observed; only survivors can be lost), and complete case is weaker than what the unblinded analysis did; letting the simulation choose keeps the rule that the estimator applied is the one evaluated.

### D3. Design data
Status: decided 2026-10-18 (Andrew; recommended option).

The design data include two outcome-observation indicators: emergency-department disposition recorded, and six-month status recorded. Death at disposition stays in the outcome store and enters only at `unblind()`. The blinding ledger records what the six-month indicator reveals: a patient lost at six months is known to have left the emergency department alive.

Reason: the design stage needs the observed censoring shares to calibrate the simulated censoring and to report censoring positivity (P2), and these indicators carry almost no information about the effect; death at disposition would.

### D4. Comparisons
Status: decided 2026-10-18 (Andrew; recommended option).

The comparisons are declared in order in one lock, the protocol's comparison first, each with its reason. A prespecified rule that reads no outcome picks the comparison that carries the confirmatory claim, from the overlap grade after declared restrictions, the share of each arm retained, and the negative-control verdict. The expected pick is Rescue.Co against other ambulance operators among direct arrivals; if the rule picks otherwise, that is reported. The rule's exact form and thresholds are written in the P3/P4 package spec, fixed and tested on simulated data before the Rescue.Co design stage is rerun.

Reason: the unblinded analysis changed the confirmatory comparison while estimates existed (inventory row 1); a rule fixed in advance makes the same choice without them, and reports it even if it disagrees with the expected one.

### D5. Arrival state
Status: decided 2026-10-18 (Andrew; recommended option). Amends decision 3 of the 2026-10-08 spec, which added arrival physiology only in a sensitivity analysis.

Primary: the treatment model and the first regressions use covariates fixed before transport; arrival state (the Glasgow Coma Scale components, blood pressure, heart rate, pulse examined) is a node after transport, used only in the censoring models and the later regressions of the longitudinal TMLE. The effect of transport stays a total effect, while missing at random becomes more plausible. Sensitivity: the unblinded analysis's choice, arrival state as a stand-in for severity at the scene, declared under P1 with its reason. The dossier shows the outcome-free evidence for each: how well transport predicts arrival state given the before-transport set, and the balance of arrival consciousness after weighting.

Reason: arrival state is recorded after transport and may be affected by it, so adjusting for it as a confounder can remove part of the effect; used after treatment, it still informs who is lost to follow-up.

### D6. Deaths before arrival
Status: decided 2026-10-18 (Andrew; recommended option). Closes the open item of the 2026-10-08 spec.

The target-trial estimand is named: the effect of transport in injured patients for whom a transport choice was made. The registry identifies it only if transport does not change survival to arrival, and the paper says so; `methods.qmd:170`, which defines the population after treatment, is replaced (Phase 6). Before unblinding, two analyses are prespecified and registered outside the package (D1):
- a two-sided quantitative bias analysis, adding hypothetical deaths drawn from the severe tail of each arm, with its bias parameters locked and a plausible range taken from outside data; if no outside source gives a range, the tipping point is reported and the paper says that no range was available;
- a restriction to patients who would almost surely have arrived alive by any route, defined from design data before unblinding.

Reason: entry to the registry requires arriving alive, after the transport being compared (inventory I1); the unblinded analysis's bias analysis added deaths to one arm at the crude risk and was chosen with estimates available.

### D7. Tolerances
Status: decided 2026-10-18 (Andrew; recommended option).

- The bias tolerance for each outcome is at most half the smallest difference that would change a decision, taken from the protocol's power calculation (protocol section on sample size); the plan states the difference and its source.
- A tolerance on the expected CI half-width is required, not optional.
- Negative-control null bands are set relative to each control's prevalence rather than at ±0.02 absolute.

P4 adds these, and plan validation warns when a bias tolerance exceeds half the smallest important difference.

Reason: the uniform 0.02 tolerance is about the size of the mortality effect of interest, so a feasible verdict could allow bias as large as the effect; an absolute band means different things for a rare and a common control.

### D8. Population gate
Status: decided 2026-10-18 (Andrew; recommended option).

Any estimand below the ATE in the ladder (trimmed ATE, ATT, ATO) must retain a declared share of each arm and stay within a declared standardized shift on named covariates. The ATO counts only if the plan declares it acceptable. An estimand that fails the gate cannot be primary, whatever its feasibility. The shares and shift limits are fixed in the plan and tested on simulated data before the Rescue.Co design stage is rerun.

Reason: in 0.4.1 the ladder takes the first feasible estimand, and the ATO is almost always feasible, so the ladder can land on a population that has moved far from the one the question is about without any check.

### D9. Approval
Status: decided 2026-10-18 (Andrew; recommended option).

The dossier is approved before unblinding by an approver from outside the team if one can be found; otherwise the blinding ledger records that every approver had seen the unblinded results. The lock hash and the plan, including the D6 bias analysis, are registered on OSF Registries, time-stamped, before unblinding. Andrew makes the registration; nothing is posted on his behalf without his word.

Reason: every member of the team has seen the unblinded estimates, so an approval by the team cannot show independence from them; a registration outside the package and outside the team's control dates the lock.

### D10. Cohort
Status: decided 2026-10-18 (Andrew; recommended option). Amends decision 5 of the 2026-10-08 spec.

The design stage starts from the protocol's cohort. The unblinded analysis's restrictions are declared as rungs and evaluated without outcomes: direct arrivals (transfers and prior care removed); known operator; the 14 Rescue.Co-classified patients with no dispatch record; and, as a sensitivity analysis, walk-ins excluded. The 14 are excluded, as `methods.qmd:165` says and as decided on 2026-10-02. `methods.qmd:192` and decision 5 of the 2026-10-08 spec are rewritten to mean the probable matches (tier 3 of Rescue.co's cross-reference): kept in the primary analysis and excluded in a sensitivity analysis.

Reason: the restrictions were chosen in the unblinded analysis while estimates existed (inventory rows 1, 3, 5, 7); declaring each as a rung lets the design stage show, without outcomes, what each does to overlap and to the share of each arm retained.

### D11. Clean-up
Status: decided 2026-10-18 (Andrew; recommended option).

Andrew makes both changes himself in `rescueCo_analysis`, so the rule that this work never writes to that repository stays intact: archive `clean_room/` with a pointer to `clean-room-sim/rescueCo/`, and remove the `rm -rf clean_room/results/*` allowance from `.claude/settings.local.json` (line 25).

Reason: `clean_room/` is an unmaintained fork written for cleanTMLE 0.1.1 that nothing live reads, and the allowance pre-approves deleting its results.

## Decisions after P1 (Andrew, 2026-10-19)

- P1 merged into cleanTMLE main with nc-exclude (52b8c07); main is pushed to GitHub, and merged local branches are deleted.
- A refusal of a column recorded before the attribute it describes is fixed is not added. Reason: Andrew refused it; it is outside the approved P1 design.
- `design_report()` re-runs the dictionary checks, as `.check_nc_not_covariate()` is re-run, so a lock re-hashed by hand cannot carry a bad dictionary into the dossier. To be built in P2.
- Next item: P2 (how the outcome is observed).

## Release scope (Andrew, 2026-10-20)

The first public release of cleanTMLE is main as it stands, with item P1, as version 0.5.0. Items P2 to P7 and the features that would match causalRisk move to `cleanTMLE/WISHLIST.md`. This supersedes the package part of D1 (which put P1 to P4 in 0.5.0). Reason: the package was taking too long to share with collaborators and was becoming too complex to show cleanly why it is useful; a first release about the size of causalRisk, with a wishlist, comes first. The decisions that depended on P2 to P4 (D2, D4, D5, D7, D8, and the comparison rungs of D10) are reopened for the case study and the paper.

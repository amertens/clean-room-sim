# Rescue.Co design stage: results so far (26 October 2026)

Aggregate, outcome-free results from the design reports of the cleanTMLE 0.5.0 re-execution (spec `docs/superpowers/specs/2026-10-20-rescueco-reexecution-0.5.0.md`). No outcome has been joined; every number below comes from the design data, the propensity score and the outcome-blind simulation. The design reports and their tables stay under `rescueCo/results/reexec_050/`, which git does not carry. The analyst had seen the outcomes in the parallel unblinded analysis, so this is a re-execution with the outcome withheld from the software, not a blinded analysis by someone who never saw it.

## Status

| Run | Lock | Status |
|---|---|---|
| First four months | L3 | done (laptop, 21 October) |
| First four months | L1 | to run |
| Final | L3 | done (laptop, 22 October) |
| Final | L1, L2, L4, L5 | done (desktop, 24 and 25 October) |
| Final | L3, death within 180 days (script 22b) | done (desktop, 25 October) |
| Final | S1 | done (desktop, 26 October) |
| Final | S4 | running |
| Final | S2, S3 | to run |
| Final | registration bundle (script 23, final) | after the above; cut-off 1 November 2026, 13:10 PDT |

Script 23 (pick) chose **L3**, Rescue.Co against other ambulance operators, as the confirmatory comparison, by the rule in the plan (`rescueCo/registration/confirmatory.csv`), and saved and hashed the longitudinal TMLE's cohort (`lmtp_cohort.csv`). The build is fingerprinted in `build_fingerprint.csv`; scripts 22b and 25 stop if it changes.

The finished locks of the laptop run rebuilt with the same hashes on the desktop, and final L1, run on 3 workers after a failed first attempt on 7, reproduced that attempt's lock hash and simulation metrics exactly.

## Results by lock

Primary outcome: GOSE 4 or below at six months (bias tolerance 0.0275, half-width tolerance 0.055). Death within 180 days: bias 0.01, half-width 0.03. Coverage 0.90 throughout. Overlap: share of patients with a propensity score outside 0.05 to 0.95.

| Lock | Patients | Overlap | Feasible estimands | Primary (candidate) | Population gate | Negative controls |
|---|---|---|---|---|---|---|
| First 4 months L3 | 226 | good (0%) | none | none | | FLAG |
| Final L3 | 1,204 | good (0.2%) | ATE, trimmed ATE, ATO | ATE (main terms) | passes | FLAG |
| Final L3, death within 180 days | 1,204 | good (0.2%) | ATE, trimmed ATE, ATO | ATE (super learner) | passes | FLAG |
| Final S1 (L3 + arrival state) | 1,204 | good (0%) | all four | ATE (main terms) | passes | FLAG |
| Final L1 (protocol) | 6,494 | poor (42%) | trimmed ATE, ATT, ATO | trimmed ATE (main terms) | fails: keeps 54% of comparators | FLAG |
| Final L2 | 4,207 | poor (30%) | none | none | | |
| Final L4 | 3,739 | poor (30%) | ATT | ATT (main terms) | fails: 0.40 SD shift | FLAG |
| Final L5 | 4,207 | moderate (11%) | ATT | ATT (main terms) | fails: 0.33 SD shift | FLAG |

L1, L2 and L4 ran with the main-terms candidate only and 100 repetitions without extension, and L5 with both candidates without extension (compute decisions of 23 October); cells left unresolved there are reported as unresolved.

## What the results show

1. **Poor overlap: the ATE fails, and the fallback estimands describe other populations.** In L1, L2, L4 and L5 the ATE misses the bias or coverage tolerance (coverage as low as 0.29 in L1). Estimands for the overlap population are often feasible, but the population gate rejects every one: they drop about half of a comparison arm or shift a named covariate by 0.33 to 0.44 SD. Only the comparison between ambulance operators can carry a claim.
2. **The four-month comparison fails on precision, not positivity.** No patient lay outside the overlap band, and the super learner met bias and coverage for every estimand on every surface; every estimand failed because its expected CI half-width (0.12 to 0.15) exceeded the smallest important difference (0.055). Main-terms TMLE was badly undercovered at this size (coverage 0.53 to 0.57).
3. **The simulation chooses the estimator, and the choice depends on the outcome.** For the primary outcome it selected main-terms TMLE in every lock, because the super learner candidate undercovered (0.86 to 0.88 in final L3). For death within 180 days, on the same patients and the same propensity score, it selected the super learner, because main-terms TMLE missed the half-width tolerance.
4. **Several verdicts sit at a tolerance.** Final L3's ATT failed on half-width by 0.0001 (0.0551 against 0.055); the 180-day ATE passed at 0.0298 against 0.03; L5's trimmed ATE failed on one surface's coverage (0.85) in a cell left unresolved because L5 had no extension. Final L3's decision rests on 12 cells unresolved at 300 repetitions.
5. **Negative controls read FLAG in every lock.** In the L3-sized cohorts no control's interval fits inside the ±0.02 band, so GO is out of reach. In the larger cohorts some controls read inside, but each lock had at least one inconclusive control, and a domain with one inconclusive control is inconclusive. No control read outside, so no lock got STOP. Insulin-treated diabetes cannot be estimated in any lock; urban household is near-constant (96.5% urban) and cannot be estimated at four months.

## The simulation's verdicts against balance diagnostics

A conventional check reads balance after weighting: an estimand passes when the largest absolute standardized mean difference after weighting is at most 0.1 (Austin 2009). Read from the same design reports (weights from the package's cross-fitted propensity score; matching was not run), against the simulation's verdicts for the primary outcome:

| Lock | ATE | Trimmed ATE | ATT | ATO |
|---|---|---|---|---|
| First 4 months L3 | 0.40 / infeasible | 0.40 / infeasible | 0.31 / infeasible | 0.36 / infeasible |
| Final L3 | 0.12 / **feasible** | 0.11 / **feasible** | 0.12 / infeasible | 0.103 / **feasible** |
| L1 | 0.21 / infeasible | 0.12 / **feasible** | 0.06 / feasible | 0.07 / feasible |
| L2 | 0.38 / infeasible | 0.101 / infeasible | 0.06 / **infeasible** | 0.05 / **infeasible** |
| L4 | 0.97 / infeasible | 0.12 / infeasible | 0.09 / feasible | 0.07 / **infeasible** |
| L5 | 0.17 / infeasible | 0.07 / **infeasible** | 0.06 / feasible | 0.04 / **infeasible** |

Each cell: largest weighted SMD / simulation verdict; bold where the two disagree. They agree on 15 of 24 estimand-comparison pairs. The simulation is stricter in 5 (an estimand balanced after weighting that TMLE could not estimate within tolerance; 2 of these rest on resolved cells) and more permissive in 4 (TMLE within tolerance despite a weighted SMD just above 0.1). The primary estimand chosen down the same order would differ in 4 of 6 comparisons at the 0.1 threshold, and in 3 of 6 at 0.15.

This shows the two rules decide differently; it cannot show which is right, because the case study has no true outcome model. The confirmatory choice (L3) is the same under either. The simulation study's comparators for claim 1 are overlap diagnostics, not balance after weighting.

## Next

S4, S2, S3 and the first four months' L1, then script 23 (final) writes the bundle and the approvals template. After approval: commit and push `rescueCo/registration/`; script 24 (steward outcomes) and script 25 (estimation) on the desktop, whose build the fingerprint records, with `R_LIBS` set to the design library (`%LOCALAPPDATA%\R\rescueCo-design-lib-4.4.2`) and no package updates in the R 4.4.2 library until script 25 has run. Script 24 needs the registry and follow-up files, which are not on the desktop.

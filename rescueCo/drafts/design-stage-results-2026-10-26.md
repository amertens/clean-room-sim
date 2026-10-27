# Rescue.Co design stage: results (final, 27 October 2026)

Aggregate, outcome-free results from the design reports of the cleanTMLE 0.5.0 re-execution (spec `docs/superpowers/specs/2026-10-20-rescueco-reexecution-0.5.0.md`). No outcome has been joined; every number below comes from the design data, the propensity score and the outcome-blind simulation. The design reports and their tables stay under `rescueCo/results/reexec_050/`, which git does not carry. The analyst had seen the outcomes in the parallel unblinded analysis, so this is a re-execution with the outcome withheld from the software, not a blinded analysis by someone who never saw it.

## Status

| Run | Lock | Status |
|---|---|---|
| First four months | L3 | done (laptop, 21 October) |
| First four months | L1 | done (desktop, 27 October) |
| Final | L3 | done (laptop, 22 October) |
| Final | L1, L2, L4, L5 | done (desktop, 24 and 25 October) |
| Final | L3, death within 180 days (script 22b) | done (desktop, 25 October) |
| Final | S1, S4, S2, S3 | done (desktop, 26 and 27 October) |
| Final | registration bundle (script 23, final) | written 27 October, 00:52 PDT, before the cut-off (1 November 2026, 13:10 PDT) |

Script 23 (pick) chose **L3**, Rescue.Co against other ambulance operators, as the confirmatory comparison, by the rule in the plan (`rescueCo/registration/confirmatory.csv`), and saved and hashed the longitudinal TMLE's cohort (`lmtp_cohort.csv`); script 23 (final) confirmed the pick. The build is fingerprinted in `build_fingerprint.csv`; scripts 22b and 25 stop if it changes.

The finished locks of the laptop run rebuilt with the same hashes on the desktop, and final L1, run on 3 workers after a failed first attempt on 7, reproduced that attempt's lock hash and simulation metrics exactly. Every lock in the bundle passes `verify_lock()`, and the bundle's lock and design-report hashes match the saved design reports. The bundle covers the final run; the four-month locks are not in it.

## Results by lock

Primary outcome: GOSE 4 or below at six months (bias tolerance 0.0275, half-width tolerance 0.055). Death within 180 days: bias 0.01, half-width 0.03. Coverage 0.90 throughout. Overlap: share of patients with a propensity score outside 0.05 to 0.95.

| Lock | Patients | Overlap | Feasible estimands | Primary (candidate) | Population gate | Negative controls |
|---|---|---|---|---|---|---|
| First 4 months L3 | 226 | good (0%) | none | none | | FLAG |
| First 4 months L1 (protocol) | 1,661 | poor (49%) | none | none | | |
| Final L3 | 1,204 | good (0.2%) | ATE, trimmed ATE, ATO | ATE (main terms) | passes | FLAG |
| Final L3, death within 180 days | 1,204 | good (0.2%) | ATE, trimmed ATE, ATO | ATE (super learner) | passes | FLAG |
| Final S1 (L3 + arrival state) | 1,204 | good (0%) | all four | ATE (main terms) | passes | FLAG |
| Final S4 (L3 without probable tier-3 dispatch matches) | 1,195 | good (0%) | all four | ATE (main terms) | passes | FLAG |
| Final L1 (protocol) | 6,494 | poor (42%) | trimmed ATE, ATT, ATO | trimmed ATE (main terms) | fails: keeps 54% of comparators; 0.44 SD shift | FLAG |
| Final L2 | 4,207 | poor (30%) | none | none | | |
| Final L4 | 3,739 | poor (30%) | ATT | ATT (main terms) | fails: 0.40 SD shift | FLAG |
| Final S2 (L4 without walk-ins) | 3,478 | poor (26%) | trimmed ATE, ATT | trimmed ATE (main terms) | fails: keeps 67% of comparators; 0.39 SD shift | FLAG |
| Final L5 | 4,207 | moderate (11%) | ATT | ATT (main terms) | fails: 0.33 SD shift | FLAG |
| Final S3 (L5 without walk-ins) | 3,946 | moderate (8%) | ATT | ATT (main terms) | fails: 0.31 SD shift | FLAG |

Four-month L1, final L1, L2, L4, S2 and S3 ran with the main-terms candidate only and 100 repetitions without extension, and L5 with both candidates without extension (compute decisions of 23 October); cells left unresolved there are reported as unresolved (8 in S2, 10 in S3). Final L3 has 12 unresolved cells at 300 repetitions, S4 has 3.

## What the results show

1. **Poor overlap: the ATE fails, and the fallback estimands describe other populations.** In L1, L2, L4, L5, S2 and S3 the ATE misses the bias or coverage tolerance (coverage as low as 0.29 in L1). Estimands for the overlap population are often feasible, but the population gate rejects every one: they drop a third to a half of the comparators or shift a named covariate by 0.31 to 0.44 SD. Leaving out walk-ins (S2, S3) narrows the overlap problem a little (26% and 8% outside the band, against 30% and 11%) and makes the trimmed ATE feasible in S2, but no estimand passes the gate. Only the comparison between ambulance operators, and its sensitivity versions, can carry a claim.
2. **At four months, the comparison between ambulance operators fails on precision and the protocol's comparison on both.** In L3 no patient lay outside the overlap band, and the super learner met bias and coverage for every estimand on every surface; every estimand failed because its expected CI half-width (0.12 to 0.15) exceeded the smallest important difference (0.055). Main-terms TMLE was badly undercovered at this size (coverage 0.53 to 0.57). In the protocol's comparison (L1), 49% of patients lay outside the band; every estimand failed the half-width tolerance (0.085 to 0.115) and also failed bias or coverage on at least one surface.
3. **The simulation chooses the estimator, and the choice depends on the outcome.** For the primary outcome it selected main-terms TMLE for every lock's primary estimand, because the super learner candidate undercovered where both ran (0.86 to 0.88 in final L3; 0.86 to 0.89 in S1 and S4). For death within 180 days, on the same patients and the same propensity score, it selected the super learner, because main-terms TMLE missed the half-width tolerance. In S1 and S4 the super learner carried the ATT: in S4 because main-terms TMLE missed the half-width tolerance (0.0559), in S1 because it undercovered on one surface (0.89).
4. **Several verdicts sit at a tolerance.** Final L3's ATT failed on half-width by 0.0001 (0.0551 against 0.055); the 180-day ATE passed at 0.0298 against 0.03; S4's ATE passed with coverage 0.903 against 0.90. L5's trimmed ATE failed on one surface's coverage (0.85) and S3's on one surface's coverage (0.87), each in a cell left unresolved because the lock had no extension; S3's trimmed ATE would have passed the population gate (it keeps 99% and 89% of the arms and shifts no named covariate by more than 0.16 SD).
5. **Negative controls read FLAG in every lock where they were judged** (all but L2 and four-month L1). In the L3-sized cohorts no control's interval fits inside the ±0.02 band, so GO is out of reach. In the larger cohorts some controls read inside, but each lock had at least one inconclusive control, and a domain with one inconclusive control is inconclusive. No control read outside, so no lock got STOP. In S3, cooking with wood or charcoal read 0.018 (0.004 to 0.033): the interval excludes zero but crosses the band edge, so it reads inconclusive, not null. Insulin-treated diabetes cannot be estimated in any lock; urban household is near-constant (96.5% urban) and cannot be estimated at four months.

## The simulation's verdicts against balance diagnostics

A conventional check reads balance after weighting: an estimand passes when the largest absolute standardized mean difference after weighting is at most 0.1 (Austin 2009). Read from the same design reports (weights from the package's cross-fitted propensity score; matching was not run), against the simulation's verdicts for the primary outcome:

| Lock | ATE | Trimmed ATE | ATT | ATO |
|---|---|---|---|---|
| First 4 months L3 | 0.40 / infeasible | 0.40 / infeasible | 0.31 / infeasible | 0.36 / infeasible |
| First 4 months L1 | 0.47 / infeasible | 0.36 / infeasible | 0.18 / infeasible | 0.18 / infeasible |
| Final L3 | 0.12 / **feasible** | 0.11 / **feasible** | 0.12 / infeasible | 0.103 / **feasible** |
| S1 | 0.13 / **feasible** | 0.13 / **feasible** | 0.12 / **feasible** | 0.105 / **feasible** |
| S4 | 0.108 / **feasible** | 0.108 / **feasible** | 0.12 / **feasible** | 0.101 / **feasible** |
| L1 | 0.21 / infeasible | 0.12 / **feasible** | 0.06 / feasible | 0.07 / feasible |
| L2 | 0.38 / infeasible | 0.101 / infeasible | 0.06 / **infeasible** | 0.05 / **infeasible** |
| L4 | 0.97 / infeasible | 0.12 / infeasible | 0.09 / feasible | 0.07 / **infeasible** |
| S2 | 0.36 / infeasible | 0.105 / **feasible** | 0.08 / feasible | 0.06 / **infeasible** |
| L5 | 0.17 / infeasible | 0.07 / **infeasible** | 0.06 / feasible | 0.04 / **infeasible** |
| S3 | 0.108 / infeasible | 0.08 / **infeasible** | 0.08 / feasible | 0.05 / **infeasible** |

Each cell: largest weighted SMD / simulation verdict; bold where the two disagree. They agree on 23 of 44 estimand-comparison pairs. The simulation is stricter in 8 (an estimand balanced after weighting that TMLE could not estimate within tolerance; 2 of these rest on resolved cells) and more permissive in 13 (TMLE within tolerance although weighting leaves a largest SMD of 0.101 to 0.130; 11 of the 13 are in L3 and its sensitivity versions). The primary estimand chosen down the same order would differ in 8 of 11 comparisons at the 0.1 threshold, and in 4 of 11 at 0.15.

The confirmatory comparison depends on the threshold. Applied in place of the simulation in the confirmatory rule, a balance check at 0.1 finds no estimand for L3, whose largest weighted SMDs are 0.103 to 0.122, and the rule moves on to L5 (any ambulance against arrival by other means), whose trimmed ATE is balanced (0.07) and passes the population gate (it keeps 99% and 85% of the arms and shifts no named covariate by more than 0.20 SD). L5's negative controls were run on its ATT, where none read outside; they were not run on its trimmed ATE. At 0.15 the balance check keeps L3 and its ATE, as the simulation does.

This shows the two rules decide differently, including on which comparison carries the claim; it cannot show which is right, because the case study has no true outcome model. Every disagreement where the simulation is more permissive sits within 0.03 of the 0.1 threshold. The simulation study's comparators for claim 1 are overlap diagnostics, not balance after weighting.

## Next

Andrew reviews each design report and completes `approvals_final_template.csv` as `approvals_final.csv`; then commit and push `rescueCo/registration/` before script 24 builds the outcome store (decision 4). Script 24 (steward outcomes) and script 25 (estimation) run on the desktop, whose build the fingerprint records, with `R_LIBS` set to the design library (`%LOCALAPPDATA%\R\rescueCo-design-lib-4.4.2`) and no package updates in the R 4.4.2 library until script 25 has run. Script 24 needs the registry and follow-up files, which are not on the desktop.

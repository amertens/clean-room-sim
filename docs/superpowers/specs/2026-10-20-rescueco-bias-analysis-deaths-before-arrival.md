# Rescue.Co: prespecified bias analysis for deaths before arrival

Approved by Andrew, 2026-10-20 (calls 1 to 10, the sources below and the Methods wording). This is decision D6 of `2026-10-18-outcome-blind-design-stage.md`, made concrete for the re-execution on cleanTMLE 0.5.0 (`2026-10-20-rescueco-reexecution-0.5.0.md`). Once approved, it is part of the analysis plan, written before unblinding: "registered" in D6 now means written into the plan document before unblinding ("Governance in the demonstration", 2026-10-20). The commit that fixes this document is recorded in `rescueCo/BLINDING_LEDGER.md` before script 24 builds the outcome store, and its hash is reported in the supplement with the lock and dossier hashes **[call 10]**. Choices made without Andrew are marked **[call]** and listed in section 10.

Nothing here reads an outcome. Every rule and every number fixed before unblinding comes from design data (arm, cohort flags, injury severity score, head injury severity, penetrating wound) or from published sources found by logged searches (`writing/lit/dba_*`).

## 1. What the registry can and cannot identify

Time zero is the transport decision, and the registry enrols a patient only on arrival alive at hospital, after the transport being compared. Write $S(a) = 1$ if a patient would arrive alive with transport $a$, where $a = 1$ is the first-named arm of the comparison (Rescue.Co, or any ambulance in L5) and $a = 0$ the comparator. Patients for whom a transport choice was made fall into four groups:

| Group | $S(1)$ | $S(0)$ | In the registry's arm 1 | In the registry's arm 0 |
|---|---|---|---|---|
| Arrive alive by either transport | 1 | 1 | yes | yes |
| Die before arrival by either transport | 0 | 0 | no | no |
| Arrive alive only with arm 1 | 1 | 0 | yes | no (died on the way) |
| Arrive alive only with arm 0 | 0 | 1 | no (died on the way) | yes |

If transport does not change survival to arrival, the last two groups are empty and, with exchangeability given the covariates fixed before transport, the registry estimates the effect among patients who would arrive alive by either transport. If they are not empty, each registry arm lacks the deaths of the group that only the other transport would have saved, and those are the most severely injured patients of that arm.

Two consequences shape the analysis.

- Patients who would die before arrival by either transport have the outcome (death, and GOSE 1) under both transports, so they carry no effect. On the risk-difference scale, the effect among all patients for whom a transport choice was made equals the effect among the remaining patients times their share. The registry cannot estimate that share, but it is close to 1 when deaths before arrival are rare among patients who were transported. The bias analysis therefore targets patients who would arrive alive with at least one of the two transports, and the paper should say so (section 9, item 2).
- The bias analysis supplies the missing deaths of the last two groups, one arm at a time. The restriction (section 6) removes the patients for whom those groups could plausibly be non-empty.

## 2. Bias parameters

**Parameter.** For each arm $a$, $\delta_a$ is the number of patients transported by arm $a$ who died before arrival but would have arrived alive with the other transport, per 1,000 patients of arm $a$ who arrived alive. This is the scale on which outside data report deaths before arrival (patients brought in dead per 1,000 patients brought to hospital), so it is the scale on which the range is locked. **[call 1]**

**Where the deaths are placed.** The hypothetical deaths are drawn from the arm's severe tail, $T = 1$ when the injury severity score is recorded and 16 or more (major trauma), or the head injury is graded severe (design-data columns `eiss_ais` with its missing indicator, and `head_injury_severity` = 3). Deaths before arrival concentrate in head injuries and multisystem injuries in the closest outside data (head injury OR 11.8 for death in transit, Reid 2018; head, chest and abdominal injuries RR about 1.3 against extremity injuries, Purcell 2020). If the severe tail holds fewer than 20 patients of either arm in a lock's complete-case cohort, it is widened to an injury severity score of 9 or more or a moderate or severe head injury, and the widening is recorded in the plan before unblinding. **[call 2]**

**From $\delta_a$ to a share of the severe tail.** With $N_a$ the arm's arrivals in the lock's cohort and $M_a$ those of them with $T = 1$, both counted before the complete-case restriction because deaths are never lost to follow-up, the added deaths are $k_a = \delta_a N_a / 1000$ and their share of the arm's severely injured patients is

$$q_a = \frac{k_a}{M_a + k_a}.$$

$N_a$ and $M_a$ are design data, so script 23 computes $k_a$ and $q_a$ for every lock and every grid value before unblinding, and the plan carries them as numbers (an appendix table written before unblinding).

**Each arm in turn.** The analysis sets $\delta_0 = 0$ while varying $\delta_1$, then the reverse. The adjustment is additive in the two arms to first order (section 3), so any combination can be read off the two one-arm results; a joint grid is not run.

## 3. Bias-adjusted estimate and its interval

**Estimand.** For arm $a$, the bias-adjusted risk at covariates $X$ is $(1 - q_a T)\,Q_a(X) + q_a T$, where $Q_a(X)$ is the risk among arm-$a$ arrivals with covariates $X$. In words: a share $q_a$ of the severely injured patients sent by arm $a$ died before arrival, and the rest have the risk seen among arm-$a$ arrivals like them. The bias-adjusted effect standardizes these risks over the same population as the lock's primary estimand (ATE, trimmed ATE, ATT or ATO, whichever the design stage selected).

**Estimator.** The bias-adjusted risk is the expected value of a transformed outcome,

$$\tilde Y_i = Y_i + q_{A_i} T_i (1 - Y_i),$$

whose regression on treatment and covariates is $\tilde Q_a(X) = Q_a(X) + q_a T (1 - Q_a(X))$. So the bias-adjusted effect is the lock's primary estimand with $\tilde Y$ in place of $Y$, and it is estimated by the primary analysis's own targeting step with three inputs changed: $\tilde Y$ for $Y$, and $\tilde Q_1$, $\tilde Q_0$ for the initial outcome predictions. The propensity score and the outcome regression are not refit on $\tilde Y$; they are the selected candidate's cross-fitted predictions, transformed. At $\delta_a = 0$ the transformation is the identity and the primary estimate is returned exactly. $\tilde Y$ and $\tilde Q_a$ are linear in $(q_1, q_0)$, so the plug-in estimate is additive in the two arms and the targeting step changes that only slightly; this is why the arms are run one at a time (section 2).

**Interval.** The 95% CI at each $\delta$ comes from the influence curve of that targeting step, recomputed at each grid value, not shifted. Three properties of the added deaths set its width.

- They are not observations, so they do not add to $n$. (Adding them as records would count them as observed patients and narrow the interval.)
- They are certain events, so they add no outcome variance. For a severely injured patient of the receiving arm, the residual becomes $(1 - q_a)(Y - Q_a(X))$, which shrinks the residual variance in the severe tail by the factor $(1 - q_a)^2$.
- They raise the arm's predicted risk in the severe tail, which changes the plug-in part of the influence curve.

The bias parameter is fixed, so it contributes no variance of its own; uncertainty about it is shown by the plausible range (section 4), not folded into the interval. **[call 8]**

**Why this replaces the unblinded analysis's version.** Scripts 61 and 69 of `rescueCo_analysis` added $k$ deaths to one arm at its crude risk $r$, so the arm's risk moved by $k(1 - r)/(N + k)$, and shifted the adjusted estimate and both CI limits by that amount. Placed in the severe tail, where the risk $m_T$ among arrivals is already high, the same deaths move the arm's risk by about $(k/N)(1 - m_T)$, which is less, so the earlier version overstated the effect of each added death. It also kept the interval's width fixed and added deaths to one arm only.

**Approximation.** The standardization is over the arrived cohort's covariates, while the population that includes the added deaths has slightly more severely injured patients ($k_a$ more, at most 1.8% of the arm within the locked ranges). The difference is of order $q_a$ times the severe tail's share times the gap between the severe tail's effect and the overall effect. The simulation check (section 8) measures it.

**Implementation (script `25b_bias_deaths_before_arrival.R`, after script 25).** For each lock in scope, rebuild $X$, $A$ and the complete-case $Y$ as `estimate_effect()` does; make the folds with the dossier's $K$ and fold seed; refit the selected candidate's propensity score (seed = fold seed) and outcome regression (seed = fold seed + 1) with the package's internal fitters, exactly as `.fit_and_target()` does; check that the targeting step on $Y$ returns the package's primary estimate to within $10^{-6}$ (if a learner is not reproducible, report both estimates and use the refit throughout, so the $\delta = 0$ row is the refit's); then run the targeting step on $\tilde Y$, $\tilde Q_1$, $\tilde Q_0$ at every grid value. Cost: one more nuisance fit per lock. The script reads only the locks, dossiers, outcome store and the prespecified table of $q_a$ values; it writes aggregate tables only.

## 4. Plausible range, with sources

No source reports the quantity directly: deaths in transit among patients for whom a transport choice was made, by transport, in Nairobi. The range is built from the closest outside data, all found by the logged searches in `writing/lit/dba_share`, `dba_transport`, `dba_kenya` and `dba_mode`.

| Source | Setting | What it gives | Used for |
|---|---|---|---|
| Reid et al. 2018, World J Surg, doi:10.1007/s00268-018-4577-y (in bib: `reid2018`) | Kamuzu Central Hospital, Lilongwe, Malawi, 2008 to 2015; no prehospital system | 701 of 85,806 trauma patients presenting died in transit (0.8%); head injury and gunshot wounds most associated | total rate (about 8 per 1,000); severe-tail rule |
| Purcell et al. 2020, World J Surg, doi:10.1007/s00268-020-05470-w (in bib: `purcell2020`) | same hospital, adults, 2008 to 2018 | 1,130 deaths before arrival among 131,020 adults presenting (about 9 per 1,000), 56% of the hospital's trauma deaths; police transport RR 1.97 (1.52 to 2.55) against minibus, adjusted | total rate; relative risk between two kinds of transport |
| Mulima et al. 2021, Afr J Emerg Med, doi:10.1016/j.afjem.2021.03.011 (in bib: `mulima2021`) | same hospital, all ages, 2008 to 2019 | 1,889 of 106,198 trauma records brought in dead (1.8%); mostly carried by police (58%) or private vehicles (27%) | upper total rate (about 18 per 1,000) |
| Henry and Reingold 2012, J Trauma Acute Care Surg, doi:10.1097/TA.0b013e31824bde1e (in bib: `henry2012`) | meta-analysis of prehospital trauma systems in low- and middle-income countries | RR of death 0.75 (0.66 to 0.85); urban 0.79 | lower anchor for the effect of an ambulance |
| Gathecha et al. 2017, BMC Public Health, doi:10.1186/s12889-016-3958-0 (in bib: `gathecha2017`) | Nairobi death certificates, January to June 2014 | 51.4% of injury deaths occurred before hospital | that Nairobi's share of injury deaths before hospital resembles Lilongwe's |
| Saidi and Oduor 2013, J Forensic Leg Med, doi:10.1016/j.jflm.2013.04.016 (in bib: `saidi2013`) | Nairobi city mortuary autopsies, 2009 to 2010 | deaths before hospital more often intentional, at night, and by gunshot than hospital deaths | context for the severe-tail rule |
| Mock et al. 1998, J Trauma, doi:10.1097/00005373-199805000-00011 (in bib: `mock1998`) | Kumasi, Monterrey, Seattle | among seriously injured patients (injury severity score 9 or more, or dead), 51%, 40% and 21% died in the field with no, basic and advanced EMS | direction of the effect of EMS; the threshold of the restriction (section 6) |

Not used: the two US analyses of patients dead on arrival by injury severity score (Elkbuli et al. 2021 and 2022), which would have anchored the rate within the severe tail, because the 2022 paper carries an expression of concern (Am Surg 2025) and both come from the same group; a cross-sectional study from Addis Ababa (low-quality venue); conference abstracts.

**Construction.** Among patients brought to the Lilongwe referral hospital, where nearly all arrive without an ambulance, deaths before arrival were about 8 to 18 per 1,000 patients. If the transport with the higher death rate before arrival has rate $r_w$ and its relative risk against the other transport is $\rho$, the excess is $r_w(1 - 1/\rho)$. At $r_w = 18$ per 1,000:

- $\rho = 1/0.75$, the effect of a prehospital system on all trauma deaths (Henry and Reingold; a lower anchor, since a prehospital system acts mostly on deaths before arrival): excess 4.5 per 1,000;
- $\rho = 2$, the adjusted relative risk between police and minibus transport (Purcell et al.): excess 9 per 1,000;
- $\rho \to \infty$, every death before arrival in that arm prevented by the other transport: excess 18 per 1,000.

**Locked ranges [call 3, call 4]:**

| Comparisons | Range of $\delta$ for the arm receiving deaths | Marked values |
|---|---|---|
| Ambulance against arrival by other means (L4, L5), and Rescue.Co against all other transport (L1, L2), whose comparator is mostly arrival by other means | 0 to 18 per 1,000, in each direction | 4.5, 9, 18 |
| Rescue.Co against other ambulance operators (L3) | 0 to 9 per 1,000, in each direction | 4.5, 9 |

The death-within-180-days lock takes the range of its comparison.

The range between ambulance operators stops at the twofold relative risk at the highest total rate, because both arms travel with a crew and a vehicle and no source compares operators. The same range is used in both directions because the outside data favor the ambulance in low- and middle-income settings (Henry and Reingold; Mock et al.) but come from settings without Nairobi's traffic and ambulance response times, so they do not exclude the other direction.

If Andrew rejects the construction, D6's fallback applies: the tipping points are reported and the paper says that no range was available.

## 5. Grid, tipping points and the prespecified reading

**Grid.** $\delta$ from 0 to 40 per 1,000 in steps of 0.5, for each arm in turn, each lock in scope and each outcome. If a tipping point is not reached by 40, the grid continues to 100 in steps of 2; beyond that, or when $q_a$ would exceed 0.9, it is reported as beyond the grid.

**Tipping points.** For each arm in turn: (i) the smallest $\delta$ at which the bias-adjusted estimate reaches no difference; (ii) the smallest $\delta$ at which the 95% CI changes between excluding and including no difference, in either direction. Each is reported per 1,000 arrivals, as $q_a$ (share of the arm's severely injured patients), and as $k_a$ (number of deaths). Adding deaths to the arm with the lower risk moves the estimate toward no difference; adding them to the other arm moves it away, and the tipping point there is the one at which a null result would become a difference.

**Reading, fixed now.** The confirmatory conclusion is called robust to deaths before arrival if no tipping point of kind (ii), for either arm and the primary outcome, falls within the locked range for that comparison. Otherwise the paper says that the conclusion could be changed by deaths before arrival within the plausible range, and names the arm and the value.

**Displays.** Main text: one figure for the confirmatory comparison, the bias-adjusted risk difference and its 95% CI against $\delta$, one panel per arm receiving deaths, rows for the primary outcome and death within 180 days, with the locked range shaded and the marked values and tipping points labelled. Supplement: a table for every lock in scope with the estimate and CI at $\delta = 0$, at the marked values and at the end of the range, the three forms of each tipping point, and $N_a$, $M_a$ and $q_a$.

**Scope [call 5].** The confirmatory comparison's two locks (primary outcome and death within 180 days) in the main text; L1 to L5 for the primary outcome in the supplement. The sensitivity locks S1 to S4 and the longitudinal TMLE sensitivity analysis are not bias-adjusted.

## 6. Restriction to patients who would almost surely have arrived alive

**Rule.** Injury severity score recorded and 8 or less, and no penetrating wound. An injury severity score below 9 means no injury graded serious or worse (abbreviated injury scale 3 or more) in any region, the threshold below which Mock et al. did not count patients as seriously injured. Penetrating wounds are excluded because gunshot wounds were the mechanism most associated with death in transit in Lilongwe (OR 38, Reid et al.) and are over-represented among Nairobi's deaths before hospital (Saidi and Oduor). Patients with a missing injury severity score are excluded from the restriction.

Head injury severity is not used, although D6 named it as an example. Whether the injury-form checkboxes were graded from the arrival Glasgow Coma Scale is an open question to co-authors (plan, call 3); if they were, restricting on them would select on a state that transport could change. An injury severity score of 8 or less already excludes head injuries graded serious or worse. **[call 6]**

**Lock R1 [call 7].** The confirmatory comparison picked by script 23, restricted by the rule, primary outcome only, with the same plan as the confirmatory lock. It is built and run through the design stage after run 2 finishes and before unblinding, so its lock and dossier hashes are reported with the others; it runs alone, as one R job, and is checkpointed like the others. Before locking, its dossier reports the share of each arm the restriction keeps. Its feasibility verdict and population gate are applied as for every lock, and the result, whatever it is, is reported next to the confirmatory estimate. It answers a different question (the effect among less severely injured patients, for whom transport may matter less), so a smaller effect there is expected and is not read as evidence of bias.

## 7. What the plan fixes before unblinding

The severe-tail rule and its fallback; the parameter and its locked ranges by comparison, with the construction and sources of section 4; the table of $N_a$, $M_a$, $k_a$ and $q_a$ for every lock in scope and every grid value (design data, from script 23); the estimator, the interval and the implementation of section 3; the grid, the tipping points and the reading of section 5; the restriction rule and R1's lock hash.

## 8. Check on simulated data before unblinding

Run after the design-stage job and before unblinding, alone as an R job, as a script `25c_check_bias_analysis.R` that calls the same function script 25b calls.

**Claim.** When deaths before arrival happen as the bias model assumes, among the severely injured patients of one arm at the locked share, the prespecified bias-adjusted estimate recovers the risk difference among patients who would arrive alive by at least one transport, and its CI covers at close to the nominal rate; the earlier crude-risk shift does not. When the deaths come from outside the severe tail, the prespecified estimate under-adjusts, and the check shows by how much.

**Data-generating process.** Each choice is set so the check runs at the registry's own scale or at the protocol's assumptions. The arm sizes and the severe tail's share in each arm are L3's, from design data (aggregate counts only). A right-skewed severity score defines $T$, and two binary confounders stand for the rest of the adjustment set; treatment follows a logistic model in all three, with propensity scores kept between 0.1 and 0.9 so that overlap does not enter the check. The outcome follows a logistic model in which severity raises the risk steeply, with the overall comparator risk at the protocol's 0.22 and a true effect equal to the protocol's smallest important difference, $-0.055$. Patients for whom a transport choice was made are generated first; deaths before arrival are then removed from the receiving arm. Scenarios: no deaths before arrival; deaths in arm 1's severe tail at $\delta_1 = 9$ and at 18 per 1,000; deaths in arm 0's severe tail at $\delta_0 = 18$; and, as the boundary, 18 per 1,000 removed from arm 1 in proportion to the outcome risk across all severity levels, so that some come from outside the severe tail. The prespecified analysis is given the true $\delta$.

**Estimands.** The risk difference among all patients who would arrive alive with at least one transport (the target), computed from one draw of $10^6$; and, for the approximation of section 3, the same effect standardized to the arrived cohort.

**Methods.** The unadjusted estimate on arrivals; the crude-risk shift of scripts 61 and 69; the prespecified bias-adjusted estimate; and an oracle, the same TMLE run on the full population with the deaths before arrival observed. All use main-terms logistic regression for both nuisance functions (correctly specified, since the check is about the adjustment and its interval, not the learners) with five-fold cross-fitting, the same targeting code as the package, and propensity truncation 0.01.

**Performance and display.** One table: rows the four methods; columns the five scenarios; cells the bias in percentage points with its Monte Carlo SE, and the coverage of the 95% CI. Expected: the prespecified estimate and the oracle have bias near zero and coverage near 0.95 in the first four scenarios; the unadjusted estimate is biased by roughly $q_a$ times the severe tail's share times $(1 - m_T)$; the shift method overcorrects and covers below 0.95; in the boundary scenario the prespecified estimate under-adjusts. A second line reports the gap between the two estimands of the approximation.

**Pilot and size.** Two repetitions to check the code; one large draw to check the severe tail's share, the risks and the propensity range; then 200 repetitions to size the run with the skill's rule, three Monte Carlo SEs for the coverage gap between the prespecified and shift methods and for each bias contrast. Expected size 1,000 to 2,000 repetitions (Monte Carlo SE of coverage about 0.007 at 1,000), minutes of compute with logistic regression. The check also confirms that the targeting step accepts an outcome between 0 and 1 (`tmle` with `family = "binomial"`, as the package calls it; if it does not, `family = "gaussian"`, which `tmle` fits with the same logistic fluctuation for an outcome bounded by 0 and 1), and that $\delta = 0$ returns the unadjusted estimate exactly.

If the check shows the prespecified estimate biased or under-covering where the claim says it should not, the computation is fixed before unblinding and the change is recorded here.

## 9. Change to `manuscript/methods.qmd` (approved 2026-10-20 and made)

1. Line 178, the placeholder. The sentence that said the two analyses were registered is replaced by four sentences that say they were prespecified in the analysis plan, citing `foxApplyingQuantitativeBias2021`, `reid2018`, `purcell2020`, `mulima2021`, `henry2012` and `mock1998`.
2. Line 177, optional and not approved, so not made. The registry identifies the effect among patients who would arrive alive by either transport; the effect among all patients for whom a transport choice was made differs from it by the share of patients who would die before arrival by any transport (section 1).

## 10. Calls for Andrew

1. **Parameter scale:** excess deaths before arrival per 1,000 arrivals in the arm receiving them, converted to a share of the arm's severely injured patients from design data.
2. **Severe tail:** injury severity score 16 or more, or severe head injury; widened to 9 or more, or moderate or severe head injury, if an arm has fewer than 20 such patients.
3. **Ranges:** 0 to 18 per 1,000 for comparisons involving arrival by other means (L1, L2, L4, L5, S2, S3); 0 to 9 per 1,000 between ambulance operators (L3, S1, S4); marked values 4.5, 9 and 18.
4. **Both directions** use the same range.
5. **Scope:** confirmatory comparison (both outcomes) in the main text; L1 to L5, primary outcome, in the supplement; no bias analysis for S1 to S4 or the longitudinal TMLE.
6. **Restriction rule:** injury severity score recorded and 8 or less, and no penetrating wound; head injury severity not used.
7. **Lock R1:** one more design-stage lock, the confirmatory comparison restricted, primary outcome only, run after run 2 and before unblinding.
8. **Deterministic analysis, not probabilistic:** the sources give the ends of a range, not a distribution over it, so a probabilistic bias analysis would rest on a shape the data do not support.
9. **Not pursued:** Rescue.Co's dispatch records may note deaths at the scene or in transit, which would bound the Rescue.Co arm only (time-zero note, 2026-10-07). Asking Rescue.Co for them would replace a judgment with data for one arm; it needs Andrew's word and does not touch the registry's outcomes.
10. **Dating the plan:** with no OSF registration, the commit that fixes this document is entered in the blinding ledger before script 24 runs, and its hash is reported in the supplement with the lock and dossier hashes.

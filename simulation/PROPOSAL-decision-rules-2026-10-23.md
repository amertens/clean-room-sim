# Proposed decision rules for the simulation study, to fix before the contrast pilot is read

Proposal for Andrew, 23 October 2026 (audit recommendation 2, with recommendations 23 and 24). Nothing in `simulation/PLAN.md` has changed. Section 1 is a new subsection for PLAN.md, to go after "Revisions after the Methods review". Section 2 gives before and after text for the lines of PLAN.md it supersedes. Section 3 lists the code changes the rules need, and section 4 the choices that are Andrew's. The margins and levels below are design choices, not results.

## 1. New subsection for PLAN.md

> ### Decision rules, fixed before the contrast pilot is summarized (2026-10-23)
>
> **What had been seen when these rules were written.** Pilot steps 1 and 2 (calibration, DGP diagnostics, and smoke runs on two datasets in each of two scenarios), the 40-dataset check of the answer key in three scenarios (`pilot/key_check_output.txt`), and the timing runs. The contrast pilot's verdicts and label fits are cached, but its summary (`pilot.R` step 3a) has not been run, and no agreement, contrast, claim 3 or claim 4 result from it has been seen. Every rule below is fixed before step 3a is run. A rule changed after that is recorded here with its date and with the outputs seen before the change.
>
> **Claim 1a: the verdict against overlap diagnostics.**
>
> - *Measure.* Balanced accuracy, the mean of sensitivity and specificity, for each estimand and sample size, over the five overlap levels and the truths inside or near the family, computed per dataset index and paired across verdict sources on the same datasets. Its Monte Carlo SE is the SD of the per-index difference over the square root of R.
> - *Comparators.* Three, in this order. First, the best tuned overlap diagnostic, whose threshold is chosen against the answer key in the same cell; it is generous to the comparator, since no analyst knows the answer key. Second, the same diagnostics with the threshold tuned at the other sample size, which an analyst could carry over from a published simulation. Third, the package's overlap grade (good below 5% outside the band, moderate below 20%), read as calling the ATE and ATT estimable when the grade is good and the trimmed ATE and ATO estimable when it is good or moderate.
> - *Pass rule, per estimand and sample size.* The verdict is **better** than a comparator when the lower bound of the 95% interval for the difference in balanced accuracy is above 0. It is **no worse** than the best tuned diagnostic when that lower bound is above -0.05. Otherwise the comparison is **unresolved**, and is reported as such.
> - *Claim 1a holds* in its stronger form if the verdict is better than all three comparators for every informative estimand and sample size, and in its weaker form if it is no worse than the best tuned diagnostic and better than the second and third comparators in every one. The weaker form is stated as "as accurate as a diagnostic tuned on the true outcome model, without a threshold that no analyst could tune". An estimand whose label does not vary across levels and truths at a sample size tests nothing; it is reported with its agreement rate and named as uninformative, and it does not count against the claim.
> - *Margin.* 0.05 in balanced accuracy, one estimand in twenty misgraded on average over the two label classes. Differences smaller than this are below the misgrading that the verdict's own noise near a tolerance produces at 200 repetitions per verdict (audit 03, section 3.3, binomial arithmetic).
>
> **Claim 1b: agreement outside the family.** The measure is the verdict's balanced accuracy for the truths outside the family (the outside model, and the nonlinear-treatment truth if it is added) minus that for the truths inside or near it, paired by dataset index. Agreement **falls** when the upper bound of the 95% interval for this difference is below 0, and is **not shown to fall** when the interval lies within ±0.05. Otherwise the comparison is unresolved. Because the verdict never sees the outcome, this difference comes only from labels that differ across truths, and the paper says so.
>
> **Claim 3: the outcome-informed comparators.** Both rules are fixed now and both are reported, so that neither can be chosen after the other is seen. Each chooses among the estimands that the dataset's own verdict judged feasible, the set the prespecified order chooses from: the estimand with the smallest P value, which selects on the estimate, and the estimand with the narrowest CI, which selects on the width. The measures are coverage of the 95% CI for the chosen estimand's true value (effect scenarios) and the rejection rate at 0.05 (null scenarios), paired by dataset with the prespecified choice. They are computed on the datasets where the verdict judged at least two estimands feasible, since elsewhere the choices coincide, and the share of datasets with two, three or four feasible estimands is reported.
>
> - Selection on the estimate **inflates the type I error** when the lower bound of the 95% interval for its rejection rate minus the prespecified choice's is above 0.
> - Selection on the width is **reported as not changing** coverage or the rejection rate when the 95% intervals of both differences lie within ±0.03. Otherwise its estimates and intervals are reported without a claim.
> - The claim is scoped to what these show. No claim is made that the gap grows with the number of feasible estimands unless a trend over two, three and four feasible estimands has a 95% interval that excludes zero.
>
> **The baseline-risk factor.** Overlap is the only factor in the main grid that changes whether an estimand is estimable, and every comparator diagnostic reads overlap, so the design gives the verdict no condition that a summary of the propensity score cannot see. The outcome's baseline risk is such a condition, because it sets the number of events.
>
> - *Truths.* The inside and outside models with control-arm risk 0.03, 0.15 (the main grid) and 0.30, b0 recalibrated for each on the population draw, at the three middle overlap levels (gamma 1.5, 1.75 and 2.0) and n = 500. Covariates and treatment are shared with the main grid, so these add label fits only.
> - *Family.* The verdict is run with the family's baseline risk at 0.03, 0.10, 0.15 and 0.30 on the same datasets. The run at 0.15 is the main grid's. The runs at 0.10 and 0.30 are the misspecified variants approved on 7 October, here at these three levels and n = 500. The run at 0.03 is new. Each family setting is crossed with each truth, so a matched family (an analyst who knows the baseline risk) and a misspecified one are both shown.
> - *A fair comparator.* An analyst who looks only at the propensity score cannot see the baseline risk, so against those diagnostics the verdict wins on this factor by construction. The comparator for this factor is therefore the best overlap diagnostic combined with an events rule an analyst could apply before the outcome: expected events in the arm with fewer, from the baseline risk the analyst gave the family, at least 10 per covariate column (`peduzziSimulationStudyNumber1996`). The package already reports these expected events before the outcome (`simulate_design()`, `expected_events`).
> - *Reporting.* The factor is reported in its own table, not pooled into Table 1. The claim it can support is what the simulation adds over overlap and an events rule judged separately, using claim 1a's rule and margin.
>
> **Dropping a truth.** A truth leaves the main display, and moves to the supplement, only if its estimability labels in the full run equal those of a truth kept in the display in every overlap level, estimand and sample size. Truths are kept in the order inside, outside, nonlinear treatment, near, so the later one moves. No truth is moved or redefined for any other reason, in particular not because the verdict agrees with it poorly.
>
> **Contrasts, levels and sizing.** Step 3a is used to size the full run: R and M are set by the largest requirement among the contrasts above, at k = 3 for each "better" or "falls" comparison and by the equivalence bound (`pilot_equivalence()`) for each margin. A contrast that the affordable R cannot resolve stays in its table, marked unresolved, and the claim it serves is rescoped to match. It is not removed. No overlap level, truth or comparator is added or removed after step 3a is read, except by the rule for dropping a truth.

## 2. Lines of PLAN.md this supersedes

### Line 24 (Goal 1, "Sharpened")

Before:

> - **Sharpened:** across scenarios crossing overlap (good, moderate, severe) with the true outcome model (inside, near, outside the declared family) at n = 500 and 2,000, the generate-treatment support map has higher sensitivity and specificity for estimability than sample-treatment and than every best-threshold overlap diagnostic when the truth is inside or near the family, and its agreement falls when the truth is outside it.

After:

> - **Sharpened:** across the five overlap levels crossed with the true outcome models at n = 500 and 2,000, the generate-treatment support map's balanced accuracy for estimability is higher than each comparator's, or no more than 0.05 below the best tuned diagnostic's, for each informative estimand when the true outcome model lies inside or near the family; and its agreement falls when the model lies outside it. The pass rules are in "Decision rules, fixed before the contrast pilot is summarized".

### Lines 29 and 30 (Goal 2)

Before:

> - On the same R datasets per scenario, with the real outcome: the prespecified estimate (the ladder's primary from that dataset's map) versus an outcome-informed choice (among the estimands, the one with the smallest p-value). Measures: coverage of the 95% CI for the chosen estimand's true value; under a null effect (separate null versions of each scenario), the rejection rate at 0.05.
> - **Sharpened:** the prespecified choice covers at the nominal rate (within Monte Carlo error) where the map says the estimand is feasible, while the outcome-informed choice undercovers and rejects a true null more often than 5%, with the gap growing as more estimands are available.

After:

> - On the same R datasets per scenario, with the real outcome: the prespecified estimate (the ladder's primary from that dataset's map) versus two outcome-informed choices among the estimands the dataset's map judged feasible, the one with the smallest p-value and the one with the narrowest CI. Measures: coverage of the 95% CI for the chosen estimand's true value; under a null effect (separate null versions of each scenario), the rejection rate at 0.05; both on the datasets where at least two estimands are feasible.
> - **Sharpened:** choosing by the smallest p-value rejects a true null more often than the prespecified choice, while choosing by the narrowest CI changes neither coverage nor the rejection rate by more than 0.03.

### Table S1 (lines 65 to 70)

The caption reads "...versus after seeing outcomes (smallest p-value or narrowest CI), on the datasets where at least two estimands are feasible...", and a third row is added: "Outcome-informed (narrowest CI)".

### Line 116 (misspecified baseline risk)

Before: "The generate-treatment map is run twice more on the same datasets with the family's baseline risk set to 0.10 (the package default, an underestimate) and 0.30 (an overestimate), to show the cost of misjudging the baseline risk. Table 1 gets a row for each." After: "Superseded on 2026-10-23 by the baseline-risk factor in 'Decision rules, fixed before the contrast pilot is summarized', which keeps these two runs, at the three middle overlap levels and n = 500, and reports them in their own table."

### Line 119 (claim 3 comparator)

Before: "The outcome-informed choice is the feasible estimand with the narrowest CI (a milder rule than the smallest P value)." After: "Superseded on 2026-10-23: both the narrowest-CI rule and the smallest-P rule are reported, each among the feasible estimands, as set out in 'Decision rules, fixed before the contrast pilot is summarized'."

## 3. Code changes the rules need

None of these runs until the design stage has finished, and step 3a runs first, on the cached 0.4.1 results, after the rules are approved. Step 3a's summary code then needs four changes before it is read: balanced accuracy per dataset index with its paired difference and 95% interval (`R/summarize.R`, beside `agreement_table()` and `table1_contrasts()`, which today report sensitivity and specificity and pick the best diagnostic by their sum); the two analyst-usable comparators (a threshold carried over from the other sample size, and the grade rule) beside the tuned one in `add_diagnostic_calls()`; the claim 1b difference; and both claim 3 rules restricted to the feasible estimands (`R/key.R:26-30` takes the smallest P value over all four estimands, and `claim3_rows()` reads it). Because the contrast pilot ran at n = 500 only, the threshold carried over from the other sample size can be sized but not computed at step 3a; it is computed in the full run.

The baseline-risk factor needs the truth's control-arm risk and the family's baseline risk as arguments (`R/dgp.R`, `R/map.R:21`), the expected events kept from each verdict run (`map.R` stores the verdict rows but not `expected_events`), the events rule, and its own sizing run, since the contrast pilot did not include it. The extra cost is one verdict run per dataset at three levels for the family at 0.03, label fits for the two new truths, and the two approved variants at three levels instead of five.

## 4. Choices for Andrew

The margins are the main judgments here. A margin of 0.05 in balanced accuracy is the one the audit proposed; a tighter 0.03 costs more repetitions. As an illustration only, if the per-index difference in balanced accuracy has an SD of 0.15 (an assumption, which step 3a replaces with its estimate), equivalence within 0.05 at 95% needs about 144 datasets per cell, and within 0.03 about 400 (`pilot_equivalence()`'s bound, (2 × 2 × SD / margin) squared). The claim 3 margin of 0.03 on coverage and the rejection rate is likewise a judgment; 0.02 would roughly double the datasets needed.

Three proposals reopen decisions of 7 October, and each is his call. Reporting the smallest-P rule beside the narrowest-CI rule reverses the choice of the narrowest CI as the only comparator; the reason predates any claim 3 result (the width of an interval is nearly independent of its estimate, so the narrowest-CI rule is expected to show no harm). Running the misspecified-baseline variants at three levels and n = 500 narrows a revision that named no levels. And the grade rule's reading (good for the ATE and ATT; good or moderate for the trimmed ATE and ATO) is this proposal's own, and he may prefer another.

One limit should be stated wherever claim 1a is reported: inside the family the agreement is partly definitional, because the inside truth is built from the family and the label uses the verdict's own tolerances. The informative cells for claim 1a are the near truth and the baseline-risk factor, and for claim 1b the outside truths.

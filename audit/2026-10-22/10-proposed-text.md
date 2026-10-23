# Proposed text for the main vignette and the paper's statement of contribution

Proposals for Andrew, 22 October 2026, from the audit in this folder. No file outside this folder was changed. Each place gives the text as it stands, quoted exactly with its lines, then the proposed text, then a short paragraph on the change and the recommendation of `REPORT.md` it serves. Numbers in the proposed text come from `writing/registry.json` by key, through inline R, or from the vignette's own data at render time. Citation keys are keys of `writing/references.bib`, each checked there. A `[cite: ...]` marks a source that sci-cite must add first, with the identifier verified in `08-literature.md`.

## 1. Opening of the main vignette

### Before

`cleanTMLE/vignettes/cleanTMLE.Rmd`, lines 21 to 81, from the opening sentence to the last line before the first chunk that builds a plan (`plan`, line 83; line 82 is blank). Quoted from the working tree on branch `vignette-limits`. Its one uncommitted edit, shown by `git -C cleanTMLE diff vignettes/cleanTMLE.Rmd`, touches lines `290-301` and `367-378` and none of these lines, so the working tree and the last commit agree here.

````markdown
cleanTMLE runs an observational comparison in stages, so that the design is
settled before anyone sees the outcome. This is the clean-room idea of
@muntner2024staging. The package adds one question to the design stage:
can these data deliver the estimand we planned, and if not, which one can
they deliver?

In this vignette we take one simulated cohort through every stage. It takes
a few minutes to build. Each step says what it decides and what to look for
in its output.

## Why decide before seeing the outcome

Two problems motivate the stages. The first is analytic flexibility: when
the person choosing the model can see how each choice moves the estimate,
the final estimate carries those choices, however honest the person. The
second is positivity. An average treatment effect (ATE) needs every kind of
patient to have a real chance of either treatment. Where a group of
patients almost never receives one treatment, an estimator extrapolates,
and it can return a precise-looking estimate that is far from the truth
[@petersen2012positivity].

Keeping the outcome away from the design solves the first problem.
Checking, still without the outcome, which estimands the design can support
solves the second. The other estimands describe other populations: the
average treatment effect among the treated (ATT); the trimmed ATE, the
effect among patients whose propensity score lies inside a band
[@crump2009trimming]; and the overlap-weighted effect (ATO), which weights
each patient by how uncertain their treatment was [@li2018overlap]. In ICH
E9(R1) terms, moving between them changes the population attribute of the
estimand; the treatment, the outcome and the summary measure stay the same.

cleanTMLE estimates each of them with targeted maximum likelihood
estimation (TMLE), which combines a model for the outcome with a model for
the propensity score and stays consistent if either one is correct
[@vanderLaan2006tmle].

## The example cohort

`cr_example` is a simulated cohort of 1,500 patients. It comes as two
separate data frames, so that the outcome stays out of the design stage.

```{r data}
library(cleanTMLE)
str(cr_example$design)
str(cr_example$outcomes)
```

The question is the effect of treatment `A` on death within 30 days
(`death_30d`), measured as a risk difference: the difference in the
proportion of patients who die within 30 days. A risk difference describes
the whole target population directly and, unlike a hazard ratio from a Cox
model, does not depend on a proportional-hazards assumption. Treatment
depends strongly on `severity`, so few severe patients go untreated and few
mild patients are treated. Because the data are simulated, we know the true
effect for each estimand; it is in `cr_example$truth`, and we compare our
estimate with it at the end.

## Step 1: declare the plan

Everything about the analysis is set once, in the plan. Nothing can be
added after the next step.
````

### After

This replaces lines `21-61`. Lines `62-81` (the `data` chunk, the paragraph on the question and the lead of Step 1) stay as they are.

````markdown
cleanTMLE helps you settle the design of an observational comparison of
two treatments before anyone sees the outcome. You write an analysis plan,
lock the design data without the outcome, and run a simulation on your
real covariates with outcomes that the package generates. The simulation
judges which estimands your data can support, and which estimator to use
for each. After a review, you add the outcome, and the package estimates
the effect with the estimator the simulation tested. The stages follow
the clean-room approach of @muntner2024.

You supply a plan, design data, a variable dictionary and, at the end,
the outcome. The plan lists the estimands you want in order of
preference, the tolerances for bias and coverage that an estimator must
meet, and the candidate estimators to compare. The design data have one
row per patient, with an id, the treatment, the covariates and any
negative-control outcomes, and no outcome. The dictionary says, for each
column of the design data, where it was recorded and when, relative to
time zero. Time zero is the moment at which treatment is assigned,
eligibility is met and follow-up starts, as in the protocol of a target
trial [@hernan2016], and the lock refuses a covariate fixed after it. The
outcome data, the id and the outcome, stay with whoever holds them until
the review is done.

Before the outcome is added, you get back a design report: a hashed
record of how well the treatment groups overlap, which population each
estimand describes, which estimands the simulation found the design can
deliver and with which estimator, and a check with negative-control
outcomes. After the review, you get the estimate from the selected
estimator and a log of every decision and who made it.

The approach has limits. The simulation's verdict holds only over the
family of outcome models written in the plan: if the real outcome depends
on the covariates in a way that none of those models resembles, an
estimand the simulation accepted can still be estimated with more bias
than you allowed. The lock refuses design data that contain the outcome
column, but it cannot tell whether another column carries information
about the outcome, such as a record that a patient was reached at
follow-up. Whether the design data are free of the outcome therefore
rests on how they are built, which in a clean room is the job of the data
steward, the person who holds the outcome.

In this vignette we take one simulated cohort through every stage. It
takes a few minutes to build; a real plan with the default settings
takes much longer, and `simulate_design()` prints an estimate of its run
time after the first repetition. Each step says what it decides and what
to look for in its output.

| Step | Function | What it decides or returns |
|---|---|---|
| Step 1 | `analysis_plan()` | The estimands in order of preference, the tolerances, the candidate estimators and the negative controls |
| Step 2 | `create_analysis_lock()` | Fixes the plan, the design data and the dictionary under one hash |
| Step 3 | `assess_design()` | How well the groups overlap, and whom each estimand describes |
| Step 4 | `simulate_design()` | Which estimands the design can deliver, and with which candidate |
| Step 5 | `negative_control_ladder()` | Whether a negative-control outcome signals residual confounding |
| Step 6 | `design_report()`, `unblind()` | The design report for review, then the outcome once it is approved |
| Step 7 | `estimate_effect()` | The estimate, from the candidate the simulation selected |
| Step 8 | `export_design_log()` | The record of every decision |

## Why decide before seeing the outcome

The first reason is analytic flexibility. When the person choosing the
model can see how each choice moves the estimate, the final estimate
carries those choices, however honest the person. The second is
positivity. An average treatment effect (ATE) needs every kind of patient
to have a real chance of either treatment. Where a group of patients
almost never receives one treatment, an estimator extrapolates, and it
can return a precise-looking estimate that is far from the truth
[@petersenDiagnosingRespondingViolations2012].

Keeping the outcome out of the design stage limits the first problem,
because no design choice can respond to the estimate. It does not close
every route to a biased answer: when the outcome can be predicted from
the covariates, choices made on the covariates alone can still bias the
estimate [cite: Fisher 2020, hiding outcomes does not prevent bias when
outcomes are predictable from covariates, arXiv:2007.02514]. Checking,
still without the outcome, which estimands the design can support tells
you whether the second problem is present, and which estimand the data
can support instead. The other estimands describe other populations: the
average treatment effect among the treated (ATT); the trimmed ATE, the
effect among patients whose propensity score lies inside a band
[@crump2009]; and the overlap-weighted effect (ATO), which weights each
patient by how uncertain their treatment was [@li2017]. In the terms of
the estimand framework
[@internationalcouncilforharmonisationoftechnicalrequirementsforpharmaceuticalsforhumanuse2019],
moving between them changes the population; the treatment, the outcome
and the summary measure stay the same.

cleanTMLE estimates each of them with targeted maximum likelihood
estimation (TMLE), which combines a model for the outcome with a model
for the propensity score and stays consistent if either one is correct
[@vanderlaan2006]. The package fits both models itself, with its own
cross-fitted super learner (Step 1 explains both terms), so that the
simulation and the final analysis run the same fitter with the same
folds, learners and truncation. It hands the last step of TMLE, the
targeting step, to the tmle package [@gruber2012] for the ATE, the
trimmed ATE and the ATT, and runs that step itself for the ATO.

## The example cohort

`cr_example` is a simulated cohort of `r nrow(cr_example$design)`
patients. It comes as two separate data frames, so that the outcome stays
out of the design stage, and with the variable dictionary for the design
data in `cr_example$dictionary`. In this cohort time zero is treatment
assignment, and every covariate was recorded before it.
````

### What changed and why

The new opening tells a first-time user what to supply (a plan, design data, a variable dictionary and, at the end, the outcome) and what comes back (a design report before the outcome is added, and an estimate and a log of decisions after it), defines time zero at its first use with the target-trial reading that `analysis_plan()` documents, and lays out the steps in a table, as pair A of `05-writing.md` proposes; this is recommendation 15. It states once, in the opening, that the verdict holds only over the outcome models written in the plan, and it says that the design data's freedom from the outcome rests on how the data steward builds them, because the lock refuses only the outcome column and keeps columns such as the outcome-observation indicators (recommendations 4 and 15). The paragraph on TMLE now says which parts the package does itself (both nuisance fits, the cross-fitting and the targeting for the ATO) and which it hands to tmle (the targeting for the ATE, the trimmed ATE and the ATT), without calling the package a wrapper (recommendation 4). The sentence on what withholding the outcome does not guarantee is the first of the four literature items of recommendation 17 and waits for Fisher 2020 to be added through sci-cite. "Solves" is gone (pair D), the typed count of patients becomes inline R, and "design report" follows the term set recommended under recommendation 15, which awaits Andrew. The citation keys are those of the project bibliography; the vignette's own bibliography uses other keys for five of these works (`muntner2024staging`, `petersen2012positivity`, `crump2009trimming`, `li2018overlap`, `vanderLaan2006tmle`), so applying the text means copying into `cleanTMLE/vignettes/references.bib` the entries for `muntner2024`, `hernan2016`, `petersenDiagnosingRespondingViolations2012`, `crump2009`, `li2017`, the estimand addendum, `vanderlaan2006` and `gruber2012`. The other parts of pair A that lie outside these lines (two candidates in the example plan, one refusal of the lock shown, the dictionary's fields and roles) are in step 11 of `09-spec-package-simplification.md`, and the vignette is rendered only after the Rescue.Co design stage has finished.

## 2. The paper's statement of contribution

### Before, introduction

`manuscript/introduction.qmd`, lines 20 to 31, the third and fourth paragraphs. Quoted from the working tree, which matches the last commit for this file.

````markdown
Simulations built on a study's own covariates can evaluate estimators without the outcome.
Plasmode simulations keep the real covariates and generate outcomes from a chosen model [@franklin2014], and simulations of this kind are used in the causal roadmap to choose the estimator before the real analysis is run [@nanceCausalRoadmapSimulations2024a; @williamsonApplicationCausalRoadmap2023].
How such a simulation generates treatment matters: keeping each patient's observed treatment can build a positivity violation into the simulated data [@shaw2026], and a plasmode simulation can show that an estimator fails but cannot certify that it will succeed in the population the data came from [cite: Karim 2026 plasmode refutes but cannot certify].
For choosing the estimand itself, design-based metrics can weigh the statistical bias of an estimator against departure from the population of interest [cite: Barnard 2026 estimand selection under positivity violations], and prespecified objective diagnostics, read before the results are unblinded, can screen out analyses likely to be biased [cite: Conover 2025 objective study validity diagnostics].
These approaches read the design through the propensity score, the weights or negative controls, while under poor overlap the bias and coverage of an estimator also depend on the outcome regression.
Whether a simulation that never uses the outcome can say which estimand a design can deliver, and what choosing the estimand at that stage does to inference, are the questions this paper takes up.

Here, we propose choosing the estimand at the design stage with an outcome-blind plasmode simulation, inside a staged workflow implemented in the R package cleanTMLE.
The simulation draws new treatments from the propensity score estimated on the study's covariates and outcomes from a prespecified family of outcome models, judges an estimand feasible when some candidate estimator meets prespecified tolerances for bias and coverage on every model in the family, and makes the first feasible estimand in a prespecified order the primary estimand, all before the outcome is joined to the data.
Our objectives were to estimate, in a simulation study, the sensitivity and specificity for estimability on the true outcome model of this verdict, of the same simulation with each patient's observed treatment kept, and of overlap diagnostics, by estimand and by how closely the true outcome model resembles the simulated family; the coverage and type I error of the estimate when the estimand is chosen before rather than after the outcome is seen; and the agreement between the confidence-interval half-width forecast before the outcome and the half-width realized after it.
We also show the workflow in a re-execution of a trauma registry analysis from Nairobi, Kenya, with the outcome withheld from the software.
The aims are methodological; the application estimates a causal effect, whose estimand and identifying assumptions are given in the Methods.
````

### After, introduction

This replaces lines 20 to 31 with four paragraphs. The inline R uses the `reg` object the drafting workflow loads in the setup chunk.

````markdown
Simulations built on a study's own covariates can evaluate estimators without the outcome.
Plasmode simulations keep the real covariates and generate outcomes from a chosen model [@franklin2014], and the causal roadmap uses simulations of this kind, after the data are collected and before any effect is estimated, to choose the estimator and complete the statistical analysis plan [@nanceCausalRoadmapSimulations2024a; @williamsonApplicationCausalRoadmap2023].
The same approach has been used to prespecify a data-adaptive estimator in randomized trials [cite: Balzer 2026, prespecified data-adaptive targeted learning tested by parametric and plasmode simulations, DOI 10.1177/17407745261417227], and a software tool compares estimators for a prespecified ATE by plasmode simulation on the analyst's own data [cite: Meng and Huang, REFINE2 plasmode comparison of estimators for a prespecified ATE, DOI 10.1093/aje/kwaf195].
In each of these uses the estimand is fixed before the simulation is run, so the simulation chooses among estimators and not among estimands.
For choosing the estimand, design-based metrics can weigh the statistical bias of an estimator against departure from the population of interest [cite: Barnard 2026, estimand selection under positivity violations, DOI 10.1093/biomtc/ujag014], and objective diagnostics, prespecified and read before the results are unblinded, can screen out analyses likely to be biased [cite: Conover 2025, objective study validity diagnostics, DOI 10.1093/jamia/ocae317]; a large drug-safety program did not unblind the evaluations that failed them [cite: Haynes 2026, diagnostics applied before unblinding in a safety program, DOI 10.1002/pds.70435].
These gates withhold the effect estimate but not the outcome, because the minimum detectable relative risk is computed from the observed number of outcome events and the expected systematic error from the estimates for negative-control outcomes [cite: Conover 2025, DOI 10.1093/jamia/ocae317].
These metrics and gates judge the design from the propensity score, the weights, outcome counts or negative controls, while under poor overlap the bias and coverage of an estimator also depend on the outcome regression.
Whether a simulation that never uses the outcome can say which estimand a design can deliver, and what choosing the estimand at that stage does to inference, are the questions this paper takes up.

Here, we extend outcome-blind plasmode simulation from the choice of the estimator to the choice of the estimand.
Before the outcome is joined to the data, the simulation draws treatments from the propensity score estimated on the study's covariates and outcomes from a prespecified family of outcome models, and it judges that the design can deliver an estimand when some candidate estimator meets prespecified tolerances for bias and coverage on every model in the family.
An order of estimands fixed in advance makes this verdict a decision, because the first estimand the design can deliver becomes the primary estimand.
The setting is a staged workflow that follows the causal roadmap [@dang2023], estimates each effect by targeted maximum likelihood estimation (TMLE) [@vanderlaan2006], and keeps the outcome out of the design stage until the statistical analysis plan is approved, as in a clean room [@muntner2024].
A flexible estimator belongs in such a workflow because its learners, tuning and truncation are choices that an analyst who has seen the outcome could steer [@kahanHowDesignPrespecified2020].
We implemented the workflow in the R package cleanTMLE, which builds on the tmle package [@gruber2012], records the plan, the design data and each design decision under a cryptographic hash, and after unblinding applies the estimator that the simulation evaluated, with the same fitter, folds, learners and truncation.
The design stage therefore produces, before any outcome is seen, the assessment of feasibility and the anticipated operating characteristics that recent guidance asks an analysis plan to state with its estimand [cite: Food and Drug Administration 2024, draft guidance on non-interventional studies, fda.gov/media/177128/download], and unlike the diagnostic gates in use today it reads no outcome count.

Our objectives were to estimate, in a simulation study, how often the verdict agrees with whether some candidate meets the same tolerances on the true outcome model, by estimand and by how closely that model resembles the simulated family, and how this agreement compares with that of overlap diagnostics and of the same simulation with each patient's observed treatment kept, which can build a positivity violation into the simulated data [@shaw2026].
We also estimated the coverage and type I error of the estimate when the estimand is chosen before rather than after the outcome is seen, and how well the confidence interval half-width forecast before the outcome matches the half-width realized after it.
To show the workflow on real data, we re-analyzed a trauma registry from Nairobi, Kenya, with the outcome withheld from the software; the aims are methodological, and the application's estimand and identifying assumptions are given in the Methods.

The verdict holds only over the family of simulated outcome models, and like any plasmode simulation it can show that an estimator fails but cannot certify that it will succeed [cite: Karim 2026, a plasmode benchmark can refute an estimator but not certify it, arXiv:2609.26186].
The hashes make a change to the plan, the design data or a design decision detectable without preventing analysis of the outcome outside the software, and nothing in the workflow removes bias from unmeasured confounding.
Withholding the outcome does not by itself prevent bias, because when the outcome can be predicted from the covariates, design choices made on those covariates can still favor one answer [cite: Fisher 2020, hiding outcomes does not prevent bias when outcomes are predictable from covariates, arXiv:2007.02514].
The simulation study has far fewer covariates than the `r reg$final_count_covariate_columns$text` covariate columns of the application, and none of its sample sizes is as small as the application's four-month comparison between ambulance operators (`r reg$first_four_months_operators_n$text` patients, against `r reg$final_operators_n$text` in the final cohort).
The application shows the software at work, not the effect of blinding, because its analysts had seen the results of an earlier, unblinded analysis of the same registry.
````

### Before, brief

`writing/brief.md`, line 14 (the message) and line 21 (key claim 5). Quoted from the working tree; the file has no uncommitted change.

````markdown
(proposed) Before any outcome is seen, an outcome-blind plasmode simulation can tell which causal estimands a design can deliver, more accurately than overlap diagnostics, when the true outcome model resembles the simulated ones, and choosing the estimand at that stage keeps inference valid where choosing it after seeing outcomes does not.

5. (proposed) In the Rescue.Co Kenya trauma registry, the outcome-blind design stage, run on the first four months of enrolment and on the final cohort, flagged before any outcome was joined the identification problems and design choices that the unblinded analysis changed while its estimates were available, and chose the confirmatory comparison, the estimand and the estimator by rules fixed in advance. (To be rewritten to what the re-execution shows, once it has run.)
````

### After, brief

````markdown
(proposed) Outcome-blind simulation, which has been used to choose the estimator for an estimand fixed in advance, can also judge before any outcome is seen which estimands a design can deliver, and an order of estimands fixed in advance turns that judgment into the choice of the primary estimand. The simulation study measures how often the judgment is right, when the true outcome model resembles the simulated ones and when it does not, and what choosing the estimand before rather than after seeing outcomes does to coverage and type I error; cleanTMLE makes the sequence checkable and applies after unblinding the estimator that the simulation evaluated.

5. (proposed) In a re-analysis of the Rescue.Co trauma registry in Nairobi, Kenya, with the outcome withheld from the software, the design stage, run on the first four months of enrollment and on the final cohort, recorded before any outcome was joined which estimands each comparison could deliver and which estimator the simulation selected for each, and the estimate for the confirmatory comparison came from the estimator the simulation selected. At four months it found no estimand that the comparison between ambulance operators could deliver, although no patient in that comparison lay outside the overlap band. The analysts wrote the decision rules after seeing the results of an earlier, unblinded analysis of the same registry, so the application shows that the rules can be applied without the outcome, not that they would have made the same choices without that knowledge. (Confirm the estimate clause once the estimates exist, and align this claim with the case-study claims of recommendation 19 of the 22 October audit.)
````

### What changed and why

Recommendation 5 sets the headline as an extension of outcome-blind simulation from the choice of the estimator to the choice of the estimand, with the causal roadmap, TMLE and the clean room as its setting, and this proposal follows it rather than pair B of `05-writing.md`, which offered the joined workflow and the choice of both estimand and estimator as the advance. The literature survey found a published source for each part of that workflow (`08-literature.md`, section 1), and the review of the framing refuted it as the headline (`07-check-framing.md`, recommendation 5). The simulation paragraph takes place 1 of `08-literature.md`, section 6, and ends on the gap the survey supports, that in the uses found the estimand is fixed before the simulation runs; the narrower fallback, that no check made before the outcome tests whether an estimator delivers the estimand, is not used, because the roadmap simulations, the trial version and REFINE2 contradict it. The sentences on design-based metrics and gates take place 2, so that the reader learns before the contribution that the gates in use read outcome counts. The contribution paragraph takes place 3 and the four parts the report names (the extension, the order fixed in advance, the test of how often the verdict is right, and software in which the estimator applied after unblinding is the one the simulation evaluated), keeps from pair B the reason a flexible estimator belongs in an outcome-blind workflow, and argues the advance on the two points the report gives: the draft guidance of the U.S. Food and Drug Administration asks an analysis plan to state the estimand with an assessment of its feasibility and its anticipated operating characteristics, which the design-stage tables provide before any outcome exists, and the gates in use today read outcome counts while cleanTMLE's design objects hold no outcome. The guideline of the International Council for Harmonisation that the survey lists (`M14`) is not cited for that sentence, because the feasibility it asks about concerns the data source (`08-literature.md`, section 4c). The last paragraph carries every limit the report wants on the first page, with the three sizes inserted from the registry (`final_count_covariate_columns`, `first_four_months_operators_n` and `final_operators_n`) and the word blinding, which the brief asks for. In the brief, the message drops "more accurately than overlap diagnostics" and "keeps inference valid", which are claims 1 and 3 still awaiting evidence (`05-writing.md`, section 8), and claim 5 is cut back to what the case study can carry, with the provenance sentence that pair H proposes for the Methods and the four-month result that the report calls the clearest case of the simulation seeing what an overlap diagnostic does not.

Four points need Andrew. The headline is his decision (report, decisions table), and the text stays a proposal until the near work has been read and added through sci-cite. He co-authored the roadmap papers, so he can confirm that their simulations fix the estimand first, and whether `williamsonApplicationCausalRoadmap2023`, kept from the current text, describes simulation, which audit 05 could not tell from its abstract. The sentence that cleanTMLE applies the estimator the simulation evaluated holds under the same code and package versions, so it should go together with the Methods sentence the report proposes for `manuscript/methods.qmd:71`, once the build record of recommendations 1 and 7 exists. And the introduction now has six paragraphs, one more than the introduction rules allow; joining the last two restores five at the cost of a paragraph of about 300 words, while moving the limits to the Discussion would take them off the first page, where the report wants them.

## Writing check

The check was run from the project root in Git Bash with `--min-level warn`, on this file and, as the task asks, on each after passage alone, written to a temporary `.md` file in this folder and deleted after the run. The before and after texts sit in fenced blocks in this file, which the check does not read, so the runs on the passages alone are the check of the proposed text.

On the vignette's after text the check reports one citation placeholder (Fisher 2020), whose preprint identifier it also reads as an identifier and an untraced number, and three names in backticks (`simulate_design()`, `cr_example` and `cr_example$dictionary`), which are the package's interface. On the introduction's after text it reports the nine citation placeholders, the identifiers inside them read as untraced numbers or labels, and one abbreviation inside them, the one for the digital object identifier; all wait for sci-cite. On the brief's after text it reports nothing. No dash, banned word, stacked hedge, glossary label or typed number remains in any after text, and the three numbers in the introduction come from inline R with registry keys. The first run on the passages showed that the checker splits a sentence at "U.S.", so the placeholder for the draft guidance went undetected; it was reworded so that the check and sci-cite both see it.

On this file the check reports only code formatting in prose (paths, registry keys and citation keys in backticks in the explanatory paragraphs), after line numbers and a guideline name in those paragraphs were put in backticks. Every run also warns that the overlap source in the unblinded analysis's repository is not found; I did not look for it, so overlap with that manuscript was not checked.

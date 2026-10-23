# Writing audit: vignettes, introduction, methods and proposed results

Audit of 2026-10-22. This is a document of proposals. It changes no vignette, manuscript file, specification or package code, and it is the only file I wrote.

This audit covers the two vignettes in `cleanTMLE/vignettes/` (as they stand in the working tree and as rendered in `docs/vignettes-0.5.0/`), `manuscript/introduction.qmd`, `manuscript/methods.qmd`, `writing/proposed-results-2026-10-21.qmd` and `writing/proposed-changes-2026-10-21.md`. The short answer is that the writing is accurate and honest but does not yet make the case Andrew wants. The vignettes let a newcomer reproduce the example, not analyze their own data. The introduction states the contribution as a method for choosing the estimand, not as the combination of the causal roadmap, TMLE and a clean room, and it leaves out the boundaries that the methods supply. The methods overclaim in places, describe a super learner where a simpler candidate was selected, and omit items the methods checklist asks for. The proposed results match the registry and the exported tables but are too dense to read. The largest risk to the argument is not the prose. The registry holds no entry for the simulation study, four of the closest papers are not in the bibliography, and the case study's rules were written by people who had seen the unblinded estimates. I read code and documents only. I ran no R process, opened no Rescue.Co outcome or data file, and read no file from the earlier unblinded analysis.

## How I worked

I loaded the review, introduction, methods and discussion skills and applied the three reviewer definitions as lenses myself: the methods lens (estimand, time zero, sensitivity analyses, limitations), the claims lens (causal wording, spin, numbers against the registry) and the reader lens (what the audience named in the brief would not follow). No subagents were used.

The writing check accepted the two `.Rmd` files, the three manuscript-side files and the proposed-changes memo. It rejects `.html`, so I read the rendered vignettes by extracting their text into my scratchpad (`C:/Users/andre/AppData/Local/Temp/claude/C--Users-andre-OneDrive-Documents-clean-room-sim/8ccd3c26-7d8d-477b-8e2a-fa0832cd1031/scratchpad/`) and reviewing the printed output by eye against `style/patterns.yml`. I also left one small bibliography-search helper in the system temporary folder. The check reports that the overlap source `rescueCo_analysis/reports/manuscript_rescueco.qmd` is not found, so overlap with the unblinded analysis's manuscript is not being checked at all. I did not look for that file, because the rules bar me from opening that repository. Running the check also caches source text under the writing-skills repository, which is how the tool behaves and not something I chose.

Numbers come from `writing/registry.json` by key. Every key used in the proposed results exists in the registry (I checked all of them), and the categorical statements about overlap grades and the population criterion match `rescueCo/results/reexec_050/*/design/design_summary.csv` and `writing/results-export/population_gate.csv`. Where I cite a literature figure I could not verify, I say so.

## 1. Main vignette, `cleanTMLE/vignettes/cleanTMLE.Rmd`

A reader who copies the code runs the whole workflow in a few minutes and gets readable output at each step. Each step says what it decides and what to look for, and the closing section says what the version does not do. That is a good tutorial skeleton. The problem is that the reader can run the example but not their own data, and they would not learn from this vignette how the package chooses among estimators, which is its central promise.

The table lists where a first-time user gets lost (reader lens).

| Where | What is missing or confusing | Evidence |
|---|---|---|
| Opening | One sentence on what the package does; nothing on what you must supply, what you get back, or how long a real run takes. | `cleanTMLE.Rmd:21-29` |
| Motivation before code | Two sections of explanation come before the first code chunk, so the first output comes late. | `cleanTMLE.Rmd:21-61`, first chunk at `:62` |
| Time zero | `time_zero` is a required argument but is not defined in this vignette; the definition is in the second vignette. | `cleanTMLE.Rmd:87`, `plan-and-dossier.Rmd:43-59` |
| One candidate | The plan lists a single candidate, so the `selected` column shows no choice. The package default crosses two libraries with two truncation levels. | `cleanTMLE.Rmd:90`, `R/plan.R:169-171` |
| Dictionary | Used at the lock and described in one sentence. Neither vignette lists its nine field names or eight role values (only `reason` and `sensitivity` are named), and the second vignette prints five of the nine columns. A newcomer cannot write one for their own data. | `cleanTMLE.Rmd:129-139`, `R/dictionary.R:8-12`, `plan-and-dossier.Rmd:208` |
| Unexplained arguments | `roles`, `seed`, `K`, `reps` and `max_reps` appear before they are explained; `roles` is explained only in the second vignette. | `cleanTMLE.Rmd:94-95`, `plan-and-dossier.Rmd:183-187` |
| Unreadable print | The identification table prints as a wide data frame. The flag and response columns wrap onto separate blocks and are empty in this example. | `cleanTMLE.Rmd:295`; rendered `cleanTMLE.html` after "The dossier also carries the identification table" |
| No resolved run | Every check is unresolved, and the vignette approves anyway "to show the next step". The reader never sees a resolved run or learns how many repetitions it takes. | `cleanTMLE.Rmd:235-239`, `:320-325` |
| No refusal shown | The lock's refusals are where a newcomer will stall, and the vignette shows none. | `cleanTMLE.Rmd:134-139` |

Three statements go further than the evidence (claims lens).

The first is "Keeping the outcome away from the design solves the first problem" and "solves the second" (`cleanTMLE.Rmd:42-44`). The vignette itself says later that the lock cannot stop anyone from analyzing the outcome elsewhere (`:144-146`), and the feasibility check holds only over the simulated outcome models, so neither problem is solved. The second is "which tells you, before seeing any outcome, how precise the answer will be" (`:232-233`). It is a forecast that depends on the simulated outcome models resembling the real one, and the paper's fourth simulation claim, which tests it, has not been run. The third is "A real analysis keeps the defaults" (`:119`). The defaults are under review (the learner-defaults item in `docs/superpowers/specs/2026-10-18-outcome-blind-design-stage.md:131-149`), the final-cohort design stage selected main-terms logistic regression over the default super learner (`writing/digest.md:23`), and tolerances are the exception that the second vignette tells the reader to set from the study question (`plan-and-dossier.Rmd:20-22`).

The family boundary is also missing where it matters. "Feasible only over this family" appears in the second vignette (`plan-and-dossier.Rmd:111-112`) but not where `simulate_design()` is introduced (`cleanTMLE.Rmd:209-233`).

The uncommitted edit (`git diff` on `cleanTMLE.Rmd`) adds the identification table after the report object is printed and a closing section, "What this version does not do". Both are useful and accurate. Two problems: the new paragraph at `:290-299` repeats what the next paragraph says at `:302-307` (the rendered report "also shows the plan, an identification section" listing time zero and the timing of each variable), and the closing section compresses four limitations into one paragraph (`:367-377`). The complete-case warning is the one most likely to change a user's analysis, and it sits in the middle of that paragraph. It should lead. The warning that the table lists what was stated and does not show identification now appears in this vignette (`:298-300`), in the second vignette (`plan-and-dossier.Rmd:238-241`), in the methods (`methods.qmd:45`) and in the release notes (`cleanTMLE/NEWS.md`). Say it once per document.

The writing check returned two hits on a word its word list flags, the verb meaning to state a plan (`cleanTMLE.Rmd:299`, `:304`; the Step 1 heading at `:78` uses it too), and a typed count of patients (`:59`, which should come from `nrow(cr_example$design)`). Otherwise it returned only code identifiers in backticks and glossary terms. "State" or "write" would do in place of the flagged verb. The backticked identifiers are the package's own argument and function names and are expected in a package vignette. The glossary hits are real findings and are in the vocabulary table in section 9.

## 2. Second vignette, `cleanTMLE/vignettes/plan-and-dossier.Rmd`

The structure is right: plan choices first, each with the reason to choose it, then a reading guide to the report by question. "Choose each part with the study question in front of you, not the data" (`:20-22`) is the best sentence in either vignette.

The gaps are specific.

`candidates` is described only in words (`:117-125`), with no code. The reader cannot see how to list two libraries and two truncation levels, which is exactly what the Rescue.Co plans did. A three-line example belongs here.

The dictionary paragraph (`:189-205`) contains a sentence that cannot be parsed on first reading: "A stand-in covariate, an outcome-observation indicator, a follow-up negative control and a restriction variable fixed after time zero, whatever its role, each need a written `reason` in their row, and a stand-in also needs a `sensitivity`". The fields and roles are not listed. Section 9, pair E, proposes a rewrite.

The heading "Identification" (`:230`) promises more than the section gives, and the section then spends its length saying that it is not a verdict (`:238-241`). The table's content is the target trial's specification of time zero, eligibility and variable timing. The methods already say the table follows the components of a target trial [@cashinTransparentReportingObservational2025] (`methods.qmd:44`). A heading such as "Time zero and the timing of each variable" says what is in the section and needs no disclaimer.

The reading guide says to "ask for more repetitions before approving" when cells are unresolved (`:216-221`). In the Rescue.Co final cohort every cell was still unresolved at the cap (`final_operators_sim_main_terms_ate_cells` and the Results at `proposed-results-2026-10-21.qmd:67`). The guide says nothing about what to do when more repetitions are not affordable, or when no estimand is feasible. Both can happen: the cap was reached in the Rescue.Co final cohort, and no estimand was feasible in the first four months (`proposed-results-2026-10-21.qmd:55`). The second is a finding, not a failure. The guide also does not warn against the obvious temptation, which is to re-lock with looser tolerances until something passes. That is adjusting the design until the claim survives, and a vignette is the place to say so. This is a policy statement, so it is a decision for Andrew (recommendation 8).

Two smaller points. `approved_by` is checked only for being non-empty text (`R/unblind.R:85-86`), so it is a label and not an authentication; the vignette should say that, given the demonstration approves its own report. And the bias advice, "set it well below the smallest effect that would matter clinically" (`:75-77`), differs from the case-study rule of at most half the smallest important difference (`docs/superpowers/specs/2026-10-18-outcome-blind-design-stage.md:67`) and from `methods.qmd:76`. Pick one wording for all three.

The writing check hits were two uses of the same flagged verb (`plan-and-dossier.Rmd:204`, `:233`; the heading at `:230` has it as well) and backticked identifiers.

## 3. Rendered 0.5.0 vignettes, `docs/vignettes-0.5.0/*.html`

The rendered copies contain the uncommitted identification chunk, so they were built from the working tree and not from the merged commit (inference: git shows the chunk as uncommitted while the rendered text includes it). The vignette branch and the rendered files should be reconciled before the release is tagged.

The printed output reads well. The first run ends with the trimmed ATE as primary, a `FLAG` verdict from the negative control whose interval is "much wider than the null band", an estimate and a true value that the interval contains, and a decision log. The tone is honest.

Four defects show only in the rendered files.

The citation to van der Laan and Rubin prints as "Laan and Rubin 2006", and the reference list sorts it under "Laan" (rendered `cleanTMLE.html`, introduction and References). The vignette bibliography entry gives the author as `van der Laan, Mark J.` without braces; the project bibliography has `{van der Laan}`. The Shaw citation prints as a 2025 preprint, while the project bibliography holds the Statistics in Medicine 2026 version (`shaw2026`), which the manuscript cites. The identification table prints with empty flag and response columns, wrapped across three blocks. The example dictionary prints "enrolment form", a British spelling that the American-English rule asks to change in the package data.

## 4. Introduction, `manuscript/introduction.qmd`

Against the introduction checklist: the problem is stated, what is known is organized by claim, the three background paragraphs end on what is unknown, the gap is specific, no "first" or "no study" claim appears, no results appear, and the aim type is stated ("The aims are methodological", `:31`). The overlap check found no reuse from the listed sources. That is a solid draft.

Seven things need work.

First, five citation placeholders remain: Song 2026 (`:8`), Liu 2025 (`:15`), Karim 2026 (`:22`), Barnard 2026 and Conover 2025 (`:23`). One resolves now, because `song2026` is in the bibliography. I could not verify "48% of studies" and "most of the changes favoured outcomes with P values below 0.05", because the bibliography entry has no abstract, so check both against the paper before replacing the placeholder, and decide how an external number is recorded (the check flags it as typed). The other four are not in `references.bib`. The literature searches in `writing/lit/` found Barnard, Huling and Wolfson 2026 (a framework for causal estimand selection under positivity violations, `writing/lit/liu_estimand2.jsonl`), Conover and colleagues 2025 (objective study validity diagnostics, interpreted before unblinding, `writing/lit/conover_gate.jsonl`), Karim 2026 (`writing/lit/karim_plasmode.jsonl`) and Liu and colleagues 2025 (a tutorial on propensity score weighting under positivity violations, `writing/lit/liu_estimand2.jsonl`). I read only their abstracts in those logs. They are the nearest prior work to the paper's claim (section 8), so they must be read and added through sci-cite before the gap is stated.

Second, the gap sentence at `:10-11` understates staging. The staging paper's abstract says staging assesses "the population size available and potential bias" before the comparative analysis, and its example uses covariate balance and negative controls [@muntner2024]. The gap is narrower and more defensible: staging and the roadmap do not test whether an estimator can deliver the estimand under poor overlap. Section 9, pair I, proposes wording.

Third, `:14` is a 75-word sentence that says the ATO "needs no overlap at all" and that the ATT "needs it only for the treated". Both are imprecise. The ATO needs no strict positivity, but it still needs some patients with a chance of either treatment; the ATT needs untreated patients who resemble the treated. The methods have the same compression ("and nothing for the ATO", `methods.qmd:28`). Pair C rewrites it.

Fourth, the contribution at `:27-31` is a proposal followed by a 90-word objectives sentence (`:29`). It does not say what the proposal joins, what the reader gains, what the evidence is, or where it stops. Its one sentence of method (`:28`) runs to sixty-six words. Pair B rewrites it, and section 8 explains why.

Fifth, TMLE is never named as a reason for the design. The strongest reason to put TMLE inside a clean room is in the simulation proposal (`simulation/PROPOSAL-tmle-comparison-2026-10-20.md:11`): a flexible estimator has library, tuning and truncation choices that an analyst who sees the outcome could steer, so it is defensible only when fixed in advance. The introduction never says this. TMLE appears only inside the package name and is never defined.

Sixth, three sources in the bibliography bear directly on "whether the data can deliver the estimand" and are not used for that point: Walker's feasibility tool is cited only for equipoise (`:17`), `gattoStructuredProcessIdentify2022` (a data feasibility assessment framework) and `desaiProcessGuideInferential2024` (a process guide) are not cited, and `dang2023`, the roadmap paper itself, is not cited although the roadmap is named in `:9`. `gruber2022` (a targeted learning statistical analysis plan) is the closest earlier work on prespecifying a TMLE analysis. Also check whether `williamsonApplicationCausalRoadmap2023` supports "simulations ... used in the causal roadmap to choose the estimator" (`:21`): its abstract in the bibliography does not mention simulation (inference from the abstract only).

Seventh, style: "read the design through the propensity score" (`:24`) is figurative where the style rules ask for literal words; "favoured" and "analysed" are British (the American-English rule); the median sentence is about a quarter longer than the twenty-seven words measured in `style/voice.md:6`, and more than a quarter of sentences run over fifty words (my rough count, with inline code counted as one word).

The check returned five placeholder errors (expected until sci-cite resolves them), a typed-number warning for 48% and 0.05, and the missing overlap source.

## 5. Methods, `manuscript/methods.qmd`

Against the methods checklist the draft is strong on time zero, the estimand's attributes for the protocol's comparison, the identification statement, the plasmode construction, the feasibility rule with its Monte Carlo handling, negative controls, the deviations list, and honest disclosure (the pilot-driven overlap levels at `:134`, the arbitrary coverage tolerance at `:77`, the favorable baseline risk at `:153`, the analyst's prior view of outcomes at `:200`). Findings follow, grouped by lens.

### Methods lens

The estimand sentence at `:176-177` describes the protocol's comparison ("injured patients for whom a transport choice was made"), but the comparison that carries the confirmatory claim is Rescue.Co against other ambulance operators among direct arrivals. For that comparison the target population, the contrast and the identifying condition all change: the registry identifies the effect only if the choice of operator does not change survival to arrival. The estimand should be stated for the confirmatory comparison (pair H).

Covariate missing data are not described. The re-execution plan sets rules for dropping, imputing and pooling covariates (`docs/superpowers/specs/2026-10-20-rescueco-reexecution-0.5.0.md:27`), and none reaches the methods (checklist item on missing data, and the methods rule on missing data). Multiplicity is not mentioned although there are four comparisons and two outcomes; saying that only the confirmatory comparison carries a confirmatory claim and the rest are descriptive would suffice (the methods rule on multiplicity). The sensitivity analyses that exclude probable dispatch matches and walk-ins (`:208`) do not name the assumption they test (the methods rule on sensitivity analyses). The risk ratio reported in the results (`proposed-results-2026-10-21.qmd:78`) has no description in the methods.

Covariates are listed (`:188`) but without the reason for each or the diagram, which is still `Figure S[X]` (`:191`).

The tolerance and null band rules are incompletely stated. `:198` says the null band "was set from the controls' prevalence". The plan's rule is a quarter of the median prevalence with a floor (`docs/superpowers/specs/2026-10-20-rescueco-reexecution-0.5.0.md:79`), and the registry shows the same band for both runs (`final_operators_nc_null_band_upper`, `first_four_months_operators_nc_null_band_upper`), which suggests the floor applied and the prevalence term never did (inference). None of the three estimable controls had an interval inside the band (`final_operators_nc_household_urban_ate`, `final_operators_nc_cooking_wood_charcoal_ate`, `final_operators_nc_hypertension_ate`); their intervals are about as wide as the band itself, so a `GO` verdict needed a narrower interval than this cohort could give (inference from those intervals). The methods should state the rule with its floor, and the results should say that the controls were close to zero but too imprecise to exclude bias as large as the band.

The simulation study covers sample sizes of 500 and 2,000 (`:140`). The first-four-months cohort (`first_four_months_operators_n`) is smaller than both, and it is the case-study result that most depends on the simulation (nothing feasible). The final operators cohort (`final_operators_n`) is inside the range. The verdict's reliability at the smaller size is untested, and the paper should say so or add a scenario. The simulation grid uses main-terms logistic regression for both nuisance functions, with the super learner on a subset (`:158`), while the application plans use the super learner for the propensity score (`docs/superpowers/specs/2026-10-20-rescueco-reexecution-0.5.0.md:77`). The evidence for the verdict therefore comes from a different configuration from the application's.

### Claims lens

`:39`: "keeps the outcome out of every design choice" contradicts `:47`, where the lock "cannot prevent an analyst from analysing the outcome outside the package". Change "every design choice" to "the design stage" (pair H).

`:117-119` describes the estimation as a super learner, while the Results report that the candidate selected for the confirmatory comparison was main-terms logistic regression for the outcome regression (`proposed-results-2026-10-21.qmd:65`). The methods should say that each estimand is estimated with the candidate the design stage selected, that the super learner is one of the two candidate libraries and the propensity-score learner, and what the selected library was. The paragraph also depends on an open item: the neural-net decay grid sat on its edge in most simulated fits, and Andrew has asked for the grid and the treatment penalty to be fixed (`docs/superpowers/specs/2026-10-18-outcome-blind-design-stage.md:141-147`). Hold `:117-119` and the edge-share placeholder at `:119` until that item is decided.

`:174` says the confirmatory rule was "written into the registered plan" (the governance note removes "registered"; item `M2` of the proposed changes fixes it) and does not say who wrote the thresholds. The spec states the expected pick before the run (`docs/superpowers/specs/2026-10-18-outcome-blind-design-stage.md:44`), and the thresholds (80% of each arm, 0.25 standard deviations, five named covariates, the ATO allowed only between operators) were written after the unblinded analysis was known to the analyst. The application shows that the rule can be applied mechanically without the outcome. It does not show the rule would have picked the same comparison without that knowledge. Add one sentence (pair H).

`:201` says the time zero and variable timing checks were added after the unblinded analysis found these problems, and that they were guarded against over-fitting. That is a good disclosure, and it limits what the identification table can claim as a finding: the design stage flags these problems because the checks were built from them (section 8).

Typed design constants appear in the methods (tolerances, the band, 80%, 0.25, 22%, 0.75, the added deaths per 1,000 at `:180`, the dates). They are plan parameters and not results, but the writing rules ask that every number trace to the registry or to a cited document. Add them to the registry from the plan document, or cite the plan with its path and line.

### Reader lens

The median sentence is about two-fifths longer than the twenty-seven words in `style/voice.md:6`, and about a fifth of sentences run over fifty words (my rough count, with inline code counted as one word). The longest run to ninety-five, ninety-two, eighty-six, eighty-two and seventy-eight words (`:162`, `:180`, `:152`, `:188` and `:209`). These are hard to follow, and `voice.md` allows splitting when that is so. Terms that need a plain replacement are in section 9: "cell" (`:73`), "estimable in our sense" (`:75`), "stand-in" (`:42`, `:182`), `GO`, `STOP` and `FLAG` (`:104`), "re-executes" (`:200`), and "population gate" (introduced by the proposed text at `proposed-changes-2026-10-21.md:11`). British spellings occur at `:39`, `:47`, `:89`, `:102`, `:120`, `:167`, `:177`, `:184`, `:188`, `:193` and `:209`. "Rescue.co" is mistyped at `:208` (the proposed change `M6` fixes it).

The check returned typed-number errors for three simulation constants (the sample size, the number of datasets and the number of simulated repetitions), warnings for typed constants that match the registry, undefined all-capital words (the three verdict labels, the registration placeholder and the archive placeholder), and one flagged word, `emulation`. The flagged word is a standard term in the target trial literature; I recommend adding it to `allow_words` in `writing/config.yml` (Andrew's call, since I do not edit that file).

## 6. Proposed results, `writing/proposed-results-2026-10-21.qmd`

The numbers are sound. All keys exist. The overlap grades and the population-criterion outcomes match the exported tables: poor for the two versions of the protocol's comparison and for Rescue.Co against arrival without an ambulance, moderate for any ambulance against arrival without one, good between operators; no estimand below the ATE passes in the three poor comparisons; only the trimmed ATE passes for any ambulance; all three pass between operators (`population_gate.csv`). The text is also candid that every cell stayed unresolved at the cap (`:67`) and that the negative controls were inconclusive.

Five things need work.

The paragraphs restate Table 1 in prose. `:32`, `:42-43` and `:46-48` carry more than fifty inline numbers between them, in sentences of seventy to eighty-five words. The global rules say to state each thing once and ask for the main text to be a plain report with detail in the supplement. Pairs F and G rewrite two of them.

`:35-36` give twenty-four numbers for loss to follow-up. The message is that about a quarter of each arm lacked six-month status and the arms differ by a few percentage points (`final_operators_treated_six_month_not_recorded`, `final_operators_comparator_six_month_not_recorded`). Put the full set in a supplementary table.

`:46` builds a range from two hand-picked keys (`final_protocol_direct_ato_max_shift` to `final_protocol_all_arrivals_trimmed_ate_max_shift`). It is correct today (I checked it against `population_gate.csv`), but it will silently stop being correct if any underlying value changes. Add `_min` and `_max` registry keys for ranges, as elsewhere.

The most important sentence is buried. The design stage found that the protocol's own comparison, on which the study was powered, has poor overlap and that no estimand below the ATE passes the population criterion, and the confirmatory claim therefore moved to a different, narrower question: Rescue.Co against other ambulance operators. That is the finding a reader should take away, and it appears only as the reason a rule picked the third comparison (`:74`). It needs a topic sentence of its own. It should be written as "poor overlap and no acceptable estimand below the ATE" and not as "infeasible" until the protocol comparison's simulation has finished (`writing/digest.md:29` lists it as pending).

The two simulation-based results point in opposite directions, and the Results give no sign that this is expected. At the first four months the main-terms candidates were biased with variance estimates wrong in both directions while the super learner met the tolerances but was too wide (`:57-59`). In the final cohort the ranking reversed by a smaller margin: the main-terms candidates met the tolerances, and the super learner's coverage fell to the low end (`final_operators_sim_super_learner_ate_coverage_min`, `:65-66`). The Discussion has to explain this, and the Results should at least keep the two paragraphs parallel so that the reader can compare them.

Smaller points: "four eligibility criteria" (`:38`) is true for four of the five locks and three for the protocol's comparison over all arrivals (`rescueCo/results/reexec_050/final/design/L1_identification.csv`), and the outline says three (`writing/outline-results.md:19`); the count is typed as a word and has no registry key. The simulation section is entirely placeholders (`:17-25`), so the proposed results are, for now, a design-stage report. If the simulation results arrive in the same proportion, the application will outweigh the evidence for the paper's four claims.

The check returned no errors and two undefined-acronym warnings (the two verdict labels are first used here and defined in the methods).

## 7. Proposed changes memo, `writing/proposed-changes-2026-10-21.md`

The form is right: before and after text, numbers by registry key with the printed values for checking, and the open question kept for Andrew. Line numbers in the memo match the current `methods.qmd` (I checked `:171`, `:174`, `:196`, `:199`, `:202`, `:204` and `:208`).

Item `M3` needs a correction. It says the lock hashes "were kept with the SAP rather than registered publicly", but the governance note says the lock and report hashes are reported in the supplement (`docs/superpowers/specs/2026-10-18-outcome-blind-design-stage.md:129`). `M3` also uses `demonstrates` (a flagged word), and neither `M3` nor `M4` says plainly what the same note requires, that every approver had seen the unblinded results; `M4` relies on `methods.qmd:200`, which says the analyst had. State it once, in `M3`.

Item `M2` keeps "population gate" as a label and runs to more than a hundred words in three sentences. Replace the label with a description of the criterion (pair G uses "population criterion" only for the first mention and then describes it).

The memo omits changes that this audit finds: the estimand for the confirmatory comparison, the null band rule, the estimation paragraph against the selected candidate, covariate missing data, multiplicity, and the provenance of the rule's thresholds (section 5). The bold lead-ins and line numbers the check flags are normal in a memo and need no action.

## 8. The contribution as an advance

### What is claimed now

The brief's message is that an outcome-blind plasmode simulation can tell which estimands a design can deliver, more accurately than overlap diagnostics, and that choosing the estimand at that stage "keeps inference valid where choosing it after seeing outcomes does not" (`writing/brief.md:14`). The introduction proposes "choosing the estimand at the design stage with an outcome-blind plasmode simulation, inside a staged workflow" (`introduction.qmd:27`). The package `README` says it "separates the design of an observational analysis from its outcome" (`cleanTMLE/README.md:3`).

Andrew's claim is different and larger: that running the causal roadmap with TMLE inside a clean room, with the estimand and the estimator both chosen before the outcome is seen, is a scientific advance. That is a claim about a workflow. The drafts present a claim about one component of it.

### Do the drafts make that case?

Not yet, for three reasons.

The statement of contribution is narrower than the claim. The causal roadmap, TMLE and the clean room appear in the introduction only as prior work, not as the thing that is joined. The estimator half of the claim, and the fact that the estimator applied after unblinding is the one the simulation evaluated (`methods.qmd:71`), is stated once in the methods and not in the introduction. That property is the clearest difference from outcome-blind simulation as the roadmap describes it [@nanceCausalRoadmapSimulations2024a], and it is undersold.

The evidence for the four simulation claims does not exist yet. The registry holds several hundred keys and every one begins with `final` or `first_four_months`; none concerns the simulation study. The Results hold only placeholders for it (`proposed-results-2026-10-21.qmd:17-25`), and the TMLE-against-weighting comparison that would speak to the TMLE half of the claim is a proposal awaiting Andrew's decisions (`simulation/PROPOSAL-tmle-comparison-2026-10-20.md:66-70`).

The nearest prior work is missing from the bibliography, so the gap cannot yet be stated with confidence. Barnard and colleagues propose a design-based procedure for choosing an estimand under positivity violations, Conover and colleagues propose gating unblinding on pre-specified diagnostics, and Liu and colleagues give guidance on choosing a primary estimand. If the first two are as close as their abstracts suggest, "no pre-outcome check judges estimability itself" (`writing/outline-introduction.md:22`) will need to become "no pre-outcome check tests whether an estimator delivers the estimand", which is narrower.

### What each part of the case rests on

| Part of the case | Prior work to position against | Evidence in the drafts | Status |
|---|---|---|---|
| Pre-specification limits analytic flexibility, and outcome switching is common | [@kahanHowDesignPrespecified2020]; [@rubin2008]; [@song2026] | Literature only | Usable once the Song figures are verified |
| The estimator can be chosen before the outcome by outcome-blind simulation | [@nanceCausalRoadmapSimulations2024a]; [@nanceApplyingCausalRoadmap2023]; [@gruber2022] (a targeted learning statistical analysis plan) | The estimator evaluated is the estimator applied (`methods.qmd:71`); the forecast of interval width (claim 4, not run); the final-cohort selection, which rests on unresolved cells (`proposed-results-2026-10-21.qmd:63-68`) | Partly |
| The estimand can be judged before the outcome | [@petersenDiagnosingRespondingViolations2012]; [@walker2013]; [@gattoStructuredProcessIdentify2022]; [cite: Barnard 2026 estimand selection under positivity violations]; [cite: Liu 2025 weighting under positivity violations] | Claims 1 and 2 (not run) | Nothing yet |
| Choosing the estimand before the outcome preserves inference | [@kahanHowDesignPrespecified2020]; [@song2026]; [cite: Conover 2025 objective study validity diagnostics] | Claim 3 (supplement, not run); close to a known result, so present it as quantification | Nothing yet |
| TMLE belongs inside a clean room, because a flexible estimator is defensible only if fixed before the outcome | [@gruber2012]; [@vanderlaanSuperLearner2007]; [@zheng2011]; [@wiederkehr2026] | The proposed comparison (not decided); the case study has no known truth, and in the final cohort the design stage chose the simpler candidate | Argued, not tested |
| A software-enforced clean room is tamper-evident | [@muntner2024] (its audit trail); [@desaiProcessGuideInferential2024] | Package behavior; no evidence of bias prevented | Argued, with a stated boundary (`methods.qmd:47`) |
| The workflow runs on a real registry and records what a reader needs | [@dang2023]; [@cashinTransparentReportingObservational2025]; [@hernan2016] | Design-stage results for the confirmatory comparison and the first four months; other comparisons pending (`writing/digest.md:29`) | Demonstration only |

### Boundaries the first page must carry

The simulation's verdict holds over the family of simulated outcome models and the estimated treatment mechanism, not for any true outcome model (`methods.qmd:65`). The lock makes changes detectable and does not stop analysis of the outcome elsewhere (`:47`). Nothing in the workflow corrects for unmeasured confounding, and a `FLAG` verdict is not evidence of its absence (`:107`; `writing/brief.md:30`). The case study is a re-analysis in which the analyst had seen the outcomes of the same registry (`:200`), and its decision rule and its time zero checks were written after those results were known. The paper covers version 0.5.0 only, with complete case as the primary handling of missing outcomes. The simulation covers two sample sizes. The introduction carries none of these today.

### Where the drafts overclaim and undersell

Overclaim. The brief's message ("keeps inference valid", `writing/brief.md:14`) is absolute, while the evidence will be coverage within Monte Carlo error in the scenarios simulated. Brief claim 5 ("chose the confirmatory comparison, the estimand and the estimator by rules fixed in advance", `:21`) is true mechanically, but the rules were fixed by people who knew the unblinded results, and the final-cohort estimator choice rests on cells that never resolved. The methods (`:39`) and the vignette (`cleanTMLE.Rmd:42-44`) overclaim in the ways listed in the first sections. The package `README` and `DESCRIPTION` promise to learn "which estimands the design can deliver" without the family boundary.

Undersell. The closed loop between the estimator evaluated and the estimator run. The reporting of Monte Carlo uncertainty in the design stage's own verdict (the "unresolved" cells), which the abstracts I read for the prior work do not mention (inference). The refusals taken from the target trial, such as time zero, eligibility, and covariates fixed after treatment. The first-four-months result, where the overlap diagnostics gave no warning (`first_four_months_operators_n_outside` and `first_four_months_operators_ate_max_weight`: no patient outside the band and a small largest weight in the operators comparison) and yet no estimand was feasible. This is the single place in the application where the simulation tells the analyst something that overlap diagnostics did not, and it is not framed as such. And the final-cohort design stage selecting a simpler candidate than the default, which shows the workflow does not push toward flexible machinery. Do not spin that last point until the learner-defaults item has ruled out a package defect, because a defect would make the selection partly an artifact.

### What the case study can and cannot carry

It can show that the workflow runs on real data, what it records, and what it said at two time points. It cannot show that staging prevents bias, because the analyst had seen the outcomes. It cannot show that the design stage would have found the identification problems blind, because the checks were added after those problems were found (`methods.qmd:201`); it shows that the package now records them. It cannot show TMLE's advantage, because there is no known truth. Brief claim 5 should be rewritten to this scope.

### If results disagree with the claim

The project rule applies: a result that contradicts a claim means rescoping the claim, never adjusting the design until the claim survives. The simulation plan already says what to do if a best-threshold overlap diagnostic matches the simulation (`simulation/PLAN.md:51`). Write that fallback paragraph now, before results exist: the advance would then rest on what diagnostics cannot give (the estimator selection, the forecast of interval width, the record), and the claim of better agreement would be dropped. The same applies to the weighting-and-matching comparison, which the proposal already says loses a subclaim, and not the design, if a contrast is unresolved (`simulation/PROPOSAL-tmle-comparison-2026-10-20.md:64`).

### A conclusion that fits the claim

The Discussion's conclusion should answer the objective in one or two sentences and hold under the stated assumptions. A skeleton that cannot outrun the evidence, to be completed only from the registry once the simulation has run:

> An outcome-blind plasmode simulation inside a locked workflow agreed with estimability on the true outcome model [R: more or less often] than overlap diagnostics when the true model lay inside or near the simulated family, and [R: result] when it lay outside. Choosing the estimand before the outcome kept coverage [R: result]; choosing it afterward [R: result]. In the Rescue.Co registry the design stage [R: what it recommended and what it could not decide]. These results hold over the simulated family and the sample sizes studied, and they do not address unmeasured confounding.

## 9. Before and after proposals

All proposed text uses only registry keys, keys in `writing/references.bib`, or the placeholders the drafts already use. Where a vignette needs a citation the vignette bibliography lacks, I name the project-bibliography key to copy. "Before" blocks are verbatim.

### A. Opening of the main vignette

Before (`cleanTMLE.Rmd:21-29`):

```text
cleanTMLE runs an observational comparison in stages, so that the design is
settled before anyone sees the outcome. This is the clean-room idea of
@muntner2024staging. The package adds one question to the design stage:
can these data deliver the estimand we planned, and if not, which one can
they deliver?

In this vignette we take one simulated cohort through every stage. It takes
a few minutes to build. Each step says what it decides and what to look for
in its output.
```

After:

> cleanTMLE helps you settle the design of an observational comparison before anyone looks at the outcome. You write a plan, lock the design data without the outcome, and run a simulation on the real covariates with made-up outcomes to see which estimands the data can support and which targeted maximum likelihood estimator (TMLE) to use for each. After a review, you add the outcome and estimate the effect with the estimator the simulation evaluated. The staging follows the clean-room approach of @muntner2024staging, and the estimation builds on the tmle package [@gruber2012].
>
> The simulation answers one question: can these data deliver the estimand we planned, and if not, which one can they deliver? It answers over the family of outcome models you write into the plan. It does not show that confounding is controlled, and the lock makes changes detectable without stopping anyone from analyzing the outcome elsewhere.
>
> You need three things: design data (an id, the treatment, the covariates and any negative control outcomes, with no outcome), a second data frame holding the outcome, and a variable dictionary that says when each design column was recorded and when the characteristic it describes was fixed. `cr_example` supplies all three for a simulated cohort of `r nrow(cr_example$design)` patients. The vignette takes a few minutes to run. A real plan with the default settings takes far longer, and `simulate_design()` prints an estimate after its first repetition. Each step says what it decides and what to look for in its output.

Then add a table of the steps (the function reference of `causalRisk` groups its functions the same way):

> | Step | Function | What it decides |
> |---|---|---|
> | 1 | `analysis_plan()` | Estimands in order of preference, tolerances, candidates, negative controls |
> | 2 | `create_analysis_lock()` | Fixes the plan, the design data and the dictionary under one hash |
> | 3 | `assess_design()` | How good the overlap is, and whom each estimand describes |
> | 4 | `simulate_design()` | Which estimands are feasible, with which candidate |
> | 5 | `negative_control_ladder()` | Whether a negative control shows residual confounding |
> | 6 | `design_report()`, `unblind()` | The reviewed decision, then admission of the outcome |
> | 7 | `estimate_effect()` | The estimate, with the candidate the simulation selected |
> | 8 | `export_design_log()` | The record of every decision |

Also add one sentence to "Why decide before seeing the outcome" tying the steps to the roadmap, whose steps the Nance abstract names as defining the quantity of interest, evaluating the assumptions, conducting estimation and interpreting results [@nanceCausalRoadmapSimulations2024a] (copy that key into the vignette bibliography):

> The steps follow the causal roadmap, which asks you to define the quantity of interest, evaluate the assumptions it needs, prespecify the estimation, and interpret the results [@nanceCausalRoadmapSimulations2024a]. Steps 1 and 2 define the question and fix the plan; steps 3 to 5 check the assumptions and the estimators before any outcome is seen; steps 6 and 7 admit the outcome and estimate with the estimator that was tested; step 8 records each decision.

### B. Statement of contribution in the introduction

Before (`introduction.qmd:27-31`):

```text
Here, we propose choosing the estimand at the design stage with an outcome-blind plasmode simulation, inside a staged workflow implemented in the R package cleanTMLE.
The simulation draws new treatments from the propensity score estimated on the study's covariates and outcomes from a prespecified family of outcome models, judges an estimand feasible when some candidate estimator meets prespecified tolerances for bias and coverage on every model in the family, and makes the first feasible estimand in a prespecified order the primary estimand, all before the outcome is joined to the data.
Our objectives were to estimate, in a simulation study, the sensitivity and specificity for estimability on the true outcome model of this verdict, of the same simulation with each patient's observed treatment kept, and of overlap diagnostics, by estimand and by how closely the true outcome model resembles the simulated family; the coverage and type I error of the estimate when the estimand is chosen before rather than after the outcome is seen; and the agreement between the confidence-interval half-width forecast before the outcome and the half-width realized after it.
We also show the workflow in a re-execution of a trauma registry analysis from Nairobi, Kenya, with the outcome withheld from the software.
The aims are methodological; the application estimates a causal effect, whose estimand and identifying assumptions are given in the Methods.
```

After:

> Here, we propose that the estimand and the estimator both be chosen before the outcome is seen, from an outcome-blind plasmode simulation run inside a staged, locked workflow.
> The analyst fixes a plan and the design data; the simulation judges which estimands the design can deliver and which targeted maximum likelihood estimation (TMLE) candidate to use for each; a review precedes outcome access; and the estimator applied afterward is the one the simulation evaluated.
> The workflow joins the causal roadmap's use of outcome-blind simulation to prespecify an estimator [@nanceCausalRoadmapSimulations2024a] with the clean room's restriction of access to the outcome [@muntner2024]. It extends the first from choosing the estimator to judging whether the estimand can be delivered. Fixing a flexible estimator in advance matters because its library, tuning and truncation are choices that an analyst who sees the outcome could steer [@kahanHowDesignPrespecified2020].
> We implemented the workflow in the R package cleanTMLE, which builds on the tmle package [@gruber2012].
> Our objectives were to estimate, in a simulation study, how often the simulation's verdict agrees with whether some candidate meets the prespecified tolerances for bias and coverage on the true outcome model, and to compare that agreement with the agreement of overlap diagnostics and of the same simulation with each patient's observed treatment kept.
> We also estimated the coverage and type I error of the estimate when the estimand is chosen before rather than after the outcome is seen, and the agreement between the confidence interval half-width forecast before the outcome and the half-width realized after it.
> To show the workflow on real data, we re-analyzed a trauma registry from Nairobi, Kenya, with the outcome withheld from the software.
> The aims are methodological. The verdict holds only over the family of simulated outcome models, the lock makes changes detectable without preventing analysis of the outcome elsewhere, and neither removes unmeasured confounding. The application estimates a causal effect, whose estimand and identifying assumptions are given in the Methods.

Notes: the word "advance" is deliberately absent. The statement says what is joined, what is extended and where it stops, and the Discussion's conclusion should carry the claim once results exist. Name Barnard, Conover and Liu in the sentence on what is extended only after those papers are read; until then that sentence is provisional (section 8). "Re-analyzed" replaces "re-execution", which is project jargon.

### C. The estimands and positivity (`introduction.qmd:14`)

Before:

```text
The average treatment effect (ATE) needs every patient to have some chance of receiving either treatment; the ATE in patients whose propensity score lies within a band, the trimmed ATE, needs this only inside the band [@crump2009]; the average treatment effect among the treated (ATT) needs it only for the treated; and the overlap-weighted effect (ATO) weights each patient by the variance of treatment given the covariates, so it needs no overlap at all [@li2017].
```

After:

> The average treatment effect (ATE) needs every patient to have some chance of receiving either treatment. The trimmed ATE, the effect in patients whose propensity score lies within a band, needs this only inside the band [@crump2009]. The average treatment effect among the treated (ATT) needs untreated patients who resemble the treated. The overlap-weighted effect (ATO) weights each patient by the variance of treatment given the covariates, which is zero where treatment is certain, so it does not require strict positivity but gives most weight to the patients whose treatment was least predictable [@li2017].

The same correction applies to `methods.qmd:28`: replace "and nothing for the ATO" with "and no strict positivity for the ATO, whose weights vanish where treatment is deterministic".

### D. Overclaims in the main vignette (`cleanTMLE.Rmd:42-44`, `:232-233`)

Before:

```text
Keeping the outcome away from the design solves the first problem.
Checking, still without the outcome, which estimands the design can support
solves the second.
```

```text
`ci_halfwidth`
is the expected half-width of the 95% confidence interval, which tells you,
before seeing any outcome, how precise the answer will be.
```

After:

> Keeping the outcome out of the design stage limits the first problem, because no design choice can respond to the estimate. Checking, still without the outcome, which estimands the design can support tells you whether the second problem is present and which estimand the data can support.

> `ci_halfwidth` is the expected half-width of the 95% confidence interval. It forecasts, before any outcome is seen, how precise the answer will be if the simulated outcome models resemble the real one.

And at `:117-124`, replace "A real analysis keeps the defaults" with: "A real analysis starts from the defaults and lets the design stage choose between them and a simpler candidate. Tolerances are the exception: set them from the study question, as the second vignette explains."

### E. The variable dictionary (`plan-and-dossier.Rmd:189-205`)

Before (the unparseable sentence in context):

```text
The lock stops
on a covariate fixed after time zero, on a treatment not fixed at time zero,
and on a column from the follow-up source that is neither an
outcome-observation indicator nor a negative control. A stand-in covariate,
an outcome-observation indicator, a follow-up negative control and a
restriction variable fixed after time zero, whatever its role, each need a
written `reason` in their row, and a stand-in also needs a `sensitivity`, the
analysis that leaves it out.
```

After:

> The dictionary has one row for each column the lock keeps and nine fields: `column`, `variable`, `source`, `follow_up`, `recorded`, `fixed`, `role`, `reason` and `sensitivity`. `recorded` and `fixed` each take the values before, at or after time zero. Age is recorded after injury but fixed before it, and a heart rate on arrival is fixed after treatment is assigned. `follow_up` is `TRUE` for a column that exists only for patients reached at follow-up. `role` is one of `id`, `treatment`, `covariate`, `stand_in`, `after_treatment`, `observation`, `negative_control` and `restriction`.
>
> The lock stops on three things: a covariate fixed after time zero, a treatment not fixed at time zero, and a follow-up column that is neither an outcome-observation indicator nor a negative control. It asks for a written `reason` in four cases: a stand-in (a post-treatment covariate that stands in for a characteristic fixed before time zero), an outcome-observation indicator, a negative control taken from follow-up, and a restriction variable fixed after time zero. A stand-in also needs a `sensitivity`, which names the analysis that leaves it out. The lock checks that a response is present, not that it is adequate, and it cannot check that the timings you state are true.

Then show one refusal, so the reader sees what a stall looks like: a short chunk (`tryCatch`) that locks the example with `severity` set to fixed after time zero and prints the message.

### F. Loss to follow-up in the Results (`proposed-results-2026-10-21.qmd:32-36`)

Before: the three paragraphs at `:32`, `:35` and `:36` (about 175 words with thirty-five inline numbers; not repeated here).

After:

> The number of patients in each comparison, after restriction to patients whose six-month status was recorded, is given in @tbl-design.
> Among direct arrivals, six-month status was not recorded for `r reg$final_operators_treated_six_month_not_recorded$text` of the patients transported by Rescue.Co and `r reg$final_operators_comparator_six_month_not_recorded$text` of those transported by other operators, and disposition from the emergency department was not recorded for `r reg$final_operators_treated_disposition_not_recorded$text` and `r reg$final_operators_comparator_disposition_not_recorded$text`.
> The design stage saw six-month status only through the data steward's indicator, and the primary analysis uses the patients whose six-month status was recorded. The shares for the other comparisons are in Supplementary Table S[X].

### G. The population criterion in the Results (`proposed-results-2026-10-21.qmd:45-48`)

Before: the paragraphs at `:45-48` (an 83-word sentence with fourteen inline numbers, followed by two sentences carrying nine and three more).

After:

> The plan applies a population criterion to the primary estimand: an estimand below the ATE is acceptable only if it keeps at least the stated share of each arm and shifts none of five named covariates by more than the stated amount. Applied here to every estimand below the ATE, for description, it was met by none in the three comparisons with poor overlap, because trimming kept too little of the comparison arm, the weighted population differed from the full cohort on a named covariate, or both (Table 1). For any ambulance against arrival without an ambulance only the trimmed ATE met it, and between operators all three did.

The share and the amount are in the Methods (item `M2` of the proposed changes), so the Results need not repeat them. Also change the table caption and the cell label "Passing the gate" to "Meeting the criterion".

### H. Three sentences in the methods

First sentence, `methods.qmd:39`. Before: "Separating the design from the outcome keeps the outcome out of every design choice, and any later departure from the SAP, such as reporting an estimand other than the primary, is labelled as such in the results and the log." After: "Separating the design from the outcome keeps the outcome out of the design stage, and any later departure from the SAP, such as reporting an estimand other than the primary, is labeled as such in the results and the log."

Second sentence, `methods.qmd:174`, after the sentence stating the rule, add: "The thresholds in this rule were written after the unblinded analysis of the same registry was known to the analyst, so the application shows that the rule can be applied without the outcome, not that it would have chosen the same comparison without that knowledge."

Third sentence, `methods.qmd:176-177`. After the sentence on the protocol's comparison, add: "For the comparison that carried the confirmatory claim, the target population is injured patients who arrived directly by ambulance, the contrast is transport by Rescue.Co against transport by another operator, and the registry identifies the effect only if the operator does not change survival to arrival."

### I. The gap sentence (`introduction.qmd:10-11`)

Before: "These frameworks leave one decision that is hard to make in advance: whether the data can deliver the estimand the protocol asks for."

After: "Staging assesses the available population size, covariate balance and negative controls [@muntner2024], and the causal roadmap uses outcome-blind simulation to choose the estimator for an estimand that is already fixed [@nanceCausalRoadmapSimulations2024a]. Neither tests, before the outcome is seen, whether an estimator can deliver the estimand the protocol asks for."

Check the second sentence against Barnard and Conover once they are read; if either does test this, narrow it.

### Plain replacements for coined terms

The glossary (`writing/glossary.yml`) already maps several of these; the check enforces them. Backticks mark the term being replaced. Names in the package interface (arguments, functions, columns) stay as they are in the code; the change is to the prose.

| Where | Term | Proposed plain term |
|---|---|---|
| Vignettes (`cleanTMLE.Rmd:98`, `:230`, `:375`; `plan-and-dossier.Rmd:26`, `:39`) | `ladder`, `rung` | "the estimands in order of preference", "the next estimand in the list" |
| Vignettes (`cleanTMLE.Rmd:215`; `plan-and-dossier.Rmd:84-88`) | `outcome surfaces` | "simulated outcome models" (the `surfaces` argument keeps its name) |
| Vignettes (heading `cleanTMLE.Rmd:269`, and throughout) and spec | `dossier` | one term shared with the paper: the SAP in the paper is the same object, so pick either "design report" (matches `design_report()`) or SAP, and use it everywhere |
| `plan-and-dossier.Rmd:230` | the heading "Identification" with its subtitle | "Time zero and the timing of each variable" |
| `plan-and-dossier.Rmd:235` | `the nodes after treatment` | "variables recorded after treatment" |
| `methods.qmd:42`, `:182`; both vignettes | `stand-in` | "proxy covariate" or, on first use, "a post-treatment covariate used as a proxy for a characteristic fixed before time zero" |
| `methods.qmd:73` | `cell` | "combination of estimand, candidate and simulated outcome model" on first use, then "condition" (the term used in Morris and colleagues' simulation reporting [@morris2019]) |
| `methods.qmd:75`, `:128`, `:146` and the brief | `estimable in our sense`, `estimability` | "meets the tolerances on the true outcome model"; "estimable" already has a technical meaning for identifiable parameters |
| `methods.qmd:104`; vignettes | `GO`, `STOP`, `FLAG` | describe the verdict first ("every control inside the band", "any control outside", "otherwise inconclusive"), then give the label once as the package's output |
| `proposed-changes-2026-10-21.md:11`; Results `:45`, `:52`, `:86` | `population gate`, `Passing the gate` | "population criterion", "meeting the criterion" |
| Results `:38`, methods `:45` | `identification table` | "the table of time zero and variable timing" (the methods already tie it to the target trial components [@cashinTransparentReportingObservational2025]); the name `identification` stays in the code |
| `methods.qmd:200`; `introduction.qmd:30` | `re-executes`, `re-execution` | "re-analyzes", "re-analysis with the outcome withheld from the software" |
| Vignettes, `methods.qmd` | "design log" | keep, and say once that it is the audit trail of the clean room described by Muntner and colleagues [@muntner2024] |
| `introduction.qmd:24` | "read the design through the propensity score" | "use only the propensity score, the weights or negative controls" |
| `cleanTMLE.Rmd:194-195` | "An answer about the treated is an answer about sicker patients." | "The ATT therefore describes sicker patients than the cohort does." |
| `cleanTMLE.Rmd:172` | "The plot shows why." | State the finding: "Near a propensity score of 0 almost every patient is a control, and near 1 most are treated (plot below)." |
| Introduction, methods, results | `favoured`, `analysed`, `analysing`, `labelled`, `enrolment`, `enrols`, "Rescue.co" | "favored", "analyzed", "analyzing", "labeled", "enrollment", "enrolls", "Rescue.Co" (American-English rule) |

The prose has no hits for the usual stock vocabulary of machine-written text (the check found only the verb for stating a plan and one standard term, which `allow_words` can admit). The texture problem is structural: very long sentences, numbers stacked in a sentence, the same warning repeated, and a few aphorisms. Those are listed above.

## 10. Recommendations in priority order

| # | Recommendation | Effort | Risk | Decision Andrew must make |
|---|---|---|---|---|
| 1 | State the contribution as a bounded advance. Rewrite the brief's message and key claims, the introduction's contribution (pair B), the vignette's "why" section, and pre-write the Discussion's conclusion and the fallback paragraph for each simulation claim. | Medium: about a day of writing, no new analysis | If claims 1 to 4 do not hold, the statement must be rescoped; the fallback paragraph is written now so that nobody tunes the design | Headline: the integrated workflow (roadmap, TMLE, clean room, estimator evaluated equals estimator applied), with the simulation as evidence, or the current estimand-feasibility message. Whether the TMLE-against-weighting-and-matching comparison becomes a main claim. |
| 2 | Close the prior-work gap. Read and add Barnard 2026, Conover 2025, Liu 2025 and Karim 2026 through sci-cite; replace the Song placeholder with `song2026` after checking the figures; cite `dang2023`, `gruber2022`, `gattoStructuredProcessIdentify2022` and `desaiProcessGuideInferential2024` where they bear on the gap. | Medium | Reading Barnard and Conover may narrow the novelty claim | Accept a narrower gap if the papers require it, and say who reads them. |
| 3 | Rescope the case-study claims and disclose provenance: brief claim 5, `methods.qmd:174` (who wrote the thresholds), `:201` and `:200` (checks built after the problems were found; analyst had seen outcomes), and the statement of the estimand for the confirmatory comparison. | Low | Reads less dramatic, but is defensible | Approve wording for the second and third sentences in pair H and for the rewritten claim 5. |
| 4 | Make the main vignette runnable on a newcomer's own data: pair A, a "what you need" paragraph, the steps table, `time_zero` defined at first use, the dictionary fields and roles (pair E), two candidates in the example plan, a readable print of the identification table, one refusal shown. | Medium: half a day to write, one re-render (needs R, so after the design stage finishes) | Low | Keep the tutorial register ("you" and "we")? Move the motivation after the first run? |
| 5 | Align the methods with the selected candidate, the open learner-defaults item and the null band rule. Hold `methods.qmd:117-119` until the item is decided; state the band rule with its floor; state that the controls could not return a `GO` verdict at this precision. | Low to medium | The item may change the super learner, so wording written now could be wrong | Hold the estimation paragraph, or write it for the released 0.5.0 and revise after a patch. |
| 6 | Adopt one vocabulary across vignettes, spec and paper (section 9 table), changing prose and keeping the interface names. | Medium: many files | Renaming arguments or roles would be a breaking change; do not | The term set: SAP or "design report" for the dossier; "proxy covariate"; "meets the tolerances" for estimable; whether the `GO`, `STOP` and `FLAG` labels stay in the paper. |
| 7 | Fix the overclaims: pair D, `methods.qmd:39`, pair C, the brief's "keeps inference valid", and the `README`'s missing family boundary. | Low | Low | None beyond approving the wording. |
| 8 | Add guidance to the second vignette for unresolved and infeasible results, `approved_by` as a label, and one wording for the bias tolerance (at most half the smallest important difference, or "well below"). | Low | A policy statement is binding on users | Your position: is re-locking with different tolerances after a negative result allowed, and how is it recorded? What is the advice at the repetition cap? |
| 9 | Restructure the proposed results: tables carry numbers, a topic sentence states what the design stage found for the protocol's comparison, ranges come from `_min` and `_max` keys, and detail moves to the supplement (pairs F and G). | Medium | Low | What stays in the main text. |
| 10 | State the simulation's scope: sample sizes cover 500 and 2,000 and not the first-four-months cohort; the grid uses main-terms nuisance models while the application uses a super learner for the propensity score. | Low to rescope, high to extend | Reviewers will ask; adding a scenario after seeing results risks looking tuned | State the limit, or add a small-sample scenario now, before results exist. |
| 11 | Close the methods checklist gaps: covariate missing data, multiplicity, the assumption tested by the probable-match and walk-in analyses, the risk ratio, reasons for covariates, and a registry or cited source for typed design constants. | Low | Low | None. |
| 12 | Small fixes: braced `van der Laan` and the published Shaw entry in the vignette bibliography, the duplicated identification paragraph, `1,500` as inline R, British spellings, "Rescue.co", `emulation` in `allow_words`, and the missing `overlap_sources` path. | Low | Low | The path of the unblinded analysis's manuscript for the overlap check (I did not look for it). |

## 11. Writing check on this file

The check was run on this file from the project root in Git Bash with `--min-level warn`: `python "C:/Users/andre/OneDrive/Documents/scientific-writing-skills/checks/check.py" "audit/2026-10-22/05-writing.md" --min-level warn`.

It reports no pattern, number, placeholder or acronym findings after the fixes made while writing. One warning remains, the missing overlap source in `writing/config.yml` (it concerns the project configuration, not this file, and I did not look for the file). Every other finding is one error of a single kind, "looks like code formatting in prose". Each is a backticked code identifier in the sense the audit rules allow me to leave: a file path with a line number, a registry key, a bibliography key, an argument or function name of the package, or a coined term that I quote so that it can be replaced. They are expected in an internal audit that cites evidence as path and line. Nothing else remains. The proposed replacement text in section 9 is written to pass the check as prose, and its only backticked items are the package's own argument, function and column names, which stay as they are in the code.

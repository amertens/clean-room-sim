# Proposed changes, 20 October 2026 (not applied)

Proposed edits to `manuscript/methods.qmd`, `writing/brief.md` and the cleanTMLE "Get started" vignette, following the decisions of 20 October (spec `docs/superpowers/specs/2026-10-18-outcome-blind-design-stage.md`, section "Case study and paper on 0.5.0", and `2026-10-20-rescueco-reexecution-0.5.0.md`). Nothing has been edited. Placeholders: `[R: ...]` is a registry number still to come; `[cite: ...]` is a source still to find through sci-cite; `[TO FILL: ...]` is a fact only Andrew can give.

## A. manuscript/methods.qmd

**A1. The staged workflow: a new paragraph after line 39.**

After: The plan also declares time zero, the moment at which, in the target trial being emulated, treatment is assigned, eligibility is met and follow-up starts [@hernan2016], and the lock takes a variable dictionary that records, for each column of the design data, the form it was recorded on, when it was recorded and when the characteristic it describes was fixed relative to time zero, and how the analysis uses it. Because adjusting for a characteristic that treatment can change can remove part of the effect being estimated, the lock refuses a covariate fixed after time zero unless it is declared as a stand-in for a characteristic fixed before time zero, with its reason and a sensitivity analysis that leaves it out. It also refuses design columns recorded only for patients reached at follow-up, and it requires a written response to each eligibility criterion assessed after time zero and to follow-up that starts at another event. The SAP lists these declarations, each flag and its response in a table that follows the components of a target trial emulation [@cashinTransparentReportingObservational2025]. The table records what was declared: the simulation assumes identification, and the table, the negative controls and the bias analyses are what bear on whether it holds.

**A2. Line 41.**

Before: any later change to the plan, the design data or the decision changes a hash

After: any later change to the plan, the variable dictionary, the design data or the decision changes a hash

**A3. Lines 165 to 169 (who was compared).**

Before: Transport by Rescue.Co was confirmed against Rescue.Co's dispatch records, and patients recorded as transported by Rescue.Co with no dispatch record were excluded [R: number excluded]. We compared transport by Rescue.Co with other modes of transport to hospital in five contrasts specified before outcome access. The primary contrast, the protocol's comparison, was Rescue.Co transport versus all other transport. The others were Rescue.Co transport versus arrival without an ambulance; any ambulance versus arrival without an ambulance; and Rescue.Co transport versus other ambulances, with transfers included and, in the primary cohort, with transfers and patients who had received prior care excluded. The cohort sizes are given in the Results.

After: Transport by Rescue.Co was confirmed against Rescue.Co's dispatch records; patients recorded as transported by Rescue.Co with no dispatch record [R: number], and ambulance arrivals with no recorded operator [R: number], were excluded. We declared four comparisons in order before outcome access, each with its own plan and lock: the protocol's comparison of transport by Rescue.Co with all other transport to hospital; transport by Rescue.Co versus other ambulance operators; transport by Rescue.Co versus arrival without an ambulance; and any ambulance versus arrival without an ambulance. The protocol's comparison was run on all arrivals and on patients who arrived directly from the scene, and the other three on direct arrivals, that is, patients who had not been transferred from another facility or sought care elsewhere before arrival. The comparison that carried the confirmatory claim was chosen by a rule written into the registered plan that used no outcome: the first comparison, in the declared order, whose overlap was graded good or moderate, whose primary estimand kept at least 80% of each arm and shifted none of five named covariates by more than 0.25 standard deviations, and whose negative-control verdict was not STOP. The cohort sizes are given in the Results.

**A4. Line 170 (the estimand).**

Before: Because the registry enrols patients on arrival at hospital, after the transport being compared, patients who died before arrival are largely absent, and the estimand is the effect of transport mode among patients who reach a participating hospital alive; this equals the effect in all injured patients only if transport mode does not change survival to arrival.

After: The target population is injured patients for whom a transport choice was made, and time zero is the transport decision, measured as the time of injury because patients who arrived by other means have no dispatch time; the six-month horizon is counted from arrival, which follows the injury within hours. The registry enrols patients only on arrival alive at a participating hospital, after the transport being compared, so it identifies the effect in this population only if transport mode does not change survival to arrival. Before outcome access we therefore registered two analyses of this assumption: a quantitative bias analysis that adds hypothetical deaths before arrival to each arm in turn, drawn from the most severely injured patients, with its bias parameters and their plausible range fixed in advance [@foxApplyingQuantitativeBias2021; cite: source of the plausible range of deaths before arrival]; and a restriction to patients who would almost surely have arrived alive by any route.

**A5. Lines 171 and 172.**

Before: ...arrival physiology, which transport could affect, was added in a sensitivity analysis. The design stage was run on the complete-case cohort the analysis uses, because which patients lack follow-up is known before the outcome is seen and is not the outcome itself.

After: ...arrival state (heart rate, blood pressure, the components of the Glasgow Coma Scale and palpable pulse), which transport could affect, was declared a stand-in for severity at the scene and added in a sensitivity analysis. The design stage was run on the complete-case cohort the analysis uses, because which patients lack follow-up is known before the outcome is seen and is not the outcome itself; a data steward supplied two indicators for this, whether disposition and whether six-month status were recorded, and the SAP reports the share of each arm lost. To show what the design stage would have found at the protocol's own outcome-blind check, which used the first four months of enrolment, we also ran it on the patients enrolled in those months, before any six-month follow-up existed.

**A6. Line 177.**

Before: Payment method, recorded at discharge, and the distance from the injury to the receiving hospital, which depends on where transport took the patient, were not adjusted for.

After: Payment method, recorded at discharge, and the hospital reached and the distance to it, which depend on where transport took the patient, were not adjusted for.

**A7. Line 184.** Delete the placeholder `[TO FILL: the bias analysis for deaths before arrival, ...]`; A4 now describes it.

**A8. Line 185 (the plan).**

Before: For each contrast, the plan specified the estimands in the order ATE, trimmed ATE, ATT and ATO, the candidates, the tolerances and the negative-control outcomes with their cohort restrictions, and the outcome was joined only after the SAP was approved.

After: For each comparison, the plan specified the estimands in the order ATE, trimmed ATE, ATT and ATO, the candidates, the negative controls, and tolerances tied to the protocol's power calculation, which assumed a 22% risk of severe disability or death in the comparison arm and a relative risk of 0.75 with Rescue.Co: an absolute bias of at most half the resulting difference of 5.5 percentage points, and an expected CI half-width of at most that difference. For death within 180 days, for which the protocol gave no power calculation, the bias tolerance was 1 percentage point. The negative-control null band was set from the controls' prevalence. The outcome was joined only after the SAP was approved by [TO FILL: approver] and the lock hashes were registered [TO FILL: OSF registration].

**A9. Line 186 (disclosure), two sentences added.**

After (added): The checks of time zero and variable timing were added to the package after the unblinded analysis of this registry had found these problems; to guard against fitting them to this study, each runs in every study, has no threshold tuned to these data, and was tested on simulated data before this design stage was run. [TO FILL: whether the approver had seen results of the unblinded analysis.]

**A10. Lines 188 to 191 (missing outcomes), replaced.**

After: Outcomes missing at follow-up were handled by complete case analysis [R: number and share of missing outcomes by arm], the estimator the design stage evaluated. Because a patient who died in the emergency department always has a recorded outcome while only survivors can be lost, complete case drops lost survivors but keeps every early death, so its risk in an arm is biased upward unless the patients lost resemble those reached, and more so in the arm with more losses. The protocol (version 1, 9 May 2023; Supplementary Material) specified inverse probability of censoring weighting, which weights followed-up patients, including those who died in the emergency department, to stand in for the patients lost and so has the same problem by construction. We therefore report, as a sensitivity analysis, a longitudinal TMLE [@williams2020; @diaz2021] that treats deaths at emergency-department disposition as observed for every patient and assumes follow-up missing at random only among patients who left the emergency department alive; it was not evaluated in the design stage.

**A11. Line 192.**

Before: Patients recorded as transported by Rescue.Co whose identity Rescue.Co had not yet confirmed [R: number] were kept in the primary analysis and excluded in a sensitivity analysis.

After: Patients recorded as transported by Rescue.Co whom Rescue.co's dispatch records matched only probably [R: number] were kept in the primary analysis and excluded in a sensitivity analysis, and patients who walked in were excluded, in a further sensitivity analysis, from the comparisons with arrival without an ambulance.

**A12. Line 193.**

Before: Other deviations from the protocol, each with its reason, are listed in Supplementary Table S[X]: ...

After: Other deviations from the protocol, each with its reason and whether estimates from the unblinded analysis existed when it was made, are listed in Supplementary Table S[X]: ... (and "the exclusion of transfers and patients with prior care from the primary cohort" becomes "the exclusion of transfers and patients with prior care from every comparison except the protocol's comparison over all arrivals")

**A13. Line 197 (software).** "cleanTMLE package version 0.4.1" becomes "version 0.5.0", and lmtp is added with its version [R: lmtp version].

## B. writing/brief.md

**B1. Claim 5.** Before: "In the Rescue.Co Kenya trauma registry, the design stage found the protocol's ATE not estimable and moved to estimands the data could support, before any outcome was joined." After: "(proposed) In the Rescue.Co Kenya trauma registry, the outcome-blind design stage, run on the first four months of enrolment and on the final cohort, flagged before any outcome was joined the identification problems and design choices that the unblinded analysis changed while its estimates were available, and chose the confirmatory comparison, the estimand and the estimator by rules fixed in advance." (To be rewritten to what the re-execution shows, once it has run.)

**B2. What not to say.** Replace "The Rescue.Co case study's missing outcomes are handled by complete case analysis, a change from the protocol's IPCW; say so and report IPCW as a sensitivity analysis." with "The Rescue.Co case study's missing outcomes are handled by complete case, the estimator the design stage evaluated; say why it is weaker here than a longitudinal TMLE (early deaths are always observed, only survivors can be lost), and report the longitudinal TMLE as a sensitivity analysis not evaluated in the design stage." Add: "Not that the identification table shows identification holds: it lists what was declared."

## C. cleanTMLE "Get started" vignette

**C1. Step 6, after the dossier is printed:** a chunk printing `dossier$identification`, with two sentences: the table lists what the plan and the dictionary declare about time zero and the timing of each variable, with each flag and its response; here nothing is flagged, because time zero is treatment assignment and every covariate is fixed before it.

**C2. A short section before the references, "What this version does not do":** a lock compares two treatments, so a study with several comparisons uses one lock for each; missing outcomes are handled by complete case (`estimate_effect()` drops patients without an observed outcome while the simulation assumes every outcome is observed, so when outcomes are missing the lock should hold the patients whose outcome is observed); restrictions say where the negative controls are checked, not target populations for the estimand ladder. Planned additions are listed in WISHLIST.md.

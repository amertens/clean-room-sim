# Clean-room, staged-analysis, and outcome-blind literature: findings for the cleanTMLE revision

Compiled September 2026. Abstracts and open-access sources only; paywalled full texts (Wiley, Dove, JClinEpi) were summarized from abstracts, indexing services, and secondary sources. Each entry gives the citation, a short summary, and a "Borrow:" line saying what cleanTMLE could take from it. The synthesis section ranks the ten most actionable borrowings.

## Papers

### Topic 1. The staging and clean room anchor paper and its citation network

**Muntner P, Hernandez RK, Kent ST, Browning JE, Gilbertson DT, Hurwitz KE, Jick SS, Lai EC, Lash TL, Monda KL, Rothman KJ, Bradbury BD, Brookhart MA (2024). Staging and clean room: constructs designed to facilitate transparency and reduce bias in comparative analyses of real-world data. Pharmacoepidemiology and Drug Safety 33(3):e5770. PMID 38419140.**
The paper defines the two constructs cleanTMLE operationalizes. Staging is "performing sequential preliminary analyses and evaluating the population size available and potential bias before conducting comparative analyses." The clean room is an environment with "restricted access to data and preliminary results, policies governing exploratory analyses and protocol deviations, and audit trail." Their staged sequence runs from protocol development and registration under restricted data access, to attempting covariate balance through propensity score modeling, to evaluating residual confounding with negative control outcomes, and only then to final comparative analyses. At each checkpoint a masked review team recommends one of three actions, to proceed, to conduct additional analyses, or to terminate the study. The authorship spans Amgen's Center for Observational Research, academic epidemiology (Lash, Rothman, Jick, Brookhart), and NoviSci/Target RWE, so the construct already has both industry and methods credibility.
Borrow: the checkpoint vocabulary. Muntner's three-way recommendation (proceed, conduct additional analyses, terminate) is the published analog of the PASS/FLAG/SEVERE gate; adopting their terms in the paper and mapping cleanTMLE's grades onto them anchors the package in the literature rather than in invented jargon. Also borrow the three clean-room ingredients (restricted access, exploratory-analysis policy, audit trail) as an explicit checklist that cleanTMLE's config and logging implement.

**Levintow SN, Orroth KK, Breskin A, Park AS, Flores-Arredondo JH, Dluzniewski P, Navar AM, Sorensen HT, Brookhart MA (2022). Use of negative control outcomes to assess the comparability of patients initiating lipid-lowering therapies. Pharmacoepidemiology and Drug Safety 31(4). PMID 34894377.**
The methodological precursor to the clean-room paper from the same group. Before any effectiveness analysis of PCSK9 inhibitors versus ezetimibe or high-intensity statins, they estimated one-year risks of negative control outcomes chosen to detect two specific bias channels, frailty (pressure ulcers, accidents, fractures) and health-seeking behavior (screening visits). PCSK9 initiators had systematically lower frailty-related risks, telling the team, before outcomes were touched, that unadjusted or even IPW-adjusted comparisons would be confounded by frailty.
Borrow: structuring the negative-control panel by named bias mechanism (frailty channel, health-seeking channel) rather than as an undifferentiated list, and reporting NCO effect estimates per mechanism in the pre-outcome decision report.

**Zou Y, et al. (2026). Use of negative control outcomes to assess cohort comparability among evolocumab initiators: a real-world database study. Clinical Epidemiology, DOI 10.2147/CLEP.S608532. PMID 42609774. Cites Muntner 2024.**
Direct follow-up application in a Chinese healthcare database. Residual bias was evaluated by estimating hazard ratios for negative control outcomes under two candidate study designs with inverse probability weighting, and the NCO results were used to judge which design produced comparable cohorts before comparative effectiveness analyses.
Borrow: running the NCO panel under each candidate design (for cleanTMLE, under each candidate estimand and weighting scheme) and letting the comparison inform the pre-outcome design choice, not just a binary proceed decision.

**Horner ME, Ogdie A, Orroth KK, et al. (2025). Implementing negative control outcomes to assess comparability of treatments for psoriasis and psoriatic arthritis. Pharmacoepidemiology and Drug Safety 34(5):e70156. PMID 40387023.**
Another Amgen-network implementation. NCOs detected healthy-user bias in the psoriasis comparisons; adding an eligibility criterion (prior topical therapy) plus weighting reduced the residual confounding to acceptable levels in psoriatic arthritis analyses. The paper shows the iterate-then-recheck loop in practice, design change, rerun NCOs, reassess.
Borrow: the loop semantics for FLAG outcomes, a FLAG should be resolvable by a documented design amendment (new eligibility criterion, different weighting) followed by a rerun of the same diagnostics, with both rounds retained in the audit trail.

**Johnson M, et al. (2024). Real-world comparative effectiveness of sotorasib versus docetaxel in second line and beyond among patients with advanced NSCLC. Lung Cancer, DOI 10.1016/j.lungcan.2024.107960. PMID 39369609. Karim N, et al. (2026). Same comparison in England national data. Lung Cancer, DOI 10.1016/j.lungcan.2026.109469. PMID 42229337. Both cite Muntner 2024.**
Applied oncology studies citing the clean-room construct. Notably, the 2024 study balanced treatment groups via overlap weighting propensity score methods, that is, the design team committed to the ATO estimand at the design stage, exactly the estimand-selection-before-outcomes move cleanTMLE is formalizing. The 2026 England replication used propensity score weighting on national registry data.
Borrow: cite these as evidence that clean-room-styled applied studies already choose overlap weighting at design time; cleanTMLE's contribution is making that choice principled, graded, and documented rather than tacit.

**Pottegard A, et al. (2025). Ten must-read papers on transparency and reproducibility in pharmacoepidemiology. Pharmacoepidemiology and Drug Safety 34:e70119. PMID 40012278. Wang SV, et al. (2025). Advancing research transparency and reproducibility in pharmacoepidemiology. Pharmacoepidemiology and Drug Safety 34:e70096. PMID 39948333.**
The Muntner paper was selected among ten must-read papers on transparency and reproducibility, with an accompanying editorial in the same journal. This confirms the clean-room construct has been picked up by the field's transparency agenda rather than remaining a one-off proposal.
Borrow: framing. cleanTMLE's introduction can position the package inside the pharmacoepidemiology transparency-and-reproducibility movement, citing the editorial pair for the movement and Muntner for the construct.

**Hernandez RK, Critchlow CW, Dreyer N, Lash TL, Reynolds RF, Sorensen HT, Lai EC, Wang SV, Bradbury BD, Brookhart MA (2025). Advancing principled pharmacoepidemiologic research to support regulatory and healthcare decision making: the era of real-world evidence. Clinical Pharmacology and Therapeutics 117(5). PMID 39807817. PMC11924150.**
A review from the clean-room authors situating staging and clean rooms inside a stepwise framework for designing and transparently executing RWE studies, alongside negative control outcomes and exposures for probing conditional exchangeability, protocol registration, and RECORD-PE reporting.
Borrow: the phrase "stepwise framework for designing and transparently executing RWE studies" as literature-anchored language for what the package pipeline is.

### Topic 2. Blinded and outcome-blind analysis in epidemiology and RWE

**MacCoun R, Perlmutter S (2015). Blind analysis: hide results to seek the truth. Nature 526:187-189. PMID 26450040.**
The canonical cross-disciplinary argument for blind analysis, from particle physics (Perlmutter) and social science (MacCoun). Data are perturbed, masked, or partially hidden so that all analytic decisions are finalized before anyone can see how those decisions move the result. They catalog concrete blinding techniques, adding noise or hidden offsets to estimates, scrambling group labels, masking outcome cells, and argue confirmation bias operates even on honest analysts, which is why procedural blinding beats good intentions.
Borrow: two things. First, the motivating argument (decisions made while results are visible are contaminated even without bad faith) for the manuscript introduction. Second, the technique of label or offset masking as a lightweight alternative when a physical two-team separation is impossible, which is the single-analyst fallback mode cleanTMLE could support.

**Dutilh G, Sarafoglou A, Wagenmakers EJ (2021). Flexible yet fair: blinding analyses in experimental psychology. Synthese 198:5745-5772.**
Develops the practicalities of blinded analysis for fields without physics-style pipelines, including the two-person protocol where a data manager creates a blinded dataset (shuffled or masked key variables) on which the analyst develops and freezes the full analysis script, which is then run once on the real data.
Borrow: the two-person operational protocol maps directly onto the design-team/outcome-team split; useful citation for the claim that the workflow generalizes beyond pharma settings to any two-person research group.

**Franklin JM, Pawar A, Martin D, Glynn RJ, Levenson M, Temple R, Schneeweiss S (2020). Nonrandomized real-world evidence to support regulatory decision making: process for a randomized trial replication project. Clinical Pharmacology and Therapeutics 107(4):817-826. PMID 31541454.**
The RCT-DUPLICATE process paper. Each trial emulation was designed and registered (clinicaltrials.gov plus a public protocol) with all design and analysis decisions locked before outcome results were generated, and prespecified agreement metrics against the index RCT. The process explicitly separates a design phase using feasibility counts and balance checks from the outcome analysis.
Borrow: the practice of registering each analysis with its decision rules and agreement criteria before results exist, and the idea of prespecified success metrics for the pipeline itself, not only the study.

**Franklin JM, Patorno E, Desai RJ, et al. (2021). Emulating randomized clinical trials with nonrandomized real-world evidence studies: first results from the RCT DUPLICATE initiative. Circulation 143(10):1002-1013. PMID 33327727. Wang SV, Schneeweiss S, et al. (2023). Emulation of randomized clinical trials with nonrandomized database analyses: results of 32 clinical trials. JAMA 329(16):1376-1385. PMID 37097356.**
The results papers. Across 32 emulations with new-user cohorts and propensity score matching in three claims databases, prespecified agreement criteria (regulatory agreement, estimate agreement, standardized difference agreement) were met in most but not all pairs, and closeness of emulation predicted agreement.
Borrow: the three-tier agreement metrics are a template for cleanTMLE's plasmode validation reporting, and the finding that design fidelity drives agreement is quantitative support for spending effort at the design stage.

**Wang SV, Pinheiro S, Hua W, Arlett P, Uyama Y, Berlin JA, Bartels DB, Kahler KH, Bessette LG, Schneeweiss S (2021). STaRT-RWE: structured template for planning and reporting on the implementation of real world evidence studies. BMJ 372:m4856. PMID 33436424. PMC8489282.**
A structured tabular template (design diagram, PICOT parameter tables, operational definitions with code lists, day-zero anchoring) intended to replace ambiguous prose in RWE protocols and reports, so that a study can be reproduced from the tables alone.
Borrow: the design diagram and parameter-table format for the design packet that crosses cleanTMLE's firewall. The packet the outcome team receives should be a STaRT-RWE-style table set, not prose.

**Wang SV, Pottegard A, Crown W, et al. (2023). HARmonized Protocol Template to Enhance Reproducibility of hypothesis evaluating real-world evidence studies on treatment effects: a good practices report of a joint ISPE/ISPOR task force. Pharmacoepidemiology and Drug Safety 32(1):44-55 (co-published in Value in Health). PMID 36215113.**
HARPER merges existing templates into a single protocol skeleton with a shared level of operational detail, covering rationale, estimand-relevant parameters, data provenance, design, analysis, and fitness-for-purpose assessment, meant to be the backbone from protocol through registration to reporting.
Borrow: name cleanTMLE's exported design report sections after HARPER sections, so the package output drops into a HARPER protocol with minimal editing. This is a low-cost, high-visibility interoperability feature.

**FDA (2023). Considerations for the use of real-world data and real-world evidence to support regulatory decision-making for drug and biological products. Final guidance, August 2023. FDA (2024). Real-world evidence: considerations regarding non-interventional studies for drug and biological products. Draft guidance, March 2024.**
The 2023 final guidance requires sponsors to ensure transparency about data access and analyses conducted before protocol finalization, and to document when the data were accessed relative to when the protocol and SAP were finalized, precisely to prevent design shopping on visible outcomes. The 2024 draft extends this to non-interventional studies, recommending a finalized protocol and detailed SAP (design schema, causal diagram, defined source population, eligibility, variables), documented data-source rationale (reliability, completeness, timing of key elements), and handling of confounding, missingness, and misclassification, before the study analyses.
Borrow: the audit-relevant facts cleanTMLE should stamp into its logs are exactly what FDA asks sponsors to attest: the date of first data access, the date of protocol/SAP freeze, which analyses ran before the freeze, and any deviations after it. A machine-generated attestation table would be a distinctive package feature.

**ICH (2025). M14: general principles on planning, designing, analysing, and reporting of non-interventional studies that utilise real-world data for safety assessment of medicines. Step 4 adopted September 2025; adopted as final FDA guidance March 2026.**
The first global (ICH) guideline for non-interventional RWD studies. It formalizes the feasibility assessment as a distinct early phase whose goal is to describe and compare the relevance and reliability of candidate data sources for the research question before protocol finalization, promotes fit-for-purpose data determinations, and takes a "show your work" stance on documentation from planning through reporting.
Borrow: vocabulary and legitimacy. "Feasibility assessment," "fit-for-purpose," and "show your work" are now ICH-anchored terms; cleanTMLE's stage 1-4 outputs can be described as the feasibility assessment deliverables in M14's sense.

**ENCePP (ongoing). Guide on methodological standards in pharmacoepidemiology (revision 12) and the ENCePP Checklist for Study Protocols; HMA-EMA Catalogues of RWD sources and studies (successor to the EU PAS Register, February 2024).**
The EU infrastructure for prespecification, protocols registered publicly before study start, a methodological checklist, and a code of conduct emphasizing scientific independence.
Borrow: the practice of registering the frozen design (cleanTMLE's pre-outcome decision report) in a public registry, with the HMA-EMA catalogue or OSF as the venue, and citing ENCePP for why registration precedes outcome analysis.

**Clinical-trial dry-run practice (statistical analysis plans; see also Leveraging Synthetic Clinical Data for Validation and Operational Readiness in Clinical Trials, Healthcare (MDPI) 2026;14(17):2870).**
In regulated trials, the analysis pipeline is routinely exercised before database lock using a dummy randomization list or synthetic data, producing all tables, figures, and listings under blinding so programming and algorithms are frozen before unblinding. This is standard operating procedure rather than a research method, which is exactly its value as a citation, outcome-blind rehearsal of the full pipeline is normal science in trials.
Borrow: the term "dry run" for cleanTMLE's stage that executes the entire TMLE workflow on simulated or permuted outcomes, and the requirement that dry-run outputs be the same artifact set as the final analysis (identical tables and figures, dummy numbers).

### Topic 3. Two-team and firewall designs

**Rubin DB (2007). The design versus the analysis of observational studies for causal effects: parallels with the design of randomized trials. Statistics in Medicine 26(1):20-36. Rubin DB (2008). For objective causal inference, design trumps analysis. Annals of Applied Statistics 2(3):808-840.**
The founding argument for outcome-free design. Observational studies should be designed to approximate randomized experiments "without examining any final outcome data," with all balancing, matching, and eligibility work completed and frozen first. Critically for cleanTMLE, Rubin is explicit that the design phase can end in rejection, a candidate dataset may have to be "rejected as inadequate" because key covariates are missing or because of "lack of overlap in the distributions of key covariates between treatment and control groups, often revealed by careful propensity score analyses."
Borrow: the rejection outcome. Rubin 2008 is the citation that a no-go decision (SEVERE) is a legitimate, indeed intended, product of outcome-free design, not a failure of the study team. Also the phrase "outcome-free design," which the device literature adopted.

**Yue LQ (2012). Regulatory considerations in the design of comparative observational studies using propensity scores. Journal of Biopharmaceutical Statistics 22(6). Yue LQ, Lu N, Xu Y (2014). Designing premarket observational comparative studies using existing data as controls: challenges and opportunities. Journal of Biopharmaceutical Statistics 24(5). PMID 25013971. Li H, Mukhi V, Lu N, Xu YL, Yue LQ (2016). A note on good practice of objective propensity score design for premarket nonrandomized medical device studies with an example. Statistics in Biopharmaceutical Research 8(3):282-286. Li H, Yue LQ (2023). Propensity score-based methods for causal inference and external data leveraging in regulatory settings: from basic ideas to implementation. Pharmaceutical Statistics 22(4).**
The FDA CDRH line of work that turned Rubin's principle into a regulatory two-stage design, used in actual device submissions since around 2013. Stage 1 (design) is conducted outcome-free by a statistician firewalled from outcome data: propensity score model building, balance assessment, definition of the analysis population, and sample size or power assessment are all completed and documented in a design report that is locked before any outcome data are linked. Stage 2 (analysis) executes the locked plan. The 2016 note gives good-practice details and a worked example; the 2023 paper provides step-by-step procedure templates for real proposals and extends the machinery to external-data leveraging (propensity-score-integrated composite likelihood).
Borrow: three mechanics. First, the "design report" as a formal, signed, frozen artifact, cleanTMLE's decision report should explicitly be one. Second, the independent-statistician firewall role, named in the config (who holds outcomes, who signs the freeze). Third, doing power and precision assessment inside stage 1 using only design-side quantities, which supports a minimum-detectable-effect gate before outcomes.

**Major-Pedersen A, et al. (2021). A joint industry-sponsored data monitoring committee model for observational, retrospective drug safety studies in the real-world setting. Pharmacoepidemiology and Drug Safety 30(1):9-16. PMC8247341.**
Documents governance for a four-sponsor observational safety program (GLP-1 receptor agonists and medullary thyroid cancer). A DMC saw unblinded drug-level data in closed sessions; sponsors, FDA, and the steering committee received only blinded class-level results through a CRO firewall; the DMC issued standardized graded recommendations (continue unaltered, continue with modifications, report to FDA early).
Borrow: the standardized recommendation grammar and the closed-session/open-session reporting split, a concrete precedent for cleanTMLE's masked review team producing a graded recommendation while most stakeholders see only the blinded summary.

**Platt R, Brown JS, Robb M, et al. (2018). The FDA Sentinel Initiative: an evolving national resource. New England Journal of Medicine 379:2091-2093. Desai RJ, Weberpals J, et al. (2025). Strengthening inferential studies in the FDA Sentinel initiative: results from a methodological demonstration project. npj Digital Medicine, DOI 10.1038/s41746-025-02234-5. PMID 41398062.**
Sentinel's routine querying tools are reusable, pre-tested, parameterized programs run against a common data model, so an analysis is specified by filling in parameters rather than writing new code, which is prespecification enforced at the software level. The 2025 demonstration project shows Sentinel's newer linked EHR-claims enterprise executing a target trial emulation with a validated computable phenotype and large-scale propensity score fine stratification.
Borrow: the reusable-module philosophy, cleanTMLE stages should be parameterized, versioned functions whose parameter files fully determine the analysis, so the audit trail is a diff of config files rather than a code review.

### Topic 4. Outcome-blind simulation and estimator selection

**Franklin JM, Schneeweiss S, Polinski JM, Rassen JA (2014). Plasmode simulation for the evaluation of pharmacoepidemiologic methods in complex healthcare databases. Computational Statistics and Data Analysis 72:219-226. PMID 24587587.**
The founding plasmode paper. Resample real cohort covariates and exposure (preserving their joint distribution and confounding structure), then simulate outcomes from a model fitted to the real data with an injected known treatment effect. Methods are then evaluated against the known truth in data that retain realistic complexity.
Borrow: the core engine, which cleanTMLE already uses; the citation matters mainly for provenance and for the design choice of preserving the covariate-exposure structure exactly while only the outcome is synthetic, which is what makes the simulation outcome-blind compatible.

**Schreck N, Slynko A, Saadati M, Benner A (2024). Statistical plasmode simulations: potentials, challenges and recommendations. Statistics in Medicine (arXiv 2305.06028).**
A general methods-and-reporting piece on plasmode simulation, giving a stepwise procedure for generation, implementation, and reporting, with attention to which parts of the data-generating process are known versus inherited from the real data.
Borrow: the reporting checklist for the plasmode stage, cleanTMLE's simulation report should state what was preserved, what was synthesized, what truth was injected, and how many replicates were run, in a fixed format.

**Shaw PA, Gruber S, Williamson BD, Desai R, Shortreed SM, Krakauer C, Nelson JC, van der Laan MJ (2025). A cautionary note for plasmode simulation studies in the setting of causal inference. arXiv 2504.11740.**
Compares two popular plasmode frameworks and shows one of them can make sound estimators appear biased with below-nominal coverage, artifacts that persist at large sample sizes and even with trial data. The mechanism is that the outcome-generation shortcut distorts the estimand and the positivity structure of the plasmode world relative to the analysis model, so estimator rankings from the simulation are misleading.
Borrow: this is the central methodological caution for cleanTMLE's outcome-blind simulation stage. The package should implement (and the paper should state) the framework Shaw et al. recommend, and should verify that the plasmode world reproduces the observed propensity structure, since our whole use case is severe practical positivity violation, exactly where the wrong plasmode recipe lies most.

**Nance N, Petersen ML, van der Laan M, Balzer LB (2024). The causal roadmap and simulations to improve the rigor and reproducibility of real-data applications. Epidemiology 35(6):791-800. PMID 39087681. PMC11444352.**
The most operational outcome-blind simulation paper. Simulations conducted after data collection but before effect estimation preserve real covariate relationships while simulating outcomes, keeping analysts blinded to the true exposure-outcome association. Candidate estimators (Super Learner library variants, with and without covariate prescreening, with and without propensity truncation; influence-curve versus bootstrap variance) are compared over about 1000 replicates with a two-step selection rule, first choose the nuisance-estimation approach minimizing empirical variance subject to adequate oracle coverage (intervals built from the true sampling variance), then choose the variance estimator minimizing estimated variance subject to nominal 95 percent coverage. The winning configuration is written into the SAP, which is finalized before outcomes are used. Case studies include an SGLT2i observational analysis with positivity violations.
Borrow: the two-step selection criterion verbatim, including oracle coverage as the intermediate metric, and the framing that the deliverable of the simulation stage is a finalized SAP naming one estimator configuration. This is the strongest available template for cleanTMLE's estimator-selection stage.

**Dang LE, Gruber S, Lee H, et al. (2023). A causal roadmap for generating high-quality real-world evidence. Journal of Clinical and Translational Science 7(1):e212.**
Report from FIORD (Forum on the Integration of Observational and Randomized Data, Washington DC, November 2022), co-sponsored by Berkeley's Forum for Collaborative Research, the Center for Targeted Machine Learning, and the Joint Initiative for Causal Inference. Lays out the seven-step causal roadmap (causal question and estimand, observed data, identifiability, statistical estimand, estimator choice, sensitivity analysis, design comparison) and endorses outcome-blind simulations that compare estimators "without information on the observed treatment-outcome association" against prespecified benchmarks (type I error, coverage, bias, precision), with all roadmap steps prespecified.
Borrow: the term "outcome-blind simulation" with this citation as its definitional source, the prespecified benchmark list, and roadmap step 3 (assess identifiability before choosing the statistical estimand), which is precisely cleanTMLE's estimability gate placed inside an established framework.

**Petersen ML, Porter KE, Gruber S, Wang Y, van der Laan MJ (2012). Diagnosing and responding to violations in the positivity assumption. Statistical Methods in Medical Research 21(1):31-54. PMID 21030422.**
Presents the parametric bootstrap as a diagnostic for positivity-driven sparsity, simulate from an estimated data-generating distribution in which the true effect is known, re-run the estimator, and read off its finite-sample bias, which quantifies how badly sparsity plus model extrapolation hurt the actual analysis. Then catalogs the response options, restrict the adjustment set, change the projection function (working MSM), restrict the sample (trimming), or modify the target intervention, and frames all of them as "trading off proximity to the initial target of inference for identifiability," advocating that the tradeoff be approached systematically.
Borrow: two things. The parametric bootstrap bias estimate is a natural severity meter behind the PASS/FLAG/SEVERE gate (a quantitative complement to PoRT's subgroup detection). And the sentence-level framing "trading off proximity to the initial target for identifiability, approached systematically" is the estimand ladder's intellectual anchor, fourteen years before we needed it.

**Phillips RV, van der Laan MJ, Lee H, Gruber S (2023). Practical considerations for specifying a super learner. International Journal of Epidemiology 52(4):1276-1285.**
Guidance for prespecifying a Super Learner (library composition, cross-validation scheme, screening) so that the machine learning component of an analysis is itself fixed in the SAP rather than tuned on visible results.
Borrow: cleanTMLE's SAP export should serialize the full SL specification chosen during the outcome-blind stage, this paper is the citation for why and how.

**Xu Y, Gruber S, van der Laan MJ (2026). Investigating targeting strategies and truncation in TMLE for the average treatment effect under practical positivity violations. arXiv 2604.20059.**
Simulation study of TMLE variants under exactly cleanTMLE's problem setting. Loss-weighted targeting can induce substantial systematic bias relative to clever-covariate-scaled targeting; fixed truncation at c/(sqrt(n) log n) with c of 5 or 6 is a robust practical default; a Lepski-type adaptive truncation with a brake mechanism stabilizes data-adaptive tuning; and targeted bootstrap variance estimation is reliable across truncation levels where influence-curve variance is not.
Borrow: concrete defaults for the outcome stage, clever-covariate targeting rather than loss-weighted, default truncation c/(sqrt(n) log n) with c=5, and targeted bootstrap variance whenever the gate has flagged positivity, plus the caveat table for the documentation.

### Topic 5. Positivity and estimability assessment methods

**Danelian G, Foucher Y, Leger M, Le Borgne F, Chatton A (2023). Identification of in-sample positivity violations using regression trees: the PoRT algorithm. Journal of Causal Inference 11(1):20220032. Implemented as port() in the RISCA R package.**
PoRT runs a sequence of shallow regression trees on the covariates (not on the fitted propensity score) to find interpretable subgroups with extreme exposure prevalence, returning for each violating subgroup a covariate-defined rule plus its exposure prevalence, size, and share of the sample. Hyperparameters control minimal subgroup size (alpha), the prevalence threshold (beta), and tree depth (gamma). Because the output is clinical rules rather than a score histogram, it supports eligibility redefinition, the most defensible response to structural violations. A 2026 extension (dePoRT, arXiv 2607.25072) covers mediation settings.
Borrow: run PoRT in the estimability stage and print its subgroup strings directly in the decision report, giving the graded gate a clinically nameable basis ("males over 75 with GCS under 6 are never treated") instead of only percentile diagnostics. RISCA::port is CRAN-available, so integration is cheap.

**Oberst M, Johansson FD, Wei D, Gao T, Brat G, Sontag D, Varshney KR (2020). Characterization of overlap in observational studies. AISTATS 2020, PMLR 108. arXiv 1907.04138. Code: github.com/clinicalml/overlap-code.**
OverRule formalizes overlap estimation as finding minimum-volume sets under coverage constraints, solved with Boolean rule classifiers. It learns two rule sets, a support rule (where the population actually lives) and propensity overlap rules, and labels a patient as in the overlap set if the support rule applies and no exclusion rule fires. Output is again interpretable rules.
Borrow: the two-layer distinction between population support and treatment overlap. cleanTMLE's report can separate "regions with no data at all" from "regions with data but one-sided treatment," which matters for interpreting a trimmed estimand. OverRule is the machine learning-flavored alternative to PoRT worth citing even if only PoRT is implemented.

**Crump RK, Hotz VJ, Imbens GW, Mitnik OA (2009). Dealing with limited overlap in estimation of average treatment effects. Biometrika 96(1):187-199. Precursor working paper: Moving the goalposts: addressing limited overlap in the estimation of average treatment effects by changing the estimand (2006, NBER technical WP 330 / IZA DP 2347).**
Derives the variance-optimal trimmed subpopulation, discard units with extreme propensity scores, and shows a rule of thumb, keep units with propensity in [0.1, 0.9], approximates the optimum across many designs. The working-paper title says the quiet part loudly, trimming is changing the estimand, and the authors present that as the honest response to limited overlap.
Borrow: the 0.1/0.9 default for the trimmed-ATE rung of the ladder, the optimal-trimming variant as an option, and the "moving the goalposts" phrase (with citation) to make explicit in the documentation that each rung of the ladder is an estimand change, not a technical fix.

**Sturmer T, Rothman KJ, Avorn J, Glynn RJ (2010). Treatment effects in the presence of unmeasured confounding: dealing with observations in the tails of the propensity score distribution, a simulation study. American Journal of Epidemiology 172(7):843-854. Sturmer T, Webster-Clark M, et al. (2021). Propensity score weighting and trimming strategies for reducing variance and bias of treatment effect estimates: a simulation study. American Journal of Epidemiology 190(8):1659-1670.**
Asymmetric trimming (drop treated below a low quantile of the treated propensity distribution and untreated above a high quantile of the untreated distribution) targets unmeasured confounding concentrated in the tails, such as frailty-driven treatment withholding. The 2021 simulations found untrimmed estimates biased under unmeasured confounding regardless of weighting method, with Sturmer and Walker trimming the only strategies that consistently reduced bias.
Borrow: offer Sturmer asymmetric trimming alongside Crump symmetric trimming, and document the distinction, Crump trims for variance and identifiability, Sturmer trims for tail-concentrated unmeasured confounding. In a trauma registry where non-treatment of the sickest is the positivity problem, the Sturmer rationale is the clinically apt one.

**Walker AM, Patrick AR, Lauer MS, Hornbrook MC, Marin MG, Platt R, Roger VL, Stang P, Schneeweiss S (2013). A tool for assessing the feasibility of comparative effectiveness research. Comparative Effectiveness Research 3:11-20.**
Introduces empirical equipoise, comparative effectiveness is feasible where alternative therapies are used as-if interchangeably. The tool converts the propensity score to a preference score (adjusting for market share) and declares empirical equipoise when at least half of patients in each treatment group have preference scores between 0.3 and 0.7. This is a published, widely adopted (OHDSI) go/no-go feasibility gate computed entirely pre-outcome.
Borrow: implement the preference-score transformation and the 0.3-0.7 with 50 percent coverage rule as the first-line PASS/FLAG threshold at the propensity stage, citing Walker for the thresholds so the gate is anchored rather than invented.

**Li F, Morgan KL, Zaslavsky AM (2018). Balancing covariates via propensity score weighting. JASA 113(521):390-400. Li F, Thomas LE, Li F (2019). Addressing extreme propensity scores via the overlap weights. American Journal of Epidemiology 188(1):250-257. Thomas LE, Li F, Pencina MJ (2020). Overlap weighting: a propensity score method that mimics attributes of a randomized clinical trial. JAMA 323(23):2417-2418. Zhou T, Tong G, Li F, Thomas LE, Li F (2022). PSweight: an R package for propensity score weighting analysis. R Journal 14(1):282-300 (arXiv 2010.08893).**
The overlap-weights program. Within the balancing-weights family (each choice of tilting function defines a target population), overlap weights h(x) = e(x)(1-e(x)) target the population with clinical equipoise, are bounded, achieve exact mean balance on covariates in the propensity model when the score is fit by logistic regression, and minimize asymptotic variance among balancing weights. The AJE paper positions ATO as the principled response to extreme scores; the JAMA piece is the two-page clinical explainer; PSweight implements ATE/ATT/ATO/matching/entropy weights with balance and effective sample size reporting.
Borrow: the ATO rung's estimation machinery and its justification language ("population in whom clinical equipoise holds"), the exact-balance property as a reportable diagnostic, and PSweight's convention of reporting effective sample size per estimand, which makes the cost of each ladder rung visible in one table.

**Mao H, Li L, Greene T (2019). Propensity score weighting analysis and treatment effect discovery. Statistical Methods in Medical Research 28(8):2439-2454.**
Studies modified IPW estimators for overlap-type populations with analytic variance formulas that account for propensity estimation, and augments them with outcome models for a double-robustness-like efficiency property. Demonstrates substantial power gains under poor overlap relative to conventional IPW.
Borrow: the augmented (outcome-model-assisted) overlap estimator as the efficient analysis-stage counterpart for the ATO rung, cleanTMLE's outcome team can use an augmented estimator matched to the estimand the design team selected.

**Matsouaka RA, Zhou Y (2024). A framework for causal inference in the presence of extreme inverse probability weights: the role of overlap weights. Biometrical Journal 66(2) (arXiv 2011.01388). Zhou Y, Matsouaka RA, Thomas L (2020). Propensity score weighting under limited overlap and model misspecification. Statistical Methods in Medical Research 29(12). Liu Y, Li H, Zhou Y, Matsouaka RA (2024). Average treatment effect on the treated under lack of positivity. Statistical Methods in Medical Research 33(10) (DOI 10.1177/09622802241269646).**
The theory line for weighted estimands under positivity failure. The 2024 framework paper studies the family of weighted ATEs indexed by tilting functions (overlap, matching, entropy, beta weights), gives large-sample theory including when positivity fails, and characterizes overlap weights as the variance-minimizing member; the estimands "completely avoid positivity." The 2020 paper adds robustness under model misspecification; the ATT paper builds trimmed and weighted ATT analogs when the treated population itself lacks support.
Borrow: the tilting-function formalism as the unifying notation for cleanTMLE's ladder (ATE, trimmed ATE, ATT, ATO are all choices of h), and the weighted-ATT variants for cases where the trauma question is intrinsically about the treated.

**Liu Y, Wang Y, Gao Y, Poteat T, Matsouaka RA (2025). A tutorial for propensity score weighting methods under violations of the positivity assumption. arXiv 2511.10077. R package ChiPS.**
A tutorial synthesizing the above: positivity violations render ATE/ATT/ATC unidentifiable, weighted estimands (WATE/WATT/WATC) restore identifiability, with guidance on choosing the primary target estimand, implementing the estimators, and post-weighting diagnostics, illustrated on NHANES and HIV case studies.
Borrow: the post-weighting diagnostic set for the outcome-stage report, and the tutorial's estimand-selection guidance as the closest published analog to cleanTMLE's decision layer, worth engaging with explicitly in the manuscript since it is contemporaneous work solving the estimation half but not the governance half.

**Kennedy EH (2019). Nonparametric causal effects based on incremental propensity score interventions. JASA 114(526):645-656 (arXiv 1704.00211). Review: Bonvini M, McClean A, Kennedy EH (2021+). Incremental causal effects: an introduction and review. arXiv 2110.10532.**
Replaces deterministic interventions with stochastic ones that multiply each subject's odds of treatment by delta. The resulting effect curve in delta is identified with no positivity conditions at all, because no one is pushed to a treatment probability of 0 or 1, and has clean nonparametric estimators even longitudinally.
Borrow: the incremental effect as the ladder's final rung or companion sensitivity analysis when even ATO is strained, reporting the effect curve over a delta range answers "would treating more of these patients help" without pretending everyone could be treated.

**Cole SR, Hernan MA (2008). Constructing inverse probability weights for marginal structural models. American Journal of Epidemiology 168(6):656-664. PMID 18682488.**
The standard weight-diagnostics reference. Stabilized weights should have mean near one (a drifting mean flags misspecification or positivity problems), extreme weights should be examined, and progressive truncation traces out the bias-variance tradeoff, with the advice that diagnostics and sensitivity analyses are essential rather than optional.
Borrow: the mean-of-stabilized-weights check and a truncation-path plot (estimate and CI as a function of truncation percentile) as cheap standard panels in the propensity-stage report.

**Gruber S, Phillips RV, Lee H, van der Laan MJ (2022). Data-adaptive selection of the propensity score truncation level for inverse-probability-weighted and targeted maximum likelihood estimators of marginal point treatment effects. American Journal of Epidemiology 191(9):1640-1651. Also: Ju C, Schwab J, van der Laan MJ (2019). On adaptive propensity score truncation in causal inference. Statistical Methods in Medical Research 28(6) (arXiv 1707.05861).**
Both make truncation itself a prespecifiable, data-adaptive procedure (collaborative-TMLE flavored selection of the truncation bound that optimizes estimated MSE of the target) rather than an arbitrary analyst choice.
Borrow: expose a data-adaptive truncation option whose selection rule is fixed in the SAP, so the gate's "trim/truncate" response is algorithmic and auditable rather than discretionary.

**Zivich PN, Cole SR, Westreich D (2022). Positivity: identifiability and estimability. arXiv 2207.05010.**
Separates deterministic from stochastic positivity and, in parallel, identifiability from estimability, a parameter can be identified yet practically inestimable in the sample at hand, and data-adaptive (machine learning) estimators cannot manufacture support where none exists, they only interpolate or extrapolate less visibly.
Borrow: the vocabulary. cleanTMLE's central question is estimability, and this paper supplies the precise distinction (identification is about the population and assumptions, estimability is about the sample) that the decision report should use, plus the argument for why flexible nuisance models do not rescue a SEVERE gate.

**D'Amour A, Ding P, Feller A, Lei L, Sekhon J (2021). Overlap in observational studies with high-dimensional covariates. Journal of Econometrics 221(2):644-654 (arXiv 1711.02582).**
Strict overlap implies explicit bounds on covariate-mean imbalance that tighten as dimension grows, so demanding overlap on many covariates is increasingly restrictive, and estimated propensity overlap in high dimensions can be an artifact of regularization.
Borrow: justification for keeping the design covariate set parsimonious and clinically curated in the estimability assessment, and a caution against reading smooth estimated-propensity histograms as evidence of overlap when the covariate space is large.

**Ben-Michael E, Feller A, Hirshberg DA, Zubizarreta JR (2021). The balancing act in causal inference. arXiv 2110.14831. Ben-Michael E, Keele L (2022+). Using balancing weights to target the treatment effect on the treated when overlap is poor. arXiv 2210.01763 (published in Statistics in Medicine).**
The balancing-weights synthesis: finding weights that achieve approximate covariate balance is shrinkage estimation of inverse propensity weights, unifying the modeling and balancing traditions. The ATT paper applies this under poor overlap and, together with related applied guidance, supports a two-track workflow, estimate the ATT with balancing weights and the ATO with overlap weights in parallel; agreement is evidence the ATT is estimable, divergence is evidence only the ATO is.
Borrow: the parallel ATT/ATO agreement check as an operational rule inside the graded decision, it converts "which estimand can the data support" from judgment into a reportable comparison, executable outcome-blind on the plasmode side and prespecifiable for the outcome side.

**Nethery RC, Mealli F, Dominici F (2019). Estimating population average causal effects in the presence of non-overlap: the effect of natural gas compressor station exposure on cancer mortality. Annals of Applied Statistics 13(2):1242-1267.**
Defines the overlap region data-adaptively, fits BART inside the overlap region, and extrapolates effects into the non-overlap region with a restricted cubic spline carrying inflated uncertainty, so a population-average estimate is retained but its non-overlap component is visibly model-dependent.
Borrow: the reporting idea more than the estimator, when the gate keeps the ATE, decompose the estimate into an overlap-region component (data-supported) and a non-overlap component (extrapolated, wider intervals), so readers see exactly what part of the answer rests on modeling.

**Zhu Y, Hubbard RA, Chubak J, Roy J, Mitra N (2021). Core concepts in pharmacoepidemiology: violations of the positivity assumption in the causal analysis of observational data: consequences and statistical approaches. Pharmacoepidemiology and Drug Safety 30(11):1471-1485. PMC8492528.**
The field-facing review. Distinguishes structural violations (treatment impossible, unfixable by sample size) from practical violations (possible but unobserved, a finite-sample problem), catalogs consequences (unstable weights, silent extrapolation) and remedies (Crump trimming, Sturmer asymmetric trimming, cardinality matching, overlap weights, BART-plus-spline extrapolation), and closes with a decision-guidance figure, determine violation type, define overlap, match method to violation type, and communicate the resulting target population.
Borrow: the structural-versus-practical taxonomy as a required field in cleanTMLE's gate report (PoRT subgroups plus clinical review decide which type each violation is), and their Figure 4 flow as the published skeleton the package's decision logic instantiates. Also the ideal citation to open the manuscript's positivity section.

**Leger M, Chatton A, Le Borgne F, Pirracchio R, Lenain R, Foucher Y (2022). Causal inference in case of near-violation of positivity: comparison of methods. Biometrical Journal 64(8):1389-1403.**
Simulation comparison of g-computation, IPW, truncated IPW, TMLE, and truncated TMLE under near-violations (chance sparsity from low exposure prevalence or small samples). Near-violation degraded all methods; g-computation and TMLE-based methods were most robust; truncation limited bias under violation but introduced bias when positivity actually held.
Borrow: evidence for the package's defaults, TMLE (or g-computation) as the base estimator under FLAG conditions, with truncation applied conditionally on diagnosed violation rather than universally, since it costs bias when support is fine.

**Also noted for topic 5.** A diagnostic for positivity with continuous treatments (arXiv 2502.11820) and data-adaptive strategies for positivity in continuous interventions (arXiv 2502.14566) extend the toolbox beyond binary exposure; sample size and power calculation methods for observational causal inference (arXiv 2501.11181) support a design-stage minimum-detectable-effect computation.

### Topic 6. Graded decision rules, estimand ladders, and estimands language

**Conover MM, Ryan PB, Chen Y, Suchard MA, Hripcsak G, Schuemie MJ (2025). Objective study validity diagnostics: a framework requiring pre-specified, empirical verification to increase trust in the reliability of real-world evidence. JAMIA 32(3):518-525. DOI 10.1093/jamia/ocae317.**
The closest published system to cleanTMLE's gate, from OHDSI. Five prespecified, objective diagnostics are computed while effect estimates remain blinded: minimum detectable relative risk (threshold MDRR under 10), empirical equipoise (preference score in 0.3 to 0.7 for more than half of patients), covariate balance (maximum standardized difference at most 0.10), generalizability (standardized difference between analytic and target population under 0.25), and expected absolute systematic error from negative controls (EASE under 0.25). Results are unblinded only for analyses passing all thresholds. Applied to 11,716 negative-control analyses from LEGEND-HTN, gating reduced EASE from 0.38 to essentially zero and the share of null-excluding intervals from 15.2 to 3.9 percent, with equipoise the single most effective gate. The authors reframe the cost of gating as converting false positives into "inestimable" findings, an explicit third category between significant and null.
Borrow: nearly everything. The five diagnostics with their exact thresholds are an off-the-shelf quantitative core for PASS/FLAG/SEVERE; "unblind only if all prespecified diagnostics pass" is the operative sentence for the workflow; and the word "inestimable" as a formal result category is precisely the vocabulary cleanTMLE needs for comparisons that fail the gate.

**Schuemie MJ, Ryan PB, Pratt N, et al. (2020). Principles of large-scale evidence generation and evaluation across a network of databases (LEGEND). JAMIA 27(8):1331-1337.**
Ten principles including prespecified analysis design, dissemination of all results regardless of direction or significance, empirical evaluation via negative and positive control questions, and empirical calibration of estimates and confidence intervals using the control distribution.
Borrow: the dissemination principle (a SEVERE or inestimable verdict is itself a published result, preventing gate-driven publication bias) and empirical calibration of intervals using the negative-control stage cleanTMLE already runs.

**FDA Sentinel ARIA sufficiency determinations (Sentinel Initiative; see also: Six years of the US FDA's postmarket Active Risk Identification and Analysis system, PMID 37391385).**
Under FDAAA, before requiring a postmarket study FDA must first determine whether its ARIA surveillance system is "sufficient" for the safety question. Sufficiency is assessed on three axes, adequate data (exposure, outcome, covariates), appropriate methods, and satisfactory precision, and an insufficiency determination (197 made to date) triggers a different pathway.
Borrow: regulatory precedent for a formal, criterion-based, pre-analysis sufficiency verdict with defined consequences. The data/methods/precision triple is a tidy top-level structure for cleanTMLE's decision report, and "sufficiency determination" is respectable vocabulary for the gate's output.

**Gatto NM, Reynolds RF, Campbell UB (2019). A structured preapproval and postapproval comparative study design framework to generate valid and transparent real-world evidence (SPACE). Clinical Pharmacology and Therapeutics 106(1):103-115. Gatto NM, Campbell UB, et al. (2022). The structured process to identify fit-for-purpose data: a data feasibility assessment framework (SPIFD). Clinical Pharmacology and Therapeutics 111(1):122-134. Gatto NM, et al. (2023). A structured process to identify fit-for-purpose study design and data (SPIFD2). Clinical Pharmacology and Therapeutics 113(6):1235-1239.**
Stepwise templates that force articulation of the design (SPACE), then data feasibility with minimal acceptance criteria per data element (SPIFD), unified in SPIFD2 with explicit target-trial articulation and bias anticipation. The templates document a reasoned go/no-go on design and data before any analysis, aimed at "decision-grade" evidence.
Borrow: the notion of minimal criteria declared per design element in a reusable template, cleanTMLE's config file can mirror SPIFD2 fields so that each gate threshold is entered as a declared minimal criterion with a rationale string, making the audit trail self-documenting.

**ICH E9(R1) (2019). Addendum on estimands and sensitivity analysis in clinical trials. Primer: Kahan BC, Hindley J, Edwards M, Cro S, Morris TP (2024). The estimands framework: a primer on the ICH E9(R1) addendum. BMJ 384:e076316.**
The estimand framework defines an estimand by five attributes, population, treatment, endpoint (variable), intercurrent-event handling, and population-level summary, and separates sensitivity analyses (same estimand, different assumptions) from supplementary analyses (different estimand). The addendum notes the principles extend to observational studies.
Borrow: state every rung of cleanTMLE's ladder as a full five-attribute estimand, making explicit that moving from ATE to trimmed ATE to ATT to ATO changes exactly one attribute, the population, while the other four are held fixed. The sensitivity-versus-supplementary distinction also cleanly separates "same estimand, vary truncation" from "different estimand, moved down the ladder" in the output tables.

**Chen J, Scharfstein D, Wang H, Yu B, Song Y, He W, Scott J, Lin X, Lee H (2023+). Estimands in real-world evidence studies. arXiv 2307.00190 (ASA BIOP RWE working group; also a Springer chapter, Estimand in real-world evidence study: from frameworks to application).**
Extends E9(R1) estimand thinking to RWE, arguing RWD studies need additional considerations for population heterogeneity, complex treatment regimes, different intercurrent-event patterns, and endpoint complexities.
Borrow: the citation that estimand language is the correct frame for observational studies, plus its treatment of intercurrent events in RWD for the survival outcome stage of the trauma case study.

**Greifer N, Stuart EA (2021, revised 2023). Choosing the causal estimand for propensity score analysis of observational studies. arXiv 2106.10577.**
Practical guidance mapping each estimand (ATE, ATT, ATU, ATO) to its target population, assumptions, interpretation, and matching/weighting methods, written for medical researchers. The core message is that the estimand is chosen by the scientific question and its feasibility, and each carries distinct positivity requirements (ATE needs two-sided positivity, ATT one-sided, ATO essentially none).
Borrow: the estimand-by-assumptions table as the printed legend for cleanTMLE's ladder, in particular the one-sided-positivity requirement for ATT, which is the technically correct middle rung for the trauma case (treated patients all have support, untreated sickest do not).

**Rizk JG (2025). When and why to use overlap weighting: clarifying its role, assumptions, and estimand in real-world studies. Journal of Clinical Epidemiology, DOI 10.1016/j.jclinepi.2025.111942. PMID 40850393.**
A commentary warning against exactly the abuse cleanTMLE's governance prevents, overlap weighting delivers stable weights and exact balance but targets a statistically defined population that is hard to characterize clinically, so researchers should "define their target estimand before choosing a method" and not adopt ATO as a workaround when estimation gets hard, since that silently distorts the question.
Borrow: the counterargument the manuscript must answer. cleanTMLE's answer is procedural, the estimand switch is made before outcome access, on prespecified graded criteria, by a team that cannot see how the switch moves the answer, and is reported as an estimand change with the weighted population characterized in a table. Citing Rizk and answering it directly will strengthen the paper.

**Crump RK, Hotz VJ, Imbens GW, Mitnik OA (2006). Moving the goalposts (working paper, cited above) and Petersen et al. 2012 (topic 4).**
Together these are the two published articulations of estimand movement as a legitimate, systematic response to failed estimability, the econometric version (change the estimand by trimming) and the epidemiologic version (trade proximity to the initial target for identifiability, systematically).
Borrow: the ladder's two foundational citations, so "estimand ladder" can be introduced as a formalization of existing practice rather than a new invention.

## Synthesis: the 10 most actionable borrowings for cleanTMLE

1. Adopt the OHDSI objective-diagnostics gate as the quantitative core of PASS/FLAG/SEVERE (Conover et al. 2025 JAMIA). This is the single highest-value borrowing because it converts cleanTMLE's graded decision from a bespoke rubric into an instance of a published, empirically validated framework. Concretely, compute their five diagnostics at the pre-outcome stage, minimum detectable relative risk under 10, preference score in 0.3 to 0.7 for over half of each arm, maximum absolute standardized mean difference at most 0.10 after weighting, generalizability SMD under 0.25 against the target population, and EASE under 0.25 from the negative-control stage, and define PASS as all five within thresholds, FLAG as recoverable failures (balance or equipoise failures that a ladder move or design amendment could fix), and SEVERE as failures no estimand on the ladder repairs. Their empirical result, that gating removed essentially all systematic error in 11,716 negative-control analyses at the cost of declaring 86 percent of analyses inestimable, is also the best available evidence that gates of this kind work.

2. Anchor the gate's verdict grammar in Muntner's checkpoint language and the ARIA sufficiency precedent. Muntner et al. 2024 already published the three-way recommendation, proceed, conduct additional analyses, or terminate, issued by a masked review team at each stage, and FDA's ARIA process shows a regulator making formal sufficiency determinations on data, methods, and precision before committing to a study. cleanTMLE should present PASS/FLAG/SEVERE as an implementation of these published verdict structures, use "sufficiency" and "inestimable" as the formal output vocabulary, and structure the decision report's top level as the ARIA triple (adequate data, appropriate methods, satisfactory precision). This costs a paragraph of writing and buys the package a lineage.

3. Formalize the estimand ladder in ICH E9(R1) attributes with Petersen's tradeoff sentence as its principle. Each rung, ATE, trimmed ATE, ATT, ATO, and optionally Kennedy's incremental effect as the terminal rung, should be written out as a five-attribute estimand in which only the population attribute moves, citing Kahan et al. 2024 for the framework, Greifer and Stuart for the estimand-to-assumptions mapping (two-sided positivity for ATE, one-sided for ATT, none for ATO), Crump's "moving the goalposts" for the econometric origin, and Petersen et al. 2012 for the principle that the analyst trades proximity to the initial target for identifiability and should do so systematically. The ladder then needs its published counterweight answered, Rizk 2025 warns against ATO as a workaround, and cleanTMLE's reply is that the switch happens outcome-blind, under prespecified criteria, and is reported as an estimand change with the implied population characterized in a weighted Table 1.

4. Implement PoRT as the named-subgroup estimability diagnostic feeding the gate. Danelian et al. 2023 is CRAN-ready (RISCA::port) and returns human-readable covariate rules with exposure prevalence and subgroup size, which is what a masked review team can actually deliberate on, and what distinguishes structural from practical violations per Zhu et al. 2021's taxonomy (structural violations warrant eligibility changes or SEVERE, practical ones warrant trimming, truncation, or a ladder move). The decision report should print PoRT rules alongside the standard propensity histograms, tag each rule structural or practical after clinical review, and record the tag in the audit trail. OverRule (Oberst et al. 2020) is the citation for the general idea of rule-based overlap characterization and for separating population support from treatment overlap.

5. Rebuild the outcome-blind simulation stage on Nance et al. 2024 and guard it with Shaw et al. 2025. Nance et al. give the complete published recipe cleanTMLE's estimator-selection stage should follow, plasmode worlds preserving real covariate and exposure structure, candidate estimator grid over Super Learner libraries, truncation choices, and variance estimators, roughly 1000 replicates, and the two-step selection rule (minimize empirical variance subject to oracle coverage, then choose the variance estimator keeping nominal coverage), with the winner frozen into the SAP. Shaw et al. 2025 then constrains the simulator itself, the plasmode framework must preserve the estimand and the positivity structure or it will misrank estimators, which is maximally relevant when the whole point is a severe practical positivity violation. Citing the FIORD report (Dang et al. 2023) supplies the term outcome-blind simulation and the prespecified benchmark list (type I error, coverage, bias, precision).

6. Freeze the design as a two-stage outcome-free design report, borrowing the FDA CDRH mechanics. The device literature (Yue 2012; Yue, Lu, Xu 2014; Li et al. 2016; Li and Yue 2023) has run outcome-free stage-1 design with a firewalled independent statistician for over a decade inside actual regulatory submissions, with the design report signed and locked before outcome linkage, and stage-1 power assessment done on design-side quantities only. cleanTMLE should name the firewall roles in its config (who holds outcomes, who signs the freeze), emit a lockable design report artifact with a hash and timestamp, and cite Rubin 2008 for the principle that this report may legitimately conclude the dataset is inadequate, the SEVERE outcome has been a sanctioned endpoint of outcome-free design since its founding paper.

7. Use Walker's empirical equipoise as the first-line feasibility gate and report the truncation path per Cole and Hernan. Walker et al. 2013 supplies the field's oldest quantitative estimability gate (preference score 0.3 to 0.7 covering at least half of each arm), already adopted by OHDSI and now embedded in the Conover framework, and it is computable the moment the propensity stage finishes. Alongside it, the propensity report should include the two standard weight diagnostics from Cole and Hernan 2008, mean stabilized weight near one and a truncation-path plot showing the estimate's sensitivity to progressive truncation, since these are the diagnostics reviewers will look for first.

8. Make truncation algorithmic, with the Xu, Gruber, van der Laan defaults. Under practical positivity violations the outcome stage should default to clever-covariate-scaled targeting (loss-weighted targeting can be badly biased), truncation at c/(sqrt(n) log n) with c of 5, and targeted bootstrap variance when the gate has flagged positivity (arXiv 2604.20059), while offering the Gruber et al. 2022 data-adaptive truncation selector as the prespecified adaptive option and citing Leger et al. 2022 for applying truncation conditionally on diagnosed violation, because it costs bias when support is actually fine. This turns the most discretionary knob in TMLE practice into a documented, prespecified reaction.

9. Adopt the parallel ATT/ATO agreement check as the ladder's operational trigger. Ben-Michael and Keele's workflow, run balancing-weight ATT and overlap-weight ATO side by side, read agreement as evidence the ATT is estimable and divergence as evidence only the ATO is, gives the gate a concrete, prespecifiable comparison to act on, and it can be rehearsed outcome-blind in the plasmode stage before being executed once on real outcomes. Pair it with the PSweight-style convention of reporting effective sample size per candidate estimand, so the precision cost of each rung is visible in one table, and with Mao, Li, Greene's augmented estimator as the efficient analysis-stage counterpart once ATO is selected.

10. Package the firewall-crossing deliverable in STaRT-RWE/HARPER form and stamp FDA-style provenance. The design packet the outcome team receives should be a STaRT-RWE-style table set (design diagram, PICOT parameter tables, operational definitions) organized under HARPER section names, so it drops into a registerable protocol unchanged, and the package should auto-generate the attestation facts FDA's 2023 guidance and ICH M14 ask for, date of first data access, date of design freeze, analyses run before the freeze, deviations after it, all from its own audit log. Registering the frozen report (HMA-EMA catalogue or OSF) before outcome access completes the chain, and LEGEND's dissemination principle, publish the verdict even when it is inestimable, closes the loop against gate-driven publication bias.

## URLs consulted

- https://onlinelibrary.wiley.com/doi/10.1002/pds.5770
- https://api.semanticscholar.org/graph/v1/paper/DOI:10.1002/pds.5770 (and /citations)
- https://www.ebi.ac.uk/europepmc/webservices/rest/search (queries for DOIs 10.1002/pds.5770, 10.1002/pds.5396, 10.2147/CLEP.S608532, 10.1016/j.lungcan.2024.107960, 10.1002/pds.70096, 10.1002/cpt.3563, 10.1002/pds.70119, 10.1016/j.lungcan.2026.109469, 10.1080/10543406.2014.926367, and title searches for Petersen 2012, Cole and Hernan 2008, Kahan 2024, Rizk 2025, Desai npj 2025)
- https://pubmed.ncbi.nlm.nih.gov/38419140/
- https://pubmed.ncbi.nlm.nih.gov/34894377/
- https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=40387023
- https://www.dovepress.com/use-of-negative-control-outcomes-to-assess-cohort-comparability-among--peer-reviewed-fulltext-article-CLEP
- https://pmc.ncbi.nlm.nih.gov/articles/PMC11191559/
- https://pmc.ncbi.nlm.nih.gov/articles/PMC11924150/
- https://onlinelibrary.wiley.com/doi/10.1002/pds.70119
- https://www.nature.com/articles/526187a
- https://pubmed.ncbi.nlm.nih.gov/26450040/
- https://link.springer.com/article/10.1007/s11229-019-02456-7
- https://pubmed.ncbi.nlm.nih.gov/31541454/
- https://pubmed.ncbi.nlm.nih.gov/33327727/
- https://pubmed.ncbi.nlm.nih.gov/37097356/
- https://pmc.ncbi.nlm.nih.gov/articles/PMC8489282/
- https://pubmed.ncbi.nlm.nih.gov/33436424/
- https://onlinelibrary.wiley.com/doi/full/10.1002/pds.5507
- https://pubmed.ncbi.nlm.nih.gov/36215113/
- https://www.fda.gov/media/171667/download
- https://www.fda.gov/regulatory-information/search-fda-guidance-documents/real-world-evidence-considerations-regarding-non-interventional-studies-drug-and-biological-products
- https://www.fda.gov/media/177128/download
- https://www.federalregister.gov/documents/2024/03/21/2024-05969/real-world-evidence-considerations-regarding-non-interventional-studies-for-drug-and-biological
- https://www.foley.com/insights/publications/2024/04/fda-new-guidance-studies-drug-safety/
- https://www.federalregister.gov/documents/2026/03/04/2026-04253/m14-general-principles-on-planning-designing-analyzing-and-reporting-of-non-interventional-studies
- https://www.fda.gov/regulatory-information/search-fda-guidance-documents/m14-general-principles-planning-designing-analyzing-and-reporting-non-interventional-studies-utilize
- https://database.ich.org/sites/default/files/ICH_M14_Step4_Final_Guideline_2025_0905.pdf
- https://www.ema.europa.eu/en/ich-m14-guideline-general-principles-plan-design-analysis-pharmacoepidemiological-studies-utilize-real-world-data-safety-assessment-medicines-scientific-guideline
- https://pmc.ncbi.nlm.nih.gov/articles/PMC5157751/
- https://encepp.europa.eu/index_en
- https://www.ema.europa.eu/en/news/encepp-guide-methodological-standards-pharmacoepidemiology-revised
- https://www.mdpi.com/2227-9032/14/17/2870
- https://arxiv.org/pdf/0811.1640
- https://projecteuclid.org/journals/annals-of-applied-statistics/volume-2/issue-3/For-objective-causal-inference-design-trumps/10.1214/08-AOAS187.pdf
- https://doi.org/10.1080/10543406.2012.715111
- https://www.tandfonline.com/doi/full/10.1080/19466315.2016.1148071
- https://onlinelibrary.wiley.com/doi/10.1002/pst.2294
- https://api.semanticscholar.org/graph/v1/paper/DOI:10.1002/pst.2294
- https://onlinelibrary.wiley.com/doi/10.1002/pst.2295
- https://www.fda.gov/media/169060/download
- https://pmc.ncbi.nlm.nih.gov/articles/PMC8247341/
- https://www.nature.com/articles/s41746-025-02234-5
- https://pubmed.ncbi.nlm.nih.gov/24587587/
- https://www.sciencedirect.com/science/article/abs/pii/S0167947313003721
- https://arxiv.org/abs/2305.06028
- https://arxiv.org/abs/2504.11740
- https://pubmed.ncbi.nlm.nih.gov/39087681/
- https://pmc.ncbi.nlm.nih.gov/articles/PMC11444352/
- https://www.cambridge.org/core/journals/journal-of-clinical-and-translational-science/article/causal-roadmap-for-generating-highquality-realworld-evidence/3F30968E70E7A13EE7FC41A46A8C3AAD
- https://arxiv.org/abs/2205.08643
- https://arxiv.org/pdf/2305.07564
- https://pubmed.ncbi.nlm.nih.gov/21030422/
- https://arxiv.org/abs/2604.20059
- https://www.degruyterbrill.com/document/doi/10.1515/jci-2022-0032/html
- https://rdrr.io/cran/RISCA/man/port.html
- https://proceedings.mlr.press/v108/oberst20a.html
- https://arxiv.org/abs/1907.04138
- https://github.com/clinicalml/overlap-code
- https://academic.oup.com/biomet/article-abstract/96/1/187/235329
- https://www.iza.org/publications/dp/2347/moving-the-goalposts-addressing-limited-overlap-in-estimation-of-average-treatment-effects-by-changing-the-estimand
- https://academic.oup.com/aje/article-abstract/172/7/843/86816
- https://academic.oup.com/aje/article/190/8/1659/6146006
- https://www.dovepress.com/a-tool-for-assessing-the-feasibility-of-comparative-effectiveness-rese-peer-reviewed-fulltext-article-CER
- https://arxiv.org/pdf/1404.1785
- https://academic.oup.com/aje/article/188/1/250/5090958
- https://www.feinberg.northwestern.edu/sites/firstdailylife/docs/JAMA_Overlap_Weighting_A_Propensity_Score_Method_That_Mimics_Attributes_of_a_Randomized_Clinical_Trial.pdf
- https://arxiv.org/pdf/2010.08893
- https://journals.sagepub.com/doi/10.1177/0962280218781171
- https://arxiv.org/abs/2011.01388
- https://onlinelibrary.wiley.com/doi/10.1002/bimj.202300156
- https://journals.sagepub.com/doi/10.1177/0962280220940334
- https://doi.org/10.1177/09622802241269646
- https://arxiv.org/abs/2511.10077
- https://www.tandfonline.com/doi/abs/10.1080/01621459.2017.1422737
- https://arxiv.org/abs/1704.00211
- https://arxiv.org/pdf/2110.10532
- https://pubmed.ncbi.nlm.nih.gov/18682488/
- https://academic.oup.com/aje/article/191/9/1640/6580570
- https://arxiv.org/pdf/1707.05861
- https://arxiv.org/abs/2207.05010
- https://arxiv.org/abs/1711.02582
- https://www.sciencedirect.com/science/article/pii/S0304407620302694
- https://arxiv.org/abs/2110.14831
- https://arxiv.org/pdf/2210.01763
- https://projecteuclid.org/journals/annals-of-applied-statistics/volume-13/issue-2/Estimating-population-average-causal-effects-in-the-presence-of-non/10.1214/18-AOAS1231.full
- https://pmc.ncbi.nlm.nih.gov/articles/PMC8492528/
- https://onlinelibrary.wiley.com/doi/abs/10.1002/bimj.202000323
- https://arxiv.org/html/2502.11820v1
- https://arxiv.org/html/2502.14566v1
- https://arxiv.org/pdf/2501.11181
- https://academic.oup.com/jamia/article/32/3/518/7950905
- https://academic.oup.com/jamia/article/27/8/1331/5895561
- https://www.sentinelinitiative.org/drugs/ongoing-aria-assessments
- https://pubmed.ncbi.nlm.nih.gov/37391385/
- https://ascpt.onlinelibrary.wiley.com/doi/10.1002/cpt.1480
- https://ascpt.onlinelibrary.wiley.com/doi/10.1002/cpt.2466
- https://ascpt.onlinelibrary.wiley.com/doi/10.1002/cpt.2883
- https://pubmed.ncbi.nlm.nih.gov/38262663/
- https://arxiv.org/abs/2307.00190
- https://link.springer.com/chapter/10.1007/978-3-031-26328-6_9
- https://arxiv.org/abs/2106.10577
- https://www.jclinepi.com/article/S0895-4356(25)00275-6/fulltext
- https://pubmed.ncbi.nlm.nih.gov/40850393/

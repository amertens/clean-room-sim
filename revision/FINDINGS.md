# Findings: reported results that are wrong, unverifiable, or internally inconsistent

Stage 0 register, per operating rule 6. Each entry records the file and location, what is claimed, what the evidence shows, and the proposed fix. Nothing here has been corrected in the prose; corrections are assigned to work packages. Line numbers refer to the current sources at commit `863efc3` (cleanTMLE submodule `9f22d50`).

## F1. The Scenario C STOP is asserted, not produced by the executed pipeline

* **Where:** `reports/manuscript_outcome_blind_dq.qmd` lines 176-181 (abstract: "the prespecified pre-outcome decision correctly returned STOP"); 508-522 (Section 1.5: "the staged decision rule returns STOP and the locked primary analysis does not run ... both a realised-bias oracle and the outcome-blind screen flag it"); 2121-2131 ("The pre-outcome decision returns STOP for Scenario C"); 2396 and 2400-2401 (@tbl-workflow-contrast: "**STOP** under the prespecified unmeasured-U DQ scenario; the locked analysis does not run as planned", "The workflow blocks publication of the misleading null", "the workflow's prespecified rule refuses to authorise the Stage 4 analysis once the unmeasured-U row of the DQ stress test exceeds the locked bias threshold"). Presentation lines 217-229 ("the prespecified rule returns stop before any outcome is read"; "analyst likely publishes a misleading null ... Stop.").
* **Evidence:** the recorded pipeline artifacts show the executed gate returned GO and authorised outcome access in Scenario C: `results_new/audit_unmeasured_conf.csv` (Check Points 1-3 GO, Composite Gate GO), `results_new/decision_log_unmeasured_conf.csv` ("Outcome unblinding authorised by gate. All required checkpoints passed."), `done_unmeasured_conf.rds` meta. In `run_simulation.R` the gate is `gate_all(cp1, cp2, cp3)` (line 664-665); the DQ results never enter it and `gate_dq()` is never called. The single negative control cannot detect U by construction (`dgp_scenarios.R:64-66`: `nc_outcome` has no U or treatment term) and was not flagged (estimate -0.016, p = 0.44). The bias and coverage quoted as evidence (`scC_bias` 0.037, `scC_cov` 0.42) are realised Monte Carlo quantities against the true effect, unavailable pre-outcome.
* **Fix (WP2 item 2):** reframe every STOP claim as the prespecified tipping-point statement the abstract already half-makes; delete "The workflow blocks publication of the misleading null", "the prespecified pre-outcome decision correctly returned STOP", and the tbl-workflow-contrast Scenario C cell's execution claims; state in limitations that the executed composite gate in the recorded run authorised all four scenarios and that the DQ verdict was not wired into it. WP4 items 2 and 5 for the slides.

## F2. The same locked rule returns STOP for Scenario A

* **Where:** manuscript lines 2090-2104 (locks `max_abs_bias = 0.02`, `min_coverage = 0.90`, then presents Scenario A as the clean GO benchmark); 2744-2786 (E-value section prints Scenario A's unmeasured-U bias ladder without comment).
* **Evidence:** recorded DQ rows for the selected candidates (`results_new/dq_stress_good_overlap.rds`, `dq_stress_unmeasured_conf.rds`): Scenario A bias 0.0264 (OR 3.0) and 0.0417 (OR 4.0), coverage 0.60 and 0.30; Scenario C 0.0325 and 0.0534, coverage 0.433 and 0.267. Both scenarios breach the locked bounds from OR 3.0 upward; at OR 2.0 both pass on bias (0.0103 vs 0.0135). No restriction of the locked grid separates them on bias; the only differing cell is OR 2.0 coverage (0.90 vs 0.833), two replicates out of thirty. The scenarios also carry different selected candidates (aggressive vs robust).
* **Fix (WP2 item 2, run-or-reframe):** the tipping-point framing must state the symmetric result: under the locked grid, the planned analysis crosses the bias tolerance between OR 2.0 and OR 3.0 in every scenario, including good overlap. If a GO/STOP contrast between A and C is retained at all, it must come from a corrected, prespecified rule executed on both (WP1 item 3), accepting whatever verdict that rule yields for both.

## F3. Blinding claims for Stage 1b are contradicted by the code

* **Where:** `cleanTMLE/vignettes/cleanTMLE-staged-analysis.qmd:594-595` ("reads only marginal event counts (never the treatment-outcome association)"); manuscript line 860 (stage table: Stage 1b "Events only") and 1048 (dossier table: "Marginal Y only") and 1715 ("(marginal event counts only)"); `cleanTMLE/vignettes/cleanTMLE-functions.qmd:181-182` (labels `summarize_event_support()` "Marginal Y only" while describing "Event counts by arm").
* **Evidence:** `estimate_design_precision()` computes per-arm events and crude per-arm rates and returns and prints them (`cleanTMLE/R/stages.R:963-973, 1024`); `summarize_event_support()` returns per-arm rates (stages.R:1084-1097). Live probe output in STAGE0_AUDIT Section 2.2.1. Per-arm event rates are the crude treatment-outcome association.
* **Fix (WP1 item 1):** restrict `estimate_design_precision()` to totals and the marginal rate; move per-arm counts behind an explicit `event_support_by_arm()` with a logged, warned access path; correct the vignette and manuscript sentences to match whichever behaviour is implemented.

## F4. The manuscript claims plasmode DGP modes that do not exist

* **Where:** manuscript lines 1010-1013 ("the package supports both modes, recording which mode was used"); 1637-1642 ("Depending on configured mode, the plasmode generator may use ... or it may take E[Y|W] from external pilot data").
* **Evidence:** `run_plasmode_feasibility()` and `run_plasmode_dq_stress()` have no mode argument and no mode recording; the only knob is `q0_library` (learner, not data source); both always fit Q0 on the lock's real outcome and refuse masked locks. The manuscript's own Future extensions item 5 (lines 3779-3784) lists the DGP-spectrum modes as not yet built. `simulate_support()` has `q0_source` options, none of which is an external-pilot mode for the candidate-selection loop.
* **Fix:** either implement `dgp_mode` with a locked field (WP1 item 2) or rewrite the two passages to describe only what exists (run-or-reframe; record the branch in DECISIONS_PENDING).

## F5. Scenario C is described with the pre-recalibration parameters

* **Where:** manuscript lines 1973-1975 ("prevalence 0.20 and odds ratios 2.0 for both A and Y"); 2764 (E-value table row label "True RD (U present, OR = 2.0)").
* **Evidence:** `run_simulation.R:146-153` and every recorded Scenario C artifact use `U_prevalence = 0.40, U_trt_OR = 3.0, U_out_OR = 3.0` (recalibrated July 2026; `results_new/_regen_scenarioC.log`, `sandbox/scenarioC_recalibration/`). Truth in the artifacts: -0.0419, versus -0.0298 in A/B.
* **Fix (WP2):** correct the design prose and the E-value row label from the artifact values; regenerate any dependent statements.

## F6. Section 12.6's candidate table and DQ numbers cannot be reproduced from the repository

* **Where:** manuscript lines 3210-3223 (Stage 2b feasibility table: glm_t01 +0.00001 / 0.0134 / 0.92; ensemble_t05 +0.00076 / 0.0135 / 0.92; glmnet_t01 +0.00019 / 0.0135 / 0.92; ensemble_t01 +0.00082 / 0.0136 / 0.92, captioned "50 replicates per candidate"); 3228-3231 (worst-case RMSE ratios 1.82 / 1.73 / 1.53 / 1.27, "the largest contributor for each candidate is the near-positivity threat at slope x4"); 3387-3388 ("Stage 2b feasibility coverage of 0.94 with bias below 0.001 across 50 replicates"); reproducibility table line 3961 (rescueco_dq_full.R -> plasmode_dq_stress.rds, plasmode_dq_degradation.csv).
* **Evidence:** `rescueCo/results/plasmode_metrics.csv` (50 reps) contains two candidates (glm_t01, glm_t05) with different values (bias 0.00016, RMSE 0.01306, coverage 0.94); `plasmode_summary.csv`, `stage2b_dq_stress.rds`, and `plasmode_dq_stress.rds` all contain a four-candidate 3-replicate smoke run (mtimes 2026-07-29) whose values differ again. The four-candidate 0.92-coverage table matches no file in the repository. The RMSE ratios match only `_dq_degradation_50rep_backup.csv` (1.818 / 1.731 / 1.525 / 1.268; mtime 2026-06-07), a file the reproducibility table never mentions; the file it does name now holds 3-rep values (worst ratios 1.14 to 1.22). The backup also shows ensemble_t05's worst case at slope x3.0, not x4. `rescueCo/logs/pipeline.log` shows a 200-replicate rerun launched 2026-06-07 with no completion line, and the 3-rep smoke run of 2026-07-29 overwriting the primary artifacts. The prose coverage 0.94 and the table coverage 0.92 disagree with each other and trace to different artifacts.
* **Fix (WP2 item 7, run-or-reframe):** rerun `rescueco_dq_full.R` at its default 200 replicates (estimated about six hours; requires the restricted registry data, present on this machine) and rebuild Section 12.6 from the fresh outputs via chunks; or remove the 50-replicate numbers and state that the case study omits the DQ stage. In either branch, the smoke-run artifacts must stop shadowing the reported ones (separate file names or a `reps` stamp in the file name).

## F7. The released case-study metadata contradicts the version claim

* **Where:** manuscript lines 2913-2916 ("All stages were run on cleanTMLE 0.2.0, with the package version and the lock fingerprints recorded in the released case-study metadata").
* **Evidence:** the only tracked metadata file, `rescueCo/results/manuscript_artifacts/case_study_metadata.json`, records `"cleanTMLE_version": "0.1.5"` (rendered 2026-07-29, the earlier single-contrast run). The tracked multi-arm CSVs (`lock_summary.csv` and the rest) record neither a package version nor a lock fingerprint. The multi-arm run of 2026-09-13 has no released metadata file. Additionally, Section 12.6's content comes from the 0.1.5-era run.
* **Fix (WP2/WP3):** emit a multi-arm metadata file (version, submodule commit, lock hashes) from script 10 or 14, track it, and make the sentence point at it; scope the 0.2.0 claim to the stages that were actually rerun, or rerun 12.6 under 0.2.0 (see F6).

## F8. "Only bounded-weight estimands" for time to arrival is contradicted by the released ladder table

* **Where:** manuscript lines 3339-3342 ("time to arrival ... is estimated here only through bounded-weight estimands"). Compare slide "The disaster outcome, re-analysed" (presentation lines 371-390), which prints the ladder rows and explicitly narrates "the override-only ATE row shows its strain in the interval width".
* **Evidence:** `rescueCo/results/multiarm/ladder_estimates.csv` contains ATE rows for `time_to_arrival_minutes` on C1 (complete case; the logged reconciliation override), C3, and C4; only C2 lacks one. The slide matches the artifact; the manuscript sentence does not. (For six-month mortality the caption at line 3293 is accurate: C2 has no ATE row and C1's is the override.)
* **Fix (WP2/WP4 item 3):** rewrite the sentence to say the C1 time-to-arrival ATE appears only as the logged, SEVERE-labelled reconciliation override while the substantive rows are the bounded-weight estimands; keep manuscript and slide aligned to the CSV.

## F9. The marginal-overlap variance slide contradicts the variance study it summarises

* **Where:** `reports/cleanTMLE_presentation.qmd` lines 211-215: "IPTW confidence intervals come out too wide ... Matched TMLE has the same problem, roughly doubling the true spread. Both are fixable with a bootstrap variance, and the workflow says so before unblinding."
* **Evidence:** `results_new/bootstrap_variance.csv`: under marginal overlap Match_TMLE SE/SD (IF) is 1.041 with coverage 0.97 (approximately calibrated, the manuscript's own reading at lines 2204-2210); IPTW SE/SD is 1.229 with bootstrap 1.220 (the bootstrap does not repair it; manuscript lines 2192-2203 and 2300-2303 say so); the anti-conservative Match_TMLE regime is good overlap (SE/SD 0.931, coverage 0.89, bootstrap raising it to about 0.95). "Roughly doubling" matches nothing in the variance study; the four-scenario run's marginal Match_TMLE SE/SD is 1.452, a different analysis and still not a doubling.
* **Fix (WP4 item 1):** rewrite the slide from `bootstrap_variance.csv`: matched TMLE undercovers under good overlap and is repaired by the re-matching bootstrap; it is calibrated under marginal overlap; IPTW over-covers under marginal overlap and the bootstrap does not fix it; plain TMLE and cross-fitted TMLE undercover under marginal overlap.

## F10. Kent (2026) is an unpublished conference presentation cited for quantitative results

* **Where:** `reports/references.bib:56-62` (`kent2026fiord`: "Conference presentation, 4th Annual Forum on the Integration of Observational and Randomized Data (FIORD), Day 1 Lunch session"); cited at manuscript lines 470, 858, 933, and 3595-3604, the last with quantitative results ("the final risk difference at four years was reported as roughly -4%").
* **Evidence:** the entry is not a citable publication. A peer-review-track source for the evolocumab example exists and resolves: Jin R, Hurwitz K, Kent ST, et al. Abstract 4363106, Circulation 2025;152(Suppl_3) (ahajournals.org/doi/10.1161/circ.152.suppl_3.4363106). Whether the abstract carries the staged-decision details the manuscript attributes to the example (negative-control flag in the health-seeking domain; review team proceeding conditional on an E-value QBA plan) is not yet verified: [[VERIFY]] before swapping.
* **Fix (WP2 item 4):** cite the Circulation abstract for the study and its result; attribute the staging-framework narrative to Muntner et al. and, if the details survive verification only in the talk, cite the presentation separately and only for what a presentation can support.

## F11. The paper's own reproducibility map points at overwritten or misplaced artifacts

* **Where:** manuscript reproducibility table lines 3950-3963; `run_gate_operating_characteristics.R:40`.
* **Evidence:** the gate-OC script writes `results/gate_oc.rds` while the manuscript reads `results_new/gate_oc.rds` (the artifact, mtime 2026-05-26, was moved by hand); the DQ row points at files now holding smoke-run values (F6); `run_simulation.R` in the working tree (modified 2026-09-13) postdates the artifacts it is credited with (Scenarios A/B/D June 8-9, C July 11), which were produced by cleanTMLE 0.1.x, before the 0.2.0 API the paper documents.
* **Fix (WP0/WP2):** point the script at `results_new/`, regenerate or re-date the provenance claims, and either rerun the four-scenario study under 0.2.0 or state the version that produced it.

## F12. Version and default inconsistencies in the methods prose

* **Where:** manuscript line 3686 ("cleanTMLE 0.2.1 adds a q0_library argument") versus DESCRIPTION Version 0.2.0 and line 433 ("produced under version 0.2.0"); line 1582 (decision thresholds "recorded ... as part of the lock fingerprint") versus `.compute_lock_hash()` covering no thresholds (a separate `thresholds_hash` exists only on the superseded path); lines 1573-1576 naming an "ESS floor at 30% of n" default that appears nowhere in the code (the superseded balance checkpoint uses 50 %).
* **Evidence:** `cleanTMLE/DESCRIPTION:3`; `cleanTMLE/R/cleanroom.R:186-195, 433-443`; no match for a 30 % ESS threshold in `cleanTMLE/R/`.
* **Fix (WP2):** correct the version reference, describe the fingerprint contents accurately (and extend the hash to cover data values and thresholds in WP3), and state the thresholds that exist.

## F13. The negative-control ladder can silently substitute unadjusted estimates for TMLE

* **Where:** `cleanTMLE/R/design_tools.R:180-198`; case-study text (manuscript lines 3153-3176) presents the unadjusted ladder as primary with the TMLE variant in the supplement.
* **Evidence:** with `method = "tmle"` (the default), any per-rung TMLE failure falls back to the unadjusted two-proportion estimate with `status = "estimated"` and no per-row record of the method actually used; a mixed table is indistinguishable from a pure one. Separately, no locked decision criteria exist for the ladder (per-rung p < alpha only), the gap WP1 item 5's `nc_criteria` addresses, and the primary/supplement assignment inverts Muntner Section 11 (adjusted population, primary estimator).
* **Fix (WP1 items 5; WP2 item 5):** record the realised method per row and never fall back silently; default the ladder to the locked primary estimator; declare `nc_criteria` on the lock; swap the primary and descriptive roles in Sections 6.6 and 12.4.

## F14. The stress test scores the full-cohort ATE even after the ladder has moved

* **Where:** `cleanTMLE/R/plasmode_dq.R` (`truth <- mean(p1_sim) - mean(p0_sim)` on the lock's cohort); `tutorials/clean_room_targeted_learning_tutorial.qmd`, Stage 2b.
* **Evidence:** on the tutorial's SEVERE design every candidate reads STOP with the tipping point at the first grid level; the cached object shows baseline coverage 0.60 to 0.63 and mean SE at 43 to 52 percent of the empirical SD before any threat is applied, while bias under the threats never exceeds 0.020. The verdict restates the support failure of an estimand the ladder has already abandoned; nothing lets the stress test score the trimmed ATE the ladder moved to.
* **Fix (done 2026-09-19):** the tutorial says what the full-cohort verdict measures and then runs `stress_test()` a second time on the common-support cohort the ladder selected (the trimmed ATE is the ATE on that cohort): baseline coverage 0.93 against 0.67 to 0.70, GO for every candidate with the worst coverage on the 0.90 floor, no tipping point inside the declared grid; the candidate is locked from that run and the design report reads it. A package-level `estimand` argument on the stress test is a feature beyond the WPs and is not proposed here.

## F15. Hybrid-mode leakage warning silenced in the tutorial; blinding contract overstated

* **Where:** `tutorials/clean_room_targeted_learning_tutorial.qmd` (`warning = FALSE` in the setup chunk; the sentence on the blinding test suite); `cleanTMLE/tests/testthat/test-blinding.R` (identity on masked and permuted locks tested only for `dgp_mode = "external_pilot"`).
* **Evidence:** the tutorial cohort's propensity c-statistic is 0.864, above the 0.80 threshold at which `.warn_hybrid_separation()` warns that the covariate-only Q0 approximates the crude association; the warning fired and was suppressed. The tutorial and the manuscript both said the suite asserts identical output for every design verb; the hybrid stress test reads the real outcome for its baseline surface.
* **Fix (done 2026-09-19):** both documents state the qualification, and the tutorial now declares `dgp_mode = "external_pilot"` on both locks with a `pilot_q0` fitted on a separately simulated pilot cohort (n = 2000, a different seed), so every design verb reads no outcome from the study cohort.

## F16. Re-declaring the ladder to attach the selected candidate duplicates the design-log entry

* **Where:** `cleanTMLE/R/estimand_ladder.R` (`declare_estimand_ladder(candidate = )` logs a full `estimand_ladder` entry on every call); the tutorial's exported log shows the ladder twice.
* **Evidence:** `lock_primary_tmle_spec()` is internal in 0.3.0, so the second declaration is the sanctioned route, and the Muntner-format export carries two identical ladder rows.
* **Fix (deferred, outside the WPs):** log the candidate attachment as an amendment entry; the tutorial explains the duplicate for now.

## F17. Manuscript statements that must be revisited when the reruns land

* **Where:** `reports/manuscript_outcome_blind_dq.qmd`: "The case-study estimates in @sec-casestudy were produced under version 0.2.0" (Current scope); "All stages were run on cleanTMLE 0.2.0" (case study); the reproducibility map's `results_new/` rows; the 50-replicate case-study DQ subsection.
* **Evidence:** the WP2 prose pass moved the manuscript's description of the package to 0.3.0 without touching any result-bearing claim; these sentences remain true of the recorded artifacts and false of the eventual 0.3.0 reruns.
* **Fix (WP2, after the reruns):** update the version statements and the replicate counts from the rerun outputs through the source map; `UNMATCHED_NUMBERS.csv` stays the work list.

## F18. `plot.plasmode_dq_results()` uses the deprecated `aes_string()`

* **Where:** `cleanTMLE/R/plasmode_dq.R`, `plot.plasmode_dq_results()`.
* **Evidence:** every call emits the ggplot2 3.0.0 deprecation warning for `aes_string()`, which the tutorial's `warning = FALSE` hides.
* **Fix (done 2026-09-19):** replaced with `aes(.data$...)`, matching the package's other plots; all three metrics build warning-free, `test-plasmode_dq.R` passes (27 assertions), reduced check clean.

# Decisions pending

Stage 0 register of decisions that belong to the author, with the information needed to make each one. Work packages do not proceed past these points without a recorded choice; run-or-reframe branches taken later will be recorded here per operating rule 3.

## D1. Paper split (WP2 item 1)

Choose between: (a) two papers, with Paper 1 carrying the estimability layer, lock and design log, Scenarios A, B, D, and case-study Sections 12.1-12.5 and 12.7, and Paper 2 carrying the DQ stress test, min_max_rmse, the gate calibration (Section 11), the FIORD two-stage selector, and the variance-method selection; or (b) one paper with Sections 7, 8, 10.3-10.5, and 11 moved to a supplement and summarised in at most 600 words. Relevant fact from the audit: the parsed body is 22,207 words against a 9,000-word target, so option (b) requires removing roughly 60 percent of the body; option (a) is the easier fit.

## D2. Scenario C reframing and its symmetric consequence (WP2 item 2; FINDINGS F1, F2)

The reframing to a tipping-point sensitivity screen is prescribed, but one consequence needs an explicit decision: under the locked grid the tipping point is crossed between OR 2.0 and OR 3.0 in Scenario A as well. The manuscript can either (a) present the screen with no GO/STOP asymmetry between A and C (the verdict is a property of the locked grid, stated identically for both), or (b) implement the corrected pure gate (WP1 item 3), run it on both scenario fixtures, and report whatever it returns, accepting that this may be STOP for both. The prompt's default is (b) with no grid tuning; confirm.

## D3. Case-study DQ stage: rerun or remove (WP2 item 7; FINDINGS F6)

Branch (run): execute `rescueCo/scripts/rescueco_dq_full.R` at its default 200 replicates, estimated about six hours on this machine, requiring the restricted registry Stage 1/2 artifacts already on disk; rebuild Section 12.6 from the fresh outputs via chunks. Branch (reframe): remove the 50-replicate numbers and state that the case study omits the DQ stage. The 50-replicate values currently survive only in `_dq_degradation_50rep_backup.csv`, and the Stage 2b candidate table matches no file at all, so leaving Section 12.6 as it stands is not an option.

## D4. Four-scenario simulation provenance (FINDINGS F11)

`results_new/` was produced in June and July 2026 under cleanTMLE 0.1.x; the working-tree `run_simulation.R` postdates it. Either (a) rerun the four scenarios under 0.2.0 (the summaries print an approximately twelve-minute per-scenario runtime claim, but the recorded June run took several hours end to end; budget half a day), or (b) keep the artifacts and state the producing version and script commit in the reproducibility section. This interacts with D2(b): a corrected gate needs fixtures consistent with whatever run is reported.

## D5. Ethics and consent placeholders (WP2 item 7)

Manuscript lines 2928-2931 need the approving ethics committee, approval number, and the consent basis (consent or waiver). Nothing in the repository supplies these. Provide the values, or direct that the placeholders become `[[FILL]]` markers and the case study be described accordingly.

## D6. Kent citation replacement (WP2 item 4; FINDINGS F10)

Swap `kent2026fiord` for the Circulation 2025 abstract (Jin, Hurwitz, Kent, et al., A4363106) for the evolocumab example. Pending verification: whether the staged-decision details (negative-control flag, conditional proceed with an E-value QBA plan) appear in the abstract text; if not, decide whether to keep a separate presentation citation for those details or drop them.

## D7. cleanroomGov fold-back (WP3 item 2)

cleanroomGov exports nine functions, is used by one script family, and is invisible in the manuscript. The WP3 criterion (fewer than ten exports) says fold it back into cleanTMLE. Confirm, including whether `attach_estimand` and `declare_sensitivity_plan` become exported cleanTMLE verbs or lock fields set by `create_analysis_lock()` arguments.

## D8. Plasmode DGP modes: implement or reframe (WP1 item 2; FINDINGS F4)

The manuscript twice claims investigator-specified and external-pilot Q0 modes exist for the candidate-selection and DQ loop; the code has only real-Y Q0 with a learner knob. WP1 item 2 says implement the three modes as locked fields. Confirm implementation (and its scope: both `run_plasmode_feasibility` and `run_plasmode_dq_stress`), or direct a reframe of the two passages to match the code.

## D9. Superseded governance layer removal (WP3 item 2)

The four-scenario simulation script and its recorded artifacts depend on the superseded checkpoint/audit layer via `:::`. Removing the layer (WP3) breaks re-running `run_simulation.R` as written. Decide the order: rebuild the simulation on the current verbs first (folds into D4a), or keep the superseded layer through 0.3.0 with the simulation pinned to it.

## D10. UNMATCHED_NUMBERS semantics for WP0 (operating rule 2)

The Stage 0 matcher matches at value level and therefore reports zero unmatched tokens while FINDINGS F6 documents genuinely unverifiable numbers. Confirm the WP0 `build_source_map.R` specification: claim-level assertions (each reported table cell and inline statistic mapped to a named file and field, failing the build when absent), with the Stage 0 value-level map retained only as a discovery aid.

## D11. Per-arm event counts access path (WP1 item 1)

WP1 prescribes `event_support_by_arm()` gated by `enforce = FALSE` or `reviewer_role = "data_manager"`, a design-log entry, and a warning. Confirm the gate style (argument-based versus role-based), since WP1 item 7's `roles` field is records-only and provides nothing to check against.

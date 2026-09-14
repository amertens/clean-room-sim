# Decisions register

Stage 0 raised D1-D11; all eleven were decided by the author on 2026-09-14 at Stage 0 approval. Each entry below records the question in brief and the recorded decision. Later run-or-reframe branches will be appended here as they are taken (operating rule 3). Execution order for the work packages, as directed: WP0, WP1, WP3, rerun all scripts, WP2, WP4.

## D1. Paper split. Decided: (a) two papers

Paper 1 carries the estimability layer, the lock and design log, Scenarios A, B, D, and case-study Sections 12.1-12.5 and 12.7. Paper 2 carries the DQ stress test, `min_max_rmse`, the gate calibration (Section 11), the FIORD two-stage selector, and the variance-method selection. **Scenario C moves to Paper 2 entirely.**

## D2. Scenario C gate. Decided: (b) pure gate, tipping point primary

Implement the pure gate (WP1 item 3), run it on all four rerun fixtures, and report **the tipping-point odds ratio per scenario as the primary output**, with the verdict secondary. No grid tuning; whatever the corrected rule returns for each scenario is what is reported.

## D3. Case-study DQ stage. Decided: run branch, 200 replicates, after WP3

Rerun at 200 replicates under 0.3.0, after WP3. Every output filename carries the replicate count. The 3-replicate smoke artifacts are renamed so they can never shadow a reported run.

## D4. Four-scenario simulation provenance. Decided: (a) rerun under 0.3.0

Rerun the four scenarios, the gate operating characteristics, and the variance study under 0.3.0, from scripts rebuilt on the current verbs. All design prose describing the configurations is machine-written from the run configuration (no hand-typed DGP parameters).

## D5. Ethics placeholders. Decided: [[FILL]] markers

Insert `[[FILL]]` markers for the approving committee, approval number, and consent basis; record the case study as unsubmittable until the values are supplied.

## D6. Kent citation. Decided: replace with two published sources

Cite Jin et al., Circulation 2025;152(Suppl_3):A4363106 for the Stage 4 estimate, and Dhalwani et al., Epidemiology 2024;35(4):579-588 for the Stage 3 negative-control-outcome feasibility. Drop the review-team conditional-proceed details. (Resolve and verify both citations at WP2 per the DOI rule.)

## D7. cleanroomGov. Decided: fold back

Fold cleanroomGov into cleanTMLE. `attach_estimand` and `declare_sensitivity_plan` become `create_analysis_lock()` arguments rather than exported verbs.

## D8. Plasmode DGP modes. Decided: implement two modes

Implement `dgp_mode` with values `"hybrid"` and `"external_pilot"` in both `run_plasmode_feasibility()` and `run_plasmode_dq_stress()`, as locked fields. Reframe `"investigator_specified"` as future work in the prose.

## D9. Superseded governance layer and old artifacts. Decided: rebuild and archive

Rebuild the simulation on the current verbs. Move the old scripts and the June/July artifacts to `archive/` with a README; delete nothing.

## D10. Source map. Decided: claim-level matcher

`build_source_map.R` asserts each reported value against a named field of a named file (implemented in WP0 with `revision/source_map_claims.R`). An empty UNMATCHED_NUMBERS.csv is never reported until it is true; the file currently lists every uncovered literal and is the WP2 work list. `--strict` (zero unmatched) is the WP2 completion gate.

## D11. Per-arm event counts. Decided: argument-based access

`event_support_by_arm(lock, reason = )` with a mandatory `reason`, a design-log entry recording the access, and a warning that arm-specific counts reveal the crude association.

## Additions to WP1 (directed at Stage 0 approval)

1. The lock fingerprint hashes design-data content, not only dimensions and column names. (Design data per WP3 item 3: covariates, treatment, missingness indicators, negative controls; the primary outcome stays out of the content hash so the fingerprint remains outcome-blind.)
2. The negative-control ladder must error or return an explicit adjusted/unadjusted column per row; it never falls back silently (finding F13).
3. Remove the internal `set.seed(1)` in `estimand_feasibility()`.
4. DQ tables report every candidate, with the selected candidate flagged.

## Additions to WP2 (directed at Stage 0 approval)

1. The reconciliation table states the parallel pipeline's data cut and package version and reports all verdict counts (agree, benign_difference, missing_in_cleanroom, missing_in_main, not_comparable).
2. Section 12.1 carries the analyst disclosure (the case-study analyst previously conducted the parallel unblinded analysis; the case study demonstrates re-execution and reconciliation under data-level masking, not personnel blinding).

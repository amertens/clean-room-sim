# Archive

Nothing in this directory is deleted; it is retired. Every file here
was superseded during the September 2026 revision (cleanTMLE 0.3.0)
and is kept verbatim for provenance. Do not run these scripts against
the current package and do not cite these result files in the
manuscripts; the reported runs live in `results_new/` (rebuilt under
cleanTMLE 0.3.0) and `rescueCo/results/`.

## Contents

- `scripts_0.1x/` - the original simulation drivers
  (`run_simulation.R`, `run_simulation_full.R`, `dgp_scenarios.R`,
  `run_gate_operating_characteristics.R`, `summarize_results.R`,
  `check_progress.R`, `_bootstrap_variance.R`, `TMLE_Demo.R`). They
  call the 0.1.x API (checkpoints, audit logs, gate tokens) that
  cleanTMLE 0.3.0 removed. The revision decision record (revision
  decision D4, option a) replaced them with drivers rebuilt on the
  sixteen-verb surface.

- `results_new_0.1x/` - the May-July 2026 simulation artifacts
  produced by those drivers under cleanTMLE 0.1.x (four scenario
  runs, gate operating characteristics, bootstrap variance study,
  and the `_bak_scenarioC_preRecalib/` backup). Retired under
  revision decisions D4 (rerun all scenarios under 0.3.0) and D9
  (move, do not delete). The Stage 0 audit finding that the recorded
  per-scenario decision logs say GO while the manuscript claimed
  STOP is documented in `revision/FINDINGS.md`; these files are the
  evidence and must stay unmodified.

- `cleanroomGov/` - the short-lived companion package from the 0.2.0
  estimation-vs-governance split. Folded back into cleanTMLE in
  0.3.0 (revision decision D7): `attach_estimand` and
  `declare_sensitivity_plan` became `create_analysis_lock()`
  arguments, and the governance-note formatters returned as cleanTMLE
  exports. `build_stage_manifest` and `summarize_stage_path` consumed
  the deleted audit-log objects and did not return.

- `PRIMARY, C1/`, `C2, C3, C4/` - stray pipeline logs from early
  multi-arm case-study runs.

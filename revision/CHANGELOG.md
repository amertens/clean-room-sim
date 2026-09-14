# Revision changelog

One entry per work package: what changed, what was verified by running code, and what remains pending. Nothing is claimed as complete that was not run. Commits are named per repository (superproject clean-room-sim; submodule cleanTMLE).

## Stage 0 (audit; superproject commit 7010e5d)

Read-only audit producing `revision/STAGE0_AUDIT.md`, `FINDINGS.md` (F1-F13), and `DECISIONS_PENDING.md`. Verified by running: as-cran check (2 ERRORs, 6 WARNINGs, 4 NOTEs), `devtools::test()` (662 pass, 0 fail), coverage (68.3 percent), clean-session vignette renders (5 min 14 s and 1 min 9 s), live blinding probes, and inspection of every recorded artifact behind the manuscript's gate claims. Approved with decisions D1-D11 on 2026-09-14.

## WP0 (tooling; superproject commit e7e2dc3, submodule commit 3d0b2e8)

* `revision/style_check.R`: failure-class and warning-class token counts by file with a scoped allowlist; `--strict` gate. Run: baseline recorded (121 "synthetic" across sources at the time; zero em-dashes and contractions in the qmds).
* `revision/build_source_map.R` + `revision/source_map_claims.R`: claim-level source map (decision D10). Run: 68 of 68 seeded claims verified against artifacts; 27 chunk-read rows; honest `UNMATCHED_NUMBERS.csv` with 183 uncovered literals (148 manuscript, 35 presentation), the WP2 work list.
* `revision/render_all.R`: renders manuscript, vignettes, slides. Run: validated on the slides target (12.5 s) to a scratch directory; the manuscript target will be exercised by WP2, the vignette targets were shown to render in Stage 0.
* `cleanTMLE/tests/testthat/test-blinding.R`: the blinding contract as tests. Run: 22 assertions passing at WP0 time, with the two known violations and the missing generator mode encoded as three explicit skips for WP1.
* Decisions D1-D11 recorded in `revision/DECISIONS_PENDING.md`.

## WP1 (blinding and gate corrections; submodule commit pending in this entry's commit, superproject pointer follows)

Package changes (all verified by the full test suite: 759 pass, 0 fail, 3 package-presence skips, 11 pre-existing experimental-scope warnings; the three WP0 skips in `test-blinding.R` are now live passing assertions):

* Item 1 and D11: `estimate_design_precision()` and `summarize_event_support()` are marginal-only; `checkpoint_cohort_adequacy()` (superseded) loses its per-arm rows; new `event_support_by_arm(lock, reason = )` with mandatory reason, warning, and design-log entry. Vignette sentences corrected.
* Item 2 and D8: `dgp_mode` (`"hybrid"`, `"external_pilot"`) as a locked, fingerprinted field read by both plasmode functions; `pilot_q0` contract; hybrid-mode c-statistic warning above 0.80 with the leakage statement in the docs; plasmode results no longer carry the lock (hash and mode instead). Verified: external-pilot runs are identical on reference, masked, and permuted locks.
* Item 3 and D2: `dq_thresholds` as a locked field; `dq_locked_verdict()` and `dq_tipping_points()` as pure functions of (locked thresholds, declared grid, metrics). Verified on the recorded Scenario A and C fixtures: the same locked grid (bias 0.02, coverage 0.90, RMSE ratio 1.5) reads STOP for both; tipping points OR 3.0 (A, bias) and OR 2.0 (C, coverage floor). No grid tuning.
* Item 4: `refine_ps_after_nco()` removed (decision recorded in NEWS; the vignette row replaced with the removal note).
* Item 5 and F13: `nc_criteria` as a locked field; `nc_ladder_verdict()` grades rungs per domain from the locked criteria only; the ladder never falls back silently (failed TMLE rows stay failed TMLE rows; per-row `method` and `domain` columns); `run_negative_control_tmle()` errors instead of substituting the IPTW helper and defaults to the locked candidate's Q library.
* Item 6: `design_report()` stores summary statistics only (per-patient g and A vectors stripped, marker recorded) and takes `dq = ` so the locked DQ verdict, tipping point, and Check Point 3 reading enter the recommendation.
* Item 7: `roles` recorded on the lock and printed in the report header.
* Item 8: `export_design_log(lock, format = "muntner")` writes the Muntner Table S1 columns; ladder switches and overrides fill the structured fields.
* Additions: the lock fingerprint now digests design-data content (outcome excluded; masked locks validate; tampered covariates fail; legacy locks validated against the legacy fingerprint with a message); `set.seed(1)` removed from `estimand_feasibility()` (deterministic sorted matching order); `summarize_dq_degradation(selected = )` flags the selected candidate while always reporting all.
* Style: the package tree (R sources, tests, vignettes, DESCRIPTION, NEWS) is clean of "synthetic" (now "simulated"), em-dashes, and contractions; `plots.R`'s `highlight` variable is an allowlisted identifier. `man/` regenerated; version bumped to 0.2.0.9000 (0.3.0 at WP3).

Pending after WP1: WP3 (API surface, seeds sweep, parallel path, pkgdown, as-cran zero-zero, 0.3.0), the reruns (D3, D4), WP2, WP4.

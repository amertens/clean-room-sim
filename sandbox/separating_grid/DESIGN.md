# DESIGN: realistic separating-candidate grid for Section 9.4

Status: **executing (2026-06-08).** The two-arm spec below was approved at Stage 1;
the sections that follow are the original design. The executed configuration was
then narrowed by the decisions in the revision log directly below, which is the
source of truth for what actually ran. The config is fingerprinted and written to
disk *before* the sweep, and the run reports whatever it produces, including a null.

## Revision log (executed configuration, supersedes the two-arm spec below)

Discoveries during build forced changes to the approved two-arm design. Each was
decided with the user; nothing about the severity-locking discipline changed.

1. **Arm B dropped (run Arm A only).** The SuperLearner arm is not just slow, it is
   numerically unstable: `SL.gam` *and* `SL.glmnet` hard-crash the R session (exit
   127, not a catchable error) after roughly a hundred fits in one process,
   confirmed by isolation tests (`_debug3.R`, `_debug5.R`, `_timecells.R`). The
   `min_max_rmse`-flips-the-estimator claim is fundamentally about truncation under
   positivity stress, which is pure GLM propensity weighting and needs no
   SuperLearner. So the study runs **Arm A only** (marginal overlap, GLM propensity,
   truncation lever); the library-richness story (Arm B) is explicitly out of scope.
2. **Arm A grid widened to five GLM truncations** {0.005, 0.01, 0.025, 0.05, 0.10},
   `SL.glm` propensity and outcome, no cross-fitting. This spans the efficient-but-
   fragile to robust-but-less-efficient ends under marginal overlap.
3. **Replication: 5 batches x 500 reps** (up from 200) for tighter between-batch
   Monte Carlo SE on the worst-case RMSE. The GLM path makes this affordable.
4. **Architecture: controller + per-batch worker processes.** Given the session
   crash history, the full run is split so each batch runs in a fresh R process that
   writes its cell to `cells/` and exits (memory freed, crash-isolated, resumable).
   The controller locks the config first, the five batch workers run concurrently,
   then an aggregate phase combines the cells. Driver:
   `separating_grid_study.R` (phases: `lock` / worker / `aggregate`).
5. **Near-positivity, calibrated:** for Arm A the calibration drops *both* candidate
   slopes (1.5 and 2.0 push the already-marginal PS to effectively deterministic
   treatment, outside the study's own overlap envelope), so the locked Arm A threat
   set is `regulatory_standard` minus near-positivity = covariate missingness (MCAR,
   3 levels), treatment misclassification (2), outcome misclassification (2), and
   unmeasured confounding (paired, 2): **9 degraded cells plus baseline**. (The
   unmeasured-confounding levels are enumerated element-wise, `plasmode_dq.R:797`,
   not as a 2x2 grid.)
6. **Package bug fixed along the way:** the NA-safety crash in
   `run_plasmode_dq_stress` (see `docs/revision/package_bugfixes.md`).

Locked-config fingerprint for the executed run: `d6b2bde2…` (full value in
`locked_config.rds` / `.yaml`).

## Purpose

Section 9.4 currently establishes only existence: the `min_max_rmse` rule changes
the locked estimator on a constructed grid whose severities were "chosen to
exhibit the mechanism" (`reports/manuscript_outcome_blind_dq.qmd:2235-2241`). This
study replaces that with a grid an analyst would plausibly prespecify and threat
severities anchored to external evidence and locked before the run. The result may
be a divergence between the selection rules or a null in which they agree. We
report whichever happens and never tune severities upward to manufacture a flip.

## Design logic: where separation can come from (and where it must not)

The separation must be a property of a realistic study regime, not of the threat
severities. There are two independent structural levers, but a verified property
of the generator forces them into **two separate arms**: the `misspec = TRUE`
branch (`run_simulation.R:136-151`) ignores `overlap_strength` and uses its own
fixed propensity structure that "does not pile the PS at 0/1." A probe
(`_probe_overlap.R`, fitted-GLM PS on n = 8000) confirms the two regimes are
mutually exclusive:

| config | induced PS (q01, q99; min, max) | outcome surface | true RD |
|---|---|---|---|
| `misspec=FALSE, overlap=1.5` | (0.105, 0.963; 0.015, 0.994) -> marginal overlap | linear | -0.030 |
| `misspec=TRUE,  overlap=1.5` | (0.184, 0.854; 0.080, 0.939) -> good overlap | nonlinear (quadratic age, sex x biomarker, sex effect modification) | -0.080 |

The task asked for both marginal overlap (so truncation matters) and a nonlinear
surface (so library richness matters). Because the generator cannot deliver both
at once and we may not modify it, the study runs **two arms**, each isolating one
lever, scored under the **same** locked threats and the **same** 5-candidate grid:

- **Arm A (truncation lever):** `misspec = FALSE, overlap_strength = 1.5`,
  generator `q0_library = NULL` (the truth is linear, so a linear generator is the
  honest choice). Marginal overlap makes aggressive vs conservative truncation
  trade variance against positivity bias. The library candidates are included but,
  under a linear truth, are expected to be non-separating; we report that plainly.
- **Arm B (library/cross-fitting lever):** `misspec = TRUE`, generator
  `q0_library = c("SL.glm","SL.gam","SL.glm.interaction","SL.mean")` so `p_base`
  captures the nonlinear truth (`cleanroom.R:1269-1309`, `:1271-1275`). Good
  overlap makes truncation roughly idle; outcome-library richness and
  cross-fitting are the live axis.

Each arm's lever is fixed by its regime, before any threat is applied. The threats
are then mild, externally anchored, and locked, identically across both arms.

## 1. Candidate grid (5 candidates)

Built with `tmle_candidate()` (`cleanTMLE/R/cleanroom.R:929-999`). Axes a real
prespecified grid contains: outcome/PS library richness, cross-fitting, and PS
truncation. No single-knob extremes (no truncation 0.001, no odds ratio 8).

| id | g_library (PS) | q_library (outcome) | cv_scheme | truncation |
|---|---|---|---|---|
| `glm_t01`    | `SL.glm` | `SL.glm` | none | 0.01 |
| `glm_t05`    | `SL.glm` | `SL.glm` | none | 0.05 |
| `ens_t01`    | `SL.glm, SL.glmnet, SL.gam` | same | none | 0.01 |
| `ens_t025`   | `SL.glm, SL.glmnet, SL.gam` | same | none | 0.025 |
| `ens_cf_t01` | `SL.glm, SL.glmnet, SL.gam` | same | `cv_tmle` (V = 5) | 0.01 |

Justification: `glm_t01` is the efficiency-first default many analysts start from;
`glm_t05` is the same learner with a more conservative truncation; the `ens_*`
pair is a standard three-learner SuperLearner that an analyst adds when they
suspect nonlinearity, at two realistic truncations; `ens_cf_t01` adds
cross-fitting, the textbook remedy for over-fit nuisance bias. Truncation spans a
realistic 0.01 to 0.05.

Dependency note: `SL.glmnet` needs `glmnet`, `SL.gam` needs `gam`, `cv_tmle` is a
candidate field consumed by the package. All are in `Suggests`. The Stage-2 script
will assert these are installed and stop with a clear message otherwise, rather
than silently dropping a learner (which would quietly change the grid).

## 2. Data-generating process (exact, per arm)

Both arms: `n = 4000`, `effect_size = -0.05`, `seed = <batch seed>`; truth via
`compute_truth(n_truth = 100000, ...)` (the engine's own large-sample Monte Carlo,
`run_simulation.R:192-226`). The two DGP functions are extracted from
`run_simulation.R` by parsing the file and evaluating only the `generate_data` and
`compute_truth` assignments, so the engine's main loop never runs (verified).

- **Arm A:** `generate_data(n=4000, overlap_strength=1.5, misspec=FALSE, ...)`;
  generator `q0_library = NULL` (linear truth -> linear generator).
- **Arm B:** `generate_data(n=4000, overlap_strength=1.5, misspec=TRUE, ...)`;
  generator `q0_library = c("SL.glm","SL.gam","SL.glm.interaction","SL.mean")`
  (matches the manuscript's Scenario D library, `run_simulation.R:260`).

`n = 4000` is a tunable choice that keeps SuperLearner + cross-fitting plasmode
runtime feasible; it is part of the locked config.

## 3. Threats: externally anchored, locked severities

We use the package's pre-existing literature-anchored preset,
`default_dq_scenarios("regulatory_standard")` (`cleanTMLE/R/plasmode_dq.R:43-52`),
which already encodes mid-range, not extreme, severities and maps onto the
manuscript's own calibration sources in @sec-calibration (`manuscript:1165-1187`).
The locked severities are therefore:

| threat | locked severities | external anchor (manuscript's own) | note |
|---|---|---|---|
| Outcome misclassification | (sens, spec) = (0.95, 0.99), (0.90, 0.95) | phenotype/algorithm-validation literature [@fda2024considerations] | real test |
| Unmeasured confounding | U OR in {1.5, 2.0}, prevalence 0.20 | routine E-value range [@vanderWeele2017evalue], [@diaz2013sensitivity] | real test; deliberately NOT the 3-8 strong-U range of the current demo |
| Covariate missingness (MCAR) | fractions {0.05, 0.10, 0.20} | observed-rate / data-profiling [@kahn2016data; @fda2024ehrclaims] | **near-null by construction**: median imputation under MCAR is unbiased (`manuscript:2082`); reported as such, not co-equal |
| Near-positivity | slope candidates {1.5, 2.0}, each locked per arm only if it passes the calibration rule below | the manuscript's realistic overlap floor (PS to ~0.04-0.99, `manuscript:1645`) | calibrated, not assumed (see rule) |

Near-positivity calibration rule (transparent, data-driven, run before locking).
The threat amplifies the fitted PS about its logit mean,
`ps_pos = plogis(mean(lp) + slope * (lp - mean(lp)))` (`plasmode_dq.R:915-918`).
For each arm the script fits the same covariate GLM PS used in the probe, applies
each candidate slope, and **keeps a slope only if the induced PS keeps its 1st-99th
percentile within [0.01, 0.99] and its min/max non-deterministic (>= 0.002,
<= 0.998)**. Slopes that push PS to effectively deterministic treatment are dropped
and recorded as dropped. This is expected to bite hardest in Arm A, which is
already at the overlap edge: if both slopes push Arm A beyond the envelope, Arm A's
near-positivity threat is dropped and the report says so, rather than dialing
positivity past the study's own realistic range. The kept slopes per arm are
written into the locked config before the scoring sweep.

Two honesty points written into the design:
- The MCAR missingness threat is flagged as a near-null test, not presented as a
  co-equal fourth threat.
- Near-positivity is anchored to the study's own overlap, not to deterministic
  treatment. Stage 2 prints the induced PS range for each slope and locks only the
  slopes that remain inside the marginal-overlap envelope. This calibration step
  happens before the scoring sweep and its outcome is recorded in the config.

Citation TODO for you: the bracketed keys above are the manuscript's existing
bibliography keys. If any anchored range needs a more specific primary citation
than the manuscript currently carries, that is left as a TODO for you to fill; the
study will not invent a citation.

## 4. Decision thresholds (locked)

`gate_dq()` package defaults (`plasmode_dq.R:1162`): `max_abs_bias = 0.02`,
`min_coverage = 0.85`, `max_rmse_ratio = 1.5`. Effective sample size from the PS
fit is recorded as a diagnostic where the fit surfaces it; if it is not separately
surfaced as a gate input, the config says so plainly rather than implying an ESS
gate that is not enforced.

## 5. Selection rules compared

`select_tmle_candidate(sim_results, rule, dq_results)` (`cleanroom.R:1567-1717`):
- `min_rmse`: lowest baseline plasmode RMSE (from `run_plasmode_feasibility`).
- `min_max_rmse`: lowest worst-case RMSE across the locked non-baseline threats
  (consumes `dq_results`; `cleanroom.R:1593-1604`).
- `fiord_two_stage`: baseline coverage within tolerance of 0.95, then smallest
  mean SE.

The headline question is whether `min_rmse` and `min_max_rmse` select the same
candidate; `fiord_two_stage` is reported alongside.

## 6. Replication and Monte Carlo SE

- 5 independent batches x 200 plasmode reps per cell (1000 total per cell), batch
  seeds spaced by 1e6, matching the existing demo design
  (`sandbox/candidate_divergence/REPORT.md`).
- `effect_sizes = c(0.05)` (single plausible target RD; keeps worst-case RMSE
  interpretation clean; tunable in the config).
- Worst-case RMSE per candidate is computed within each batch, then averaged
  across batches; MC SE = between-batch SD / sqrt(5). A reported divergence is
  tested against this MC SE (separation must exceed roughly 2 combined MC SE to be
  called real).
- Runtime guard: Stage 2 runs a smoke test (1 batch x 20 reps) first and prints an
  estimated full-run time. If SL.gam + cv_tmle make 5 x 200 prohibitive, the
  documented fallback is 5 x 100; the choice is recorded in the config before the
  full run.

## 7. Locking discipline

A single R list `config` (grid spec, DGP args, locked severities including the
calibrated near-positivity slopes, thresholds, rules, effect sizes, reps, batch
seeds) is serialized to `sandbox/separating_grid/locked_config.rds` and a
human-readable `sandbox/separating_grid/locked_config.yaml`, and fingerprinted with
`digest::digest()`, **before** the scoring sweep. RESULTS.md records the
fingerprint and confirms (by file mtime) that the config predates the results file.
The near-positivity PS-range calibration is the one step that runs before locking;
its result is written into the config that is then fingerprinted.

## 8. Outputs (all new, under sandbox/separating_grid/)

- `locked_config.rds`, `locked_config.yaml` (written before the sweep)
- `separating_grid_results.rds`, `separating_grid_results.csv` (tidy)
- `RESULTS.md` (Stage 3: computed numbers only)
- No existing manuscript output or `results_new/` file is touched.

## 9. Stage 3 reporting contract (honest outcome)

RESULTS.md contains only computed numbers and:
- A table: candidate, axes, baseline RMSE (MC SE), worst-case RMSE (MC SE),
  binding threat, and which rule selects which candidate.
- A one-line verdict: do `min_rmse` and `min_max_rmse` diverge at the anchored
  severities?
- If divergence: the worst-case separation between the two selected candidates and
  its combined MC SE, whether it exceeds noise, and confirmation the config was
  locked before the run.
- If null: stated plainly, with no severity escalation. Then, quantitatively from
  the run, how far the binding threat would have to move for a flip, and whether
  that severity is inside or beyond the externally anchored range, with a
  recommendation on whether Section 9.4 should stay framed as an existence
  demonstration.

## Open choices for you before Stage 2

1. Grid: accept the 5 candidates, or adjust the library/truncation/cv axes.
2. Replication budget: 5 x 200 (default) versus 5 x 100 fallback, pending the
   Stage-2 runtime smoke test.
3. Sample size `n = 4000`: accept or change.
4. Whether to keep `effect_sizes = c(0.05)` or add a second target (e.g. 0.03).

# Scenario C recalibration: making detection and fragility align

## The problem being fixed

The manuscript's Scenario C reports a realised TMLE bias of ~0.007 under a true
unmeasured confounder of OR ≈ 2.0, then presents the workflow's STOP as
"blocking a biased estimate the conventional pipeline would publish." But 0.007
is **below** the paper's own 0.02 bias threshold, with ~0.93 coverage, so on the
realised metrics the analysis passes. Under an oracle defined on the *realised*
estimator's bias (as the contribution text specifies), Scenario C is therefore a
**false positive**: the oracle returns GO while the screen returns STOP. The STOP
is produced only by the top rung of the DQ severity sweep (confounding stronger
than the truth), i.e. it is a *fragility* signal dressed up as a *detection*.

Detection and fragility align only when, at the designated true confounder
strength, the realised analysis is *materially* biased (oracle STOP) **and** the
outcome-blind screen also flags it (screen STOP).

## Method

`sandbox/scenarioC_recalibration/scenarioC_smoke.R` reuses the DGP and the
oracle-vs-screen machinery of `run_gate_operating_characteristics.R` and adds the
missing piece: **realised-TMLE Wald coverage**. It sweeps confounder prevalence
`U_prev` × strength `s` (the OR of `U` on both treatment and outcome), computing,
per cell:

- `bias_real`, `cover_real`  → oracle STOP if `|bias_real| > 0.02` (detection)
- `bias_screen_sl`           → screen STOP if `|bias_screen| > 0.02` (fragility)

`s` maps to an RR-scale E-value via `EV = s + sqrt(s(s-1))`, so each strength has
a defensible "plausible residual confounding" reading.

Run: `n = 2000`, `reps = 60`, linear surface, threshold 0.02, seed 2026.

## Results (60 reps, linear surface)

| U_prev |  s  | E-value | true_rd | bias_real | coverage | screen_sl | oracle | screen |
|:------:|:---:|:-------:|:-------:|:---------:|:--------:|:---------:|:------:|:------:|
| 0.3 | 1.0 | 1.00 | — | ~0.000 | ~0.95 | ~0.00 | 0 | 0 |
| 0.3 | 1.5 | 2.37 | — | <0.02 | ~0.93 | <0.02 | 0 | 0 |
| 0.3 | **2.0** | **3.41** | 0.132 | **0.0228** | 0.82 | 0.0204 | 1 | 1 |
| 0.3 | **2.5** | **4.44** | 0.131 | **0.0410** | 0.57 | 0.0320 | 1 | 1 |
| 0.3 | 3.0 | 5.45 | 0.129 | 0.0549 | 0.30 | 0.0521 | 1 | 1 |
| 0.3 | 3.5 | 6.46 | 0.128 | 0.0668 | 0.15 | 0.0597 | 1 | 1 |
| 0.5 | 2.0 | 3.41 | 0.133 | 0.0248 | 0.77 | 0.0265 | 1 | 1 |
| 0.5 | 2.5 | 4.44 | 0.131 | 0.0496 | 0.37 | 0.0378 | 1 | 1 |
| 0.5 | 3.0 | 5.45 | 0.128 | 0.0671 | 0.15 | 0.0595 | 1 | 1 |
| 0.5 | 3.5 | 6.46 | 0.126 | 0.0797 | 0.07 | 0.0644 | 1 | 1 |

(s = 1.0 and 1.5 cells: realised bias below threshold, oracle GO — correctly *not*
aligned.)

## Recommendation

**Use U_prev = 0.3, s = 2.5 (E-value 4.44) as the revised true Scenario C.**

At this cell the realised analysis is genuinely broken — bias 0.041 (2× the
threshold), coverage 0.57 — so the conventional pipeline publishes a materially
biased, under-covering estimate, and both the oracle and the outcome-blind screen
return STOP. Detection and fragility agree, and the E-value of 4.44 is a
defensible "plausible residual confounding" magnitude to justify from a
fit-for-purpose review.

Why not the weakest aligned cell (s = 2.0)? At s = 2.0 the realised bias (0.0228)
and screen bias (0.0204) sit right on the 0.02 threshold; the verdict is
Monte-Carlo-fragile and could flip at higher reps. s = 2.5 clears the threshold
decisively. (Use s = 3.0 for an even more emphatic demonstration: bias 0.055,
coverage 0.30.)

Notably, at the manuscript's stated OR ≈ 2.0 the realised bias in this DGP is
already ~0.023, i.e. right at the boundary — the original "0.007" must have come
from weaker effective confounding (different `U_prev`/coefficients in
`run_simulation.R`). Either way, moving the designated true strength to s ≥ 2.5
converts Scenario C from a false-positive illustration into an honest true
positive.

## Accompanying gate-OC reframing

1. **Define the oracle on the realised estimator's bias** (as the contribution
   text already says). With the revised Scenario C at s = 2.5, `stop_oracle = 1`,
   so the screen's STOP is a **true positive**, not a disagreement.
2. **Report the screen's operating characteristics as literal numbers.** Extend
   `run_gate_operating_characteristics.R` to print sensitivity, specificity, and
   false-STOP rate against the realised-bias oracle across the full s-grid,
   separating the sub-threshold cells (s ≤ 1.5, oracle GO — where a STOP is a
   false positive, the *cost* of screening) from the supra-threshold cells
   (s ≥ 2.5, oracle STOP — where a STOP is power). This gives the reader the
   power/false-positive trade-off the current text omits.
3. **Key the DQ severity range to an external E-value argument that brackets the
   truth.** Sweep s over, say, {1.5, 2.0, 2.5, 3.0} (E-values 2.4–5.5) and state
   that the range is the fit-for-purpose review's plausibility bound for residual
   confounding. The screen STOPs because a rung *within the plausible range*
   breaches the threshold, not because an arbitrary top rung does.
4. **Drop the "conventional pipeline publishes bias 0.007" framing** and replace
   it with the s = 2.5 contrast: the conventional pipeline publishes RD-bias 0.041
   at 57% coverage; the staged workflow returns STOP before the outcome is read.

## Caveats / next steps

- 60 reps gives coverage MCSE ≈ 0.06; confirm the recommended cell at 200 reps
  before it goes in the manuscript (`SC_REPS=200`).
- Only the **linear** surface was swept (`SC_NL=0`), which is the right base case
  for the alignment argument (GLM and SL screens agree). The nonlinear surface
  (`SC_NL=1`) is where the GLM screen under-predicts and the SL screen is needed;
  that is a separate point about screen *fidelity*, not detection/fragility
  alignment, and can be reported alongside.
- The DGP here is the gate-OC script's (`W1` continuous, `W2` binary), not
  `run_simulation.R`'s five-covariate DGP. To put these numbers in the manuscript
  as "Scenario C," port the recommended `U_prev`/`s` into `run_simulation.R`'s
  unmeasured-confounding generator and re-confirm bias/coverage there.

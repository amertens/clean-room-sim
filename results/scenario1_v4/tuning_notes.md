# Scenario 1 v4 — Tuning Notes

## Goal

Demonstrate that an overlap flag can be a **false stop** for TMLE *inferential*
performance (not just point estimation) in the weak-W4 mode, while being a
**true stop** in the strong-W4 mode.

## Key Design Choices

### 1. Mode-Specific Overlap Severity

v3 used the same `s_bad = 2.0` for both modes, which made the overlap
violation equally severe regardless of whether W4 mattered for outcomes.
v4 introduces mode-specific `s_bad`:

- **Weak:** `s_bad = 1.50` — moderate PS separation. Enough to trigger the
  overlap flag, but not so severe that even irrelevant-W4 settings break
  inference.
- **Strong:** `s_bad = 2.0` — severe PS separation, matching v3.

**Tuning range:** If weak-mode false-stop fails, decrease `s_bad_weak` toward
1.35. If it succeeds too easily, increase toward 1.60.

### 2. Moderated Nonlinearity

v3 used `h1 = 0.6*W1 + 0.4*W1^2` and `h3 = 0.5*|W3|`. v4 reduces these to:

- `h1 = 0.50*W1 + 0.25*W1^2`
- `h3 = 0.35*|W3|`

This keeps the misspecification signal (GLM vs truth) present but reduces the
difficulty of learning the outcome surface, helping TMLE-ML achieve good
coverage.

### 3. Reduced Treatment Heterogeneity (Weak Mode)

v3 weak mode: `tau_int1 = 0.25`, `tau_int2 = 0.10`
v4 weak mode: `tau_int1 = 0.15`, `tau_int2 = 0.05`

Less effect modification means the ATE is more homogeneous across covariate
strata, reducing the difficulty of the estimation problem in the false-stop
scenario.

### 4. NP Bootstrap as Primary Inference

The v3 multiplier bootstrap perturbs the EIF without refitting nuisance models.
At N=200-400, this can produce under-coverage because nuisance estimation
uncertainty is non-negligible. v4 uses full NP bootstrap (resample + refit
Q and g + retarget) as the primary inferential check.

- **Fast preset:** B_np = 75 (development)
- **Eval preset:** B_np = 150 (final evaluation)

### 5. Sample Size and Replicates

| Preset | N    | Reps | B_np | B_mult | V_cf |
|--------|------|------|------|--------|------|
| fast   | 400  | 120  | 75   | 200    | 3    |
| eval   | 600  | 300  | 150  | 300    | 5    |

### 6. Overlap Thresholds (Slightly Tighter than v3)

| Threshold      | v3   | v4   | Rationale |
|----------------|------|------|-----------|
| extreme_prop   | 0.35 | 0.30 | More sensitive flag |
| ess_frac       | 0.25 | 0.30 | More sensitive flag |
| max_w          | 50   | 40   | More sensitive flag |

Tighter thresholds make the false-stop scenario more common (more replicates
flagged) but also more realistic — a practical analyst would use conservative
thresholds.

## Expected Results

### Weak W4 (False Stop)

- **Overlap flag rate:** ~40-60% of replicates flagged as "bad"
- **TMLE-ML bias (bad overlap):** |bias| < 0.01
- **TMLE-ML NP bootstrap coverage (bad overlap):** ≥ 0.93 (target: ~0.93-0.96)
- **Gate result:** acceptable = TRUE → demonstrates false stop

### Strong W4 (True Stop)

- **Overlap flag rate:** ~40-60% of replicates flagged
- **TMLE-ML bias (bad overlap):** may be > 0.01
- **TMLE-ML NP bootstrap coverage (bad overlap):** likely < 0.90
- **Gate result:** acceptable = FALSE → demonstrates true stop

## Coverage MCSE Reference

For the fast preset (120 reps, ~60 per stratum):

- True coverage 0.95 → MCSE ≈ sqrt(0.95*0.05/60) ≈ 0.028
- True coverage 0.90 → MCSE ≈ sqrt(0.90*0.10/60) ≈ 0.039

For the eval preset (300 reps, ~150 per stratum):

- True coverage 0.95 → MCSE ≈ sqrt(0.95*0.05/150) ≈ 0.018
- True coverage 0.90 → MCSE ≈ sqrt(0.90*0.10/150) ≈ 0.024

## If Tuning Fails

If the weak-mode false-stop doesn't materialize (coverage in bad-overlap
stratum is too low):

1. **Decrease `s_bad_weak`** from 1.50 toward 1.35 (milder overlap violation)
2. **Increase N** from 400 toward 500-600 (more data helps TMLE)
3. **Decrease `beta4_Q`** below 0.04 (make W4 even more irrelevant to Y)
4. **Add more SL learners** (e.g., SL.gam, SL.randomForest)

If the strong-mode true-stop is too easy (coverage is already terrible):

1. **Decrease `s_bad_strong`** from 2.0 toward 1.8
2. **Increase `tau_int3`** above 0.20 (stronger A×W4 interaction)

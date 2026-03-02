# Analysis Description: HCV Treatment and Acute Kidney Injury (AKI) Clean-Room Simulation

## Overview

This repository implements a **clean-room simulation study** that estimates the causal effect of **sofosbuvir-based hepatitis C virus (HCV) treatment** on the risk of **acute kidney injury (AKI)**. It uses **Targeted Maximum Likelihood Estimation (TMLE)** within a clean-room design framework, where distinct teams handle data generation, method specification, quality assurance, and reporting in isolation to reduce researcher degrees of freedom and bias.

The simulation is written in R and structured as a realistic pharmacoepidemiologic study with time-to-event outcomes, informative censoring, treatment switching, and optionally non-proportional hazards.

---

## Repository Structure

```
clean-room-sim/
├── R/
│   ├── dgp/
│   │   └── DGP.R                  # Data-generating process (core simulation engine)
│   └── tmle/
│       └── tmle_pipeline.R        # TMLE estimation pipeline (placeholder)
├── simulation/
│   └── run_simulations.R          # Simulation driver (placeholder)
├── analysis/
│   └── analysis_plan.qmd          # Analysis plan integration (placeholder)
├── reports/
│   ├── README.md
│   └── TL-SAP.qmd                # Targeted Learning Statistical Analysis Plan
├── tests/
│   └── test_sim.R                 # Tests (placeholder)
├── tutorial code/
│   ├── tmle_clean_room_tutorial.qmd   # Quarto tutorial: TMLE in clean-room design
│   └── frengression_example.Rmd       # R Markdown: calling Python frengression from R
├── data/                          # Generated simulation datasets (CSV)
└── clean-room-sim.Rproj           # RStudio project file
```

---

## Data-Generating Process (DGP)

The core of the analysis is in `R/dgp/DGP.R`, which defines the function `generate_hcv_data()`. This function simulates a realistic observational cohort study of HCV patients, producing a dataset that mirrors the structure of administrative claims or EHR data.

### Population

- **Default cohort size**: N = 125,000 individuals
- **Demographics**: age (normal, mean 48, SD 13, minimum 18), sex (58% male), race (5 categories), region (4 US regions), enrollment duration (Poisson, mean 420 days)

### Baseline Covariates (34 binary/continuous variables)

The simulated patients have a rich set of clinical characteristics:

| Category | Variables |
|---|---|
| Renal/hepatic | `ckd`, `prior_aki`, `cirrhosis`, `portal_htn`, `esld` |
| Cardiovascular | `heart_failure`, `hypertension` |
| Infections | `hiv`, `sepsis` |
| Metabolic | `diabetes`, `bmi`, `overweight_obese` |
| Behavioral | `smoking`, `alcohol`, `substance_abuse` |
| Cancer | `cancer`, `chemo` |
| Nephrotoxic medications | `nsaid`, `acearb`, `diuretic`, `aminoglycoside`, `contrast` |
| Other medications | `statin`, `aspirin`, `beta_blocker`, `ccb`, `art` |
| Treatment history | `prior_sof`, `prior_nonsof` |
| Other | `dehydration`, `obstruction` |

### Cohort Eligibility Criteria

Patients are excluded if:
- Enrollment < 365 days
- Age < 18
- Prior AKI event
- Prior sofosbuvir or non-sofosbuvir HCV treatment

### Treatment Assignment (Propensity Score Model)

- **Treatment**: binary indicator for sofosbuvir-based regimen (target prevalence ~36%)
- **Simple model**: logistic function of age, cirrhosis, CKD, HIV, diabetes, cancer, plus noise
- **Complex model** (when `complexity = TRUE`): adds non-linear terms — `bmi^2`, `sin(bmi)`, `(age/50)^3`, `ckd*cancer` interaction, `hiv*log(age)` interaction
- Propensity scores are clipped to [0.05, 0.95] to maintain positivity
- Supports counterfactual overrides: `"all_treated"` or `"all_control"` for computing true causal effects

### Outcome Model (AKI Event Times)

- **Outcome**: time to first acute kidney injury event
- **Baseline hazard**: `h0 * exp(linear predictor)` where the linear predictor depends on age, CKD, cirrhosis, heart failure, NSAIDs, contrast agents
- **Simple model**: proportional hazards with a single HR for treatment
- **Complex model** (when `complexity = TRUE`): adds `age^2`, `bmi^2`, `sin(bmi)`, `heart_failure*acearb` interaction, `nsaid*treatment` interaction, `contrast*log(age)` interaction

#### Non-Proportional Hazards (when `np_hazard = TRUE`)

The treatment effect varies over time using a piecewise exponential model with a change-point at `tau = 45` days:
- **Early period** (0–45 days): HR = 1.25 (elevated risk, possibly due to drug initiation effects)
- **Late period** (>45 days): HR = 0.70 (protective effect)

This creates a realistic scenario where the treatment effect is not constant over follow-up, challenging standard Cox regression assumptions.

### Censoring Mechanisms

#### Administrative Censoring

- **Independent** (when `dep_censor = FALSE`): exponential with rate 1/100 days
- **Dependent** (when `dep_censor = TRUE`): censoring rate depends on baseline risk (linear predictor) and treatment, violating independent censoring assumptions
- Hard administrative cutoff at 720 days (max follow-up)

#### Treatment Switching Censoring (when `switch_on = TRUE`)

- Treatment switching is modeled as a survival process with its own hazard function
- Switch hazard depends on treatment assignment (`gamma_A = 0.60`) and CKD status (`gamma_ckd = 0.40`)
- Patients are censored at switch time + 30-day risk window
- Creates informative censoring tied to both treatment and patient characteristics

### Observed Data Construction

- **Follow-up time** = min(event time, administrative censoring time, switch censoring time)
- **Event indicator** = 1 if event occurred before any censoring, 0 otherwise
- Internal latent variables (true event times, censoring times) are dropped from the analysis dataset

### Optional Missingness

- Can introduce missing data in `region` (5% MCAR) and `ckd` (10% MCAR)
- Optional imputation via `missForest` (random forest-based multiple imputation)

### Two Pre-Configured Simulation Scenarios

1. **Simple scenario** (`sim_hcv_aki.csv`): proportional hazards, independent censoring, linear models, `h0 = 3e-4`
2. **Complex scenario** (`sim_hcv_aki_complex.csv`): non-proportional hazards, dependent censoring, non-linear models with interactions, `h0 = 5e-5`

---

## Estimation: Targeted Maximum Likelihood Estimation (TMLE)

### Target Parameter

The **average treatment effect (ATE)**: the difference in mean potential outcomes under treatment vs. control:

```
ATE = E[Y(1)] - E[Y(0)]
```

where `Y(a)` is the potential outcome under treatment assignment `a`.

### TMLE Procedure (from tutorial)

The tutorial in `tutorial code/tmle_clean_room_tutorial.qmd` demonstrates TMLE using the `tmle` R package:

1. **Estimate the outcome regression** Q(A, W) = E[Y | A, W]
2. **Estimate the propensity score** g(W) = P(A = 1 | W)
3. **Compute the clever covariate** H(A, W) = A/g(W) - (1-A)/(1-g(W))
4. **Fluctuate** Q along the clever covariate direction (targeting step)
5. **Substitute** to get the ATE estimate, confidence intervals, and p-value

The tutorial uses simple main-effects models for both Q and g, but notes these can be replaced with Super Learner ensembles.

### Super Learner Libraries (from TL-SAP)

The Statistical Analysis Plan (`reports/TL-SAP.qmd`) specifies ensemble machine learning via Super Learner for both Q and g models:

- **GLM** (generalized linear model)
- **GAM** (generalized additive model)
- **Random Forest**
- **XGBoost**
- **HAL** (Highly Adaptive Lasso, optional)
- **SL.mean** (baseline fallback)

Cross-validated risk metric: Bernoulli log-likelihood.

### Truncation and Diagnostics (from TL-SAP)

- g-value truncation at [0.01, 0.99]
- Effective sample size under weights
- Positivity violation checks
- Influence-curve based standard errors and 95% CIs

### Current Implementation Status

- `R/tmle/tmle_pipeline.R`: **placeholder** — the full TMLE pipeline for the realistic DGP has not yet been implemented
- `simulation/run_simulations.R`: **placeholder** — the simulation driver that would run TMLE across multiple seeds/scenarios has not yet been implemented
- `analysis/analysis_plan.qmd`: **placeholder** — analysis plan integration
- `tests/test_sim.R`: **placeholder** — unit and integration tests
- The tutorial demonstrates a working TMLE on a simple simulated dataset (n = 1000, 2 covariates, continuous outcome) as a proof of concept

---

## Clean-Room Design Framework

The analysis follows a clean-room study design with separated team responsibilities:

| Team | Role |
|---|---|
| **Data Providers** | Deliver de-identified data extracts or synthetic data; document fields and privacy rules |
| **Simulation & QA** | Generate realistic synthetic data; validate distributional assumptions; supply quality gates |
| **Methods** | Choose estimators (TMLE); set model specifications; review causal assumptions |
| **Platform** | Maintain the clean-room runtime; manage packages; ensure reproducibility |
| **Privacy & Compliance** | Verify outputs comply with disclosure controls |
| **Reporting** | Convert outputs into stakeholder-friendly narratives and visualizations |

### Release Checklist

Before results leave the clean room:
1. Methods confirms model choices and sensitivity analyses are documented
2. Simulation & QA verifies reproducibility and design targets are met
3. Privacy & Compliance applies disclosure controls (rounding, minimum cell sizes)
4. Platform packages scripts, session info, and seeds into version control
5. Reporting prepares final summaries with any approved aggregation

---

## Planned Simulation Evaluation Design (from TL-SAP)

The TL-SAP (`reports/TL-SAP.qmd`) specifies a comprehensive simulation evaluation framework:

### DGP Scenarios

1. **Unconfounded scenario** — benchmark where treatment is randomized
2. **Baseline confounding only** — confounding through measured covariates, no censoring complications
3. **Time-varying confounding + informative censoring** — censoring depends on covariates and treatment
4. **Treatment switching scenario** — informative switching tied to treatment and CKD status
5. **Positivity-stress scenario** — rare SOF or rare non-SOF strata to test near-violations of positivity

### Simulation Parameters

- Default sample size: n = 50,000 per replicate
- 500 simulation replicates per scenario
- Fixed seeds for reproducibility

### Evaluation Metrics

| Metric | Definition |
|---|---|
| Bias | estimated ATE minus true ATE |
| RMSE | root mean squared error across replicates |
| Empirical SE vs IC-based SE | agreement between simulation variance and influence-curve SE estimates |
| 95% CI Coverage | proportion of replicates whose CI contains the true ATE |

### Comparator Estimators

- **TMLE** (primary)
- **IPTW** (inverse probability of treatment weighting)
- **Unadjusted** (naive comparison)
- **Cox regression** (standard pharmacoepi approach)

### Planned Reporting Outputs

For each scenario:
- Risk curves under SOF and non-SOF
- Risk difference (RD) and risk ratio (RR) at each time point
- Convergence and diagnostic summaries
- TMLE vs comparator estimator performance tables

---

## Key Causal Assumptions

For TMLE to identify the ATE, the following assumptions are required:

1. **Consistency**: the observed outcome under treatment `a` equals the potential outcome Y(a)
2. **No unmeasured confounding (exchangeability)**: Y(a) ⊥ A | W — treatment is as-if random conditional on measured covariates
3. **Positivity**: 0 < P(A = 1 | W) < 1 for all covariate strata in the population (enforced in the DGP by clipping propensity scores to [0.05, 0.95])
4. **No interference**: one patient's treatment does not affect another's outcome

The DGP is constructed so that assumptions 1–3 hold by design in the simple scenario. In the complex scenario, dependent censoring and non-proportional hazards introduce additional challenges that require appropriate handling (e.g., inverse probability of censoring weights, flexible outcome models).

---

## Challenges Built into the Simulation

The DGP deliberately introduces several realistic analytic challenges:

| Challenge | DGP Toggle | Why It Matters |
|---|---|---|
| Non-proportional hazards | `np_hazard = TRUE` | Treatment effect changes direction over time; standard Cox models will be biased |
| Dependent censoring | `dep_censor = TRUE` | Censoring depends on covariates and treatment; naive Kaplan-Meier and unadjusted analyses are biased |
| Informative treatment switching | `switch_on = TRUE` | Switching depends on treatment and CKD; creates confounding in the censoring mechanism |
| Non-linear confounding | `complexity = TRUE` | Quadratic, trigonometric, and interaction terms in both propensity and outcome models; linear regression / logistic regression will be misspecified |
| High-dimensional covariate space | Always on | 30+ covariates with varying prevalences; variable selection and regularization matter |
| Positivity near-violations | Always on | Some covariate strata may have very low treatment probability, even with clipping |

---

## Dependencies

### Currently Used
- **R packages**: `tidyverse`, `tmle`, `ggplot2`, `dplyr`, `here`
- **Optional**: `missForest` (for imputation), `reticulate` (for Python interop in the frengression example)
- **Python** (optional): `frengression` package (accessed via `reticulate`)

### Planned (from TL-SAP)
- `ltmle` — longitudinal TMLE implementation
- `lmtp` — alternative TMLE/g-estimation for modified treatment policies
- `SuperLearner` — ensemble machine learning for Q and g models
- `data.table` — fast data manipulation for simulation loops
- `simcausal` — alternative causal simulation framework (optional)

---

## What Remains to Be Implemented

1. **Full TMLE pipeline** (`R/tmle/tmle_pipeline.R`): Apply TMLE to the realistic HCV-AKI simulated data, including Super Learner for flexible estimation, inverse probability of censoring weights for dependent censoring, and handling of the time-to-event outcome structure
2. **Simulation driver** (`simulation/run_simulations.R`): Run the analysis across multiple simulation replications to assess estimator performance (bias, coverage, MSE)
3. **Comparison estimators**: Benchmark TMLE against naive approaches (unadjusted, regression-adjusted, IPW) to demonstrate the value of doubly-robust estimation under model misspecification
4. **Diagnostics**: Propensity score overlap checks, influence curve diagnostics, sensitivity analyses for unmeasured confounding
5. **Reporting outputs**: Formatted tables and figures summarizing simulation results across scenarios

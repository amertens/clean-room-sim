# HCV-AKI Clean-Room Simulation: Targeted Learning in Pharmacoepidemiology

A reproducible, end-to-end R package-style repository demonstrating
**Targeted Maximum Likelihood Estimation (TMLE)** applied to a simulated
pharmacoepidemiology case study under a **clean-room staged analysis
framework**.

## Purpose

This project implements:

1. **A clean-room (staged) workflow** for applying TMLE to an HCV-AKI
   pharmacoepidemiology case study, with explicit checkpoints, decision
   logging, and outcome-blind governance at each stage.

2. **A two-phase simulation study** with an **outcome-blind selection phase**
   followed by an **unblinded analysis**, comparing TMLE, cross-fitted TMLE,
   AIPW, IPTW, G-computation, and Cox PH across five scenarios with
   increasing complexity (nonlinear confounding, informative censoring,
   treatment switching, non-proportional hazards).

3. **A draft manuscript** (`reports/clean_room_tmle_manuscript.qmd`) that
   documents the methodology, simulation results, and findings.

## Clean-Room Stages

| Stage | Purpose | Allowed Outputs | Checkpoint |
|-------|---------|-----------------|------------|
| **1. Cohort Build** | Feasibility assessment | Cohort size, missingness, overall event rate (NO treatment comparisons) | Min N, min events, max missingness |
| **2. Design Checks** | PS overlap and balance | PS overlap, SMDs, ESS, weight distributions (NO outcome models) | Max SMD, min ESS, overlap |
| **3. Estimation** | TMLE + comparators | Risk differences at 90/180 days, diagnostics, Cox HR (secondary) | Requires Stage 2 PASS |
| **4. Reporting** | Locked report | Compiled report with all artifacts | Requires Stage 3 complete |

All decisions are logged to `outputs/decision_log.csv` with timestamps,
justifications, and outcome-blindness confirmation.

## Repository Structure

```
clean-room-sim/
├── DESCRIPTION              # R package metadata
├── config/
│   └── default.yml          # All analysis parameters
├── R/
│   ├── dgp/
│   │   └── DGP.R            # Data-generating process
│   ├── tmle/
│   │   └── tmle_pipeline.R  # TMLE implementations
│   ├── estimators/
│   │   ├── cox_ph_estimator.R
│   │   ├── iptw_survival.R
│   │   ├── aipw_survival.R  # Augmented IPW (doubly-robust)
│   │   └── gcomp_risk.R
│   ├── staging/
│   │   ├── stage1_cohort.R  # Stage 1: Cohort build
│   │   ├── stage2_design.R  # Stage 2: Design checks
│   │   ├── stage3_estimation.R  # Stage 3: Estimation
│   │   ├── stage4_reporting.R   # Stage 4: Reporting
│   │   ├── checkpoints.R       # Checkpoint read/write
│   │   └── decision_log.R      # Decision logging
│   └── utils/
│       ├── config.R         # Config loader
│       └── helpers.R        # Shared utilities
├── simulation/
│   └── run_simulations.R    # Full simulation driver
├── analysis/
│   ├── case_study_clean_room.qmd  # End-to-end case study
│   └── simulation_study.qmd      # Simulation results
├── tests/
│   └── testthat/            # Unit tests
├── outputs/                 # Generated outputs (gitignored)
├── run_stage1.R             # Stage 1 runner
├── run_stage2.R             # Stage 2 runner
└── run_stage3_report.R      # Stage 3+4 runner
```

## Quickstart

### Prerequisites

Install required R packages:

```r
install.packages(c("yaml", "jsonlite", "ggplot2", "dplyr", "data.table",
                   "survival", "SuperLearner", "tmle", "gam", "glmnet",
                   "ranger"))
```

### Run the Staged Pipeline

```bash
# Stage 1: Build cohort and check feasibility
Rscript run_stage1.R

# Stage 2: Design adequacy checks (PS overlap, balance)
Rscript run_stage2.R

# Stage 3 + 4: Estimation and locked report
Rscript run_stage3_report.R
```

### Run the Simulation Study

```bash
Rscript simulation/run_simulations.R
```

This runs a two-phase simulation:
1. **Phase 1 (Outcome-blind)**: Assesses PS overlap, ESS, convergence, and
   runtime for each estimator *without* using outcome data.
2. **Phase 2 (Unblinded)**: Runs all six estimators (TMLE, TMLE-CF, AIPW,
   IPTW, G-computation, Cox PH) across five scenarios with 200 replicates.
Results are saved to `outputs/simulation/`.

### Run Tests

```bash
Rscript -e "testthat::test_dir('tests/testthat')"
```

### Render Vignettes and Manuscript

```bash
quarto render analysis/case_study_clean_room.qmd
quarto render analysis/simulation_study.qmd
quarto render reports/clean_room_tmle_manuscript.qmd
```

## Configuration

All parameters are controlled via `config/default.yml`:

- **DGP parameters**: sample size, hazard rates, follow-up
- **Scenario toggles**: np_hazard, dep_censor, complexity, switch_on
- **Simulation settings**: replicates, time points, parallelization
- **TMLE settings**: truncation bounds, SuperLearner libraries
- **Stage gate thresholds**: min N, max SMD, min ESS

Override by editing the config or passing a custom path:

```bash
Rscript run_stage1.R config/my_custom.yml nonlinear
```

## Data Policy

No large datasets are committed. All data is generated deterministically
from seeds specified in the config. The DGP (`R/dgp/DGP.R`) produces
cohorts on-the-fly. Only minimal artifacts (JSON checkpoints, CSV
summaries) are saved to `outputs/`.

## Estimand

**Primary (PH-free)**: Risk difference and risk ratio at clinically
meaningful time points (90, 180 days):

$$\psi(t) = E[I(T^{a=1} \le t)] - E[I(T^{a=0} \le t)]$$

**Secondary (assumes PH)**: Hazard ratio from Cox regression, with
Schoenfeld test for PH diagnostic.

## Simulation Scenarios

| Scenario | Nonlinear | Dep. Censoring | Switching | Non-PH |
|----------|-----------|----------------|-----------|--------|
| simple | - | - | - | - |
| nonlinear | Yes | - | - | - |
| dep_censor | Yes | Yes | - | - |
| switching | Yes | - | Yes | - |
| np_hazard | Yes | - | - | Yes |

## Reproducibility

- All seeds are set via `set.seed()` with `L'Ecuyer-CMRG` for
  parallel-safe reproducibility.
- Session info is captured in all outputs.
- Config is stored alongside results for auditability.

## License

MIT

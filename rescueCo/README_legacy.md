# Clean-Room TMLE Workflow: RescueCo Kenya Trauma Registry

## Purpose

This workflow implements a **clean-room causal inference analysis** of the RescueCo Kenya Trauma Registry, comparing Rescue.Co EMS transport with non-EMS transport on trauma outcomes. It demonstrates how clean-room design principles (staged, outcome-blind decision-making) can be combined with modern targeted learning methods.

The workflow supports two outcome frameworks simultaneously:
- **Binary:** Good functional outcome (GOSE > 4)
- **Time-to-event:** Survival analysis through 30, 90, and 180 days

## Difference from Main Analysis

This workflow is **completely separate** from the existing `analysis/` pipeline:
- No scripts in `analysis/` are modified
- All outputs go to `rescueCo/results/` and `rescueCo/reports/`
- Uses its own configuration (`rescueCo/config/clean_room_config.yml`)
- Reads the same input data but processes it independently
- Enforces strict staged clean-room protocol (no outcome peeking)

## Key Design Principles

- **Pre-specification:** All variable selections, SuperLearner libraries, truncation rules, and gating criteria are declared in `clean_room_config.yml` before analysis begins
- **Decision logging:** Every analytic choice is recorded with timestamp, type, and rationale in a structured decision log (CSV)
- **Double robustness:** TMLE provides consistent estimates if either the outcome or treatment model is correctly specified
- **Semiparametric efficiency:** Influence-function-based inference achieves the efficiency bound among regular estimators
- **Negative controls:** Pre-treatment covariates are tested for zero effect to detect residual confounding before outcomes are examined
- **GO/FLAG/STOP gating:** Formal decision criteria determine whether to proceed, investigate further, or halt based on diagnostics

## Estimators Compared

| ID | Estimator | Outcome | Cohort |
|----|-----------|---------|--------|
| A | PS-matched logistic regression | Binary GOSE | Matched |
| B | PS-matched TMLE | Binary GOSE | Matched |
| C | Outcome-blind simulation | Synthetic | Both |
| D | Full-cohort TMLE | Binary GOSE | Full |
| E | Cox PH regression | Survival | Full |
| F | Cox PH in matched cohort | Survival | Matched |
| G | Discrete-time survival TMLE | Survival | Full |
| H | Discrete-time survival TMLE | Survival | Matched |

## How to Run

### Prerequisites

R packages required:
```r
install.packages(c("yaml", "SuperLearner", "tmle", "MatchIt",
                   "survival", "ggplot2", "randomForest", "glmnet",
                   "gam", "quarto"))
```

### Full Pipeline

From the **project root** directory:
```r
source("rescueCo/scripts/00_run_all.R")
```

### Individual Stages

Run stages sequentially (each depends on prior stage outputs):
```r
source("rescueCo/scripts/01_stage1_build_cohort.R")
source("rescueCo/scripts/02_stage2_design_diagnostics.R")
source("rescueCo/scripts/02b_stage2b_negative_controls.R")
source("rescueCo/scripts/03_stage3_outcome_blind_simulation.R")
source("rescueCo/scripts/04_stage4_binary_outcome_analysis.R")
source("rescueCo/scripts/05_stage5_survival_outcome_analysis.R")
source("rescueCo/scripts/06_stage6_render_report.R")
```

## Output Locations

| Type | Location |
|------|----------|
| Intermediate RDS | `rescueCo/results/*.rds` |
| Summary CSVs | `rescueCo/results/*.csv` |
| Diagnostic plots | `rescueCo/results/*.png` |
| Pipeline log | `rescueCo/logs/pipeline.log` |
| HTML report | `rescueCo/reports/rescueco_clean_room_case_study.html` |

## Directory Structure

```
rescueCo/
  config/
    clean_room_config.yml       # All parameters and variable mappings
  R/
    utils.R                     # Utilities, logging, imputation, SMD
    ps_matching.R               # Propensity score and matching
    tmle_binary.R               # Binary outcome TMLE functions
    tmle_survival.R             # Discrete-time survival TMLE
    simulation.R                # Outcome-blind simulation study
    negative_controls.R         # NC TMLE, plasmode sim, GO/FLAG/STOP
    plotting.R                  # All visualization functions
  scripts/
    00_run_all.R                # Master pipeline runner
    01_stage1_build_cohort.R    # Load data, build cohort
    02_stage2_design_diagnostics.R  # PS estimation, matching, balance, IPTW diagnostics
    02b_stage2b_negative_controls.R # Negative control TMLE + plasmode + gating
    03_stage3_outcome_blind_simulation.R  # Simulation study
    04_stage4_binary_outcome_analysis.R   # GOSE analysis
    05_stage5_survival_outcome_analysis.R # Survival analysis
    06_stage6_render_report.R   # Compile report and decision log
  reports/
    rescueco_clean_room_case_study.qmd  # Quarto report template
  results/                      # Generated outputs (RDS, CSV, PNG)
  logs/                         # Pipeline execution logs
```

## Staged Clean-Room Protocol

| Stage | What happens | Outcome data used? |
|-------|-------------|-------------------|
| 1 | Cohort build, covariate prep, missingness | No |
| 1b | Exploratory data analysis (Quarto report) | No |
| 2 | PS estimation, matching, balance, IPTW diagnostics | No |
| 2b | Negative control TMLE + plasmode simulation + GO/FLAG/STOP | No |
| 3 | Simulation with synthetic outcomes | No (synthetic only) |
| 4 | Binary GOSE analysis + IC/clever-covariate diagnostics | Yes |
| 5 | Survival analysis + risk curves with CI bands | Yes |
| 6 | Report compilation, decision log | Summary only |

## Stage 2b: Negative Control Analysis

Stage 2b applies TMLE to pre-treatment covariates that ambulance transport should not causally affect. If a non-zero effect is found, it signals residual confounding or model misspecification.

**Negative control outcomes** (configured in `clean_room_config.yml`):
- `chronic_hypertension`, `chronic_diabetes_insulin`, `chronic_hiv_art` (pre-existing conditions)
- `household_urban`, `fuel_wood` (SES indicators)

**Plasmode simulation** resamples treatment from the fitted PS and generates outcomes under a known zero ATE. This compares the bias, RMSE, and coverage of TMLE vs. PS matching vs. g-computation.

**GO/FLAG/STOP gating:**
- **GO:** All NC outcomes pass (CI includes 0) and plasmode coverage >= 0.90
- **FLAG:** 1-2 NC failures or coverage 0.80-0.90
- **STOP:** 3+ NC failures or coverage < 0.80

## Assumptions and User Actions

### GOSE outcome
- Variable: `gose` (integer 1-8) from 6-month follow-up
- Binary threshold: GOSE > 4 = good outcome
- ~35% missing due to loss to follow-up
- If `gose` column name differs, update `binary_outcome.gose_variable` in config

### Survival outcome
- **Constructed as a proxy** from available mortality indicators:
  - In-hospital death: event at day 1
  - Death by 6 months: event at day 90 (midpoint)
  - Alive at follow-up: censored at day 180
  - Lost to follow-up: censored at day 30
- **Users should verify/replace** these proxy survival times if more precise timing is available
- Update `survival_outcome` section in config for different variable names

### To plug in exact variable names
1. Edit `rescueCo/config/clean_room_config.yml`
2. Update the `binary_outcome`, `survival_outcome`, and `treatment` sections
3. Re-run the pipeline from Stage 1

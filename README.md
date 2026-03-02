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

---

## Repository Structure

```
clean-room-sim/
├── DESCRIPTION                          # R package metadata
├── NAMESPACE                            # R package exports
├── LICENSE                              # MIT
├── ANALYSIS.md                          # Detailed analysis description
│
├── config/
│   └── default.yml                      # All analysis parameters (DGP, scenarios,
│                                        #   TMLE, stage gates, output dirs)
│
├── R/                                   # ── Core R source code ──
│   ├── dgp/
│   │   └── DGP.R                        # Data-generating process (generate_hcv_data)
│   ├── tmle/
│   │   └── tmle_pipeline.R              # TMLE + cross-fitted TMLE for survival
│   ├── estimators/
│   │   ├── aipw_survival.R              # Augmented IPW (doubly-robust)
│   │   ├── cox_ph_estimator.R           # Cox PH + Schoenfeld PH test
│   │   ├── gcomp_risk.R                 # G-computation plug-in estimator
│   │   └── iptw_survival.R              # Inverse probability weighted KM
│   ├── staging/
│   │   ├── stage1_cohort.R              # Stage 1: Cohort build & feasibility
│   │   ├── stage2_design.R              # Stage 2: PS overlap, balance, ESS
│   │   ├── stage3_estimation.R          # Stage 3: TMLE + comparator estimation
│   │   ├── stage4_reporting.R           # Stage 4: Locked report compilation
│   │   ├── checkpoints.R               # Read/write/enforce stage checkpoints
│   │   └── decision_log.R              # Append-only decision audit trail
│   └── utils/
│       ├── config.R                     # YAML config loader & validator
│       └── helpers.R                    # Shared utilities (SMD, ESS, seeds, etc.)
│
├── run_stage1.R                         # ── Stage runner: cohort build ──
├── run_stage2.R                         # ── Stage runner: design checks ──
├── run_stage3_report.R                  # ── Stage runner: estimation + report ──
│
├── simulation/
│   └── run_simulations.R                # Two-phase simulation driver
│
├── analysis/                            # ── Interactive analysis vignettes ──
│   ├── case_study_clean_room.qmd        # End-to-end 4-stage case study demo
│   ├── simulation_study.qmd             # Simulation results display & plots
│   └── analysis_plan.qmd               # TL-SAP integration (placeholder)
│
├── reports/                             # ── Manuscript & statistical plan ──
│   ├── clean_room_tmle_manuscript.qmd   # Draft manuscript (methods + results)
│   └── TL-SAP.qmd                       # Targeted Learning Statistical Analysis Plan
│
├── tutorial code/                       # ── Learning materials ──
│   ├── tmle_clean_room_tutorial.qmd     # Simple TMLE tutorial
│   └── frengression_example.Rmd         # Python interop demo
│
├── tests/
│   ├── testthat.R                       # Test runner
│   └── testthat/
│       ├── test-dgp.R                   # DGP unit tests
│       ├── test-tmle.R                  # TMLE estimator tests
│       ├── test-staging.R               # Staging pipeline tests
│       └── test-simulation.R            # Simulation integration tests
│
├── outputs/                             # Generated outputs (gitignored)
│   ├── stage1/                          # Stage 1 artifacts
│   ├── stage2/                          # Stage 2 artifacts + plots
│   ├── stage3/                          # Stage 3 estimates + diagnostics
│   ├── report/                          # Locked final report (all artifacts)
│   ├── simulation/                      # Simulation CSV results
│   ├── checkpoint_1.json               # Stage 1 gate (PASS/FAIL)
│   ├── checkpoint_2.json               # Stage 2 gate (PASS/FAIL)
│   └── decision_log.csv                # Full audit trail
│
└── .github/workflows/
    └── check.yml                        # CI/CD (R CMD check + tests)
```

---

## How to Run (Execution Order)

### Prerequisites

Install required R packages:

```r
install.packages(c("yaml", "jsonlite", "ggplot2", "dplyr", "data.table",
                   "survival", "SuperLearner", "tmle", "gam", "glmnet",
                   "ranger"))
```

### A. Staged Clean-Room Pipeline

Run these **in order** — each stage gates on the previous checkpoint:

```bash
# Stage 1: Build cohort, assess feasibility (outcome-blind)
Rscript run_stage1.R                          # default: config/default.yml, "simple"
Rscript run_stage1.R config/default.yml nonlinear  # or specify scenario

# Stage 2: PS overlap, covariate balance, ESS (outcome-blind)
#   Requires: checkpoint_1.json = PASS
Rscript run_stage2.R

# Stage 3 + 4: Estimation (TMLE + comparators) and locked report
#   Requires: checkpoint_2.json = PASS
Rscript run_stage3_report.R
```

**What each stage produces** — see [Outputs Reference](#outputs-reference) below.

### B. Simulation Study

```bash
Rscript simulation/run_simulations.R
```

This is self-contained and independent of the staged pipeline. It runs a
two-phase simulation:

1. **Phase 1 (Outcome-blind)**: Assesses PS overlap, ESS, convergence, and
   runtime for each estimator *without* using outcome data.
2. **Phase 2 (Unblinded)**: Computes ground-truth risks, then runs all six
   estimators (TMLE, TMLE-CF, AIPW, IPTW, G-computation, Cox PH) across
   five scenarios.

Set `QUICK_TEST <- TRUE` at the top of the script for a fast smoke test
(10 reps, N=500) vs production settings (200 reps, N=5000).

### C. Run Tests

```bash
Rscript -e "testthat::test_dir('tests/testthat')"
```

### D. Render Reports

After running the simulation and/or staged pipeline:

```bash
# Case study vignette (runs the 4-stage pipeline interactively)
quarto render analysis/case_study_clean_room.qmd

# Simulation results (reads saved CSVs from outputs/simulation/)
quarto render analysis/simulation_study.qmd

# Draft manuscript (reads saved CSVs from outputs/simulation/)
quarto render reports/clean_room_tmle_manuscript.qmd

# Statistical analysis plan (standalone, no data needed)
quarto render reports/TL-SAP.qmd
```

---

## Outputs Reference

All outputs are written to `outputs/` (gitignored). The directory is created
automatically by each script.

### Staged Pipeline Outputs

| Stage | File | Format | Description |
|-------|------|--------|-------------|
| **1** | `outputs/stage1/stage1_report.json` | JSON | Cohort size, event rate, missingness, covariates |
| **1** | `outputs/checkpoint_1.json` | JSON | PASS/FAIL with criteria (min N, treated frac, missingness) |
| **2** | `outputs/stage2/stage2_diagnostics.json` | JSON | PS summary, weighted SMD, ESS, overlap, truncation |
| **2** | `outputs/stage2/ess_table.csv` | CSV | Effective sample size by treatment group |
| **2** | `outputs/stage2/ps_overlap.png` | PNG | PS density plot (treated vs control) |
| **2** | `outputs/stage2/smd_love_plot.png` | PNG | Covariate balance Love plot |
| **2** | `outputs/stage2/weight_distribution.png` | PNG | IPW weight histogram |
| **2** | `outputs/checkpoint_2.json` | JSON | PASS/FAIL with criteria (max SMD, min ESS, overlap) |
| **3** | `outputs/stage3/stage3_estimates.csv` | CSV | Risk estimates: all methods x time points (RD, RR, SE, CI) |
| **3** | `outputs/stage3/cox_hr_secondary.csv` | CSV | Cox PH hazard ratio, CI, PH test p-value |
| **3** | `outputs/stage3/tmle_diagnostics_t{t}.json` | JSON | TMLE clever covariate, IC, positivity, eps per time point |
| **4** | `outputs/report/` | mixed | Locked report: copies all stage artifacts + `config_used.yml`, `session_info.txt`, `report_metadata.json` |
| all | `outputs/decision_log.csv` | CSV | Append-only audit trail (meeting ID, timestamp, stage, decision, outcome-blind flag) |

### Simulation Outputs

| Phase | File | Format | Description |
|-------|------|--------|-------------|
| **1** | `outputs/simulation/phase1_outcome_blind.csv` | CSV | Per-replicate PS diagnostics (convergence, ESS, SMD, overlap) |
| **1** | `outputs/simulation/phase1_summary.csv` | CSV | Aggregated Phase 1 diagnostics by scenario |
| **2** | `outputs/simulation/true_values.csv` | CSV | Monte Carlo ground-truth risks per scenario x time point |
| **2** | `outputs/simulation/simulation_results.csv` | CSV | Full results: all reps x estimators x time points (estimate, SE, CI, bias, coverage, runtime) |
| **2** | `outputs/simulation/simulation_summary.csv` | CSV | Aggregated performance: bias, RMSE, empirical SD, mean SE, coverage, runtime per method |
| | `outputs/simulation/session_info.txt` | TXT | R session and package versions |

### Reports and Vignettes

| Document | Reads From | Renders To | Description |
|----------|------------|------------|-------------|
| `analysis/case_study_clean_room.qmd` | Runs pipeline live (sources `R/`) | HTML | Interactive walkthrough of all 4 clean-room stages on one dataset |
| `analysis/simulation_study.qmd` | `outputs/simulation/*.csv` | HTML | Simulation results: performance tables, bias/RMSE/coverage plots, SE calibration, runtime comparison, key findings |
| `reports/clean_room_tmle_manuscript.qmd` | `outputs/simulation/*.csv` | HTML/PDF | Draft manuscript with Introduction, Methods, Results (Phase 1 + Phase 2 tables and figures), Discussion, Conclusion |
| `reports/TL-SAP.qmd` | None (standalone) | PDF | Targeted Learning Statistical Analysis Plan: estimand, DGP, estimation strategy, evaluation metrics |

---

## Clean-Room Stages

| Stage | Purpose | Allowed Outputs | Checkpoint |
|-------|---------|-----------------|------------|
| **1. Cohort Build** | Feasibility assessment | Cohort size, missingness, overall event rate (NO treatment comparisons) | Min N, min events, max missingness |
| **2. Design Checks** | PS overlap and balance | PS overlap, SMDs, ESS, weight distributions (NO outcome models) | Max SMD, min ESS, overlap |
| **3. Estimation** | TMLE + comparators | Risk differences at 90/180 days, diagnostics, Cox HR (secondary) | Requires Stage 2 PASS |
| **4. Reporting** | Locked report | Compiled report with all artifacts | Requires Stage 3 complete |

All decisions are logged to `outputs/decision_log.csv` with timestamps,
justifications, and outcome-blindness confirmation.

---

## Simulation Scenarios

| Scenario | Nonlinear | Dep. Censoring | Switching | Non-PH |
|----------|-----------|----------------|-----------|--------|
| simple | - | - | - | - |
| nonlinear | Yes | - | - | - |
| dep_censor | Yes | Yes | - | - |
| switching | Yes | - | Yes | - |
| np_hazard | Yes | - | - | Yes |

## Estimand

**Primary (PH-free)**: Risk difference and risk ratio at clinically
meaningful time points (90, 180 days):

$$\psi(t) = E[I(T^{a=1} \le t)] - E[I(T^{a=0} \le t)]$$

**Secondary (assumes PH)**: Hazard ratio from Cox regression, with
Schoenfeld test for PH diagnostic.

---

## Configuration

All parameters are controlled via `config/default.yml`:

| Section | Key Parameters | Controls |
|---------|----------------|----------|
| `dgp` | N, p_sof, h0, HR_early, HR_late, tau, max_follow | Data-generating process |
| `scenarios` | simple, nonlinear, dep_censor, switching, np_hazard | Scenario toggle flags |
| `simulation` | n_replicates, time_points, parallel, n_cores | Simulation execution |
| `tmle` | truncation_lower/upper, sl_library_Q, sl_library_g | TMLE algorithm tuning |
| `stage1` | min_N, min_treated_frac, max_missingness | Stage 1 gate thresholds |
| `stage2` | max_smd, min_ess_frac, min_ps_overlap | Stage 2 gate thresholds |
| `output` | base_dir, stage1_dir, ..., simulation_dir | Output directory paths |

Override by editing the config or passing a custom path:

```bash
Rscript run_stage1.R config/my_custom.yml nonlinear
```

## Reproducibility

- All seeds are set via `set.seed()` with `L'Ecuyer-CMRG` for
  parallel-safe reproducibility.
- Session info is captured in all outputs.
- Config is stored alongside results for auditability.

## License

MIT

# ============================================================
# Stage 2b/3: Outcome-Blind Plasmode + Data-Quality Stress Test
# ============================================================
# Modernised: uses the cleanTMLE plasmode API
#   - expand_tmle_candidate_grid() to enumerate candidates
#   - validate_tmle_candidates() to reject ill-specified ones
#   - run_plasmode_feasibility() for baseline RMSE/coverage
#   - run_plasmode_dq_stress() for the DQ stress curve (the package's
#     distinguishing contribution: covariate missingness, treatment
#     misclassification, outcome misclassification, unmeasured confounding)
#   - summarize_dq_degradation() to flag fragile candidates
#   - select_tmle_candidate() / lock_primary_tmle_spec() to lock a spec
#   - gate_check() for a structured GO/FLAG/STOP from plasmode metrics
#
# NO real GOSE or survival outcomes are used at this stage.
# ============================================================

library(SuperLearner)

source(file.path(if (dir.exists("rescueCo")) "rescueCo/R" else "R",
                  "bootstrap.R"))
source("rescueCo/R/utils.R")

cfg <- load_cr_config()
cr_log("=== Stage 2b/3: Outcome-Blind Plasmode + DQ Stress (cleanTMLE) ===")

# Load Stage 2 outputs
dat          <- load_stage_output("stage1_cohort.rds")
W_matrix     <- load_stage_output("stage1_W_matrix.rds")
match_result <- load_stage_output("stage2_match_result.rds")
decisions    <- tryCatch(load_stage_output("stage2_decisions.rds"),
                         error = function(e) load_stage_output("stage1_decisions.rds"))
lock         <- tryCatch(load_stage_output("stage2_lock.rds"),
                         error = function(e) load_stage_output("stage1_lock.rds"))
audit        <- tryCatch(load_stage_output("stage2_audit.rds"),
                         error = function(e) load_stage_output("stage1_audit.rds"))
results_dir  <- cfg$paths$results

# Need an unmasked lock for plasmode (plasmode generates synthetic outcomes
# from a Q0 model fit on the real covariate distribution; it does NOT use
# the real outcome column for selection)
lock_for_plasmode <- tryCatch(
  cleanTMLE::unmask_outcome(lock, load_stage_output("stage1_lock_unmasked.rds")),
  error = function(e) lock)

# ============================================================
# STEP 1: Define candidate TMLE specifications
# ============================================================
cr_log("Defining TMLE candidate specifications...")

# Hand-built grid: vary truncation (positivity-cutoff) and learner library
candidates <- list(
  cleanTMLE::tmle_candidate(
    candidate_id = "glm_t01",
    label        = "GLM only, truncation 0.01",
    g_library    = c("SL.glm"),
    q_library    = c("SL.glm"),
    truncation   = 0.01),
  cleanTMLE::tmle_candidate(
    candidate_id = "glmnet_t01",
    label        = "glmnet only, truncation 0.01",
    g_library    = c("SL.glmnet"),
    q_library    = c("SL.glmnet"),
    truncation   = 0.01),
  cleanTMLE::tmle_candidate(
    candidate_id = "ensemble_t01",
    label        = "GLM + glmnet + Mean, truncation 0.01",
    g_library    = c("SL.glm", "SL.glmnet", "SL.mean"),
    q_library    = c("SL.glm", "SL.glmnet", "SL.mean"),
    truncation   = 0.01),
  cleanTMLE::tmle_candidate(
    candidate_id = "ensemble_t05",
    label        = "GLM + glmnet + Mean, truncation 0.05",
    g_library    = c("SL.glm", "SL.glmnet", "SL.mean"),
    q_library    = c("SL.glm", "SL.glmnet", "SL.mean"),
    truncation   = 0.05))

# Validate
cleanTMLE::validate_tmle_candidates(candidates)
cr_log(paste("Defined", length(candidates), "candidate TMLE specifications"))

# Also accept the package's default expansion as a sanity check
default_grid <- tryCatch(
  cleanTMLE::expand_tmle_candidate_grid(
    g_libraries  = list(c("SL.glm"), c("SL.glm", "SL.glmnet", "SL.mean")),
    Q_libraries  = list(c("SL.glm"), c("SL.glm", "SL.glmnet", "SL.mean")),
    truncations  = c(0.01, 0.05)),
  error = function(e) { cr_log(paste("expand_tmle_candidate_grid failed:", e$message)); NULL })
if (!is.null(default_grid))
  cr_log(paste("Default candidate-grid size:", length(default_grid)))

# ============================================================
# STEP 2: Plasmode feasibility — baseline performance per candidate
# ============================================================
cr_log("Running plasmode feasibility (this may take several minutes)...")

n_reps <- cfg$simulation$n_sims %||% 30L
plas <- tryCatch(
  cleanTMLE::run_plasmode_feasibility(
    lock_for_plasmode,
    tmle_candidates = candidates,
    effect_sizes    = c(0.05, 0.10),
    reps            = as.integer(n_reps),
    verbose         = FALSE),
  error = function(e) { cr_log(paste("run_plasmode_feasibility failed:", e$message)); NULL })

if (!is.null(plas)) {
  # summarize_plasmode_results() returns a `plasmode_results` object whose
  # metrics live in $metrics. write.csv can't coerce it directly.
  cleanTMLE::summarize_plasmode_results(plas)  # prints to console
  if (!is.null(plas$metrics)) {
    write.csv(as.data.frame(plas$metrics),
      file.path(results_dir, "plasmode_summary.csv"),
      row.names = FALSE)
  }
  cr_log("cleanTMLE plasmode feasibility complete.")
} else {
  # Fallback to the local plasmode-simulation helper. Same idea: synthetic
  # outcomes generated from a Q0 model fit on the real (X, A); evaluate
  # bias / RMSE / coverage of each candidate. Output is on the same scale
  # as the cleanTMLE summary so downstream code can consume either.
  cr_log("Falling back to local plasmode simulation (rescueCo/R/simulation.R)...")
  source("rescueCo/R/simulation.R", local = TRUE)
  source("rescueCo/R/utils.R",      local = TRUE)
  match_result <- load_stage_output("stage2_match_result.rds")
  ps_trunc     <- load_stage_output("stage2_ps_truncated.rds")
  W_sim <- sanitize_W(W_matrix, max_cols = NULL, A = dat$A)
  plas_local <- tryCatch(
    run_plasmode_simulation(
      A          = dat$A,
      W          = W_sim,
      ps         = ps_trunc,
      matched_idx = match_result$matched_idx,
      sl_lib     = "SL.glm",
      n_sims     = as.integer(min(n_reps, 30L)),
      true_ate   = 0.05,
      seed       = cfg$seed),
    error = function(e) { cr_log(paste("Local plasmode failed:", e$message)); NULL })
  if (!is.null(plas_local) && !is.null(plas_local$metrics)) {
    plas <- plas_local
    write.csv(plas_local$metrics,
      file.path(results_dir, "plasmode_summary.csv"), row.names = FALSE)
    cr_log("Local plasmode complete (fallback).")
  } else {
    cr_log("Plasmode skipped — pipeline will continue without selection")
  }
}

# ============================================================
# STEP 3: DQ stress test — the marquee feature
# ============================================================
cr_log("Running plasmode DQ stress test (multiple data-quality threats)...")

dq <- tryCatch(
  cleanTMLE::run_plasmode_dq_stress(
    lock_for_plasmode,
    tmle_candidates = candidates,
    effect_sizes    = c(0.05),
    reps            = as.integer(min(n_reps, 5L)),
    data_quality_scenarios = list(
      covariate_missingness = list(fractions = c(0.10, 0.20)),
      treatment_misclass    = list(sensitivity = c(0.95),
                                    specificity = c(0.95)),
      outcome_misclass      = list(sensitivity = c(0.90),
                                    specificity = c(0.95)),
      unmeasured_confounding = list(U_prevalence = 0.20,
                                     U_treatment_OR = c(2.0),
                                     U_outcome_OR   = c(2.0))),
    verbose = FALSE),
  error = function(e) { cr_log(paste("run_plasmode_dq_stress failed:", e$message)); NULL })

if (!is.null(dq)) {
  saveRDS(dq, file.path(results_dir, "plasmode_dq_stress.rds"))
  dq_summary <- tryCatch(cleanTMLE::summarize_dq_degradation(dq),
    error = function(e) NULL)
  if (!is.null(dq_summary)) {
    write.csv(dq_summary, file.path(results_dir, "plasmode_dq_degradation.csv"),
              row.names = FALSE)
    cr_log(paste("DQ stress: scored",
                 length(unique(dq_summary$scenario)), "scenarios across",
                 length(unique(dq_summary$candidate_id)), "candidates"))
  }
}

# ============================================================
# STEP 4: Select primary TMLE candidate (rule-based)
# ============================================================
cr_log("Selecting primary TMLE candidate (rule: min_max_rmse — manuscript primary)...")
# `min_max_rmse` (added in cleanTMLE 0.1.1) selects the candidate with the
# minimum WORST-CASE RMSE across the DQ-stress scenarios, *not* the candidate
# with minimum baseline RMSE. This is the rule the methods paper recommends:
# it favours candidates that stay robust across the plausible data-quality
# threats rather than the one that simply wins under perfect conditions.
selected <- if (!is.null(plas)) {
  tryCatch({
    if (!is.null(dq))
      cleanTMLE::select_tmle_candidate(plas, rule = "min_max_rmse",
                                         dq_results = dq)
    else
      cleanTMLE::select_tmle_candidate(plas, rule = "min_rmse")
  }, error = function(e) {
    cr_log(paste("min_max_rmse failed:", e$message,
                  " — falling back to min_rmse"))
    tryCatch(cleanTMLE::select_tmle_candidate(plas, rule = "min_rmse"),
             error = function(e2) NULL)
  })
} else NULL

if (!is.null(selected)) {
  lock <- tryCatch(cleanTMLE::lock_primary_tmle_spec(lock, selected),
                   error = function(e) { cr_log(paste("lock_primary_tmle_spec failed:", e$message)); lock })
  cr_log(paste("Locked primary TMLE spec:", selected$candidate_id,
               " (truncation =", selected$truncation, ")"))
  decisions <- log_decision(decisions, "stage2b",
    paste("Locked TMLE spec:", selected$candidate_id),
    "Minimum-RMSE rule on plasmode feasibility (pre-outcome)",
    type = "pre-specified")
}

# ============================================================
# STEP 5: GATE CHECK — composite GO/FLAG/STOP from plasmode metrics
# ============================================================
gate <- tryCatch(
  cleanTMLE::gate_check(plas,
    rmse_threshold     = 0.05,
    coverage_threshold = 0.85),
  error = function(e) { cr_log(paste("gate_check failed:", e$message)); NULL })
if (!is.null(gate) && inherits(gate, "cleantmle_checkpoint")) {
  audit <- cleanTMLE::record_checkpoint(audit, gate)
  cr_log(paste("Gate decision (Stage 2b plasmode):", gate$decision))
} else if (!is.null(gate)) {
  cr_log(paste("Gate result returned but not a cleantmle_checkpoint; skipping audit append."))
}

audit <- cleanTMLE::record_stage(audit, "Stage 2b",
  if (!is.null(selected)) paste("Locked TMLE candidate:", selected$candidate_id) else "Plasmode complete")

# ============================================================
# STEP 6: Save Stage 2b outputs
# ============================================================
save_stage_output(plas,        "stage2b_plasmode.rds")
save_stage_output(dq,          "stage2b_dq_stress.rds")
save_stage_output(selected,    "stage2b_selected_candidate.rds")
save_stage_output(candidates,  "stage2b_candidates.rds")
save_stage_output(decisions,   "stage3_decisions.rds")
save_stage_output(lock,        "stage3_lock.rds")
save_stage_output(audit,       "stage3_audit.rds")

# Backwards-compat names
save_stage_output(list(plasmode = plas, dq_stress = dq,
                        selected_candidate = selected),
                  "stage3_simulation_results.rds")

cr_log("Stage 2b/3 complete. NO real outcomes were used.")

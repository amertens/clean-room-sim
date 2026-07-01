# ============================================================
# rescueCo case study: FULL data-quality stress run (reported result)
# ============================================================
# Regenerates the case-study DQ stress artifacts at the reported replicate
# count, now including the near-positivity threat (consistent with the main
# simulation's five-threat sweep). Candidates use SL.glmnet.bounded rather
# than SL.glmnet so the regularised learners cannot run away on the
# near-positivity designs.
#
# Usage:  Rscript rescueCo/scripts/rescueco_dq_full.R [dq_reps] [feas_reps]
#   defaults: dq_reps = 200, feas_reps = 50   (pass "3 3" for a smoke test)
#
# Reuses the staged Stage 1/2 outputs already on disk; produces no real
# outcome estimates. Outputs:
#   rescueCo/results/plasmode_dq_stress.rds
#   rescueCo/results/plasmode_dq_degradation.csv
#   rescueCo/results/stage2b_dq_stress.rds   (back-compat name)
# ============================================================

library(SuperLearner)

args     <- commandArgs(trailingOnly = TRUE)
dq_reps  <- if (length(args) >= 1) as.integer(args[[1]]) else 200L
feas_reps <- if (length(args) >= 2) as.integer(args[[2]]) else 50L

source(file.path(if (dir.exists("rescueCo")) "rescueCo/R" else "R", "bootstrap.R"))
source("rescueCo/R/utils.R")

cfg <- load_cr_config()
cr_log(sprintf("=== rescueCo FULL DQ stress (dq_reps=%d, feas_reps=%d) ===",
               dq_reps, feas_reps))

lock        <- tryCatch(load_stage_output("stage2_lock.rds"),
                        error = function(e) load_stage_output("stage1_lock.rds"))
results_dir <- cfg$paths$results

lock_for_plasmode <- tryCatch(
  cleanTMLE::unmask_outcome(lock, load_stage_output("stage1_lock_unmasked.rds"),
                            allow_unauthorized = TRUE),
  error = function(e) lock)

# Candidate grid: same truncation/library design as the staged pipeline, but
# SL.glmnet -> SL.glmnet.bounded so near-positivity cannot trigger runaway.
candidates <- list(
  cleanTMLE::tmle_candidate(
    candidate_id = "glm_t01", label = "GLM only, truncation 0.01",
    g_library = c("SL.glm"), q_library = c("SL.glm"), truncation = 0.01),
  cleanTMLE::tmle_candidate(
    candidate_id = "glmnet_t01", label = "glmnet(bounded) only, truncation 0.01",
    g_library = c("SL.glmnet.bounded"), q_library = c("SL.glmnet.bounded"),
    truncation = 0.01),
  cleanTMLE::tmle_candidate(
    candidate_id = "ensemble_t01", label = "GLM + glmnet(bounded) + Mean, trunc 0.01",
    g_library = c("SL.glm", "SL.glmnet.bounded", "SL.mean"),
    q_library = c("SL.glm", "SL.glmnet.bounded", "SL.mean"), truncation = 0.01),
  cleanTMLE::tmle_candidate(
    candidate_id = "ensemble_t05", label = "GLM + glmnet(bounded) + Mean, trunc 0.05",
    g_library = c("SL.glm", "SL.glmnet.bounded", "SL.mean"),
    q_library = c("SL.glm", "SL.glmnet.bounded", "SL.mean"), truncation = 0.05))
cleanTMLE::validate_tmle_candidates(candidates)

# Baseline feasibility (needed for the minimax selection).
cr_log("Plasmode feasibility...")
plas <- cleanTMLE::run_plasmode_feasibility(
  lock_for_plasmode, tmle_candidates = candidates,
  effect_sizes = c(0.05, 0.10), reps = feas_reps, verbose = FALSE)

# DQ stress: original four threats at their existing levels PLUS near-positivity.
cr_log("Plasmode DQ stress (five threats incl. near-positivity)...")
dq <- cleanTMLE::run_plasmode_dq_stress(
  lock_for_plasmode, tmle_candidates = candidates,
  effect_sizes = c(0.05), reps = dq_reps,
  data_quality_scenarios = list(
    near_positivity            = list(slopes = c(2.0, 3.0, 4.0)),
    covariate_missingness      = list(fractions = c(0.10, 0.20)),
    covariate_missingness_mar  = list(fractions = c(0.10, 0.20), treatment_OR = 3),
    covariate_missingness_mnar = list(fractions = c(0.10, 0.20), strength = 1.5),
    treatment_misclass     = list(sensitivity = c(0.95), specificity = c(0.95)),
    outcome_misclass       = list(sensitivity = c(0.90), specificity = c(0.95)),
    unmeasured_confounding = list(U_prevalence = 0.20,
                                  U_treatment_OR = c(2.0), U_outcome_OR = c(2.0))),
  verbose = TRUE)

saveRDS(dq, file.path(results_dir, "plasmode_dq_stress.rds"))
save_stage_output(dq, "stage2b_dq_stress.rds")

dq_summary <- cleanTMLE::summarize_dq_degradation(dq)
write.csv(dq_summary, file.path(results_dir, "plasmode_dq_degradation.csv"),
          row.names = FALSE)
cr_log(sprintf("DQ stress: scored %d scenarios across %d candidates",
               length(unique(dq_summary$scenario)),
               length(unique(dq_summary$candidate))))

# Minimax selection (manuscript primary rule), reported for prose.
selected <- tryCatch(
  cleanTMLE::select_tmle_candidate(plas, rule = "min_max_rmse", dq_results = dq),
  error = function(e) { cr_log(paste("min_max_rmse failed:", e$message)); NULL })
if (!is.null(selected))
  cr_log(sprintf("Minimax-selected candidate: %s (truncation = %s)",
                 selected$candidate_id, selected$truncation))

cr_log("rescueCo FULL DQ stress complete. NO real outcomes were used.")
cat("RESCUECO_DQ_FULL_DONE\n")

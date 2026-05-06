# ============================================================
# Clean-Room Workflow: Propensity Score & Matching Functions
# ============================================================

#' Estimate propensity scores using SuperLearner
#'
#' @param A Binary treatment vector
#' @param W Covariate matrix
#' @param sl_lib Character vector of SuperLearner libraries
#' @param cv_folds Number of CV folds
#' @param seed Random seed
#' @return List with ps_scores, sl_fit, and diagnostics
estimate_propensity_score <- function(A, W, sl_lib, cv_folds = 5, seed = 42,
                                      cl = NULL) {
  require(SuperLearner)
  set.seed(seed)


  # Remove observations with missing treatment
  complete <- !is.na(A)
  A_cc <- A[complete]
  W_cc <- as.matrix(W[complete, , drop = FALSE])

  # Sanitize: replace NA/NaN/Inf, drop zero-var cols, fix column names
  for (j in seq_len(ncol(W_cc))) {
    bad <- is.na(W_cc[, j]) | is.nan(W_cc[, j]) | is.infinite(W_cc[, j])
    if (any(bad)) {
      med_j <- median(W_cc[!bad, j], na.rm = TRUE)
      if (is.na(med_j)) med_j <- 0
      W_cc[bad, j] <- med_j
    }
  }
  col_var <- apply(W_cc, 2, var, na.rm = TRUE)
  keep_cols <- !is.na(col_var) & col_var > 1e-6
  W_cc <- W_cc[, keep_cols, drop = FALSE]
  colnames(W_cc) <- make.names(colnames(W_cc), unique = TRUE)

  # Final complete-cases filter: SuperLearner's model.frame will drop rows
  # with any NA, causing length mismatches between Y and X in CV folds
  W_cc <- as.data.frame(W_cc)
  cc2 <- complete.cases(W_cc)
  if (!all(cc2)) {
    message("PS estimation: dropping ", sum(!cc2), " rows with residual NAs in W")
    W_cc <- W_cc[cc2, , drop = FALSE]
    A_cc <- A_cc[cc2]
    # Update the complete-case map
    idx_full <- which(complete)
    complete <- rep(FALSE, length(A))
    complete[idx_full[cc2]] <- TRUE
  }

  # Fit SuperLearner (with optional parallel cluster)
  message("  Fitting PS SuperLearner (", cv_folds, "-fold CV, ",
          length(sl_lib), " learners)...")
  sl_fit <- tryCatch({
    run_sl(Y = A_cc, X = W_cc, family = binomial(),
           SL.library = sl_lib,
           cvControl = list(V = cv_folds), cl = cl)
  }, error = function(e) {
    warning("SuperLearner failed for PS, falling back to GLM: ", e$message)
    SuperLearner(
      Y = A_cc, X = W_cc, family = binomial(),
      SL.library = "SL.glm",
      cvControl = list(V = cv_folds)
    )
  })

  ps <- as.numeric(sl_fit$SL.predict)

  # Full-length PS vector (NA for missing treatment)
  ps_full <- rep(NA_real_, length(A))
  ps_full[complete] <- ps

  list(
    ps_scores = ps_full,
    sl_fit    = sl_fit,
    cv_risk   = sl_fit$cvRisk,
    coef      = sl_fit$coef,
    n_complete = sum(complete),
    n_missing  = sum(!complete)
  )
}

#' Perform 1:1 nearest-neighbor matching on logit PS
#'
#' @param A Treatment vector
#' @param ps Propensity score vector
#' @param caliper_sd_mult Caliper as multiplier of SD(logit PS)
#' @return List with matched indices, match object, diagnostics
perform_ps_matching <- function(A, ps, caliper_sd_mult = 0.2) {
  require(MatchIt)

  # Filter to non-missing

  keep <- !is.na(A) & !is.na(ps) & ps > 0 & ps < 1
  df <- data.frame(A = A[keep], ps = ps[keep], idx = which(keep))

  logit_ps <- log(df$ps / (1 - df$ps))
  caliper <- caliper_sd_mult * sd(logit_ps, na.rm = TRUE)

  # MatchIt

  m_out <- matchit(
    A ~ ps,
    data = df,
    method = "nearest",
    distance = df$ps,
    caliper = caliper,
    ratio = 1,
    replace = FALSE
  )

  md <- match.data(m_out)

  # Map back to original indices

  matched_original_idx <- df$idx[as.integer(rownames(md))]

  list(
    match_object     = m_out,
    match_data       = md,
    matched_idx      = matched_original_idx,
    caliper          = caliper,
    n_treated_matched = sum(md$A == 1),
    n_control_matched = sum(md$A == 0),
    n_treated_total  = sum(df$A == 1),
    frac_treated_matched = sum(md$A == 1) / sum(df$A == 1)
  )
}

#' Generate PS overlap density plot data
ps_overlap_data <- function(ps, A) {
  data.frame(
    ps = ps[!is.na(ps) & !is.na(A)],
    group = ifelse(A[!is.na(ps) & !is.na(A)] == 1, "Treated", "Control")
  )
}

#' Generate Love plot data (SMD before and after matching)
love_plot_data <- function(W, A, matched_idx) {
  smd_pre  <- compute_all_smds(W, A)
  smd_post <- compute_all_smds(W, A, matched_idx)

  merge(
    smd_pre,  smd_post,
    by = "covariate", suffixes = c("_pre", "_post")
  )
}

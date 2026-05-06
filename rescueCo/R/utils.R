# ============================================================
# Clean-Room Workflow: Utility Functions
# ============================================================

#' Locate the rescueCo/ folder relative to the working directory.
#' Allows the bundle to be lifted into any project — scripts work whether
#' invoked from the project root or from inside rescueCo/.
.cr_root <- function() {
  for (cand in c(".", "rescueCo", "..", "../rescueCo")) {
    if (dir.exists(file.path(cand, "config")) &&
        file.exists(file.path(cand, "config", "clean_room_config.yml"))) {
      return(normalizePath(cand, mustWork = TRUE))
    }
  }
  if (file.exists("config/clean_room_config.yml")) return(normalizePath("."))
  if (file.exists("rescueCo/config/clean_room_config.yml"))
    return(normalizePath("rescueCo"))
  stop("Cannot locate rescueCo/ folder. Run from project root or inside rescueCo/.")
}

#' Load clean-room config (auto-locates relative to working directory).
load_cr_config <- function(config_path = NULL) {
  if (!requireNamespace("yaml", quietly = TRUE)) stop("Package 'yaml' required")
  if (is.null(config_path)) {
    cr <- .cr_root()
    # Find the right working layout: when invoked from project root, paths
    # should remain "rescueCo/..."; when invoked from inside clean_room,
    # config paths reference the same files but the working directory is
    # already rescueCo/. Read from whichever exists.
    if (file.exists("rescueCo/config/clean_room_config.yml")) {
      config_path <- "rescueCo/config/clean_room_config.yml"
    } else {
      config_path <- file.path(cr, "config", "clean_room_config.yml")
    }
  }
  cfg <- yaml::read_yaml(config_path)
  # Resolve data paths: if running from inside rescueCo/, strip the leading
  # "rescueCo/" so paths still find the .dta files.
  fix_path <- function(p) {
    if (is.null(p) || !nzchar(p)) return(p)
    if (file.exists(p)) return(p)
    p_strip <- sub("^rescueCo/", "", p)
    if (file.exists(p_strip)) return(p_strip)
    p
  }
  if (!is.null(cfg$data)) {
    for (k in names(cfg$data)) cfg$data[[k]] <- fix_path(cfg$data[[k]])
  }
  cfg
}

#' Log a decision with stage, type, and rationale
log_decision <- function(decision_log, stage, decision, rationale,
                         type = c("pre-specified", "design-stage",
                                  "simulation-selected", "final analysis")) {
  type <- match.arg(type)

  new_row <- data.frame(
    timestamp = Sys.time(),
    stage     = stage,
    type      = type,
    decision  = decision,
    rationale = rationale,
    stringsAsFactors = FALSE
  )

  rbind(decision_log, new_row)
}

#' Initialize empty decision log
init_decision_log <- function() {
  data.frame(
    timestamp = character(0),
    stage     = character(0),
    type      = character(0),
    decision  = character(0),
    rationale = character(0),
    stringsAsFactors = FALSE
  )
}

#' Save intermediate results as RDS
save_stage_output <- function(obj, filename, results_dir = "rescueCo/results") {
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(results_dir, filename)
  saveRDS(obj, path)
  message("Saved: ", path)
  invisible(path)
}

#' Load intermediate results
load_stage_output <- function(filename, results_dir = "rescueCo/results") {
  path <- file.path(results_dir, filename)
  if (!file.exists(path)) stop("Stage output not found: ", path)
  readRDS(path)
}

#' Write a log message to file and console
cr_log <- function(msg, log_dir = "rescueCo/logs") {
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  log_file <- file.path(log_dir, "pipeline.log")
  entry <- paste0("[", Sys.time(), "] ", msg)
  message(entry)
  cat(entry, "\n", file = log_file, append = TRUE)
}

#' Impute missing values in a data frame
#' Numeric: median, Categorical: explicit "Missing" level
impute_missing <- function(df) {
  for (col in names(df)) {
    if (is.numeric(df[[col]])) {
      med <- median(df[[col]], na.rm = TRUE)
      df[[col]][is.na(df[[col]])] <- med
    } else if (is.factor(df[[col]]) || is.character(df[[col]])) {
      df[[col]] <- as.character(df[[col]])
      df[[col]][is.na(df[[col]]) | df[[col]] == ""] <- "Missing"
      df[[col]] <- as.factor(df[[col]])
    }
  }
  df
}

#' Summarize missingness across columns
missingness_summary <- function(df) {
  count_missing <- function(x) {
    n_na <- sum(is.na(x))
    if (is.character(x)) n_na <- n_na + sum(x == "", na.rm = TRUE)
    n_na
  }
  data.frame(
    variable   = names(df),
    n_total    = nrow(df),
    n_missing  = sapply(df, count_missing),
    pct_missing = sapply(df, function(x) round(count_missing(x) / length(x) * 100, 1)),
    stringsAsFactors = FALSE
  )
}

#' Truncate propensity scores
truncate_ps <- function(ps, lower = 0.01, upper = 0.99) {
  pmin(pmax(ps, lower), upper)
}

#' Compute standardized mean difference
compute_smd <- function(x, group, weighted = FALSE, weights = NULL) {
  g1 <- x[group == 1]
  g0 <- x[group == 0]
  if (weighted && !is.null(weights)) {
    w1 <- weights[group == 1]
    w0 <- weights[group == 0]
    m1 <- weighted.mean(g1, w1, na.rm = TRUE)
    m0 <- weighted.mean(g0, w0, na.rm = TRUE)
  } else {
    m1 <- mean(g1, na.rm = TRUE)
    m0 <- mean(g0, na.rm = TRUE)
  }
  s1 <- sd(g1, na.rm = TRUE)
  s0 <- sd(g0, na.rm = TRUE)
  pooled_sd <- sqrt((s1^2 + s0^2) / 2)
  if (pooled_sd == 0) return(0)
  (m1 - m0) / pooled_sd
}

#' Compute SMDs for all covariates
compute_all_smds <- function(W, A, matched_idx = NULL) {
  if (!is.null(matched_idx)) {
    W <- W[matched_idx, , drop = FALSE]
    A <- A[matched_idx]
  }
  smds <- sapply(seq_len(ncol(W)), function(j) {
    compute_smd(W[, j], A)
  })
  data.frame(
    covariate = colnames(W),
    smd = smds,
    stringsAsFactors = FALSE
  )
}

#' Effective sample size under weighting
effective_sample_size <- function(weights) {
  sum(weights)^2 / sum(weights^2)
}

# ============================================================
# IPTW Weight Diagnostics
# ============================================================

#' Summarize IPTW weights by treatment group
make_weight_summary_table <- function(weights, A) {
  summarize_group <- function(w) {
    data.frame(
      n      = length(w),
      min    = min(w),
      Q1     = quantile(w, 0.25, names = FALSE),
      median = median(w),
      Q3     = quantile(w, 0.75, names = FALSE),
      max    = max(w),
      mean   = mean(w),
      SD     = sd(w),
      stringsAsFactors = FALSE
    )
  }
  rbind(
    cbind(group = "Treated (A=1)", summarize_group(weights[A == 1])),
    cbind(group = "Control (A=0)", summarize_group(weights[A == 0]))
  )
}

#' Detect extreme IPTW weights
#' @return list with n_extreme, pct_extreme, max_min_ratio, extreme_idx, flag
detect_extreme_weights <- function(weights, A,
                                   pctl = 0.99,
                                   ratio_threshold = 50) {
  cutoff <- quantile(weights, pctl, names = FALSE)
  extreme_idx <- which(weights > cutoff)
  max_min_ratio <- max(weights) / max(min(weights), 1e-10)

  list(
    n_extreme     = length(extreme_idx),
    pct_extreme   = length(extreme_idx) / length(weights) * 100,
    threshold     = cutoff,
    max_min_ratio = max_min_ratio,
    extreme_idx   = extreme_idx,
    flag          = (length(extreme_idx) / length(weights) > 0.05) ||
                    (max_min_ratio > ratio_threshold)
  )
}

#' Compute IPTW-weighted SMDs for all covariates
compute_weighted_smds <- function(W, A, weights) {
  smds <- sapply(seq_len(ncol(W)), function(j) {
    compute_smd(W[, j], A, weighted = TRUE, weights = weights)
  })
  data.frame(
    covariate = colnames(W),
    smd = smds,
    stringsAsFactors = FALSE
  )
}

# ============================================================
# Progress bar helper
# ============================================================

#' Create a simple text progress bar that works in console and Rscript
#' @param total Total iterations
#' @param prefix Label prefix
#' @return A list with $tick() and $done() methods
cr_progress <- function(total, prefix = "") {
  env <- new.env(parent = emptyenv())
  env$i <- 0L
  env$total <- total
  env$prefix <- prefix
  env$start <- Sys.time()
  env$width <- 30L

  tick <- function() {
    env$i <- env$i + 1L
    pct  <- env$i / env$total
    filled <- round(pct * env$width)
    bar <- paste0(
      "\r", env$prefix, " [",
      strrep("=", filled), strrep(" ", env$width - filled),
      "] ", env$i, "/", env$total,
      " (", round(pct * 100), "%)"
    )
    if (env$i == env$total) {
      elapsed <- as.numeric(difftime(Sys.time(), env$start, units = "secs"))
      bar <- paste0(bar, "  done in ", round(elapsed, 1), "s\n")
    }
    cat(bar)
    flush.console()
  }

  done <- function() {
    if (env$i < env$total) {
      elapsed <- as.numeric(difftime(Sys.time(), env$start, units = "secs"))
      cat("\n", env$prefix, "finished in", round(elapsed, 1), "s\n")
    }
  }

  list(tick = tick, done = done)
}

# ============================================================
# Parallel SuperLearner helper (Windows-safe via PSOCK cluster)
# ============================================================

#' Set up or return a PSOCK cluster for SuperLearner on Windows
#' @param n_cores Number of cores
#' @return cluster object (or NULL if parallelization fails)
setup_sl_cluster <- function(n_cores = 2) {
  if (!requireNamespace("parallel", quietly = TRUE)) return(NULL)
  tryCatch({
    cl <- parallel::makeCluster(n_cores, type = "PSOCK")
    # Export libraries needed by SL learners
    parallel::clusterEvalQ(cl, {
      suppressMessages({
        library(SuperLearner)
        if (requireNamespace("glmnet", quietly = TRUE)) library(glmnet)
        if (requireNamespace("gam", quietly = TRUE)) library(gam)
      })
    })
    cl
  }, error = function(e) {
    message("Parallel cluster setup failed: ", e$message, ". Running sequentially.")
    NULL
  })
}

#' Stop a cluster safely
stop_sl_cluster <- function(cl) {
  if (!is.null(cl)) {
    tryCatch(parallel::stopCluster(cl), error = function(e) NULL)
  }
}

#' Run SuperLearner with optional parallelization
#' @param Y Outcome vector
#' @param X Covariate data.frame
#' @param family family argument
#' @param SL.library learner library
#' @param cvControl cv control list
#' @param cl Optional parallel cluster
#' @return SuperLearner fit
run_sl <- function(Y, X, family, SL.library, cvControl, cl = NULL) {
  if (!is.null(cl)) {
    # mcSuperLearner doesn't work on Windows; use snowSuperLearner with PSOCK
    tryCatch(
      snowSuperLearner(
        Y = Y, X = X, family = family,
        SL.library = SL.library,
        cvControl = cvControl,
        cluster = cl
      ),
      error = function(e) {
        message("snowSuperLearner failed, falling back to sequential: ", e$message)
        SuperLearner(
          Y = Y, X = X, family = family,
          SL.library = SL.library,
          cvControl = cvControl
        )
      }
    )
  } else {
    SuperLearner(
      Y = Y, X = X, family = family,
      SL.library = SL.library,
      cvControl = cvControl
    )
  }
}

#' Sanitize a covariate matrix for use with SuperLearner / tmle
#'
#' Fixes common problems: NA/NaN/Inf values, zero-variance columns,
#' column names that break formula-based learners (SL.gam, SL.glm),
#' and optionally reduces to a manageable number of columns.
#'
#' @param W numeric matrix or data.frame of covariates
#' @param max_cols Maximum columns to retain (NULL = keep all). Columns are
#'   selected by absolute correlation with \code{A} (if supplied) or variance.
#' @param A Optional treatment vector used to rank columns by relevance.
#' @return A clean numeric matrix ready for SuperLearner
sanitize_W <- function(W, max_cols = NULL, A = NULL) {
  W <- as.matrix(W)


  ## 1. Replace Inf / NaN / NA with column medians (or 0)
  for (j in seq_len(ncol(W))) {
    bad <- is.na(W[, j]) | is.nan(W[, j]) | is.infinite(W[, j])
    if (any(bad)) {
      med_j <- median(W[!bad, j], na.rm = TRUE)
      if (is.na(med_j)) med_j <- 0
      W[bad, j] <- med_j
    }
  }

  ## 2. Drop zero / near-zero variance columns
  col_var <- apply(W, 2, var, na.rm = TRUE)
  keep <- !is.na(col_var) & col_var > 1e-6
  W <- W[, keep, drop = FALSE]

  ## 3. Make column names safe for R formulas
  colnames(W) <- make.names(colnames(W), unique = TRUE)

  ## 4. Optionally select top columns by |cor(W_j, A)|

  if (!is.null(max_cols) && ncol(W) > max_cols) {
    if (!is.null(A) && length(A) == nrow(W)) {
      cors <- abs(apply(W, 2, function(x) {
        tryCatch(cor(x, A, use = "complete.obs"), error = function(e) 0)
      }))
      top_idx <- order(cors, decreasing = TRUE)[seq_len(max_cols)]
    } else {
      # fallback: highest-variance columns
      col_var2 <- apply(W, 2, var, na.rm = TRUE)
      top_idx <- order(col_var2, decreasing = TRUE)[seq_len(max_cols)]
    }
    W <- W[, top_idx, drop = FALSE]
  }

  ## 5. Final assertion
  stopifnot(!anyNA(W))
  W
}

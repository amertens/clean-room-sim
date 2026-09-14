# ============================================================
# Multi-arm Stage 3: outcome-blind support simulation
# ============================================================
# Runs simulate_support() on each lock: the generate-treatment plasmode over
# the prespecified confounding-by-modification grid, with per-estimand truths
# and the support map the gate reads. On the primary cohort and the protocol
# contrast it additionally runs the deprecated sample-treatment design on a
# reduced grid, so the Shaw et al. (2025) artifact is shown side by side on
# this study's own data.
#
# Outcome-blind: the locks carry no outcome column; Q0 is the synthetic
# family anchored at a prespecified 3 percent baseline rate (the registry's
# marginal six-month mortality is public knowledge from the protocol; no
# treatment-outcome association is used).
#
# Checkpointed per grid cell via simulate_support(checkpoint_file = ...).
# Env: MULTIARM_SIM_ONLY = comma-separated contrast names (default: all),
#      MULTIARM_SIM_REPS (default 200).
# ============================================================

suppressMessages(library(cleanTMLE))

source(file.path(if (dir.exists("rescueCo")) "rescueCo/R" else "R",
                 "bootstrap.R"))
source("rescueCo/R/utils.R")
cr_log("=== Multi-arm Stage 3: outcome-blind support simulation ===")

OUT <- "rescueCo/results/multiarm"
locks <- readRDS(file.path(OUT, "multiarm_locks.rds"))

REPS <- suppressWarnings(as.integer(Sys.getenv("MULTIARM_SIM_REPS", "200")))
if (is.na(REPS)) REPS <- 200L
only <- trimws(strsplit(Sys.getenv("MULTIARM_SIM_ONLY", ""), ",")[[1]])
run_set <- if (length(only[nzchar(only)])) only[nzchar(only)] else names(locks)

# Prespecified surface family: additive mortality effect of 2 percentage
# points at the centre, baseline rate 3 percent, confounding up to 2 log-odds
# per SD of the propensity direction, modification up to full (effect
# vanishing one SD into the control-typical region).
SURFACE <- cleanTMLE:::support_surfaces(confounding = c(0, 1, 2),
                            modification = c(0, 1),
                            complexity = "linear",
                            effect = 0.02, base_rate = 0.03)

map_rows <- list(); cmp_rows <- list()

for (cn in intersect(names(locks), run_set)) {
  lk <- locks[[cn]]
  lite_file <- file.path(OUT, paste0("design_lite_", cn, ".rds"))
  full_file <- file.path(OUT, paste0("design_", cn, ".rds"))
  if (!file.exists(lite_file) && file.exists(full_file)) {
    st_full <- readRDS(full_file)
    saveRDS(list(ps_raw = st_full$ps$ps_raw %||% st_full$ps$ps,
                 ps = st_full$ps$ps, support = st_full$support,
                 feasibility = st_full$feasibility), lite_file)
    rm(st_full); gc()
  }
  psf <- if (file.exists(lite_file)) {
    lt <- readRDS(lite_file)
    cleanTMLE:::wrap_ps_fit(lk, ps_scores = lt$ps_raw)
  } else NULL
  ck  <- file.path(OUT, paste0("supportsim_", cn, "_ckpt.rds"))
  res_file <- file.path(OUT, paste0("supportsim_", cn, ".rds"))

  if (file.exists(res_file) && !nzchar(Sys.getenv("MULTIARM_FORCE"))) {
    cr_log(paste("[", cn, "] support simulation exists; loading"))
    sim <- readRDS(res_file)
  } else {
    cr_log(paste("[", cn, "] simulate_support:", REPS, "reps, n =",
                 nrow(lk$data)))
    t0 <- Sys.time()
    sim <- simulate_support(lk, ps_fit = psf, surface = SURFACE,
                            reps = REPS, band = c(0.05, 0.95),
                            q0_source = "synthetic",
                            checkpoint_file = ck, verbose = TRUE)
    cr_log(paste("[", cn, "] done in",
                 round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1),
                 "min"))
    saveRDS(sim, res_file)
    p <- plot(sim)
    ggplot2::ggsave(file.path(OUT, paste0("fig_supportmap_", cn, ".png")), p,
                    width = 9, height = 6, dpi = 150)
  }
  map_rows[[cn]] <- cbind(contrast = cn, sim$metrics,
                          stringsAsFactors = FALSE)
  cr_log(paste("[", cn, "] feasible per map:",
               paste(names(sim$feasible)[sim$feasible], collapse = ", ")))

  # The design comparison on the required pair: generate vs sample treatment
  # at a reduced grid, same replicate count.
  if (cn %in% c("PRIMARY", "C1")) {
    cmp_file <- file.path(OUT, paste0("designcmp_", cn, ".rds"))
    if (file.exists(cmp_file) && !nzchar(Sys.getenv("MULTIARM_FORCE"))) {
      cmp <- readRDS(cmp_file)
    } else {
      surf2 <- cleanTMLE:::support_surfaces(confounding = c(0, 2), modification = 0,
                                complexity = "linear",
                                effect = 0.02, base_rate = 0.03)
      cr_log(paste("[", cn, "] design comparison (generate vs sample),",
                   REPS, "reps"))
      gen <- simulate_support(lk, ps_fit = psf, surface = surf2, reps = REPS,
                              checkpoint_file = file.path(
                                OUT, paste0("designcmp_gen_", cn, ".rds")),
                              verbose = FALSE)
      smp <- suppressWarnings(simulate_support(
        lk, ps_fit = psf, surface = surf2, reps = REPS,
        design = "sample_treatment",
        checkpoint_file = file.path(OUT,
                                    paste0("designcmp_smp_", cn, ".rds")),
        verbose = FALSE))
      cmp <- rbind(cbind(design = "generate_treatment", gen$metrics),
                   cbind(design = "sample_treatment", smp$metrics))
      saveRDS(cmp, cmp_file)
    }
    cmp_rows[[cn]] <- cbind(contrast = cn, cmp, stringsAsFactors = FALSE)
  }
}

if (length(map_rows))
  utils::write.csv(do.call(rbind, map_rows),
                   file.path(OUT, "support_map.csv"), row.names = FALSE)
if (length(cmp_rows))
  utils::write.csv(do.call(rbind, cmp_rows),
                   file.path(OUT, "plasmode_design_comparison.csv"),
                   row.names = FALSE)
cr_log("Multi-arm Stage 3 complete for:", paste(run_set, collapse = ", "))

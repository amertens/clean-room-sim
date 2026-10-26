# ============================================================
# 22b. The confirmatory comparison's lock for death within 180 days
# ============================================================
# Spec section 4, run 2: "then the confirmatory comparison's lock for death
# within 180 days". Script 22 never built it, and script 22 cannot change while
# the supervisor relaunches it for each lock, so this script builds the lock
# from the confirmatory comparison's final lock: the same complete-case rows,
# design data, covariates, dictionary, seed and plan, with only the outcome,
# its baseline risk and its tolerances changed (spec section 5 and call 7:
# baseline risk 0.05, bias 0.01, half-width 0.03). Reads design data only; no
# outcome exists in this process.
#
# Usage, after script 23 (pick) has recorded the confirmatory comparison, one
# R job at a time, and before the sensitivity locks so it meets the cut-off:
#   Rscript rescueCo/scripts/22b_secondary_lock.R check            # build the lock only
#   Rscript rescueCo/scripts/22b_secondary_lock.R run <workers>
# Decision 5 of 2026-10-23.
# ============================================================

source(file.path(if (dir.exists("rescueCo")) "rescueCo/R" else "R", "bootstrap.R"))
source("rescueCo/R/utils.R")
source("rescueCo/R/reexec_guards.R")
.libPaths(c("rescueCo/lib", .libPaths()))
suppressMessages(library(cleanTMLE))
rx_check_fingerprint()
qp <- c("C:/Program Files/RStudio/resources/app/bin/quarto/bin/quarto.exe",
        file.path(Sys.getenv("LOCALAPPDATA"), "Programs/Quarto/bin/quarto.exe"))   # RStudio's, or a user install
if (any(file.exists(qp))) Sys.setenv(QUARTO_PATH = qp[file.exists(qp)][1])
cfg <- load_cr_config()
args <- commandArgs(trailingOnly = TRUE)
MODE <- args[1]
WORKERS <- if (length(args) >= 2) as.integer(args[2]) else 1L
stopifnot(MODE %in% c("check", "run"))

conf <- utils::read.csv(file.path(REG_DIR, "confirmatory.csv"), stringsAsFactors = FALSE)
NM <- conf$lock[1]
if (is.na(NM) || !NM %in% paste0("L", 1:5))
  stop("No confirmatory comparison is recorded; run script 23 in pick mode first.", call. = FALSE)
OUT <- file.path(cfg$data$design_dir, "locks", "final")
RES <- file.path(cfg$paths$results, "reexec_050", "final")
NEW <- paste0(NM, "_death180")

base <- load_lock(file.path(OUT, paste0(NM, "_lock.rds")))
p <- base$plan
plan <- analysis_plan(
  outcome = "death_by_6mo", estimands = p$estimands, time_zero = p$time_zero,
  eligibility = p$eligibility, ps_library = p$ps_library, candidates = p$candidates,
  tolerance = list(bias = 0.01, coverage = p$tolerance$coverage, ci_halfwidth = 0.03),
  surfaces = utils::modifyList(p$surfaces, list(baseline_risk = 0.05)),
  trim_band = p$trim_band, negative_controls = p$negative_controls,
  nc_exclude = p$nc_exclude, nc_criteria = p$nc_criteria, roles = p$roles,
  K = p$K, V = p$V, reps = p$reps, max_reps = p$max_reps, seed = p$seed)

# Everything but the outcome, its tolerances and its surfaces must equal the
# confirmatory lock's plan.
same <- setdiff(names(p), c("outcome", "tolerance", "surfaces"))
off <- same[!vapply(same, function(k) identical(plan[[k]], p[[k]]), logical(1))]
if (length(off))
  stop("The rebuilt plan differs from ", NM, "'s plan in: ", paste(off, collapse = ", "), ".",
       call. = FALSE)
lk <- create_analysis_lock(base$data, treatment = base$treatment, covariates = base$covariates,
                           plan = plan, id = base$id, dictionary = base$dictionary)
stopifnot(identical(lk$data, base$data), identical(lk$dictionary, base$dictionary))
rx_log(sprintf("22b %s: built from %s (lock %s); n = %d; plan differs only in outcome, baseline risk and tolerances",
               NEW, NM, substr(base$lock_hash, 1, 12), nrow(lk$data)))
if (MODE == "check") quit(save = "no")

done <- file.path(OUT, paste0(NEW, "_dossier.rds"))
if (file.exists(done)) stop(NEW, " already has a design report.", call. = FALSE)
if (WORKERS > 1L) future::plan(future::multisession, workers = WORKERS)
t0 <- Sys.time()
save_lock(lk, file.path(OUT, paste0(NEW, "_lock.rds")))
ds <- assess_design(lk); saveRDS(ds, file.path(OUT, paste0(NEW, "_design.rds")))
sm <- simulate_design(lk, ds); saveRDS(sm, file.path(OUT, paste0(NEW, "_simulation.rds")))
nc <- negative_control_ladder(lk, ds, sm); saveRDS(nc, file.path(OUT, paste0(NEW, "_nc.rds")))
dos <- design_report(lk, ds, sm, nc, file = file.path(RES, paste0(NEW, "_dossier.html")))
saveRDS(dos, done)
utils::write.csv(dos$identification, file.path(RES, paste0(NEW, "_identification.csv")), row.names = FALSE)
rx_log(sprintf("22b %s: primary %s; NC %s; lock %s; %.1f min", NEW, dos$decision$primary,
               dos$decision$nc_verdict, substr(lk$lock_hash, 1, 12),
               as.numeric(Sys.time() - t0, units = "mins")))

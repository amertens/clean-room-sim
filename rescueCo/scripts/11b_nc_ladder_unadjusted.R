# ============================================================
# Multi-arm Stage 2b: negative-control ladder, unadjusted
# ============================================================
# The main pipeline's negative-control step (analysis/18) is an unadjusted
# two-proportion difference and frames itself as a balance check. The
# TMLE-adjusted ladder from script 11 conditions on the design columns,
# which for the SES controls include near-copies of the control itself
# (cooking_fuel_* for the wood/charcoal control), so its restricted-rung
# estimates are not comparable to the main pipeline's. This script reruns
# the ladder unadjusted on the two ambulance locks and overwrites
# nc_ladder.csv; the adjusted version is kept beside it as
# nc_ladder_tmle.csv with a caveat column.
# ============================================================

suppressMessages(library(cleanTMLE))
source(file.path(if (dir.exists("rescueCo")) "rescueCo/R" else "R",
                 "bootstrap.R"))
source("rescueCo/R/utils.R")
cr_log("=== Stage 2b: unadjusted negative-control ladder ===")

OUT <- "rescueCo/results/multiarm"
locks <- readRDS(file.path(OUT, "multiarm_locks.rds"))

rows <- list()
for (cn in c("C4", "PRIMARY", "C1", "C2", "C3")) {
  lk <- locks[[cn]]
  restr <- if (cn == "C4") list(
    `transfers excluded` = lk$data$is_transfer == 0,
    `transfers and prior care excluded` = lk$data$pre_hospital_care == 0)
  else list()
  lad <- run_negative_control_ladder(lk, restrictions = restr,
                                     method = "unadjusted", verbose = FALSE)
  rows[[cn]] <- cbind(contrast = cn, lad$table, stringsAsFactors = FALSE)
  for (s in lad$turned_null) cr_log(paste("[", cn, "]", s))
}
tab <- do.call(rbind, rows)

# Preserve the adjusted version for the record before overwriting.
old <- file.path(OUT, "nc_ladder.csv")
if (file.exists(old)) {
  adj <- utils::read.csv(old, stringsAsFactors = FALSE)
  if (!"method" %in% names(adj)) {
    adj$method <- "tmle_adjusted"
    adj$caveat <- paste("Adjusted fit conditions on design columns that are",
                        "near-copies of the SES controls; restricted-rung",
                        "estimates are not comparable to the main pipeline.")
    utils::write.csv(adj, file.path(OUT, "nc_ladder_tmle.csv"),
                     row.names = FALSE)
  }
}
tab$method <- "unadjusted"
utils::write.csv(tab, old, row.names = FALSE)
cr_log(paste("Wrote", old, "(", nrow(tab), "rows )"))

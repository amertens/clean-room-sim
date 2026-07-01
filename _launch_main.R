# Wrapper: run main simulation overnight (logs to results_new/_overnight_main.log)
setwd("C:/Users/andre/OneDrive/Documents/clean-room-sim")
log_path <- "results_new/_overnight_main.log"
con <- file(log_path, open = "wt")
sink(con, type = "output"); sink(con, type = "message")
cat(sprintf("=== overnight run started: %s ===\n", Sys.time()))
tryCatch(
  source("run_simulation.R", echo = FALSE),
  error = function(e) cat("FATAL ERROR:", conditionMessage(e), "\n")
)
cat(sprintf("=== finished: %s ===\n", Sys.time()))
sink(type = "message"); sink(type = "output"); close(con)

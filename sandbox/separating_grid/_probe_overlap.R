# Diagnostic probe (not part of the locked study). Verifies:
#  (1) required packages are available,
#  (2) the two DGP functions can be extracted from run_simulation.R without
#      running its main loop,
#  (3) the induced propensity-score range under the candidate DGP configs,
#      to test whether marginal overlap and a nonlinear surface can co-occur.

pk <- c("cleanTMLE","SuperLearner","glmnet","gam","digest","yaml","sandwich","nnet")
cat("== package availability ==\n")
for (p in pk) cat(sprintf("  %-14s %s\n", p, ifelse(requireNamespace(p, quietly=TRUE), "OK", "MISSING")))

# Extract ONLY generate_data + compute_truth from run_simulation.R (no side effects).
sim_path <- "C:/Users/andre/OneDrive/Documents/clean-room-sim/run_simulation.R"
exprs <- parse(file = sim_path)
got <- character(0)
for (e in exprs) {
  if (is.call(e) && (identical(e[[1]], as.name("<-")) || identical(e[[1]], as.name("=")))) {
    lhs <- e[[2]]
    if (is.name(lhs) && as.character(lhs) %in% c("generate_data","compute_truth")) {
      eval(e, envir = globalenv()); got <- c(got, as.character(lhs))
    }
  }
}
cat("\n== extracted functions ==\n  ", paste(got, collapse=", "), "\n")

ps_report <- function(label, df_args) {
  d <- do.call(generate_data, c(list(n = 8000, seed = 101), df_args))
  ps <- predict(glm(treatment ~ age + sex + biomarker + comorbidity + ckd,
                    data = d, family = binomial()), type = "response")
  cat(sprintf("\n%-42s PS range [%.3f, %.3f]  q01=%.3f q99=%.3f  treated=%.2f  Yrate=%.3f\n",
              label, min(ps), max(ps), quantile(ps,.01), quantile(ps,.99),
              mean(d$treatment), mean(d$event_24)))
}

cat("\n== induced propensity-score ranges ==")
ps_report("misspec=FALSE, overlap_strength=0.5", list(overlap_strength=0.5, misspec=FALSE))
ps_report("misspec=FALSE, overlap_strength=1.5", list(overlap_strength=1.5, misspec=FALSE))
ps_report("misspec=TRUE,  overlap_strength=1.5", list(overlap_strength=1.5, misspec=TRUE))
ps_report("misspec=TRUE,  overlap_strength=0.5", list(overlap_strength=0.5, misspec=TRUE))

cat("\n== true RD ==\n")
cat(sprintf("  misspec=FALSE overlap=1.5: RD=%.5f\n",
            compute_truth(1e5, 1.5, -0.05, 7, misspec=FALSE)$RD))
cat(sprintf("  misspec=TRUE  overlap=1.5: RD=%.5f\n",
            compute_truth(1e5, 1.5, -0.05, 7, misspec=TRUE)$RD))
cat("\nDONE\n")

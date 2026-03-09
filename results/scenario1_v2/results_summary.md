
> cat("# Scenario 1 v2: Bad Overlap, TMLE OK — Results Summary\n\n")
# Scenario 1 v2: Bad Overlap, TMLE OK — Results Summary


> cat(sprintf("**True RD:** %.6f\n\n", truth$RD_true))
**True RD:** 0.077603


> cat(sprintf("**Replicates:** %d (N = %d)\n\n", params$reps, params$n_main))
**Replicates:** 200 (N = 500)


> cat("## Estimator Performance (all replicates)\n\n")
## Estimator Performance (all replicates)


> cat("| Method | Bias | RMSE | Emp SD | Mean SE | Coverage | MCSE Cov |\n")
| Method | Bias | RMSE | Emp SD | Mean SE | Coverage | MCSE Cov |

> cat("|--------|------|------|--------|---------|----------|----------|\n")
|--------|------|------|--------|---------|----------|----------|

> for (i in seq_len(nrow(metrics_all))) {
+   m <- metrics_all[i, ]
+   cat(sprintf("| %s | %.4f | %.4f | %.4f | %s | %s | %s |\n",
+               m$ .... [TRUNCATED] 
| Crude | 0.1084 | 0.1143 | 0.0364 | 0.0381 | 17.00% | 0.0266 |
| Adj. Regression | 0.0005 | 0.0528 | 0.0529 | 0.0013 | 4.50% | 0.0147 |
| IPTW | 0.0224 | 0.1097 | 0.1077 | 0.0858 | 85.50% | 0.0249 |
| PS-Matched Reg. | 0.0026 | 0.0652 | 0.0654 | 0.0024 | 5.00% | 0.0154 |
| TMLE (IC) | -0.0077 | 0.1043 | 0.1043 | 0.0707 | 80.00% | 0.0283 |
| TMLE (tboot) | -0.0077 | 0.1043 | 0.1043 | 0.0708 | 74.00% | 0.0310 |
| TMLE (npboot) | -0.0077 | 0.1043 | 0.1043 | — | 0.00% | 0.0000 |

> cat("\n## Stop/Go Gate\n\n")

## Stop/Go Gate


> cat(sprintf("Replicates flagged 'bad' overlap: %d/%d (%.1f%%)\n\n",
+             n_flagged, params$reps, pct_flagged))
Replicates flagged 'bad' overlap: 200/200 (100.0%)


> cat("| Flag | N | TMLE Bias | IC Cov | tboot Cov | Acceptable |\n")
| Flag | N | TMLE Bias | IC Cov | tboot Cov | Acceptable |

> cat("|------|---|-----------|--------|-----------|------------|\n")
|------|---|-----------|--------|-----------|------------|

> for (i in seq_len(nrow(gate_results))) {
+   g <- gate_results[i, ]
+   cat(sprintf("| %s | %s | %s | %s | %s | %s |\n",
+               g$overlap_f .... [TRUNCATED] 
| bad | 200 | -0.00766 | 0.800 | 0.740 | FALSE |
| ok | 0 | — | — | — | — |

> sink()

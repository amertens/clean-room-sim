# source_map_claims.R -- claim registry for build_source_map.R (WP0).
#
# Each claim binds ONE literal number printed in a document to ONE field
# of ONE result file, and build_source_map.R verifies the binding by
# recomputing it. Columns of a claim:
#   document : "manuscript" or "presentation" (extend as documents split)
#   anchor   : regex that must match exactly one source line, the line
#              carrying the printed value
#   value    : the number exactly as printed (may carry +, %, commas)
#   source   : result file path, relative to the repository root
#   extract  : one-argument function; receives the loaded object (a data
#              frame for .csv, the object for .rds) and returns a single
#              numeric, on the printed scale
#   note     : short provenance remark (optional)
#
# Verification rule: |extracted - printed| <= 0.5 * 10^-d + 1e-9, where d
# is the number of printed decimals. Inline `r` expressions never need a
# claim (they are computed at render); claims cover literals only.

claim <- function(document, anchor, value, source, extract, note = "") {
  list(document = document, anchor = anchor, value = value,
       source = source, extract = extract, note = note)
}

.ma  <- "rescueCo/results/multiarm"
.rr  <- "rescueCo/results"

claims <- list(

  # -- Walkthrough and abstract: C1 support headline ------------------------
  claim("manuscript", "outside \\[0\\.05, 0\\.95\\]; max weight 257", "40.6",
        file.path(.ma, "support_by_contrast.csv"),
        function(df) df$pct_outside_band[df$contrast == "C1"]),
  claim("manuscript", "outside \\[0\\.05, 0\\.95\\]; max weight 257", "257",
        file.path(.ma, "support_by_contrast.csv"),
        function(df) df$max_iptw_weight[df$contrast == "C1"]),
  claim("manuscript", "placed 41 percent of the sample", "41",
        file.path(.ma, "support_by_contrast.csv"),
        function(df) df$pct_outside_band[df$contrast == "C1"],
        "abstract rounds 40.57 to 41"),
  claim("manuscript", "five transport\\s*$|on 8,323 patients", "8,323",
        file.path(.ma, "lock_summary.csv"),
        function(df) df$n[df$contrast == "C1"]),

  # -- C1 violation regions (tree search) -----------------------------------
  claim("manuscript", "injured at home \\(1,386 patients", "1,386",
        file.path(.ma, "violation_regions.csv"),
        function(df) df$n[df$contrast == "C1"][1]),
  claim("manuscript", "injured at home \\(1,386 patients", "1.4",
        file.path(.ma, "violation_regions.csv"),
        function(df) 100 * df$p_treated[df$contrast == "C1"][1]),
  claim("manuscript", "\\(4.2% treated\\)", "4.2",
        file.path(.ma, "violation_regions.csv"),
        function(df) 100 * df$p_treated[df$contrast == "C1" & df$n == 1062]),

  # -- C1 trimmed-population profile ---------------------------------------
  claim("manuscript", "retains 4,946 patients after removing 3,377", "4,946",
        file.path(.ma, "who_is_unsupported.csv"),
        function(df) df$n_kept[df$contrast == "C1" & df$population == "all"][1]),
  claim("manuscript", "retains 4,946 patients after removing 3,377", "3,377",
        file.path(.ma, "who_is_unsupported.csv"),
        function(df) df$n_removed[df$contrast == "C1" & df$population == "all"][1]),
  claim("manuscript", "injured at home \\(40% versus 1%", "40",
        file.path(.ma, "who_is_unsupported.csv"),
        function(df) 100 * df$mean_removed[df$contrast == "C1" &
          df$population == "all" & grepl("Private.House", df$variable)][1]),
  claim("manuscript", "injured at home \\(40% versus 1%", "1.10",
        file.path(.ma, "who_is_unsupported.csv"),
        function(df) df$smd[df$contrast == "C1" & df$population == "all" &
          grepl("Private.House", df$variable)][1]),
  claim("manuscript", "assaults \\(39% versus", "39",
        file.path(.ma, "who_is_unsupported.csv"),
        function(df) 100 * df$mean_removed[df$contrast == "C1" &
          df$population == "all" & grepl("assault", df$variable)][1]),
  claim("manuscript", "assaults \\(39% versus", "10",
        file.path(.ma, "who_is_unsupported.csv"),
        function(df) 100 * df$mean_kept[df$contrast == "C1" &
          df$population == "all" & grepl("assault", df$variable)][1]),

  # -- Negative-control ladder (case study) ---------------------------------
  claim("manuscript", "urban\\s*$|urban household \\+0.068", "+0.068",
        file.path(.ma, "nc_ladder.csv"),
        function(df) df$estimate[df$contrast == "C4" &
          df$negative_control == "nc_household_urban" &
          df$cohort == "full cohort"]),
  claim("manuscript", "95% CI 0.049 to 0.086", "0.049",
        file.path(.ma, "nc_ladder.csv"),
        function(df) df$ci_lower[df$contrast == "C4" &
          df$negative_control == "nc_household_urban" &
          df$cohort == "full cohort"]),
  claim("manuscript", "95% CI 0.049 to 0.086", "0.086",
        file.path(.ma, "nc_ladder.csv"),
        function(df) df$ci_upper[df$contrast == "C4" &
          df$negative_control == "nc_household_urban" &
          df$cohort == "full cohort"]),
  claim("manuscript", "0.052, CI \\$-\\$0.074 to \\$-\\$0.029|0\\.052, CI",
        "0.052",
        file.path(.ma, "nc_ladder.csv"),
        function(df) abs(df$estimate[df$contrast == "C4" &
          df$negative_control == "nc_cooking_fuel_wood_charcoal" &
          df$cohort == "full cohort"])),
  claim("manuscript", "household \\+0.001", "+0.001",
        file.path(.ma, "nc_ladder.csv"),
        function(df) df$estimate[df$contrast == "C4" &
          df$negative_control == "nc_household_urban" &
          df$cohort == "transfers excluded"]),
  claim("manuscript", "urban-household control fails on the full\\s*$|cohort \\(\\+0.013, CI 0.003 to 0.024\\)",
        "+0.013",
        file.path(.ma, "nc_ladder.csv"),
        function(df) df$estimate[df$contrast == "C1" &
          df$negative_control == "nc_household_urban" &
          df$cohort == "full cohort"]),

  # -- Generate- versus sample-treatment comparison (designcmp) -------------
  claim("manuscript", "ATT 0.018 versus", "0.018",
        file.path(.ma, "designcmp_C1.rds"),
        function(x) x$bias[x$confounding == 2 & x$modification == 0 &
          x$design == "sample_treatment" & x$estimand == "ATT"]),
  claim("manuscript", "ATT 0.018 versus", "0.013",
        file.path(.ma, "designcmp_C1.rds"),
        function(x) x$bias[x$confounding == 2 & x$modification == 0 &
          x$design == "generate_treatment" & x$estimand == "ATT"]),
  claim("manuscript", "ATO 0.016 versus 0.010", "0.016",
        file.path(.ma, "designcmp_C1.rds"),
        function(x) x$bias[x$confounding == 2 & x$modification == 0 &
          x$design == "sample_treatment" & x$estimand == "ATO"]),
  claim("manuscript", "ATO 0.016 versus 0.010", "0.010",
        file.path(.ma, "designcmp_C1.rds"),
        function(x) x$bias[x$confounding == 2 & x$modification == 0 &
          x$design == "generate_treatment" & x$estimand == "ATO"]),
  claim("manuscript", "0.69 versus 0.78", "0.69",
        file.path(.ma, "designcmp_C1.rds"),
        function(x) x$coverage[x$confounding == 2 & x$modification == 0 &
          x$design == "sample_treatment" & x$estimand == "ATT"]),
  claim("manuscript", "0.69 versus 0.78", "0.78",
        file.path(.ma, "designcmp_C1.rds"),
        function(x) x$coverage[x$confounding == 2 & x$modification == 0 &
          x$design == "generate_treatment" & x$estimand == "ATT"]),
  claim("manuscript", "ATO 0.64 versus 0.79", "0.64",
        file.path(.ma, "designcmp_C1.rds"),
        function(x) x$coverage[x$confounding == 2 & x$modification == 0 &
          x$design == "sample_treatment" & x$estimand == "ATO"]),
  claim("manuscript", "ATO 0.64 versus 0.79", "0.79",
        file.path(.ma, "designcmp_C1.rds"),
        function(x) x$coverage[x$confounding == 2 & x$modification == 0 &
          x$design == "generate_treatment" & x$estimand == "ATO"]),
  claim("manuscript", "matched ATT\\s*$|0.63 versus\\s*$|0.86\\), while the ATE", "0.63",
        file.path(.ma, "designcmp_C1.rds"),
        function(x) x$coverage[x$confounding == 2 & x$modification == 0 &
          x$design == "sample_treatment" & x$estimand == "matched_ATT"]),
  claim("manuscript", "0.86\\), while the ATE", "0.86",
        file.path(.ma, "designcmp_C1.rds"),
        function(x) x$coverage[x$confounding == 2 & x$modification == 0 &
          x$design == "generate_treatment" & x$estimand == "matched_ATT"]),
  claim("manuscript", "balloons more than fivefold \\(0.066", "0.066",
        file.path(.ma, "designcmp_C1.rds"),
        function(x) x$rmse[x$confounding == 2 & x$modification == 0 &
          x$design == "sample_treatment" & x$estimand == "ATE"]),
  claim("manuscript", "balloons more than fivefold \\(0.066", "0.012",
        file.path(.ma, "designcmp_C1.rds"),
        function(x) x$rmse[x$confounding == 2 & x$modification == 0 &
          x$design == "generate_treatment" & x$estimand == "ATE"]),

  # -- Section 12.6 worst-case RMSE ratios ----------------------------------
  # These verify against the surviving 50-replicate BACKUP csv; the files
  # the reproducibility table names hold 3-replicate smoke values (Stage 0
  # finding F6). The binding is recorded here so the defect stays visible
  # until decision D3's 200-replicate rerun replaces the section.
  claim("manuscript", "worst-case RMSE ratio relative to baseline is 1.82", "1.82",
        file.path(.rr, "_dq_degradation_50rep_backup.csv"),
        function(df) max(df$rmse_ratio[df$candidate == "glm_t01"]),
        "F6: only the 50-rep backup reproduces this"),
  claim("manuscript", "1.73 for `glmnet_t01`", "1.73",
        file.path(.rr, "_dq_degradation_50rep_backup.csv"),
        function(df) max(df$rmse_ratio[df$candidate == "glmnet_t01"]),
        "F6"),
  claim("manuscript", "1.53 for `ensemble_t01`", "1.53",
        file.path(.rr, "_dq_degradation_50rep_backup.csv"),
        function(df) max(df$rmse_ratio[df$candidate == "ensemble_t01"]),
        "F6"),
  claim("manuscript", "1.27 for\\s*$|`ensemble_t05`; the largest", "1.27",
        file.path(.rr, "_dq_degradation_50rep_backup.csv"),
        function(df) max(df$rmse_ratio[df$candidate == "ensemble_t05"]),
        "F6"),

  # -- Mortality and time-to-arrival narrative ------------------------------
  claim("manuscript", "C4 primary row is \\$-2.9\\$", "-2.9",
        file.path(.ma, "ladder_estimates.csv"),
        function(df) 100 * df$estimate[df$contrast == "C4" &
          df$outcome == "death_by_6mo" & df$estimand == "ATE"]),
  claim("manuscript", "95% CI \\$-4.9\\$ to", "-4.9",
        file.path(.ma, "ladder_estimates.csv"),
        function(df) 100 * df$ci_lower[df$contrast == "C4" &
          df$outcome == "death_by_6mo" & df$estimand == "ATE"]),
  claim("manuscript", "\\$-0.9\\$\\), and on the primary", "-0.9",
        file.path(.ma, "ladder_estimates.csv"),
        function(df) 100 * df$ci_upper[df$contrast == "C4" &
          df$outcome == "death_by_6mo" & df$estimand == "ATE"]),
  claim("manuscript", "treated \\(ATT \\$-2.8\\$", "-2.8",
        file.path(.ma, "ladder_estimates.csv"),
        function(df) 100 * df$estimate[df$contrast == "PRIMARY" &
          df$outcome == "death_by_6mo" & df$estimand == "ATT"]),
  claim("manuscript", "ATE \\$-1.7\\$", "-1.7",
        file.path(.ma, "ladder_estimates.csv"),
        function(df) 100 * df$estimate[df$contrast == "PRIMARY" &
          df$outcome == "death_by_6mo" & df$estimand == "ATE"]),
  claim("manuscript", "higher mortality \\(\\+1.8 points", "+1.8",
        file.path(.ma, "ladder_estimates.csv"),
        function(df) 100 * df$estimate[df$contrast == "C3" &
          df$outcome == "death_by_6mo" & df$estimand == "ATE"]),
  claim("manuscript", "CI 0.8 to 2.8", "0.8",
        file.path(.ma, "ladder_estimates.csv"),
        function(df) 100 * df$ci_lower[df$contrast == "C3" &
          df$outcome == "death_by_6mo" & df$estimand == "ATE"]),
  claim("manuscript", "CI 0.8 to 2.8", "2.8",
        file.path(.ma, "ladder_estimates.csv"),
        function(df) 100 * df$ci_upper[df$contrast == "C3" &
          df$outcome == "death_by_6mo" & df$estimand == "ATE"]),
  claim("manuscript", "C3\\s*$|\\$-245\\$ minutes", "-245",
        file.path(.ma, "ladder_estimates.csv"),
        function(df) df$estimate[df$contrast == "C3" &
          df$outcome == "time_to_arrival_minutes" & df$estimand == "ATE"]),
  claim("manuscript", "CI \\$-320\\$ to \\$-169\\$", "-320",
        file.path(.ma, "ladder_estimates.csv"),
        function(df) df$ci_lower[df$contrast == "C3" &
          df$outcome == "time_to_arrival_minutes" & df$estimand == "ATE"]),
  claim("manuscript", "CI \\$-320\\$ to \\$-169\\$", "-169",
        file.path(.ma, "ladder_estimates.csv"),
        function(df) df$ci_upper[df$contrast == "C3" &
          df$outcome == "time_to_arrival_minutes" & df$estimand == "ATE"]),
  claim("manuscript", "null \\$-0.1\\$ minutes", "-0.1",
        file.path(.ma, "ladder_estimates.csv"),
        function(df) df$estimate[df$contrast == "PRIMARY" &
          df$outcome == "time_to_arrival_minutes" & df$estimand == "ATE"]),
  claim("manuscript", "\\$-32\\$ to \\$31\\$", "-32",
        file.path(.ma, "ladder_estimates.csv"),
        function(df) df$ci_lower[df$contrast == "PRIMARY" &
          df$outcome == "time_to_arrival_minutes" & df$estimand == "ATE"]),
  claim("manuscript", "\\$-32\\$ to \\$31\\$", "31",
        file.path(.ma, "ladder_estimates.csv"),
        function(df) df$ci_upper[df$contrast == "PRIMARY" &
          df$outcome == "time_to_arrival_minutes" & df$estimand == "ATE"]),

  # -- Presentation: lock table and estimand table --------------------------
  claim("presentation", "^\\| C1 \\| Rescue.Co vs everyone else", "8,323",
        file.path(.ma, "lock_summary.csv"),
        function(df) df$n[df$contrast == "C1"]),
  claim("presentation", "^\\| C1 \\| Rescue.Co vs everyone else", "1,039",
        file.path(.ma, "lock_summary.csv"),
        function(df) df$n_treated[df$contrast == "C1"]),
  claim("presentation", "^\\| C2 \\|", "7,165",
        file.path(.ma, "lock_summary.csv"),
        function(df) df$n[df$contrast == "C2"]),
  claim("presentation", "^\\| C3 \\|", "2,197",
        file.path(.ma, "lock_summary.csv"),
        function(df) df$n_treated[df$contrast == "C3"]),
  claim("presentation", "^\\| C4 \\|", "2,197",
        file.path(.ma, "lock_summary.csv"),
        function(df) df$n[df$contrast == "C4"]),
  claim("presentation", "^\\| Primary \\|", "1,616",
        file.path(.ma, "lock_summary.csv"),
        function(df) df$n[df$contrast == "PRIMARY"]),
  claim("presentation", "^\\| Primary \\|", "1,003",
        file.path(.ma, "lock_summary.csv"),
        function(df) df$n_treated[df$contrast == "PRIMARY"]),
  claim("presentation", "^\\| ATE \\| 257 \\|", "257",
        file.path(.ma, "estimand_feasibility.csv"),
        function(df) df$max_weight[df$contrast == "C1" & df$estimand == "ATE"]),
  claim("presentation", "^\\| Trimmed ATE \\| 20 \\|", "20",
        file.path(.ma, "estimand_feasibility.csv"),
        function(df) df$max_weight[df$contrast == "C1" &
          grepl("^trimmed_ATE", df$estimand)][1]),
  claim("presentation", "^\\| ATT \\| 5.0 \\|", "5.0",
        file.path(.ma, "estimand_feasibility.csv"),
        function(df) df$max_weight[df$contrast == "C1" & df$estimand == "ATT"]),
  claim("presentation", "^\\| ATO \\| 1.0 \\|", "1.0",
        file.path(.ma, "estimand_feasibility.csv"),
        function(df) df$max_weight[df$contrast == "C1" & df$estimand == "ATO"]),
  claim("presentation", "8,323 trauma patients with known transport", "8,323",
        file.path(.ma, "lock_summary.csv"),
        function(df) df$n[df$contrast == "C1"]),
  claim("presentation", "injured at a private home, 1,386 patients", "1,386",
        file.path(.ma, "violation_regions.csv"),
        function(df) df$n[df$contrast == "C1"][1]),
  claim("presentation", "injured at a private home, 1,386 patients", "1.4",
        file.path(.ma, "violation_regions.csv"),
        function(df) 100 * df$p_treated[df$contrast == "C1"][1]),
  claim("presentation", "other non-street injury places, 1,438", "1,438",
        file.path(.ma, "violation_regions.csv"),
        function(df) df$n[df$contrast == "C1" & df$n == 1438]),
  claim("presentation", "other non-street injury places, 1,438", "2.5",
        file.path(.ma, "violation_regions.csv"),
        function(df) 100 * df$p_treated[df$contrast == "C1" & df$n == 1438]),
  claim("presentation", "assault victims elsewhere, 1,062", "1,062",
        file.path(.ma, "violation_regions.csv"),
        function(df) df$n[df$contrast == "C1" & df$n == 1062]),
  claim("presentation", "assault victims elsewhere, 1,062", "4.2",
        file.path(.ma, "violation_regions.csv"),
        function(df) 100 * df$p_treated[df$contrast == "C1" & df$n == 1062]),
  claim("presentation", "A trimmed analysis removes 3,377 people", "3,377",
        file.path(.ma, "who_is_unsupported.csv"),
        function(df) df$n_removed[df$contrast == "C1" &
          df$population == "all"][1])
)

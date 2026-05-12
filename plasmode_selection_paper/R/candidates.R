# Candidate-grid construction for the plasmode-selection paper.
#
# Three nuisance-library tiers (parametric / parametric+ML / ML+screener)
# crossed with two PS truncation levels and two screener choices, plus
# a "discrete vs convex" SL axis. The grid is paper-sized but bounded
# below 32 candidates per analysis so it stays tractable on a workstation.

library(cleanTMLE)

# --- 1. Nuisance libraries -------------------------------------------

#' The three library tiers used in the paper. Each returns a character
#' vector compatible with SuperLearner::SuperLearner.
NUISANCE_LIBRARIES <- list(
  # Tier 1: parametric only. Represents a traditional analyst who
  # picks a GLM by reflex and does not invoke any ML learner.
  parametric = c("SL.glm", "SL.mean"),

  # Tier 2: parametric + glmnet. The minimum-effort ML library;
  # representative of a "good practice" applied analyst.
  parametric_plus_glmnet = c("SL.glm", "SL.glmnet", "SL.mean"),

  # Tier 3: full SL library with smoothers and trees. Representative
  # of a workshop-recommended default. Adds GAM and ranger.
  rich = c("SL.glm", "SL.glmnet", "SL.gam", "SL.ranger", "SL.mean"),

  # Tier 4: rich library wrapped with a correlation screener.
  rich_screener = list(
    c("SL.glm", "screen.corP"),
    c("SL.glmnet", "screen.corP"),
    c("SL.gam", "screen.corP"),
    c("SL.ranger", "screen.corP"),
    c("SL.mean", "All")
  )
)

# --- 2. The candidate grid --------------------------------------------

#' Build the paper's full candidate grid.
#'
#' Six candidates per analysis, capturing the contrast the paper makes:
#' Tier 1 vs Tier 2 vs Tier 3 vs Tier 4, with two truncation levels on
#' the rich and rich+screener tiers.
build_candidate_grid <- function() {
  cands <- list()

  cands[["param_t01"]] <- tmle_candidate(
    candidate_id = "param_t01",
    label        = "Parametric, trunc = 0.01",
    g_library    = NUISANCE_LIBRARIES$parametric,
    q_library    = NUISANCE_LIBRARIES$parametric,
    truncation   = 0.01
  )

  cands[["paramplus_t01"]] <- tmle_candidate(
    candidate_id = "paramplus_t01",
    label        = "Parametric + glmnet, trunc = 0.01",
    g_library    = NUISANCE_LIBRARIES$parametric_plus_glmnet,
    q_library    = NUISANCE_LIBRARIES$parametric_plus_glmnet,
    truncation   = 0.01
  )

  cands[["rich_t01"]] <- tmle_candidate(
    candidate_id = "rich_t01",
    label        = "Rich SL library, trunc = 0.01",
    g_library    = NUISANCE_LIBRARIES$rich,
    q_library    = NUISANCE_LIBRARIES$rich,
    truncation   = 0.01
  )

  cands[["rich_t05"]] <- tmle_candidate(
    candidate_id = "rich_t05",
    label        = "Rich SL library, trunc = 0.05",
    g_library    = NUISANCE_LIBRARIES$rich,
    q_library    = NUISANCE_LIBRARIES$rich,
    truncation   = 0.05
  )

  cands[["rich_screener_t01"]] <- tmle_candidate(
    candidate_id = "rich_screener_t01",
    label        = "Rich + corP screener, trunc = 0.01",
    g_library    = NUISANCE_LIBRARIES$rich_screener,
    q_library    = NUISANCE_LIBRARIES$rich_screener,
    truncation   = 0.01,
    screener     = "corP"
  )

  cands[["rich_screener_t05"]] <- tmle_candidate(
    candidate_id = "rich_screener_t05",
    label        = "Rich + corP screener, trunc = 0.05",
    g_library    = NUISANCE_LIBRARIES$rich_screener,
    q_library    = NUISANCE_LIBRARIES$rich_screener,
    truncation   = 0.05,
    screener     = "corP"
  )

  cands
}

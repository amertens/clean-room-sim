# build_source_map.R -- claim-level source map for the revision (WP0).
#
# Regenerates revision/SOURCE_MAP.csv and revision/UNMATCHED_NUMBERS.csv.
#
# Model (decision D10): every number a reader sees must be either
#   (a) computed at render time (an inline `r` expression or a chunk that
#       reads a named result file): recorded here at chunk granularity;
#   (b) a literal covered by a VERIFIED claim in
#       revision/source_map_claims.R, which binds the printed value to a
#       named field of a named result file and is recomputed on every run
#       (a claim that stops matching its source FAILS the build); or
#   (c) a classified non-result token (year, citation, section reference,
#       version, small structural count, declared grid parameter).
# Every remaining literal lands in UNMATCHED_NUMBERS.csv.
#
# Exit status: any FAILED claim exits 1 unconditionally. --strict also
# exits 1 while UNMATCHED_NUMBERS.csv is non-empty (the WP2 completion
# gate; operating rule: the file must be empty at the end of WP2).
#
# Usage:
#   Rscript revision/build_source_map.R            # report and write CSVs
#   Rscript revision/build_source_map.R --strict   # WP2 completion gate

options(width = 220, stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
strict <- "--strict" %in% args

root <- local({
  d <- getwd()
  while (!file.exists(file.path(d, "clean-room-sim.Rproj")) && dirname(d) != d)
    d <- dirname(d)
  d
})
setwd(root)
commit <- tryCatch(system("git rev-parse --short HEAD", intern = TRUE),
                   error = function(e) NA_character_)

documents <- c(
  manuscript   = "reports/manuscript_outcome_blind_dq.qmd",
  presentation = "reports/cleanTMLE_presentation.qmd"
)
# Producing scripts for the result files (from the manuscript's own
# script-to-output map, corrected by the Stage 0 audit).
script_of <- function(f) {
  rules <- c(
    "results_new/(simulation_results|simulation_config|summary_|plasmode_|dq_stress_|audit_|decision_log_|interim_|done_|lock_)" = "run_simulation.R",
    "results_new/gate_oc"            = "run_gate_operating_characteristics.R",
    "results_new/bootstrap_variance" = "_bootstrap_variance.R",
    "sandbox/validation"             = "sandbox/validation/validate_vs_tmle.R",
    "sandbox/candidate_divergence"   = "sandbox/candidate_divergence/divergence_study.R",
    "multiarm/(multiarm_|lock_summary)" = "rescueCo/scripts/10_multiarm_build.R",
    "multiarm/(support_by_contrast|estimand_feasibility|who_is|violation|nc_ladder|collider|near_deterministic|design_)" = "rescueCo/scripts/11_multiarm_design.R",
    "multiarm/(supportsim|designcmp|support_map)" = "rescueCo/scripts/12_multiarm_support_sim.R",
    "multiarm/(ladder_|sensitivity_evalues|sl_propensity)" = "rescueCo/scripts/13_multiarm_estimation.R",
    "rescueCo/results/(plasmode_dq|_dq_degradation|stage2b_dq)" = "rescueCo/scripts/rescueco_dq_full.R",
    "rescueCo/results/plasmode_fidelity" = "sandbox/rescueco_fidelity/plasmode_fidelity.R",
    "rescueCo/results/reconciliation" = "rescueCo/scripts/build_reconciliation_table.R"
  )
  for (re in names(rules)) if (grepl(re, f)) return(unname(rules[re]))
  ""
}

# ---- load claims ----------------------------------------------------------
source("revision/source_map_claims.R")  # defines claims (a list) and claim()

.load_cache <- new.env()
load_source <- function(path) {
  key <- path
  if (!is.null(.load_cache[[key]])) return(.load_cache[[key]])
  if (!file.exists(path)) stop("claim source file missing: ", path)
  obj <- if (grepl("\\.csv$", path)) read.csv(path, check.names = FALSE)
         else readRDS(path)
  .load_cache[[key]] <- obj
  obj
}

norm_num <- function(s) suppressWarnings(as.numeric(gsub("[,+%]", "", s)))
decimals_of <- function(s) {
  s <- gsub("[,+%]", "", s)
  if (grepl("\\.", s)) nchar(sub("^-?[0-9]*\\.", "", s)) else 0L
}

doc_lines <- lapply(documents, readLines, warn = FALSE, encoding = "UTF-8")

verify_claims <- function() {
  rows <- list()
  for (i in seq_along(claims)) {
    cl <- claims[[i]]
    lines <- doc_lines[[cl$document]]
    hit <- grep(cl$anchor, lines)
    status <- "OK"; msg <- ""; got <- NA_real_
    if (!length(hit)) {
      status <- "ANCHOR_MISSING"
      msg <- "anchor matches no line"
    } else {
      got <- tryCatch({
        v <- cl$extract(load_source(cl$source))
        if (length(v) != 1L || !is.finite(v))
          stop("extractor returned ", length(v), " value(s)")
        as.numeric(v)
      }, error = function(e) { status <<- "EXTRACT_ERROR"
                               msg <<- conditionMessage(e); NA_real_ })
      if (status == "OK") {
        want <- norm_num(cl$value)
        tol <- 0.5 * 10^(-decimals_of(cl$value)) + 1e-9
        if (!is.finite(want)) { status <- "BAD_VALUE"; msg <- "unparseable printed value" }
        else if (abs(got - want) > tol) {
          status <- "FAILED"
          msg <- sprintf("printed %s but source gives %.6g (tol %.4g)",
                         cl$value, got, tol)
        }
      }
    }
    rows[[i]] <- data.frame(
      idx = i, document = cl$document,
      line = if (length(hit)) hit[1] else NA_integer_,
      value = cl$value, source_file = cl$source,
      source_field = paste(deparse(body(cl$extract)), collapse = " "),
      note = cl$note, status = status, computed = got, message = msg)
  }
  do.call(rbind, rows)
}
cv <- verify_claims()

# ---- token census ---------------------------------------------------------
extract_tokens <- function(lines, docname) {
  n <- length(lines); drop <- rep(FALSE, n)
  if (lines[1] == "---") { e <- which(lines[-1] == "---")[1] + 1L; drop[1:e] <- TRUE }
  in_chunk <- FALSE
  for (i in seq_len(n)) {
    if (grepl("^```", lines[i])) { in_chunk <- !in_chunk; drop[i] <- TRUE; next }
    if (in_chunk) drop[i] <- TRUE
  }
  out <- list()
  for (i in which(!drop)) {
    ln2 <- gsub("`r[^`]*`", " INLINE ", lines[i])
    ln2 <- gsub("@[A-Za-z0-9_:-]+", " ", ln2)
    ms <- gregexpr("[+-]?\\d{1,3}(,\\d{3})+(\\.\\d+)?|[+-]?\\d*\\.\\d+|[+-]?\\d+",
                   ln2, perl = TRUE)[[1]]
    if (ms[1] == -1) next
    toks <- regmatches(ln2, list(ms))[[1]]
    for (k in seq_along(toks)) {
      st <- ms[k]
      out[[length(out) + 1L]] <- data.frame(
        document = docname, line = i, token = toks[k],
        context = substr(ln2, max(1, st - 45),
                         min(nchar(ln2), st + attr(ms, "match.length")[k] + 45)))
    }
  }
  do.call(rbind, out)
}
toks <- do.call(rbind, Map(extract_tokens, doc_lines, names(documents)))
toks$numeric <- norm_num(toks$token)

classify <- function(toks) {
  cls <- rep("candidate", nrow(toks)); ctx <- toks$context; tk <- toks$token
  num <- toks$numeric
  is_year <- num >= 1900 & num <= 2100 & !grepl("\\.", tk)
  cls[is_year] <- "year_or_citation"
  cls[grepl("E9\\(R1\\)|ICH|M14", ctx) & num %in% c(9, 1, 14)] <- "spec_id"
  cls[grepl("0\\.1\\.|0\\.2\\.|version|v0\\.", ctx, ignore.case = TRUE) & num < 10] <- "version"
  cls[grepl("Section|@sec|Stage|Step|Check Point|slide|Appendix", ctx, ignore.case = TRUE) &
        num < 30 & !grepl("\\.", tk)] <- "structural_ref"
  cls[num %in% c(0, 1) & !grepl("\\.", tk) & cls == "candidate"] <- "small_ordinal"
  cls[num %in% 2:12 & !grepl("\\.", tk) & cls == "candidate" &
        grepl("five|four|three|two|six|seven|eight|steps?|verbs?|threat|famil|scenario|estimator|fold|arm|hospital|decimal|stage|section|candidate|learner|column|item|part|reason|way|percentage points|contrast",
              ctx, ignore.case = TRUE)] <- "count_smallint"
  cls[cls == "candidate" &
        grepl("threshold|band|\\[0\\.05|0\\.95\\]|trunc|tolerance|floor|window|alpha|prevalence|OR|sens|spec|slope|fraction|seed|reps|replicat|n = |N = |B = ",
              ctx, ignore.case = TRUE)] <- "declared_parameter"
  cls
}
toks$class <- classify(toks)

# Coverage by verified claims: same document, same normalised value, within
# two lines of the claim's anchor line.
ok <- cv[cv$status == "OK", ]
covered <- rep(FALSE, nrow(toks))
for (j in seq_len(nrow(ok))) {
  same <- toks$document == ok$document[j] &
    abs(toks$numeric - norm_num(ok$value[j])) < 1e-9 &
    !is.na(ok$line[j]) & abs(toks$line - ok$line[j]) <= 2
  covered <- covered | same
}
toks$covered <- covered

# ---- computed content (chunk file reads) ---------------------------------
chunk_reads <- function(lines, docname) {
  # Any quoted .csv/.rds path inside a code chunk counts as a render-time
  # read (covers read.csv, readRDS, file.path pieces, and helper wrappers).
  in_chunk <- FALSE
  rows <- list()
  for (i in seq_along(lines)) {
    if (grepl("^```", lines[i])) { in_chunk <- !in_chunk; next }
    if (!in_chunk) next
    paths <- unlist(regmatches(lines[i],
                               gregexpr('"[^"]*\\.(csv|rds)"', lines[i])))
    for (p in gsub('"', "", paths))
      rows[[length(rows) + 1L]] <- data.frame(
        document = docname, line = i, source_file = p)
  }
  do.call(rbind, rows)
}
reads <- do.call(rbind, Map(chunk_reads, doc_lines, names(documents)))
n_inline <- vapply(doc_lines, function(l)
  sum(lengths(regmatches(l, gregexpr("`r [^`]+`", l)))), integer(1))

# ---- outputs --------------------------------------------------------------
sm_claim <- data.frame(
  document = ok$document, location = paste0("line ", ok$line),
  value = ok$value, source_file = ok$source_file,
  source_field = ok$source_field,
  script = vapply(ok$source_file, script_of, character(1)),
  commit = commit, kind = "verified_claim", note = ok$note)
sm_comp <- if (!is.null(reads)) data.frame(
  document = reads$document, location = paste0("line ", reads$line),
  value = "(computed at render)",
  source_file = sub("^\\.\\./", "", reads$source_file),
  source_field = "(chunk reads this file)",
  script = vapply(sub("^\\.\\./", "", reads$source_file), script_of,
                  character(1)),
  commit = commit, kind = "computed_chunk", note = "") else NULL
source_map <- rbind(sm_claim, sm_comp)
write.csv(source_map, "revision/SOURCE_MAP.csv", row.names = FALSE)

unm <- toks[toks$class == "candidate" & !toks$covered,
            c("document", "line", "token", "context")]
names(unm) <- c("document", "location", "value", "context")
unm$location <- paste0("line ", unm$location)
write.csv(unm, "revision/UNMATCHED_NUMBERS.csv", row.names = FALSE)

decl <- toks[toks$class == "declared_parameter" & !toks$covered, ]

cat("== build_source_map ==\n")
cat(sprintf("Claims: %d total; %d verified; %d failed; %d anchor missing; %d extract errors\n",
            nrow(cv), sum(cv$status == "OK"), sum(cv$status == "FAILED"),
            sum(cv$status == "ANCHOR_MISSING"), sum(cv$status == "EXTRACT_ERROR")))
bad <- cv[cv$status != "OK", ]
if (nrow(bad)) {
  cat("\n-- claims needing attention --\n")
  print(bad[, c("idx", "document", "value", "status", "message")], row.names = FALSE)
}
cat(sprintf("\nInline `r` expressions (computed at render): manuscript %d, presentation %d\n",
            n_inline[1], n_inline[2]))
cat(sprintf("Chunk-level file-read rows recorded: %d\n",
            if (is.null(sm_comp)) 0L else nrow(sm_comp)))
cat("\nToken census by class and coverage:\n")
print(table(toks$document, toks$class))
cat(sprintf("\nUNMATCHED literals (candidate class, no verified claim): %d (manuscript %d, presentation %d)\n",
            nrow(unm), sum(unm$document == "manuscript"),
            sum(unm$document == "presentation")))
cat(sprintf("Declared parameters without a claim (reported, not failing): %d\n",
            nrow(decl)))
cat("\nWrote revision/SOURCE_MAP.csv (", nrow(source_map), " rows) and ",
    "revision/UNMATCHED_NUMBERS.csv (", nrow(unm), " rows)\n", sep = "")

n_claim_fail <- sum(cv$status != "OK")
if (n_claim_fail > 0) {
  cat("\nbuild_source_map: FAILING (", n_claim_fail, " claim problem(s))\n", sep = "")
  quit(status = 1)
}
if (strict && nrow(unm) > 0) {
  cat("\nbuild_source_map --strict: FAILING (", nrow(unm),
      " unmatched literal(s))\n", sep = "")
  quit(status = 1)
}
invisible(NULL)

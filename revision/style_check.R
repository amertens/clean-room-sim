# style_check.R -- prose-style gate for the cleanTMLE revision (WP0).
#
# Scans the manuscript, vignettes, presentation, package prose (roxygen,
# DESCRIPTION, NEWS), and the revision documents for the banned style
# tokens of the revision protocol, and reports counts by file.
#
# Failure-class tokens (violations; --strict exits 1 when any remain
# after the allowlist): em-dashes; contractions; the banned vocabulary
# "leverage", "delve", "crucial", "landscape", "notably", "importantly",
# "highlight", "underscore" (as a verb); and "synthetic" (global
# terminology change: "simulated" replaces "synthetic").
#
# Warning-class tokens (reported, never fatal): "robust" outside the
# technical collocations (doubly robust, robust variance, robust
# standard error, robustness); en-dashes; and the per-section counts of
# the corrective constructions "rather than" and "not X, but Y" (the
# protocol allows at most one per section).
#
# Usage:
#   Rscript revision/style_check.R                 # report on the default set
#   Rscript revision/style_check.R --strict        # exit 1 on any violation
#   Rscript revision/style_check.R --strict FILE.. # check specific files
#
# Run before every commit; a commit that touches a document must leave
# that document with zero failure-class tokens.

options(width = 200)
args <- commandArgs(trailingOnly = TRUE)
strict <- "--strict" %in% args
files_arg <- setdiff(args, "--strict")

root <- local({
  d <- getwd()
  while (!file.exists(file.path(d, "clean-room-sim.Rproj")) && dirname(d) != d)
    d <- dirname(d)
  d
})
setwd(root)

default_files <- c(
  Sys.glob("reports/*.qmd"),
  Sys.glob("cleanTMLE/vignettes/*.qmd"),
  Sys.glob("cleanTMLE/R/*.R"),
  "cleanTMLE/DESCRIPTION",
  "cleanTMLE/NEWS.md",
  Sys.glob("revision/*.md")
)
files <- if (length(files_arg)) files_arg else default_files
files <- files[file.exists(files)]

# Allowlist: exemptions with a file regex, a line regex, and optionally a
# token name the exemption is scoped to. Keep this short and justified;
# every entry is a documented exception, not a loophole.
allowlist <- list(
  list(file = "manuscript_outcome_blind_dq\\.qmd$",
       line = "underscore-prefixed scripts", token = "underscore",
       why  = "describes filenames beginning with an underscore"),
  list(file = "revision/(STAGE0_AUDIT|FINDINGS|DECISIONS_PENDING|CHANGELOG)\\.md$",
       line = ".", token = "synthetic",
       why  = "audit documents quote and count the banned token by design"),
  list(file = "revision/(STAGE0_AUDIT|FINDINGS|DECISIONS_PENDING|CHANGELOG)\\.md$",
       line = ".", token = "underscore",
       why  = "audit documents name the token and underscore-prefixed files"),
  list(file = "revision/(STAGE0_AUDIT|FINDINGS|DECISIONS_PENDING|CHANGELOG)\\.md$",
       line = ".", token = "highlight",
       why  = "audit documents quote and count the banned token by design"),
  list(file = "cleanTMLE/R/plots\\.R$",
       line = "highlight <-|\\[highlight", token = "highlight",
       why  = "local variable name in plot code, not prose")
)

contraction_re <- paste0(
  "\\b(don't|doesn't|didn't|can't|won't|wouldn't|couldn't|shouldn't|",
  "isn't|aren't|wasn't|weren't|hasn't|haven't|hadn't|it's|that's|",
  "there's|here's|what's|who's|let's|we're|we've|we'll|they're|",
  "they've|you're|you've|I'm|I've|I'll|I'd)\\b")

fail_tokens <- list(
  em_dash      = "—",
  contraction  = contraction_re,
  leverage     = "\\bleverag(e|es|ed|ing)\\b",
  delve        = "\\bdelv(e|es|ed|ing)\\b",
  crucial      = "\\bcrucial(ly)?\\b",
  landscape    = "\\blandscape(s)?\\b",
  notably      = "\\bnotably\\b",
  importantly  = "\\bimportantly\\b",
  highlight    = "\\bhighlight(s|ed|ing)?\\b",
  underscore   = "\\bunderscor(e|es|ed|ing)\\b",
  synthetic    = "\\bsynthetic(ally)?\\b"
)
warn_tokens <- list(
  en_dash      = "–",
  robust_fill  = paste0("\\brobust\\b",
                        "(?<!doubly robust)(?<!doubly-robust)"),
  rather_than  = "\\brather than\\b",
  not_x_but    = "\\bnot\\b[^.;:]{0,60}, but\\b|, not\\b"
)
robust_technical <- paste0(
  "doubly[ -]robust|robust variance|robust standard error|robustness|",
  "robust_|\"robust\"|`robust`|'robust'")

allowed_line <- function(path, line_text, token = NULL) {
  for (a in allowlist) {
    tok_ok <- is.null(a$token) || (!is.null(token) && identical(a$token, token))
    if (tok_ok && grepl(a$file, path) && grepl(a$line, line_text))
      return(TRUE)
  }
  FALSE
}

count_hits <- function(lines, re, perl = TRUE) {
  m <- gregexpr(re, lines, ignore.case = TRUE, perl = perl)
  vapply(m, function(x) if (x[1] == -1) 0L else length(x), integer(1))
}

fail_rows <- list(); warn_rows <- list(); detail <- list()
for (f in files) {
  lines <- readLines(f, warn = FALSE, encoding = "UTF-8")
  fr <- c(file = f)
  for (nm in names(fail_tokens)) {
    hits <- count_hits(lines, fail_tokens[[nm]])
    excl <- vapply(seq_along(lines), function(i)
      hits[i] > 0 && allowed_line(f, lines[i], nm), logical(1))
    hits[excl] <- 0L
    fr[nm] <- sum(hits)
    if (sum(hits) > 0)
      detail[[length(detail) + 1L]] <- data.frame(
        file = f, token = nm, line = which(hits > 0),
        text = substr(lines[hits > 0], 1, 100))
  }
  fail_rows[[f]] <- as.data.frame(as.list(fr), check.names = FALSE)

  # robust outside technical collocations
  rb <- 0L
  for (i in seq_along(lines)) {
    l <- lines[i]
    n_all  <- count_hits(l, "\\brobust(ly)?\\b")
    n_tech <- count_hits(l, robust_technical)
    rb <- rb + max(0L, n_all - n_tech)
  }
  wr <- c(file = f,
          en_dash = sum(count_hits(lines, warn_tokens$en_dash)),
          robust_nontechnical = rb,
          rather_than = sum(count_hits(lines, warn_tokens$rather_than)))
  warn_rows[[f]] <- as.data.frame(as.list(wr), check.names = FALSE)
}

fail_tab <- do.call(rbind, fail_rows); rownames(fail_tab) <- NULL
warn_tab <- do.call(rbind, warn_rows); rownames(warn_tab) <- NULL
num_cols <- setdiff(names(fail_tab), "file")
fail_tab[num_cols] <- lapply(fail_tab[num_cols], as.integer)
fail_tab$TOTAL <- rowSums(fail_tab[num_cols])

cat("== style_check: failure-class token counts by file ==\n")
print(fail_tab[order(-fail_tab$TOTAL), ], row.names = FALSE)
cat("\n== warning-class counts by file ==\n")
print(warn_tab, row.names = FALSE)

# Per-section corrective-construction counts for .qmd files (warn only).
cat("\n== sections with more than one corrective construction (warn) ==\n")
any_sec <- FALSE
for (f in files[grepl("\\.qmd$", files)]) {
  lines <- readLines(f, warn = FALSE, encoding = "UTF-8")
  heads <- grep("^#{1,2} ", lines)
  if (!length(heads)) next
  bounds <- c(heads, length(lines) + 1L)
  for (k in seq_along(heads)) {
    sec <- lines[bounds[k]:(bounds[k + 1L] - 1L)]
    n_rt <- sum(count_hits(sec, warn_tokens$rather_than))
    n_nb <- sum(count_hits(sec, "\\bnot\\b[^.;:]{0,60}, but\\b"))
    if (n_rt + n_nb > 1L) {
      any_sec <- TRUE
      cat(sprintf("  %s :: %-60s rather_than=%d not_but=%d\n",
                  basename(f), substr(lines[heads[k]], 1, 60), n_rt, n_nb))
    }
  }
}
if (!any_sec) cat("  none\n")

n_fail <- sum(fail_tab$TOTAL)
if (n_fail > 0) {
  cat("\n== offending lines (failure class) ==\n")
  dt <- do.call(rbind, detail)
  print(utils::head(dt[order(dt$file, dt$line), ], 60), row.names = FALSE)
  if (nrow(dt) > 60) cat("  ...", nrow(dt) - 60, "more rows\n")
}
cat(sprintf("\nstyle_check: %d failure-class token(s) in %d file(s).\n",
            n_fail, sum(fail_tab$TOTAL > 0)))
if (strict && n_fail > 0) quit(status = 1)
invisible(NULL)

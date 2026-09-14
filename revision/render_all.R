# render_all.R -- render every revision deliverable from source (WP0).
#
# Renders, in order: the manuscript (html and docx), the two package
# vignettes (html), and the presentation (pptx), each from its Quarto
# source. Rendered documents are build outputs; prose edits go through
# the .qmd sources only (operating rule 4).
#
# Usage:
#   Rscript revision/render_all.R                    # render everything in place
#   Rscript revision/render_all.R manuscript slides  # subset of targets
#   Rscript revision/render_all.R --outdir DIR ...   # render into DIR instead
#
# Caveat for --outdir validation runs: Quarto renders beside the source
# and then relocates the output, so a tracked render at the default
# location is still displaced. Restore it afterwards with
# `git checkout -- <render>` and delete any reports/.gitignore that
# Quarto creates.
#
# Targets: manuscript, vignettes, slides. Fails loudly on the first
# render error and prints per-target wall times.

args <- commandArgs(trailingOnly = TRUE)
outdir <- NULL
if (length(w <- which(args == "--outdir"))) {
  if (length(args) < w + 1L) stop("--outdir needs a directory argument")
  outdir <- normalizePath(args[w + 1L], mustWork = FALSE)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  args <- args[-c(w, w + 1L)]
}
targets <- if (length(args)) args else c("manuscript", "vignettes", "slides")
bad <- setdiff(targets, c("manuscript", "vignettes", "slides"))
if (length(bad)) stop("Unknown target(s): ", paste(bad, collapse = ", "))

root <- local({
  d <- getwd()
  while (!file.exists(file.path(d, "clean-room-sim.Rproj")) && dirname(d) != d)
    d <- dirname(d)
  d
})
setwd(root)

quarto_bin <- "C:/Program Files/RStudio/resources/app/bin/quarto/bin/quarto.exe"
if (!file.exists(quarto_bin))
  quarto_bin <- Sys.which("quarto")
if (!nzchar(quarto_bin) || !file.exists(quarto_bin))
  stop("quarto executable not found; install Quarto or adjust quarto_bin.")

render_one <- function(qmd, to = NULL, label = qmd) {
  stopifnot(file.exists(qmd))
  args <- c("render", qmd)
  if (!is.null(to)) args <- c(args, "--to", to)
  if (!is.null(outdir)) args <- c(args, "--output-dir", outdir)
  t0 <- Sys.time()
  status <- system2(quarto_bin, shQuote(args, type = "cmd"))
  secs <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  if (status != 0)
    stop(sprintf("Render FAILED (exit %d) after %.1f s: %s", status, secs,
                 label), call. = FALSE)
  cat(sprintf("Rendered %-55s %8.1f s\n", label, secs))
  invisible(secs)
}

t_all <- Sys.time()
if ("manuscript" %in% targets) {
  render_one("reports/manuscript_outcome_blind_dq.qmd", to = "html",
             label = "manuscript (html)")
  render_one("reports/manuscript_outcome_blind_dq.qmd", to = "docx",
             label = "manuscript (docx)")
}
if ("vignettes" %in% targets) {
  render_one("cleanTMLE/vignettes/cleanTMLE-staged-analysis.qmd", to = "html",
             label = "vignette: staged analysis")
  render_one("cleanTMLE/vignettes/cleanTMLE-functions.qmd", to = "html",
             label = "vignette: functions")
}
if ("slides" %in% targets) {
  render_one("reports/cleanTMLE_presentation.qmd", to = "pptx",
             label = "presentation (pptx)")
}
cat(sprintf("\nrender_all: done in %.1f min%s\n",
            as.numeric(difftime(Sys.time(), t_all, units = "mins")),
            if (is.null(outdir)) "" else paste0(" (outputs in ", outdir, ")")))

# Generate static PNG flowcharts to replace mermaid diagrams (which do not
# render in docx output) for the manuscript appendix. Base R graphics only:
#   Rscript make_figures.R
# Produces reports/figures/fig-roadmap.png and fig-governance.png, the two
# appendix figures the manuscript includes.

# Resolve an output directory next to this script (reports/figures), so the
# figures land in the right place regardless of the working directory.
.args <- commandArgs(FALSE)
.sf <- sub("^--file=", "", .args[grep("^--file=", .args)])
.base <- if (length(.sf)) dirname(normalizePath(.sf)) else getwd()
figdir <- file.path(.base, "figures")
dir.create(figdir, showWarnings = FALSE)

## ---- helpers ---------------------------------------------------------------

# Draw a rounded-ish rectangle (plain rect) with centred, wrapped text.
draw_box <- function(xc, yc, w, h, label, fill = "white",
                     border = "grey30", cex = 1, wrap = 28) {
  rect(xc - w/2, yc - h/2, xc + w/2, yc + h/2,
       col = fill, border = border, lwd = 1.4)
  lines <- strwrap(label, width = wrap)
  n <- length(lines)
  ys <- yc + (rev(seq_len(n)) - (n + 1)/2) * strheight("Ag", cex = cex) * 1.35
  text(xc, ys, lines, cex = cex)
}

varrow <- function(x, y0, y1) arrows(x, y0, x, y1, length = 0.10, lwd = 1.4, col = "grey30")
harrow <- function(x0, x1, y) arrows(x0, y, x1, y, length = 0.10, lwd = 1.4, col = "grey30")

## ---- 1. cleanTMLE inside the causal roadmap --------------------------------

roadmap <- list(
  list("1. Causal question and target trial", "white"),
  list("2. Observed data and source characterisation", "white"),
  list("3. Identification assumptions", "white"),
  list("4. Statistical estimand", "white"),
  list("5. Lock candidate estimators and thresholds (create_analysis_lock, lock_primary_tmle_spec)", "grey85"),
  list("Pre-outcome checks (checkpoints, plasmode, DQ stress)", "grey85"),
  list("Pre-outcome decision (gate_check, authorize_outcome_analysis)", "grey85"),
  list("Unblind and run locked primary analysis (four-step TMLE, IPCW-TMLE)", "grey85"),
  list("At-outcome diagnostics: negative controls, balance", "white"),
  list("Post-outcome sensitivity: E-value, QBA, causal-gap", "white"),
  list("Interpretation and transportability", "white")
)
png(file.path(figdir, "fig-roadmap.png"), width = 1500, height = 2400, res = 200)
op <- par(mar = c(0.5, 0.5, 0.5, 0.5)); on.exit(par(op), add = TRUE)
n <- length(roadmap)
plot.new(); plot.window(xlim = c(0, 10), ylim = c(0, n + 0.5))
bw <- 9; bh <- 0.72; xc <- 5
for (i in seq_len(n)) {
  yc <- n - i + 0.5
  draw_box(xc, yc, bw, bh, roadmap[[i]][[1]], fill = roadmap[[i]][[2]], cex = 0.85, wrap = 46)
  if (i < n) varrow(xc, yc - bh/2, (n - i - 0.5) + bh/2)
}
invisible(dev.off())

## ---- 2. Governance layers --------------------------------------------------

png(file.path(figdir, "fig-governance.png"), width = 1700, height = 1500, res = 200)
op <- par(mar = c(0.5, 0.5, 0.5, 0.5)); on.exit(par(op), add = TRUE)
plot.new(); plot.window(xlim = c(0, 12), ylim = c(0, 9))
main <- c("Protocol and estimand", "Analysis-lock object", "Audit trail",
          "Decision log", "Pre-outcome decision point", "Primary analysis")
xc <- 8.5; bw <- 6; bh <- 0.8
for (i in seq_along(main)) {
  yc <- 9 - i*1.3 + 0.2
  draw_box(xc, yc, bw, bh, main[i], fill = "white", cex = 0.9, wrap = 32)
  assign(paste0("y", i), yc)
}
for (i in 1:5) varrow(xc, get(paste0("y", i)) - bh/2, get(paste0("y", i+1)) + bh/2)
# External inputs feeding the pre-outcome decision point (box 5)
ext <- c("External validation studies", "Independent review and role governance",
         "Data standards and provenance")
exc <- 2.4; ebw <- 3.6
ey <- c(y5 + 1.5, y5, y5 - 1.5)
for (j in seq_along(ext)) {
  draw_box(exc, ey[j], ebw, 0.95, ext[j], fill = "grey92", cex = 0.82, wrap = 22)
  harrow(exc + ebw/2, xc - bw/2, ey[j])
}
invisible(dev.off())

cat("Wrote figures/fig-roadmap.png and fig-governance.png
")

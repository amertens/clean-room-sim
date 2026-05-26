# Generate static PNG flowcharts to replace mermaid diagrams (which do not
# render in docx output). Uses only base R graphics. Run from reports/:
#   Rscript make_figures.R
# Produces reports/figures/fig-roadmap.png, fig-governance.png, fig-cohort-flow.png

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

## ---- 3. Rescue.Co cohort flow ---------------------------------------------

png(file.path(figdir, "fig-cohort-flow.png"), width = 1600, height = 2000, res = 200)
op <- par(mar = c(0.5, 0.5, 0.5, 0.5)); on.exit(par(op), add = TRUE)
plot.new(); plot.window(xlim = c(0, 12), ylim = c(0, 12))
flow <- c("Registry observations (all 2018-2024 trauma transports)",
          "Ambulance arrivals at participating hospitals",
          "Direct-from-scene cohort (inter-facility transfers excluded)",
          "Analytic cohort: n = 1,693 (1,013 Rescue.Co + 680 other)",
          "Complete 6-month GOSE: n = 1,277",
          "Matched subset: n = 615 (matching) / 1,230 (matched-TMLE)")
xc <- 4.5; bw <- 7.5; bh <- 0.9
ys <- numeric(length(flow))
for (i in seq_along(flow)) {
  yc <- 12 - i*1.7 + 0.3; ys[i] <- yc
  draw_box(xc, yc, bw, bh, flow[i], fill = "white", cex = 0.85, wrap = 40)
  if (i > 1) varrow(xc, ys[i-1] - bh/2, yc + bh/2)
}
# Side note: missing GOSE handled by IPCW, branching off the analytic cohort (box 4)
draw_box(10.3, ys[4], 3, 1.0, "Missing GOSE n = 416 (24.6%): handled by IPCW",
         fill = "grey92", cex = 0.78, wrap = 20)
arrows(xc + bw/2, ys[4], 10.3 - 1.5, ys[4], length = 0.10, lwd = 1.2,
       lty = 2, col = "grey40")
invisible(dev.off())

## ---- 4. Negative-control attrition ---------------------------------------

png(file.path(figdir, "fig-nc-attrition.png"), width = 1600, height = 1500, res = 200)
op <- par(mar = c(0.5, 0.5, 0.5, 0.5)); on.exit(par(op), add = TRUE)
plot.new(); plot.window(xlim = c(0, 12), ylim = c(0, 10))
xc <- 4.5; bw <- 7.5
yA <- 9; yB <- 6.7; yD <- 4.2; yE <- 1.9
draw_box(xc, yA, bw, 1.4, "5 prespecified NCs: chronic_hypertension, chronic_diabetes_insulin, chronic_hiv_art, household_urban, fuel_wood", cex = 0.82, wrap = 40)
draw_box(xc, yB, bw, 0.8, "NZV filter at data-sanitisation step", cex = 0.85, wrap = 40)
draw_box(xc, yD, bw, 1.0, "2 analysed: chronic_hypertension, household_urban", cex = 0.85, wrap = 40)
draw_box(xc, yE, bw, 0.8, "Both pass NC checks", cex = 0.85, wrap = 40)
varrow(xc, yA - 0.7, yB + 0.4)
varrow(xc, yB - 0.4, yD + 0.5)
varrow(xc, yD - 0.5, yE + 0.4)
# Side branch: 3 dropped
draw_box(10.3, yB, 3, 1.2, "3 dropped: chronic_diabetes_insulin, chronic_hiv_art, fuel_wood",
         fill = "grey92", cex = 0.78, wrap = 20)
arrows(xc + bw/2, yB, 10.3 - 1.5, yB, length = 0.10, lwd = 1.2, lty = 2, col = "grey40")
invisible(dev.off())

cat("Wrote figures/fig-roadmap.png, fig-governance.png, fig-cohort-flow.png, fig-nc-attrition.png\n")

# Guards for the steps after the design stage (scripts 22b, 22c, 23, 24 and 25):
# the build fingerprint, the approvals file and the check that it is pushed, a
# log of their own, and the rows of a primary comparison before the
# complete-case restriction. Decisions of 2026-10-23 on report recommendation 1.

REG_DIR <- "rescueCo/registration"   # tracked: build, hashes, approvals
RX_LOG <- "rescueCo/logs/reexec_050_after_design.log"
FP_PKGS <- c("cleanTMLE", "tmle", "glmnet", "xgboost", "earth", "nnet",
             "SuperLearner", "lmtp", "generics", "digest", "future", "future.apply")

# cr_log() appends to rescueCo/logs/pipeline.log, whose first lines hold
# outcome counts from the unblinded analysis and which no one opens (ledger,
# section 5), so these scripts log here instead.
rx_log <- function(msg) {
  entry <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", msg)
  message(entry)
  cat(entry, "\n", file = RX_LOG, append = TRUE)
}

rx_md5 <- function(files) unname(tools::md5sum(files))

# One md5 for an installed package folder, from every file's path and md5.
.rx_folder_md5 <- function(path) {
  f <- sort(list.files(path, recursive = TRUE, full.names = TRUE))
  rel <- substring(f, nchar(path) + 2L)
  digest::digest(paste(rel, tools::md5sum(f), collapse = "\n"), algo = "md5", serialize = FALSE)
}

# R's version, and for each package that builds, fits or targets the models
# its version, its library and the md5 of its installed folder. Run after
# .libPaths(c("rescueCo/lib", .libPaths())), as every script here does.
rx_fingerprint <- function() {
  rows <- lapply(FP_PKGS, function(p) {
    path <- find.package(p, quiet = TRUE)
    if (!length(path))
      return(data.frame(package = p, version = NA, library = NA, folder_md5 = NA))
    lib <- if (grepl("rescueCo/lib/", normalizePath(path, winslash = "/"), fixed = TRUE))
      "rescueCo/lib" else "user library"
    data.frame(package = p, version = as.character(utils::packageVersion(p, lib.loc = dirname(path))),
               library = lib, folder_md5 = .rx_folder_md5(path))
  })
  rbind(data.frame(package = "R", version = R.version.string, library = NA, folder_md5 = NA),
        do.call(rbind, rows))
}

# Stops unless the installed build is the one script 22c recorded, so the
# analysis runs the fitter the simulation evaluated (report recommendations 1
# and 7; the lock itself stores only cleanTMLE's version string).
rx_check_fingerprint <- function(file = file.path(REG_DIR, "build_fingerprint.csv")) {
  if (!file.exists(file))
    stop("No build fingerprint at ", file, "; run script 22c first.", call. = FALSE)
  rec <- utils::read.csv(file, colClasses = "character")
  now <- rx_fingerprint()
  key <- function(d) stats::setNames(paste(d$version, d$library, d$folder_md5), d$package)
  a <- key(rec); b <- key(now)
  bad <- union(setdiff(names(a), names(b)), names(a)[names(a) %in% names(b) & a != b[names(a)]])
  if (length(bad))
    stop("The installed build differs from the recorded one for: ", paste(bad, collapse = ", "),
         ". Nothing is estimated on a build the design stage did not evaluate.", call. = FALSE)
  invisible(TRUE)
}

# The approvals file and everything else in rescueCo/registration must be
# committed and on a remote branch before an outcome is scored (decision 4),
# so a dated copy of the hashes exists off this machine.
rx_check_pushed <- function(paths = REG_DIR) {
  git <- Sys.which("git")
  if (!nzchar(git))
    stop("git is not on the PATH, so the check that the approvals are pushed cannot run.", call. = FALSE)
  st <- system2(git, c("status", "--porcelain", "--", paths), stdout = TRUE)
  if (length(st))
    stop(paths, " has uncommitted changes:\n", paste(st, collapse = "\n"), call. = FALSE)
  head <- system2(git, c("log", "-1", "--format=%H", "--", paths), stdout = TRUE)
  if (!length(head) || !nzchar(head[1]))
    stop(paths, " is not committed.", call. = FALSE)
  rem <- system2(git, c("branch", "-r", "--contains", head[1]), stdout = TRUE)
  if (!length(rem))
    stop("The commit holding the approvals (", substr(head[1], 1, 12),
         ") is on no remote branch; push it before the outcome store is built.", call. = FALSE)
  invisible(head[1])
}

# The approvals file: one row per final lock, written by script 23 (final) as
# a template and completed by Andrew, who fills approved_by and approved_on
# for each design report he approves (decision 3). Returns the rows to
# unblind, after checking each hash against the design report on disk.
rx_read_approvals <- function(lockd, file = file.path(REG_DIR, "approvals_final.csv"),
                              require_pushed = TRUE) {
  if (!file.exists(file))
    stop("No approvals at ", file, ". Complete the template from script 23 and commit it.",
         call. = FALSE)
  ap <- utils::read.csv(file, colClasses = "character")
  need <- c("lock", "role", "lock_hash", "dossier_hash", "approved_by", "approved_on")
  miss <- setdiff(need, names(ap))
  if (length(miss)) stop("The approvals file lacks: ", paste(miss, collapse = ", "), ".", call. = FALSE)
  run <- ap[ap$role != "not run", , drop = FALSE]
  blank <- is.na(run$approved_by) | !nzchar(trimws(run$approved_by)) |
    is.na(suppressWarnings(as.Date(run$approved_on)))
  if (any(blank))
    stop("Every design report to be unblinded needs approved_by and a date in approved_on; missing for: ",
         paste(run$lock[blank], collapse = ", "), ".", call. = FALSE)
  for (i in seq_len(nrow(run))) {
    f <- file.path(lockd, paste0(run$lock[i], "_dossier.rds"))
    if (!file.exists(f)) stop("No design report for approved lock ", run$lock[i], ".", call. = FALSE)
    d <- readRDS(f)
    if (!identical(d$dossier_hash, run$dossier_hash[i]) || !identical(d$lock_hash, run$lock_hash[i]))
      stop("The design report for ", run$lock[i], " does not match the approved hashes.", call. = FALSE)
  }
  if (require_pushed) rx_check_pushed()
  run
}

# The rows of a primary comparison before the complete-case restriction, as
# script 22 builds them (its LOCKS, L1 to L5, and lock_rows()). Script 23 checks
# that the complete-case part of these rows is exactly the lock's patients, so
# this copy cannot drift from script 22 unnoticed.
RX_PRIMARY <- list(
  L1 = list(treated = "rescueco", control = c("other_operator", "other_means"), direct = FALSE),
  L2 = list(treated = "rescueco", control = c("other_operator", "other_means"), direct = TRUE),
  L3 = list(treated = "rescueco", control = "other_operator", direct = TRUE),
  L4 = list(treated = "rescueco", control = "other_means", direct = TRUE),
  L5 = list(treated = c("rescueco", "other_operator"), control = "other_means", direct = TRUE))

rx_primary_rows <- function(D, nm) {
  s <- RX_PRIMARY[[nm]]
  keep <- D$arm %in% c(s$treated, s$control) & !(D$crossref_tier %in% 5L)
  if (isTRUE(s$direct)) keep <- keep & D$direct_arrival == 1
  keep
}

# Binds data frames whose columns differ, filling the missing ones with NA.
rx_bind_fill <- function(l) {
  cols <- unique(unlist(lapply(l, names)))
  do.call(rbind, lapply(l, function(x) { x[setdiff(cols, names(x))] <- NA; x[cols] }))
}

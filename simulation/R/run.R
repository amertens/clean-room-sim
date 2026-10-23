# Running and caching. A cell is (stage, unit, n, idx). The unit of the map
# stages is an overlap level (a row of overlap_levels()): map runs never see the
# outcome, so one run serves every outcome model and effect at that level. The
# unit of the key stage is a scenario (a row of scenario_table()). Each cell's
# rows are saved to cache/<stage>/<unit id>_n<n>_i<idx>.rds together with the
# key that determined them: the stage, the unit row, n, idx, reps (map stages),
# the cleanTMLE version and the text of the code files the stage runs. A cached
# file is used only when its key is identical, so a hit equals recomputation.

STAGES <- c("map_generate", "map_sample", "key")

sim_dir <- function() getOption("sim.dir", ".")

stage_files <- function(stage)
  file.path(sim_dir(), "R", c("dgp.R", if (stage == "key") c("map.R", "key.R") else "map.R"))

code_hash <- function(stage) {
  txt <- vapply(stage_files(stage), function(f) gsub("\r", "", paste(readLines(f, warn = FALSE),
                                                                    collapse = "\n")),
                character(1))
  digest::digest(unname(txt), algo = "sha256")
}

cell_key <- function(stage, unit, n, idx, reps) {
  list(stage = stage, unit = as.list(unit), n = as.integer(n), idx = as.integer(idx),
       reps = if (stage == "key") NA_integer_ else as.integer(reps),
       cleanTMLE = as.character(utils::packageVersion("cleanTMLE")),
       code = code_hash(stage))
}

cell_path <- function(stage, unit, n, idx, cache_dir)
  file.path(cache_dir, stage, sprintf("%s_n%d_i%d.rds", unit$id, as.integer(n),
                                      as.integer(idx)))

run_cell <- function(stage, unit, n, idx, reps, cache_dir) {
  path <- cell_path(stage, unit, n, idx, cache_dir)
  key <- cell_key(stage, unit, n, idx, reps)
  if (file.exists(path)) {
    obj <- readRDS(path)
    if (identical(obj$key, key)) return(obj$rows)
  }
  rows <- switch(stage,
    map_generate = run_map(unit, n, idx, "generate", reps),
    map_sample = run_map(unit, n, idx, "sample", reps),
    key = run_key(unit, n, idx))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp", Sys.getpid())
  saveRDS(list(key = key, rows = rows, metrics = attr(rows, "metrics")), tmp)
  file.rename(tmp, path)  # atomic: a killed run leaves no truncated entry
  rows
}

# Runs every (unit, n, idx) of `stage` not already cached. `units` is a subset
# of rows of overlap_levels() (map stages) or of scenario_table() (key stage).
# Sequential by default; workers > 1 uses future.apply with that many
# background R sessions, each started with simulation/lib first on its library
# path (so it loads the project's cleanTMLE) and sourcing the code.
run_stage <- function(stage, units, ns, idx, reps = 200L, workers = 1L,
                      cache_dir = file.path(sim_dir(), "cache")) {
  stage <- match.arg(stage, STAGES)
  # Force every argument here: a promise sent to a worker would be evaluated
  # there, where the caller's variables do not exist.
  force(units); force(reps); force(cache_dir)
  cells <- expand.grid(i = seq_len(nrow(units)), n = ns, idx = idx)
  if (workers <= 1L)
    return(do.call(rbind, lapply(seq_len(nrow(cells)), function(j)
      run_cell(stage, units[cells$i[j], ], cells$n[j], cells$idx[j], reps, cache_dir))))
  dir <- sim_dir()
  lib <- normalizePath(file.path(dir, "lib"), "/")
  oplan <- future::plan(future::multisession, workers = workers,
                        rscript_libs = c(lib, "*"))
  on.exit(future::plan(oplan), add = TRUE)
  rows <- future.apply::future_lapply(seq_len(nrow(cells)), function(j) {
    .libPaths(unique(c(lib, .libPaths())))
    options(sim.dir = dir)
    for (f in c("dgp.R", "map.R", "key.R", "run.R")) source(file.path(dir, "R", f))
    if (dirname(find.package("cleanTMLE")) != lib)
      stop("Worker loaded cleanTMLE from ", find.package("cleanTMLE"))
    run_cell(stage, units[cells$i[j], ], cells$n[j], cells$idx[j], reps, cache_dir)
  }, future.seed = NULL, future.chunk.size = 1L, future.globals = c("stage", "units", "cells", "reps", "cache_dir",
                                             "dir", "lib"),
     future.packages = character())
  do.call(rbind, rows)
}

# All cached rows of a stage in one data frame (what = "rows" or "metrics").
collect <- function(stage, what = c("rows", "metrics"),
                    cache_dir = file.path(sim_dir(), "cache")) {
  what <- match.arg(what)
  fs <- list.files(file.path(cache_dir, stage), pattern = "\\.rds$", full.names = TRUE)
  out <- do.call(rbind, lapply(fs, function(f) readRDS(f)[[what]]))
  rownames(out) <- NULL
  out
}

# Running and caching. A cell is (stage, scenario, n, idx). Each cell's rows are
# saved to cache/<stage>/<scenario>_n<n>_i<idx>.rds together with the key that
# determined them: the stage, the scenario row, n, idx, reps (map stages), the
# cleanTMLE version and the text of the code files the stage runs. A cached
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

cell_key <- function(stage, scenario, n, idx, reps) {
  list(stage = stage, scenario = as.list(scenario), n = as.integer(n), idx = as.integer(idx),
       reps = if (stage == "key") NA_integer_ else as.integer(reps),
       cleanTMLE = as.character(utils::packageVersion("cleanTMLE")),
       code = code_hash(stage))
}

cell_path <- function(stage, scenario, n, idx, cache_dir)
  file.path(cache_dir, stage, sprintf("%s_n%d_i%d.rds", scenario$id, as.integer(n),
                                      as.integer(idx)))

run_cell <- function(stage, scenario, n, idx, reps, cache_dir) {
  path <- cell_path(stage, scenario, n, idx, cache_dir)
  key <- cell_key(stage, scenario, n, idx, reps)
  if (file.exists(path)) {
    obj <- readRDS(path)
    if (identical(obj$key, key)) return(obj$rows)
  }
  rows <- switch(stage,
    map_generate = run_map(scenario, n, idx, "generate", reps),
    map_sample = run_map(scenario, n, idx, "sample", reps),
    key = run_key(scenario, n, idx))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp", Sys.getpid())
  saveRDS(list(key = key, rows = rows, metrics = attr(rows, "metrics")), tmp)
  file.rename(tmp, path)  # atomic: a killed run leaves no truncated entry
  rows
}

# Runs every (scenario, n, idx) of `stage` not already cached. `scenarios` is
# a subset of rows of the scenario table. Sequential by default; workers > 1
# uses future.apply with that many background R sessions (each sources the
# code and puts simulation/lib first on its library path).
run_stage <- function(stage, scenarios, ns, idx, reps = 200L, workers = 1L,
                      cache_dir = file.path(sim_dir(), "cache")) {
  stage <- match.arg(stage, STAGES)
  cells <- expand.grid(i = seq_len(nrow(scenarios)), n = ns, idx = idx)
  one <- function(j) {
    c <- cells[j, ]
    run_cell(stage, scenarios[c$i, ], c$n, c$idx, reps, cache_dir)
  }
  rows <- if (workers > 1L) {
    oplan <- future::plan(future::multisession, workers = workers)
    on.exit(future::plan(oplan), add = TRUE)
    dir <- sim_dir()
    future.apply::future_lapply(seq_len(nrow(cells)), function(j) {
      .libPaths(c(file.path(dir, "lib"), .libPaths()))
      options(sim.dir = dir)
      for (f in c("dgp.R", "map.R", "key.R", "run.R")) source(file.path(dir, "R", f))
      one(j)
    }, future.seed = NULL)
  } else lapply(seq_len(nrow(cells)), one)
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

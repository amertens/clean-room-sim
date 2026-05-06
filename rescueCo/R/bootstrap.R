# bootstrap.R — auto-detect the layout so scripts work whether invoked
# from the project root (rescueCo/ as child) OR from inside rescueCo/.
#
# Usage at the top of any script:
#   source(file.path(if (dir.exists("clean_room")) "rescueCo/R" else "R",
#                    "bootstrap.R"))
#   cr_source("utils.R")
#   cr_source("ps_matching.R")
#
# After cr_bootstrap() runs, the working directory is *always* the parent
# of rescueCo/ (i.e., the "project root"), so paths like
# "rescueCo/results/..." are unambiguous.

cr_bootstrap <- function() {
  # Already at the right place
  if (dir.exists("rescueCo") &&
      file.exists("rescueCo/config/clean_room_config.yml")) {
    return(invisible(getwd()))
  }
  # We're inside rescueCo/ — go up one level
  if (file.exists("config/clean_room_config.yml") &&
      dir.exists("R")) {
    setwd("..")
    return(invisible(getwd()))
  }
  # Walk up the tree to find a rescueCo/ folder
  for (depth in 1:5) {
    candidate <- normalizePath(do.call(file.path, as.list(c(getwd(), rep("..", depth)))))
    if (dir.exists(file.path(candidate, "rescueCo"))) {
      setwd(candidate)
      return(invisible(getwd()))
    }
  }
  stop("Cannot locate rescueCo/ — run from project root or inside rescueCo/.")
}

# Source any helper file from rescueCo/R/ regardless of cwd
cr_source <- function(file) {
  for (p in c(file.path("rescueCo/R", file),
              file.path("R", file))) {
    if (file.exists(p)) {
      source(p, local = parent.frame())
      return(invisible(p))
    }
  }
  stop("Cannot locate ", file, " in rescueCo/R/")
}

# Auto-run the layout detection on source
cr_bootstrap()

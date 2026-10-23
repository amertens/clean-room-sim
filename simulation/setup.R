# Install cleanTMLE into the project library simulation/lib (git-ignored), so the
# simulation does not depend on whichever cleanTMLE is in the user library.
# The source is the cleanTMLE submodule's `main` branch. When the submodule's
# working tree is on another branch, export main first, for example
#   git -C cleanTMLE archive main --prefix=cleanTMLE/ | tar -x -C <dir>
# and pass <dir>/cleanTMLE as the first argument.
# Run from the clean-room-sim folder: Rscript simulation/setup.R [source]

args <- commandArgs(trailingOnly = TRUE)
src <- if (length(args)) args[1] else "cleanTMLE"
lib <- "simulation/lib"
dir.create(lib, showWarnings = FALSE)
install.packages(src, repos = NULL, type = "source", lib = lib)
.libPaths(c(lib, .libPaths()))
cat("cleanTMLE", as.character(utils::packageVersion("cleanTMLE", lib.loc = lib)),
    "installed in", normalizePath(lib), "\n")

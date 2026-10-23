# ============================================================
# 22c. Record the build the design stage ran on
# ============================================================
# Writes rescueCo/registration/build_fingerprint.csv: R's version and, for
# cleanTMLE (from rescueCo/lib) and each package that fits or targets the
# models, its version, its library and an md5 of its installed folder. Scripts
# 22b and 25 stop if the installed build differs. Run once, as soon as the
# design stage allows a short R job and before anything is installed into
# rescueCo/lib, then commit the file. Reads no data.
#
# The fingerprint shows that nothing changed after it was recorded; that the
# locks already built ran on the same build rests on the lock's version string
# and on nothing having been installed into rescueCo/lib since 20 October.
#
# Usage: Rscript rescueCo/scripts/22c_record_build.R
# Decision of 2026-10-23 (report recommendations 1 and 7).
# ============================================================

source(file.path(if (dir.exists("rescueCo")) "rescueCo/R" else "R", "bootstrap.R"))
source("rescueCo/R/utils.R")
source("rescueCo/R/reexec_guards.R")
.libPaths(c("rescueCo/lib", .libPaths()))

fp <- rx_fingerprint()
ct <- fp[fp$package == "cleanTMLE", ]
if (is.na(ct$version) || !identical(ct$library, "rescueCo/lib"))
  stop("cleanTMLE was not found in rescueCo/lib.", call. = FALSE)
f <- file.path(REG_DIR, "build_fingerprint.csv")
if (file.exists(f))
  stop("A fingerprint is already recorded at ", f, "; replacing it needs a ledger entry first.",
       call. = FALSE)
dir.create(REG_DIR, showWarnings = FALSE)
utils::write.csv(fp, f, row.names = FALSE)
rx_log(paste("22c: recorded the build; cleanTMLE", ct$version, "folder md5", ct$folder_md5,
             "; commit", f))

# Entry point for testthat
library(testthat)

# Source all R modules
r_files <- list.files("R", pattern = "\\.R$", recursive = TRUE,
                      full.names = TRUE)
for (f in r_files) source(f, local = FALSE)

test_dir("tests/testthat")

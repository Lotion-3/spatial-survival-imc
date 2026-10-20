# Run the R test files; exits non-zero on any failure.
#   Rscript tests/R/run_tests.R          (from the repository root)
root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "..", ".."))
Sys.setenv(SPATIALSURV_ROOT = root)
res <- testthat::test_dir(file.path(root, "tests", "R"), reporter = "summary", stop_on_failure = TRUE)

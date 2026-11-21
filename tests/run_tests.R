# Run the test suite from the repository root:
#   Rscript tests/run_tests.R
# Needs the app's R packages (shiny, yaml, QDNAseq, DNAcopy) plus testthat.
# No external tools, model files or sequencing data are needed.
library(testthat)
args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grepl("^--file=", args)])
tests_dir <- if (length(script) == 1) file.path(dirname(normalizePath(script)), "testthat") else "tests/testthat"
test_dir(tests_dir, reporter = ProgressReporter$new(show_praise = FALSE), stop_on_failure = TRUE)

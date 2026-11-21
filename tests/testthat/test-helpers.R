test_that("pad_barcode zero-pads single digits only", {
  expect_equal(pad_barcode(1), "01")
  expect_equal(pad_barcode("7"), "07")
  expect_equal(pad_barcode("12"), "12")
})

test_that("apply_extra_path prepends directories and ignores blanks", {
  old <- Sys.getenv("PATH")
  on.exit(Sys.setenv(PATH = old))
  apply_extra_path(c("", "/opt/a/bin", "/opt/b/bin"))
  expect_true(startsWith(Sys.getenv("PATH"), "/opt/a/bin:/opt/b/bin:"))
  apply_extra_path(character(0))
  expect_true(startsWith(Sys.getenv("PATH"), "/opt/a/bin:/opt/b/bin:"))
})

test_that("run_cmd returns output and raises an error with stderr on failure", {
  quiet <- function(...) invisible(NULL)
  expect_equal(run_cmd("echo hello", log_fun = quiet), "hello")
  # system() also warns about the non-zero exit; only the error matters here
  err <- expect_error(suppressWarnings(run_cmd("sh -c 'echo boom >&2; exit 3'", log_fun = quiet)), "exit 3")
  expect_match(conditionMessage(err), "boom")
})

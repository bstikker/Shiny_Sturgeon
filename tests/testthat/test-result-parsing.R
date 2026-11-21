test_that("v1 result CSV parses into class/score, without the number_probes column", {
  res <- parse_sturgeon_result(fixture("sturgeon_v1_result.csv"))
  expect_equal(colnames(res), c("class", "score"))
  expect_equal(nrow(res), 5)
  expect_false("number_probes" %in% res$class)
  expect_equal(as.numeric(res$score[res$class == "Glioblastoma - GBM - RTK II"]), 0.91)
})

test_that("earlier sturgeon-v2 output (comma-separated, probes last) parses", {
  res <- parse_sturgeon_v2_result(fixture("sturgeon_v2_result_legacy.csv"))
  expect_equal(nrow(res), 5)
  expect_false("probes" %in% res$class)
  expect_equal(as.numeric(res$score[res$class == "GBM_RTK2"]), 0.88)
})

test_that("public sturgeon-v2 result.txt (tab-separated, one row per sampling) is averaged", {
  res <- parse_sturgeon_v2_result(fixture("sturgeon_v2_result.txt"))
  expect_equal(nrow(res), 5)
  expect_false("probes" %in% res$class)
  expect_equal(as.numeric(res$score[res$class == "GBM_RTK2"]), mean(c(0.9, 0.8, 0.7)))
})

test_that("v1 results with more than one row are rejected", {
  f <- tempfile(fileext = ".csv")
  writeLines(c("number_probes,A,B", "10,0.5,0.5", "10,0.4,0.6"), f)
  expect_error(parse_sturgeon_result(f), "expected exactly one")
})

test_that("a result with only the metadata column is rejected", {
  f <- tempfile(fileext = ".csv")
  writeLines(c("number_probes", "10"), f)
  expect_error(parse_sturgeon_result(f), "no class columns")
})

test_that("find_single_file finds exactly one match and errors otherwise", {
  d <- tempfile("pred_"); dir.create(d)
  expect_error(find_single_file(d, "\\.csv$"), "No file found")
  file.create(file.path(d, "merged_probes_methyl_calls_general.csv"))
  expect_equal(basename(find_single_file(d, "^merged_probes_methyl_calls_.*\\.csv$")),
               "merged_probes_methyl_calls_general.csv")
  file.create(file.path(d, "merged_probes_methyl_calls_other.csv"))
  expect_error(find_single_file(d, "^merged_probes_methyl_calls_.*\\.csv$"), "exactly one")
})

test_that("resolve_sturgeon_v2_output handles both a result file and a result directory", {
  f <- tempfile(fileext = ".csv"); file.create(f)
  expect_equal(resolve_sturgeon_v2_output(f), f)

  d <- tempfile("v2_out_"); dir.create(d)
  expect_error(resolve_sturgeon_v2_output(d), "without result.txt")
  file.create(file.path(d, "result.txt"))
  expect_equal(resolve_sturgeon_v2_output(d), file.path(d, "result.txt"))

  expect_error(resolve_sturgeon_v2_output(tempfile()), "did not produce")
})

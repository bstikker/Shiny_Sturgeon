iter <- function(scores) data.frame(class = names(scores), score = unname(scores), stringsAsFactors = FALSE)
t0 <- as.POSIXct("2026-01-01 10:00:00")

test_that("update_progress starts a table and adds a TIME row", {
  p <- update_progress(data.frame(), iter(c(A = 0.6, B = 0.4)), "iteration_1", t0)
  expect_equal(colnames(p), c("class", "iteration_1"))
  expect_equal(p$iteration_1[p$class == "TIME"], "10:00")
})

test_that("update_progress adds one column per iteration and keeps classes aligned", {
  p <- update_progress(data.frame(), iter(c(A = 0.6, B = 0.4)), "iteration_1", t0)
  p <- update_progress(p, iter(c(B = 0.7, A = 0.3)), "iteration_2", t0 + 300)
  expect_equal(colnames(p), c("class", "iteration_1", "iteration_2"))
  expect_equal(p$iteration_2[p$class == "B"], "0.7")
  expect_equal(p$iteration_1[p$class == "B"], "0.4")
  expect_equal(p$iteration_2[p$class == "TIME"], "10:05")
})

test_that("the progress table survives a write.csv/read.csv round trip (resume)", {
  p <- update_progress(data.frame(), iter(c("Glioma IDH - A IDH - HG" = 0.6, B = 0.4)), "iteration_1", t0)
  p <- update_progress(p, iter(c("Glioma IDH - A IDH - HG" = 0.8, B = 0.2)), "iteration_2", t0)
  f <- tempfile(fileext = ".txt")
  write.csv(p, f, row.names = FALSE)
  back <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  expect_true("class" %in% colnames(back))
  expect_true("Glioma IDH - A IDH - HG" %in% back$class)
  expect_equal(latest_scores_table(back)$class[1], "Glioma IDH - A IDH - HG")
})

test_that("an old space-separated progress file has no class column after read.csv", {
  back <- read.csv(fixture("classifier_progress_old_format.txt"), stringsAsFactors = FALSE, check.names = FALSE)
  expect_false("class" %in% colnames(back))
  expect_equal(nrow(latest_scores_table(back)), 0)
})

test_that("latest_scores_table sorts the last iteration and drops TIME", {
  p <- update_progress(data.frame(), iter(c(A = 0.2, B = 0.5, C = 0.3)), "iteration_1", t0)
  tbl <- latest_scores_table(p)
  expect_equal(tbl$class, c("B", "C", "A"))
  expect_false("TIME" %in% tbl$class)
  expect_equal(nrow(latest_scores_table(data.frame())), 0)
})

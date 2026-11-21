ct <- list("Glioblastoma - GBM - RTK II" = "#123456", "Glioma IDH - A IDH - HG" = "#abcdef")

test_that("class_color uses color_translation, then the v2 mapping, then the fallback", {
  mapping <- load_class_mapping(file.path(repo_root, "class_mapping_v2.csv"))
  expect_equal(class_color("Glioblastoma - GBM - RTK II", ct), "#123456")
  expect_equal(class_color("GBM_RTK2", ct, mapping), "#123456")
  expect_equal(class_color("GBM_RTK2", ct), fallback_class_color("GBM_RTK2"))
  expect_equal(class_color("Unknown", ct, mapping), fallback_class_color("Unknown"))
})

test_that("class_color works with a plain named vector and NULL color_translation", {
  vec <- c("Glioma IDH - A IDH - HG" = "#abcdef")
  expect_equal(class_color("Glioma IDH - A IDH - HG", vec), "#abcdef")
  expect_equal(class_color("missing", vec), fallback_class_color("missing"))
  expect_equal(class_color("missing", NULL), fallback_class_color("missing"))
})

test_that("fallback colors are stable and valid", {
  expect_equal(fallback_class_color("A_IDH_HG"), fallback_class_color("A_IDH_HG"))
  expect_silent(grDevices::col2rgb(fallback_class_color("A_IDH_HG")))
})

test_that("the bundled class mapping file loads", {
  mapping <- load_class_mapping(file.path(repo_root, "class_mapping_v2.csv"))
  expect_true(length(mapping) >= 1)
  expect_equal(mapping[["A_IDH_HG"]], "Glioma IDH - A IDH - HG")
  expect_null(load_class_mapping(""))
})

test_that("load_class_mapping rejects bad files", {
  f <- tempfile(fileext = ".csv")
  writeLines(c("v2,v1", "A,B"), f)
  expect_error(load_class_mapping(f), "missing column")
  writeLines(c("v2_class,v1_class", "A,B", "A,C"), f)
  expect_error(load_class_mapping(f), "more than once")
})

test_that("plot_confidence_over_time draws empty, short and long progress tables", {
  pdf(NULL); on.exit(dev.off())
  expect_silent(plot_confidence_over_time(data.frame(), ct))
  expect_silent(plot_confidence_over_time(data.frame(x = 1), ct))  # no class column

  t0 <- as.POSIXct("2026-01-01 10:00:00")
  p <- data.frame()
  for (i in 1:30) {
    p <- update_progress(p, data.frame(class = c("GBM_RTK2", "A_IDH_HG"), score = c(i / 30, 1 - i / 30)),
                         paste0("iteration_", i), t0 + 60 * i)
  }
  mapping <- load_class_mapping(file.path(repo_root, "class_mapping_v2.csv"))
  expect_silent(plot_confidence_over_time(p, ct, class_mapping = mapping))
})

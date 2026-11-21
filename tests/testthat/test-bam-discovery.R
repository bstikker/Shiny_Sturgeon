make_bam_pass <- function() {
  root <- tempfile("bam_pass_")
  dir.create(file.path(root, "barcode01", "subdir"), recursive = TRUE)
  dir.create(file.path(root, "barcode02"))
  dir.create(file.path(root, "unclassified"))
  file.create(file.path(root, "barcode01", c("run_0.bam", "run_1.bam", "run_1.bam.bai")))
  file.create(file.path(root, "barcode02", "other.bam"))
  file.create(file.path(root, "unclassified", "uncl_0.bam"))
  root
}

test_that("find_new_bams lists the barcode's BAMs, skipping indexes and subfolders", {
  root <- make_bam_pass()
  found <- find_new_bams(root, "01", use_unclassified = FALSE, already_processed = character(0))
  expect_setequal(basename(found), c("run_0.bam", "run_1.bam"))
})

test_that("find_new_bams adds unclassified reads only when asked", {
  root <- make_bam_pass()
  found <- find_new_bams(root, "01", use_unclassified = TRUE, already_processed = character(0))
  expect_setequal(basename(found), c("run_0.bam", "run_1.bam", "uncl_0.bam"))
})

test_that("find_new_bams leaves out files that were already processed", {
  root <- make_bam_pass()
  done <- file.path(root, "barcode01", "run_0.bam")
  found <- find_new_bams(root, "01", use_unclassified = FALSE, already_processed = done)
  expect_equal(basename(found), "run_1.bam")
})

test_that("find_new_bams returns nothing for a barcode folder that doesn't exist yet", {
  root <- make_bam_pass()
  expect_length(find_new_bams(root, "05", FALSE, character(0)), 0)
})

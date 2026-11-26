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

test_that("the 'no new BAMs' log message names the watched folders", {
  bam_pass <- tempfile("bam_pass_"); dir.create(file.path(bam_pass, "barcode03"), recursive = TRUE)
  msg <- no_new_bams_message(bam_pass, "03", FALSE)
  expect_length(msg, 1)
  expect_match(msg, file.path(bam_pass, "barcode03"), fixed = TRUE)

  msg <- no_new_bams_message(bam_pass, "03", TRUE)
  expect_match(msg[1], "unclassified", fixed = TRUE)
  expect_match(msg[2], "WARNING: .*unclassified does not exist")
})

test_that("the 'no new BAMs' log message warns when the barcode or MinKNOW folder is missing", {
  bam_pass <- tempfile("bam_pass_"); dir.create(bam_pass)
  msg <- no_new_bams_message(bam_pass, "05", FALSE)
  expect_match(msg[2], "barcode05 does not exist")

  msg <- no_new_bams_message(file.path(bam_pass, "typo"), "05", TRUE)
  expect_length(msg, 2)
  expect_match(msg[2], "MinKNOW folder .*typo does not exist")
})

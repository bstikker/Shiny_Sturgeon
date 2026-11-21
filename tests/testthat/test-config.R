valid_cfg <- function() {
  dir <- tempfile("cfg_"); dir.create(dir)
  bam_pass <- file.path(dir, "bam_pass"); dir.create(bam_pass)
  files <- file.path(dir, c("probes.bed", "model.zip", "bins.rds", "ref.bed", "colors.rds"))
  file.create(files)
  list(bam_pass_folder = bam_pass, output_folder = file.path(dir, "case"),
       probes_file = files[1], model_file_v1 = files[2], bins_rds = files[3],
       cnvplot_reffile = files[4], color_translation_rds = files[5],
       model_file_v2 = "", sturgeon_v2_bin = "sturgeon-v2", sturgeon_v2_mode = "single",
       sturgeon_v2_file_type = "bed")
}

test_that("the bundled default config loads", {
  cfg <- load_config(repo_root)
  expect_equal(cfg$barcode, "01")
  expect_equal(cfg$sturgeon_v2_file_type, "bed")
  expect_true("class_mapping_v2" %in% names(cfg))
})

test_that("a complete config validates", {
  expect_length(validate_config(valid_cfg()), 0)
})

test_that("missing required fields are reported", {
  cfg <- valid_cfg(); cfg$model_file_v1 <- ""
  expect_match(validate_config(cfg), "'model_file_v1' is empty", all = FALSE)
})

test_that("files that don't exist are reported, including the optional ones", {
  cfg <- valid_cfg(); cfg$class_mapping_v2 <- "/nonexistent/mapping.csv"
  expect_match(validate_config(cfg), "class_mapping_v2", all = FALSE)
})

test_that("relative output and bam_pass folders are rejected", {
  cfg <- valid_cfg(); cfg$output_folder <- "relative/case"
  expect_match(validate_config(cfg), "absolute path", all = FALSE)
})

test_that("v2 parameters become required once a v2 model is set", {
  cfg <- valid_cfg()
  cfg$model_file_v2 <- cfg$model_file_v1
  cfg$sturgeon_v2_file_type <- ""
  cfg$sturgeon_v2_bin <- NA
  problems <- validate_config(cfg)
  expect_match(problems, "sturgeon_v2_file_type", all = FALSE)
  expect_match(problems, "sturgeon_v2_bin", all = FALSE)
})

test_that("normalize_cfg_paths expands ~ in path fields only", {
  cfg <- normalize_cfg_paths(list(output_folder = "~/case", samtools_bin = "samtools",
                                  barcode = "~01", extra_path_dirs = c("~/bin", "")))
  expect_false(startsWith(cfg$output_folder, "~"))
  expect_equal(cfg$samtools_bin, "samtools")
  expect_equal(cfg$barcode, "~01")
  expect_false(startsWith(unname(cfg$extra_path_dirs[1]), "~"))
})

test_that("save_config writes a config.yml that loads back", {
  cfg <- valid_cfg()
  path <- save_config(cfg, cfg$output_folder)
  expect_true(file.exists(path))
  expect_equal(yaml::read_yaml(path)$probes_file, cfg$probes_file)
})

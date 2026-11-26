valid_cfg <- function() {
  dir <- tempfile("cfg_"); dir.create(dir)
  bam_pass <- file.path(dir, "bam_pass"); dir.create(bam_pass)
  files <- file.path(dir, c("probes.bed", "model.zip", "bins.rds", "ref.bed", "colors.rds", "model_v2.zip"))
  file.create(files)
  list(bam_pass_folder = bam_pass, output_folder = file.path(dir, "case"),
       probes_file = files[1], model_file_v1 = files[2], bins_rds = files[3],
       cnvplot_reffile = files[4], color_translation_rds = files[5],
       model_file_v2 = files[6], sturgeon_v2_bin = "sturgeon-v2", sturgeon_v2_mode = "single",
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

test_that("the v2 model is required, like the v1 model", {
  cfg <- valid_cfg(); cfg$model_file_v2 <- ""
  expect_match(validate_config(cfg), "'model_file_v2' is empty", all = FALSE)
})

test_that("the v2 parameters are required, except num_cpgs and num_samplings", {
  cfg <- valid_cfg()
  cfg$sturgeon_v2_num_cpgs <- ""
  cfg$sturgeon_v2_num_samplings <- NA
  expect_length(validate_config(cfg), 0)
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

test_that("config_local.yml (written by install.sh) overrides the default, and a saved case config overrides both", {
  dir <- tempfile("app_"); dir.create(dir)
  file.copy(file.path(repo_root, "config_default.yml"), dir)
  yaml::write_yaml(list(samtools_bin = "/opt/env/bin/samtools", probes_file = "/data/probes.bed"),
                   file.path(dir, "config_local.yml"))
  cfg <- load_config(dir)
  expect_equal(cfg$samtools_bin, "/opt/env/bin/samtools")
  expect_equal(cfg$probes_file, "/data/probes.bed")
  expect_equal(cfg$barcode, "01")

  saved <- file.path(dir, "saved.yml")
  yaml::write_yaml(list(probes_file = "/case/probes.bed"), saved)
  cfg <- load_config(dir, saved_path = saved)
  expect_equal(cfg$probes_file, "/case/probes.bed")
  expect_equal(cfg$samtools_bin, "/opt/env/bin/samtools")
})

test_that("load_config works without a config_local.yml", {
  expect_equal(load_config(repo_root, local_path = tempfile())$samtools_bin, "samtools")
})

test_that("resolve_tool finds tools by full path, in extra PATH directories, and on PATH", {
  bin_dir <- tempfile("bin_"); dir.create(bin_dir)
  tool <- fake_tool(bin_dir, "my-tool")
  expect_equal(resolve_tool(tool), tool)
  expect_equal(resolve_tool("my-tool", bin_dir), tool)
  expect_equal(resolve_tool("my-tool"), "")
  expect_equal(resolve_tool(file.path(bin_dir, "absent")), "")
  expect_equal(resolve_tool(""), "")
  expect_equal(resolve_tool(NA), "")
  expect_true(nzchar(resolve_tool("sh")))
  expect_true(nzchar(resolve_tool("sh", list())))  # extra_path_dirs: [] read from YAML
  expect_true(nzchar(resolve_tool("sh", NULL)))

  not_exec <- file.path(bin_dir, "not-exec"); file.create(not_exec)
  expect_equal(resolve_tool(not_exec), "")
})

setup_cfg <- function() {
  cfg <- valid_cfg()
  bin_dir <- tempfile("bin_"); dir.create(bin_dir)
  for (t in c("samtools", "modkit", "sturgeon", "sturgeon-v2")) fake_tool(bin_dir, t)
  dir.create(file.path(cfg$bam_pass_folder, "barcode01"))
  c(cfg, list(samtools_bin = "samtools", modkit_bin = "modkit", sturgeon_bin = "sturgeon",
              extra_path_dirs = bin_dir, barcode = "01", use_unclassified = FALSE,
              relevant_genes_bed = "", class_mapping_v2 = "", probes_file_v2 = ""))
}
status_of <- function(res, item) res$status[res$item == item]

test_that("check_setup reports a complete setup as ok", {
  res <- check_setup(setup_cfg())
  expect_false(any(res$status %in% c("problem", "warning")))
  expect_equal(status_of(res, "samtools"), "ok")
  expect_equal(status_of(res, "sturgeon-v2"), "ok")
  expect_equal(status_of(res, "Sturgeon v2 model"), "ok")
  expect_equal(status_of(res, "CNV gene labels"), "not used")
  expect_match(res$detail[res$item == "Output folder"], "will be created")
})

test_that("check_setup reports missing tools, files and folders", {
  cfg <- setup_cfg()
  cfg$modkit_bin <- "/nonexistent/modkit"
  cfg$bins_rds <- "/nonexistent/bins.rds"
  cfg$color_translation_rds <- ""
  cfg$bam_pass_folder <- "/nonexistent/bam_pass"
  res <- check_setup(cfg)
  expect_equal(status_of(res, "modkit"), "problem")
  expect_equal(status_of(res, "CNV bins"), "problem")
  expect_equal(status_of(res, "Class colors"), "problem")
  expect_equal(status_of(res, "MinKNOW folder"), "problem")
  expect_equal(status_of(res, "samtools"), "ok")
})

test_that("check_setup requires sturgeon-v2 and the v2 model", {
  cfg <- setup_cfg()
  cfg$sturgeon_v2_bin <- "sturgeon-v2-absent"
  cfg$model_file_v2 <- ""
  res <- check_setup(cfg)
  expect_equal(status_of(res, "sturgeon-v2"), "problem")
  expect_equal(status_of(res, "Sturgeon v2 model"), "problem")
  fake_tool(cfg$extra_path_dirs, "sturgeon-v2-absent")
  expect_equal(status_of(check_setup(cfg), "sturgeon-v2"), "ok")
})

test_that("check_setup warns about a barcode folder that doesn't exist yet and an output folder to resume", {
  cfg <- setup_cfg()
  cfg$barcode <- "07"
  cfg$use_unclassified <- TRUE
  dir.create(cfg$output_folder)
  file.create(file.path(cfg$output_folder, "processed_files.txt"))
  res <- check_setup(cfg)
  watched <- res[res$item == "Watched folder", ]
  expect_equal(nrow(watched), 2)
  expect_true(all(watched$status == "warning"))
  expect_match(watched$detail[1], "barcode07")
  expect_equal(status_of(res, "Output folder"), "warning")
  expect_false(any(res$status == "problem"))
})

test_that("check_setup rejects a relative output folder", {
  cfg <- setup_cfg(); cfg$output_folder <- "case"
  expect_equal(status_of(check_setup(cfg), "Output folder"), "problem")
})

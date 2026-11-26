# =============================================================================
# make_demo_case.R -- build a demo case folder that opens in the app with
# populated confidence-over-time plots, without running the pipeline. It is
# what the README screenshot shows. To test an installation for real, use
# tools/make_test_case.R instead.
#
#   Rscript tools/make_demo_case.R /tmp/shiny_sturgeon_demo
#
# Then start the app, enter <folder>/case as the output folder on the Setup
# tab, click "Load config.yml from output folder", and click "Start
# monitoring" on the Live Run tab. The app resumes the "previous run" it
# finds in the folder and draws its plots.
#
# The scores are made-up curves and bam_pass/ is empty, so starting the demo
# never runs samtools, modkit or sturgeon. Tool and reference paths come
# from this machine's installation (config_local.yml, written by
# install.sh); reference files that aren't installed are replaced by empty
# placeholders. "Start monitoring" checks that the tools exist, so it only
# works after install.sh. There is no CNV plot (that needs BAM files).
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript tools/make_demo_case.R <empty folder to create>")
root <- normalizePath(args[1], mustWork = FALSE)
if (dir.exists(root) && length(list.files(root, all.files = TRUE, no.. = TRUE)) > 0) {
  stop(root, " already exists and is not empty -- pick a new folder")
}

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
app_dir <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), ".."))

case_dir <- file.path(root, "case")
ref_dir <- file.path(root, "reference")
bam_pass <- file.path(root, "bam_pass")
for (d in c(case_dir, ref_dir, file.path(bam_pass, "barcode01"), file.path(bam_pass, "unclassified"))) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ---- installed paths, or placeholders where nothing is installed ------------
source(file.path(app_dir, "config.R"))
cfg <- load_config(app_dir)  # includes config_local.yml from install.sh, if any
installed <- function(field) {
  v <- cfg[[field]]
  !is.null(v) && length(v) == 1 && nzchar(v) && file.exists(v)
}
placeholder <- function(field, name, writer = writeLines) {
  if (installed(field)) return(cfg[[field]])
  path <- file.path(ref_dir, name)
  writer("synthetic demo placeholder -- not a real reference file", path)
  path
}
cfg$probes_file <- placeholder("probes_file", "probes_PLACEHOLDER.bed")
cfg$model_file_v1 <- placeholder("model_file_v1", "v1_model_PLACEHOLDER.zip")
cfg$model_file_v2 <- placeholder("model_file_v2", "v2_model_PLACEHOLDER.zip")
cfg$cnvplot_reffile <- placeholder("cnvplot_reffile", "cnv_reference_PLACEHOLDER.bed")
cfg$bins_rds <- placeholder("bins_rds", "bins_PLACEHOLDER.rds", function(x, f) saveRDS(x, f))

v1_classes <- c("Glioblastoma - GBM - RTK II", "Glioblastoma - GBM - RTK I",
                "Glioma IDH - A IDH - HG", "Glioma IDH - O IDH - O IDH",
                "Glioblastoma - DMG - K27")
color_translation <- setNames(as.list(c("#c0392b", "#e67e22", "#2471a3", "#17a589", "#7d3c98")),
                              v1_classes)
# the real color_translation.rds (from install.sh) uses these same class names
if (installed("color_translation_rds")) {
  color_rds <- cfg$color_translation_rds
} else {
  color_rds <- file.path(ref_dir, "color_translation_demo.rds")
  saveRDS(color_translation, color_rds)
}

# ---- synthetic confidence-over-time curves ----------------------------------
# The top class climbs as reads accumulate, the others share what is left.
n_iter <- 14
start_time <- as.POSIXct("2026-01-01 10:02:00", tz = "UTC")

make_progress <- function(classes, top_curve, seed) {
  set.seed(seed)
  others <- matrix(runif((length(classes) - 1) * n_iter), ncol = n_iter)
  others[1, ] <- others[1, ] + 1.5  # a clear runner-up
  progress <- data.frame()
  for (i in seq_len(n_iter)) {
    rest <- others[, i] / sum(others[, i]) * (1 - top_curve[i])
    scores <- signif(c(top_curve[i], rest), 4)
    it <- data.frame(class = classes, score = as.character(scores), stringsAsFactors = FALSE)
    it <- rbind(it, data.frame(class = "TIME", score = format(start_time + (i - 1) * 300, "%H:%M")))
    colnames(it)[2] <- paste0("iteration_", i)
    progress <- if (nrow(progress) == 0) it else merge(progress, it, by = "class", all = TRUE)
  }
  progress
}

x <- seq_len(n_iter)
top_v1 <- pmin(0.98, 0.25 + 0.75 / (1 + exp(-(x - 6) / 1.6)))
top_v2 <- pmin(0.99, 0.30 + 0.70 / (1 + exp(-(x - 5) / 1.4)))

write.csv(make_progress(v1_classes, top_v1, seed = 1),
          file.path(case_dir, "classifier_progress_v1.txt"), row.names = FALSE)
write.csv(make_progress(c("GBM_RTK2", "GBM_RTK1", "A_IDH_HG", "O_IDH", "DMG_K27"), top_v2, seed = 2),
          file.path(case_dir, "classifier_progress_v2.txt"), row.names = FALSE)

# ---- what a resumed run expects to find --------------------------------------
for (i in seq_len(n_iter)) dir.create(file.path(case_dir, paste0("iteration", i)), showWarnings = FALSE)
writeLines(file.path(bam_pass, "barcode01", sprintf("synthetic_demo_%02d.bam", seq_len(n_iter))),
           file.path(case_dir, "processed_files.txt"))

cfg$bam_pass_folder <- bam_pass
cfg$output_folder <- case_dir
cfg$poll_interval_sec <- 60L
cfg$color_translation_rds <- color_rds
if (!nzchar(if (is.null(cfg$class_mapping_v2)) "" else cfg$class_mapping_v2)) cfg$class_mapping_v2 <- file.path(app_dir, "class_mapping_v2.csv")
invisible(save_config(cfg, case_dir))

tools_found <- all(vapply(c("samtools_bin", "modkit_bin", "sturgeon_bin", "sturgeon_v2_bin"),
                          function(f) nzchar(resolve_tool(cfg[[f]], unlist(cfg$extra_path_dirs))), logical(1)))
cat("Demo case written to", root, "\n",
    "In the app: Setup tab -> output folder =", case_dir, "\n",
    "-> 'Load config.yml from output folder' -> Live Run tab -> 'Start monitoring'\n")
if (!tools_found) {
  cat("Note: samtools, modkit, sturgeon or sturgeon-v2 isn't installed here, so 'Start monitoring'\n",
      "will refuse to start. Run ./install.sh first.\n")
}

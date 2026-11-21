# =============================================================================
# make_demo_case.R -- build a demo case folder that opens in the app with
# populated confidence-over-time plots, without sequencing data, model files
# or external tools.
#
#   Rscript tools/make_demo_case.R /tmp/shiny_sturgeon_demo
#
# Then start the app, enter <folder>/case as the output folder on the Setup
# tab, click "Load config.yml from output folder", and click "Start
# monitoring" on the Live Run tab. The app resumes the "previous run" it
# finds in the folder and draws its plots.
#
# Everything here is SYNTHETIC: the scores are made-up curves, the
# reference/model files are empty placeholders, and bam_pass/ is empty.
# Starting the demo never runs samtools, modkit or sturgeon, because there
# are no new BAM files to process. There is no CNV plot (that needs real
# QDNAseq bins and BAM files).
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

# ---- placeholder reference files (they only need to exist) ------------------
placeholders <- file.path(ref_dir, c("probes_PLACEHOLDER.bed", "v1_model_PLACEHOLDER.zip",
                                     "v2_model_PLACEHOLDER.zip", "cnv_reference_PLACEHOLDER.bed"))
for (f in placeholders) writeLines("synthetic demo placeholder -- not a real reference file", f)
bins_rds <- file.path(ref_dir, "bins_PLACEHOLDER.rds")
saveRDS("synthetic demo placeholder -- not real QDNAseq bins", bins_rds)

v1_classes <- c("Glioblastoma - GBM - RTK II", "Glioblastoma - GBM - RTK I",
                "Glioma IDH - A IDH - HG", "Glioma IDH - O IDH - O IDH",
                "Glioblastoma - DMG - K27")
color_translation <- setNames(as.list(c("#c0392b", "#e67e22", "#2471a3", "#17a589", "#7d3c98")),
                              v1_classes)
color_rds <- file.path(ref_dir, "color_translation_demo.rds")
saveRDS(color_translation, color_rds)

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

cfg <- yaml::read_yaml(file.path(app_dir, "config_default.yml"))
cfg$bam_pass_folder <- bam_pass
cfg$output_folder <- case_dir
cfg$poll_interval_sec <- 60
cfg$samtools_bin <- "demo-not-used"
cfg$modkit_bin <- "demo-not-used"
cfg$sturgeon_bin <- "demo-not-used"
cfg$sturgeon_v2_bin <- "demo-not-used"
cfg$probes_file <- placeholders[1]
cfg$model_file_v1 <- placeholders[2]
cfg$model_file_v2 <- placeholders[3]
cfg$cnvplot_reffile <- placeholders[4]
cfg$bins_rds <- bins_rds
cfg$color_translation_rds <- color_rds
cfg$class_mapping_v2 <- file.path(app_dir, "class_mapping_v2.csv")
yaml::write_yaml(cfg, file.path(case_dir, "config.yml"))

cat("Demo case written to", root, "\n",
    "In the app: Setup tab -> output folder =", case_dir, "\n",
    "-> 'Load config.yml from output folder' -> Live Run tab -> 'Start monitoring'\n")

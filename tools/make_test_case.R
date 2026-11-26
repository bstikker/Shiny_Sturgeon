# =============================================================================
# make_test_case.R -- test an installation end to end with synthetic reads.
#
#   Rscript tools/make_test_case.R ~/sturgeon_test
#
# Builds <folder>/bam_pass/ (a stand-in for a MinKNOW run folder) and
# <folder>/case/config.yml, which uses this machine's tools and reference
# files (from config_local.yml, written by install.sh) and points at that
# bam_pass/. It then writes synthetic BAM files into bam_pass/barcode01/, one
# every 30 seconds, like MinKNOW does during a run.
#
# While it runs, start the app (./run_app.sh) and on the Setup tab enter
# <folder>/case as the output folder, click "Load config.yml from output
# folder", "Check installation", and then "Start monitoring" on the Live Run
# tab. Every new file goes through samtools, modkit, Sturgeon v1 and v2, and
# every 3rd iteration the CNV plot is rebuilt (the reads include a fake chr7
# gain and chr10 loss).
#
# Options:
#   --files N        BAM files to write for barcode01 (default 6)
#   --reads N        reads per file (default 1000)
#   --interval SEC   seconds between files (default 30; 0 writes all at once)
#
# The methylation calls are RANDOM: this tests that the pipeline runs, and
# any classification made from it is meaningless.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
usage <- "Usage: Rscript tools/make_test_case.R <new folder> [--files N] [--reads N] [--interval SEC]"
opt <- list(files = 6, reads = 1000, interval = 30)
folder <- NULL
i <- 1
while (i <= length(args)) {
  a <- args[i]
  if (a %in% c("--files", "--reads", "--interval")) {
    if (i == length(args)) stop(usage)
    value <- suppressWarnings(as.numeric(args[i + 1]))
    if (is.na(value) || value < 0) stop(a, " needs a number. ", usage)
    opt[[sub("^--", "", a)]] <- value
    i <- i + 2
  } else if (is.null(folder) && !startsWith(a, "--")) {
    folder <- a
    i <- i + 1
  } else {
    stop(usage)
  }
}
if (is.null(folder)) stop(usage)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
app_dir <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), ".."))
source(file.path(app_dir, "config.R"))

root <- normalizePath(path.expand(folder), mustWork = FALSE)
if (dir.exists(root) && length(list.files(root, all.files = TRUE, no.. = TRUE)) > 0) {
  stop(root, " already exists and is not empty -- pick a new folder")
}

# ---- this machine's tools and reference files --------------------------------
cfg <- load_config(app_dir)
if (!file.exists(file.path(app_dir, "config_local.yml"))) {
  message("Note: no config_local.yml (written by install.sh); using the paths in config_default.yml.")
}
extra <- unlist(cfg$extra_path_dirs)
samtools <- resolve_tool(cfg$samtools_bin, extra)
if (!nzchar(samtools)) stop("samtools not found (", cfg$samtools_bin, "). Run ./install.sh first.")
if (is.null(cfg$probes_file) || !nzchar(cfg$probes_file) || !file.exists(cfg$probes_file)) {
  stop("Probes file not found (", cfg$probes_file, "). Run ./install.sh first.")
}
python <- resolve_tool("python3", extra)
if (!nzchar(python)) stop("python3 not found")

# ---- the test case ---------------------------------------------------------------
bam_pass <- file.path(root, "bam_pass")
case_dir <- file.path(root, "case")
dir.create(bam_pass, recursive = TRUE)

cfg$bam_pass_folder <- bam_pass
cfg$output_folder <- case_dir
cfg$barcode <- "01"
cfg$use_unclassified <- TRUE
cfg$poll_interval_sec <- 15L
cfg$cnv_every_n_iterations <- 3L
config_path <- save_config(cfg, case_dir)

cat("\nTest case created in ", root, "\n\n",
    "Now, while this keeps running:\n",
    "  1. start the app:  ", file.path(app_dir, "run_app.sh"), "  (in another terminal)\n",
    "  2. Setup tab: output folder = ", case_dir, "\n",
    "     click 'Load config.yml from output folder', then 'Check installation'\n",
    "  3. Live Run tab: click 'Start monitoring'\n\n",
    sprintf("Writing %d synthetic BAM files into %s%s. Press Ctrl+C to stop early.\n\n",
            opt$files, file.path(bam_pass, "barcode01"),
            if (opt$interval > 0) sprintf(", one every %g seconds", opt$interval) else ""),
    sep = "")

status <- system2(python, c("-I", shQuote(file.path(app_dir, "tools", "make_synthetic_bams.py")),
                            "--out", shQuote(bam_pass),
                            "--n-files", opt$files, "--reads-per-file", opt$reads,
                            "--unclassified-files", 1,
                            "--probes", shQuote(cfg$probes_file),
                            "--gain", 7, "--loss", 10,
                            "--samtools", shQuote(samtools),
                            "--live-interval", opt$interval))
if (!identical(status, 0L)) stop("make_synthetic_bams.py failed")
cat("\nAll files written. The app processes the last ones on its next check.\n")

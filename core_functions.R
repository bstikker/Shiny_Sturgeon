# =============================================================================
# core_functions.R
#
# Cleaned, Linux-native pipeline logic for live intraoperative Sturgeon
# classification.
#
# This replaces:
#   - live_processing_windows_minknow_basecalling.R  (Windows/WSL polling loop)
#   - peroperative_seq_functions_shorter.R (function library)
#   - sturgeon_guppy_modkit_SLURM.sh                  (Sturgeon CLI wrapper)
#
# Dropped on purpose (dead weight for the current MinKNOW-does-basecalling
# workflow):
#   - wrappert_guppy_R10_guppy6.5 / _V2 / _V2_modkit   (old fast5->guppy path,
#     obsolete now that MinKNOW basecalls + aligns in real time)
#   - functioncheck()                                  (old megalodon/guppy
#     self-test with hardcoded cluster paths)
#   - all Windows-drive-letter <-> /mnt/c/... path translation and the
#     "wsl <cmd>" prefixing done via shell() (Windows-only R function)
#   - the duplicate "_windows" copies of add_and_plot()/confidence_over_time_plot()
#
# Every external path (samtools, sturgeon, reference files) is now passed in
# explicitly via the config list produced by config.R, instead of being
# picked up from loose global variables.
# =============================================================================

library(DNAcopy)
library(QDNAseq)

# ---- small helpers ----------------------------------------------------------

#' Zero-pad a barcode number/string to two digits, e.g. 1 -> "01", "7" -> "07"
pad_barcode <- function(barcode) {
  barcode <- as.character(barcode)
  if (nchar(barcode) == 1) barcode <- paste0("0", barcode)
  barcode
}

#' Prepend extra directories (e.g. a conda env's bin/) to PATH for this R
#' session, so samtools/sturgeon/modkit resolve correctly regardless of how
#' the app itself was launched.
apply_extra_path <- function(extra_path_dirs) {
  extra_path_dirs <- extra_path_dirs[nzchar(extra_path_dirs)]
  if (length(extra_path_dirs) == 0) return(invisible(NULL))
  current <- Sys.getenv("PATH")
  Sys.setenv(PATH = paste(paste(extra_path_dirs, collapse = ":"), current, sep = ":"))
  invisible(NULL)
}

#' Run a shell command, capture output, and raise a real R error (with the
#' command's stderr/stdout attached) if it fails instead of silently
#' continuing -- the original scripts used shell()/system() and ignored
#' exit codes, so a failed samtools/sturgeon call could pass through
#' unnoticed.
#'
#' stderr is explicitly redirected into the captured stream (`2>&1`).
#' `system(..., intern = TRUE)` on its own only captures stdout -- stderr
#' still goes to the R process's own stderr, invisibly, wherever that
#' happens to be connected (a terminal Shiny isn't attached to, a log file
#' nobody's watching, or nowhere). Since Python CLIs like sturgeon-v2 print
#' their real tracebacks to stderr, that gap meant a failure only ever
#' surfaced here as its generic top-level "an exception occurred" message,
#' with the actual cause silently discarded instead of showing up in the
#' app's log.
run_cmd <- function(cmd, log_fun = message) {
  log_fun(paste0("$ ", cmd))
  out <- system(paste0(cmd, " 2>&1"), intern = TRUE)
  status <- attr(out, "status")
  if (!is.null(status) && status != 0) {
    stop(sprintf("Command failed (exit %s): %s\n%s", status, cmd, paste(out, collapse = "\n")))
  }
  invisible(out)
}

# ---- BAM discovery & merging ------------------------------------------------

#' Find (and optionally unclassified) BAM files for a barcode that have not
#' been processed yet.
#'
#' Matches the current live-analysis script's folder layout: MinKNOW's
#' realtime BAM output folder (named e.g. "bam_pass/", NOT "pod5/pass/" --
#' this app used to assume the older WSL-era layout with an extra "pass/"
#' level, which was wrong for the real Linux setup) contains "barcode<XX>/"
#' directly, and "unclassified/" as a sibling of the barcode folders.
#'
#' @param bam_pass_folder MinKNOW's realtime BAM output folder (the one
#'   containing barcode01/, barcode02/, ..., unclassified/)
#' @param barcode     zero-padded barcode string, e.g. "01"
#' @param use_unclassified include reads MinKNOW couldn't assign a barcode to
#' @param already_processed character vector of file paths already handled
find_new_bams <- function(bam_pass_folder, barcode, use_unclassified, already_processed) {
  # no filename pattern filter (e.g. "bam_runid") -- the real script just
  # takes every file in the barcode folder except .bai indexes, since
  # MinKNOW's realtime bam output filenames aren't guaranteed to match a
  # fixed substring
  list_bams <- function(dir) {
    files <- list.files(dir, full.names = TRUE)
    files <- files[!grepl("\\.bai$", files)]
    files[!file.info(files)$isdir %in% TRUE]  # skip any subdirectories, if present
  }

  barcode_dir <- file.path(bam_pass_folder, paste0("barcode", barcode))
  bams <- list_bams(barcode_dir)

  if (isTRUE(use_unclassified)) {
    unclassified_dir <- file.path(bam_pass_folder, "unclassified")
    bams <- c(bams, list_bams(unclassified_dir))
  }

  setdiff(bams, already_processed)
}

#' The folders find_new_bams() looks in for a barcode, and which of them do
#' not exist (yet). Used for the log: MinKNOW only creates barcode<NN>/ once
#' the first read for that barcode is written, so a missing folder early in a
#' run is normal, but one that never appears usually means a wrong path or
#' barcode.
describe_watch_folders <- function(bam_pass_folder, barcode, use_unclassified) {
  folders <- file.path(bam_pass_folder, paste0("barcode", barcode))
  if (isTRUE(use_unclassified)) folders <- c(folders, file.path(bam_pass_folder, "unclassified"))
  list(folders = folders, missing = folders[!dir.exists(folders)])
}

#' Log line(s) for a poll that found nothing new: names the folders being
#' watched, and warns about each one that does not exist.
no_new_bams_message <- function(bam_pass_folder, barcode, use_unclassified) {
  watch <- describe_watch_folders(bam_pass_folder, barcode, use_unclassified)
  msg <- sprintf("No new BAM files in %s%s, waiting...", watch$folders[1],
                 if (isTRUE(use_unclassified)) " (or unclassified/)" else "")
  if (!dir.exists(bam_pass_folder)) {
    return(c(msg, sprintf("WARNING: the MinKNOW folder %s does not exist. Check 'bam_pass_folder' on the Setup tab.",
                          bam_pass_folder)))
  }
  for (d in watch$missing) {
    msg <- c(msg, sprintf("WARNING: %s does not exist (yet). MinKNOW creates it when the first reads for it arrive; if it never appears, check the barcode and MinKNOW folder on the Setup tab.", d))
  }
  msg
}

#' Merge + coordinate-sort a set of BAMs into one file for this iteration.
#' Returns the path to the sorted BAM.
merge_and_sort_bams <- function(bam_files, iteration_dir, samtools_bin = "samtools",
                                 threads = 4, log_fun = message) {
  if (!dir.exists(iteration_dir)) dir.create(iteration_dir, recursive = TRUE)

  merged_bam <- file.path(iteration_dir, "merged.bam")
  sorted_bam <- file.path(iteration_dir, "merged_srt.bam")

  run_cmd(sprintf('%s merge -f -@ %d -o "%s" %s',
                   samtools_bin, threads, merged_bam,
                   paste(shQuote(bam_files), collapse = " ")), log_fun)

  run_cmd(sprintf('%s sort -@ %d -o "%s" "%s"',
                   samtools_bin, threads, sorted_bam, merged_bam), log_fun)

  sorted_bam
}

# ---- modkit preprocessing ----------------------------------------------------
#
# Two explicit modkit
# steps happen first, per iteration, on that iteration's own new reads only:
#   1. adjust-mods (5hmC "h" calls folded into 5mC "m") -- always required
#   2. extract full -- writes a per-iteration modkit .txt file into a
#      CUMULATIVE directory that just keeps growing, one file per iteration,
#      nothing ever removed
# Only step 2's *output* accumulates -- step 1 and 2 themselves only ever
# process the current iteration's new reads, so this stays cheap as a case
# goes on. It's `sturgeon inputtobed` (below) that then aggregates every
# modkit .txt file written so far into one combined bed reflecting all reads
# processed to date -- which is what makes "confidence over time" meaningful
# (climbing as more cumulative data comes in), rather than each iteration
# being an independent, noisy classification of just its own small new batch.

#' modkit adjust-mods --convert h m: fold 5hmC calls into 5mC before
#' extraction. Returns the path to the adjusted BAM.
run_modkit_adjust_mods <- function(in_bam, out_bam, modkit_bin = "modkit", log_fun = message) {
  run_cmd(sprintf('%s adjust-mods --convert h m "%s" "%s"', modkit_bin, in_bam, out_bam), log_fun)
  if (!file.exists(out_bam)) stop("modkit adjust-mods did not produce: ", out_bam)
  out_bam
}

#' modkit extract full: pull per-read methylation calls out of the adjusted
#' BAM into a text file. Returns the path to that text file.
run_modkit_extract_full <- function(in_bam, out_txt, modkit_bin = "modkit", threads = 10, log_fun = message) {
  out_dir <- dirname(out_txt)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  run_cmd(sprintf('%s extract full -t %d "%s" "%s"', modkit_bin, threads, in_bam, out_txt), log_fun)
  if (!file.exists(out_txt)) stop("modkit extract full did not produce: ", out_txt)
  out_txt
}

# ---- Sturgeon classification ------------------------------------------------
#
# Split into two steps (inputtobed, predict) instead of one combined
# run_sturgeon() so a case can run two classifier models (e.g. Sturgeon v1
# and v2) against the same cumulative data. inputtobed (the methylation
# aggregation step) is the one both models should share when they use the
# same probes file, so running a second model is cheap, not a doubled
# pipeline.

#' Aggregate every modkit-extracted .txt file in `modkit_txt_dir` (cumulative
#' across all iterations so far, see above) into one combined bed file in
#' `bed_out_dir`. Re-running this after a new iteration's modkit.txt has been
#' added naturally re-aggregates everything, including the new data.
#'
#' `bed_out_dir` is wiped and recreated on every call. This is deliberate:
#' whether sturgeon names its bed output the same way every run, or derives
#' the name from its (growing, ever-changing) set of input files isn't
#' something to assume either way -- so any stale bed left over from a
#' previous call is removed first, guaranteeing exactly one (fresh, fully
#' up to date) bed file is ever present for find_single_bed_file() to find.
run_sturgeon_inputtobed <- function(modkit_txt_dir, bed_out_dir, probes_file,
                                     sturgeon_bin = "sturgeon", margin = 50, log_fun = message) {
  if (dir.exists(bed_out_dir)) unlink(bed_out_dir, recursive = TRUE)
  dir.create(bed_out_dir, recursive = TRUE)
  run_cmd(sprintf('%s inputtobed --margin %d -i "%s" -o "%s" -s modkit --probes-file "%s"',
                   sturgeon_bin, margin, modkit_txt_dir, bed_out_dir, probes_file), log_fun)
  invisible(bed_out_dir)
}

#' Find exactly one file matching `pattern` in `dir`. Used instead of
#' hardcoding output filenames, which turned out not to be safe: v1's result
#' file is actually named after the model (`merged_probes_methyl_calls_general.csv`
#' for the "general" model, confirmed against a real run), not a fixed
#' "..._model.csv" like the original scripts assumed. Errors loudly (rather
#' than guessing) if the match isn't unique -- a live case is the wrong place
#' to silently pick the wrong file.
find_single_file <- function(dir, pattern, what = "file") {
  matches <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(matches) == 0) {
    stop("No ", what, " found in ", dir, " matching ", pattern)
  }
  if (length(matches) > 1) {
    stop("Expected exactly one ", what, " in ", dir, " but found ", length(matches),
         ": ", paste(basename(matches), collapse = ", "))
  }
  matches[1]
}

#' Run Sturgeon prediction against an already-extracted bed folder for one
#' model, writing results into `predict_out_dir`. Two models can each get
#' their own predict_out_dir so their result CSVs don't overwrite each other.
#' The result filename is discovered rather than assumed (see find_single_file).
run_sturgeon_predict <- function(bed_dir, model_file, predict_out_dir,
                                  sturgeon_bin = "sturgeon", log_fun = message) {
  if (!dir.exists(predict_out_dir)) dir.create(predict_out_dir, recursive = TRUE)
  run_cmd(sprintf('%s predict -p --i "%s" -o "%s" -m "%s"',
                   sturgeon_bin, bed_dir, predict_out_dir, model_file), log_fun)

  find_single_file(predict_out_dir, "^merged_probes_methyl_calls_.*\\.csv$", "Sturgeon v1 result CSV")
}

# ---- Sturgeon v2 classification ----------------------------------------------
#
# sturgeon-v2 has a different CLI shape than v1: it's a separate program
# (`sturgeon-v2 single/multi ...`, not a subcommand of `sturgeon`), it takes
# a single bed FILE via -i rather than a directory, and it writes a single
# OUTPUT_FILE rather than writing into a directory. 

#' Find the single .bed file inputtobed wrote into `bed_dir`.
find_single_bed_file <- function(bed_dir) {
  find_single_file(bed_dir, "\\.bed$", "bed file")
}

#' Run sturgeon-v2 prediction against the bed file in `bed_dir`, and return
#' the path of the result file it wrote.
#'
#' `output_file` is passed to sturgeon-v2 as `-o`. Two sturgeon-v2 versions
#' are in use, and they treat `-o` differently:
#'   - earlier versions write a single comma-separated result FILE at that path;
#'   - the public release (github.com/deRidderLab/sturgeon-v2) treats `-o` as an
#'     output DIRECTORY and writes a tab-separated `result.txt` inside it.
#' resolve_sturgeon_v2_output() finds the result in either case, so the same
#' config works with both.
#'
#' NOTE: file_type is passed straight through from config. mode is "single"
#' or "multi" per sturgeon-v2's usage; "single" (one sample at a time) is
#' what this app is built around, since each iteration is one case's
#' accumulated reads.
#'
#' num_cpgs (-c) and num_samplings (-s) are OPTIONAL per sturgeon-v2's own
#' --help: left unset, sturgeon-v2 itself defaults to using all available
#' CpG sites and 1 random sampling. Confirmed (2026-09) that "-c all" is NOT
#' accepted -- -c takes a number only, and omitting the flag entirely is how
#' you get "all sites". So these are only added to the command when the
#' user has actually supplied a value; otherwise the flags are left out and
#' sturgeon-v2's own defaults apply.
run_sturgeon_v2_predict <- function(bed_dir, model_file, output_file, sturgeon_v2_bin = "sturgeon-v2",
                                     mode = "single", num_cpgs = NULL, num_samplings = NULL, file_type,
                                     log_fun = message) {
  input_bed <- find_single_bed_file(bed_dir)

  out_dir <- dirname(output_file)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  cmd <- sprintf('%s %s -i "%s" -m "%s" -o "%s" -f %s',
                  sturgeon_v2_bin, mode, input_bed, model_file, output_file, file_type)

  has_value <- function(x) !is.null(x) && !(length(x) == 1 && is.na(x)) && nzchar(trimws(as.character(x)))
  if (has_value(num_cpgs)) cmd <- paste0(cmd, sprintf(' -c %s', trimws(as.character(num_cpgs))))
  if (has_value(num_samplings)) cmd <- paste0(cmd, sprintf(' -s %s', trimws(as.character(num_samplings))))

  run_cmd(cmd, log_fun)

  resolve_sturgeon_v2_output(output_file)
}

#' Return the path of the result file sturgeon-v2 wrote for `-o output_path`:
#' the path itself if it is a file (earlier sturgeon-v2 versions), or
#' `output_path/result.txt` if it is a directory (the public release).
resolve_sturgeon_v2_output <- function(output_path) {
  if (dir.exists(output_path)) {
    result_txt <- file.path(output_path, "result.txt")
    if (!file.exists(result_txt)) {
      stop("sturgeon-v2 wrote an output directory without result.txt: ", output_path)
    }
    return(result_txt)
  }
  if (!file.exists(output_path)) {
    stop("sturgeon-v2 did not produce the expected result file: ", output_path)
  }
  output_path
}

# ---- shared result parsing ---------------------------------------------------
#
# Confirmed against real output from both models: v1 and v2 write the same
# *shape* of CSV -- one header row of class names plus one metadata column,
# one data row of scores plus one metadata value -- just with the metadata
# column in a different place and named differently:
#   v1: "number_probes,<class1>,<class2>,...\n19404,4.05e-10,...\n"
#   v2: "<class1>,<class2>,...,probes\n0.0,0.00067,...,16178\n"
# One shared parser handles both; only the metadata column name differs.

#' Parse a wide-format Sturgeon result file (one column per class, plus one
#' metadata column) into a tidy class/score data.frame.
#'
#' The separator is detected from the header line: the public sturgeon-v2
#' release writes tab-separated `result.txt`, everything else is
#' comma-separated.
#'
#' @param metadata_column the one column to exclude (not a class) -- v1 uses
#'   "number_probes", v2 uses "probes"
#' @param average_rows if TRUE, more than one data row is allowed and the
#'   score per class is the mean over rows. The public sturgeon-v2 writes one
#'   row per random sampling (`-s`), and its own summary is that mean. If
#'   FALSE (v1), anything other than exactly one row is an error.
parse_wide_sturgeon_csv <- function(result_csv, metadata_column, average_rows = FALSE) {
  header <- tryCatch(readLines(result_csv, n = 1, warn = FALSE),
                     error = function(e) stop("Couldn't read ", result_csv, ": ", conditionMessage(e)))
  sep <- if (length(header) == 1 && grepl("\t", header, fixed = TRUE)) "\t" else ","

  raw <- tryCatch(
    read.table(result_csv, sep = sep, header = TRUE, stringsAsFactors = FALSE,
               check.names = FALSE, quote = "\"", comment.char = ""),
    error = function(e) stop("Couldn't read ", result_csv, " as a table: ", conditionMessage(e))
  )
  if (nrow(raw) < 1) stop("No data rows in ", result_csv)
  if (nrow(raw) > 1 && !isTRUE(average_rows)) {
    stop(result_csv, " has ", nrow(raw), " data rows -- expected exactly one (one sample per iteration). ",
         "If multi-row output is now expected, update parse_wide_sturgeon_csv().")
  }

  class_cols <- setdiff(colnames(raw), metadata_column)
  if (length(class_cols) == 0) {
    stop(result_csv, " has no class columns after excluding '", metadata_column,
         "' -- columns found: ", paste(colnames(raw), collapse = ", "))
  }

  scores <- if (nrow(raw) == 1) {
    as.character(unlist(raw[1, class_cols], use.names = FALSE))
  } else {
    as.character(vapply(class_cols, function(cl) mean(as.numeric(raw[[cl]])), numeric(1)))
  }

  data.frame(class = class_cols, score = scores, stringsAsFactors = FALSE)
}

#' Parse one iteration's Sturgeon v1 result CSV into a tidy two-column
#' data.frame: class, score.
parse_sturgeon_result <- function(result_csv) {
  parse_wide_sturgeon_csv(result_csv, metadata_column = "number_probes")
}

#' Parse one iteration's sturgeon-v2 result file the same way.
parse_sturgeon_v2_result <- function(result_file) {
  parse_wide_sturgeon_csv(result_file, metadata_column = "probes", average_rows = TRUE)
}

# ---- progress tracking (data only -- no plotting side effects here) --------

#' Fold one iteration's result into the running wide-format progress table.
#' This is deliberately *pure data* -- plotting is a separate function
#' (plot_confidence_over_time) so the app can re-render the plot on demand
#' without re-running the pipeline.
#'
#' @param progress existing wide data.frame (class, iteration_1, iteration_2, ...)
#'                 or an empty data.frame() for the first call
#' @param iteration_result data.frame(class, score) from parse_sturgeon_result()
#' @param iteration_label column name to use for this iteration, e.g. "iteration_3"
#' @param timestamp POSIXct time to record for this iteration (defaults to now)
update_progress <- function(progress, iteration_result, iteration_label, timestamp = Sys.time()) {
  iteration_result$score <- as.character(iteration_result$score)
  tm <- format(as.POSIXct(timestamp), format = "%H:%M")
  iteration_result <- rbind(iteration_result, data.frame(class = "TIME", score = tm))
  colnames(iteration_result)[2] <- iteration_label

  if (is.null(progress) || nrow(progress) == 0) {
    return(iteration_result)
  }
  merge(progress, iteration_result, by = "class", all = TRUE)
}

# ---- plotting ----------------------------------------------------------------

#' Deterministic fallback color for a class name that isn't in
#' color_translation.rds -- e.g. every Sturgeon v2 class, since v2's
#' abbreviated names (e.g. "A_IDH_HG") don't match v1's full taxonomy names
#' (e.g. "Glioma IDH - A IDH - HG") that color_translation.rds is keyed by.
#' Same class name always gets the same color (stable across iterations and
#' across the v1/v2 plots), and different classes get visually distinct
#' colors, unlike a flat grey fallback. v2 classes listed in the optional
#' class mapping file (see load_class_mapping()) get their v1 color instead.
fallback_class_color <- function(cls) {
  hue <- sum(utf8ToInt(cls)) %% 360
  grDevices::hcl(h = hue, c = 65, l = 50)
}

#' Read the optional v2 -> v1 class-name mapping: a CSV with columns
#' `v2_class` and `v1_class`. Returns a named character vector
#' (names = v2 class, values = v1 class), or NULL if `path` is blank.
#' Only used to pick plot colors -- class names shown in plots and tables
#' are never changed.
load_class_mapping <- function(path) {
  if (is.null(path) || !nzchar(path)) return(NULL)
  df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, strip.white = TRUE)
  missing_cols <- setdiff(c("v2_class", "v1_class"), colnames(df))
  if (length(missing_cols) > 0) {
    stop(path, " is missing column(s): ", paste(missing_cols, collapse = ", "),
         " -- expected a CSV with columns v2_class,v1_class")
  }
  df <- df[nzchar(df$v2_class) & nzchar(df$v1_class), , drop = FALSE]
  dups <- unique(df$v2_class[duplicated(df$v2_class)])
  if (length(dups) > 0) {
    stop(path, " maps these v2 classes more than once: ", paste(dups, collapse = ", "))
  }
  stats::setNames(df$v1_class, df$v2_class)
}

#' Plot color for one class: its own entry in color_translation if it has
#' one, else the color of the v1 class it maps to (class_mapping), else the
#' deterministic fallback color.
class_color <- function(cls, color_translation, class_mapping = NULL) {
  lookup <- function(name) {
    # Look up by name membership first, not a bare `[[name]]` -- if
    # color_translation is a plain named vector (not a list), `[[` on a
    # missing name throws "subscript out of bounds" instead of returning
    # NULL like a list would, which crashed this plot the first time a
    # class name wasn't present in color_translation.rds.
    if (is.null(color_translation) || is.null(name) || is.na(name) ||
        !(name %in% names(color_translation))) {
      return(NULL)
    }
    clr <- color_translation[[name]]
    if (is.null(clr) || is.na(clr) || !nzchar(clr)) NULL else clr
  }
  clr <- lookup(cls)
  if (is.null(clr) && !is.null(class_mapping) && cls %in% names(class_mapping)) {
    clr <- lookup(class_mapping[[cls]])
  }
  if (is.null(clr)) clr <- fallback_class_color(cls)
  clr
}

#' Render the confidence-over-time plot (base R graphics, matches the
#' original look) for the current progress table. Call this inside a Shiny
#' renderPlot(); it draws to whatever graphics device is active.
#' @param class_mapping optional named vector from load_class_mapping(), so
#'   v2 classes are drawn in the color of their v1 counterpart
plot_confidence_over_time <- function(progress, color_translation, plot_title = "Confidence over time",
                                      class_mapping = NULL) {
  # Guard against a malformed/unexpected progress table -- e.g. one restored
  # by resume from a classifier_progress_v1/v2.txt written in an older
  # format (space-separated, pre-dating the write.csv fix), which read.csv()
  # misparses as a single oddly-named column with no "class" column at all.
  # Previously this fell through to nrow()/seq_len() on that broken shape
  # and threw "argument must be coercible to non-negative integer" instead
  # of degrading to the normal "waiting for first iteration" placeholder.
  if (is.null(progress) || !is.data.frame(progress) || nrow(progress) == 0 ||
      !("class" %in% colnames(progress)) || ncol(progress) < 2) {
    plot.new()
    title(main = paste0(plot_title, " -- waiting for first iteration..."))
    return(invisible(NULL))
  }

  mgd <- progress
  n_iter <- ncol(mgd) - 1  # minus the "class" column
  time_row <- mgd[mgd$class == "TIME", , drop = FALSE]
  data_rows <- mgd[mgd$class != "TIME", , drop = FALSE]

  plot('', xlim = c(0, max(n_iter, 5) + 1), ylim = c(0, 1),
       main = plot_title, xlab = "Iteration", ylab = "Confidence",
       xaxt = "n")

  for (i in seq_len(nrow(data_rows))) {
    cls <- data_rows[i, "class"]
    clr <- class_color(cls, color_translation, class_mapping)
    y <-suppressWarnings(as.numeric(data_rows[i, 2:ncol(data_rows)]))
    lines(x = seq_len(n_iter), y = y, col = clr, lwd = 2)
    last_y <- y[length(y)]
    if (!is.na(last_y) && last_y > 0.5) {
      text(x = n_iter + 0.3, y = last_y, labels = cls, adj = 0, cex = 0.8)
    }
  }

  # Cap the number of x-axis labels regardless of how many iterations have
  # piled up (e.g. an overnight run) -- one tick per iteration is fine for a
  # typical few-dozen-iteration case, but at a few hundred+ iterations the
  # labels start overlapping into unreadable clutter, and this plot redraws
  # on every new iteration, so more labels also means every future redraw
  # gets progressively slower. All data points still plot at full
  # resolution via lines() above; only the axis labeling is thinned.
  max_axis_labels <- 20
  if (n_iter <= max_axis_labels) {
    label_idx <- seq_len(n_iter)
  } else {
    step <- ceiling(n_iter / max_axis_labels)
    label_idx <- unique(c(seq(1, n_iter, by = step), n_iter))
  }

  axis(side = 1, at = label_idx, labels = label_idx)
  if (nrow(time_row) == 1) {
    axis(side = 1, at = label_idx,
         labels = as.character(time_row[1, label_idx + 1]),  # +1: col 1 is "class"
         line = 1, tick = FALSE, cex.axis = 0.7)
  }
  abline(h = 0.95, col = "red", lty = 2)
  abline(h = 0.80, col = "orange", lty = 2)
}

#' Latest score per class, sorted descending -- used for the "current top
#' calls" table in the UI.
latest_scores_table <- function(progress) {
  # Same defensive guard as plot_confidence_over_time() -- a progress table
  # missing "class" (e.g. from misparsing an old-format resume file) used to
  # reach `c("class", colnames(progress)[last_col])`, which references a
  # column that doesn't exist and throws "undefined columns selected"
  # instead of just showing an empty table.
  if (is.null(progress) || !is.data.frame(progress) || nrow(progress) == 0 ||
      !("class" %in% colnames(progress)) || ncol(progress) < 2) {
    return(data.frame(class = character(), score = character()))
  }
  last_col <- ncol(progress)
  tbl <- progress[progress$class != "TIME", c("class", colnames(progress)[last_col]), drop = FALSE]
  colnames(tbl) <- c("class", "score")
  tbl$score <- suppressWarnings(as.numeric(tbl$score))
  tbl[order(-tbl$score), ]
}

# ---- CNV plotting (DNAcopy), adapted from plot_cnv_from_bam_DNAcopy --------

#' Generate a CNV plot from a BAM file using QDNAseq bin counts + DNAcopy
#' segmentation. Same algorithm as the original plot_cnv_from_bam_DNAcopy(),
#' but every reference file is an explicit argument instead of a global.
plot_cnv_from_bam <- function(bam_file, bins, cnvplot_reffile, relevant_genes_cnvplot,
                               tmp_bed = tempfile(fileext = ".bed"),
                               makeplot = TRUE, lines_only = FALSE, binsize = 1e6) {
  readCounts <- binReadCounts(bins, bamfiles = bam_file, isNotPassingQualityControls = NA,
                               minMapq = 2, isDuplicate = NA, isSecondaryAlignment = NA)
  exportBins(readCounts, file = tmp_bed, format = "bed", filter = FALSE, logTransform = FALSE)
  on.exit(unlink(tmp_bed), add = TRUE)

  reference <- read.table(cnvplot_reffile, skip = 1)
  colnames(reference) <- c("chrom", "start", "end", "name", "coverage", "orientation")
  sample <- read.table(tmp_bed, skip = 1)
  colnames(sample) <- c("chrom", "start", "end", "name", "coverage", "orientation")
  total_reads_title <- sum(sample$coverage)

  sample$coverage <- as.numeric(sample$coverage) + 0.001
  sample$refcov <- reference$coverage
  sample <- sample[sample$chrom != "Y" & sample$chrom != "X", ]
  sample$reltoref <- sample$coverage / sample$refcov
  blacklist <- sample[sample$refcov < 500, ]
  sample <- sample[!row.names(sample) %in% row.names(blacklist), ]

  mn <- mean(sample$reltoref)
  sample$logtr <- log(sample$reltoref / mn, 2)
  sample$logtr <- pmin(pmax(sample$logtr, -3), 3)
  sample$chromname <- ifelse(nchar(sample$chrom) == 1, paste0("chr0", sample$chrom), paste0("chr", sample$chrom))

  cna_obj <- CNA(sample$logtr, chrom = sample$chromname, maploc = sample$start + binsize,
                 data.type = "logratio", sampleid = "cnvplot")
  smoothed <- smooth.CNA(cna_obj)
  segmented <- segment(smoothed, verbose = 1)

  point_colors <- data.frame(pointx = seq_along(cna_obj$cnvplot), pointy = cna_obj$cnvplot,
                              pointcolor = "black")
  for (x in seq_len(nrow(segmented$segRows))) {
    rows <- segmented$segRows[x, 1]:segmented$segRows[x, 2]
    if (segmented$output[x, "seg.mean"] > 0.5) point_colors[rows, "pointcolor"] <- "green"
    if (segmented$output[x, "seg.mean"] < -0.5) point_colors[rows, "pointcolor"] <- "blue"
  }

  if (isTRUE(makeplot)) {
    plot("", xlim = c(0, length(cna_obj$cnvplot)), ylim = c(-3, 3),
         main = paste0("CNVs, n reads = ", total_reads_title), xaxt = "n",
         xlab = "genomic position", ylab = "log ratio")
    if (!lines_only) {
      points(x = point_colors$pointx, y = point_colors$pointy, pch = 16, cex = 0.5,
             col = point_colors$pointcolor)
    }
    for (x in seq_len(nrow(segmented$segRows))) {
      segments(x0 = segmented$segRows[x, 1], x1 = segmented$segRows[x, 2],
                y0 = segmented$output[x, "seg.mean"], y1 = segmented$output[x, "seg.mean"],
                col = "red", lwd = 2)
    }

    fr <- data.frame(chrom = smoothed$chrom, loc = smoothed$maploc, val = smoothed$cnvplot)
    nr <- 1
    for (i in seq_len(nrow(fr) - 1)) {
      if (fr$chrom[i] != fr$chrom[i + 1]) {
        abline(v = as.numeric(row.names(fr)[i]) + 0.5, col = "grey")
        chrname <- gsub("chr", "", fr[i, "chrom"])
        ypos <- ifelse(nr %% 2 == 1, 2, 1.8)
        nr <- nr + 1
        text(x = i, adj = c(1, 0), labels = chrname, y = ypos)
      }
    }

    if (!lines_only && !is.null(relevant_genes_cnvplot) && file.exists(relevant_genes_cnvplot)) {
      relevant_genes <- read.table(relevant_genes_cnvplot, header = TRUE)
      fr$start <- fr$loc - binsize
      for (i in seq_len(nrow(relevant_genes))) {
        genestart <- relevant_genes[i, "start"]
        genechrom <- strsplit(relevant_genes[i, "chrom"], split = "hr")[[1]][2]
        genechrom <- ifelse(nchar(genechrom) == 1, paste0("chr0", genechrom), paste0("chr", genechrom))
        binline <- fr[fr$chrom == genechrom & fr$start < genestart & (fr$start + binsize * 2) > genestart, ]
        if (nrow(binline) == 0) next
        ycoord <- binline$val
        text(x = as.numeric(row.names(binline)), y = ycoord - 0.1 + relevant_genes[i, "yoffset"],
             cex = 0.5, srt = -90, adj = c(0.5, 1), labels = relevant_genes[i, "name"], pos = 1)
        points(x = as.numeric(row.names(binline)), y = ycoord, col = "red")
      }
    }
  }

  invisible(list(pointdata = point_colors, segdata = segmented))
}

#' Merge every merged_srt.bam produced so far (across all iteration folders)
#' into one BAM and CNV-plot it. Used for the periodic "every N iterations"
#' CNV snapshot.
merge_and_plot_cnv <- function(iteration_bams, cnv_dir, cnv_tag, samtools_bin,
                                bins, cnvplot_reffile, relevant_genes_cnvplot,
                                threads = 4, log_fun = message) {
  if (!dir.exists(cnv_dir)) dir.create(cnv_dir, recursive = TRUE)
  merged <- file.path(cnv_dir, paste0("cnv_", cnv_tag, ".bam"))
  run_cmd(sprintf('%s merge -f -@ %d -o "%s" %s',
                   samtools_bin, threads, merged, paste(shQuote(iteration_bams), collapse = " ")), log_fun)
  plot_cnv_from_bam(merged, bins = bins, cnvplot_reffile = cnvplot_reffile,
                     relevant_genes_cnvplot = relevant_genes_cnvplot)
}

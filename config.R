# =============================================================================
# config.R -- load/save helpers for the app's YAML configuration.
# =============================================================================

library(yaml)

#' Fields in the config that represent filesystem paths -- expanded so a "~"
#' typed into the Setup tab (or written by hand into a saved config.yml)
#' resolves correctly. This matters because pipeline commands are built with
#' each path wrapped in double quotes for `system()` (e.g. `-i "~/foo"`), and
#' shells do NOT tilde-expand inside double quotes -- an un-expanded "~"
#' gets passed to the external tool literally (as a path component named
#' "~") and fails to resolve, instead of falling back to $HOME. Expanding
#' once here means every path field is safe regardless of whether it reached
#' the app via the Setup tab or a hand-edited YAML file.
PATH_FIELDS <- c("bam_pass_folder", "output_folder", "samtools_bin", "modkit_bin",
                  "sturgeon_bin", "probes_file", "model_file_v1", "bins_rds",
                  "cnvplot_reffile", "relevant_genes_bed", "color_translation_rds",
                  "class_mapping_v2", "model_file_v2", "probes_file_v2", "sturgeon_v2_bin")

#' path.expand() every path-like field in a config list (see PATH_FIELDS).
#' Bare command names with no "~" (e.g. "samtools") pass through unchanged.
normalize_cfg_paths <- function(cfg) {
  for (field in PATH_FIELDS) {
    if (!is.null(cfg[[field]]) && nzchar(cfg[[field]])) {
      cfg[[field]] <- path.expand(cfg[[field]])
    }
  }
  if (!is.null(cfg$extra_path_dirs) && length(cfg$extra_path_dirs) > 0) {
    cfg$extra_path_dirs <- vapply(cfg$extra_path_dirs,
                                   function(d) if (nzchar(d)) path.expand(d) else d,
                                   character(1))
  }
  cfg
}

#' Load config in three layers, later ones overriding earlier ones:
#'   1. config_default.yml -- the bundled template
#'   2. config_local.yml   -- written by install.sh with this machine's tool
#'      and reference-file paths (not in git; absent on a manual install)
#'   3. <output_folder>/config.yml, if `saved_path` is given and exists, so
#'      each case keeps its own saved settings
load_config <- function(app_dir, saved_path = NULL,
                        local_path = file.path(app_dir, "config_local.yml")) {
  default_path <- file.path(app_dir, "config_default.yml")
  cfg <- yaml::read_yaml(default_path)

  if (!is.null(local_path) && file.exists(local_path)) {
    local_cfg <- yaml::read_yaml(local_path)
    if (length(local_cfg) > 0) cfg[names(local_cfg)] <- local_cfg
  }

  if (!is.null(saved_path) && file.exists(saved_path)) {
    saved <- yaml::read_yaml(saved_path)
    cfg[names(saved)] <- saved
  }
  normalize_cfg_paths(cfg)
}

#' Save the current config (as a named list) to <output_folder>/config.yml.
save_config <- function(cfg, output_folder) {
  if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
  path <- file.path(output_folder, "config.yml")
  yaml::write_yaml(cfg, path)
  path
}

#' Minimal sanity check before starting a run -- catches the most common
#' "forgot to fill in a path" mistakes up front instead of failing deep into
#' the first iteration during a live case.
validate_config <- function(cfg) {
  problems <- character()

  required_paths <- c("bam_pass_folder", "output_folder", "probes_file", "model_file_v1",
                       "bins_rds", "cnvplot_reffile", "color_translation_rds")
  for (field in required_paths) {
    val <- cfg[[field]]
    if (is.null(val) || !nzchar(val)) {
      problems <- c(problems, sprintf("'%s' is empty.", field))
    }
  }

  # model_file_v2 / probes_file_v2 / class_mapping_v2 are optional -- only checked for existence
  # if the user actually filled them in (blank means "v2 not used" / "v2
  # reuses the v1 probes file").
  file_fields <- c("probes_file", "model_file_v1", "bins_rds", "cnvplot_reffile",
                    "color_translation_rds", "class_mapping_v2", "model_file_v2", "probes_file_v2")
  for (field in file_fields) {
    val <- cfg[[field]]
    if (!is.null(val) && nzchar(val) && !file.exists(val)) {
      problems <- c(problems, sprintf("'%s' points to a file that doesn't exist: %s", field, val))
    }
  }

  if (!is.null(cfg$bam_pass_folder) && nzchar(cfg$bam_pass_folder) && !dir.exists(cfg$bam_pass_folder)) {
    problems <- c(problems, sprintf("bam_pass_folder does not exist: %s", cfg$bam_pass_folder))
  }

  # Catch relative paths early -- a relative output_folder/bam_pass_folder
  # gets created relative to wherever the app process's working directory
  # happens to be (not necessarily your home directory), which has already
  # caused confusing "files ended up nested inside the app's own install
  # folder" surprises. "~" is fine here -- normalize_cfg_paths() already
  # expands it to an absolute path before validate_config() ever sees it.
  for (field in c("bam_pass_folder", "output_folder")) {
    val <- cfg[[field]]
    if (!is.null(val) && nzchar(val) && !startsWith(val, "/")) {
      problems <- c(problems, sprintf(
        "'%s' should be a full absolute path starting with '/' -- got '%s', which is a relative path and would be created relative to the app's working directory.",
        field, val))
    }
  }

  # if v2 is enabled (model_file_v2 set), its executable/mode/file-type still
  # need real values -- but num_cpgs (-c) and num_samplings (-s) are
  # confirmed OPTIONAL in sturgeon-v2's own CLI (left unset, it defaults to
  # "all available CpG sites" / "1 random sampling" itself), so they are
  # deliberately NOT required here. run_sturgeon_v2_predict() only adds
  # -c/-s to the command when the user has actually filled them in.
  if (!is.null(cfg$model_file_v2) && nzchar(cfg$model_file_v2)) {
    v2_required <- c("sturgeon_v2_bin", "sturgeon_v2_mode", "sturgeon_v2_file_type")
    for (field in v2_required) {
      val <- cfg[[field]]
      # A blank Shiny numericInput comes through as R's NA, not "".
      # as.character(NA) is the *string* "NA", which nzchar() sees as
      # non-empty -- so check is.na() explicitly rather than relying on
      # nzchar(as.character(val)) alone, in case a numeric field is ever
      # added back to this list.
      is_missing <- is.null(val) || (length(val) == 1 && is.na(val)) ||
        !nzchar(trimws(as.character(val)))
      if (is_missing) {
        problems <- c(problems, sprintf("'%s' is empty (required because a Sturgeon v2 model is configured).", field))
      }
    }
  }

  problems
}

# ---- installation check ------------------------------------------------------

#' Find an executable the way the app will run it: a value containing "/" is
#' used as a path; a bare name is looked up in `extra_path_dirs` first (the
#' app prepends those to PATH before a run) and then on PATH. Returns the
#' full path, or "" if not found.
resolve_tool <- function(bin, extra_path_dirs = character(0)) {
  if (is.null(bin) || length(bin) != 1 || is.na(bin) || !nzchar(trimws(bin))) return("")
  bin <- path.expand(trimws(bin))
  if (grepl("/", bin, fixed = TRUE)) {
    ok <- file.exists(bin) && !dir.exists(bin) && file.access(bin, 1) == 0
    return(if (ok) bin else "")
  }
  extra_path_dirs <- extra_path_dirs[!is.na(extra_path_dirs) & nzchar(extra_path_dirs)]
  for (d in path.expand(extra_path_dirs)) {
    candidate <- file.path(d, bin)
    if (file.exists(candidate) && !dir.exists(candidate) && file.access(candidate, 1) == 0) return(candidate)
  }
  unname(Sys.which(bin))
}

#' Check that every tool, reference file and folder the configuration points
#' to can be found, before a run is started. Returns a data frame with one
#' row per item and a `status` of "ok", "problem" (the run will fail),
#' "warning" (worth a look, but a run can start) or "not used" (an optional
#' item left blank).
check_setup <- function(cfg) {
  rows <- list()
  add <- function(item, status, detail) {
    rows[[length(rows) + 1]] <<- data.frame(item = item, status = status, detail = detail,
                                             stringsAsFactors = FALSE)
  }
  val <- function(field) {
    v <- cfg[[field]]
    if (is.null(v) || length(v) == 0 || (length(v) == 1 && is.na(v))) return("")
    trimws(as.character(v[[1]]))
  }
  v2_used <- nzchar(val("model_file_v2"))
  extra <- if (is.null(cfg$extra_path_dirs)) character(0) else unlist(cfg$extra_path_dirs)

  for (d in extra[!is.na(extra) & nzchar(extra)]) {
    if (!dir.exists(path.expand(d))) add(paste("Extra PATH directory", d), "warning", "folder does not exist")
  }

  tools <- list(
    list("samtools", "samtools_bin", TRUE),
    list("modkit", "modkit_bin", TRUE),
    list("sturgeon (v1)", "sturgeon_bin", TRUE),
    list("sturgeon-v2", "sturgeon_v2_bin", v2_used)
  )
  for (t in tools) {
    label <- t[[1]]; bin <- val(t[[2]])
    if (!isTRUE(t[[3]])) { add(label, "not used", "no Sturgeon v2 model set"); next }
    if (!nzchar(bin)) { add(label, "problem", sprintf("'%s' is empty", t[[2]])); next }
    found <- resolve_tool(bin, extra)
    if (nzchar(found)) add(label, "ok", found)
    else add(label, "problem", sprintf("not found: %s (set the full path, or add its folder to 'Extra PATH directories')", bin))
  }

  files <- list(
    list("Probes file", "probes_file", TRUE),
    list("Sturgeon v1 model", "model_file_v1", TRUE),
    list("CNV bins", "bins_rds", TRUE),
    list("CNV reference coverage", "cnvplot_reffile", TRUE),
    list("Class colors", "color_translation_rds", TRUE),
    list("CNV gene labels", "relevant_genes_bed", FALSE),
    list("v2-to-v1 class mapping", "class_mapping_v2", FALSE),
    list("Sturgeon v2 model", "model_file_v2", FALSE),
    list("Probes file for v2", "probes_file_v2", FALSE)
  )
  for (f in files) {
    label <- f[[1]]; path <- path.expand(val(f[[2]]))
    if (!nzchar(path)) {
      if (isTRUE(f[[3]])) add(label, "problem", sprintf("'%s' is empty", f[[2]]))
      else add(label, "not used", "optional, left blank")
    } else if (!file.exists(path) || dir.exists(path)) {
      add(label, "problem", paste("file not found:", path))
    } else {
      add(label, "ok", path)
    }
  }

  bam_pass <- path.expand(val("bam_pass_folder"))
  if (!nzchar(bam_pass)) {
    add("MinKNOW folder", "problem", "'bam_pass_folder' is empty")
  } else if (!dir.exists(bam_pass)) {
    add("MinKNOW folder", "problem", paste("folder does not exist:", bam_pass))
  } else {
    add("MinKNOW folder", "ok", bam_pass)
    watch <- describe_watch_folders(bam_pass, pad_barcode(val("barcode")), isTRUE(cfg$use_unclassified))
    for (d in watch$folders) {
      if (d %in% watch$missing) {
        add("Watched folder", "warning",
            paste(d, "does not exist yet (MinKNOW creates it when the first reads arrive; check the barcode if it never appears)"))
      } else {
        add("Watched folder", "ok", d)
      }
    }
  }

  out <- path.expand(val("output_folder"))
  if (!nzchar(out)) {
    add("Output folder", "problem", "'output_folder' is empty")
  } else if (!startsWith(out, "/")) {
    add("Output folder", "problem", paste("not a full path starting with '/':", out))
  } else if (dir.exists(out)) {
    if (file.access(out, 2) != 0) add("Output folder", "problem", paste("not writable:", out))
    else if (file.exists(file.path(out, "processed_files.txt")))
      add("Output folder", "warning", paste(out, "already holds a run; Start will resume it"))
    else add("Output folder", "ok", out)
  } else {
    parent <- dirname(out)
    while (!dir.exists(parent) && parent != dirname(parent)) parent <- dirname(parent)
    if (file.access(parent, 2) == 0) add("Output folder", "ok", paste(out, "(will be created)"))
    else add("Output folder", "problem", paste("cannot be created; no write access to", parent))
  }

  do.call(rbind, rows)
}

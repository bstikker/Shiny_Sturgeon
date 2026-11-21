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

#' Load config: prefer <output_folder>/config.yml if it exists (so each case
#' keeps its own saved settings), otherwise fall back to the bundled default.
load_config <- function(app_dir, saved_path = NULL) {
  default_path <- file.path(app_dir, "config_default.yml")
  cfg <- yaml::read_yaml(default_path)

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

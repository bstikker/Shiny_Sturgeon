# =============================================================================
# Sturgeon Live Classifier -- local Shiny GUI
#
# Run with:  shiny::runApp("path/to/sturgeon_shiny_app")
# (or open this file in RStudio and click "Run App")
#
# Replaces the SLURM/Windows/WSL-era scripts with a single local app:
#   - Setup tab: fill in paths once, save them for this case
#   - Live Run tab: Start/Stop watches the MinKNOW output folder, merges new
#     BAMs, runs Sturgeon, updates the confidence-over-time plot, and
#     refreshes a CNV plot every N iterations (or on demand)
#
# Design notes (why it looks the way it does):
#   - The original R script used a blocking `while(TRUE)` loop with
#     Sys.sleep(). In Shiny that would freeze the UI and make Stop
#     unresponsive. This app uses the standard non-blocking Shiny polling
#     pattern (observe + invalidateLater) instead: the UI stays responsive
#     at all times, and Stop takes effect immediately.
#   - Processing (samtools/sturgeon calls) still runs synchronously within
#     each poll tick, which briefly pauses the UI while an iteration is
#     processed. That's a deliberate simplicity/robustness trade-off over
#     fully async processing (which would need the `future`/`promises`
#     packages) -- fewer moving parts to break mid-case. See README.md.
#   - Every iteration and every error is written to the on-screen log AND
#     kept as plain files on disk (classifier_progress.txt, processed
#     file list, per-iteration folders), so nothing is lost if the app
#     window is closed and reopened mid-case.
# =============================================================================

library(shiny)
library(yaml)

app_dir <- getwd()
source(file.path(app_dir, "core_functions.R"))
source(file.path(app_dir, "config.R"))

`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && !nzchar(as.character(a)))) b else a

default_cfg <- load_config(app_dir)

# ---- UI ----------------------------------------------------------------

ui <- fluidPage(
  titlePanel("Sturgeon Live Intraoperative Classifier"),

  tabsetPanel(
    id = "main_tabs",

    tabPanel("Setup",
      br(),
      fluidRow(
        column(6,
          h4("Sequencing run"),
          textInput("bam_pass_folder", "MinKNOW realtime BAM output folder (e.g. .../bam_pass -- contains barcode01/, barcode02/, ..., unclassified/)",
                     value = default_cfg$bam_pass_folder, width = "100%"),
          textInput("barcode", "Barcode", value = default_cfg$barcode),
          checkboxInput("use_unclassified", "Include unclassified reads",
                         value = isTRUE(default_cfg$use_unclassified)),

          h4("Output"),
          textInput("output_folder", "Output folder for this case (will be created)",
                     value = default_cfg$output_folder, width = "100%"),

          h4("Timing"),
          numericInput("poll_interval_sec", "Check for new files every (seconds)",
                        value = default_cfg$poll_interval_sec, min = 5),
          numericInput("cnv_every_n_iterations", "Auto-update CNV plot every N iterations",
                        value = default_cfg$cnv_every_n_iterations, min = 1)
        ),
        column(6,
          h4("Tools"),
          textInput("samtools_bin", "samtools path (or just 'samtools' if on PATH)",
                     value = default_cfg$samtools_bin, width = "100%"),
          textInput("modkit_bin", "modkit path (or just 'modkit' if on PATH)",
                     value = default_cfg$modkit_bin, width = "100%"),
          textInput("sturgeon_bin", "sturgeon path (point at your venv's sturgeon if not on PATH)",
                     value = default_cfg$sturgeon_bin, width = "100%"),
          textInput("extra_path_dirs", "Extra PATH directories (comma-separated, e.g. conda env bin/)",
                     value = paste(default_cfg$extra_path_dirs, collapse = ", "), width = "100%"),

          h4("Reference files"),
          textInput("probes_file", "Probes file (probelocs_chm13.bed)",
                     value = default_cfg$probes_file, width = "100%"),
          textInput("model_file_v1", "Sturgeon v1 model (.zip)",
                     value = default_cfg$model_file_v1, width = "100%"),
          textInput("bins_rds", "CNV bins .rds",
                     value = default_cfg$bins_rds, width = "100%"),
          textInput("cnvplot_reffile", "CNV reference coverage .bed",
                     value = default_cfg$cnvplot_reffile, width = "100%"),
          textInput("relevant_genes_bed", "Relevant genes .bed (optional, CNV plot labels)",
                     value = default_cfg$relevant_genes_bed, width = "100%"),
          textInput("color_translation_rds", "Class-to-color .rds",
                     value = default_cfg$color_translation_rds, width = "100%"),

          h4("Sturgeon v2 (optional)"),
          helpText("Leave the model blank to run v1 only. Fill it in and the app runs",
                    "both models every iteration, with two separate confidence-over-time plots.",
                    "v2 has its own CLI with its own required parameters (below) -- all of them",
                    "become mandatory once a v2 model is set."),
          textInput("model_file_v2", "Sturgeon v2 model file -- optional",
                     value = default_cfg$model_file_v2, width = "100%"),
          textInput("probes_file_v2", "Probes file for v2 -- only if different from v1 above",
                     value = default_cfg$probes_file_v2, width = "100%"),
          textInput("sturgeon_v2_bin", "sturgeon-v2 executable (path, if not on PATH)",
                     value = default_cfg$sturgeon_v2_bin, width = "100%"),
          selectInput("sturgeon_v2_mode", "Mode", choices = c("single", "multi"),
                       selected = default_cfg$sturgeon_v2_mode %||% "single"),
          textInput("sturgeon_v2_num_cpgs", "NUM_CPGS -- optional, a number; leave blank to use all available CpG sites (sturgeon-v2's own default)",
                     value = default_cfg$sturgeon_v2_num_cpgs %||% "", width = "100%"),
          numericInput("sturgeon_v2_num_samplings", "NUM_SAMPLINGS -- optional; leave blank for sturgeon-v2's default (1 random sampling)",
                        value = suppressWarnings(as.numeric(default_cfg$sturgeon_v2_num_samplings))),
          textInput("sturgeon_v2_file_type", "FILE_TYPE flag for -f (best guess: \"bed\" -- verify)",
                     value = default_cfg$sturgeon_v2_file_type, width = "100%")
        )
      ),
      br(),
      fluidRow(
        column(12,
          actionButton("save_cfg_btn", "Save configuration", class = "btn-primary"),
          actionButton("load_cfg_btn", "Load config.yml from output folder"),
          br(), br(),
          verbatimTextOutput("configStatus")
        )
      )
    ),

    tabPanel("Live Run",
      br(),
      fluidRow(
        column(12,
          actionButton("start_btn", "Start monitoring", class = "btn-success btn-lg"),
          actionButton("stop_btn", "Stop", class = "btn-danger btn-lg"),
          actionButton("cnv_now_btn", "Generate CNV plot now"),
          actionButton("reset_btn", "Reset for new case"),
          br(), br(),
          h4(textOutput("statusText")),
          br()
        )
      ),
      fluidRow(
        column(6,
          h4("Confidence over time -- Sturgeon v1"),
          plotOutput("confidencePlotV1", height = "380px"),
          downloadButton("downloadConfidenceV1", "Download plot (PDF)"),
          downloadButton("downloadProgressV1", "Download data (CSV)"),
          br(), br(),
          h5("Current top calls (v1)"),
          tableOutput("scoresTableV1")
        ),
        column(6,
          h4("Confidence over time -- Sturgeon v2"),
          plotOutput("confidencePlotV2", height = "380px"),
          downloadButton("downloadConfidenceV2", "Download plot (PDF)"),
          downloadButton("downloadProgressV2", "Download data (CSV)"),
          br(), br(),
          h5("Current top calls (v2)"),
          tableOutput("scoresTableV2")
        )
      ),
      br(),
      fluidRow(
        column(7,
          h4("CNV plot"),
          helpText("Shared between both models -- CNV is computed straight from the BAMs, independent of which classifier ran."),
          plotOutput("cnvPlot", height = "380px"),
          downloadButton("downloadCNV", "Download CNV plot (PDF)")
        ),
        column(5,
          h4("Log"),
          div(style = "height: 420px; overflow-y: scroll; background: #111; color: #0f0;
                        font-family: monospace; font-size: 12px; padding: 8px;",
              verbatimTextOutput("logOutput"))
        )
      )
    )
  )
)

# ---- SERVER --------------------------------------------------------------

server <- function(input, output, session) {

  rv <- reactiveValues(
    cfg = NULL,
    running = FALSE,
    iteration = 1,
    processed_files = character(0),
    iteration_bams = character(0),
    progress_v1 = data.frame(),
    progress_v2 = data.frame(),
    log_lines = character(0),
    color_translation = NULL,
    bins = NULL,
    cnv_counter = 0,
    latest_cnv_bam = NULL,
    cnv_trigger = 0
  )

  add_log <- function(msg) {
    rv$log_lines <- c(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", msg), rv$log_lines)
  }

  config_status_msg <- reactiveVal("")
  output$configStatus <- renderText(config_status_msg())

  gather_cfg_from_inputs <- function() {
    normalize_cfg_paths(list(
      bam_pass_folder = input$bam_pass_folder,
      barcode = pad_barcode(input$barcode),
      use_unclassified = input$use_unclassified,
      output_folder = input$output_folder,
      poll_interval_sec = input$poll_interval_sec,
      cnv_every_n_iterations = input$cnv_every_n_iterations,
      samtools_bin = input$samtools_bin,
      modkit_bin = input$modkit_bin,
      sturgeon_bin = input$sturgeon_bin,
      extra_path_dirs = trimws(strsplit(input$extra_path_dirs, ",")[[1]]),
      probes_file = input$probes_file,
      model_file_v1 = input$model_file_v1,
      bins_rds = input$bins_rds,
      cnvplot_reffile = input$cnvplot_reffile,
      relevant_genes_bed = input$relevant_genes_bed,
      color_translation_rds = input$color_translation_rds,
      model_file_v2 = input$model_file_v2,
      probes_file_v2 = input$probes_file_v2,
      sturgeon_v2_bin = input$sturgeon_v2_bin,
      sturgeon_v2_mode = input$sturgeon_v2_mode,
      sturgeon_v2_num_cpgs = input$sturgeon_v2_num_cpgs,
      sturgeon_v2_num_samplings = input$sturgeon_v2_num_samplings,
      sturgeon_v2_file_type = input$sturgeon_v2_file_type
    ))
  }

  apply_cfg_to_inputs <- function(cfg) {
    updateTextInput(session, "bam_pass_folder", value = cfg$bam_pass_folder %||% "")
    updateTextInput(session, "barcode", value = cfg$barcode %||% "")
    updateCheckboxInput(session, "use_unclassified", value = isTRUE(cfg$use_unclassified))
    updateTextInput(session, "output_folder", value = cfg$output_folder %||% "")
    updateNumericInput(session, "poll_interval_sec", value = cfg$poll_interval_sec %||% 15)
    updateNumericInput(session, "cnv_every_n_iterations", value = cfg$cnv_every_n_iterations %||% 5)
    updateTextInput(session, "samtools_bin", value = cfg$samtools_bin %||% "samtools")
    updateTextInput(session, "modkit_bin", value = cfg$modkit_bin %||% "modkit")
    updateTextInput(session, "sturgeon_bin", value = cfg$sturgeon_bin %||% "sturgeon")
    updateTextInput(session, "extra_path_dirs", value = paste(cfg$extra_path_dirs %||% character(0), collapse = ", "))
    updateTextInput(session, "probes_file", value = cfg$probes_file %||% "")
    updateTextInput(session, "model_file_v1", value = cfg$model_file_v1 %||% "")
    updateTextInput(session, "bins_rds", value = cfg$bins_rds %||% "")
    updateTextInput(session, "cnvplot_reffile", value = cfg$cnvplot_reffile %||% "")
    updateTextInput(session, "relevant_genes_bed", value = cfg$relevant_genes_bed %||% "")
    updateTextInput(session, "color_translation_rds", value = cfg$color_translation_rds %||% "")
    updateTextInput(session, "model_file_v2", value = cfg$model_file_v2 %||% "")
    updateTextInput(session, "probes_file_v2", value = cfg$probes_file_v2 %||% "")
    updateTextInput(session, "sturgeon_v2_bin", value = cfg$sturgeon_v2_bin %||% "sturgeon-v2")
    updateSelectInput(session, "sturgeon_v2_mode", selected = cfg$sturgeon_v2_mode %||% "single")
    updateTextInput(session, "sturgeon_v2_num_cpgs", value = cfg$sturgeon_v2_num_cpgs %||% "")
    updateNumericInput(session, "sturgeon_v2_num_samplings", value = suppressWarnings(as.numeric(cfg$sturgeon_v2_num_samplings)))
    updateTextInput(session, "sturgeon_v2_file_type", value = cfg$sturgeon_v2_file_type %||% "bed")
  }

  # (%||% is defined once at file scope, above, so it's available both here
  # and while building `ui` at app-load time.)

  # ---- Setup tab -------------------------------------------------------

  observeEvent(input$save_cfg_btn, {
    cfg <- gather_cfg_from_inputs()
    problems <- validate_config(cfg)
    if (length(problems) > 0) {
      config_status_msg(paste("Not saved -- fix these first:\n",
                               paste("-", problems, collapse = "\n")))
    } else {
      path <- save_config(cfg, cfg$output_folder)
      config_status_msg(paste("Saved to", path))
    }
  })

  observeEvent(input$load_cfg_btn, {
    path <- file.path(input$output_folder, "config.yml")
    if (file.exists(path)) {
      cfg <- yaml::read_yaml(path)
      apply_cfg_to_inputs(cfg)
      config_status_msg(paste("Loaded", path))
    } else {
      config_status_msg(paste("No config.yml found in", input$output_folder))
    }
  })

  # ---- Live Run tab: start / stop / reset ------------------------------

  observeEvent(input$start_btn, {
    cfg <- gather_cfg_from_inputs()
    problems <- validate_config(cfg)
    if (length(problems) > 0) {
      showModal(modalDialog(
        title = "Can't start -- fix configuration first",
        tags$ul(lapply(problems, tags$li)),
        easyClose = TRUE
      ))
      return(invisible(NULL))
    }

    rv$cfg <- cfg
    save_config(cfg, cfg$output_folder)
    apply_extra_path(cfg$extra_path_dirs)

    tryCatch({
      rv$bins <- readRDS(cfg$bins_rds)
      rv$color_translation <- readRDS(cfg$color_translation_rds)
    }, error = function(e) {
      showModal(modalDialog(title = "Failed to load reference files",
                             conditionMessage(e), easyClose = TRUE))
    })

    # resume support: pick up where a previous run in this output folder left
    # off. This matters beyond just "don't redo work" -- if the browser tab
    # disconnects (laptop sleep, network drop, an overnight run outliving the
    # session) Shiny starts a brand new session with empty reactiveValues,
    # so the confidence-over-time plots and CNV plot would otherwise come
    # back blank even though every iteration's data is safely on disk.
    processed_path <- file.path(cfg$output_folder, "processed_files.txt")
    if (file.exists(processed_path) && length(rv$processed_files) == 0) {
      rv$processed_files <- readLines(processed_path)
      existing_its <- list.files(cfg$output_folder, pattern = "^iteration[0-9]+$")
      n_done <- 0
      if (length(existing_its) > 0) {
        nums <- as.integer(gsub("iteration", "", existing_its))
        n_done <- max(nums, na.rm = TRUE)
        rv$iteration <- n_done + 1

        # rebuild the cumulative CNV bam list (in-memory only, never
        # persisted directly) from each iteration's modkit-adjusted bam, so
        # the next CNV snapshot still covers the whole case, not just
        # whatever happens after this reconnect.
        candidate_bams <- file.path(cfg$output_folder, paste0("iteration", sort(nums)), "merged_srt_adjust.bam")
        rv$iteration_bams <- candidate_bams[file.exists(candidate_bams)]
        rv$cnv_counter <- n_done %% max(cfg$cnv_every_n_iterations, 1)
      }

      # restore the confidence-over-time history for both models, if present.
      # If the file predates the write.csv fix (still old space-separated
      # write.table format) or is otherwise malformed, read.csv() won't
      # error -- it'll just produce a data frame with no "class" column.
      # Detect that here and skip it with a clear log line, rather than
      # handing plot_confidence_over_time()/latest_scores_table() a shape
      # they don't expect.
      load_progress_csv <- function(path, label) {
        if (!file.exists(path)) return(NULL)
        df <- tryCatch(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
                        error = function(e) NULL)
        if (is.null(df) || !("class" %in% colnames(df))) {
          add_log(paste0("Could not restore ", label, " history from ", path,
                          " (unexpected format, likely from an older app version) -- ",
                          "starting that plot fresh; the file will be overwritten by the next iteration."))
          return(NULL)
        }
        df
      }
      restored_v1 <- load_progress_csv(file.path(cfg$output_folder, "classifier_progress_v1.txt"), "v1")
      if (!is.null(restored_v1)) rv$progress_v1 <- restored_v1
      restored_v2 <- load_progress_csv(file.path(cfg$output_folder, "classifier_progress_v2.txt"), "v2")
      if (!is.null(restored_v2)) rv$progress_v2 <- restored_v2

      # restore the most recent CNV snapshot so the CNV plot isn't blank
      # until the next auto-trigger fires
      cnv_dir <- file.path(cfg$output_folder, "CNV_data")
      if (dir.exists(cnv_dir)) {
        existing_cnv <- list.files(cnv_dir, pattern = "^cnv_it[0-9]+\\.bam$", full.names = TRUE)
        if (length(existing_cnv) > 0) {
          cnv_nums <- as.integer(gsub("[^0-9]", "", basename(existing_cnv)))
          rv$latest_cnv_bam <- existing_cnv[which.max(cnv_nums)]
        }
      }

      add_log(paste0("Resumed previous run: ", length(rv$processed_files),
                      " files already processed, starting at iteration ", rv$iteration))
    }

    rv$running <- TRUE
    add_log("Monitoring started.")
  })

  observeEvent(input$stop_btn, {
    rv$running <- FALSE
    add_log("Monitoring stopped.")
  })

  observeEvent(input$reset_btn, {
    rv$running <- FALSE
    rv$iteration <- 1
    rv$processed_files <- character(0)
    rv$iteration_bams <- character(0)
    rv$progress_v1 <- data.frame()
    rv$progress_v2 <- data.frame()
    rv$cnv_counter <- 0
    rv$latest_cnv_bam <- NULL
    rv$cnv_trigger <- 0
    rv$log_lines <- character(0)
    add_log("Reset for a new case. Update Setup tab (output folder etc.) before starting.")
  })

  output$statusText <- renderText({
    if (isTRUE(rv$running)) {
      sprintf("Running -- iteration %d, %d file(s) processed so far",
              rv$iteration, length(rv$processed_files))
    } else {
      "Idle"
    }
  })

  # ---- polling loop (non-blocking) -------------------------------------

  poll_once <- function() {
    cfg <- rv$cfg
    new_bams <- tryCatch(
      find_new_bams(cfg$bam_pass_folder, cfg$barcode, cfg$use_unclassified, rv$processed_files),
      error = function(e) {
        add_log(paste("ERROR scanning for BAM files:", conditionMessage(e)))
        character(0)
      }
    )

    if (length(new_bams) == 0) {
      add_log("No new BAM files, waiting...")
      return(invisible(NULL))
    }

    add_log(sprintf("Found %d new BAM file(s) -- waiting 5s in case a file is still being written...",
                     length(new_bams)))
    Sys.sleep(5)  # settle delay, confirmed still needed on the real setup

    add_log(sprintf("Processing iteration %d", rv$iteration))
    iteration_dir <- file.path(cfg$output_folder, paste0("iteration", rv$iteration))

    # cumulative across the whole case, NOT per-iteration: each iteration's
    # modkit.txt is added here and never removed, so inputtobed below always
    # aggregates every read processed so far -- this is what makes
    # "confidence over time" mean something (climbing as data accumulates)
    # rather than each iteration being an independent, noisy read of just
    # its own small new batch.
    modkit_txt_dir <- file.path(cfg$output_folder, "modkit_extract")
    bed_dir_v1 <- file.path(cfg$output_folder, "bed_v1")

    ok <- TRUE
    adjusted_bam <- NULL
    tryCatch({
      sorted_bam <- merge_and_sort_bams(new_bams, iteration_dir, cfg$samtools_bin, log_fun = add_log)

      # modkit preprocessing on just this iteration's new reads (cheap --
      # never redone on already-processed reads)
      adjusted_bam <- run_modkit_adjust_mods(sorted_bam, file.path(iteration_dir, "merged_srt_adjust.bam"),
                                              cfg$modkit_bin, log_fun = add_log)
      run_modkit_extract_full(adjusted_bam, file.path(modkit_txt_dir, paste0("iteration_", rv$iteration, "_modkit.txt")),
                               cfg$modkit_bin, log_fun = add_log)

      # aggregate ALL modkit.txt so far (cumulative) into one fresh bed, then predict from it
      run_sturgeon_inputtobed(modkit_txt_dir, bed_dir_v1, cfg$probes_file, cfg$sturgeon_bin, log_fun = add_log)

      result_csv_v1 <- run_sturgeon_predict(bed_dir_v1, cfg$model_file_v1,
                                             file.path(iteration_dir, "v1"), cfg$sturgeon_bin, log_fun = add_log)
      parsed_v1 <- parse_sturgeon_result(result_csv_v1)
      rv$progress_v1 <- update_progress(rv$progress_v1, parsed_v1, paste0("iteration_", rv$iteration))
      # write.csv, not write.table's default space-separated format -- class
      # names can contain spaces (e.g. "Glioma IDH - A IDH - HG"), which
      # would misalign columns on a later read.csv() round-trip (used by the
      # resume logic above to restore this plot's history after a dropped
      # browser session).
      write.csv(rv$progress_v1, file = file.path(cfg$output_folder, "classifier_progress_v1.txt"),
                row.names = FALSE)

      if (nzchar(cfg$model_file_v2 %||% "")) {
        # reuse v1's cumulative bed unless a different probes file was given for v2
        if (nzchar(cfg$probes_file_v2 %||% "")) {
          bed_dir_v2 <- file.path(cfg$output_folder, "bed_v2")
          run_sturgeon_inputtobed(modkit_txt_dir, bed_dir_v2, cfg$probes_file_v2, cfg$sturgeon_bin, log_fun = add_log)
        } else {
          bed_dir_v2 <- bed_dir_v1
        }
        result_csv_v2 <- run_sturgeon_v2_predict(
          bed_dir_v2, cfg$model_file_v2,
          output_file = file.path(iteration_dir, "v2", "sturgeon_v2_result.csv"),
          sturgeon_v2_bin = cfg$sturgeon_v2_bin, mode = cfg$sturgeon_v2_mode,
          num_cpgs = cfg$sturgeon_v2_num_cpgs, num_samplings = cfg$sturgeon_v2_num_samplings,
          file_type = cfg$sturgeon_v2_file_type, log_fun = add_log
        )
        parsed_v2 <- parse_sturgeon_v2_result(result_csv_v2)
        rv$progress_v2 <- update_progress(rv$progress_v2, parsed_v2, paste0("iteration_", rv$iteration))
        write.csv(rv$progress_v2, file = file.path(cfg$output_folder, "classifier_progress_v2.txt"),
                  row.names = FALSE)
      }

      add_log(sprintf("Iteration %d complete.", rv$iteration))
    }, error = function(e) {
      ok <<- FALSE
      add_log(sprintf("ERROR in iteration %d: %s -- will retry with next new file(s).",
                       rv$iteration, conditionMessage(e)))
    })

    # mark these files processed regardless of success, so a permanently broken
    # file can't wedge the loop retrying it forever; failures are visible in the log
    rv$processed_files <- c(rv$processed_files, new_bams)
    write(new_bams, file = file.path(cfg$output_folder, "processed_files.txt"), append = TRUE)

    if (isTRUE(ok) && !is.null(adjusted_bam)) {
      # CNV uses the modkit-adjusted bam (matches the real script moving
      # merged_srt_adjust.bam, not the plain sorted one, into its cumulative pool)
      rv$iteration_bams <- c(rv$iteration_bams, adjusted_bam)
      rv$cnv_counter <- rv$cnv_counter + 1
      rv$iteration <- rv$iteration + 1

      if (rv$cnv_counter >= cfg$cnv_every_n_iterations) {
        rv$cnv_counter <- 0
        generate_cnv_snapshot()
      }
    }
  }

  generate_cnv_snapshot <- function() {
    req(rv$cfg)
    if (length(rv$iteration_bams) == 0) {
      add_log("No processed iterations yet -- nothing to build a CNV plot from.")
      return(invisible(NULL))
    }
    add_log("Building CNV snapshot (merging iterations so far)...")
    cfg <- rv$cfg
    tryCatch({
      cnv_dir <- file.path(cfg$output_folder, "CNV_data")
      if (!dir.exists(cnv_dir)) dir.create(cnv_dir, recursive = TRUE)
      merged <- file.path(cnv_dir, paste0("cnv_it", rv$iteration, ".bam"))
      run_cmd(sprintf('%s merge -f -@ 4 -o "%s" %s',
                       cfg$samtools_bin, merged, paste(shQuote(rv$iteration_bams), collapse = " ")),
              add_log)
      rv$latest_cnv_bam <- merged
      rv$cnv_trigger <- rv$cnv_trigger + 1
      add_log("CNV plot updated.")
    }, error = function(e) {
      add_log(paste("ERROR building CNV plot:", conditionMessage(e)))
    })
  }

  observeEvent(input$cnv_now_btn, {
    generate_cnv_snapshot()
  })

  observe({
    if (!isTRUE(rv$running)) return(NULL)
    invalidateLater(max(rv$cfg$poll_interval_sec, 5) * 1000, session)
    isolate(poll_once())
  })

  # ---- outputs -----------------------------------------------------------

  output$confidencePlotV1 <- renderPlot({
    req(rv$color_translation)
    plot_confidence_over_time(rv$progress_v1, rv$color_translation, plot_title = "Confidence over time -- Sturgeon v1")
  })

  output$scoresTableV1 <- renderTable({
    latest_scores_table(rv$progress_v1)
  }, digits = 3)

  output$confidencePlotV2 <- renderPlot({
    if (is.null(rv$cfg) || !nzchar(rv$cfg$model_file_v2 %||% "")) {
      plot.new()
      title(main = "Sturgeon v2 not configured (see Setup tab)")
      return(invisible(NULL))
    }
    req(rv$color_translation)
    plot_confidence_over_time(rv$progress_v2, rv$color_translation, plot_title = "Confidence over time -- Sturgeon v2")
  })

  output$scoresTableV2 <- renderTable({
    latest_scores_table(rv$progress_v2)
  }, digits = 3)

  output$cnvPlot <- renderPlot({
    rv$cnv_trigger
    req(rv$latest_cnv_bam, rv$bins)
    plot_cnv_from_bam(rv$latest_cnv_bam, bins = rv$bins, cnvplot_reffile = rv$cfg$cnvplot_reffile,
                       relevant_genes_cnvplot = rv$cfg$relevant_genes_bed)
  })

  output$logOutput <- renderText({
    paste(rv$log_lines, collapse = "\n")
  })

  output$downloadConfidenceV1 <- downloadHandler(
    filename = function() paste0("confidence_over_time_v1_", format(Sys.time(), "%Y%m%d_%H%M"), ".pdf"),
    content = function(file) {
      pdf(file)
      plot_confidence_over_time(rv$progress_v1, rv$color_translation, plot_title = "Confidence over time -- Sturgeon v1")
      dev.off()
    }
  )

  output$downloadProgressV1 <- downloadHandler(
    filename = function() paste0("classifier_progress_v1_", format(Sys.time(), "%Y%m%d_%H%M"), ".csv"),
    content = function(file) write.csv(rv$progress_v1, file, row.names = FALSE)
  )

  output$downloadConfidenceV2 <- downloadHandler(
    filename = function() paste0("confidence_over_time_v2_", format(Sys.time(), "%Y%m%d_%H%M"), ".pdf"),
    content = function(file) {
      req(rv$cfg, cancelOutput = TRUE)
      validate(need(nzchar(rv$cfg$model_file_v2 %||% ""), "Sturgeon v2 is not configured for this case."))
      pdf(file)
      plot_confidence_over_time(rv$progress_v2, rv$color_translation, plot_title = "Confidence over time -- Sturgeon v2")
      dev.off()
    }
  )

  output$downloadProgressV2 <- downloadHandler(
    filename = function() paste0("classifier_progress_v2_", format(Sys.time(), "%Y%m%d_%H%M"), ".csv"),
    content = function(file) write.csv(rv$progress_v2, file, row.names = FALSE)
  )

  output$downloadCNV <- downloadHandler(
    filename = function() paste0("cnv_plot_", format(Sys.time(), "%Y%m%d_%H%M"), ".pdf"),
    content = function(file) {
      req(rv$latest_cnv_bam, rv$bins)
      pdf(file, useDingbats = FALSE)
      plot_cnv_from_bam(rv$latest_cnv_bam, bins = rv$bins, cnvplot_reffile = rv$cfg$cnvplot_reffile,
                         relevant_genes_cnvplot = rv$cfg$relevant_genes_bed)
      dev.off()
    }
  )
}

shinyApp(ui = ui, server = server)

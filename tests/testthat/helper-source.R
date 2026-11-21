# Load the app's functions (not the Shiny app itself) into the test environment.
repo_root <- normalizePath(file.path(testthat::test_path(), "..", ".."))
fixture <- function(...) file.path(repo_root, "tests", "fixtures", ...)

suppressPackageStartupMessages({
  source(file.path(repo_root, "core_functions.R"), local = TRUE)
  source(file.path(repo_root, "config.R"), local = TRUE)
})

# Write an executable shell script that stands in for an external tool
# (samtools, modkit, sturgeon, sturgeon-v2). It logs its arguments to
# `<dir>/<name>.args` and then runs `body`.
fake_tool <- function(dir, name, body = "exit 0") {
  path <- file.path(dir, name)
  writeLines(c("#!/bin/sh",
               sprintf('printf "%%s\\n" "$*" >> "%s"', file.path(dir, paste0(name, ".args"))),
               body), path)
  Sys.chmod(path, "755")
  path
}

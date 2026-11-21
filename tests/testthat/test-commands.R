# The external tools are replaced by small shell scripts (fake_tool()) that
# record their arguments and write the files the real tool would write.
# This checks the commands the app builds and how it handles their output.
quiet <- function(...) invisible(NULL)

test_that("merge_and_sort_bams calls samtools merge then sort", {
  bin <- tempfile("bin_"); dir.create(bin)
  samtools <- fake_tool(bin, "samtools")
  out <- merge_and_sort_bams(c("/data/a b.bam", "/data/c.bam"), file.path(bin, "it1"), samtools, log_fun = quiet)
  args <- readLines(file.path(bin, "samtools.args"))
  expect_match(args[1], "^merge -f -@ 4 -o .*merged.bam /data/a b.bam /data/c.bam$")
  expect_match(args[2], "^sort -@ 4 -o .*merged_srt.bam .*merged.bam$")
  expect_equal(basename(out), "merged_srt.bam")
})

test_that("modkit steps pass the right arguments and check their outputs", {
  bin <- tempfile("bin_"); dir.create(bin)
  # writes its last argument, like modkit writing its output file
  modkit <- fake_tool(bin, "modkit", 'for last; do :; done; touch "$last"')
  out_bam <- run_modkit_adjust_mods("in.bam", file.path(bin, "adj.bam"), modkit, log_fun = quiet)
  out_txt <- run_modkit_extract_full(out_bam, file.path(bin, "extract", "it1.txt"), modkit, log_fun = quiet)
  args <- readLines(file.path(bin, "modkit.args"))
  expect_match(args[1], "^adjust-mods --convert h m in.bam ")
  expect_match(args[2], "^extract full -t 10 ")
  expect_true(file.exists(out_txt))

  broken <- fake_tool(bin, "broken_modkit")
  expect_error(run_modkit_adjust_mods("in.bam", file.path(bin, "none.bam"), broken, log_fun = quiet),
               "did not produce")
})

test_that("sturgeon inputtobed recreates the bed folder and predict finds the result", {
  bin <- tempfile("bin_"); dir.create(bin)
  sturgeon <- fake_tool(bin, "sturgeon", paste(
    'case "$1" in',
    '  inputtobed) out=$(echo "$*" | sed -E "s/.* -o ([^ ]+) .*/\\1/"); touch "$out/merged_probes_methyl_calls.bed" ;;',
    sprintf('  predict) out=$(echo "$*" | sed -E "s/.* -o ([^ ]+) .*/\\1/"); cp "%s" "$out/merged_probes_methyl_calls_general.csv" ;;',
            fixture("sturgeon_v1_result.csv")),
    'esac', sep = "\n"))
  bed_dir <- file.path(bin, "bed_v1")
  dir.create(bed_dir); file.create(file.path(bed_dir, "stale.bed"))

  run_sturgeon_inputtobed(file.path(bin, "modkit"), bed_dir, "/ref/probes.bed", sturgeon, log_fun = quiet)
  expect_equal(list.files(bed_dir), "merged_probes_methyl_calls.bed")

  csv <- run_sturgeon_predict(bed_dir, "/ref/general.zip", file.path(bin, "v1"), sturgeon, log_fun = quiet)
  expect_equal(basename(csv), "merged_probes_methyl_calls_general.csv")
  expect_equal(nrow(parse_sturgeon_result(csv)), 5)

  args <- readLines(file.path(bin, "sturgeon.args"))
  expect_match(args[1], "^inputtobed --margin 50 -i .* -s modkit --probes-file /ref/probes.bed$")
  expect_match(args[2], "^predict -p --i .* -m /ref/general.zip$")
})

v2_setup <- function(writes_directory) {
  bin <- tempfile("bin_"); dir.create(bin)
  bed_dir <- file.path(bin, "bed"); dir.create(bed_dir)
  file.create(file.path(bed_dir, "merged_probes_methyl_calls.bed"))
  # find the value after -o and write a result there, either as a file
  # (earlier sturgeon-v2) or as <dir>/result.txt (public sturgeon-v2)
  write_result <- if (writes_directory) {
    sprintf('mkdir -p "$out"; cp "%s" "$out/result.txt"', fixture("sturgeon_v2_result.txt"))
  } else {
    sprintf('cp "%s" "$out"', fixture("sturgeon_v2_result_legacy.csv"))
  }
  v2 <- fake_tool(bin, "sturgeon-v2", paste(
    'while [ $# -gt 0 ]; do [ "$1" = "-o" ] && out="$2"; shift; done', write_result, sep = "\n"))
  list(bin = bin, bed_dir = bed_dir, v2 = v2)
}

test_that("sturgeon-v2 works when it writes a result file at -o (earlier versions)", {
  s <- v2_setup(writes_directory = FALSE)
  out <- run_sturgeon_v2_predict(s$bed_dir, "/ref/cns-v2.zip", file.path(s$bin, "v2", "sturgeon_v2_result.csv"),
                                 s$v2, file_type = "bed", log_fun = quiet)
  expect_equal(basename(out), "sturgeon_v2_result.csv")
  expect_equal(as.numeric(parse_sturgeon_v2_result(out)$score[1]), 0.88)
  args <- readLines(file.path(s$bin, "sturgeon-v2.args"))
  expect_match(args, "^single -i .*merged_probes_methyl_calls.bed -m /ref/cns-v2.zip -o .* -f bed$")
})

test_that("sturgeon-v2 works when it writes result.txt into an -o directory (public release)", {
  s <- v2_setup(writes_directory = TRUE)
  out <- run_sturgeon_v2_predict(s$bed_dir, "/ref/cns-v2.zip", file.path(s$bin, "v2", "sturgeon_v2_result.csv"),
                                 s$v2, file_type = "bed", log_fun = quiet)
  expect_equal(basename(out), "result.txt")
  expect_equal(as.numeric(parse_sturgeon_v2_result(out)$score[1]), 0.8)
})

test_that("sturgeon-v2 -c and -s are only passed when filled in", {
  s <- v2_setup(writes_directory = FALSE)
  run_sturgeon_v2_predict(s$bed_dir, "m.zip", file.path(s$bin, "o1.csv"), s$v2, file_type = "bed",
                          num_cpgs = "", num_samplings = NA, log_fun = quiet)
  run_sturgeon_v2_predict(s$bed_dir, "m.zip", file.path(s$bin, "o2.csv"), s$v2, file_type = "bed",
                          num_cpgs = " 5000 ", num_samplings = 3, log_fun = quiet)
  args <- readLines(file.path(s$bin, "sturgeon-v2.args"))
  expect_false(grepl(" -c | -s ", args[1]))
  expect_match(args[2], " -c 5000 -s 3$")
})

# Test fixtures

Small, **synthetic** files for the tests in `tests/testthat/`. No file here comes from a patient or a real sequencing run. Class names are real Sturgeon class names, but every score is made up.

| File | Mimics |
|---|---|
| `sturgeon_v1_result.csv` | `sturgeon predict` output (`merged_probes_methyl_calls_general.csv`): comma-separated, `number_probes` first |
| `sturgeon_v2_result_legacy.csv` | earlier sturgeon-v2 output: one comma-separated file at the `-o` path, `probes` last |
| `sturgeon_v2_result.txt` | public sturgeon-v2 (deRidderLab) `result.txt`: tab-separated, one row per random sampling (`-s 3`) |
| `classifier_progress_old_format.txt` | a progress file from an older app version (space-separated `write.table`), which resume must skip |

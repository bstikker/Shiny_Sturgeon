# Changelog

All notable changes to Shiny-Sturgeon are documented here.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added
- `install.sh`: installs the system libraries, creates both conda environments (Sturgeon v1 needs Python < 3.10, sturgeon-v2 needs >= 3.10), runs `renv::restore()`, downloads the Sturgeon models and the Sturgeon v1 `R_scripts/` reference files from their original sources after the Sturgeon license is accepted, and writes `config_local.yml` with all tool and file paths. Downloads identify themselves as `Shiny-Sturgeon-installer`, because Zenodo refuses unnamed script requests from some networks.
- `config_local.yml` is loaded on top of `config_default.yml`, so after `install.sh` the Setup tab only needs the MinKNOW folder, barcode and output folder.
- `run_app.sh` starts the app and opens it in the browser.
- **Check installation** button on the Setup tab: lists every tool, reference file and folder with ok / problem / warning.
- [Step-by-step installation guide](docs/install-guide.md) with screenshots.

### Changed
- The log names the folder(s) being watched, and warns when the MinKNOW folder or the barcode folder doesn't exist.
- Long log lines wrap instead of being cut off.
- **Start monitoring** refuses to start when samtools, modkit or sturgeon can't be found, instead of failing in the first iteration.

## [1.0.0] - 2026-11-25

First tagged release. Validated by replaying a previous case with the public sturgeon-v2 (see "Tested with" in the README).

### Added
- Shiny GUI with a **Setup** tab (YAML config, saved per case to `<output_folder>/config.yml`) and a **Live Run** tab.
- Live monitoring of a MinKNOW `bam_pass/` folder for a chosen barcode, optionally including `unclassified/` reads.
- Per-iteration pipeline: merge and sort new BAMs, `modkit adjust-mods` and `modkit extract full`, cumulative `sturgeon inputtobed`, and classification.
- Sturgeon v1 classification with a confidence-over-time plot and a top-calls table.
- Optional Sturgeon v2 classification running side by side with v1, sharing the `inputtobed` step when both use the same probes file. Works with the public sturgeon-v2 release (tab-separated `result.txt` in the `-o` folder, one row per random sampling, averaged) and with the earlier single-file output.
- Optional `class_mapping_v2` setting: a CSV mapping v2 class names to v1 class names, so v2 classes are drawn in their v1 color. The bundled `class_mapping_v2.csv` covers five classes.
- Periodic and on-demand CNV plots (QDNAseq/DNAcopy) built from all reads processed so far.
- Resuming an interrupted run from an existing output folder.
- Configuration validation before a run starts (absolute paths, required files, v2 parameters).
- R package versions recorded with renv (`renv.lock`); conda environment files in `envs/` with the tested tool versions.
- Test suite (`tests/`, run with `Rscript tests/run_tests.R`) with synthetic fixtures and stand-in scripts for the external tools, and a GitHub Actions workflow that also runs synthetic BAMs through samtools and modkit 0.5.0.
- `tools/make_synthetic_bams.py` (synthetic MinKNOW-style BAMs with methylation tags) and `tools/make_demo_case.R` (demo case with made-up scores, used for the README screenshot).
- `CITATION.cff`, citing the Sturgeon v1 (Nature 2023) and Sturgeon V2 (medRxiv 2026) papers.

### Known limitations
- The v1↔v2 class mapping covers five classes; other v2 classes get fallback colors.
- Long class labels next to the confidence-plot lines can be cut off at the plot edge.

[Unreleased]: https://github.com/tgac-vumc/Shiny_Sturgeon/compare/v1.0.0...HEAD
[1.0.0]: https://github.com/tgac-vumc/Shiny_Sturgeon/releases/tag/v1.0.0

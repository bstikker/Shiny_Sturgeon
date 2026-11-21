# Changelog

All notable changes to Shiny-Sturgeon are documented here.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added
- Optional `class_mapping_v2` setting: a CSV mapping v2 class names to v1 class names, so v2 classes are drawn in their v1 color. The bundled `class_mapping_v2.csv` covers five classes (from the NEURON pipeline).
- Support for the public sturgeon-v2 release, which writes a tab-separated `result.txt` into the `-o` folder with one row per random sampling (scores are averaged over samplings). The earlier single-file output still works.
- Test suite (`tests/`, run with `Rscript tests/run_tests.R`) with synthetic fixtures and stand-in scripts for the external tools.
- GitHub Actions workflow running the tests, an R syntax check, and synthetic BAMs through samtools and modkit 0.5.0.
- `tools/make_synthetic_bams.py`: generates a synthetic MinKNOW-style `bam_pass/` folder with methylation tags.
- `tools/make_demo_case.R`: builds a demo case with made-up scores; README screenshot made from it.
- Conda environment files in `envs/` with the tool versions used by the NEURON pipeline.

### Changed
- `sturgeon_v2_file_type` is documented as `bed` (confirmed against sturgeon-v2's CLI) instead of a guess.

### Fixed
- README: QDNAseq is installed from Bioconductor, not CRAN; sturgeon-v2 links point to github.com/deRidderLab/sturgeon-v2.

## [1.0.0] - 2026-11-21

First tagged release.

### Added
- Shiny GUI with a **Setup** tab (YAML config, saved per case to `<output_folder>/config.yml`) and a **Live Run** tab.
- Live monitoring of a MinKNOW `bam_pass/` folder for a chosen barcode, optionally including `unclassified/` reads.
- Per-iteration pipeline: merge and sort new BAMs, `modkit adjust-mods` and `modkit extract full`, cumulative `sturgeon inputtobed`, and classification.
- Sturgeon v1 classification with a confidence-over-time plot and a top-calls table.
- Optional Sturgeon v2 classification running side by side with v1, sharing the `inputtobed` step when both use the same probes file.
- Periodic and on-demand CNV plots (QDNAseq/DNAcopy) built from all reads processed so far.
- Resuming an interrupted run from an existing output folder.
- Configuration validation before a run starts (absolute paths, required files, v2 parameters).

### Known limitations
- No v1↔v2 class-name mapping, so v2 classes use fallback colors instead of the shared clinical color scheme.

[Unreleased]: https://github.com/bstikker/Shiny_Sturgeon/compare/v1.0.0...HEAD
[1.0.0]: https://github.com/bstikker/Shiny_Sturgeon/releases/tag/v1.0.0
